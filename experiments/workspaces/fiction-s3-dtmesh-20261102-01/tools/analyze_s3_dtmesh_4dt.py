#!/usr/bin/env python3
"""S3 M2/M3/M4 x 4 dt (0.005/0.0025/0.00125/0.000625) 对比。

主 run 提供前三档，ref run 提供 0.000625。纯读取冻结 CSV，不插值。
"""
import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MESHES = ["M2", "M3", "M4"]
DTS = [0.005, 0.0025, 0.00125, 0.000625]
COLORS = {"M2": "#1f77b4", "M3": "#d62728", "M4": "#2ca02c"}
DT_COLORS = {0.005: "#9467bd", 0.0025: "#ff7f0e", 0.00125: "#17becf", 0.000625: "#8c564b"}
T_END = 0.5
SAMPLE = 0.005


def load_case(run, mesh, dt):
    p = Path(run) / mesh / f"dt-{dt:g}" / "steps.csv"
    if not p.exists():
        return None
    lines = [ln for ln in open(p) if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    by = {int(r["step"]): r for r in rows}
    times = [j * SAMPLE for j in range(int(round(T_END / SAMPLE)) + 1)]
    u = [0.0]
    for t in times[1:]:
        step = int(round(t / dt))
        if step not in by:
            return None
        if abs(float(by[step]["time"]) - t) > 1e-12:
            return None
        u.append(float(by[step]["tip_deflection"]))
    peak = max(rows, key=lambda r: float(r["tip_deflection"]))
    wall = sum(float(r["wall_seconds"]) for r in rows)
    rss = None
    tl = p.parent / "time.log"
    if tl.exists():
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", tl.read_text())
        if m:
            rss = int(m.group(1))
    return {
        "rows": rows, "u": np.array(u), "times": np.array(times),
        "peak_m": float(peak["tip_deflection"]), "peak_t": float(peak["time"]),
        "tip_end_m": float(rows[-1]["tip_deflection"]),
        "max_speed": max(float(r["max_speed"]) for r in rows),
        "min_J": min(float(r["min_J"]) for r in rows),
        "max_dbc": max(float(r["dbc_error"]) for r in rows),
        "wall_s": wall, "ms_per_step": 1000 * wall / len(rows),
        "peak_rss_kib": rss, "steps": len(rows),
        "all_converged": all(r["converged"] == "1" for r in rows),
        "max_iter": max(int(r["iterations"]) for r in rows),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cases = {}
    for m in MESHES:
        for dt in DTS:
            run = args.main if dt >= 0.00125 else args.ref
            c = load_case(run, m, dt)
            if c is not None:
                cases[(m, dt)] = c

    # case table
    status = []
    for (m, dt), c in sorted(cases.items(), key=lambda kv: (kv[0][0], -kv[0][1])):
        status.append({"mesh": m, "dt": dt, "steps": c["steps"],
                       "all_converged": c["all_converged"],
                       "max_iter": c["max_iter"], "min_J": c["min_J"],
                       "max_dbc": c["max_dbc"], "peak_mm": c["peak_m"] * 1e3,
                       "peak_t_s": c["peak_t"], "tip_0.5_mm": c["tip_end_m"] * 1e3,
                       "max_speed": c["max_speed"], "wall_s": c["wall_s"],
                       "ms_per_step": c["ms_per_step"],
                       "peak_rss_kib": c["peak_rss_kib"]})
    with open(out / "case-values.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(status[0].keys()))
        w.writeheader(); w.writerows(status)

    # time errors
    errors = []
    for m in MESHES:
        if not all((m, dt) in cases for dt in DTS):
            continue
        us = [cases[(m, dt)]["u"] for dt in DTS]
        e = [float(np.max(np.abs(us[i] - us[i + 1]))) for i in range(3)]
        rms = [float(np.sqrt(np.mean((us[i] - us[i + 1]) ** 2))) for i in range(3)]
        p12 = math.log2(e[0] / e[1]) if e[1] > 0 else float("nan")
        p23 = math.log2(e[1] / e[2]) if e[2] > 0 else float("nan")
        errors.append({"mesh": m,
                       "e1_0.005_to_0.0025_m": e[0],
                       "e2_0.0025_to_0.00125_m": e[1],
                       "e3_0.00125_to_0.000625_m": e[2],
                       "ratio_e1_e2": e[0] / e[1] if e[1] else None,
                       "ratio_e2_e3": e[1] / e[2] if e[2] else None,
                       "p12": p12, "p23": p23,
                       "rms1_m": rms[0], "rms2_m": rms[1], "rms3_m": rms[2],
                       "rms_ratio_e1_e2": rms[0] / rms[1] if rms[1] else None,
                       "rms_ratio_e2_e3": rms[1] / rms[2] if rms[2] else None,
                       "monotone_3": bool(e[2] < e[1] < e[0]),
                       "peak_finest_mm": cases[(m, DTS[-1])]["peak_m"] * 1e3})
    with open(out / "time-errors.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(errors[0].keys()))
        w.writeheader(); w.writerows(errors)

    summary = {"main_run": args.main, "ref_run": args.ref, "dts": DTS,
               "cases": status, "time_errors": errors}
    with open(out / "analysis-summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # common trajectories
    with open(out / "tip-common-times.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s"] + [f"{m}_dt{dt:g}" for m in MESHES for dt in DTS])
        for i, t in enumerate(cases[(MESHES[0], DTS[0])]["times"]):
            row = [f"{t:.6f}"]
            for m in MESHES:
                for dt in DTS:
                    row.append(f"{cases[(m, dt)]['u'][i]:.12g}" if (m, dt) in cases else "")
            w.writerow(row)

    # plots
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.6), sharey=True)
    for ax, m in zip(axes, MESHES):
        for dt in DTS:
            c = cases.get((m, dt))
            if c:
                ax.plot(c["times"], c["u"] * 1000, lw=1.2, color=DT_COLORS[dt], label=f"dt={dt:g} s")
        ax.set_title(m); ax.set_xlabel("time [s]"); ax.set_xlim(0, T_END)
        ax.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("tip deflection [mm]")
    fig.tight_layout(); fig.savefig(out / "tip-trajectories.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for r in errors:
        xs = np.array([0.005, 0.0025, 0.00125])
        ys = np.array([r["e1_0.005_to_0.0025_m"], r["e2_0.0025_to_0.00125_m"],
                       r["e3_0.00125_to_0.000625_m"]]) * 1000
        ax.plot(xs, ys, "o-", color=COLORS[r["mesh"]],
                label=f"{r['mesh']} (p12={r['p12']:.2f}, p23={r['p23']:.2f})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("dt [s]"); ax.set_ylabel("max tip difference to next finer dt [mm]")
    ax.set_title("Time-step difference over three halvings")
    ax.legend(frameon=False, fontsize=7); fig.tight_layout()
    fig.savefig(out / "time-error-convergence.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for m in MESHES:
        c = cases[(m, DTS[-1])]
        ax.plot([m], [c["peak_m"] * 1e3], "o", color=COLORS[m], label=f"{m} peak")
        ax.annotate(f"{c['peak_m']*1e3:.1f} mm\n@ {c['peak_t']:.3f} s",
                    (m, c["peak_m"] * 1e3), textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_ylabel("peak tip deflection [mm]"); ax.set_title("Peak amplitude, dt=0.000625 s")
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(out / "mesh-peak-amplitude.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    for m in MESHES:
        xs = DTS
        ys = [cases[(m, dt)]["ms_per_step"] for dt in DTS if (m, dt) in cases]
        ax.plot(xs, ys, "o-", color=COLORS[m], label=m)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("dt [s]"); ax.set_ylabel("wall time per step [ms]")
    ax.set_title("Cost per step"); ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(out / "cost-per-step.png", dpi=180); plt.close(fig)

    print(json.dumps({"cases": len(status), "errors": errors}, indent=2))


if __name__ == "__main__":
    main()
