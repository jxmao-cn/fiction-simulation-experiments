#!/usr/bin/env python3
"""S3 M2/M3/M4 x dt 对比：读取冻结 CSV，生成公共时刻轨迹、误差表与图。

用法：
  /home/steven/miniconda3/envs/elastogen/bin/python tools/analyze_s3_dtmesh.py \
      --run simulator/output/validation/s3-dtmesh-20261102-01 \
      --figures figures/validation/s3-dtmesh-20261102-01

约定：
  - 主量是 steps.csv 的 tip_deflection（x=1 端面平均向下 z 位移）。
  - 公共时刻 t=0,0.005,...,0.5 s；t=0 显式补 0。
  - 不插值、不做相位对齐。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MESHES = ["M2", "M3", "M4"]
DTS = ["0.005", "0.0025", "0.00125"]
COLORS = {"M2": "#1f77b4", "M3": "#d62728", "M4": "#2ca02c"}
NODES = {"M2": 525, "M3": 3321, "M4": 10309}
DT_COLORS = {"0.005": "#9467bd", "0.0025": "#ff7f0e", "0.00125": "#17becf"}
COMMON_DT = 0.005
T_END = 0.5


def project_root(start: Path) -> Path:
    p = start.resolve()
    for parent in [p, *p.parents]:
        if (parent / "fiction.project.json").exists():
            return parent
    raise SystemExit("cannot find fiction.project.json")


def read_steps(path: Path) -> tuple[list[dict], dict]:
    meta = {}
    rows = []
    with path.open() as f:
        first = f.readline()
        for part in first.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                meta[k.strip().lstrip("# ").strip()] = v.strip()
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows, meta


def case_dir(run: Path, mesh: str, dt: str) -> Path:
    return run / mesh / f"dt-{dt}"


def load_case(run: Path, mesh: str, dt: str):
    d = case_dir(run, mesh, dt)
    steps = d / "steps.csv"
    if not steps.exists():
        return None
    rows, meta = read_steps(steps)
    meta["dt"] = dt
    t = np.array([float(r["time"]) for r in rows])
    u = np.array([float(r["tip_deflection"]) for r in rows])
    return {"t": t, "u": u, "rows": rows, "meta": meta, "dir": d}


def common_curve(case, times: np.ndarray):
    """按整数步索引取公共时刻；返回 (u, max_time_error)。"""
    dt = float(case["meta"].get("dt", "nan"))
    if not math.isfinite(dt) or dt <= 0:
        # 从 steps 的第一行反推
        dt = float(case["t"][0])
    u = np.zeros_like(times)
    max_dt_err = 0.0
    for i, tt in enumerate(times):
        if tt == 0.0:
            u[i] = 0.0
            continue
        k = int(round(tt / dt))
        idx = k - 1
        if idx < 0 or idx >= len(case["u"]):
            return None, None
        max_dt_err = max(max_dt_err, abs(case["t"][idx] - tt))
        u[i] = case["u"][idx]
    return u, max_dt_err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--figures", required=True)
    args = ap.parse_args()

    root = project_root(Path.cwd())
    run = (root / args.run).resolve()
    fig = (root / args.figures).resolve()
    fig.mkdir(parents=True, exist_ok=True)

    times = np.arange(0.0, T_END + 1e-12, COMMON_DT)
    cases = {}
    for mesh in MESHES:
        for dt in DTS:
            c = load_case(run, mesh, dt)
            if c is not None:
                u, err = common_curve(c, times)
                if u is not None:
                    c["common_u"] = u
                    c["common_time_error"] = err
                    cases[(mesh, dt)] = c

    if not cases:
        raise SystemExit(f"no complete cases under {run}")

    # 逐病例安全/求解状态
    status_rows = []
    for (mesh, dt), c in sorted(cases.items()):
        rows = c["rows"]
        all_conv = all(r.get("converged") == "1" for r in rows)
        finite = all(math.isfinite(float(r[k])) for r in rows for k in
                     ("tip_deflection", "min_J", "dbc_error"))
        min_j = min(float(r["min_J"]) for r in rows)
        max_dbc = max(float(r["dbc_error"]) for r in rows)
        tip = np.array([float(r["tip_deflection"]) for r in rows])
        t = np.array([float(r["time"]) for r in rows])
        status_rows.append({"mesh": mesh, "dt": float(dt), "steps": len(rows),
                            "all_converged": all_conv, "finite": finite,
                            "min_J": min_j, "max_dbc_error": max_dbc,
                            "peak_tip_m": float(np.max(tip)),
                            "peak_time_s": float(t[int(np.argmax(tip))])})
    with (fig / "case-status.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(status_rows[0].keys()))
        w.writeheader()
        w.writerows(status_rows)

    # 逐网格时间差
    err_rows = []
    for mesh in MESHES:
        if not all((mesh, dt) in cases for dt in DTS):
            continue
        u1 = cases[(mesh, "0.005")]["common_u"]
        u2 = cases[(mesh, "0.0025")]["common_u"]
        u3 = cases[(mesh, "0.00125")]["common_u"]
        e1 = float(np.max(np.abs(u1 - u2)))
        e2 = float(np.max(np.abs(u2 - u3)))
        rms1 = float(np.sqrt(np.mean((u1 - u2) ** 2)))
        rms2 = float(np.sqrt(np.mean((u2 - u3) ** 2)))
        p = math.log2(e1 / e2) if e2 > 0 and e1 > 0 else float("nan")
        peak = float(np.max(np.abs(u3)))
        err_rows.append({
            "mesh": mesh,
            "e_max_0.005_vs_0.0025_m": e1,
            "e_max_0.0025_vs_0.00125_m": e2,
            "e_max_0.005_vs_0.0025_rel": e1 / peak if peak > 0 else float("nan"),
            "e_max_0.0025_vs_0.00125_rel": e2 / peak if peak > 0 else float("nan"),
            "rms_0.005_vs_0.0025_m": rms1,
            "rms_0.0025_vs_0.00125_m": rms2,
            "peak_finest_m": peak,
            "observed_order_p": p,
            "monotone_decrease": bool(e2 < e1),
        })

    with (fig / "error-summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(err_rows[0].keys()) if err_rows else ["mesh"])
        w.writeheader()
        w.writerows(err_rows)

    # 公共时刻轨迹
    with (fig / "tip-common-times.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_s"] + [f"{m}_dt{d}" for m in MESHES for d in DTS])
        for i, tt in enumerate(times):
            row = [f"{tt:.6f}"]
            for m in MESHES:
                for d in DTS:
                    c = cases.get((m, d))
                    row.append(f"{c['common_u'][i]:.12g}" if c else "")
            w.writerow(row)

    # 同 dt 的网格差
    mesh_rows = []
    for dt in DTS:
        if not all((m, dt) in cases for m in MESHES):
            continue
        a = cases[("M2", dt)]["common_u"]
        b = cases[("M3", dt)]["common_u"]
        c = cases[("M4", dt)]["common_u"]
        mesh_rows.append({
            "dt": float(dt),
            "max_M2_M3_m": float(np.max(np.abs(a - b))),
            "max_M3_M4_m": float(np.max(np.abs(b - c))),
            "rms_M2_M3_m": float(np.sqrt(np.mean((a - b) ** 2))),
            "rms_M3_M4_m": float(np.sqrt(np.mean((b - c) ** 2))),
        })
    with (fig / "mesh-difference.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mesh_rows[0].keys()) if mesh_rows else ["dt"])
        w.writeheader()
        w.writerows(mesh_rows)

    # 图 1：三网格的 dt 轨迹（共享 y 轴）
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})
    fig1, axes = plt.subplots(1, 3, figsize=(12.5, 3.6), sharey=True)
    for ax, mesh in zip(axes, MESHES):
        for dt in DTS:
            c = cases.get((mesh, dt))
            if c is None:
                continue
            ax.plot(c["t"], c["u"] * 1000, lw=1.3, color=DT_COLORS[dt],
                    label=f"dt={dt} s")
        ax.set_title(f"{mesh}  ({NODES[mesh]} nodes)")
        ax.set_xlabel("time [s]")
        ax.set_xlim(0, T_END)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("tip deflection [mm]")
    fig1.tight_layout()
    fig1.savefig(fig / "s3-tip-trajectories.png", dpi=180)
    fig1.savefig(fig / "s3-tip-trajectories.svg")
    plt.close(fig1)

    # 图 2：误差随 dt（log-log），并标观测阶
    fig2, ax = plt.subplots(figsize=(5.4, 3.8))
    for r in err_rows:
        xs = np.array([0.005, 0.0025])
        ys = np.array([r["e_max_0.005_vs_0.0025_m"], r["e_max_0.0025_vs_0.00125_m"]]) * 1000
        ax.plot(xs, ys, "o-", color=COLORS[r["mesh"]], label=f"{r['mesh']} (p={r['observed_order_p']:.2f})")
        if r["mesh"] == "M2":
            ax.plot(xs, ys[0] * (xs / xs[0]), "k--", lw=0.9, label="slope 1 guide")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.005, 0.0025])
    ax.set_xticklabels(["0.005", "0.0025"])
    ax.set_xlabel("dt [s]")
    ax.set_ylabel("max tip difference to next finer dt [mm]")
    ax.set_title("Time-step difference decreases with finer dt")
    ax.legend(frameon=False)
    fig2.tight_layout()
    fig2.savefig(fig / "s3-time-error-convergence.png", dpi=180)
    fig2.savefig(fig / "s3-time-error-convergence.svg")
    plt.close(fig2)

    # 图 2b：按最细 dt 峰值归一化的误差
    fig2b, ax = plt.subplots(figsize=(5.4, 3.8))
    for r in err_rows:
        xs = np.array([0.005, 0.0025])
        ys = np.array([r["e_max_0.005_vs_0.0025_rel"], r["e_max_0.0025_vs_0.00125_rel"]]) * 100
        ax.plot(xs, ys, "o-", color=COLORS[r["mesh"]], label=f"{r['mesh']} (p={r['observed_order_p']:.2f})")
        if r["mesh"] == "M2":
            ax.plot(xs, ys[0] * (xs / xs[0]), "k--", lw=0.9, label="slope 1 guide")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dt [s]")
    ax.set_ylabel("max difference / peak of finest dt [%]")
    ax.set_title("Normalized time-step difference")
    ax.legend(frameon=False)
    fig2b.tight_layout()
    fig2b.savefig(fig / "s3-time-error-relative.png", dpi=180)
    fig2b.savefig(fig / "s3-time-error-relative.svg")
    plt.close(fig2b)

    # 图 3：同 dt 的网格比较
    fig3, ax = plt.subplots(figsize=(6.0, 3.8))
    dt_plot = "0.00125"
    for mesh in MESHES:
        c = cases.get((mesh, dt_plot))
        if c is None:
            continue
        ax.plot(c["t"], c["u"] * 1000, lw=1.4, color=COLORS[mesh], label=mesh)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("tip deflection [mm]")
    ax.set_title(f"Mesh comparison at dt={dt_plot} s")
    ax.set_xlim(0, T_END)
    ax.legend(frameon=False)
    fig3.tight_layout()
    fig3.savefig(fig / "mesh-comparison.png", dpi=180)
    fig3.savefig(fig / "mesh-comparison.svg")
    plt.close(fig3)

    # 图 4：每步成本（若 time.log 有 Elapsed，另行解析；此处用 steps.csv wall 求和）
    cost_rows = []
    for mesh in MESHES:
        for dt in DTS:
            c = cases.get((mesh, dt))
            if c is None:
                continue
            wall = sum(float(r["wall_seconds"]) for r in c["rows"])
            cost_rows.append({"mesh": mesh, "dt": float(dt), "steps": len(c["rows"]),
                              "wall_s": wall, "ms_per_step": 1000 * wall / len(c["rows"])})
    with (fig / "cost-summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mesh", "dt", "steps", "wall_s", "ms_per_step"])
        w.writeheader()
        w.writerows(cost_rows)

    fig4, ax = plt.subplots(figsize=(5.0, 3.6))
    labels = [f"{r['mesh']}\n{r['steps']} steps" for r in cost_rows]
    vals = [r["ms_per_step"] for r in cost_rows]
    ax.bar(range(len(vals)), vals, color=[COLORS[r["mesh"]] for r in cost_rows])
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("wall time per step [ms]")
    ax.set_title("Single-run cost (includes I/O and diagnostics)")
    fig4.tight_layout()
    fig4.savefig(fig / "s3-cost-per-step.png", dpi=180)
    fig4.savefig(fig / "s3-cost-per-step.svg")
    plt.close(fig4)

    # 图 5：峰值尖端位移随网格（同 dt 对照）
    fig5, ax = plt.subplots(figsize=(5.2, 3.6))
    for dt in DTS:
        xs, ys = [], []
        for mi, mesh in enumerate(MESHES):
            c = cases.get((mesh, dt))
            if c is None: continue
            xs.append(mi); ys.append(float(np.max(c["u"])) * 1000)
        if xs:
            ax.plot(xs, ys, "o-", color=DT_COLORS[dt], label=f"dt={dt} s")
    ax.set_xticks(range(len(MESHES)))
    ax.set_xticklabels([f"{m}\n({NODES[m]} nodes)" for m in MESHES])
    ax.set_ylabel("peak tip deflection [mm]")
    ax.set_title("Spatial effect: peak amplitude vs mesh")
    ax.legend(frameon=False, fontsize=8)
    fig5.tight_layout()
    fig5.savefig(fig / "mesh-peak-amplitude.png", dpi=180)
    fig5.savefig(fig / "mesh-peak-amplitude.svg")
    plt.close(fig5)

    # 图 6：同 dt 下相邻网格的最大差
    if mesh_rows:
        fig6, ax = plt.subplots(figsize=(5.2, 3.6))
        ds = [r["dt"] for r in mesh_rows]
        ax.plot(ds, [r["max_M2_M3_m"] * 1000 for r in mesh_rows], "o-", color="#1f77b4", label="M2 vs M3")
        ax.plot(ds, [r["max_M3_M4_m"] * 1000 for r in mesh_rows], "s-", color="#2ca02c", label="M3 vs M4")
        ax.set_xscale("log")
        ax.set_xticks(ds)
        ax.set_xticklabels([str(d) for d in ds])
        ax.set_xlabel("dt [s]")
        ax.set_ylabel("max tip difference [mm]")
        ax.set_title("Spatial effect at fixed dt")
        ax.legend(frameon=False, fontsize=8)
        fig6.tight_layout()
        fig6.savefig(fig / "mesh-difference.png", dpi=180)
        fig6.savefig(fig / "mesh-difference.svg")
        plt.close(fig6)

    summary = {
        "run": str(run.relative_to(root)),
        "common_times": len(times),
        "cases_found": len(cases),
        "case_status": status_rows,
        "time_errors": err_rows,
        "mesh_differences": mesh_rows,
        "common_time_max_error_s": max((c["common_time_error"] for c in cases.values()), default=None),
    }
    with (fig / "analysis-summary.json").open("w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
