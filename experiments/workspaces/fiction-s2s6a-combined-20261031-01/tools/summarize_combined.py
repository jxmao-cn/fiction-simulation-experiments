#!/usr/bin/env python3
"""v0.3 停止判据候选确认汇总（run-id s2s6a-combined-20261031-01）。

输入（由 run_combined.sh 产生）：
  OUT/s2/M3/{M3.csv,M3-state.csv,run.log,time_rss.txt,exit_code.txt,trace/iterations.csv}
  OUT/s2/M6/{M6.csv,M6-state.csv,...}
  OUT/s6a/body1/{body1_steps.csv,body1_state.csv,...}
  OUT/s6a/body4/{body4_steps.csv,body4_state.csv,block_diagonal.json,trajectory_compare.csv,summary.json,...}

输出：
  OUT/<case>/config.json, OUT/<case>/summary.json
  OUT/cases.csv, OUT/judgment.json

预登记门槛（与 protocol-v0.3 草案 §4 一致，不放宽）：
  - 物理解不变：末态最大节点差 < 1e-8 m
  - S2：converged=1 且 stopReason=combined；相对误差保持 M3 25.94% / M6 4.26%
  - S6a：10/10 步 converged；每步停止时 free_gradient/gradScale < 1e-6
  - S6a 4 体：块对角跨体 triplet/entry = 0；4 体 vs 4 次独立单体 < 1e-8 m
"""
import argparse
import csv
import json
import math
import os

RUN_ID = "s2s6a-combined-20261031-01"
POS_TOL = 1e-8
GRAD_RATIO_TOL = 1e-6
ERR_TOL = 1e-6

ARCHIVE = "/home/steven/mjx/fiction/simulator/output/validation"
STOPCRIT = os.path.join(ARCHIVE, "stopcriterion-20261031-01")


def read_rows(path):
    with open(path) as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    return list(csv.DictReader(lines))


def read_state(path, cols=("x", "y", "z")):
    out = []
    for r in read_rows(path):
        out.append(tuple(float(r[c]) for c in cols))
    return out


def max_node_diff(a, b):
    if len(a) != len(b):
        raise ValueError(f"节点数不一致: {len(a)} vs {len(b)}")
    return max(math.dist(p, q) for p, q in zip(a, b))


def parse_time_rss(path):
    wall = rss = None
    if not os.path.exists(path):
        return wall, rss
    for line in open(path, errors="ignore"):
        s = line.strip()
        if s.startswith("Elapsed (wall clock) time"):
            wall = s.split("):", 1)[-1].strip()
        elif s.startswith("Maximum resident set size"):
            rss = int(s.split(":", 1)[1].strip())
    return wall, rss


def scan_finite(paths):
    """对 CSV 数值单元做有限性扫描，返回非有限值计数。"""
    bad = 0
    for p in paths:
        if not os.path.exists(p):
            bad += 1
            continue
        with open(p) as f:
            for line in f:
                if line.lstrip().startswith("#") or not line.strip():
                    continue
                for cell in line.rstrip("\n").split(","):
                    try:
                        v = float(cell)
                    except ValueError:
                        continue
                    if not math.isfinite(v):
                        bad += 1
    return bad


def s2_case(out, mesh_name, mesh_path, csv_name, ref_state, ref_err, tip_name,
            reference_label):
    d = os.path.join(out, "s2", mesh_name)
    row = read_rows(os.path.join(d, csv_name))[0]
    cur = read_state(os.path.join(d, f"{mesh_name}-state.csv"))
    ref = read_state(ref_state, ("x", "y", "z"))
    diff = max_node_diff(cur, ref)

    # 同判据复现性交叉检查（stopcriterion 的 combined 运行）
    repro_diff = None
    repro_path = os.path.join(STOPCRIT, "s2", "M3", "combined", "M3-state.csv")
    if mesh_name == "M3" and os.path.exists(repro_path):
        repro_diff = max_node_diff(cur, read_state(repro_path))

    grad_scale = float(row["grad_scale"])
    rel_grad = float(row["free_gradient_inf"]) / grad_scale
    conv = int(float(row["converged"]))
    stop_reason = row["stop_reason"]
    iters = int(float(row["iterations"]))
    ls = int(float(row["line_search_iterations"]))
    err = float(row["relative_error"])
    wall, rss = parse_time_rss(os.path.join(d, "time_rss.txt"))
    bad = scan_finite([os.path.join(d, csv_name),
                       os.path.join(d, f"{mesh_name}-state.csv")])
    exit_code = int(open(os.path.join(d, "exit_code.txt")).read().strip())

    gates = {
        "physical_unchanged_lt_1e8m": diff < POS_TOL,
        "converged": conv == 1,
        "stop_reason_combined": stop_reason == "combined",
        "no_premature_stop_grad_ratio_lt_1e6": rel_grad < GRAD_RATIO_TOL,
        "accuracy_error_unchanged": abs(err - ref_err) <= ERR_TOL,
        "all_finite": bad == 0,
        "exit_code_zero": exit_code == 0,
    }
    if repro_diff is not None:
        gates["reproduces_stopcriterion_combined"] = repro_diff == 0.0

    cfg = {
        "run_id": RUN_ID,
        "case": f"S2 {mesh_name} static cantilever (NeoHookean tet, E=5.90e8, nu=0.3, rho=1000)",
        "mesh": mesh_path,
        "stop_mode": "combined",
        "tol": 1e-10,  # cantilever_check 硬编码（静态梁，M3/M6 同一程序）
        "max_iter": 100,
        "static_solve": True,
        "grad_tol_rel": float(row["grad_tol_rel"]),
        "grad_scale_n": grad_scale,
        "actual_step_tol_m": float(row["actual_step_tol"]),
        "actual_step_count": int(float(row["actual_step_count"])),
        "reference_state": reference_label,
        "reference_relative_error": ref_err,
        "preregistered_gates": {
            "state_max_node_diff_m": f"< {POS_TOL}",
            "relative_error_abs_change": f"<= {ERR_TOL}",
        },
    }
    res = {
        "case": f"s2/{mesh_name}",
        "mesh_nodes": int(float(row["nodes"])),
        "mesh_tets": int(float(row["tets"])),
        "stop_mode": row["stop_mode"],
        "stop_reason": stop_reason,
        "iterations": iters,
        "line_search_iterations": ls,
        "converged": conv,
        "max_iter_hit": int(float(row["max_iter_hit"])),
        "line_search_failed": int(float(row["line_search_failed"])),
        "solver": row["solver"],
        "search_direction_rate": float(row["search_direction_rate"]),
        "free_gradient_inf_n": float(row["free_gradient_inf"]),
        "free_gradient_relative": rel_grad,
        "min_J": float(row["min_J"]),
        "max_J": float(row["max_J"]),
        "dbc_satisfied": int(float(row["dbc_satisfied"])),
        "dbc_total": int(float(row["dbc_total"])),
        "dbc_max_error": float(row["dbc_max_error"]),
        "numerical_pass_flag": int(float(row["numerical_pass"])),
        "tip_deflection_m": float(row[tip_name]),
        "relative_error": err,
        "reference_relative_error": ref_err,
        "reference_state": reference_label,
        "max_node_diff_vs_reference_m": diff,
        "max_node_diff_vs_stopcriterion_combined_m": repro_diff,
        "solver_wall_seconds": float(row["solve_wall_seconds"]),
        "process_wall": wall,
        "max_rss_kib": rss,
        "non_finite_cells": bad,
        "exit_code": exit_code,
        "gates": gates,
        "pass_all_gates": all(gates.values()),
    }
    write_json(os.path.join(d, "config.json"), cfg)
    write_json(os.path.join(d, "summary.json"), res)
    return res


def s6a_case(out, body, mesh_path, ref_state, grad_scale_expected):
    d = os.path.join(out, "s6a", body)
    steps_path = os.path.join(d, f"{body}_steps.csv")
    rows = read_rows(steps_path)
    state_path = os.path.join(d, f"{body}_state.csv")
    cur = read_state(state_path, ("global_x", "global_y", "global_z"))
    ref = read_state(ref_state, ("global_x", "global_y", "global_z"))
    diff = max_node_diff(cur, ref)

    grad_scale = float(grad_scale_expected)
    rel = [float(r["free_gradient_inf"]) / grad_scale for r in rows]
    iters = [int(float(r["iterations"])) for r in rows]
    lss = [int(float(r["line_search_iterations"])) for r in rows]
    reasons = [r["stop_reason"] for r in rows]
    conv = [int(float(r["converged"])) for r in rows]
    wall, rss = parse_time_rss(os.path.join(d, "time_rss.txt"))
    bad = scan_finite([steps_path, state_path])
    exit_code = int(open(os.path.join(d, "exit_code.txt")).read().strip())

    extra = {}
    if body == "body4":
        bd = json.load(open(os.path.join(d, "block_diagonal.json")))
        summ = json.load(open(os.path.join(d, "program_summary.json")))
        cmp = summ["independent_compare"]
        extra = {
            "block_diagonal": {
                "cross_body_triplets_final": bd["body4_final"]["cross_body_triplets"],
                "assembled_cross_body_entries_final": bd["body4_final"]["assembled_cross_body_entries"],
                "tets_spanning_bodies": bd["body4_final"]["tets_spanning_bodies"],
                "block_pass": bd["body4_final"]["pass"],
            },
            "independent_compare": cmp,
            "no_contact": summ["no_contact"],
            "program_summary_pass": summ["pass"],
            "program_stop_criterion": summ.get("stop_criterion"),
            "cost": summ["cost"],
        }
        # 独立单体逐步收敛情况
        indep_rows = []
        for b in range(4):
            p = os.path.join(d, f"independent_body{b}_steps.csv")
            ir = read_rows(p)
            indep_rows.append({
                "body": b,
                "all_converged": all(int(float(r["converged"])) == 1 for r in ir),
                "iterations": [int(float(r["iterations"])) for r in ir],
                "stop_reasons": [r["stop_reason"] for r in ir],
            })
        extra["independent_runs"] = indep_rows

    gates = {
        "all_10_steps_converged": len(conv) == 10 and all(c == 1 for c in conv),
        "physical_unchanged_lt_1e8m": diff < POS_TOL,
        "no_maxiter_stop": all(r != "maxIter" for r in reasons),
        "no_premature_stop_grad_ratio_lt_1e6": max(rel) < GRAD_RATIO_TOL,
        "all_finite": bad == 0,
        "exit_code_zero": exit_code == 0,
    }
    if body == "body4":
        gates.update({
            "block_diagonal_zero_cross_body": (
                extra["block_diagonal"]["cross_body_triplets_final"] == 0
                and extra["block_diagonal"]["assembled_cross_body_entries_final"] == 0
                and extra["block_diagonal"]["tets_spanning_bodies"] == 0),
            "independent_lt_1e8m": extra["independent_compare"]["max_position_diff_m"] < POS_TOL,
        })

    cfg = {
        "run_id": RUN_ID,
        "case": f"S6a {body} (sphere1K, NH E=1e5 nu=0.3 rho=1000, gravity, dt=0.005, 10 steps)",
        "mesh": mesh_path,
        "stop_mode": "combined",
        "tol": 1e-8,
        "max_iter": 100,
        "steps": 10,
        "grad_tol_rel": 1e-6,
        "grad_scale_n": grad_scale,
        "actual_step_tol_m": 1e-12,
        "actual_step_count": 2,
        "reference_state": ref_state,
        "preregistered_gates": {
            "per_step_converged": "10/10",
            "state_max_node_diff_m": f"< {POS_TOL}",
            "block_diagonal_cross_body": "0",
            "independent_vs_4body_m": f"< {POS_TOL}",
        },
    }
    res = {
        "case": f"s6a/{body}",
        "stop_mode": "combined",
        "stop_reasons": reasons,
        "iterations_per_step": iters,
        "total_iterations": sum(iters),
        "line_search_per_step": lss,
        "total_line_search": sum(lss),
        "state_change_inf_per_step": [float(r["state_change_inf"]) for r in rows],
        "free_gradient_relative_per_step": rel,
        "max_free_gradient_relative_at_stop": max(rel),
        "converged_per_step": conv,
        "min_J_per_step": [float(r["minJ"]) for r in rows],
        "max_disp_m": float(rows[-1]["max_disp"]),
        "max_node_diff_vs_reference_m": diff,
        "reference_state": ref_state,
        "process_wall": wall,
        "max_rss_kib": rss,
        "non_finite_cells": bad,
        "exit_code": exit_code,
        "gates": gates,
        "pass_all_gates": all(gates.values()),
        **extra,
    }
    write_json(os.path.join(d, "config.json"), cfg)
    write_json(os.path.join(d, "summary.json"), res)
    return res


def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = args.out

    m3 = s2_case(
        out, "M3",
        os.path.join(ARCHIVE, "s2-cantilever-20261031-01/meshes/M3.msh"),
        "M3.csv",
        os.path.join(ARCHIVE, "s2-cantilever-20261031-01/M3-state.csv"),
        0.2594127814953327,
        "tip_area_weighted",
        "s2-cantilever-20261031-01/M3-state.csv (v0.2 strict, tol=1e-10, step 判据 maxIter 末态)")
    m6 = s2_case(
        out, "M6",
        os.path.join(ARCHIVE, "s2-extension-20261031-01/meshes/M6.msh"),
        "M6.csv",
        os.path.join(ARCHIVE, "s2-extension-20261031-01/M6-tol1e8/state.csv"),
        0.042597010571187602,
        "tip_area_weighted",
        "s2-extension-20261031-01/M6-tol1e8/state.csv (v0.2 tol=1e-8)")
    b1 = s6a_case(
        out, "body1",
        "/home/steven/mjx/Archive/IPC/input/tetMeshes/sphere1K.msh",
        os.path.join(ARCHIVE, "s6a-multibody-20261031-01/body1_state.csv"),
        5086.259120998311)
    b4 = s6a_case(
        out, "body4",
        "/home/steven/mjx/Archive/IPC/input/tetMeshes/sphere1K.msh",
        os.path.join(ARCHIVE, "s6a-multibody-20261031-01/body4_state.csv"),
        20345.036483993244)

    cases = [m3, m6, b1, b4]
    with open(os.path.join(out, "cases.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "stop_mode", "stop_reason", "iterations", "line_search_iterations",
                    "converged", "total_iterations", "steps_converged", "max_free_gradient_relative_at_stop",
                    "search_direction_rate", "max_node_diff_vs_reference_m", "relative_error",
                    "relative_error_reference", "max_disp_or_tip_m", "block_diagonal_cross_body",
                    "independent_vs_4body_m", "non_finite_cells", "process_wall",
                    "max_rss_kib", "exit_code", "pass_all_gates"])
        for c in cases:
            if c["case"].startswith("s2/"):
                w.writerow([c["case"], c["stop_mode"], c["stop_reason"], c["iterations"],
                            c["line_search_iterations"], c["converged"], "", "",
                            f"{c['free_gradient_relative']:.6e}", f"{c['search_direction_rate']:.6e}",
                            f"{c['max_node_diff_vs_reference_m']:.6e}", f"{c['relative_error']:.9f}",
                            f"{c['reference_relative_error']:.9f}", f"{c['tip_deflection_m']:.12f}",
                            "", "", c["non_finite_cells"], c["process_wall"], c["max_rss_kib"],
                            c["exit_code"], c["pass_all_gates"]])
            else:
                conv_n = sum(c["converged_per_step"])
                bd = c.get("block_diagonal", {})
                indep = c.get("independent_compare", {}).get("max_position_diff_m", "")
                w.writerow([c["case"], c["stop_mode"], "|".join(c["stop_reasons"]), "",
                            0, "", c["total_iterations"], f"{conv_n}/10",
                            f"{c['max_free_gradient_relative_at_stop']:.6e}", "",
                            f"{c['max_node_diff_vs_reference_m']:.6e}", "", "",
                            f"{c['max_disp_m']:.12f}",
                            bd.get("cross_body_triplets_final", ""),
                            indep, c["non_finite_cells"], c["process_wall"], c["max_rss_kib"],
                            c["exit_code"], c["pass_all_gates"]])

    judgment = {
        "run_id": RUN_ID,
        "nature": "候选判据确认（不是协议通过）；S2 精度门槛不变；生产代码与既有 run-id 未改动",
        "preregistered_thresholds": {
            "state_max_node_diff_m": POS_TOL,
            "free_gradient_relative_at_stop": GRAD_RATIO_TOL,
            "s2_relative_error_abs_change": ERR_TOL,
            "s6a_steps_converged": "10/10",
            "block_diagonal_cross_body": 0,
            "independent_vs_4body_m": POS_TOL,
        },
        "cases": {
            "s2/M3": {
                "stop_reason": m3["stop_reason"], "iterations": m3["iterations"],
                "max_node_diff_vs_reference_m": m3["max_node_diff_vs_reference_m"],
                "reproduces_stopcriterion_combined_m": m3["max_node_diff_vs_stopcriterion_combined_m"],
                "relative_error": m3["relative_error"],
                "free_gradient_relative": m3["free_gradient_relative"],
                "gates": m3["gates"], "pass": m3["pass_all_gates"],
            },
            "s2/M6": {
                "stop_reason": m6["stop_reason"], "iterations": m6["iterations"],
                "max_node_diff_vs_reference_m": m6["max_node_diff_vs_reference_m"],
                "relative_error": m6["relative_error"],
                "free_gradient_relative": m6["free_gradient_relative"],
                "gates": m6["gates"], "pass": m6["pass_all_gates"],
            },
            "s6a/body1": {
                "stop_reasons": b1["stop_reasons"],
                "iterations_per_step": b1["iterations_per_step"],
                "total_iterations": b1["total_iterations"],
                "max_node_diff_vs_reference_m": b1["max_node_diff_vs_reference_m"],
                "max_free_gradient_relative_at_stop": b1["max_free_gradient_relative_at_stop"],
                "gates": b1["gates"], "pass": b1["pass_all_gates"],
            },
            "s6a/body4": {
                "stop_reasons": b4["stop_reasons"],
                "iterations_per_step": b4["iterations_per_step"],
                "total_iterations": b4["total_iterations"],
                "max_node_diff_vs_reference_m": b4["max_node_diff_vs_reference_m"],
                "block_diagonal": b4["block_diagonal"],
                "independent_compare": b4["independent_compare"],
                "max_free_gradient_relative_at_stop": b4["max_free_gradient_relative_at_stop"],
                "gates": b4["gates"], "pass": b4["pass_all_gates"],
            },
        },
        "false_positive_check": {
            "definition": "停止时 free_gradient_inf / gradScale >= 1e-6 即记为提前误停",
            "max_relative_gradient_any_stop": max(
                m3["free_gradient_relative"], m6["free_gradient_relative"],
                b1["max_free_gradient_relative_at_stop"], b4["max_free_gradient_relative_at_stop"]),
            "false_positive": any([
                m3["free_gradient_relative"] >= GRAD_RATIO_TOL,
                m6["free_gradient_relative"] >= GRAD_RATIO_TOL,
                b1["max_free_gradient_relative_at_stop"] >= GRAD_RATIO_TOL,
                b4["max_free_gradient_relative_at_stop"] >= GRAD_RATIO_TOL]),
        },
        "candidate_confirmed_on_s2_s6a": all(c["pass_all_gates"] for c in cases),
        "note": "即使全部通过，也只说明候选判据在 S2 M3/M6 + S6a 1/4 体上重验收通过，不构成协议 v0.3 生效。",
    }
    write_json(os.path.join(out, "judgment.json"), judgment)
    print(json.dumps({k: v["pass"] for k, v in judgment["cases"].items()}, ensure_ascii=False))
    print("ALL_PASS" if judgment["candidate_confirmed_on_s2_s6a"] else "NOT_ALL_PASS")


if __name__ == "__main__":
    main()
