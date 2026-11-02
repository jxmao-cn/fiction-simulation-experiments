#!/usr/bin/env python3
"""停止判据诊断汇总：逐 case 写 config.json，并对照基线做预登记判定。

输入（由 run_stopcriterion.sh 产生）：
  OUT/s2/M3/{step,grad,actual,combined}/M3.csv, M3-state.csv
  OUT/s6a/body1/{step,grad,actual,combined}/body1_steps.csv, body1_state.csv
输出：
  OUT/s2/M3/<mode>/config.json, summary.json
  OUT/s6a/body1/<mode>/config.json, summary.json
  OUT/cases.csv, OUT/judgment.json
"""
import argparse
import csv
import json
import os

MODES = ["step", "grad", "actual", "combined"]
POS_TOL = 1e-8      # 物理解不变门槛 [m]
GRAD_RATIO_TOL = 1e-6  # 不提前误停门槛


def read_csv_rows(path):
    rows = []
    with open(path) as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    return list(csv.DictReader(lines))


def fnum(row, key):
    return float(row[key])


def s2_case(out, mode, grad_scale):
    d = os.path.join(out, "s2", "M3", mode)
    csv_path = os.path.join(d, "M3.csv")
    state_path = os.path.join(d, "M3-state.csv")
    row = read_csv_rows(csv_path)[0]
    state = []
    with open(state_path) as f:
        for r in csv.DictReader(f):
            state.append((float(r["x"]), float(r["y"]), float(r["z"])))
    grad_scale = float(row["grad_scale"])
    res = {
        "case": "s2/M3",
        "mode": mode,
        "tol": 1e-10,
        "max_iter": 100,
        "grad_scale_n": grad_scale,
        "grad_tol_rel": fnum(row, "grad_tol_rel"),
        "actual_step_tol": fnum(row, "actual_step_tol"),
        "actual_step_count": int(float(row["actual_step_count"])),
        "iterations": int(float(row["iterations"])),
        "line_search_iterations": int(float(row["line_search_iterations"])),
        "converged": int(float(row["converged"])),
        "max_iter_hit": int(float(row["max_iter_hit"])),
        "line_search_failed": int(float(row["line_search_failed"])),
        "stop_reason": row["stop_reason"],
        "search_direction_rate": fnum(row, "search_direction_rate"),
        "free_gradient_inf": fnum(row, "free_gradient_inf"),
        "free_gradient_relative": fnum(row, "free_gradient_inf") / grad_scale,
        "tip_deflection_m": fnum(row, "tip_area_weighted"),
        "tip_min_m": fnum(row, "tip_min"),
        "tip_max_m": fnum(row, "tip_max"),
        "min_J": fnum(row, "min_J"),
        "max_J": fnum(row, "max_J"),
        "dbc_satisfied": int(float(row["dbc_satisfied"])),
        "dbc_total": int(float(row["dbc_total"])),
        "dbc_max_error": fnum(row, "dbc_max_error"),
        "solver": row["solver"],
        "wall_seconds_solver": fnum(row, "solve_wall_seconds"),
        "numerical_pass_flag": int(float(row["numerical_pass"])),
        "state": state,
    }
    return res


def s6a_case(out, mode):
    d = os.path.join(out, "s6a", "body1", mode)
    steps = read_csv_rows(os.path.join(d, "body1_steps.csv"))
    state = []
    with open(os.path.join(d, "body1_state.csv")) as f:
        for r in csv.DictReader(f):
            state.append((float(r["global_x"]), float(r["global_y"]), float(r["global_z"])))
    rows = []
    for r in steps:
        rows.append({
            "step": int(r["step"]),
            "converged": int(r["converged"]),
            "iterations": int(r["iterations"]),
            "line_search_iterations": int(r["line_search_iterations"]),
            "stop_reason": r["stop_reason"],
            "max_iter_hit": int(r["max_iter_hit"]),
            "search_direction_rate": float(r["search_direction_rate"]),
            "free_gradient_inf": float(r["free_gradient_inf"]),
            "state_change_inf": float(r["state_change_inf"]),
            "min_J": float(r["minJ"]),
            "max_disp": float(r["max_disp"]),
            "dbc_satisfied": int(r["dbc_satisfied"]),
            "dbc_total": int(r["dbc_total"]),
            "wall_seconds": float(r["wall_seconds"]),
        })
    # gradScale 从 config.json 读（含场景总重）
    with open(os.path.join(d, "config.json")) as f:
        cfg = json.load(f)
    grad_scale = cfg["grad_scale_n"]
    for r in rows:
        r["free_gradient_relative"] = r["free_gradient_inf"] / grad_scale
    res = {
        "case": "s6a/body1",
        "mode": mode,
        "tol": 1e-8,
        "max_iter": 100,
        "grad_scale_n": grad_scale,
        "steps": rows,
        "final_max_disp_m": rows[-1]["max_disp"],
        "total_iterations": sum(r["iterations"] for r in rows),
        "total_line_search_iterations": sum(r["line_search_iterations"] for r in rows),
        "total_wall_seconds": sum(r["wall_seconds"] for r in rows),
        "all_steps_converged": all(r["converged"] == 1 for r in rows),
        "max_iter_steps": [r["step"] for r in rows if r["max_iter_hit"] == 1],
        "min_J": min(r["min_J"] for r in rows),
        "state": state,
    }
    return res


def max_pos_diff(a, b):
    return max(max(abs(x[i] - y[i]) for i in range(3)) for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--s2-grad-scale", required=True)
    ap.add_argument("--s6a-grad-scale", required=True)
    args = ap.parse_args()
    out = args.out

    # ---------------- 逐 case config.json ----------------
    for mode in MODES:
        d = os.path.join(out, "s2", "M3", mode)
        cfg = {
            "run_id": "stopcriterion-20261031-01",
            "case": "S2 M3 static cantilever (NeoHookean tet, E=5.90e8, nu=0.3, rho=1000)",
            "mesh": "fiction/simulator/output/validation/s2-cantilever-20261031-01/meshes/M3.msh",
            "stop_mode": mode,
            "tol": 1e-10,
            "max_iter": 100,
            "static_solve": True,
            "grad_tol_rel": 1e-7,
            "grad_scale_n": float(args.s2_grad_scale),
            "actual_step_tol_m": 1e-12,
            "actual_step_count": 2,
            "reference_state": "fiction/simulator/output/validation/s2-cantilever-20261031-01/M3-state.csv",
        }
        with open(os.path.join(d, "config.json"), "w") as f:
            json.dump(cfg, f, indent=2)
        d = os.path.join(out, "s6a", "body1", mode)
        cfg = {
            "run_id": "stopcriterion-20261031-01",
            "case": "S6a 1-body (sphere1K, NH E=1e5, nu=0.3, rho=1000, top cap DBC)",
            "mesh": "Archive/IPC/input/tetMeshes/sphere1K.msh",
            "stop_mode": mode,
            "tol": 1e-8,
            "max_iter": 100,
            "dt": 0.005,
            "steps": 10,
            "gravity": [0.0, 0.0, -9.81],
            "dbc_k": 1e8,
            "grad_tol_rel": 1e-7,
            "grad_scale_n": float(args.s6a_grad_scale),
            "actual_step_tol_m": 1e-12,
            "actual_step_count": 2,
            "only_body1": True,
            "reference_state": "fiction/simulator/output/validation/s6a-multibody-20261031-01/body1_state.csv",
        }
        with open(os.path.join(d, "config.json"), "w") as f:
            json.dump(cfg, f, indent=2)

    s2 = {m: s2_case(out, m, args.s2_grad_scale) for m in MODES}
    s6a = {m: s6a_case(out, m) for m in MODES}

    # ---------------- 判定 ----------------
    judgment = {"run_id": "stopcriterion-20261031-01", "pos_tol_m": POS_TOL,
                "grad_ratio_tol": GRAD_RATIO_TOL, "cases": {}}

    base = s2["step"]
    s2j = {}
    for m in MODES:
        c = s2[m]
        v = {
            "stop_reason": c["stop_reason"],
            "iterations": c["iterations"],
            "line_search_iterations": c["line_search_iterations"],
            "converged": c["converged"],
            "max_iter_hit": c["max_iter_hit"],
            "free_gradient_relative_at_stop": c["free_gradient_relative"],
            "search_direction_rate_at_stop": c["search_direction_rate"],
            "tip_deflection_m": c["tip_deflection_m"],
            "max_node_pos_diff_vs_baseline_m": max_pos_diff(c["state"], base["state"]),
            "tip_diff_vs_baseline_m": abs(c["tip_deflection_m"] - base["tip_deflection_m"]),
            "wall_seconds_solver": c["wall_seconds_solver"],
            "dbc_ok": c["dbc_satisfied"] == c["dbc_total"],
            "min_J": c["min_J"],
        }
        v["physical_solution_unchanged"] = bool(v["max_node_pos_diff_vs_baseline_m"] < POS_TOL
                                                and v["tip_diff_vs_baseline_m"] < POS_TOL)
        v["no_premature_stop"] = bool(v["free_gradient_relative_at_stop"] < GRAD_RATIO_TOL)
        # 基线 maxIter 行程 -> 变体是否消除
        if base["max_iter_hit"] == 1:
            v["false_negative_removed"] = bool(c["converged"] == 1 and c["stop_reason"] != "maxIter")
        else:
            v["false_negative_removed"] = None
        s2j[m] = v
    judgment["cases"]["s2/M3"] = {
        "baseline": {"stop_reason": base["stop_reason"], "iterations": base["iterations"],
                     "converged": base["converged"], "max_iter_hit": base["max_iter_hit"],
                     "free_gradient_relative_at_stop": base["free_gradient_relative"],
                     "tip_deflection_m": base["tip_deflection_m"]},
        "modes": s2j,
    }

    b6 = s6a["step"]
    s6j = {}
    fn_steps = b6["max_iter_steps"]
    for m in MODES:
        c = s6a[m]
        rows = c["steps"]
        variant_hits = [r["step"] for r in rows if r["max_iter_hit"] == 1]
        fn_removed = True
        fn_detail = {}
        for s in fn_steps:
            r = rows[s - 1]
            ok = (r["converged"] == 1 and r["stop_reason"] != "maxIter")
            fn_detail[str(s)] = {"converged": r["converged"], "stop_reason": r["stop_reason"],
                                 "iterations": r["iterations"],
                                 "free_gradient_relative": r["free_gradient_relative"],
                                 "false_negative_removed": bool(ok)}
            fn_removed = fn_removed and ok
        v = {
            "stop_reasons": [r["stop_reason"] for r in rows],
            "iterations": [r["iterations"] for r in rows],
            "line_search_iterations": [r["line_search_iterations"] for r in rows],
            "converged_per_step": [r["converged"] for r in rows],
            "total_iterations": c["total_iterations"],
            "total_line_search_iterations": c["total_line_search_iterations"],
            "total_wall_seconds": c["total_wall_seconds"],
            "all_steps_converged": c["all_steps_converged"],
            "max_iter_steps": variant_hits,
            "max_node_pos_diff_vs_baseline_m": max_pos_diff(c["state"], b6["state"]),
            "final_max_disp_m": c["final_max_disp_m"],
            "max_free_gradient_relative_at_stop": max(r["free_gradient_relative"] for r in rows),
            "min_J": c["min_J"],
            "false_negative_removed": bool(fn_removed),
            "false_negative_detail": fn_detail,
        }
        v["physical_solution_unchanged"] = bool(v["max_node_pos_diff_vs_baseline_m"] < POS_TOL)
        v["no_premature_stop"] = bool(v["max_free_gradient_relative_at_stop"] < GRAD_RATIO_TOL)
        s6j[m] = v
    judgment["cases"]["s6a/body1"] = {
        "baseline": {"max_iter_steps": fn_steps, "total_iterations": b6["total_iterations"],
                     "final_max_disp_m": b6["final_max_disp_m"],
                     "stop_reasons": [r["stop_reason"] for r in b6["steps"]]},
        "modes": s6j,
    }

    for case in ("s2/M3", "s6a/body1"):
        for m, v in judgment["cases"][case]["modes"].items():
            v["pass_all_gates"] = bool(v["physical_solution_unchanged"] and v["no_premature_stop"]
                                        and v.get("false_negative_removed", True))

    with open(os.path.join(out, "judgment.json"), "w") as f:
        json.dump(judgment, f, indent=2)

    # ---------------- 表格 ----------------
    with open(os.path.join(out, "cases.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "mode", "stop_reason", "iterations", "line_search_iterations",
                    "converged", "max_iter_hit", "search_direction_rate",
                    "free_gradient_inf_n", "free_gradient_relative",
                    "max_node_pos_diff_vs_baseline_m", "tip_or_max_disp_m",
                    "wall_seconds", "physical_unchanged", "no_premature_stop",
                    "false_negative_removed", "pass_all_gates"])
        for m in MODES:
            v = judgment["cases"]["s2/M3"]["modes"][m]
            c = s2[m]
            w.writerow(["s2/M3", m, v["stop_reason"], v["iterations"], v["line_search_iterations"],
                        v["converged"], v["max_iter_hit"], v["search_direction_rate_at_stop"],
                        c["free_gradient_inf"], v["free_gradient_relative_at_stop"],
                        v["max_node_pos_diff_vs_baseline_m"], v["tip_deflection_m"],
                        v["wall_seconds_solver"], v["physical_solution_unchanged"],
                        v["no_premature_stop"], v["false_negative_removed"], v["pass_all_gates"]])
        for m in MODES:
            v = judgment["cases"]["s6a/body1"]["modes"][m]
            c = s6a[m]
            w.writerow(["s6a/body1", m, "|".join(v["stop_reasons"]), v["total_iterations"],
                        v["total_line_search_iterations"], int(v["all_steps_converged"]),
                        len(v["max_iter_steps"]), "|".join("%.3e" % x for x in
                        [r["search_direction_rate"] for r in c["steps"]]),
                        max(r["free_gradient_inf"] for r in c["steps"]),
                        v["max_free_gradient_relative_at_stop"],
                        v["max_node_pos_diff_vs_baseline_m"], v["final_max_disp_m"],
                        v["total_wall_seconds"], v["physical_solution_unchanged"],
                        v["no_premature_stop"], v["false_negative_removed"], v["pass_all_gates"]])

    for m in MODES:
        with open(os.path.join(out, "s2", "M3", m, "summary.json"), "w") as f:
            json.dump(s2[m], f, indent=2)
        with open(os.path.join(out, "s6a", "body1", m, "summary.json"), "w") as f:
            json.dump(s6a[m], f, indent=2)
    print("wrote judgment.json / cases.csv / summaries under", out)


if __name__ == "__main__":
    main()
