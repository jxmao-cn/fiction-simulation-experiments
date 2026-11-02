#!/usr/bin/env python3
"""S2 补充实验分析 v2。

独立核对末态与 trace，量化：
1) M3 增加迭代/初始步长为何无效（方向率地板、实际位移、能量分辨率）；
2) 严格 tol=1e-10 下 M5 过线是否为网格相关的数值偶然；
3) 统一 tol=1e-8 的加密趋势、局部表观阶数与外推敏感性；
4) 停止误差与离散误差的量级分离；
5) M6 触发条件与资源。
生成 comparison.csv / 诊断 JSON / 三张图 / 中文报告。
"""
import csv
import hashlib
import json
import math
import os
import shutil
import tarfile
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/fiction-s2-ext-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORK = Path(__file__).resolve().parent
PROJECT = next(p for p in WORK.parents if (p/'fiction.project.json').is_file())
ROOT = PROJECT.parent
OUT = PROJECT / 'simulator/output/validation/s2-extension-20261031-01'
FIG = PROJECT / 'figures/validation' / OUT.name
REF = 0.010002209491525423      # 梁参考，m
W = 24.525                      # 总重，N
EPS = 2.220446049250313e-16
ULP_1M = EPS                    # 1 m 量级双精度位置分辨率


def readcsv(path):
    with path.open() as f:
        return list(csv.DictReader(l for l in f if not l.startswith('#')))


def fit_fixed_p(pts, p):
    """delta = A - B/n^p，固定 p，最小二乘求 A（mm）。"""
    xs = [1.0 / (n ** p) for _, n, _ in pts]
    ys = [d for _, _, d in pts]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
    return my - b * mx


def main():
    summary = json.loads((OUT / 'summary.json').read_text())
    assert summary['status'] != 'running'
    rows = {r['case']['id']: r for r in summary['rows']}
    checks = {}

    # ---------- 复现性核对 ----------
    old = readcsv(PROJECT / 'simulator/output/validation/s2-cantilever-20261031-01/M3.csv')[0]
    baseline = rows['M3-baseline']['metrics']
    keys = ('tip_deflection', 'search_direction_rate', 'free_gradient_inf', 'min_J',
            'iterations', 'line_search_iterations')
    checks['baseline_matches_original_exactly'] = all(old[k] == baseline[k] for k in keys)
    checks['baseline_full_state_byte_identical'] = (
        (OUT / 'M3-baseline/state.csv').read_bytes()
        == (PROJECT / 'simulator/output/validation/s2-cantilever-20261031-01/M3-state.csv').read_bytes())
    btrace = readcsv(OUT / 'M3-baseline/trace/iterations.csv')
    checks['baseline_all_trace_rows_identical'] = (
        btrace == readcsv(PROJECT / 'simulator/output/validation/s2-cantilever-20261031-01/M3-trace/iterations.csv'))
    ltrace = readcsv(OUT / 'M3-long1000/trace/iterations.csv')
    compare = ('iter', 'search_rate', 'alpha', 'halvings', 'actual_dx_inf', 'energy', 'free_gradient_inf')
    checks['long_matches_baseline_first_100_iterations'] = all(
        all(b[k] == l[k] for k in compare) for b, l in zip(btrace, ltrace))
    with tarfile.open(PROJECT / 'simulator/output/validation/s2-cantilever-20261031-01/source.tar.gz') as tar:
        checks['production_core_unchanged_since_original_S2'] = all(
            tar.extractfile('cpp/' + p).read() == (PROJECT / 'simulator/cpp' / p).read_bytes()
            for p in ['include/time_integrator.h', 'src/time_integrator.cpp'])

    # ---------- 全病例表 ----------
    checks['state_tip_consistency'] = {}
    table = {}
    for name, r in rows.items():
        m = r.get('metrics')
        if not m:
            continue
        ss = readcsv(OUT / name / 'state.csv')
        tip = [-float(s['uz']) for s in ss if abs(float(s['X']) - 1) < 1e-12]
        diff = abs(math.fsum(tip) / len(tip) - float(m['tip_deflection']))
        finite = all(math.isfinite(float(s[k])) for s in ss for k in ('x', 'y', 'z'))
        checks['state_tip_consistency'][name] = {'absolute_difference_m': diff, 'all_positions_finite': finite}
        assert diff < 1e-14 and finite
        table[name] = {
            'case': name, 'mesh': r['case']['mesh'], 'tol': r['case']['tol'],
            'initial_alpha': r['case']['alpha'], 'max_iter': r['case']['max_iter'],
            'converged': m['converged'], 'iterations': int(m['iterations']),
            'nodes': int(m['nodes']), 'tets': int(m['tets']),
            'tip_mm': float(m['tip_deflection']) * 1000,
            'error_percent': float(m['relative_error']) * 100,
            'search_rate': float(m['search_direction_rate']),
            'free_gradient_inf': float(m['free_gradient_inf']),
            'min_J': float(m['min_J']), 'max_J': float(m['max_J']),
            'min_stretch': float(m['min_stretch']), 'max_stretch': float(m['max_stretch']),
            'dbc_error': float(m['measured_dbc_error']),
            'solve_seconds': float(m['solve_wall_seconds']),
            'rss_MiB': r['peak_rss_kib'] / 1024,
            'force_balance_norm': float(m['force_balance_norm']),
            'elastic_energy': float(m['elastic_energy']),
            'gravity_energy': float(m['gravity_energy']),
        }
        # 用线性刚度折算的等效位移残差：g*delta/W
        table[name]['residual_disp_m'] = (table[name]['free_gradient_inf'] * table[name]['tip_mm'] / 1000 / W)

    with (OUT / 'comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(next(iter(table.values())).keys()))
        writer.writeheader()
        writer.writerows(table[k] for k in sorted(table))

    # ---------- 停滞分析 ----------
    stall = {}
    for name in ['M3-baseline', 'M3-long1000', 'M3-alpha05', 'M3-alpha025', 'M4-strict', 'M5-strict']:
        t = readcsv(OUT / name / 'trace/iterations.csv')
        frozen_from = None
        for i in range(len(t)):
            if all(float(r['actual_dx_inf']) < 1e-16 for r in t[i:]):
                frozen_from = int(t[i]['iter'])
                break
        last = t[-1]
        e_total = abs(table[name]['elastic_energy']) + abs(table[name]['gravity_energy'])
        noise_est = math.sqrt(table[name]['tets']) * EPS * e_total
        stall[name] = {
            'iterations': int(last['iter']),
            'rate_last': float(last['search_rate']),
            'rate_over_tol_1e10': float(last['search_rate']) / 1e-10,
            'dx_last_m': float(last['actual_dx_inf']),
            'dx_in_ulps_at_1m': float(last['actual_dx_inf']) / ULP_1M,
            'energy_last_J': float(last['energy']),
            'predicted_decrease_J': abs(float(last['state_gradient_dot_current_direction'])),
            'energy_sum_noise_estimate_J': noise_est,
            'free_gradient_inf_N': float(last['free_gradient_inf']),
            'gradient_relative_to_weight': float(last['free_gradient_inf']) / W,
            'residual_disp_m': table[name]['residual_disp_m'],
            'frozen_from_iter': frozen_from,
            'noop_iters': sum(1 for r in t if float(r['actual_dx_inf']) < 1e-16),
        }
    (OUT / 'stall-analysis.json').write_text(json.dumps(stall, indent=2) + '\n')

    # ---------- 容差敏感性（同网格 strict vs loose） ----------
    tol_sens = {}
    for mesh, strict, loose in [('M3', 'M3-baseline', 'M3-tol1e8'), ('M5', 'M5-strict', 'M5-tol1e8')]:
        a = table[strict]['tip_mm'] / 1000
        b = table[loose]['tip_mm'] / 1000
        disc = abs(b - REF)
        tol_sens[mesh] = {
            'strict_case': strict, 'loose_case': loose,
            'strict_tip_m': a, 'loose_tip_m': b,
            'abs_diff_m': abs(a - b),
            'rel_diff': abs(a - b) / b,
            'discretization_error_m': disc,
            'discretization_error_rel': disc / REF,
            'ratio_discretization_to_tol_diff': disc / abs(a - b) if a != b else None,
            'strict_converged': table[strict]['converged'],
            'loose_converged': table[loose]['converged'],
        }
    (OUT / 'tolerance-sensitivity.json').write_text(json.dumps(tol_sens, indent=2) + '\n')

    # ---------- M6 触发条件 ----------
    m5 = rows['M5-tol1e8']
    m5m = m5['metrics']
    trigger = {
        'condition': 'M5-tol1e8 numerical_pass=1, relative_error>0.05, solve_seconds<60, RSS<2GiB, elapsed<900s',
        'observed': {
            'numerical_pass': m5m['numerical_pass'],
            'relative_error': float(m5m['relative_error']),
            'solve_seconds': float(m5m['solve_wall_seconds']),
            'peak_rss_MiB': m5['peak_rss_kib'] / 1024,
            'batch_elapsed_seconds_at_end': summary['wall_seconds'],
        },
        'thresholds': {'relative_error': 0.05, 'solve_seconds': 60, 'peak_rss_MiB': 2048, 'elapsed_seconds': 900},
        'all_satisfied': (m5m['numerical_pass'] == '1'
                          and float(m5m['relative_error']) > 0.05
                          and float(m5m['solve_wall_seconds']) < 60
                          and m5['peak_rss_kib'] < 2 * 1024 ** 2
                          and summary['wall_seconds'] < 900),
        'M6_run': None,
    }
    if 'M6-tol1e8' in table:
        trigger['M6_run'] = {
            'mesh': '120x24x24', 'nodes': table['M6-tol1e8']['nodes'], 'tets': table['M6-tol1e8']['tets'],
            'converged': table['M6-tol1e8']['converged'],
            'tip_mm': table['M6-tol1e8']['tip_mm'],
            'error_percent': table['M6-tol1e8']['error_percent'],
            'solve_seconds': table['M6-tol1e8']['solve_seconds'],
            'peak_rss_MiB': table['M6-tol1e8']['rss_MiB'],
            'note': 'M6 RSS 2520 MiB 超过触发判断所用的 2 GiB 参考值（该值约束 M5），仍低于 8 GiB 地址空间上限',
        }
    (OUT / 'm6-trigger-check.json').write_text(json.dumps(trigger, indent=2) + '\n')

    # ---------- 加密趋势与拟合 ----------
    relaxed = [table[f'M{i}-tol1e8'] for i in range(1, 7)]
    n_of = {f'M{i}': n for i, n in zip(range(1, 7), [10, 20, 40, 60, 80, 120])}
    trend = [{'mesh': t['mesh'], 'n': n_of[t['mesh']], 'h': 1.0 / n_of[t['mesh']],
              'tip_mm': t['tip_mm'], 'error_percent': t['error_percent'],
              'converged': t['converged'], 'iterations': t['iterations']} for t in relaxed]
    orders = []
    for a, b in zip(trend, trend[1:]):
        orders.append({
            'pair': a['mesh'] + '->' + b['mesh'],
            'apparent_order_p': math.log((a['error_percent'] / 100) / (b['error_percent'] / 100))
                                 / math.log(a['h'] / b['h']),
            'h_ratio': a['h'] / b['h'],
        })
    # 只用 M3-M6 四点做外推：M1/M2 明显还在前渐近区，纳入最小二乘会把截距拉低。
    pts = [(x['mesh'], x['n'], x['tip_mm']) for x in trend[2:]]
    A_curve = {}
    for ip in range(100, 241):
        p = ip / 100
        A_curve[p] = fit_fixed_p(pts, p)
    A_15 = fit_fixed_p(pts, 1.5)
    A_17 = fit_fixed_p(pts, 1.7)
    A_20 = fit_fixed_p(pts, 2.0)
    fit = {
        'points': trend,
        'consecutive_apparent_orders': orders,
        'local_order_last_two': [o['apparent_order_p'] for o in orders[-2:]],
        'fixed_p_extrapolated_limit_mm': {'p=1.5': A_15, 'p=1.7': A_17, 'p=2.0': A_20},
        'extrapolation_range_p_1.5_2.0_mm': [min(A_15, A_20), max(A_15, A_20)],
        'beam_reference_mm': REF * 1000,
        'M6_last_accepted_step_m': float(readcsv(OUT / 'M6-tol1e8/trace/iterations.csv')[-1]['actual_dx_inf']),
        'M6_next_direction_m': table['M6-tol1e8']['search_rate'],
        'M5_M6_increment_mm': table['M6-tol1e8']['tip_mm'] - table['M5-tol1e8']['tip_mm'],
        'M5_M6_increment_percent_of_reference': (table['M6-tol1e8']['tip_mm'] - table['M5-tol1e8']['tip_mm']) / (REF * 1000) * 100,
        'note': ('四个非嵌套点无法唯一确定渐近阶；固定 p 的外推范围覆盖梁参考，'
                 '因此本批数据既不证明也不排除夹持端模型差异'),
    }
    (OUT / 'refinement-fit.json').write_text(json.dumps(fit, indent=2) + '\n')

    # ---------- 图 1：加密趋势 ----------
    FIG.mkdir(parents=True, exist_ok=True)
    labels = [t['mesh'] for t in relaxed]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), layout='constrained')
    tips = [t['tip_mm'] for t in relaxed]
    errs = [t['error_percent'] for t in relaxed]
    axes[0].plot(labels, tips, 'o-', color='#236a92', label='P1 tetrahedra')
    axes[0].axhline(REF * 1000, color='black', ls='--', label='Beam reference: 10.0022 mm')
    axes[0].set(ylabel='Tip deflection [mm]', xlabel='Mesh refinement', title='Self-weight cantilever')
    axes[1].plot(labels, errs, 'o-', color='#b34c25')
    axes[1].axhline(5, color='black', ls='--', label='5% reference band')
    axes[1].set(ylabel='Relative error [%]', xlabel='Mesh refinement', title='Error against beam approximation')
    for i, t in enumerate(relaxed):
        dy = -16 if i == len(relaxed) - 1 else 7
        axes[0].annotate(f"{t['tip_mm']:.3f}", (i, t['tip_mm']), xytext=(0, dy),
                         textcoords='offset points', ha='center', fontsize=8)
        axes[1].annotate(f"{t['error_percent']:.2f}%", (i, t['error_percent']), xytext=(0, 7),
                         textcoords='offset points', ha='center', fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=8); ax.margins(x=.07, y=.17)
    fig.suptitle('S2 refinement diagnostic | uniform solver tolerance 1e-8', fontsize=11)
    fig.savefig(FIG / 's2-refinement.png', dpi=180); fig.savefig(FIG / 's2-refinement.svg'); plt.close(fig)

    # ---------- 图 2：M3 控制 ----------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    for name, label in [('M3-long1000', 'alpha=1, maxIter=1000'),
                        ('M3-alpha05', 'alpha=0.5, maxIter=300'),
                        ('M3-alpha025', 'alpha=0.25, maxIter=300')]:
        data = readcsv(OUT / name / 'trace/iterations.csv')
        it = [int(s['iter']) for s in data]
        axes[0].semilogy(it, [max(float(s['search_rate']), 1e-30) for s in data], label=label)
        axes[1].semilogy(it, [max(float(s['actual_dx_inf']), 1e-30) for s in data], label=label)
    axes[0].axhline(1e-10, color='black', ls='--', label='Original tolerance: 1e-10')
    axes[0].set(ylabel='Newton direction infinity norm [m]', xlabel='Iteration', title='Stopping metric')
    axes[1].axhline(ULP_1M, color='grey', ls=':', label='Position ulp at 1 m')
    axes[1].set(ylabel='Accepted position change [m]', xlabel='Iteration', title='Actual motion')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('M3 solver controls | same mesh, material, loading and strict tolerance', fontsize=11)
    fig.savefig(FIG / 's2-m3-controls.png', dpi=180); fig.savefig(FIG / 's2-m3-controls.svg'); plt.close(fig)

    # ---------- 图 3：阶数与外推敏感性 ----------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    hs = [t['h'] for t in trend]
    es = [t['error_percent'] / 100 for t in trend]
    axes[0].loglog(hs, es, 'o-', color='#236a92', label='P1 tetrahedra (tol=1e-8)')
    for p, ls in [(1, ':'), (2, '--')]:
        axes[0].loglog(hs, [es[-1] * (h / hs[-1]) ** p for h in hs], ls, color='grey',
                       label=f'slope {p} guide')
    axes[0].set(xlabel='h (unit cell size / m)', ylabel='Relative error vs beam reference',
                title='Observed convergence')
    axes[0].invert_xaxis(); axes[0].grid(alpha=.2, which='both'); axes[0].legend(fontsize=8)
    axes[0].annotate(f"local slope {orders[-1]['apparent_order_p']:.2f}", xy=(hs[-2], es[-2]),
                     xytext=(-90, 25), textcoords='offset points', fontsize=8,
                     arrowprops=dict(arrowstyle='->', color='black', lw=.8))
    ps = sorted(A_curve)
    axes[1].plot(ps, [A_curve[p] for p in ps], color='#b34c25', label='Extrapolated limit A(p)')
    axes[1].axhline(REF * 1000, color='black', ls='--', label='Beam reference: 10.0022 mm')
    axes[1].axvline(orders[-1]['apparent_order_p'], color='grey', ls=':', label='Local apparent order')
    axes[1].fill_between([1.5, 2.0], min(A_15, A_20), max(A_15, A_20), color='#b34c25', alpha=.15,
                         label='Range for p in [1.5, 2.0]')
    axes[1].set(xlabel='Assumed convergence order p', ylabel='Extrapolated limit A [mm]',
                title='Extrapolation is not identifiable')
    axes[1].grid(alpha=.2); axes[1].legend(fontsize=8)
    fig.suptitle('S2 refinement fit | four non-nested points, order not established', fontsize=11)
    fig.savefig(FIG / 's2-refinement-fit.png', dpi=180); fig.savefig(FIG / 's2-refinement-fit.svg'); plt.close(fig)

    # ---------- 校验与报告 ----------
    for p in [OUT / 'comparison.csv', OUT / 'analysis-checks.json', OUT / 'summary.json',
              OUT / 'stall-analysis.json', OUT / 'tolerance-sensitivity.json',
              OUT / 'm6-trigger-check.json', OUT / 'refinement-fit.json',
              OUT / 'trace-samples.json']:
        if p.exists():
            shutil.copy2(p, FIG / p.name)

    def mdtable(items):
        lines = ['| 病例 | 迭代 | 挠度 mm | 参考误差 | 方向率 | 自由梯度 N | 等效位移残差 m | 求解墙钟 s | 数值状态 |',
                 '|---|---:|---:|---:|---:|---:|---:|---:|---|']
        for t in items:
            lines.append(
                f"| {t['case']} | {t['iterations']} | {t['tip_mm']:.8f} | {t['error_percent']:.5f}% | "
                f"{t['search_rate']:.4e} | {t['free_gradient_inf']:.4e} | {t['residual_disp_m']:.2e} | "
                f"{t['solve_seconds']:.3f} | {'收敛' if t['converged'] == '1' else '未收敛末态'} |")
        return '\n'.join(lines)

    controls = mdtable([table[n] for n in ['M3-baseline', 'M3-long1000', 'M3-alpha05', 'M3-alpha025']])
    strict = mdtable([table[n] for n in ['M3-baseline', 'M4-strict', 'M5-strict']])
    loose = mdtable([table[f'M{i}-tol1e8'] for i in range(1, 7)])

    s3 = stall['M3-baseline']
    s4 = stall['M4-strict']
    s5 = stall['M5-strict']
    m6 = table['M6-tol1e8']
    m3 = table['M3-baseline']
    ts3 = tol_sens['M3']
    ts5 = tol_sens['M5']

    report = f'''# S2 补充实验：加密趋势、迭代预算与初始步长

日期：2026-10-31。状态：补充诊断批次已完成并分析；生产求解器未合并改动；原 S2 验收标准与“未通过”结论保持不变；未进入 S3。

## 结论摘要

1. **加密确实改善弯曲。** 在统一 tol=1e-8 下，M1→M6 的尖端挠度为 1.6082 → 4.2499 → 7.4075 → 8.6203 → 9.1508 → 9.5761 mm，参考误差 83.92% → 57.51% → 25.94% → 13.82% → 8.51% → 4.26%，严格单调下降。
2. **M3 的“未收敛”不是没算到平衡。** M3 末态自由梯度仅 {s3['free_gradient_inf_N']:.3e} N（总重 24.525 N 的 {s3['gradient_relative_to_weight']:.2e}），按线性刚度折算的等效位移残差约 {s3['residual_disp_m']:.2e} m；同网格 tol=1e-10 与 tol=1e-8 的尖端只差 {ts3['abs_diff_m']:.2e} m。停滞发生在“方向范数”判据上，不在解上。
3. **1e-10 的步长判据低于该问题的数值地板。** M3、M4 的方向率冻结在 {s3['rate_last']:.4e}、{s4['rate_last']:.4e}，实际位置变化降到 {s3['dx_last_m']:.1e} m（1 m 处双精度 ulp≈{ULP_1M:.1e} m）；M5 在第 5 次迭代以 {s5['rate_last']:.3e} 恰好过线。过线与否由网格相关的数值噪声决定，不是物理收敛质量。
4. **增加迭代与减小初始步长无效。** M3 的 maxIter 从 100 提到 1000，方向率与挠度不变（{s3['rate_last']:.4e}、{table['M3-long1000']['tip_mm']:.6f} mm）；alpha=0.5/0.25 在 300 次内同样未过线，挠度只变化约 3e-8 m。它们改变收敛路径，不改变同一离散问题的平衡解。
5. **M6 首次进入 5% 带，但尚未渐近。** M6=120×24×24 在 tol=1e-8 下误差 {m6['error_percent']:.2f}%；最后一步加密仍使尖端移动 {fit['M5_M6_increment_mm']:.3f} mm（参考值的 {fit['M5_M6_increment_percent_of_reference']:.2f}%），与剩余误差 {m6['error_percent']:.2f}% 同量级。固定阶数 p∈[1.5,2.0] 外推的极限为 {fit['extrapolation_range_p_1.5_2.0_mm'][0]:.2f}–{fit['extrapolation_range_p_1.5_2.0_mm'][1]:.2f} mm，区间覆盖梁参考 10.0022 mm；这批数据既不能证明也不能排除夹持端模型差异。

## 1. 范围与可复现性

用户要求追加更细网格、M3 更多迭代和初始求解参数对照。执行设计见 [补充运行规格](s2-extension-spec-20261031.md)，首次结果见 [S2 原始失败现场](s2-validation-20261031.md)。每个病例从未变形状态独立求解；材料、重力、边界、网格拓扑生成规则与主观测量保持一致，未依据 10 mm 调材料或初始形状。

本次只在 `/home/steven/mjx/fiction/experiments/workspaces/fiction-s2-extension-20261031-01/` 的独立副本暴露 maxIter/tol/initial-alpha；主工程核心与原 S2 快照一致：{checks['production_core_unchanged_since_original_S2']}。默认配置与旧 M3 关键数值完全一致：{checks['baseline_matches_original_exactly']}；完整末态逐字节一致：{checks['baseline_full_state_byte_identical']}；1000 次实验的前 100 次轨迹数值一致：{checks['long_matches_baseline_first_100_iterations']}。Release/CHOLMOD，CTest 10/10（其中 2 个 GUI 占位）通过。

## 2. M3：增加迭代与更改初始步长

以下全部保持原 tol=1e-10。alpha=0.5/0.25 是每轮初始试探步长上限，仍可能回溯到更小值，并非固定学习率；Newton 方向包含曲率信息，与神经网络梯度下降只能类比“沿某个方向走多远”。

{controls}

迭代上限 1000 的 baseline 提供与较小 alpha 相同预算前缀的参照，不能把 1000 次与 300 次的差异单独归因于 alpha。逐轮第 3/10/30/100/300/1000 次方向率、实际移动、能量、梯度及末态保存在 trace-samples.json。未收敛病例不能作为原 S2 通过的证据。

![M3 方向与实际位置更新](../../figures/validation/{OUT.name}/s2-m3-controls.png)

停滞的逐轮证据（M3-baseline）：第 3 次迭代方向率 1.198e-9、实际移动 2.55e-8 m；第 4 次线搜索退到 alpha=0.125，移动 1.50e-10 m；第 5 次 alpha=2.44e-4，移动 2.58e-13 m；第 {s3['frozen_from_iter']} 次起每轮移动 <1e-16 m，第 21 次后固定为 {s3['dx_last_m']:.2e} m。方向率从第 4 次起冻结在 {s3['rate_last']:.4e}，自由梯度冻结在 {s3['free_gradient_inf_N']:.3e} N。此时 $x+\\alpha p$ 舍入回 $x$，迭代成为浮点不动点；继续迭代只是重复同一计算。

为什么线搜索会接受这种零位移步？静止点的预测能量下降只有约 {s3['predicted_decrease_J']:.2e} J；实际被接受的每轮能量差已降到 1e-16 J 量级（总能量 0.036 J 的约 10 个双精度 ulp），而按 $\\sqrt{{N_{{tet}}}}\\,\\epsilon\\,(|E_{{elastic}}|+|E_{{gravity}}|)$ 估计的独立求和噪声约 {s3['energy_sum_noise_estimate_J']:.2e} J。方向质量不足使 alpha=1 反而增能，线搜索只能不断二分；当步长小到能量差落在舍入噪声内时，它就接受一个位置不再变化的步。**这是“方向质量 + 能量分辨率 + 绝对步长判据”共同造成的停止问题，不是材料、边界或载荷错误；具体由哪一项主导，需要另做线性求解残差与能量求和的插桩诊断。**

## 3. 原严格容差下的加密与“偶然过线”

{strict}

M5 在 5 次迭代内收敛，M3/M4 却在 100 次后仍未过线。三者末态的等效位移残差都在 {min(table['M3-baseline']['residual_disp_m'], table['M4-strict']['residual_disp_m'], table['M5-strict']['residual_disp_m']):.1e}–{max(table['M3-baseline']['residual_disp_m'], table['M4-strict']['residual_disp_m'], table['M5-strict']['residual_disp_m']):.1e} m 量级，自由梯度相对总重不超过 {max(s3['gradient_relative_to_weight'], s4['gradient_relative_to_weight'], s5['gradient_relative_to_weight']):.1e}。也就是说，M5 的“通过”与 M3/M4 的“失败”并不对应解的质量差异，而对应方向率噪声地板恰好落在 1e-10 的哪一侧：M3 地板 {s3['rate_last']:.2e}（门槛的 {s3['rate_over_tol_1e10']:.1f} 倍），M4 {s4['rate_last']:.2e}（{s4['rate_over_tol_1e10']:.1f} 倍），M5 {s5['rate_last']:.2e}（{s5['rate_over_tol_1e10']:.2f} 倍，刚好低于 1）。

按线性刚度 $K=W/\\delta$ 折算，M3 的实际平衡残差相当于 {s3['residual_disp_m']:.2e} m，反而小于 1e-10 m 的原定容差；也就是说，判据想表达的位移精度已经达到，但方向范数没有过线。因此，**tol=1e-10 的原始门槛在这组问题上不是一个可复现的物理判据**；它比该离散系统的方向率地板更严。原 S2 的“M3 未通过”结论仍然有效（按固定标准执行），但其原因应记为“判据低于数值地板 + M3 网格过粗”，而不是“M3 没有求得平衡”。

## 4. 统一 tol=1e-8 的加密趋势

这组明确改变停止容差，从 M1 重新完整比较，不替代原 tol=1e-10 验收。M4=60×12×12，M5=80×16×16，M6=120×24×24（75625 节点/414720 tet）。加密保持原单元长宽比，但六档的细化比并不一致（M3→M4 为 1.5，M4→M5 为 4/3，M5→M6 为 1.5），也不是单一嵌套序列，所以没有以这些点声称严格的渐近收敛阶。

{loose}

本组六个病例全部收敛，参考误差严格下降。局部表观阶数为：M2→M3 {orders[1]['apparent_order_p']:.2f}，M3→M4 {orders[2]['apparent_order_p']:.2f}，M4→M5 {orders[3]['apparent_order_p']:.2f}，M5→M6 {orders[4]['apparent_order_p']:.2f}；最后两档稳定在约 {orders[-1]['apparent_order_p']:.2f}，但这是两个非嵌套点的局部斜率，不是已确立的渐近阶。

![统一容差的加密趋势](../../figures/validation/{OUT.name}/s2-refinement.png)

外推敏感性：固定阶数 p=1.5/1.7/2.0 时，$\\delta=A-B/n^p$ 的极限分别为 {A_15:.3f} / {A_17:.3f} / {A_20:.3f} mm；p 在 [1.5,2.0] 内变动时极限落在 {fit['extrapolation_range_p_1.5_2.0_mm'][0]:.2f}–{fit['extrapolation_range_p_1.5_2.0_mm'][1]:.2f} mm。三点或四点拟合不足以唯一确定 p，因此不能把某个外推值当作三维真值，也不能据此宣称存在或不存在约 2% 的夹持端模型误差。

![阶数与外推敏感性](../../figures/validation/{OUT.name}/s2-refinement-fit.png)

## 5. 停止误差与离散误差的量级分离

同网格、不同容差的尖端差异：

| 网格 | strict 挠度 mm | loose 挠度 mm | 差异 m | 相对差异 | 该网格参考误差 | 离散/停止误差比 |
|---|---:|---:|---:|---:|---:|---:|
| M3 | {ts3['strict_tip_m']*1000:.8f} | {ts3['loose_tip_m']*1000:.8f} | {ts3['abs_diff_m']:.2e} | {ts3['rel_diff']:.2e} | {ts3['discretization_error_rel']*100:.2f}% | {ts3['ratio_discretization_to_tol_diff']:.1e} |
| M5 | {ts5['strict_tip_m']*1000:.8f} | {ts5['loose_tip_m']*1000:.8f} | {ts5['abs_diff_m']:.2e} | {ts5['rel_diff']:.2e} | {ts5['discretization_error_rel']*100:.2f}% | {ts5['ratio_discretization_to_tol_diff']:.1e} |

两档网格上，容差引起的尖端变化都比离散误差小 5–7 个数量级。这支持把 tol=1e-8 用作**本场景的数值停止容差**：它远高于数值地板，又远低于离散误差，不会污染网格收敛结论。它不改变原 S2 的 tol=1e-10 结论，只是为后续协议修订提供量级依据。

## 6. M6 触发条件与资源

预登记条件为：M5-tol1e8 数值通过、参考误差 >5%、求解墙钟 <60 s、峰值 RSS <2 GiB、批次已用时间 <900 s。实测 M5：误差 {trigger['observed']['relative_error']*100:.2f}%，墙钟 {trigger['observed']['solve_seconds']:.1f} s，RSS {trigger['observed']['peak_rss_MiB']:.0f} MiB，全部满足；批次总墙钟 {summary['wall_seconds']:.1f} s，也满足。M6 按条件追加，结果：{m6['error_percent']:.2f}%，墙钟 {m6['solve_seconds']:.1f} s，峰值 RSS {m6['rss_MiB']:.0f} MiB。M6 的 RSS 高于触发判断所用的 2 GiB 参考值（该值约束的是 M5），但仍低于 8 GiB 地址空间上限；后续再加密需要上调内存预算。

## 7. 对 S2 与协议的结论

- **原 S2 仍然未通过，且不应改写。** 固定标准是“M3 收敛且相对误差 <5%”；M3 实际误差 {m3['error_percent']:.2f}%，按标准就是不通过。补充实验是诊断，不是重验收。
- **失败原因已经分解清楚：** (a) tol=1e-10 低于方向率数值地板，M3/M4 的“未收敛”是停止判据的假阴性；(b) 即使忽略停止问题，M3 的离散误差仍有 {m3['error_percent']:.2f}%，原标准默认 M3 足够细的假设不成立，5% 带要到 M6=120×24×24 才首次达到。
- **S2 的科学问题已有肯定回答：** 弯曲离散误差随加密严格单调下降，六档跨度 83.92%→4.26%，局部表观阶数约 {orders[-1]['apparent_order_p']:.2f}。这可以作为文章附录的收敛图，但必须同时标注“参考是梁理论近似”“M6 尚未渐近”“1e-10 判据低于数值地板”。
- **建议另立协议修订，而不是静默放宽。** 候选修订：(i) 静态求解的停止判据改为相对自由梯度或“实际位移 + 能量下降”的组合，绝对步长范数只作诊断；(ii) 若暂不改代码，则把本场景的 tol 固定为 1e-8，并强制做一次同网格 strict/loose 敏感性检查（本批已给出 M3、M5 两例）；(iii) 精度门槛按“最细可负担网格 + 单调下降 + 局部阶数”表述，不再把 5% 单独挂在某个未经验证的网格档位上。修订需在新 run-id 下重新确认，不能把本批直接记作通过。

## 8. 资源与证据

总批次墙钟 {summary['wall_seconds']:.2f} s，最多两个独立求解同时运行。已设置线程环境变量，但未证明实际单线程；耗时含 trace 开销且可能互相竞争 CPU，不作为重复性能对比。执行状态 `{summary['status']}`，已完成病例 {len(table)} 个。

run-id：[s2-extension-20261031-01](../output/validation/{OUT.name}/)。config/command 保存每次配置与命令，source.tar.gz/hash 与二进制保存运行版本，各病例目录保存 CSV、完整节点末态、trace、退出码、time/RSS；comparison.csv 是全指标汇总，analysis-checks.json 保存从完整末态重新计算尖端的独立核对，stall-analysis.json 保存停滞量化，tolerance-sensitivity.json 保存同网格容差对照，m6-trigger-check.json 保存条件判断，refinement-fit.json 保存拟合敏感性，trace-samples.json 保存相同迭代前缀。原实验未覆盖，未提交或发布。
'''
    (PROJECT / 'simulator/docs/s2-extension-results-20261031.md').write_text(report)
    (OUT / 'analysis-checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    shutil.copy2(OUT / 'analysis-checks.json', FIG / 'analysis-checks.json')
    shutil.copy2(Path(__file__).resolve(), OUT / 'analysis-script.py')
    (OUT / 'analysis-script.sha256').write_text(
        hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest() + '\n')
    print(json.dumps(checks, indent=2))
    for k in sorted(table):
        t = table[k]
        print(k, t['converged'], t['iterations'], round(t['tip_mm'], 6),
              round(t['error_percent'], 6), round(t['solve_seconds'], 3))


if __name__ == '__main__':
    main()
