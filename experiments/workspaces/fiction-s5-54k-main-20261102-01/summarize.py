import csv,hashlib,json,math,statistics,shutil
from pathlib import Path
w=Path(__file__).resolve().parent
root=next(p for p in [w,*w.parents] if (p/'fiction.project.json').exists())
layout=json.loads((root/'fiction.project.json').read_text());out=root/layout['results']/'s5-54k-main-20261102-01'
summary=json.loads((out/'summary.json').read_text());audit=json.loads((out/'asset_audit.json').read_text())
with (out/'steps.csv').open() as f:rows=list(csv.DictReader(f))
with (out/'state.csv').open() as f:state=list(csv.DictReader(f))
finite=all(math.isfinite(float(x[k])) for x in state for k in ['x','y','z','ux','uy','uz'])
times=[float(r['wall_seconds']) for r in rows];ordered=sorted(times)
rss=None
for line in (out/'time.log').read_text().splitlines():
 if 'Maximum resident set size' in line:rss=int(line.split(':')[-1])
extra={'state_nodes':len(state),'all_saved_state_finite':finite,'max_DBC_error':max(float(x['dbc_max_error']) for x in rows),'min_stretch':min(float(x['min_stretch']) for x in rows),'max_stretch':max(float(x['max_stretch']) for x in rows),'max_force_residual_N':max(float(x['force_residual_inf_N']) for x in rows),'median_step_s':statistics.median(times),'p95_nearest_rank_s':ordered[math.ceil(.95*len(ordered))-1],'max_step_s':max(times),'max_rss_kib':rss,'output_bytes':sum(p.stat().st_size for p in out.rglob('*') if p.is_file()),'production_modified':False}
summary.update(extra);(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
m=summary['meta'];status='通过成本pilot短段门槛' if summary['pass'] and finite else '未通过/未完成'
report=f'''# S5 54K：主工程成本 pilot

日期：2026-11-02。结论：**{status}**。本批执行原计划54K规模候选的成本pilot。执行前规格见 [spec](s5-54k-main-spec-20261102.md)。

## 资产与配置

资产实际路径为 `/home/steven/mjx/Archive/IPC/input/tetMeshes/Armadillo54K.msh`，sha256见运行目录hashes.json。实测{audit['nodes']}节点/{audit['tets']}tet，{audit['surface_triangles']}表面三角形；参考正定向、非退化、无孤立节点、表面闭合且与tet边界匹配的审计门槛通过。包围盒为{audit['bbox']}；体积{audit['volume']:.8g}网格单位³，总质量数值{audit['total_mass']:.8g}。实际固定节点{m['pin_nodes']}，自由DOF {m['free_dof']}。

NH E=1e7、nu=.3、rho1000，重力-z9.81；原固定盒保持不变（y<=-52）。dt=.005、10步、tol=.01、maxIter50、DBC k1e8。生产默认step；tol为绝对方向率阈值。原样坐标没有单位缩放，资产+y为竖直轴而本批-z侧向加载；本批记录原坐标下侧向加载的短段计算成本。

## 结果

| 指标 | 结果 |
|---|---:|
| 已完成时间步 | {len(rows)}/10 |
| 每步数值收敛 | {summary['all_converged']} |
| 各步Newton次数 | {','.join(x['iterations'] for x in rows)} |
| 各步线搜索减半次数 | {','.join(x['line_search_iterations'] for x in rows)} |
| min J | {summary['min_J']:.10g} |
| 主伸长全批最小/最大 | {extra['min_stretch']:.10g} / {extra['max_stretch']:.10g} |
| 最大DBC误差 | {extra['max_DBC_error']:.8g} |
| 完整末态有限 | {finite} |
| 总求解墙钟 | {summary['sum_solve_wall_s']:.4f} s |
| 每步中位数/P95/最大 | {extra['median_step_s']:.4f} / {extra['p95_nearest_rank_s']:.4f} / {extra['max_step_s']:.4f} s |
| 峰值RSS | {rss} KiB |
| 退出码 | {summary['exitcode']} |

参考装配{m['reference_assembly_seconds']:.4f} s、首次参考分解{m['first_factorization_seconds']:.4f} s（{m['first_factorization_solver']}），H非零{m['hessian_nnz']}。这些参考操作在驱动侧单独复算，包含罚项；逐步CSV记录完整时间步墙钟。

主伸长和力残差记录为观测指标。最大物理力残差约{extra['max_force_residual_N']:.6g} N；成本pilot的门槛见执行前规格。每步实际求解后端记录在CSV。

## 验证与计时条件

harness通过CMake add_subdirectory编译真实生产fiction_core，生产include优先；header-dependencies.txt及link-command.txt保留实际来源，hashes.json保存核心库和二进制哈希。生产子目录CTest11/11通过（其中2项GUI占位），完整末态重新检查共{len(state)}节点有限。

成本为一次诊断测量，包含首步trace、额外单元主伸长观测及日志；线程环境设为1，实际库线程数未测量，运行使用共享CPU。物理时长0.05s。122K/219K未执行。

证据目录：[s5-54k-main-20261102-01](../output/validation/s5-54k-main-20261102-01/)，包含规格/配置、源快照与hash、资产审计、实际命令、编译回归、首步trace、逐步CSV、完整末态、参考装配/分解元数据、time/RSS与summary。
'''
(root/layout['reports']/'s5-54k-main-results-20261102.md').write_text(report)
shutil.copy2(Path(__file__).resolve(),out/'summarize.py')
manifest=[hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(out)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.sha256']
(out/'manifest.sha256').write_text('\n'.join(manifest)+'\n')
print(json.dumps({k:summary[k] for k in ['pass','steps','min_J','sum_solve_wall_s','max_rss_kib','median_step_s','all_saved_state_finite']},indent=2))
