import csv,json,math,os,subprocess,time,tarfile,hashlib,shutil
from pathlib import Path
w=Path(__file__).resolve().parent
root=next(p for p in w.parents if (p/'fiction.project.json').exists());layout=json.loads((root/'fiction.project.json').read_text())
prod=root/layout['simulator'];out=root/layout['results']/'s3-algebra-main-20261102-01';build=Path('/tmp/fiction-s3-algebra-main-20261102-build')
out.mkdir(parents=True,exist_ok=True);assert not (out/'config.json').exists(),'no overwrite'
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
cases=[('ref',.0003125,1600,1e-6),('ref-tight',.0003125,1600,1e-7),('refcheck',.00015625,3200,1e-6),('refcheck-tight',.00015625,3200,1e-7)]
(out/'config.json').write_text(json.dumps({'cases':cases,'mode':'step','max_iter':50,'case_seconds':300,'batch_seconds':1200,'algebraic_gate_m':1e-5},indent=2)+'\n')
shutil.copy2(prod/'docs/s3-algebra-main-spec-20261102.md',out/'spec.md')
def run(cmd,log):
 with log.open('w') as f:return subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
def read(p):
 with p.open() as f:return list(csv.DictReader(l for l in f if not l.startswith('#')))
def sample(rows,dt):
 vals=[0.]
 for j in range(1,101):
  row=rows[round(j*.005/dt)-1];assert abs(float(row['time'])-j*.005)<1e-12
  vals.append(float(row['tip_deflection']))
 return vals
def error(a,b):
 ds=[x-y for x,y in zip(a,b)];return {'max_m':max(map(abs,ds)),'rms_m':math.sqrt(math.fsum(x*x for x in ds)/len(ds))}
for cmd,name in [(['cmake','-S',str(w),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release'],'configure'),(['cmake','--build',str(build),'-j','3'],'build'),(['ctest','--test-dir',str(build/'production'),'--output-on-failure'],'ctest')]:
 if run(cmd,out/(name+'.log')):raise SystemExit(name+' failed')
paths=[p for p in (prod/'cpp').rglob('*') if p.is_file() and 'build' not in p.parts]+[p for p in w.rglob('*') if p.is_file()]
with tarfile.open(out/'source.tar.gz','w:gz') as t:
 for p in paths:t.add(p,arcname=str(p.relative_to(root)))
mesh=root/layout['results']/'s3-bend-20261031-01/M2.msh'
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
for label,p in [('binary',build/'dynamic_bar_check'),('core',build/'production/libfiction_core.a'),('mesh',mesh)]:hashes[label]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
shutil.copy2(build/'CMakeFiles/dynamic_bar_check.dir/link.txt',out/'link-command.txt');shutil.copy2(build/'CMakeFiles/dynamic_bar_check.dir/dynamic_bar_check.cpp.o.d',out/'header-dependencies.txt');shutil.copy2(build/'dynamic_bar_check',out/'dynamic_bar_check')
start=time.monotonic();results=[];samples={};status='passed'
for name,dt,n,tol in cases:
 d=out/name;d.mkdir();remaining=min(300,1200-(time.monotonic()-start))
 if remaining<=0:status='batch_timeout';break
 cmd=['prlimit','--as='+str(8*1024**3),'--','timeout','--kill-after=5',str(remaining),'/usr/bin/time','-v','-o',str(d/'time.log'),str(out/'dynamic_bar_check'),'--mesh',str(mesh),'--out-dir',str(d),'--dt',str(dt),'--steps',str(n),'--tol',str(tol),'--max-iter','50','--mode','bend','--gravity-z','-9.81']
 (d/'command.json').write_text(json.dumps(cmd,indent=2)+'\n');rc=run(cmd,d/'run.log');(d/'exitcode').write_text(str(rc)+'\n')
 rows=read(d/'steps.csv') if (d/'steps.csv').exists() else []
 good=rc==0 and len(rows)==n and all(r['converged']=='1' and float(r['min_J'])>.3 and float(r['dbc_error'])<1e-6 for r in rows)
 good=good and all(math.isfinite(float(r[k])) for r in rows for k in ['tip_deflection','min_J','max_J','search_direction_rate','free_gradient_inf'])
 result={'case':name,'dt':dt,'tol':tol,'steps':len(rows),'exitcode':rc,'numerical_pass':good};results.append(result)
 if not good:status='numerical_or_resource_failure';break
 samples[name]=sample(rows,dt)
 if name.endswith('-tight'):
  base=name[:-6];result['tip_comparison']=error(samples[name],samples[base]);a=read(d/'state.csv');b=read(out/base/'state.csv');assert len(a)==len(b)
  result['max_node_state_difference_m']=max(math.dist([float(x[k]) for k in ['x','y','z']],[float(y[k]) for k in ['x','y','z']]) for x,y in zip(a,b))
  with (d/'algebraic_comparison.csv').open('w') as f:
   writer=csv.writer(f);writer.writerow(['time','baseline_tip_m','tight_tip_m','difference_m'])
   writer.writerows((j*.005,samples[base][j],samples[name][j],samples[name][j]-samples[base][j]) for j in range(101))
  result['algebraic_pass']=result['tip_comparison']['max_m']<=1e-5 and result['max_node_state_difference_m']<=1e-5
  if not result['algebraic_pass']:status='algebraic_accuracy_not_established';break
 print(name,result,flush=True)
 (out/'summary.json').write_text(json.dumps({'status':'running','results':results},indent=2)+'\n')
summary={'status':status,'results':results,'wall_seconds':time.monotonic()-start,'original_S3_status':'unchanged: failed','production_modified':False}
if 'refcheck' in samples:summary['time_difference_diagnostic']=error(samples['ref'],samples['refcheck'])
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# S3代数求解精度资格：主工程执行结果','','日期2026-11-02。状态：**'+status+'**。原S3未通过；旧参考解资格未建立。','', '采用真实生产core的step模式，材料/加载/几何固定；运行前规格见 [规格](s3-algebra-main-spec-20261102.md)。','', '```json',json.dumps(summary,indent=2),'```','', '松紧容差门槛：101共同物理时刻尖端最大差与完整末态节点最大差均<=10 μm。数值/资源/门槛失败即停，未运行项标为未运行。两时间步差记录为诊断指标。','',f'证据：[运行目录](../output/validation/{out.name}/)，保存源码快照、哈希、生产core链接/头依赖、构建/CTest、命令/CSV/末态/资源/逐采样差值。未改生产代码。']
(prod/'docs/s3-algebra-main-results-20261102.md').write_text('\n'.join(lines)+'\n');print('FINAL',status,flush=True)
