import csv,hashlib,json,math,os,shutil,subprocess,tarfile
from pathlib import Path
w=Path(__file__).resolve().parent
root=next(p for p in [w,*w.parents] if (p/'fiction.project.json').exists())
layout=json.loads((root/'fiction.project.json').read_text());prod=root/layout['simulator'];results=root/layout['results'];out=results/'s6a-single-main-20261102-01';build=Path('/tmp/fiction-s6a-single-main-20261102-build')
assert not (out/'summary.json').exists(),'refuse overwrite'
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']})
with (out/'rebuild.log').open('w') as f:subprocess.run(['cmake','--build',str(build),'-j','3'],stdout=f,stderr=subprocess.STDOUT,check=True)
with (out/'ctest.log').open('w') as f:subprocess.run(['ctest','--test-dir',str(build/'production'),'--output-on-failure'],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
asset=root/'assets/ipc/tetMeshes/sphere1K.msh'
paths=[p for p in (prod/'cpp').rglob('*') if p.is_file() and 'build' not in p.parts]+[p for p in w.rglob('*') if p.is_file()]
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
hashes['asset']=hashlib.sha256(asset.read_bytes()).hexdigest()
(out/'hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
with tarfile.open(out/'source.tar.gz','w:gz') as t:
 for p in paths:t.add(p,arcname=str(p.relative_to(root)))
shutil.copy2(build/'CMakeFiles/s6a_single_main.dir/s6a_single_main.cpp.o.d',out/'header-dependencies.txt')
shutil.copy2(build/'CMakeFiles/s6a_single_main.dir/link.txt',out/'link-command.txt')
shutil.copy2(build/'s6a_single_main',out/'s6a_single_main')
config={'run_id':out.name,'production_core':str(prod/'cpp/src/time_integrator.cpp'),'asset':str(asset),'material':{'E':1e5,'nu':.3,'rho':1000},'gravity_z':-9.81,'dt':.005,'steps':10,'tol':1e-8,'max_iter':100,'mode':'combined','gradTolRel':1e-6,'gradScale_N':5086.259120998311,'actualStepTol_m':1e-12,'actualStepCount':2,'max_seconds':600,'address_space_bytes':8*1024**3,'fail_first':True}
(out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
cmd=['prlimit','--as='+str(8*1024**3),'--','timeout','--kill-after=5','600','/usr/bin/time','-v','-o',str(out/'time.log'),str(build/'s6a_single_main'),str(asset),str(out)]
(out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (out/'run.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
(out/'exitcode').write_text(str(r.returncode)+'\n')
def rows(p):
 with p.open() as f:return list(csv.DictReader(l for l in f if not l.startswith('#')))
summary={'exitcode':r.returncode,'pass':False}
if (out/'steps.csv').exists():
 data=rows(out/'steps.csv');state=rows(out/'state.csv');diffs={}
 for name,path in [('original',results/'s6a-multibody-20261031-01/body1_state.csv'),('corrected',results/'stopcriterion-corrected-20261031-01/body1/body1_state.csv')]:
  ref=rows(path);assert len(ref)==len(state)
  diffs[name]=max(math.dist([float(s[k]) for k in ['global_x','global_y','global_z']],[float(t[k]) for k in ['global_x','global_y','global_z']]) for s,t in zip(state,ref))
 allfinite=all(math.isfinite(float(s[k])) for s in state for k in ['global_x','global_y','global_z'])
 traces=[]
 for i in range(1,len(data)+1):
  tr=rows(out/f'iterations.csv.step-{i}.csv');traces.append({'step':i,'final_actual_dx_inf_m':tr[-1]['actual_dx_inf'],'last_two':[{k:x[k] for k in ['iter','actual_dx_inf','alpha','free_gradient_inf']} for x in tr[-2:]]})
 summary.update(steps=len(data),total_iterations=sum(int(d['iterations']) for d in data),all_step_gates=all(d['gate_pass']=='1' for d in data),all_state_finite=allfinite,state_differences_m=diffs,trace_stop_motion=traces,min_J=min(float(d['minJ']) for d in data),max_force_residual_N=max(float(d['force_residual_inf_N']) for d in data),max_relative_force_residual=max(float(d['relative_force_residual']) for d in data),max_DBC_error_m=max(float(d['dbc_error']) for d in data),stop_reasons=[d['stop_reason'] for d in data],iterations=[int(d['iterations']) for d in data])
 summary['pass']=r.returncode==0 and len(data)==10 and summary['all_step_gates'] and allfinite and all(v<1e-8 for v in diffs.values())
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
