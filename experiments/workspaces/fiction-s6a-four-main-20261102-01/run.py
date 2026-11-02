import csv,hashlib,json,math,os,shutil,subprocess,tarfile
from pathlib import Path
w=Path(__file__).resolve().parent
root=next(p for p in [w,*w.parents] if (p/'fiction.project.json').exists())
layout=json.loads((root/'fiction.project.json').read_text());prod=root/layout['simulator'];out=root/layout['results']/'s6a-four-main-20261102-01';build=Path('/tmp/fiction-s6a-four-main-20261102-build')
assert not (out/'summary.json').exists()
env=dict(os.environ);env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
with (out/'ctest.log').open('w') as f:r=subprocess.run(['ctest','--test-dir',str(build/'production'),'--output-on-failure'],stdout=f,stderr=subprocess.STDOUT,env=env)
if r.returncode:raise SystemExit('Regression failed; stopped')
asset=root/'assets/ipc/tetMeshes/sphere1K.msh'
files=[p for p in (prod/'cpp').rglob('*') if p.is_file() and 'build' not in p.parts]+[p for p in w.rglob('*') if p.is_file()]
with tarfile.open(out/'source.tar.gz','w:gz') as t:
 for p in files:t.add(p,arcname=str(p.relative_to(root)))
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
for name,p in [('asset',asset),('binary',build/'s6a_four_main'),('core',build/'production/libfiction_core.a')]:hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'hashes.json').write_text(json.dumps(hashes,indent=2))
shutil.copy2(build/'CMakeFiles/s6a_four_main.dir/link.txt',out/'link-command.txt');shutil.copy2(build/'CMakeFiles/s6a_four_main.dir/s6a_four_main.cpp.o.d',out/'header-dependencies.txt');shutil.copy2(build/'s6a_four_main',out/'s6a_four_main')
config={'mode':'combined','dt':.005,'steps':10,'tol':1e-8,'max_iter':100,'gradTolRel':1e-6,'gradScale_N':5086.259120998311,'actualStepTol_m':1e-12,'actualStepCount':2,'fail_first':True,'budget_seconds':600,'address_space_bytes':8*1024**3}
(out/'config.json').write_text(json.dumps(config,indent=2))
cmd=['prlimit','--as='+str(8*1024**3),'--','timeout','--kill-after=5','600','/usr/bin/time','-v','-o',str(out/'time.log'),str(build/'s6a_four_main'),str(asset),str(out)]
(out/'command.json').write_text(json.dumps(cmd,indent=2))
with (out/'run.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
(out/'exitcode').write_text(str(r.returncode))
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
s={'exitcode':r.returncode,'pass':False}
if r.returncode==0:
 rows=read(out/'steps.csv');state=read(out/'state.csv');old=read(root/layout['results']/'s6a-multibody-20261031-01/body4_state.csv');assert len(state)==len(old)==7040
 diff=max(math.dist([float(a[k]) for k in ('global_x','global_y','global_z')],[float(b[k]) for k in ('global_x','global_y','global_z')]) for a,b in zip(state,old))
 cmp=read(out/'trajectory_compare.csv');maxcmp=max(float(x['max_position_diff_m']) for x in cmp)
 s.update(steps=len(rows),all_gates=all(x['gate_pass']=='1' for x in rows),old_state_max_node_difference_m=diff,independent_max_difference_m=maxcmp,comparison_rows=len(cmp),min_J=min(float(x['minJ']) for x in rows),max_relative_force_residual=max(float(x['relative_force_residual']) for x in rows),max_DBC_error_m=max(float(x['dbc_error']) for x in rows),iterations=[int(x['iterations']) for x in rows],all_states_finite=all(math.isfinite(float(x[k])) for x in state for k in ('global_x','global_y','global_z')))
 s['pass']=len(rows)==10 and s['all_gates'] and s['all_states_finite'] and diff<1e-8 and maxcmp<1e-8 and len(cmp)==44
(out/'summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
raise SystemExit(0 if s['pass'] else 1)
