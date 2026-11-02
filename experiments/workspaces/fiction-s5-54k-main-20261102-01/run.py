import csv,hashlib,json,os,shutil,subprocess,tarfile
from pathlib import Path
w=Path(__file__).resolve().parent
root=next(p for p in [w,*w.parents] if (p/'fiction.project.json').exists())
layout=json.loads((root/'fiction.project.json').read_text());prod=root/layout['simulator']
out=root/layout['results']/'s5-54k-main-20261102-01';build=Path('/tmp/fiction-s5-54k-main-20261102-build')
assert not (out/'summary.json').exists(),'Refuse overwrite'
asset=root.parent/'Archive/IPC/input/tetMeshes/Armadillo54K.msh'
env=dict(os.environ);env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']})
with (out/'ctest.log').open('w') as f:
 subprocess.run(['ctest','--test-dir',str(build/'production'),'--output-on-failure'],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
spec=prod/'docs/s5-54k-main-spec-20261102.md';shutil.copy2(spec,out/'spec.md')
paths=[p for p in (prod/'cpp').rglob('*') if p.is_file() and 'build' not in p.parts]+[p for p in w.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
for name,p in [('asset',asset),('binary',build/'s5_54k_main'),('production_core',build/'production/libfiction_core.a')]:hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
with tarfile.open(out/'source.tar.gz','w:gz') as t:
 for p in paths:t.add(p,arcname=str(p.relative_to(root)))
for a,b in [('CMakeFiles/s5_54k_main.dir/link.txt','link-command.txt'),('CMakeFiles/s5_54k_main.dir/s5_scale_check.cpp.o.d','header-dependencies.txt'),('s5_54k_main','s5_54k_main')]:shutil.copy2(build/a,out/b)
args=[str(build/'s5_54k_main'),'--mesh',str(asset),'--steps','10','--pin-box','-63.52761','-54.19321','-57.689','63.53793','-52','57.70009','--csv',str(out/'steps.csv'),'--state',str(out/'state.csv'),'--audit',str(out/'asset_audit.json'),'--meta',str(out/'meta.json'),'--trace',str(out/'first-step.csv'),'--label',out.name]
cmd=['prlimit','--as='+str(32*1024**3),'--','timeout','--kill-after=5','600','/usr/bin/time','-v','-o',str(out/'time.log'),*args]
(out/'config.json').write_text(json.dumps({'run_id':out.name,'mode':'step','parameters':'spec.md','physical_seconds':.05,'thread_env':{k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']},'address_space_bytes':32*1024**3,'timeout_seconds':600},indent=2)+'\n')
(out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (out/'run.log').open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
(out/'exitcode').write_text(str(r.returncode)+'\n')
summary={'exitcode':r.returncode,'pass':False}
if (out/'steps.csv').exists():
 with (out/'steps.csv').open() as f:rows=list(csv.DictReader(f))
 summary.update(steps=len(rows),all_converged=all(x['converged']=='1' for x in rows),min_J=min(float(x['minJ']) for x in rows) if rows else None,sum_solve_wall_s=sum(float(x['wall_seconds']) for x in rows))
 summary['pass']=r.returncode==0 and len(rows)==10 and summary['all_converged']
if (out/'meta.json').exists():summary['meta']=json.loads((out/'meta.json').read_text())
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
