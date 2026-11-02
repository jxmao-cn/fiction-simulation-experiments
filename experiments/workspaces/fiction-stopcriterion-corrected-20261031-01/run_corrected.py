from pathlib import Path
import os,json,subprocess,csv,hashlib,tarfile,sys,shutil
w=Path(__file__).resolve().parent
root=next(p for p in w.parents if (p/'fiction.project.json').exists())
layout=json.loads((root/'fiction.project.json').read_text());out=root/layout['results']/'stopcriterion-corrected-20261031-01'
build=Path('/tmp/fiction-stopcriterion-corrected-build')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
config={'purpose':'isolated corrected criterion diagnostics; not formal acceptance','dynamic_residual':'raw IP gradient / h^2 [N]','normalization':'single body weight shared by 1/4 body and independent solves','grad_tol_rel':1e-6,'actual_step_tol_m':1e-12,'steps':10,'max_iter':100,'python':sys.executable}
(out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
with (out/'ctest.log').open('w') as log:
 r=subprocess.run(['ctest','--test-dir',str(build),'--output-on-failure'],env=env,stdout=log,stderr=subprocess.STDOUT)
if r.returncode:raise SystemExit('ctest failed')
with (out/'units-negative-check.log').open('w') as log:
 subprocess.run([str(build/'stopcriterion_units_check')],stdout=log,stderr=subprocess.STDOUT,check=True)
with tarfile.open(out/'source.tar.gz','w:gz') as tar:
 tar.add(w/'cpp',arcname='cpp');tar.add(Path(__file__).resolve(),arcname='run_corrected.py')
base=root/layout['workspaces']/'fiction-s2s6a-combined-20261031-01'
with (out/'source-change.diff').open('w') as f:subprocess.run(['diff','-ruN',str(base/'cpp'),str(w/'cpp')],stdout=f)
hashes={str(p.relative_to(w)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (w/'cpp').rglob('*') if p.is_file()}
for b in ['cantilever_check','s6a_multibody_check','stopcriterion_units_check']:hashes['binary/'+b]=hashlib.sha256((build/b).read_bytes()).hexdigest()
(out/'hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
args=['--stop-mode','combined','--grad-tol','1e-6','--actual-step-tol','1e-12','--actual-step-count','2']
rows=[]
for name in ['M3','body1','body4']:
 d=out/name;d.mkdir(exist_ok=False)
 if name=='M3':
  cmd=[str(build/'cantilever_check'),'--mesh',str(root/layout['results']/'s2-cantilever-20261031-01/meshes/M3.msh'),'--csv',str(d/'metrics.csv'),'--state',str(d/'state.csv'),'--trace-dir',str(d/'trace'),*args,'--grad-scale','24.525']
 else:
  cmd=[str(build/'s6a_multibody_check'),'--mesh',str(root/layout['assets']/'ipc/tetMeshes/sphere1K.msh'),'--out',str(d),'--run-id',out.name,*args]
  if name=='body1':cmd+=['--only-body1']
 (d/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
 with (d/'run.log').open('w') as log:
  r=subprocess.run(['prlimit','--as=8589934592','--','timeout','--kill-after=5','600','/usr/bin/time','-v','-o',str(d/'time.log'),*cmd],env=env,stdout=log,stderr=subprocess.STDOUT)
 (d/'exitcode').write_text(str(r.returncode)+'\n');rows.append({'case':name,'exitcode':r.returncode})
 print(rows[-1],flush=True)
 if r.returncode not in (0,1,2):break
(out/'summary.json').write_text(json.dumps({'runs':rows,'production_modified':False,'original_acceptance_changed':False},indent=2)+'\n')
