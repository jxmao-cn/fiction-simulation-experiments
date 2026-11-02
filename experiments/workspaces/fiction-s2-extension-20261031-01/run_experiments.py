#!/usr/bin/env python3
"""S2预登记独立补实验，两个有界求解并发；生产源码只读。"""
import concurrent.futures as futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import threading
import time

WORK = Path(__file__).resolve().parent
PROJECT = next(p for p in WORK.parents if (p/'fiction.project.json').is_file())
PROD = PROJECT/'simulator'
OUT = PROD/'output/validation/s2-extension-20261031-01'
BUILD = Path('/tmp/fiction-s2-extension-20261031-build')
SPECS = [('M1',10,2),('M2',20,4),('M3',40,8),('M4',60,12),('M5',80,16)]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT/'config.json').exists():
        raise SystemExit('Refusing overwrite of existing batch')
    env = dict(os.environ)
    env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']})
    cases = [dict(id='M3-long1000',mesh='M3',tol=1e-10,alpha=1.,max_iter=1000),
             dict(id='M3-baseline',mesh='M3',tol=1e-10,alpha=1.,max_iter=100),
             dict(id='M3-alpha05',mesh='M3',tol=1e-10,alpha=.5,max_iter=300),
             dict(id='M3-alpha025',mesh='M3',tol=1e-10,alpha=.25,max_iter=300),
             *[dict(id=m+'-strict',mesh=m,tol=1e-10,alpha=1.,max_iter=100) for m in ['M4','M5']],
             *[dict(id=m+'-tol1e8',mesh=m,tol=1e-8,alpha=1.,max_iter=100) for m,_,_ in SPECS]]
    config = {'run_id':OUT.name,'python':sys.executable,'cases':cases,'max_concurrent':2,
              'case_seconds':600,'batch_seconds':1500,'address_space_bytes':8*1024**3,
              'thread_env':{k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']},
              'conditional_M6':'M5-tol1e8 numerical_pass=1, error>0.05, solve_seconds<60, RSS<2GiB',
              'purpose':'diagnostic; does not replace original S2 acceptance',
              'spec':'docs/s2-extension-spec-20261031.md'}
    (OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    shutil.copy2(PROD/'docs/s2-extension-spec-20261031.md',OUT/'execution-spec.md')
    paths=[p for sub in ('cpp','tools') for p in (WORK/sub).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    paths += [Path(__file__).resolve()]
    with tarfile.open(OUT/'source.tar.gz','w:gz') as t:
        for p in paths:t.add(p,arcname=str(p.relative_to(WORK)))
    hashes={str(p.relative_to(WORK)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    shutil.copy2(BUILD/'cantilever_check',OUT/'cantilever_check')
    hashes['cantilever_check_binary']=hashlib.sha256((OUT/'cantilever_check').read_bytes()).hexdigest()
    (OUT/'source-hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    with (OUT/'ctest.log').open('w') as log:
        r=subprocess.run(['ctest','--test-dir',str(BUILD),'--output-on-failure'],env=env,stdout=log,stderr=subprocess.STDOUT)
    if r.returncode:raise SystemExit('Regression failed; no physical runs')
    for cmd,name in [(['lscpu'],'cpu.txt'),(['free','-b'],'memory-before.txt'),(['uname','-a'],'uname.txt')]:
        with (OUT/name).open('w') as f:subprocess.run(cmd,stdout=f,check=True)
    start=time.monotonic()
    unsafe=threading.Event()
    rows=[]
    lock=threading.Lock()

    def generate(m,nx,ny):
        mesh=OUT/'meshes'/f'{m}.msh'
        cmd=[sys.executable,str(WORK/'tools/make_structured_bar.py'),'--nx',str(nx),'--ny',str(ny),'--nz',str(ny),'--output',str(mesh),'--metadata',str(mesh.with_suffix('.json'))]
        with (OUT/f'{m}-generate.log').open('w') as f:
            r=subprocess.run(['timeout','--kill-after=5','120',*cmd],stdout=f,stderr=subprocess.STDOUT)
        (OUT/f'{m}-generate.exitcode').write_text(str(r.returncode)+'\n')
        if r.returncode:unsafe.set();return False
        return True

    for spec in SPECS:
        if not generate(*spec):raise SystemExit('Mesh generation failed; stopped')

    def run(c):
        if unsafe.is_set():return dict(case=c,status='not_started_after_unsafe_result')
        remain=min(600,1500-(time.monotonic()-start))
        if remain<=0:return dict(case=c,status='not_started_budget')
        d=OUT/c['id'];d.mkdir()
        cmd=['prlimit',f'--as={8*1024**3}','--','timeout','--signal=TERM','--kill-after=5',str(remain),
             '/usr/bin/time','-v','-o',str(d/'time.log'),str(OUT/'cantilever_check'),
             '--mesh',str(OUT/'meshes'/f"{c['mesh']}.msh"),'--csv',str(d/'metrics.csv'),
             '--state',str(d/'state.csv'),'--trace-dir',str(d/'trace'),
             '--tol',str(c['tol']),'--max-iter',str(c['max_iter']),'--initial-alpha',str(c['alpha'])]
        (d/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        with (d/'run.log').open('w') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        (d/'exitcode').write_text(str(r.returncode)+'\n')
        metrics=None
        if (d/'metrics.csv').exists():
            with (d/'metrics.csv').open() as f:raw=list(csv.DictReader(l for l in f if not l.startswith('#')))
            if raw:metrics=raw[-1]
        rss=None
        if (d/'time.log').exists():
            for line in (d/'time.log').read_text().splitlines():
                if 'Maximum resident set size' in line:rss=int(line.split(':')[-1])
        status='converged' if metrics and metrics['numerical_pass']=='1' else 'solver_nonconvergence'
        if r.returncode not in (0,1) or not metrics:
            status='resource_or_execution_failure';unsafe.set()
        elif any(not math.isfinite(float(metrics[k])) for k in ('min_J','max_J','tip_deflection','free_gradient_inf','min_stretch','max_stretch')) or float(metrics['min_J'])<=.9 or float(metrics['measured_dbc_error'])>=1e-6:
            status='unsafe_state';unsafe.set()
        result=dict(case=c,status=status,exitcode=r.returncode,peak_rss_kib=rss,metrics=metrics)
        (d/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        with lock:
            print(c['id'],status,'tip_mm='+str(float(metrics['tip_deflection'])*1000) if metrics else '',flush=True)
        return result

    with futures.ThreadPoolExecutor(max_workers=2) as pool:
        jobs={pool.submit(run,c):c for c in cases}
        for job in futures.as_completed(jobs):
            rows.append(job.result())
            (OUT/'summary.json').write_text(json.dumps({'status':'running','rows':rows},indent=2)+'\n')
    m5=next((r for r in rows if r['case']['id']=='M5-tol1e8'),None)
    if not unsafe.is_set() and m5 and m5.get('metrics'):
        m=m5['metrics']
        if m['numerical_pass']=='1' and float(m['relative_error'])>.05 and float(m['solve_wall_seconds'])<60 and m5['peak_rss_kib']<2*1024**2 and time.monotonic()-start<900:
            if generate('M6',120,24):rows.append(run(dict(id='M6-tol1e8',mesh='M6',tol=1e-8,alpha=1.,max_iter=100)))
    summary={'status':'diagnostic_complete' if not unsafe.is_set() else 'stopped_unsafe_or_resource',
             'rows':rows,'wall_seconds':time.monotonic()-start,'original_S2_passed':False,'S3_started':False}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary['status'], 'batch_seconds=',summary['wall_seconds'],flush=True)


if __name__=='__main__':main()
