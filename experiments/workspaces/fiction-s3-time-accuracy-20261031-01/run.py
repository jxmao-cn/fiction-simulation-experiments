#!/usr/bin/env python3
import csv,json,math,os,subprocess,sys,time,hashlib,tarfile,shutil
from pathlib import Path
W=Path(__file__).resolve().parent
R=next(p for p in W.parents if (p/'fiction.project.json').exists())
P=json.loads((R/'fiction.project.json').read_text())
O=R/P['results']/'s3-time-accuracy-20261031-01'
B=Path('/tmp/fiction-s3-time-accuracy-20261031-build')
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')

def run(cmd,log,timeout=600):
    with log.open('w') as f:return subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=ENV,timeout=timeout).returncode
def read(p):
    with p.open() as f:return list(csv.DictReader(l for l in f if not l.startswith('#')))
def sample(rows,dt):
    out=[0.]
    for j in range(1,101):
        i=round(j*.005/dt)-1
        assert abs(float(rows[i]['time'])-j*.005)<1e-12
        out.append(float(rows[i]['tip_deflection']))
    return out
def error(a,b):
    d=[x-y for x,y in zip(a,b)]
    return dict(max_m=max(map(abs,d)),rms_m=math.sqrt(math.fsum(v*v for v in d)/len(d)))

def main():
    O.mkdir(parents=True,exist_ok=True)
    if (O/'config.json').exists():raise SystemExit('refuse overwrite')
    cfg={'spec':'s3-time-accuracy-revision-20261031.md','mode':'original step','cases':[
        ['ref',.0003125,1600,1e-6],['refcheck',.00015625,3200,1e-6],
        ['ref-tight',.0003125,1600,1e-7],['refcheck-tight',.00015625,3200,1e-7],
        ['dt005',.005,100,1e-6],['dt0025',.0025,200,1e-6],['dt00125',.00125,400,1e-6],['dt000625',.000625,800,1e-6]],
        'case_seconds':180,'batch_seconds':900,'rms_gate_m':.001,'max_gate_m':.002,'reference_gate_m':.00025,'algebraic_gate_m':.00001}
    (O/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    shutil.copy2(R/P['reports']/cfg['spec'],O/'spec.md')
    # Only the header protocol label changes; numerical code remains original step implementation.
    runner=W/'cpp/examples/dynamic_bar_check.cpp'
    runner.write_text(runner.read_text().replace('protocol=fiction-validation-v0.3-S3','protocol=s3-time-accuracy-candidate-20261031'))
    with tarfile.open(O/'source.tar.gz','w:gz') as t:
        t.add(W/'cpp',arcname='cpp');t.add(Path(__file__),arcname='run.py')
    mesh=R/P['results']/'s3-bend-20261031-01/M2.msh'
    (O/'mesh.sha256').write_text(hashlib.sha256(mesh.read_bytes()).hexdigest()+'\n')
    for cmd,name in [(['cmake','-S',str(W/'cpp'),'-B',str(B),'-DCMAKE_BUILD_TYPE=Release'],'configure'),(['cmake','--build',str(B),'-j','3'],'build'),(['ctest','--test-dir',str(B),'--output-on-failure'],'ctest')]:
        if run(cmd,O/(name+'.log')):raise SystemExit(name+' failed')
    shutil.copy2(B/'dynamic_bar_check',O/'dynamic_bar_check')
    (O/'binary.sha256').write_text(hashlib.sha256((O/'dynamic_bar_check').read_bytes()).hexdigest()+'\n')
    start=time.monotonic();results=[];data={};status='running'
    for name,dt,n,tol in cfg['cases']:
        remaining=min(180,900-(time.monotonic()-start))
        if remaining<=0:status='batch_timeout';break
        d=O/name;d.mkdir()
        cmd=['prlimit','--as='+str(8*1024**3),'--','timeout','--kill-after=5',str(remaining),'/usr/bin/time','-v','-o',str(d/'time.log'),str(O/'dynamic_bar_check'),'--mesh',str(mesh),'--out-dir',str(d),'--dt',str(dt),'--steps',str(n),'--tol',str(tol),'--max-iter','50','--mode','bend','--gravity-z','-9.81']
        (d/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        rc=run(cmd,d/'run.log',int(remaining)+15)
        (d/'exitcode').write_text(str(rc)+'\n')
        row={'case':name,'dt':dt,'tol':tol,'exitcode':rc}
        rows=read(d/'steps.csv') if (d/'steps.csv').exists() else []
        row['steps']=len(rows)
        good=rc==0 and len(rows)==n and all(r['converged']=='1' and float(r['min_J'])>.3 and float(r['dbc_error'])<1e-6 for r in rows)
        row['numerical_pass']=good
        if good:
            for r in rows:
                for k in ('tip_deflection','min_J','max_J','search_direction_rate','free_gradient_inf'):
                    good=good and math.isfinite(float(r[k]))
        row['numerical_pass']=good
        results.append(row)
        if not good:status='numerical_or_resource_failure';break
        data[name]=sample(rows,dt)
        if name=='refcheck':
            row['reference_check']=error(data['ref'],data['refcheck'])
            if row['reference_check']['max_m']>.00025:status='reference_not_established';break
        if name.endswith('-tight'):
            row['algebraic_check']=error(data[name[:-6]],data[name])
            if row['algebraic_check']['max_m']>.00001:status='algebraic_accuracy_not_established';break
        if name.startswith('dt'):
            row['time_error']=error(data[name],data['ref'])
            row['time_accuracy_pass']=row['time_error']['max_m']<=.002 and row['time_error']['rms_m']<=.001
        print(name,row,flush=True)
        (O/'summary.json').write_text(json.dumps({'status':'running','results':results},indent=2)+'\n')
    else:status='candidate_passed' if all(r.get('time_accuracy_pass',True) for r in results) else 'candidate_accuracy_failed'
    s={'status':status,'results':results,'wall_seconds':time.monotonic()-start,'original_S3_status':'unchanged: failed','production_modified':False}
    (O/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
    lines=['# S3 短窗口时间精度候选：执行结果','',f"日期：2026-10-31。状态：**{status}**。原S3未通过状态不变。",'',
           '执行设计见 [候选规格](s3-time-accuracy-revision-20261031.md)。采用原step判据；只更改CSV协议标签，未采用新combined算法。生产源码与旧run未改。',
           '',f"批次数值墙钟 {s['wall_seconds']:.2f} s，实际完成 {len(results)} 个病例。每病例180s、总900s上限。",'',
           '| 病例 | dt | tol | 完成步数 | 退出码 | 数值通过 |','|---|---:|---:|---:|---:|---|']
    for r in results:lines.append(f"| {r['case']} | {r['dt']} | {r['tol']} | {r['steps']} | {r['exitcode']} | {r['numerical_pass']} |")
    lines+=['','## 资格与停止','', '```json',json.dumps(s,indent=2),'```','',
            '未执行病例不会记为通过。若数值失败，先停止并保留现场；不放宽tol、不追加迭代追绿。参考或代数资格失败时，不进行粗档精度判定。',
            '',f"证据：[s3-time-accuracy-20261031-01](../output/validation/{O.name}/)，含spec/config、源码快照、二进制及输入哈希、构建/CTest、实际命令、逐步CSV、末态、退出码及time记录。",
            '', '本结果仅适用于候选0–0.5s窗口；人工阻尼、相位和振幅差都属于时间误差，未做扣除或相位对齐。']
    (R/P['reports']/'s3-time-accuracy-results-20261031.md').write_text('\n'.join(lines)+'\n')
    print('FINAL',status,flush=True)

if __name__=='__main__':main()
