import csv, hashlib, json, math, os, shutil, subprocess, tarfile
from pathlib import Path
w=Path(__file__).resolve().parent
root=next(p for p in [w,*w.parents] if (p/'fiction.project.json').exists())
layout=json.loads((root/'fiction.project.json').read_text());prod=root/layout['simulator'];results=root/layout['results']
build=Path('/tmp/fiction-s3-visual-main-20261102-build')
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']})
font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'

def readcsv(p):
    with p.open() as f:return list(csv.DictReader(l for l in f if not l.startswith('#')))

def overlay(title, time_text):
    text=f'{title} | Physical time {time_text}\nNH E=44 MPa, nu=0.3 | No contact | Capability sample; accuracy unverified'
    # Text is constant except ffmpeg's generated timestamp expression.
    return f"drawtext=fontfile={font}:text='{text}':fontcolor=white:fontsize=19:x=18:y=16:box=1:boxcolor=black@0.7:boxborderw=8"

def main():
    all_results=[]
    for mode,run_id,old_id,steps,fps,stride in [
        ('bend','s3-visual-main-20261102-01','s3-bend-20261031-01',2000,25,8),
        ('twist','s3t-visual-main-20261102-01','s3t-twist-20261031-01',500,20,10)]:
        out=results/run_id;out.mkdir(parents=True,exist_ok=False)
        fig=root/layout['figures']/run_id;fig.mkdir(parents=True,exist_ok=False)
        asset=results/old_id/'M2.msh'
        paths=[p for p in (prod/'cpp').rglob('*') if p.is_file() and 'build' not in p.parts]+[p for p in w.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
        hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        hashes.update(asset=hashlib.sha256(asset.read_bytes()).hexdigest(),binary=hashlib.sha256((build/'dynamic_bar_visual').read_bytes()).hexdigest(),production_core=hashlib.sha256((build/'production/libfiction_core.a').read_bytes()).hexdigest())
        (out/'hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
        with tarfile.open(out/'source.tar.gz','w:gz') as tar:
            for p in paths:tar.add(p,arcname=str(p.relative_to(root)))
        shutil.copy2(build/'dynamic_bar_visual',out/'dynamic_bar_visual')
        shutil.copy2('/tmp/fiction-s3-visual-configure.log',out/'configure.log');shutil.copy2('/tmp/fiction-s3-visual-build.log',out/'build.log')
        shutil.copy2(build/'CMakeFiles/dynamic_bar_visual.dir/link.txt',out/'link-command.txt')
        shutil.copy2(build/'CMakeFiles/dynamic_bar_visual.dir/dynamic_bar_visual.cpp.o.d',out/'header-dependencies.txt')
        shutil.copy2(prod/'docs/s3-s3t-visual-spec-20261102.md',out/'spec.md')
        config={'mode':mode,'run_id':run_id,'dt':.005,'steps':steps,'tol':1e-4,'max_iter':50,'stop_mode':'step','E':4.4e7,'nu':.3,'rho':1000,'gravity_z':-9.81 if mode=='bend' else 0,'twist_deg':90,'twist_time':2,'video_physical_end_s':2 if mode=='bend' else 2.5,'video_fps':fps,'sample_every_steps':stride,'max_state_diff_m':1e-8,'camera':{'target':[.5,0,-.03],'radius':1.5,'yaw':-75,'pitch':18,'fov':45},'size':[960,540],'supersampling':2,'environment':{k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']}}
        (out/'config-visual.json').write_text(json.dumps(config,indent=2)+'\n')
        cmd=['prlimit','--as='+str(8*1024**3),'--','timeout','--kill-after=5','600','/usr/bin/time','-v','-o',str(out/'time.log'),str(out/'dynamic_bar_visual'),'--mesh',str(asset),'--out-dir',str(out),'--mode',mode,'--dt','.005','--steps',str(steps),'--tol','1e-4','--max-iter','50','--twist-deg','90','--twist-time','2','--gravity-z',str(config['gravity_z'])]
        (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        with (out/'run.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
        (out/'exitcode').write_text(str(r.returncode)+'\n')
        record={'run_id':run_id,'exitcode':r.returncode,'numerical_pass':False,'media_complete':False}
        if (out/'steps.csv').exists():record['completed_steps']=len(readcsv(out/'steps.csv'))
        if r.returncode==0:
            cur=readcsv(out/'state.csv');old=readcsv(results/old_id/'state.csv');assert len(cur)==len(old)
            diff=max(math.dist([float(x[k]) for k in ['x','y','z']],[float(y[k]) for k in ['x','y','z']]) for x,y in zip(cur,old))
            record.update(max_state_diff_m=diff,numerical_pass=diff<1e-8)
            rows=readcsv(out/'steps.csv');record.update(all_steps_converged=all(row['converged']=='1' for row in rows),min_J=min(float(row['min_J']) for row in rows))
            if record['numerical_pass']:
                media=[]
                title='S3 bend' if mode=='bend' else 'S3t twist 90 deg'
                end=2. if mode=='bend' else 2.5
                for p in sorted(out.glob('keyframe_t*.png')):
                    t=int(p.stem.split('t')[-1])/1000
                    cmdm=['ffmpeg','-y','-v','error','-i',str(p),'-vf',overlay(title,f'{t:.3f} s'),'-frames:v','1',str(fig/p.name)]
                    media.append(cmdm);z=subprocess.run(cmdm,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    if z.returncode:record['media_error']=z.stderr.decode();break
                    shutil.copy2(p.with_suffix('.obj'),fig/p.with_suffix('.obj').name)
                if 'media_error' not in record:
                    cmdm=['ffmpeg','-y','-v','error','-framerate',str(fps),'-i',str(out/'raw_frames/frame_%04d.png'),'-vf',overlay(title,'%{pts\\:hms} / '+f'{end:.3f} s'),'-c:v','libx264','-pix_fmt','yuv420p','-crf','18',str(fig/'video.mp4')]
                    media.append(cmdm);z=subprocess.run(cmdm,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    if z.returncode:record['media_error']=z.stderr.decode()
                    else:
                        probe=subprocess.run(['ffprobe','-v','error','-count_frames','-select_streams','v:0','-show_entries','stream=nb_read_frames,r_frame_rate,width,height,duration','-of','json',str(fig/'video.mp4')],stdout=subprocess.PIPE,check=True)
                        (out/'ffprobe.json').write_bytes(probe.stdout);record['video_probe']=json.loads(probe.stdout)
                        record['media_complete']=int(record['video_probe']['streams'][0]['nb_read_frames'])==51
                (out/'media-commands.json').write_text(json.dumps(media,indent=2)+'\n')
                (fig/'README.md').write_text(f'# {title}\n\nFixed-camera capability sample, NH E=44 MPa, nu=0.3, no contact. Accuracy unverified.\n\nVideo: 51 simulated frames sampled every {stride} solver steps; dt=0.005 s; {fps} fps, no interpolation. Physical interval 0–{end} s; expected container duration {51/fps:.2f} s. Numerical run covers {steps*.005} s. Red nodes mark the fixed left end.\n')
        (out/'visual-summary.json').write_text(json.dumps(record,indent=2)+'\n');all_results.append(record)
        print(json.dumps(record),flush=True)
    (w/'results.json').write_text(json.dumps(all_results,indent=2)+'\n')

if __name__=='__main__':main()
