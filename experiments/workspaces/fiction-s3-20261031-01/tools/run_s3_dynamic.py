#!/usr/bin/env python3
"""S3/S3t 单次运行封装：记录配置/命令/源码哈希/二进制哈希/时间与 RSS/退出码。

用法示例：
  run_s3_dynamic.py --run-dir DIR --mode bend --dt 0.005 --steps 2000 \
      --tol 1e-4 --max-iter 50 --mesh FILE [--tag NAME] [--twist-deg 90 --twist-time 2.0]
"""
import argparse
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = Path('/tmp/fiction-s3-build/dynamic_bar_check')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', required=True, type=Path)
    ap.add_argument('--mesh', required=True, type=Path)
    ap.add_argument('--mode', default='bend', choices=['bend', 'twist'])
    ap.add_argument('--dt', default=0.005, type=float)
    ap.add_argument('--steps', default=2000, type=int)
    ap.add_argument('--tol', default=1e-4, type=float)
    ap.add_argument('--max-iter', default=50, type=int)
    ap.add_argument('--twist-deg', default=90.0, type=float)
    ap.add_argument('--twist-time', default=2.0, type=float)
    ap.add_argument('--gravity-z', default=-9.81, type=float)
    ap.add_argument('--tag', default='')
    ap.add_argument('--timeout', default=3600, type=int)
    a = ap.parse_args()

    out = a.run_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'steps.csv').exists():
        sys.exit(f'拒绝覆盖已有结果: {out}/steps.csv')

    env = dict(os.environ)
    env.update({k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')})

    cmd = [str(BINARY), '--mesh', str(a.mesh.resolve()), '--out-dir', str(out),
           '--dt', repr(a.dt), '--steps', str(a.steps), '--tol', repr(a.tol),
           '--max-iter', str(a.max_iter), '--mode', a.mode,
           '--twist-deg', repr(a.twist_deg), '--twist-time', repr(a.twist_time),
           '--gravity-z', repr(a.gravity_z)]
    (out / 'command.txt').write_text(' '.join(cmd) + '\n')
    config = {
        'protocol': 'fiction-validation-v0.3-S3',
        'spec': 'docs/s3-s6a-execution-design-20261031.md',
        'run_id': out.name, 'tag': a.tag, 'mode': a.mode,
        'mesh': str(a.mesh.resolve()), 'mesh_sha256': sha(a.mesh),
        'binary': str(BINARY), 'binary_sha256': sha(BINARY),
        'dt': a.dt, 'steps': a.steps, 'sim_time': a.dt * a.steps,
        'tol': a.tol, 'max_iter': a.max_iter,
        'initial_dbc_stiffness': 1e8,
        'gravity_z': 0.0 if a.mode == 'twist' else a.gravity_z,
        'twist_deg': a.twist_deg if a.mode == 'twist' else 0.0,
        'twist_time': a.twist_time if a.mode == 'twist' else 0.0,
        'material': {'E': 4.40e7, 'nu': 0.3, 'rho': 1000},
        'build': '/tmp/fiction-s3-build', 'build_type': 'Release', 'solver': 'CHOLMOD',
        'environment': {k: env[k] for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                            'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')},
        'python': sys.executable,
    }
    (out / 'config.json').write_text(json.dumps(config, indent=2) + '\n')

    files = sorted([ROOT / 'cpp/CMakeLists.txt',
                    *[f for sub in ('cpp/include', 'cpp/src', 'cpp/examples', 'tools')
                      for f in (ROOT / sub).rglob('*')
                      if f.is_file() and '__pycache__' not in f.parts]])
    hashes = {str(f.relative_to(ROOT)): sha(f) for f in files}
    with tarfile.open(out / 'source.tar.gz', 'w:gz') as tar:
        for f in files:
            tar.add(f, arcname=str(f.relative_to(ROOT)))
    (out / 'source-hashes.json').write_text(json.dumps(hashes, indent=2) + '\n')
    ctest_log = Path('/tmp/fiction-s3-build/ctest.log')
    if ctest_log.exists():
        shutil.copy2(ctest_log, out / 'ctest.log')

    t0 = time.time()
    usage0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (out / 'run.log').open('w') as log:
        proc = subprocess.run(cmd, cwd='/tmp', env=env, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, timeout=a.timeout)
        log.write(proc.stdout)
    wall = time.time() - t0
    usage1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    (out / 'exitcode').write_text(f'{proc.returncode}\n')
    timeinfo = {
        'wall_seconds': wall,
        'child_user_seconds': usage1.ru_utime - usage0.ru_utime,
        'child_sys_seconds': usage1.ru_stime - usage0.ru_stime,
        'child_max_rss_kb': usage1.ru_maxrss,
        'exitcode': proc.returncode,
    }
    (out / 'time-resource.json').write_text(json.dumps(timeinfo, indent=2) + '\n')
    print(json.dumps(timeinfo))
    print(proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else '(no output)')
    sys.exit(0 if proc.returncode == 0 else 1)


if __name__ == '__main__':
    main()
