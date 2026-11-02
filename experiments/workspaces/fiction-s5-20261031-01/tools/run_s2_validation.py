#!/usr/bin/env python3
"""固定配置的 S2 验收批次；失败即停，输出不得覆盖已有数值结果。"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--build', required=True, type=Path)
    p.add_argument('--run-dir', required=True, type=Path)
    args = p.parse_args()
    build, out = args.build.resolve(), args.run_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any((out / f'{m}.csv').exists() for m in ('M1', 'M2', 'M3')):
        p.error('已有结果，拒绝覆盖；请使用新的 run-id')
    env = dict(os.environ)
    env.update({key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')})
    config = {'protocol': 'fiction-validation-v0.2-S2', 'spec': 'docs/s2-execution-spec-20261031.md',
              'run_id': out.name, 'python': sys.executable, 'build': str(build),
              'material': {'E': 5.90e8, 'nu': .3, 'rho': 1000},
              'gravity': [0, 0, -9.81], 'solver': {'static': True, 'h': 1, 'tol': 1e-10,
              'max_iter': 100, 'initial_dbc_stiffness': 1e8},
              'tip_definition': 'arithmetic mean downward displacement of x=L nodes',
              'thresholds': {'min_J': .9, 'dbc_error_m': 1e-6,
                             'strictly_decreasing_relative_error': True, 'M3_relative_error': .05},
              'per_mesh_seconds': 600, 'batch_seconds': 1800, 'address_space_bytes': 8*1024**3,
              'environment': {key: env[key] for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                              'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}}
    (out/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    source_files = sorted([ROOT/'cpp/CMakeLists.txt', *[f for sub in ('cpp/include', 'cpp/src',
                            'cpp/examples', 'tools') for f in (ROOT/sub).rglob('*')
                            if f.is_file() and '__pycache__' not in f.parts],
                           ROOT/'docs/s2-execution-spec-20261031.md'])
    with tarfile.open(out/'source.tar.gz', 'w:gz') as tar:
        for f in source_files:
            tar.add(f, arcname=str(f.relative_to(ROOT)))
    hashes = {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files}
    binary = build/'cantilever_check'
    hashes[str(binary)] = hashlib.sha256(binary.read_bytes()).hexdigest()
    hashes['protocol_v0.2'] = hashlib.sha256((next(p for p in ROOT.parents if (p/'fiction.project.json').is_file())/'planning/physics-cpp-simulator-validation-protocol-v0.2.md').read_bytes()).hexdigest()
    (out/'source-hashes.json').write_text(json.dumps(hashes, indent=2)+'\n')
    shutil.copy2(binary, out/'cantilever_check')
    for cmd, name in [(['git', 'status', '--short'], 'source-status.txt'),
                      (['git', 'diff'], 'source.patch'), (['uname', '-a'], 'uname.txt'),
                      (['lscpu'], 'cpu.txt'), (['free', '-b'], 'memory-before.txt')]:
        with (out/name).open('w') as log:
            subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    with (out/'ctest.log').open('w') as log:
        check = subprocess.run(['ctest', '--test-dir', str(build), '--output-on-failure'],
                               env=env, stdout=log, stderr=subprocess.STDOUT)
    if check.returncode:
        (out/'summary.json').write_text(json.dumps({'status': 'regression_failed', 'returncode': check.returncode}, indent=2)+'\n')
        return 1

    def memory_limit():
        resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))

    started, rows, commands = time.monotonic(), [], []
    status = 'incomplete'
    for name, nx, ny in [('M1', 10, 2), ('M2', 20, 4), ('M3', 40, 8)]:
        mesh = out/'meshes'/f'{name}.msh'
        gen = [sys.executable, str(ROOT/'tools/make_structured_bar.py'), '--nx', str(nx),
               '--ny', str(ny), '--nz', str(ny), '--output', str(mesh),
               '--metadata', str(mesh.with_suffix('.json'))]
        with (out/'mesh-generation.log').open('a') as log:
            subprocess.run(gen, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = ['/usr/bin/time', '-v', '-o', str(out/f'{name}-time.log'),
               str(out/'cantilever_check'), '--mesh', str(mesh), '--csv', str(out/f'{name}.csv'),
               '--state', str(out/f'{name}-state.csv'), '--trace-dir', str(out/f'{name}-trace')]
        commands.append({'mesh_generation': gen, 'solve': cmd})
        (out/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
        remaining = min(600, 1800-(time.monotonic()-started))
        if remaining <= 0:
            status = 'batch_timeout'
            break
        # timeout also terminates the actual solver child, not only GNU time.
        with (out/f'{name}.log').open('w') as log:
            result = subprocess.run(['timeout', '--signal=TERM', '--kill-after=5', str(remaining),
                                     *cmd], env=env, preexec_fn=memory_limit,
                                    stdout=log, stderr=subprocess.STDOUT)
        (out/f'{name}.exitcode').write_text(str(result.returncode)+'\n')
        if (out/f'{name}.csv').exists():
            with (out/f'{name}.csv').open() as f:
                data = list(csv.DictReader(line for line in f if not line.startswith('#')))
            if data:
                rows.append(dict(mesh=name, **data[0]))
        print(f'{name}: exit={result.returncode}', flush=True)
        if result.returncode:
            status = 'numerical_or_resource_failure'
            break
        if not rows or rows[-1]['mesh'] != name or rows[-1]['numerical_pass'] != '1':
            status = 'missing_or_invalid_result'
            break
    else:
        errors = [float(r['relative_error']) for r in rows]
        monotonic = errors[0] > errors[1] > errors[2]
        accuracy = errors[-1] < .05
        status = 'passed' if monotonic and accuracy else 'accuracy_failed'
    summary = {'run_id': out.name, 'status': status, 'rows': rows,
               'wall_seconds': time.monotonic()-started, 'S3_started': False}
    if len(rows) == 3:
        errors = [float(r['relative_error']) for r in rows]
        summary.update(strictly_decreasing_error=errors[0] > errors[1] > errors[2],
                       M3_below_5_percent=errors[-1] < .05)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('S2 batch:', status)
    return 0 if status == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
