#!/usr/bin/env python3
"""从冻结 S2 summary 生成英文标注精度图；失败结果显式标记。"""
import argparse
import json
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/fiction-s2-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run-dir', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    a = p.parse_args()
    summary = json.loads((a.run_dir/'summary.json').read_text())
    rows = summary.get('rows', [])
    if not rows:
        raise SystemExit('No numeric data; no figure created')
    a.output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    names = [r['mesh'] + ('\n(not converged)' if r['converged'] != '1' else '') for r in rows]
    tips = [float(r['tip_deflection'])*1000 for r in rows]
    errs = [float(r['relative_error'])*100 for r in rows]
    axes[0].plot(names, tips, 'o-', label='P1 tetrahedra')
    axes[0].axhline(float(rows[0]['reference_tip'])*1000, color='black', ls='--',
                    label='Beam reference (bending + shear)')
    axes[0].set(ylabel='Tip deflection [mm]', xlabel='Mesh refinement', title='Self-weight cantilever')
    axes[0].legend(fontsize=8)
    axes[1].plot(names, errs, 'o-', color='#b34c25')
    axes[1].axhline(5, color='black', ls='--', label='M3 acceptance threshold: 5%')
    axes[1].set(ylabel='Relative error [%]', xlabel='Mesh refinement', title='Error against beam reference')
    axes[1].legend(fontsize=8)
    for i, (tip, err) in enumerate(zip(tips, errs)):
        axes[0].annotate(f'{tip:.4f}', (i, tip), xytext=(0, 7), textcoords='offset points', ha='center', fontsize=8)
        axes[1].annotate(f'{err:.2f}%', (i, err), xytext=(0, 7), textcoords='offset points', ha='center', fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
        ax.margins(y=.18)
    title = 'Validation passed' if summary['status'] == 'passed' else 'Validation failed; inspect numerical status'
    fig.suptitle(f'S2 | {title} | tip = end-face node mean', fontsize=11)
    fig.savefig(a.output_dir/'s2-convergence.png', dpi=180)
    fig.savefig(a.output_dir/'s2-convergence.svg')
    plt.close(fig)
    (a.output_dir/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(a.output_dir/'s2-convergence.png')


if __name__ == '__main__':
    main()
