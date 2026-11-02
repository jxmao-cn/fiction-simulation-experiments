#!/usr/bin/env python3
"""S3/S3t 结果图（英文标注），全部数据来自 run-id 目录下的 steps.csv。"""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/fiction-s3-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_steps(path):
    lines = Path(path).read_text().splitlines()
    hdr = [l for l in lines if not l.startswith('#')][0].split(',')
    out = {k: [] for k in hdr}
    for l in lines[2:]:
        parts = l.split(',')
        for k, v in zip(hdr, parts):
            try:
                out[k].append(float(v))
            except ValueError:
                out[k].append(float('nan'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--s3-dir', required=True, type=Path)
    ap.add_argument('--s3t-dir', required=True, type=Path)
    ap.add_argument('--output-dir', required=True, type=Path)
    ap.add_argument('--s3t-output-dir', type=Path, default=None)
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    s3t_out = a.s3t_output_dir or a.output_dir
    s3t_out.mkdir(parents=True, exist_ok=True)
    s3 = read_steps(a.s3_dir / 'steps.csv')
    s3sum = json.loads((a.s3_dir / 'summary.json').read_text())
    s3t = read_steps(a.s3t_dir / 'steps.csv')
    s3tsum = json.loads((a.s3t_dir / 'summary.json').read_text())

    # ---------------- S3: tip deflection & speed ----------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax[0].plot(s3['time'], [v * 1000 for v in s3['tip_deflection']], lw=1.0, color='#1f4e79')
    ax[0].axhline(s3['tip_deflection'][-1] * 1000, color='black', ls='--', lw=0.8,
                  label=f"final {s3['tip_deflection'][-1]*1000:.2f} mm")
    ax[0].set(xlabel='Time [s]', ylabel='Tip deflection [mm]',
              title='S3 self-weight cantilever: x=1 face mean downward displacement')
    ax[0].legend(fontsize=8)
    ax[1].semilogy(s3['time'], s3['max_speed'], lw=1.0, label='max node speed (all nodes)')
    ax[1].semilogy(s3['time'], s3['tip_mean_speed'], lw=0.9, alpha=.8, label='mean tip speed')
    ax[1].axhline(0.05, color='black', ls='--', lw=0.8, label='terminal-speed report line 0.05 m/s')
    ax[1].set(xlabel='Time [s]', ylabel='Speed [m/s]', title='S3 speed history (dt = 5 ms, 2000 steps)')
    ax[1].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=.2)
    fig.savefig(a.output_dir / 's3-tip-deflection-and-speed.png', dpi=170)
    fig.savefig(a.output_dir / 's3-tip-deflection-and-speed.svg')
    plt.close(fig)

    # ---------------- S3: time-step comparison ----------------
    sweep = a.s3_dir / 'dt-sweep'
    comp = json.loads((sweep / 'comparison.json').read_text())
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3), constrained_layout=True)
    styles = {'0.01': ('#b34c25', '-'), '0.005': ('#1f4e79', '-'), '0.0025': ('black', '--')}
    for dt, sub in (('0.01', 'dt-0.01'), ('0.005', 'dt-0.005'), ('0.0025', 'dt-0.0025')):
        d = read_steps(sweep / sub / 'steps.csv')
        c, ls = styles[dt]
        ax[0].plot(d['time'], [v * 1000 for v in d['tip_deflection']], ls, color=c, lw=1.0,
                   label=f'dt = {dt} s')
    for t in comp['comparison_times_s']:
        ax[0].axvline(t, color='gray', ls=':', lw=0.7)
    ax[0].set(xlabel='Time [s]', ylabel='Tip deflection [mm]',
              title='S3 time-step comparison (first 2 s, free oscillation)')
    ax[0].legend(fontsize=8, title='reference: dt = 0.0025 s', title_fontsize=8)

    dts = ['0.01', '0.005']
    times = comp['comparison_times_s']
    for dt in dts:
        ys = [comp['points'][str(t)]['values'][dt]['relative_diff_vs_ref'] * 100 for t in times]
        ax[1].plot(times, ys, 'o-', color=styles[dt][0], lw=1.0, ms=5, label=f'dt = {dt} s')
    ax[1].axhline(2.0, color='black', ls='--', lw=0.8, label='pre-registered 2% threshold')
    ax[1].set_xticks(times)
    ax[1].set(xlabel='Comparison time [s]', ylabel='|tip difference| vs dt = 0.0025 s  [%]',
              title='Pre-registered pointwise criterion: FAILED' if not comp['pass']
                    else 'Pre-registered pointwise criterion: passed')
    ax[1].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=.2)
    fig.savefig(a.output_dir / 's3-dt-comparison.png', dpi=170)
    fig.savefig(a.output_dir / 's3-dt-comparison.svg')
    plt.close(fig)

    # ---------------- S3: settled diagnostic ----------------
    pairs = [('0.01', sweep / 'diag-settled-dt0.01', 'post-hoc 1000 steps'),
             ('0.005', a.s3_dir, 'main 2000 steps'),
             ('0.0025', sweep / 'diag-settled-dt0.0025', 'post-hoc 4000 steps')]
    fig, ax = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
    recs = []
    for dt, d, lab in pairs:
        st = read_steps(Path(d) / 'steps.csv')
        tail = [(t, v) for t, v in zip(st['time'], st['tip_deflection']) if t > 8.0]
        mean_tail = sum(v for _, v in tail) / len(tail) * 1000
        max_sp = max(v for t, v in zip(st['time'], st['max_speed']) if t > 8.0)
        recs.append((dt, mean_tail, st['tip_deflection'][-1] * 1000, max_sp, lab))
    ref = [r for r in recs if r[0] == '0.0025'][0][1]
    ax.axhspan(ref * 0.98, ref * 1.02, color='green', alpha=.12, label='+/-2% of dt=0.0025 s')
    for i, (dt, mean_tail, endv, max_sp, lab) in enumerate(recs):
        ax.plot([i], [mean_tail], 'o', ms=8, color='#1f4e79')
        ax.plot([i], [endv], 'o', ms=8, mfc='none', mec='#b34c25')
        ax.annotate(f'{mean_tail:.3f} mm\n({100*(mean_tail-ref)/ref:+.2f}%)', (i, mean_tail),
                    xytext=(0, 13), textcoords='offset points', ha='center', fontsize=8)
    ax.set_xticks(range(len(recs)))
    ax.set_xticklabels([f"dt = {r[0]} s\n{r[4]}\n$v_{{max}}$ = {r[3]:.1e} m/s" for r in recs], fontsize=7.5)
    ax.plot([], [], 'o', color='#1f4e79', label='mean tip over last 2 s')
    ax.plot([], [], 'o', mfc='none', mec='#b34c25', label='instantaneous value at t = 10 s')
    ax.set_xlim(-0.45, 2.45)
    ax.set(xlabel='Simulation case', ylabel='Tip deflection [mm]',
           title='Post-hoc diagnostic: quasi-static limit vs dt (not an acceptance test)')
    ax.legend(fontsize=7.5, loc='lower left')
    ax.grid(alpha=.2)
    fig.savefig(a.output_dir / 's3-settled-vs-dt.png', dpi=170)
    fig.savefig(a.output_dir / 's3-settled-vs-dt.svg')
    plt.close(fig)
    json.dump([{'dt': r[0], 'mean_tip_last2s_mm': r[1], 'end_tip_mm': r[2],
                'max_speed_last2s': r[3], 'case': r[4]} for r in recs],
              open(sweep / 'settled-diagnostic.json', 'w'), indent=2)

    # ---------------- S3t: rotation, torque, tracking error ----------------
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.0), constrained_layout=True)
    ax[0].plot(s3t['time'], s3t['target_rotation_deg'], '--', color='black', lw=0.9,
               label='prescribed target')
    ax[0].plot(s3t['time'], s3t['actual_rotation_deg'], color='#1f4e79', lw=1.2,
               label='measured x=1 face rotation')
    ax[0].set(xlabel='Time [s]', ylabel='Rotation about x [deg]',
              title='S3t driven rotation (90 deg, 2 s smoothstep)')
    ax[0].legend(fontsize=8)
    ax[1].plot(s3t['time'], s3t['reaction_torque_x'], color='#b34c25', lw=1.2,
               label='wall reaction torque')
    ax[1].plot(s3t['time'], s3t['drive_torque_x'], color='#1f4e79', lw=1.0, ls='--',
               label='drive-end material torque')
    ax[1].set(xlabel='Time [s]', ylabel='Torque about x [N m]',
              title='S3t reaction / drive torque (no gravity)')
    ax[1].legend(fontsize=8)
    err_deg = [abs(t - m) for t, m in zip(s3t['target_rotation_deg'], s3t['actual_rotation_deg'])]
    ax[2].semilogy(s3t['time'], err_deg, color='black', lw=1.0, label='|target - measured| [deg]')
    ax[2].axhline(1e-6, color='gray', ls='--', lw=0.8, label='1e-6 m DBC gate (in degrees ~ 2.3e-3)')
    ax[2].set(xlabel='Time [s]', ylabel='Angle error [deg]',
              title='S3t rotation tracking error\n(max %.2e deg)' % max(err_deg))
    ax[2].legend(fontsize=7.5)
    for a_ in ax:
        a_.grid(alpha=.2)
    fig.savefig(s3t_out / 's3t-rotation-and-torque.png', dpi=170)
    fig.savefig(s3t_out / 's3t-rotation-and-torque.svg')
    plt.close(fig)

    # ---------------- S3t: principal stretches & J ----------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), constrained_layout=True)
    ax[0].plot(s3t['time'], s3t['min_stretch'], color='#1f4e79', lw=1.0, label='min principal stretch')
    ax[0].plot(s3t['time'], s3t['max_stretch'], color='#b34c25', lw=1.0, label='max principal stretch')
    ax[0].axhline(1.0, color='black', ls=':', lw=0.8)
    ax[0].set(xlabel='Time [s]', ylabel='Singular value of F [-]',
              title='S3t principal stretches (all tets)')
    ax[0].legend(fontsize=8)
    ax[1].plot(s3t['time'], s3t['min_J'], color='#1f4e79', lw=1.0, label='min det F')
    ax[1].plot(s3t['time'], s3t['max_J'], color='#b34c25', lw=1.0, label='max det F')
    ax[1].axhline(0.3, color='black', ls='--', lw=0.8, label='material-inversion gate J > 0.3')
    ax[1].set(xlabel='Time [s]', ylabel='det F [-]',
              title='S3t volume ratio range (min %.4f)' % min(s3t['min_J']))
    ax[1].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=.2)
    ax[1].set_ylim(0.0, max(1.15, max(s3t['max_J']) + 0.05))
    fig.savefig(s3t_out / 's3t-stretch-and-J.png', dpi=170)
    fig.savefig(s3t_out / 's3t-stretch-and-J.svg')
    plt.close(fig)

    print('S3 figures:', a.output_dir)
    print('S3t figures:', s3t_out)
    print('S3 pass:', s3sum['gates']['pass'], 'dt sweep pass:', comp['pass'])
    print('S3t pass:', s3tsum['gates']['pass'])


if __name__ == '__main__':
    main()
