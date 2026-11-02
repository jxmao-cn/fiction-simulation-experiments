#!/usr/bin/env python3
"""从 dt-sweep 的 steps.csv 复算预登记 dt 对照与事后准静态诊断（不改变判据）。"""
import json
import sys
from pathlib import Path


def read_steps(path):
    lines = Path(path).read_text().splitlines()
    hdr = [l for l in lines if not l.startswith('#')][0].split(',')
    return [dict(zip(hdr, l.split(','))) for l in lines[2:]]


def main():
    sweep = Path(sys.argv[1]).resolve()
    dts = {'0.01': 'dt-0.01', '0.005': 'dt-0.005', '0.0025': 'dt-0.0025'}
    data = {k: read_steps(sweep / v / 'steps.csv') for k, v in dts.items()}
    times = [0.5, 1.0, 1.5, 2.0]
    comp = {'reference_dt': 0.0025, 'comparison_times_s': times,
            'tolerance_relative': 0.02, 'criterion': 'pointwise relative tip difference vs dt=0.0025',
            'points': {}, 'pass': True}
    for t in times:
        vals = {}
        for k, rows in data.items():
            i = min(range(len(rows)), key=lambda j: abs(float(rows[j]['time']) - t))
            vals[k] = (float(rows[i]['time']), float(rows[i]['tip_deflection']))
        ref = vals['0.0025'][1]
        rec = {'reference_tip_deflection_m': ref, 'values': {}}
        for k in ('0.01', '0.005', '0.0025'):
            tv, v = vals[k]
            rel = abs(v - ref) / abs(ref)
            rec['values'][k] = {'time': tv, 'tip_deflection_m': v,
                                'relative_diff_vs_ref': rel, 'within_2pct': rel < 0.02}
            if k != '0.0025' and rel >= 0.02:
                comp['pass'] = False
        comp['points'][str(t)] = rec
    (sweep / 'comparison.json').write_text(json.dumps(comp, indent=2) + '\n')

    pairs = [('0.01', sweep / 'diag-settled-dt0.01', 'post-hoc 1000 steps'),
             ('0.005', sweep.parent, 'main 2000 steps'),
             ('0.0025', sweep / 'diag-settled-dt0.0025', 'post-hoc 4000 steps')]
    recs = []
    for dt, d, lab in pairs:
        rows = read_steps(Path(d) / 'steps.csv')
        tail = [(float(r['time']), float(r['tip_deflection'])) for r in rows
                if float(r['time']) > 8.0]
        mean_tail = sum(v for _, v in tail) / len(tail) * 1000
        max_sp = max(float(r['max_speed']) for r in rows if float(r['time']) > 8.0)
        recs.append({'dt': dt, 'case': lab, 'mean_tip_last2s_mm': mean_tail,
                     'end_tip_mm': float(rows[-1]['tip_deflection']) * 1000,
                     'max_speed_last2s': max_sp})
    ref = [r for r in recs if r['dt'] == '0.0025'][0]['mean_tip_last2s_mm']
    for r in recs:
        r['relative_diff_of_window_mean_pct'] = 100 * (r['mean_tip_last2s_mm'] - ref) / ref
    (sweep / 'settled-diagnostic.json').write_text(json.dumps(
        {'note': 'post-hoc diagnostic, not an acceptance test', 'records': recs}, indent=2) + '\n')
    print('dt criterion pass:', comp['pass'])
    for t in times:
        print(t, {k: round(v['relative_diff_vs_ref'] * 100, 3) for k, v in comp['points'][str(t)]['values'].items()})
    for r in recs:
        print(r['dt'], round(r['mean_tip_last2s_mm'], 4), round(r['relative_diff_of_window_mean_pct'], 4))


if __name__ == '__main__':
    main()
