"""离线汇总三个组别；不运行求解器。"""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_csv(p):
    with p.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))


summary = {'groups': {}, 'same_initial_state': False}
groups = ['control', 'iteration300', 'penalty1e10']
for name in groups:
    folder = ROOT / 'results' / name
    metrics = read_csv(folder / 'metrics.csv')
    trace = read_csv(folder / 'load-1.05.csv')
    last = trace[-1]
    tail = [row for row in trace if int(row['iter']) >= 100]
    summary['groups'][name] = {
        'metrics_final': metrics[-1],
        'trace_final': last,
        'exitcode': int((folder / 'exitcode').read_text()),
        'initial_sha256': hashlib.sha256((folder / 'final-load-initial.csv').read_bytes()).hexdigest(),
        'trace_states': len(trace),
        'zero_coordinate_updates': sum(row['trial_equals_x'] == '1' for row in trace if int(row['iter']) > 0),
        'tail100_states': len(tail),
        'tail100_distinct_energy': len({row['energy'] for row in tail}),
        'tail100_dbc_error_range': max([float(row['dbc_error_max']) for row in tail], default=0)
                                    - min([float(row['dbc_error_max']) for row in tail], default=0),
    }
summary['same_initial_state'] = len({x['initial_sha256'] for x in summary['groups'].values()}) == 1
assert summary['same_initial_state']
assert summary['groups']['control']['metrics_final']['pass'] == '0'
assert summary['groups']['iteration300']['metrics_final']['pass'] == '0'
assert summary['groups']['penalty1e10']['metrics_final']['pass'] == '1'
(ROOT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

# 科学图采用英文；已满足全部节点时，未满足集合为空，不把比例0画成力完全为零。
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
labels = {'control': 'Control: 100, initial k=1e8',
          'iteration300': 'Iteration group: 300, initial k=1e8',
          'penalty1e10': 'Penalty group: 100, initial k=1e10'}
for name in groups:
    trace = read_csv(ROOT / 'results' / name / 'load-1.05.csv')
    x = [int(row['iter']) for row in trace]
    axes[0, 0].semilogy(x, [float(row['dbc_error_max']) for row in trace], label=labels[name])
    axes[0, 1].semilogy(x, [float(row['k']) for row in trace])
    axes[1, 0].semilogy(x[1:], [float(row['actual_dx_inf']) or math.nan for row in trace[1:]])
    axes[1, 1].semilogy(x, [float(row['force_balance_relative_unsatisfied'])
                                 if int(row['unsatisfied_nodes']) else math.nan for row in trace])
axes[0, 0].axhline(1e-10, color='black', linestyle=':', label='DBC tracking tolerance')
axes[0, 0].set_ylabel('Max DBC position error (m)')
axes[0, 1].set_ylabel('Adaptive penalty parameter k')
axes[1, 0].set_ylabel('Accepted coordinate update, infinity norm (m)')
axes[1, 1].set_ylabel('Relative net force on unsatisfied DBC nodes')
axes[0, 0].legend(fontsize=8)
for ax in axes.flat:
    ax.set_xlabel('Newton iteration')
    ax.grid(True, alpha=0.25)
fig.suptitle('M2, 5% stretch: identical incoming state; variant-A in all groups')
fig.savefig(ROOT / 'comparison.png', dpi=170)
plt.close(fig)
print(json.dumps({name: {'pass': row['metrics_final']['pass'],
                       'iterations': row['metrics_final']['iterations'],
                       'boundary_error': row['trace_final']['dbc_error_max'],
                       'balance_ratio': row['trace_final']['force_balance_relative_unsatisfied']}
                  for name, row in summary['groups'].items()}, indent=2))
