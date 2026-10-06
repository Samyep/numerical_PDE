"""Typeset the four new panels from archived results; does not run PDE solvers."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / 'data' / 'batchwise_extended'
F = ROOT / 'figures'
headline = json.loads((D / 'hjb_headline_extracted.json').read_text())['rows']
suite = json.loads((D / 'batchwise_extended_suite_summary.json').read_text())
methods = [('raw', 'Raw (no correction)', 'o'),
           ('samplewise', 'Samplewise IR', 's'),
           ('batch', 'Batch-IR', '^')]
dimensions = [r['d'] for r in headline]

def newfig():
    return plt.subplots(figsize=(3.25, 2.35))

def save(fig, ax, name):
    ax.tick_params(labelsize=8)
    ax.xaxis.label.set_size(8.5)
    ax.yaxis.label.set_size(8.5)
    ax.grid(True, axis='y', alpha=.22)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, frameon=False, handlelength=1.5, labelspacing=.25)
    fig.tight_layout(pad=.7)
    fig.savefig(F / f'{name}.pdf', bbox_inches='tight', metadata={'Creator':'Archived Batch-IR result typesetting'})
    fig.savefig(F / f'{name}.png', dpi=180, bbox_inches='tight')
    plt.close(fig)

fig, ax = newfig()
for k,label,m in methods:
    ax.errorbar(dimensions, [r['relative_l2'][k]['mean'] for r in headline],
                yerr=[r['relative_l2'][k]['std'] for r in headline], marker=m,
                linewidth=1.25, markersize=3.8, capsize=2.5, label=label)
ax.set(xlabel='Dimension $d$', ylabel='Mean relative $L_2$ error', xticks=dimensions, yscale='log')
save(fig, ax, 'extended_hjb_headline')

fig, ax = newfig(); x = np.arange(len(suite['funding_deep']))
for j,(k,label,_) in enumerate(methods):
    ax.bar(x+(j-1)*.25, [r[f'{k}_mae'] for r in suite['funding_deep']], width=.25, label=label)
ax.set_xticks(x, [f"$n={r['n']},M={r['M']}$" for r in suite['funding_deep']])
ax.set(xlabel='Separate depth / budget settings', ylabel='Value MAE')
save(fig, ax, 'extended_funding_deep')

for metric,ylabel,name in [
    ('generator_mse', 'Generator MSE', 'extended_generator_mse'),
    ('correction_variance', 'Averaged-correction variance', 'extended_correction_variance')]:
    fig, ax = newfig()
    for k,label,m in methods:
        ax.plot(dimensions, [r['mechanism'][metric][k] for r in headline],
                marker=m, linewidth=1.25, markersize=3.8, label=label)
    ax.set(xlabel='Dimension $d$', ylabel=ylabel, xticks=dimensions, yscale='log')
    save(fig, ax, name)
print('Four archived-data panels written; original figures left unchanged.')
