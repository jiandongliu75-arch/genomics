"""Plot the descriptive library-level signals used in the fibroblast feasibility audit."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/fibroblast-audit-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / 'results/fibroblast_reassessment'
df = pd.read_csv(OUT / 'fibroblast_marker_profiles_by_library.tsv.gz', sep='\t')
df = df[df['subset'].eq('all_reviewed_fibroblasts')]
stages = ['P0', 'P7', 'P14', 'P21']
colors = ['#2678b2', '#ed8d23', '#35934b', '#b94a60']
genes = ['Postn', 'Tnc', 'Col8a1', 'Gfpt2', 'Fn1', 'Col3a1', 'Lamc1', 'Hspg2']
fig, axes = plt.subplots(2, 4, figsize=(14, 7), constrained_layout=True)
for ax, gene in zip(axes.flat, genes):
    sub = df[df.gene.eq(gene)]
    means = sub.groupby(['layer', 'age'], observed=True).mean_log1p_CP10k.mean()
    ax.plot(range(4), [means['corrected', s] for s in stages], color='#243c50', lw=1.7,
            label='Corrected: mean of libraries')
    ax.plot(range(4), [means['raw', s] for s in stages], color='#999999', lw=1.2, ls='--',
            label='Raw: mean of libraries')
    for i, stage in enumerate(stages):
        vals = sub[sub.layer.eq('corrected') & sub.age.eq(stage)].sort_values('sample').mean_log1p_CP10k.to_numpy()
        assert len(vals) == 2
        ax.scatter(i + np.array([-.07, .07]), vals, color=colors[i], s=43, edgecolor='white', lw=.6, zorder=3)
    ax.set_title(gene, fontsize=14, fontstyle='italic')
    ax.set_xticks(range(4), stages)
    ax.set_ylabel('Mean log1p(CP10k)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', alpha=.15)
axes[0, 0].legend(fontsize=8, loc='best')
fig.suptitle('Fibroblast expression programs across postnatal stages\n'
             '14,985 reviewed nuclei; dots show the two libraries at each stage', fontsize=15)
fig.savefig(OUT / 'fibroblast_library_signal_audit.png', dpi=170)
fig.savefig(OUT / 'fibroblast_library_signal_audit.pdf')
