"""Show descriptive cross-age signature resemblance, with within-stage cells masked."""
from pathlib import Path
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/ventricular_stage_subtype_catalog'
STAGES = ['P0', 'P7', 'P14', 'P21']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
catalog = pd.read_csv(OUT / 'stage_subtype_catalog.tsv', sep='\t')
ids = [cluster for stage in STAGES for cluster in catalog.loc[catalog.stage.eq(stage), 'microenvironment_cluster']]
counts = catalog.set_index('microenvironment_cluster').nuclei
families = catalog.set_index('microenvironment_cluster').functional_family
labels = [f'{cluster}  {counts[cluster]:,}' for cluster in ids]
colors = {'Mechanical-response-associated': '#bd5b17', 'Conduction-program-enriched': '#8f3985', 'Regional electrical-associated': '#147f9c', 'Cycling': '#1f944d'}
pdf = PdfPages(OUT / 'cross_stage_signature_similarity.pdf')
for filename, title in [('relative_signature_similarity_matrix.tsv', 'Relative gene signatures'), ('interpreted_genes_excluded_similarity_matrix.tsv', 'Previously interpreted genes withheld')]:
    matrix = pd.read_csv(OUT / filename, sep='\t', index_col=0).loc[ids, ids]
    same_stage = np.array([[a.split('_')[0] == b.split('_')[0] for b in ids] for a in ids])
    data = np.ma.masked_where(same_stage, matrix.to_numpy())
    fig, ax = plt.subplots(figsize=(14, 13))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#eeeeee')
    im = ax.imshow(data, cmap=cmap, vmin=-0.7, vmax=0.7, interpolation='none')
    ax.set_xticks(range(len(ids)), ids, rotation=65, ha='right')
    ax.set_yticks(range(len(ids)), labels)
    for tick, cluster in zip(ax.get_yticklabels(), ids):
        tick.set_color(colors.get(families[cluster], '#333333'))
    offset = 0
    for stage in STAGES:
        n = sum(cluster.startswith(stage + '_') for cluster in ids)
        if offset:
            ax.axhline(offset - .5, color='black', linewidth=1)
            ax.axvline(offset - .5, color='black', linewidth=1)
        ax.text(offset + (n - 1) / 2, -1.4, stage, ha='center', weight='bold', fontsize=12)
        offset += n
    fig.colorbar(im, ax=ax, shrink=.68, pad=.025, label='Spearman correlation across relative gene signatures', extend='both')
    ax.set_title(title, pad=37, fontsize=17, weight='bold')
    fig.text(.5, .025, 'Group count-sum log1p(CPM) minus same-stage remainder; descriptive resemblance, not lineage or differential expression.\n'
             'Gray blocks mask within-stage comparisons. Numbers beside groups are retained nuclei; no forced one-to-one matches.\n'
             'Orange: Ankrd1-associated; purple: conduction-program candidates; blue: regional electrical candidates; green: cycling.\n'
             'Pooled signatures can reflect library weighting, age and stage composition. Small/marker-poor groups retain their review flags.',
             ha='center', fontsize=10, linespacing=1.45)
    fig.subplots_adjust(left=.17, bottom=.18, top=.90, right=.94)
    suffix = '' if filename.startswith('relative') else '_interpreted_genes_excluded'
    fig.savefig(OUT / f'cross_stage_signature_similarity{suffix}.png', dpi=190, facecolor='white')
    pdf.savefig(fig, facecolor='white')
    plt.close(fig)
pdf.close()
print('Saved primary/sensitivity cross-stage heatmaps and two-page PDF.')
