"""Plot raw versus supported CytoSPACE assignments and withheld-gene recovery."""
import os
import argparse
os.environ.setdefault('MPLCONFIGDIR', '/tmp/heart-mpl')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from run_cytospace_xenium import SOURCE, OUT, STAGES
from plot_xenium_ventricular import background, frame, COLORS, fine_colors


def scatter(ax, obs, field, colors, vent, descriptions=None):
    ax.scatter(obs.x_centroid, obs.y_centroid, c='#e6e6e6', s=.35, edgecolors='none', rasterized=True)
    labels = obs.loc[vent, field].fillna('Not supported by all checks')
    handles = []
    for label in sorted(labels.unique()):
        mask = vent & obs[field].eq(label)
        if label == 'Not supported by all checks':
            continue
        color = colors.get(label, '#333333')
        ax.scatter(obs.loc[mask, 'x_centroid'], obs.loc[mask, 'y_centroid'], c=[color], s=1.0, edgecolors='none', alpha=.9, rasterized=True)
        display = f'{label} ({mask.sum():,})'
        if descriptions and label in descriptions:
            display = f'{label}: {descriptions[label]}\n({mask.sum():,})'
        handles.append(Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markersize=5, label=display))
    if handles:
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False, borderaxespad=0)
    else:
        ax.text(.5, .5, 'No groups pass all checks', transform=ax.transAxes, ha='center', va='center', fontsize=13,
                bbox={'facecolor': 'white', 'alpha': .9, 'edgecolor': 'none'})


def save(fig, folder, name, pdf):
    fig.savefig(folder / f'{name}.png', dpi=180, facecolor='white')
    fig.savefig(folder / f'{name}.pdf', facecolor='white')
    pdf.savefig(fig, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stages', nargs='+', default=STAGES, choices=STAGES)
    args = parser.parse_args()
    catalog = pd.read_csv(SOURCE / 'reference_group_catalog.tsv', sep='\t')
    families = {**COLORS, 'Uncharacterized ventricular candidate': '#6e6e6e'}
    fine = fine_colors(catalog)
    group_ids = [f'P{age}_M{i}' for age, n in [(0, 8), (7, 9), (14, 5), (21, 6)] for i in range(n)]
    short = dict(zip(group_ids, [
        'cycling', 'Ddc-associated metabolic', 'Zbtb16/Ncam1 working', 'Ankrd1 remodeling A',
        'fatty-acid-associated', 'Slit2/conduction program', 'Ankrd1 remodeling B', 'DNA maintenance',
        'Ppara/Rreb1 working', 'Ankrd1 remodeling', 'Fgf14/Pdzrn4 electrical', 'cycling',
        'Ivns1abp/Myom2 working', 'DNA maintenance', 'Slit2/conduction program', 'uncharacterized A', 'uncharacterized B',
        'Slit2/Mpped2 electrical', 'Lrtm1 working', 'Ankrd1 remodeling', 'Kcnd2/Pdzrn4 reference', 'Slit2/conduction program',
        'Lrtm1/Stk39 reference', 'marker-poor working', 'Ankrd1 remodeling A', 'Slit2/Mpped2 electrical',
        'Ankrd1 remodeling B', 'Slit2/conduction program']))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
    with PdfPages(OUT / 'cytospace_spatial_maps_and_checks.pdf') as pdf:
        for stage in args.stages:
            folder = OUT / stage
            obs = pd.read_csv(folder / 'cytospace_cell_annotations.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
            vent = obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
            bg, extent, bounds = background(stage, obs)
            fig, axes = plt.subplots(2, 2, figsize=(21, 17))
            for ax, field, title, colors in [
                (axes[0, 0], 'cyto_corrected_seed42_family', 'Raw CytoSPACE families: every included cell is assigned', families),
                (axes[0, 1], 'cyto_stable_family_candidate', 'Families passing panel, seed, library and count checks', families),
                (axes[1, 0], 'cyto_corrected_seed42_group', 'Raw reference groups: exploratory forced assignments', fine),
                (axes[1, 1], 'cyto_stable_group_candidate', 'Fine groups passing all checks', fine)]:
                frame(ax, bg, extent, bounds, title)
                scatter(ax, obs, field, colors, vent, short if 'group' in field else None)
            fig.suptitle(f'{stage}: CytoSPACE assignments on native Xenium coordinates', fontsize=21, weight='bold', y=.986)
            fig.text(.5, .045, 'Broad identities are supplied as input; subtype labels/fractions are not supplied. Native assignments are compulsory, not confidence-filtered.\n'
                     'Checks: same-stage reference-panel cross-library recoverability, 3 seeds, original counts, each reference library, expression fit and measured marker support.\n'
                     'Gray cells retain their original coordinates. Kcnd2 and Irx5 are unassayed; no measured electrical gradient is shown.', ha='center', fontsize=10)
            fig.subplots_adjust(left=.055, right=.84, top=.948, bottom=.11, wspace=.74, hspace=.26)
            save(fig, folder, f'{stage}_cytospace_raw_and_supported_maps', pdf)
        metrics = pd.read_csv(OUT / 'heldout_gene_recovery.tsv', sep='\t')
        genes = metrics.gene.drop_duplicates().tolist()
        fig, axes = plt.subplots(1, 3, figsize=(16, 10))
        for ax, key, title, limits in [
            (axes[0], 'cell_spearman_rho', 'Withheld genes: cell-level Spearman r', (-.5, .5)),
            (axes[1], 'depth_area_gene_count_adjusted_r', 'Withheld genes: depth/area-adjusted r', (-.5, .5)),
            (axes[2], 'grid250um_spearman_rho', 'Withheld genes: 250 µm grid Spearman r', (-1, 1))]:
            values = metrics.pivot(index='gene', columns='stage', values=key).reindex(index=genes, columns=STAGES)
            im = ax.imshow(values, cmap='RdBu_r', vmin=limits[0], vmax=limits[1], aspect='auto')
            ax.set_xticks(np.arange(4), STAGES)
            ax.set_yticks(np.arange(len(genes)), genes)
            ax.set_title(title, fontsize=11, weight='bold')
            for i in range(len(genes)):
                for j in range(4):
                    value = values.iloc[i, j]
                    ax.text(j, i, 'NA' if pd.isna(value) else f'{value:.2f}', ha='center', va='center', fontsize=8,
                            color='white' if abs(value) > .75 * limits[1] else 'black')
            fig.colorbar(im, ax=ax, fraction=.035, pad=.025)
        fig.suptitle('Gene recovery when all listed genes are omitted from CytoSPACE matching', fontsize=16, weight='bold')
        fig.text(.5, .025, 'Evaluation is within the fixed ventricular input cohort. Broad input annotations are retained. Grid bins require ≥20 cells.\n'
                 'These are reference-derived predictions compared with measured Xenium genes; they do not validate unassayed Kcnd2 or Irx5.', ha='center', fontsize=10)
        fig.subplots_adjust(top=.91, bottom=.10, wspace=.48, left=.06, right=.97)
        save(fig, OUT, 'heldout_gene_recovery', pdf)


if __name__ == '__main__':
    main()
