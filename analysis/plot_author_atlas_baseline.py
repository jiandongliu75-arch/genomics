"""Plot retrieved published CM states and the independent parameter baseline."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import anndata as ad
import numpy as np
import pandas as pd

from acquire_author_heart_atlas import OUT, PREVIOUS, STAGES
from run_author_parameter_baseline import VENT_STATES
from analyze_author_atlas_baseline import CM_STATES, MARKERS
from plot_xenium_ventricular import background, frame, COLORS

STATE_COLORS = {'v.CM': '#7393ad', 'Ankrd1.CM': '#e67e22', 'Slit2.CM': '#22a7a5',
    'Ddc.CM': '#3275b9', 'p.CM': '#3baf52', 'a.CM': '#d56c9d', 'av.CM': '#8b62ac'}


def categorical(ax, obs, labels, colors, mask):
    ax.scatter(obs.x_centroid, obs.y_centroid, s=.3, color='#e6e6e6', edgecolors='none', rasterized=True)
    handles = []
    for label in sorted(labels[mask].unique()):
        select = mask & labels.eq(label)
        color = colors.get(label, '#555555')
        ax.scatter(obs.loc[select, 'x_centroid'], obs.loc[select, 'y_centroid'], s=.85, color=color, edgecolors='none', rasterized=True)
        handles.append(Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markersize=5, label=f'{label} ({select.sum():,})'))
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)


def save(fig, folder, name, pdf):
    fig.savefig(folder / f'{name}.png', dpi=180, facecolor='white')
    fig.savefig(folder / f'{name}.pdf', facecolor='white')
    pdf.savefig(fig, facecolor='white')
    plt.close(fig)


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
    with PdfPages(OUT / 'published_CM_maps_and_independent_baseline.pdf') as pdf:
        for stage in STAGES:
            folder = OUT / stage
            obs = pd.read_csv(folder / 'original_cells_with_author_labels.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
            pairs = pd.read_csv(folder / 'independent_default_parameters/annotated_assignments.tsv.gz', sep='\t').set_index('SpotID')
            obs['baseline_reference_family'] = pairs.reference_mapping_family.reindex(obs.index)
            cm = obs.author_atlas_cell_state.isin(CM_STATES)
            vent = obs.author_atlas_cell_state.isin(VENT_STATES)
            bg, extent, bounds = background(stage, obs)
            fig, axes = plt.subplots(2, 2, figsize=(21, 17))
            frame(axes[0, 0], bg, extent, bounds, 'Published author CM states: retrieved final annotations')
            categorical(axes[0, 0], obs, obs.author_atlas_cell_state, STATE_COLORS, cm)
            frame(axes[0, 1], bg, extent, bounds, 'Published ventricular states + chamber-unresolved p.CM')
            categorical(axes[0, 1], obs, obs.author_atlas_cell_state, STATE_COLORS, vent)
            frame(axes[1, 0], bg, extent, bounds, 'Independent reference families: raw/default 10,000-cell run')
            categorical(axes[1, 0], obs, obs.baseline_reference_family, COLORS, obs.baseline_reference_family.notna())
            frame(axes[1, 1], bg, extent, bounds, 'Measured Slit2 RNA in the selected CM cohort')
            data = ad.read_h5ad(PREVIOUS / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad', backed='r')
            counts = data[:, ['Slit2']].X.toarray().ravel()
            axes[1, 1].scatter(obs.loc[vent, 'x_centroid'], obs.loc[vent, 'y_centroid'], s=.45, color='#dddddd', edgecolors='none', rasterized=True)
            positive = vent.to_numpy() & (counts > 0)
            sc = axes[1, 1].scatter(obs.loc[positive, 'x_centroid'], obs.loc[positive, 'y_centroid'], c=np.log1p(counts[positive]),
                s=1.2, cmap='viridis', norm=Normalize(0, np.log1p(10)), edgecolors='none', rasterized=True)
            cb = fig.colorbar(sc, ax=axes[1, 1], fraction=.035, pad=.02, ticks=np.log1p([0, 1, 3, 10]))
            cb.ax.set_yticklabels(['0', '1', '3', '≥10'])
            cb.set_label('Original Slit2 molecules', fontsize=8)
            data.file.close()
            fig.suptitle(f'{stage}: published CM states and an independent CytoSPACE parameter comparison', fontsize=20, weight='bold', y=.986)
            fig.text(.5, .043, 'Top: authors’ published final labels, verified against original barcodes, centroids and transcript totals. These are not newly inferred identities.\n'
                'Bottom left: our independent snRNA reference, raw counts, default seed/downsampling and 10,000-cell partitions; no ventricular subtype labels supplied.\n'
                'The two annotation systems differ. v.CM is one working-ventricular state; p.CM includes chamber-unresolved cycling CMs.\n'
                'Kcnd2 and Irx5 are unassayed; TF_Irx5 is a derived feature. Positions and measured counts are unchanged.', ha='center', fontsize=9)
            fig.subplots_adjust(left=.055, right=.84, top=.947, bottom=.12, wspace=.74, hspace=.26)
            save(fig, folder, f'{stage}_published_CM_states_and_baseline', pdf)
        profile = pd.read_csv(OUT / 'measured_marker_profiles_by_author_CM_state.tsv', sep='\t')
        genes = [g for g in MARKERS if g in set(profile.gene)]
        profile['gene_display_z'] = profile.groupby('gene').mean_measured_log1pCP10k.transform(lambda x: (x-x.mean())/max(x.std(), 1e-8))
        fig, axes = plt.subplots(2, 2, figsize=(19, 11))
        for stage, ax in zip(STAGES, axes.ravel()):
            table = profile[profile.stage.eq(stage)]
            for row in table.itertuples():
                if row.gene not in genes:
                    continue
                x, y = genes.index(row.gene), CM_STATES.index(row.author_CM_state)
                im = ax.scatter(x, y, s=12+180*row.measured_positive_fraction, c=[row.gene_display_z], cmap='RdBu_r', vmin=-2, vmax=2, edgecolors='none')
            ax.set_xticks(range(len(genes)), genes, rotation=55, ha='right', fontsize=8)
            ax.set_yticks(range(len(CM_STATES)), CM_STATES)
            ax.set_ylim(len(CM_STATES)-.5, -.5)
            ax.set_xlim(-.5, len(genes)-.5)
            ax.set_title(stage, weight='bold')
            ax.grid(alpha=.12)
            fig.colorbar(im, ax=ax, fraction=.03, pad=.02, label='Within-gene display z')
        fig.suptitle('Directly measured Xenium marker profiles under the published CM-state labels', fontsize=17, weight='bold')
        fig.text(.5, .025, 'Dot area reflects detection fraction. Color is standardized mean log1p(CP10k) within each gene across state/stage summaries.\n'
            'Blank rows indicate no published segments in that state. These are descriptive profiles, not independent label validation or animal-replicated differential expression.', ha='center', fontsize=10)
        fig.subplots_adjust(left=.07, right=.96, top=.91, bottom=.16, wspace=.30, hspace=.68)
        save(fig, OUT, 'measured_markers_by_published_CM_state', pdf)


if __name__ == '__main__':
    main()
