"""Plot independently subclustered ventricular cells with support status kept visible.

Inputs are produced by the within-stage analysis; this script never fits or alters
an embedding, cluster, state label or support status.
"""
from pathlib import Path
import argparse
import os
import re
import textwrap

os.environ.setdefault('MPLCONFIGDIR', '/tmp/ventricular-stage-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
STAGES = ('P0', 'P7', 'P14', 'P21')
XY = ('within_stage_UMAP1', 'within_stage_UMAP2')
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00',
          '#6A3D9A', '#17BECF', '#B15928', '#E31A1C', '#526D35',
          '#4D4D4D', '#9467BD', '#1B9E77', '#A6761D', '#E7298A']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                     'pdf.fonttype': 42, 'axes.titlesize': 15})


def cluster_key(value):
    return tuple(int(v) if v.isdigit() else v.lower()
                 for v in re.split(r'(\d+)', str(value)))


def load_inputs(out):
    cells = pd.read_csv(out / 'cell_annotations.csv.gz', index_col=0)
    summary = pd.read_csv(out / 'cluster_summary.tsv', sep='\t')
    required = {'age', 'sample', 'stage_subcluster_id', *XY,
                'stage_state_label', 'stage_support_status', 'cycling_state'}
    assert required.issubset(cells.columns), required - set(cells.columns)
    assert {'stage', 'stage_subcluster_id', 'nuclei', 'stage_state_label',
            'stage_support_status'}.issubset(summary.columns)
    assert cells.index.is_unique and cells.index.notna().all()
    assert set(cells.age) == set(STAGES)
    assert np.isfinite(cells[list(XY)].to_numpy()).all()
    assert cells[list(required - set(XY))].notna().all().all()
    assert summary.stage_subcluster_id.is_unique
    assert set(summary.stage_subcluster_id) == set(cells.stage_subcluster_id)
    for cid, group in cells.groupby('stage_subcluster_id', sort=False):
        row = summary.loc[summary.stage_subcluster_id.eq(cid)].iloc[0]
        assert len(group) == int(row.nuclei)
        assert group.age.eq(row.stage).all(), cid
        assert group.stage_state_label.eq(row.stage_state_label).all(), cid
        assert group.stage_support_status.eq(row.stage_support_status).all(), cid
    assert int(summary.nuclei.sum()) == len(cells)
    return cells, summary


def palette(ids):
    if len(ids) <= len(COLORS):
        return {cid: COLORS[i] for i, cid in enumerate(ids)}
    # Keep numeric labels when many groups exceed the handpicked palette.
    return {cid: plt.get_cmap('hsv')(i / len(ids)) for i, cid in enumerate(ids)}


def anchors(cells, ids):
    span = np.maximum(np.ptp(cells[list(XY)].to_numpy(), axis=0), 1e-9)
    result = {}
    for cid in ids:
        xy = cells.loc[cells.stage_subcluster_id.eq(cid), list(XY)].to_numpy()
        center = np.median(xy, axis=0)
        # Select an observed cell, avoiding a label anchor in an empty hole.
        result[cid] = xy[np.argmin(np.sum(((xy - center) / span) ** 2, axis=1))]
    return result


def base(ax, cells, title):
    xy = cells[list(XY)].to_numpy()
    span = np.maximum(np.ptp(xy, axis=0), 1e-9)
    ax.set(xlim=(xy[:, 0].min() - .09 * span[0], xy[:, 0].max() + .09 * span[0]),
           ylim=(xy[:, 1].min() - .09 * span[1], xy[:, 1].max() + .09 * span[1]),
           xlabel='Within-stage UMAP 1', ylabel='Within-stage UMAP 2', title=title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal', adjustable='box')
    ax.spines[['top', 'right']].set_visible(False)


def scatter(ax, cells, color_column, colors):
    shuffled = cells.sample(frac=1, random_state=42)
    mapped = shuffled[color_column].map(colors)
    assert mapped.notna().all()
    ax.scatter(shuffled[XY[0]], shuffled[XY[1]],
               c=mapped.tolist(), s=4.2, alpha=.82,
               linewidths=0, rasterized=True)


def label_ids(ax, points, fontsize=10):
    for cid, xy in points.items():
        ax.text(*xy, str(cid), ha='center', va='center', fontsize=fontsize,
                weight='bold', zorder=5,
                bbox={'facecolor': 'white', 'edgecolor': '#aaaaaa',
                      'linewidth': .4, 'alpha': .92, 'pad': 2})


def draw_stage(cells, summary, stage, out, pdf):
    s = cells.loc[cells.age.eq(stage)]
    ids = sorted(s.stage_subcluster_id.unique(), key=cluster_key)
    colors = palette(ids)
    points = anchors(s, ids)
    libraries = sorted(s['sample'].unique())
    assert len(libraries) == 2, (stage, libraries)
    sample_colors = dict(zip(libraries, ['#0072B2', '#D55E00']))
    table_height = max(2.7, .62 * (len(ids) + 1))
    fig = plt.figure(figsize=(22, 10 + table_height))
    grid = fig.add_gridspec(2, 2, height_ratios=(8, table_height),
                           left=.045, right=.985, top=.87, bottom=.08,
                           hspace=.19, wspace=.15)
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    scatter(axes[0], s, 'stage_subcluster_id', colors)
    base(axes[0], s, 'Within-stage ventricular subclusters')
    scatter(axes[1], s, 'sample', sample_colors)
    base(axes[1], s, 'Same UMAP colored by library')
    for ax in axes:
        label_ids(ax, points)
    fig.legend(handles=[Line2D([], [], marker='o', linestyle='', color=c,
                                  label=f'{lib} (n={s["sample"].eq(lib).sum():,})',
                                  markersize=9)
                           for lib, c in sample_colors.items()],
                   loc='upper right', bbox_to_anchor=(.985, .955),
                   frameon=False, fontsize=11)
    table_ax = fig.add_subplot(grid[1, :])
    table_ax.axis('off')
    rows = []
    ordered = summary.set_index('stage_subcluster_id').loc[ids]
    for cid, row in ordered.iterrows():
        g = s.loc[s.stage_subcluster_id.eq(cid)]
        counts = g['sample'].value_counts()
        rows.append([cid, f'{int(row.nuclei):,}',
                     textwrap.fill(str(row.stage_state_label), 48),
                     textwrap.fill(str(row.stage_support_status), 43),
                     ' / '.join(f'{int(counts.get(lib, 0)):,}' for lib in libraries)])
    table = table_ax.table(cellText=rows,
                          colLabels=['Subcluster', 'Nuclei', 'Proposed state',
                                     'Evidence / sensitivity status', 'Library 1 / 2 nuclei'],
                          cellLoc='left', colLoc='left', loc='center',
                          colWidths=[.10, .07, .35, .32, .16], bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor('#dddddd')
        cell.set_linewidth(.5)
        if r == 0:
            cell.set_facecolor('#eeeeee')
            cell.get_text().set_weight('bold')
        elif c == 0:
            cell.get_text().set_color(colors[ids[r - 1]])
            cell.get_text().set_weight('bold')
    fig.suptitle(f'{stage} ventricular-like cardiomyocytes · {len(s):,} nuclei',
                 fontsize=23, weight='bold', y=.975)
    fig.text(.5, .035,
             'Identical coordinates and axis limits in both panels · Numeric groups are computational subclusters\n'
             'Proposed states remain RNA-based interpretations; evidence and sensitivity status are listed above.\n'
             'Library 1 / 2 order: ' + ' / '.join(libraries),
             ha='center', fontsize=11)
    stem = out / f'{stage}_ventricular_subclusters'
    fig.savefig(stem.with_suffix('.png'), dpi=180, facecolor='white')
    fig.savefig(stem.with_suffix('.pdf'), facecolor='white')
    pdf.savefig(fig)
    plt.close(fig)


def draw_overview(cells, summary, out, pdf):
    max_groups = int(summary.groupby('stage').size().max())
    fig = plt.figure(figsize=(28, max(21, 14 + .55 * max_groups)))
    outer = fig.add_gridspec(2, 2, left=.035, right=.99, bottom=.07,
                            top=.92, wspace=.10, hspace=.15)
    for slot, stage in zip(outer, STAGES):
        inner = GridSpecFromSubplotSpec(1, 2, subplot_spec=slot,
                                       width_ratios=[1.65, 1], wspace=.06)
        ax, info = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        s = cells.loc[cells.age.eq(stage)]
        ids = sorted(s.stage_subcluster_id.unique(), key=cluster_key)
        colors = palette(ids)
        scatter(ax, s, 'stage_subcluster_id', colors)
        base(ax, s, f'{stage} · {len(s):,} nuclei · {len(ids)} groups')
        label_ids(ax, anchors(s, ids), fontsize=8.5)
        info.axis('off')
        ordered = summary.set_index('stage_subcluster_id').loc[ids]
        y = .98
        gap = .93 / len(ids)
        for cid, row in ordered.iterrows():
            state = textwrap.fill(str(row.stage_state_label), 30)
            support = textwrap.fill(str(row.stage_support_status), 30)
            info.text(0, y, f'{cid} · n={int(row.nuclei):,}\n{state}\n[{support}]',
                      color=colors[cid], ha='left', va='top',
                      fontsize=min(10.5, 80 / max(1, len(ids))),
                      weight='bold', linespacing=1.15,
                      transform=info.transAxes)
            y -= gap
    fig.suptitle('Ventricular states resolved independently within each stage',
                 fontsize=26, weight='bold', y=.975)
    fig.text(.5, .945, 'Color and IDs identify within-stage computational groups; proposed labels and support are shown beside each panel.',
             ha='center', fontsize=14)
    fig.text(.5, .025,
             'P0, P7, P14 and P21 have independently fitted coordinate spaces; distances and orientation are not comparable across panels.\n'
             'These groups are proposed transcriptional states, not established anatomical or functional cell types.\n'
             'See each stage figure for matching library-colored UMAPs and per-library counts.',
             ha='center', fontsize=13)
    fig.savefig(out / 'ventricular_within_stage_overview.png', dpi=160, facecolor='white')
    fig.savefig(out / 'ventricular_within_stage_overview.pdf', facecolor='white')
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path,
                        default=ROOT / 'results/ventricular_within_stage')
    args = parser.parse_args()
    out = args.outdir
    cells, summary = load_inputs(out)
    with PdfPages(out / 'ventricular_within_stage_figures.pdf') as pdf:
        draw_overview(cells, summary, out, pdf)
        for stage in STAGES:
            draw_stage(cells, summary, stage, out, pdf)
    print(f'Saved overview and four stage comparisons for {len(cells):,} nuclei; coordinates and annotations unchanged.')


if __name__ == '__main__':
    main()
