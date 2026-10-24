"""Display measured RNA and candidate labels in native Xenium section coordinates."""
from pathlib import Path
import argparse
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/heart-mpl')
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, to_rgb
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/xenium_ventricular_mapping'
STAGES = ['P0', 'P7', 'P14', 'P21']
COLORS = {'Ankrd1-associated': '#d87523', 'Conduction-program candidate': '#8e44ad', 'Cycling-associated': '#35a853', 'Regional electrical candidate': '#1f9eab', 'DNA-maintenance-associated': '#c74d63', 'Other working CM': '#356bad', 'Unresolved ventricular state': '#b8b8b8', 'Not assigned ventricular CM': '#eeeeee'}
LINEAGE = {'Ventricular-enriched CM': '#356bad', 'Atrial-marker-enriched CM': '#e78335', 'CM; chamber unresolved': '#819bb1', 'CM/EC program overlap; review': '#9d9d9d', 'Endothelial': '#32a68a', 'Fibroblast': '#b35494', 'Pericyte': '#967348', 'Smooth muscle': '#d5ae27', 'Immune': '#cf4c4c', 'Epicardial': '#7d6bc1', 'Unresolved cell identity': '#dddddd', 'Low-information/QC flag': '#eeeeee'}
MARKER_PANELS = [('Tnnt2', ['Tnnt2']), ('Myl2', ['Myl2']), ('Myl7', ['Myl7']), ('Ankrd1', ['Ankrd1']), ('Slit2', ['Slit2']), ('Cntn2', ['Cntn2']), ('Gja5', ['Gja5']), ('Sema3a', ['Sema3a']), ('Ddc', ['Ddc']), ('Mki67 + Top2a', ['Mki67', 'Top2a']), ('Npr3', ['Npr3']), ('Upk3b + Wt1', ['Upk3b', 'Wt1'])]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})


def background(stage, obs):
    folder = OUT / stage
    image = np.asarray(Image.open(folder / 'DAPI_overview.png')) / 255
    info = json.loads((folder / 'morphology_display.json').read_text())
    x, y = obs.x_centroid.to_numpy(), obs.y_centroid.to_numpy()
    bounds = [max(0, x.min() - 100), min(info['extent_um'][1], x.max() + 100), max(0, y.min() - 100), min(info['extent_um'][2], y.max() + 100)]
    return 1 - .35 * image, info['extent_um'], bounds


def frame(ax, bg, extent, bounds, title):
    ax.imshow(bg, cmap='gray', vmin=0, vmax=1, extent=extent, origin='upper', interpolation='nearest', rasterized=True)
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[3], bounds[2])
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=13, weight='bold', pad=9)
    ax.set_xlabel('X (µm)')
    ax.set_ylabel('Y (µm)')
    ax.tick_params(labelsize=8)
    x = bounds[0] + .06 * (bounds[1] - bounds[0])
    y = bounds[3] - .06 * (bounds[3] - bounds[2])
    length = 500 if bounds[1] - bounds[0] > 1500 else 100
    ax.plot([x, x + length], [y, y], color='black', lw=3)
    ax.text(x + length / 2, y - .014 * (bounds[3] - bounds[2]), f'{length} µm', ha='center', va='bottom', fontsize=8)


def scatter_categories(ax, obs, field, colors, mask=None, point_size=.8, legend=True):
    mask = np.ones(len(obs), dtype=bool) if mask is None else np.asarray(mask)
    values = obs[field].astype(str)
    handles = []
    for value in sorted(values[mask].unique(), key=lambda value: (value.startswith('Unresolved'), value)):
        selected = mask & values.eq(value).to_numpy()
        color = colors.get(value, '#bbbbbb')
        ax.scatter(obs.loc[selected, 'x_centroid'], obs.loc[selected, 'y_centroid'], s=point_size, color=color, edgecolors='none', rasterized=True, alpha=.85)
        handles.append(Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markersize=5, label=f'{value} ({selected.sum():,})'))
    if legend:
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False, borderaxespad=0)


def save(fig, folder, name, pdf):
    fig.savefig(folder / (name + '.png'), dpi=180, facecolor='white')
    fig.savefig(folder / (name + '.pdf'), facecolor='white')
    pdf.savefig(fig, facecolor='white')
    plt.close(fig)


def fine_colors(catalog):
    result = {'Not resolved': '#b8b8b8'}
    for family, group in catalog.groupby('mapping_family', sort=False):
        base = np.array(to_rgb(COLORS.get(family, '#7d8c99')))
        for i, cluster in enumerate(group.microenvironment_cluster):
            mix = min(.40, .14 * i)
            result[cluster] = tuple(base * (1 - mix) + mix)
    return result


def maps(stage, data, catalog, pdf):
    obs = data.obs
    folder = OUT / stage
    bg, extent, bounds = background(stage, obs)
    fig, axes = plt.subplots(2, 2, figsize=(21, 17))
    frame(axes[0, 0], bg, extent, bounds, 'Measured-marker-supported lineage candidates')
    scatter_categories(axes[0, 0], obs, 'spatial_lineage_candidate', LINEAGE, point_size=.65)
    vent = obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
    frame(axes[0, 1], bg, extent, bounds, 'Ventricular program-family candidates')
    scatter_categories(axes[0, 1], obs, 'spatial_mapped_family', COLORS, mask=vent, point_size=1.0)
    frame(axes[1, 0], bg, extent, bounds, 'Fine reference groups passing transfer criteria')
    scatter_categories(axes[1, 0], obs, 'spatial_mapped_fine_group', fine_colors(catalog), mask=vent, point_size=1.0)
    frame(axes[1, 1], bg, extent, bounds, 'Measured multi-gene program coexpression')
    base = vent.to_numpy()
    axes[1, 1].scatter(obs.loc[base, 'x_centroid'], obs.loc[base, 'y_centroid'], s=.5, c='#cccccc', edgecolors='none', rasterized=True)
    handles = []
    for field, label, color in [('measured_Ankrd1_program_positive', 'Ankrd1-associated', COLORS['Ankrd1-associated']), ('measured_Slit2_program_positive', 'Slit2-associated', '#1f9eab'), ('measured_conduction_overlap_positive', 'Cntn2/Gja5 + conduction overlap', COLORS['Conduction-program candidate']), ('measured_cycling_program_positive', 'Cycling-marker coexpression', COLORS['Cycling-associated'])]:
        mask = obs[field].to_numpy(dtype=bool)
        axes[1, 1].scatter(obs.loc[mask, 'x_centroid'], obs.loc[mask, 'y_centroid'], s=1.4, c=color, edgecolors='none', rasterized=True)
        handles.append(Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markersize=5, label=f'{label} ({mask.sum():,})'))
    axes[1, 1].legend(handles=handles, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8, frameon=False)
    fig.suptitle(f'{stage}: ventricular candidates on the acquired cardiac section', fontsize=20, weight='bold')
    fig.text(.5, .025, 'Gray image: acquired DAPI; colored points: original cell centroids. Candidate labels are reference-assisted, not anatomical or functional verification.\n'
             'Fine assignments require panel/library/raw-reference agreement. Unresolved cells remain gray; measured program flags can overlap.\n'
             'Kcnd2 and Irx5 are NOT assayed in this Xenium panel. No Kcnd2 expression, imputation or transmural coordinate is shown.', ha='center', fontsize=10, linespacing=1.5)
    fig.subplots_adjust(left=.05, right=.83, bottom=.10, top=.94, wspace=.60, hspace=.25)
    save(fig, folder, f'{stage}_spatial_candidate_maps', pdf)


def marker_maps(stage, data, maxima, pdf):
    obs = data.obs
    bg, extent, bounds = background(stage, obs)
    fig, axes = plt.subplots(3, 4, figsize=(22, 18))
    for ax, (label, genes) in zip(axes.ravel(), MARKER_PANELS):
        frame(ax, bg, extent, bounds, label)
        values = np.asarray(data.X[:, data.var_names.get_indexer(genes)].sum(axis=1)).ravel()
        selected = obs.pass_expression_QC.to_numpy() & (values > 0)
        vmax = maxima[label]
        order = np.flatnonzero(selected)
        order = order[np.argsort(values[order])]
        points = ax.scatter(obs.iloc[order].x_centroid, obs.iloc[order].y_centroid, c=values[order], norm=LogNorm(vmin=1, vmax=vmax), cmap='magma_r', s=.7, edgecolors='none', rasterized=True)
        fig.colorbar(points, ax=ax, shrink=.66, pad=.02).set_label('Measured counts/cell', fontsize=8)
    fig.suptitle(f'{stage}: directly measured Xenium RNA on DAPI', fontsize=21, weight='bold')
    fig.text(.5, .025, 'All QC-pass segmented cells shown; these genes are not all CM-specific (notably Gja5). Same gene uses the same scale across stages.\n'
             'No count correction, smoothing or expression imputation. Marker counts depend on cell area, segmentation and detection efficiency.\n'
             'Npr3 and Upk3b/Wt1 offer candidate anatomical landmarks; they are not traced wall boundaries. Kcnd2/Irx5 are absent from the panel.', ha='center', fontsize=11, linespacing=1.5)
    fig.subplots_adjust(left=.04, right=.97, bottom=.10, top=.94, wspace=.32, hspace=.28)
    save(fig, OUT / stage, f'{stage}_measured_marker_maps', pdf)


def zoom(stage, data, catalog, pdf):
    obs = data.obs
    candidates = obs.spatial_mapped_family.eq('Conduction-program candidate').to_numpy()
    reason = 'density of reference-assisted conduction-program candidates'
    if candidates.sum() < 5:
        candidates = obs.measured_Slit2_program_positive.to_numpy(dtype=bool)
        reason = 'density of measured Slit2-program-positive candidates'
    if candidates.sum() < 5:
        candidates = obs.spatial_lineage_candidate.eq('Ventricular-enriched CM').to_numpy()
        reason = 'density of ventricular marker-supported candidates'
    coords = obs.loc[candidates, ['x_centroid', 'y_centroid']].to_numpy()
    if not len(coords):
        return
    density = cKDTree(coords).query_ball_point(coords, 180, return_length=True)
    center = coords[int(np.argmax(density))]
    bounds = [center[0] - 300, center[0] + 300, center[1] - 300, center[1] + 300]
    selected = obs.x_centroid.between(bounds[0], bounds[1]) & obs.y_centroid.between(bounds[2], bounds[3])
    ids = set(obs.index[selected])
    source = BASE / 'data/xenium' / str(obs.section.iloc[0]) / 'cell_boundaries.csv.gz'
    if source.exists():
        chunks = []
        for chunk in pd.read_csv(source, chunksize=250000):
            chunks.append(chunk[chunk.cell_id.isin(ids)])
        boundaries = pd.concat(chunks, ignore_index=True)
    else:
        published = json.loads((OUT / stage / 'zoom_selection.json').read_text())
        center = np.array(published['center_um'])
        bounds = published['bounds_um']
        reason = 'the published illustrative ROI (full boundary archive unavailable locally)'
        selected = obs.x_centroid.between(bounds[0], bounds[1]) & obs.y_centroid.between(bounds[2], bounds[3])
        boundaries = pd.read_csv(OUT / stage / 'zoom_original_cell_boundary_vertices.tsv.gz', sep='\t')
    boundaries.to_csv(OUT / stage / 'zoom_original_cell_boundary_vertices.tsv.gz', sep='\t', index=False)
    polygons, cell_ids = [], []
    for cell_id, group in boundaries.groupby('cell_id', sort=False):
        polygons.append(group[['vertex_x', 'vertex_y']].to_numpy())
        cell_ids.append(cell_id)
    positions = obs.index.get_indexer(cell_ids)
    assert (positions >= 0).all()
    bg, extent, _ = background(stage, obs)
    fig, axes = plt.subplots(1, 3, figsize=(20, 7.5))
    for ax, gene in zip(axes[:2], ['Slit2', 'Cntn2']):
        frame(ax, bg, extent, bounds, f'Measured {gene} + original boundaries')
        values = data.X[positions, data.var_names.get_loc(gene)].toarray().ravel()
        maximum = max(3, int(np.quantile(values[values > 0], .98))) if (values > 0).any() else 3
        colors = plt.get_cmap('magma_r')(LogNorm(vmin=1, vmax=maximum)(np.maximum(values, 1)))
        colors[values == 0] = [.96, .96, .96, .15]
        collection = PolyCollection(polygons, facecolors=colors, edgecolors='#999999', linewidths=.25, rasterized=True)
        ax.add_collection(collection)
        scalar = plt.cm.ScalarMappable(norm=LogNorm(vmin=1, vmax=maximum), cmap='magma_r')
        fig.colorbar(scalar, ax=ax, shrink=.65).set_label('Measured counts/cell')
    ax = axes[2]
    frame(ax, bg, extent, bounds, 'Reference-group candidates + original boundaries')
    palette = fine_colors(catalog)
    names = obs.iloc[positions].spatial_mapped_fine_group.astype(str)
    colors = [palette.get(value, '#dddddd') for value in names]
    ax.add_collection(PolyCollection(polygons, facecolors=colors, edgecolors='#999999', linewidths=.25, alpha=.8, rasterized=True))
    handles = [Line2D([0], [0], marker='s', color='none', markerfacecolor=palette.get(value, '#dddddd'), label=f'{value} ({(names == value).sum()})') for value in sorted(names.unique())]
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9, frameon=False)
    fig.suptitle(f'{stage}: 600 × 600 µm illustrative region', fontsize=18, weight='bold')
    fig.text(.5, .025, f'ROI selected by {reason}; not an unbiased enrichment test.\nOriginal boundaries and RNA counts retained; labels do not verify Purkinje identity. Kcnd2 is not measured.', ha='center', fontsize=10, linespacing=1.5)
    fig.subplots_adjust(left=.05, right=.88, top=.88, bottom=.17, wspace=.45)
    (OUT / stage / 'zoom_selection.json').write_text(json.dumps({'center_um': center.tolist(), 'bounds_um': bounds, 'selection': reason, 'window_um': 600, 'cells_with_centroid_in_window': int(selected.sum()), 'boundary_cells_plotted': len(cell_ids)}, indent=2))
    save(fig, OUT / stage, f'{stage}_boundary_zoom', pdf)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stages', nargs='+', default=STAGES)
    args = parser.parse_args()
    values = {label: [] for label, _ in MARKER_PANELS}
    for stage in STAGES:
        data = ad.read_h5ad(OUT / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        for label, genes in MARKER_PANELS:
            count = np.asarray(data.X[:, data.var_names.get_indexer(genes)].sum(axis=1)).ravel()
            values[label].append(count[count > 0])
    maxima = {label: max(3, float(np.quantile(np.concatenate(parts), .995))) for label, parts in values.items()}
    (OUT / 'measured_marker_color_scales.json').write_text(json.dumps(maxima, indent=2))
    catalog = pd.read_csv(OUT / 'reference_group_catalog.tsv', sep='\t')
    name = 'xenium_ventricular_spatial_figures.pdf' if args.stages == STAGES else 'xenium_' + '_'.join(args.stages) + '_preview.pdf'
    with PdfPages(OUT / name) as pdf:
        for stage in args.stages:
            data = ad.read_h5ad(OUT / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
            maps(stage, data, catalog[catalog.stage.eq(stage)], pdf)
            marker_maps(stage, data, maxima, pdf)
            zoom(stage, data, catalog[catalog.stage.eq(stage)], pdf)
            print('SPATIAL_FIGURES_READY', stage, flush=True)


if __name__ == '__main__':
    main()
