"""Discover fibroblast states within each age, then test matrix and stage effects."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/fibroblast-state-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/fibroblast-state-numba')
import argparse
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import harmonypy
from anndata.io import read_elem, sparse_dataset
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
from fibroblast_state_common import *


def label_clusters(genes, zlog, x, cells, labels, global_mu, global_sd, stage):
    rows, markers, evidence = [], [], []
    family = np.array([UNRESOLVED] * len(cells), dtype=object)
    gene_lookup = {g: i for i, g in enumerate(genes)}
    for cl in sorted(set(labels), key=int):
        mask = labels == cl
        n = int(mask.sum())
        mean = np.asarray(zlog[mask].mean(axis=0)).ravel()
        frac = np.asarray((x[mask] >= 1).mean(axis=0)).ravel()
        diff = mean - np.asarray(zlog[~mask].mean(axis=0)).ravel()
        fracrest = np.asarray((x[~mask] >= 1).mean(axis=0)).ravel()
        eligible_libraries = []
        library_diffs = []
        for lib in sorted(cells['sample'].unique()):
            a = mask & cells['sample'].eq(lib).to_numpy()
            b = ~mask & cells['sample'].eq(lib).to_numpy()
            if min(a.sum(), b.sum()) >= 20:
                eligible_libraries.append(lib)
                library_diffs.append(np.asarray(zlog[a].mean(axis=0)).ravel() - np.asarray(zlog[b].mean(axis=0)).ravel())
        positive = np.sum(np.asarray(library_diffs) > .05, axis=0) if library_diffs else np.zeros(len(genes))
        rank = np.flatnonzero((diff > .15) & (frac >= .15) & ((frac-fracrest) >= .05) & ~np.isin(genes, SEX + CM_CHECK))
        rank = rank[np.argsort(diff[rank])[::-1]][:60]
        for j in rank:
            markers.append({'stage': stage, 'cluster_id': stage+'_F'+cl, 'gene': genes[j], 'log_expression_difference': diff[j],
                'fraction_ge1': frac[j], 'rest_fraction_ge1': fracrest[j], 'positive_libraries': int(positive[j])})
        candidates = []
        for fam, supported, anchors in choose_markers(genes, x, mean, frac):
            js = [gene_lookup[g] for g in PANELS[fam] if g in gene_lookup]
            score = float(((mean[js] - global_mu[js])/global_sd[js]).mean())
            # Nuclear RNA has the full panel; use >=3 markers and an anchor.
            qualified = len(supported) >= 3 and len(anchors) >= 1
            candidates.append((fam, score, qualified, supported, anchors))
        candidates.sort(key=lambda a: a[1], reverse=True)
        qualified = [a for a in candidates if a[2]]
        winner = qualified[0] if qualified else None
        margin = winner[1] - max([a[1] for a in candidates if a[0] != winner[0]], default=winner[1]) if winner else np.nan
        chosen = winner[0] if winner and winner[1] > .05 and margin >= .1 else UNRESOLVED
        family[mask] = chosen
        nc = cells.loc[mask, 'sample'].value_counts()
        both_markers = [genes[j] for j in rank if positive[j] >= 2]
        rows.append({'stage': stage, 'cluster_id': stage+'_F'+cl, 'nuclei': n, 'family': chosen,
            'best_qualified_family': winner[0] if winner else UNRESOLVED, 'family_score_margin': margin,
            'libraries_ge20': int((nc >= 20).sum()), 'largest_library_fraction': float(nc.max()/n),
            'two_library_positive_markers': len(both_markers), 'top_markers': ';'.join(both_markers[:12]) or 'None',
            'median_ambient_fraction': float(cells.loc[mask, 'decontx_contamination'].median()),
            'cycling_fraction': float(cells.loc[mask, 'cycling_marker_count'].ge(2).mean())})
        for fam, score, qualified, supported, anchors in candidates:
            evidence.append({'stage': stage, 'cluster_id': stage+'_F'+cl, 'candidate_family': fam,
                'score': score, 'panel_qualified': qualified, 'supported_genes': ';'.join(supported) or 'None',
                'anchors': ';'.join(anchors) or 'None'})
    return family, rows, markers, evidence


def regress_covariates(pc, obs):
    cov = np.column_stack([np.log1p(obs.total_counts), obs.decontx_contamination, obs.pct_counts_mt,
                           obs.cycling_marker_count])
    cov = (cov - cov.mean(0)) / np.maximum(cov.std(0), 1e-8)
    return pc - cov @ np.linalg.lstsq(cov, pc-pc.mean(0), rcond=None)[0]


def main():
    OUT.mkdir(exist_ok=True)
    with h5py.File(BASE / 'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        ii = np.flatnonzero(obs.reviewed_identity.eq('Fibroblast'))
        obs = obs.iloc[ii].copy()
        x = sparse_dataset(f['X'])[ii].tocsr().astype(np.float32)
        raw = sparse_dataset(f['layers/counts_raw'])[ii].tocsr().astype(np.float32)
    symbols = var.gene_symbol.astype(str).to_numpy()
    genes = pd.unique(symbols)
    if len(genes) != len(symbols):
        gi = pd.Index(genes).get_indexer(symbols)
        collapse = sparse.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)),gi)),shape=(len(symbols),len(genes)))
        x = (x@collapse).tocsr().astype(np.float32)
        raw = (raw@collapse).tocsr().astype(np.float32)
        var = var.groupby('gene_symbol',observed=True)[['mt','ribo']].any().reindex(genes)
        var['gene_symbol'] = genes
    assert len(obs) == 14985 and len(set(genes)) == len(genes)
    w = ad.AnnData(x.copy(), obs=obs.copy(), var=var.copy())
    sc.pp.normalize_total(w, target_sum=1e4)
    sc.pp.log1p(w)
    zlog = w.X.tocsr()
    # Equal weight to each library prevents a large library setting the panel calibration.
    mu = np.mean([np.asarray(zlog[obs['sample'].eq(s).to_numpy()].mean(axis=0)).ravel() for s in obs['sample'].unique()], axis=0)
    sd = np.maximum(np.sqrt(np.asarray(zlog.power(2).mean(axis=0)).ravel()-np.asarray(zlog.mean(axis=0)).ravel()**2), .25)
    annotations, summaries, marker_rows, panel_rows, mixing, sensitivity, embgenes = [], [], [], [], [], [], []
    for stage in STAGES:
        print('RNA_STAGE', stage, flush=True)
        folder = OUT / stage
        folder.mkdir(exist_ok=True)
        ix = np.flatnonzero(obs.age.eq(stage))
        o = obs.iloc[ix].copy()
        for col in o.select_dtypes(include='category').columns:
            o[col] = o[col].cat.remove_unused_categories()
        work = ad.AnnData(zlog[ix].copy(), obs=o[['sample','age']].copy(), var=var.copy())
        sc.pp.highly_variable_genes(work, n_top_genes=2000, flavor='seurat', batch_key='sample')
        eligible = ~var.mt.to_numpy() & ~var.ribo.to_numpy() & ~np.isin(genes, SEX + CM_CHECK)
        hv = work.var.highly_variable.to_numpy() & eligible
        for j in np.flatnonzero(hv):embgenes.append({'stage':stage,'gene':genes[j]})
        a = work[:,hv].copy()
        sc.pp.scale(a, zero_center=False, max_value=10)
        sc.pp.pca(a, n_comps=30, random_state=42, svd_solver='arpack')
        pc = a.obsm['X_pca'].copy()
        h = harmonypy.run_harmony(pc, o[['sample']], ['sample'], theta=2, lamb=1, nclust=20,
            max_iter_harmony=15, max_iter_kmeans=10, random_state=42, ncores=2, verbose=False)
        adjusted = np.asarray(h.Z_corr)
        if adjusted.shape != pc.shape:adjusted = adjusted.T
        adjusted += pc.mean(0)-adjusted.mean(0)
        q = ad.AnnData(sparse.csr_matrix((len(o),0)), obs=o.copy())
        reps = {'primary': adjusted, 'unadjusted': pc, 'QC_sensitivity': regress_covariates(adjusted,o)}
        # Marker/program withholding tests whether the whole-transcriptome grouping survives explicit panels.
        excluded = ANNOTATION_GENES | set(sum(MATRIX_PANELS.values(), [])) | set(WITHHOLD)
        hvm = hv & ~np.isin(genes, list(excluded))
        b = work[:,hvm].copy()
        sc.pp.scale(b, zero_center=False, max_value=10)
        sc.pp.pca(b, n_comps=30, random_state=42, svd_solver='arpack')
        hb = harmonypy.run_harmony(b.obsm['X_pca'], o[['sample']], ['sample'], theta=2, lamb=1,
            nclust=20, max_iter_harmony=15, max_iter_kmeans=10, random_state=42, ncores=2, verbose=False)
        pb = np.asarray(hb.Z_corr)
        reps['panel_withheld'] = pb if pb.shape == pc.shape else pb.T
        for name, rep in reps.items():
            q.obsm['X_'+name] = rep
            sc.pp.neighbors(q, n_neighbors=20, use_rep='X_'+name, key_added=name, random_state=42)
            sc.tl.leiden(q, neighbors_key=name, resolution=.45, flavor='igraph', n_iterations=3,
                         random_state=42, key_added=name+'_cluster')
            nn = NearestNeighbors(n_neighbors=21).fit(rep).kneighbors(rep, return_distance=False)[:,1:]
            libs = o['sample'].astype(str).to_numpy()
            expected = 1-np.sum((o['sample'].value_counts()/len(o))**2)
            mixing.append({'stage':stage,'representation':name,'cross_library_neighbor_fraction':float(np.mean(libs[:,None]!=libs[nn])),
                           'random_expected':float(expected)})
        for seed in [7,123]:
            sc.tl.leiden(q, neighbors_key='primary', resolution=.45, flavor='igraph', n_iterations=3,
                         random_state=seed, key_added='seed_'+str(seed))
        for res in [.25,.7]:
            sc.tl.leiden(q, neighbors_key='primary', resolution=res, flavor='igraph', n_iterations=3,
                         random_state=42, key_added='resolution_'+str(res))
        sc.tl.umap(q, neighbors_key='primary', min_dist=.4, maxiter=250, random_state=42)
        labels = q.obs.primary_cluster.astype(str).to_numpy()
        fam, rows, markers, panels = label_clusters(genes,zlog[ix],x[ix],o,labels,mu,sd,stage)
        q.obs['fibroblast_cluster_id'] = stage+'_F'+q.obs.primary_cluster.astype(str)
        q.obs['fibroblast_family'] = fam
        for row in rows:
            mask = q.obs.fibroblast_cluster_id.eq(row['cluster_id']).to_numpy()
            checks = []
            for key in ['unadjusted_cluster','QC_sensitivity_cluster','panel_withheld_cluster','seed_7','seed_123','resolution_0.25','resolution_0.7']:
                other = q.obs[key].astype(str).to_numpy()
                jac = max(np.sum(mask&(other==c))/np.sum(mask|(other==c)) for c in set(other))
                row[key+'_best_jaccard'] = float(jac)
                if key in ['seed_7','seed_123','QC_sensitivity_cluster','panel_withheld_cluster']:checks.append(jac)
            ok = row['nuclei']>=80 and row['libraries_ge20']==2 and row['largest_library_fraction']<=.9 and row['two_library_positive_markers']>=3
            row['support'] = 'supported RNA grouping' if ok and min(checks)>=.5 else 'candidate / sensitivity-dependent'
            if row['family'] == UNRESOLVED:row['support'] = 'unresolved family; '+row['support']
            q.obs.loc[mask,'fibroblast_support'] = row['support']
        for key in ['unadjusted_cluster','QC_sensitivity_cluster','panel_withheld_cluster','seed_7','seed_123','resolution_0.25','resolution_0.7']:
            sensitivity.append({'stage':stage,'comparison':key,'clusters':q.obs[key].nunique(),
                                'adjusted_rand_index':adjusted_rand_score(labels,q.obs[key].astype(str))})
        q.obs['fibroblast_UMAP1'],q.obs['fibroblast_UMAP2'] = q.obsm['X_umap'][:,0],q.obsm['X_umap'][:,1]
        q.obs.to_csv(folder/'RNA_cell_annotations.tsv.gz', sep='\t', na_rep='NA')
        q.write_h5ad(folder/'RNA_embeddings_metadata_only.h5ad', compression='gzip')
        annotations.append(q.obs)
        summaries.extend(rows);marker_rows.extend(markers);panel_rows.extend(panels)
    ann = pd.concat(annotations).loc[obs.index]
    ann.to_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',na_rep='NA')
    save_tsv(pd.DataFrame(summaries),OUT/'RNA_cluster_summary.tsv')
    save_tsv(pd.DataFrame(marker_rows),OUT/'RNA_within_stage_markers.tsv')
    save_tsv(pd.DataFrame(panel_rows),OUT/'RNA_family_marker_evidence.tsv')
    save_tsv(pd.DataFrame(mixing),OUT/'RNA_library_mixing.tsv')
    save_tsv(pd.DataFrame(sensitivity),OUT/'RNA_clustering_sensitivity.tsv')
    save_tsv(pd.DataFrame(embgenes),OUT/'RNA_embedding_genes.tsv')
    # Matrix outcomes deliberately omit all supplied state/cycling marker genes.
    scores = pd.DataFrame(index=ann.index)
    score_rows = []
    for name, panel in MATRIX_PANELS.items():
        for version in ['full','annotation_markers_withheld']:
            gl = [g for g in panel if g in genes and (version=='full' or g not in ANNOTATION_GENES)]
            js = np.flatnonzero(np.isin(genes,gl))
            if not len(js):continue
            scores[name+'__'+version] = np.asarray(zlog[:,js].mean(axis=1)).ravel()
            score_rows.append({'program':name,'version':version,'genes':';'.join(gl),'gene_count':len(gl)})
    # Independent response-associated genes and cell-cycle/foreign-lineage diagnostics.
    extra = {'TGFbeta_feedback':'Smad7 Skil Pmepa1 Bambi Serpine1 Junb'.split(), 'CM_carryover':CM_CHECK,
             'mural_check':MURAL_CHECK,'fibroblast_check':FB_CHECK}
    for name, gl in extra.items():
        js = np.flatnonzero(np.isin(genes,gl));scores[name] = np.asarray(zlog[:,js].mean(axis=1)).ravel()
    scores.to_csv(OUT/'RNA_cell_program_scores.tsv.gz',sep='\t')
    save_tsv(pd.DataFrame(score_rows),OUT/'RNA_matrix_program_definitions.tsv')
    combined = ann[['sample','age','fibroblast_family','fibroblast_cluster_id','total_counts','decontx_contamination']].join(scores)
    grouped = combined.groupby(['sample','age','fibroblast_family'],observed=True)
    means = grouped[list(scores)].mean().reset_index()
    sizes = grouped.size().rename('nuclei').reset_index()
    means = means.merge(sizes,on=['sample','age','fibroblast_family'])
    means = means.merge(grouped.agg(median_UMI=('total_counts','median'),mean_ambient=('decontx_contamination','mean')).reset_index(),on=['sample','age','fibroblast_family'])
    save_tsv(means,OUT/'RNA_matrix_programs_by_library_family.tsv')
    # Cell type x library pseudobulks preserve replicate labels. Original and corrected values are both retained.
    meta = ann[['sample','age','fibroblast_family','fibroblast_cluster_id']].copy()
    meta['group'] = meta['sample'].astype(str)+'::'+meta.fibroblast_cluster_id.astype(str)
    cats = pd.Categorical(meta.group)
    design = sparse.csr_matrix((np.ones(len(meta)),(cats.codes,np.arange(len(meta)))),shape=(len(cats.categories),len(meta)))
    for layer, counts in [('corrected',x),('raw',raw)]:
        pb = pd.DataFrame((design@counts).toarray().T,index=genes,columns=cats.categories)
        pb.index.name='gene';pb.to_csv(OUT/('RNA_'+layer+'_pseudobulk_counts.tsv.gz'),sep='\t',float_format='%.7g')
    gm = meta.groupby('group',observed=True).first().loc[cats.categories].copy()
    gm['nuclei'] = meta.groupby('group',observed=True).size().reindex(cats.categories)
    gm.to_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t')
    write_json(OUT/'RNA_methods.json',{'nuclei':len(ann),'libraries':8,'stage_counts':ann.age.value_counts().to_dict(),
        'HVG':2000,'neighbors':20,'resolution':.45,'seed':42,'stage_independent_clustering':True,
        'batch':'Harmony between libraries within each age only; no age regression',
        'annotation':'Post hoc exact paper panels; >=3 detected panel genes, >=1 anchor, score>0.05,margin>=0.1; otherwise unresolved',
        'state_name_limit':'Kirrel3-associated is not anatomically validated valve identity; labels are provisional RNA programs',
        'matrix_outcomes':'Mean log1p CP10k; full and all annotation-marker-withheld panels; not protein or mechanics',
        'support':'n>=80,both libraries>=20,largest library<=90%,>=3 positive markers in both; seed/QC/panel-withheld best Jaccard>=0.5',
        'unchanged_original_counts':True})
    print(pd.DataFrame(summaries)[['stage','cluster_id','nuclei','family','support','top_markers']].to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
