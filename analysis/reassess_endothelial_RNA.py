"""Endothelial-specific QC, stage-resolved discovery, and sensitivity analyses."""
import argparse
import h5py
import anndata as ad
import scanpy as sc
import harmonypy
from anndata.io import read_elem, sparse_dataset
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from endothelial_common import *

def collapse_read():
    with h5py.File(RNA_INPUT) as f:
        fullobs, var = read_elem(f['obs']), read_elem(f['var'])
        ii = np.flatnonzero(fullobs.reviewed_identity.eq('Endothelial'))
        obs = fullobs.iloc[ii].copy()
        x = sparse_dataset(f['X'])[ii].tocsr().astype(np.float32)
        raw = sparse_dataset(f['layers/counts_raw'])[ii].tocsr().astype(np.float32)
    symbols = var.gene_symbol.to_numpy(dtype=str)
    genes = pd.unique(symbols)
    c = sparse.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), pd.Index(genes).get_indexer(symbols))), shape=(len(symbols), len(genes)))
    x, raw = (x @ c).astype(np.float32).tocsr(), (raw @ c).astype(np.float32).tocsr()
    var = var.groupby('gene_symbol', observed=True)[['mt', 'ribo']].any().reindex(genes)
    var.index.name = 'gene'
    np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(), obs.total_counts, rtol=2e-5)
    np.testing.assert_allclose(np.asarray(raw.sum(axis=1)).ravel(), obs.raw_total_counts, rtol=2e-5)
    assert obs.scdblfinder_class.eq('singlet').all() and obs.pass_final_qc.all()
    assert obs.n_genes_by_counts.between(750, 7500).all() and obs.total_counts.ge(1000).all()
    assert obs.pct_counts_mt.le(5).all() and obs.decontx_contamination.le(.5).all()
    return obs, var, genes, x, raw

def harmony(pc, obs):
    if obs['sample'].nunique() < 2: return pc.copy()
    h = harmonypy.run_harmony(pc, obs[['sample']], ['sample'], theta=2, lamb=1, nclust=min(20, max(2, len(obs)//40)),
        max_iter_harmony=15, max_iter_kmeans=10, random_state=42, ncores=2, verbose=False)
    out = np.asarray(h.Z_corr)
    if out.shape != pc.shape: out = out.T
    return out + pc.mean(axis=0) - out.mean(axis=0)

def pca(z, obs, var, genes, extra_exclude=()):
    a = ad.AnnData(z.copy(), obs=obs[['sample', 'age']].copy(), var=var.copy())
    for c in a.obs.select_dtypes(include='category').columns:
        a.obs[c] = a.obs[c].cat.remove_unused_categories()
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor='seurat', batch_key='sample')
    use = a.var.highly_variable.to_numpy() & ~var.mt.to_numpy() & ~var.ribo.to_numpy() & ~np.isin(genes, WITHHOLD + list(extra_exclude))
    a = a[:, use].copy()
    sc.pp.scale(a, zero_center=False, max_value=10)
    sc.pp.pca(a, n_comps=min(30, len(obs)-2, use.sum()-2), svd_solver='arpack', random_state=42)
    return a.obsm['X_pca'].copy(), use

def cluster(q, pc, key, resolution=.45, seed=42):
    q.obsm['X_' + key] = pc
    sc.pp.neighbors(q, n_neighbors=min(20, len(q)-1), use_rep='X_' + key, key_added=key, random_state=seed)
    sc.tl.leiden(q, neighbors_key=key, resolution=resolution, random_state=seed, flavor='igraph', n_iterations=3, key_added=key+'_cluster')

def annotate(stage, o, labels, z, x, genes):
    # Stage-specific calibration, with equal library weighting in the baseline.
    mu = np.mean([np.asarray(z[o['sample'].eq(s).to_numpy()].mean(axis=0)).ravel() for s in o['sample'].unique()], axis=0)
    second = np.asarray(z.power(2).mean(axis=0)).ravel()
    sd = np.maximum(np.sqrt(np.maximum(second-np.asarray(z.mean(axis=0)).ravel()**2, 0)), .25)
    rows, markers, evidence = [], [], []
    fam = np.full(len(o), UNRESOLVED, object)
    for cl in sorted(set(labels), key=int):
        m = labels == cl
        mean = np.asarray(z[m].mean(axis=0)).ravel()
        frac = np.asarray((x[m] >= 1).mean(axis=0)).ravel()
        rest = np.asarray(z[~m].mean(axis=0)).ravel()
        restfrac = np.asarray((x[~m] >= 1).mean(axis=0)).ravel()
        diffs = []
        for lib in sorted(o['sample'].unique()):
            a, b = m & o['sample'].eq(lib).to_numpy(), ~m & o['sample'].eq(lib).to_numpy()
            if min(a.sum(), b.sum()) >= 10:
                diffs.append(np.asarray(z[a].mean(axis=0)).ravel()-np.asarray(z[b].mean(axis=0)).ravel())
        positives = np.sum(np.asarray(diffs) > .05, axis=0) if diffs else np.zeros(len(genes))
        technical = pd.Series(genes).str.lower().str.startswith('mt-').to_numpy() | pd.Series(genes).str.match(r'^Rp[sl]\d').to_numpy()
        rank = np.flatnonzero((mean-rest > .15) & (frac >= .15) & (frac-restfrac >= .05) & ~np.isin(genes, WITHHOLD) & ~technical)
        rank = rank[np.argsort((mean-rest)[rank])[::-1]][:60]
        cid = stage + '_EC' + cl
        for j in rank:
            markers.append({'stage': stage, 'cluster_id': cid, 'gene': genes[j], 'mean_log_difference': mean[j]-rest[j],
                'fraction_ge1': frac[j], 'rest_fraction_ge1': restfrac[j], 'positive_libraries': positives[j]})
        candidates = []
        for f in FAMILIES:
            ix = np.flatnonzero(np.isin(genes, PANELS[f]))
            detected = genes[ix[(frac[ix] >= .15) & (mean[ix] >= .1)]]
            anchors = sorted(set(detected) & set(ANCHORS[f]))
            required = 2 if f in FAMILIES[1:] else 1
            qualified = len(detected) >= 3 and len(anchors) >= required
            if f == FAMILIES[4]: qualified &= {'Prox1','Flt4'}.issubset(detected)
            score = float(np.mean((mean[ix]-mu[ix])/sd[ix]))
            candidates.append((f, score, qualified, detected.tolist(), anchors))
        candidates.sort(key=lambda v: v[1], reverse=True)
        qual = [a for a in candidates if a[2]]
        win = qual[0] if qual else None
        margin = win[1]-max(v[1] for v in candidates if v[0] != win[0]) if win else np.nan
        name = win[0] if win and win[1] > .05 and margin >= .10 else UNRESOLVED
        fam[m] = name
        top = [genes[j] for j in rank if positives[j] == 2]
        if not top: top = genes[rank[:3]].tolist()
        alias = '/'.join(top[:2]) + ' ECs' if top else 'Marker-poor ECs'
        nc = o.loc[m, 'sample'].value_counts()
        rows.append({'stage':stage, 'cluster_id':cid, 'gene_name':alias, 'family':name, 'nuclei':int(m.sum()),
            'best_qualified_family':win[0] if win else UNRESOLVED, 'family_score_margin':margin,
            'libraries_ge10':int((nc >= 10).sum()), 'largest_library_fraction':float(nc.max()/m.sum()),
            'two_library_markers':int(np.sum(positives[rank] == 2)), 'top_markers':';'.join(top[:12]),
            'median_ambient_fraction':float(o.loc[m, 'decontx_contamination'].median()),
            'cycling_fraction':float(o.loc[m, 'cycling_marker_count'].ge(2).mean())})
        for f, score, qualified, detected, anchors in candidates:
            evidence.append({'stage':stage,'cluster_id':cid,'family':f,'score':score,'qualified':qualified,
                'detected_genes':';'.join(detected), 'anchors':';'.join(anchors)})
    return fam, rows, markers, evidence

def coherent_mural_review(obs,var,genes,x,raw):
    """Two-pass conservative group review; never call every member a doublet."""
    excluded=pd.Series(False,index=obs.index);assignments=[];rows=[]
    for stage in STAGES:
        print('PROVISIONAL_LINEAGE_QC',stage,flush=True)
        ii=np.flatnonzero(obs.age.eq(stage) & obs.primary_EC_eligible)
        o=obs.iloc[ii].copy();pc,_=pca(lognorm(x[ii]),o,var,genes)
        q=ad.AnnData(sparse.csr_matrix((len(o),0)),obs=o[['sample','age']].copy())
        cluster(q,harmony(pc,o),'provisional_QC')
        labels=q.obs.provisional_QC_cluster.to_numpy(dtype=str)
        for cl in sorted(set(labels),key=int):
            m=labels==cl;libs=o.loc[m,'sample'].value_counts();libs=libs[libs>0];shared=None
            for lib,n in libs.items():
                if n<10:continue
                ix=ii[m & o['sample'].eq(lib).to_numpy()]
                j=np.flatnonzero(np.isin(genes,LINEAGES['Mural']))
                fc=np.asarray((x[ix][:,j]>=1).mean(axis=0)).ravel()
                fr=np.asarray((raw[ix][:,j]>=1).mean(axis=0)).ravel()
                present=set(genes[j[(fc>=.15)&(fr>=.15)]])
                shared=present if shared is None else shared & present
            flagged=(len(libs)==2 and (libs>=10).all() and shared is not None and len(shared)>=3)
            if flagged:excluded.loc[o.index[m]]=True
            cid=stage+'_ECQC'+cl
            rows.append({'stage':stage,'provisional_QC_cluster':cid,'nuclei':int(m.sum()),'coherent_mural_group_review':bool(flagged),
                'mural_markers_ge15pct_both_libraries_raw_and_corrected':';'.join(sorted(shared or [])),
                'minimum_nuclei_in_library':int(libs.min()),'interpretation':'Mixed EC/mural group review, not proof each nucleus is a doublet' if flagged else 'No coherent mural-group exclusion'})
            assignments.extend({'cell_id':c,'provisional_QC_cluster':cid,'coherent_mural_group_review':bool(flagged)} for c in o.index[m])
    save(pd.DataFrame(rows),OUT/'snRNA/provisional_mural_group_review.tsv')
    save(pd.DataFrame(assignments),OUT/'snRNA/provisional_QC_cluster_assignments.tsv.gz')
    return excluded

def main():
    folder = OUT / 'snRNA'; folder.mkdir(parents=True, exist_ok=True)
    obs, var, genes, x, raw = collapse_read()
    qc = lineage_evidence(x, genes, 'snRNA'); qc.index = obs.index
    obs = obs.join(qc)
    obs['primary_EC_eligible'] = obs.EC_lineage_supported & ~obs.mixed_lineage_review
    obs['QC_disposition'] = np.select([~obs.EC_lineage_supported, obs.mixed_lineage_review],
        ['Insufficient canonical EC support; review', 'Mixed-lineage program; excluded primary, retained review'], default='Primary endothelial reference')
    obs['ambient_gt30_sensitivity_flag'] = obs.decontx_contamination.gt(.30)
    obs['coherent_mural_group_review']=coherent_mural_review(obs,var,genes,x,raw)
    obs.loc[obs.coherent_mural_group_review,'primary_EC_eligible']=False
    obs.loc[obs.coherent_mural_group_review,'QC_disposition']='Coherent EC/mural program; excluded primary, retained review'
    obs.to_csv(folder/'all_endothelial_QC_decisions.tsv.gz', sep='\t', na_rep='NA')
    # Demonstrate what the inherited ambient correction changes; keep raw counts.
    profiles = []
    for lib, ix in obs.groupby('sample', observed=True).indices.items():
        for layer, counts in [('raw',raw),('DecontX_corrected',x)]:
            e = lineage_evidence(counts[ix], genes, 'snRNA')
            for panel in LINEAGES:
                profiles.append({'sample':lib,'stage':str(obs.iloc[ix[0]].age),'layer':layer,'panel':panel,'nuclei':len(ix),
                    'median_marker_fraction':e[panel+'_fraction'].median(),'mean_marker_fraction':e[panel+'_fraction'].mean(),
                    'strong_program_fraction':e['strong_'+panel+'_program'].mean() if panel != 'Endothelial' else e.EC_lineage_supported.mean()})
    save(pd.DataFrame(profiles), folder/'ambient_correction_lineage_comparison.tsv')
    saves, summaries, markers, panels, mixing, sens = [], [], [], [], [], []
    for stage in STAGES:
        print('ENDOTHELIAL_RNA_STAGE', stage, flush=True)
        ii = np.flatnonzero(obs.age.eq(stage) & obs.primary_EC_eligible)
        o = obs.iloc[ii].copy()
        o = o.drop(columns=[c for c in o if c.startswith('leiden_')])
        z = lognorm(x[ii]); rz = lognorm(raw[ii])
        pc, use = pca(z, o, var, genes)
        adj = harmony(pc, o)
        q = ad.AnnData(sparse.csr_matrix((len(o),0)), obs=o.copy())
        cov = np.column_stack([np.log1p(o.total_counts),o.decontx_contamination,o.pct_counts_mt,o.cycling_marker_count])
        cov = (cov-cov.mean(axis=0))/np.maximum(cov.std(axis=0),1e-8)
        residual = adj-cov @ np.linalg.lstsq(cov,adj-adj.mean(axis=0),rcond=None)[0]
        wp, _ = pca(z, o, var, genes, sum(PANELS.values(),[])+CYCLING)
        rp, _ = pca(rz, o, var, genes)
        for key, rep in [('primary',adj),('unadjusted',pc),('QC_adjusted',residual),('markers_withheld',harmony(wp,o)),('raw_counts',harmony(rp,o))]:
            cluster(q, rep, key)
            nn = NearestNeighbors(n_neighbors=21).fit(rep).kneighbors(rep,return_distance=False)[:,1:]
            libs=o['sample'].to_numpy(dtype=str)
            mixing.append({'stage':stage,'representation':key,'cross_library_neighbor_fraction':np.mean(libs[:,None] != libs[nn]),
                'random_expected':1-np.sum((o['sample'].value_counts()/len(o))**2)})
        for seed in [7,123]:
            sc.tl.leiden(q, neighbors_key='primary', resolution=.45, random_state=seed, flavor='igraph', n_iterations=3, key_added='seed_'+str(seed))
        for res in [.25,.7]:
            sc.tl.leiden(q, neighbors_key='primary', resolution=res, random_state=42, flavor='igraph', n_iterations=3, key_added='resolution_'+str(res))
        sc.tl.umap(q, neighbors_key='primary', random_state=42, min_dist=.4, maxiter=300)
        lab=q.obs.primary_cluster.to_numpy(dtype=str)
        family, rows, mr, ev=annotate(stage,o,lab,z,x[ii],genes)
        q.obs['EC_cluster_id']=stage+'_EC'+q.obs.primary_cluster.astype(str)
        q.obs['EC_family']=family
        for row in rows:
            m=q.obs.EC_cluster_id.eq(row['cluster_id']).to_numpy();checks=[]
            for key in ['unadjusted_cluster','QC_adjusted_cluster','markers_withheld_cluster','raw_counts_cluster','seed_7','seed_123','resolution_0.25','resolution_0.7']:
                other=q.obs[key].to_numpy(dtype=str)
                jac=max(np.sum(m & (other==c))/np.sum(m | (other==c)) for c in set(other))
                row[key+'_best_jaccard']=float(jac)
                if key in ['QC_adjusted_cluster','markers_withheld_cluster','seed_7','seed_123']:checks.append(jac)
            ok=row['nuclei'] >= 50 and row['libraries_ge10']==2 and row['largest_library_fraction'] <= .9 and row['two_library_markers'] >= 3 and min(checks)>=.5
            row['grouping_support']='Supported RNA grouping' if ok else 'Candidate / sensitivity-dependent grouping'
            mlow=m & o.decontx_contamination.le(.30).to_numpy()
            row['low_ambient_nuclei']=int(mlow.sum())
            q.obs.loc[m,'EC_gene_name']=row['gene_name']
            q.obs.loc[m,'EC_grouping_support']=row['grouping_support']
        for key in ['unadjusted_cluster','QC_adjusted_cluster','markers_withheld_cluster','raw_counts_cluster','seed_7','seed_123','resolution_0.25','resolution_0.7']:
            sens.append({'stage':stage,'comparison':key,'clusters':q.obs[key].nunique(),'adjusted_rand_index':adjusted_rand_score(lab,q.obs[key])})
        # Rebuild the graph after excluding >30% estimated ambient contamination.
        low=o.decontx_contamination.le(.30).to_numpy()
        if low.sum()>=50:
            lp,_=pca(lognorm(x[ii[low]]),o.loc[low],var,genes)
            lq=ad.AnnData(sparse.csr_matrix((low.sum(),0)),obs=o.loc[low,['sample','age']].copy())
            cluster(lq,harmony(lp,o.loc[low]),'low_ambient')
            q.obs['low_ambient_cluster']=pd.Series(lq.obs.low_ambient_cluster.astype(str),index=lq.obs_names).reindex(q.obs_names).fillna('Excluded >30%')
            for row in rows:
                m=q.obs.EC_cluster_id.eq(row['cluster_id']).to_numpy()[low];other=lq.obs.low_ambient_cluster.to_numpy(dtype=str)
                row['low_ambient_best_jaccard']=max(np.sum(m & (other==c))/max(1,np.sum(m | (other==c))) for c in set(other))
            sens.append({'stage':stage,'comparison':'low_ambient_recluster','clusters':lq.obs.low_ambient_cluster.nunique(),'adjusted_rand_index':adjusted_rand_score(lab[low],lq.obs.low_ambient_cluster)})
        q.obs['EC_UMAP1'],q.obs['EC_UMAP2']=q.obsm['X_umap'][:,0],q.obsm['X_umap'][:,1]
        dest=folder/stage;dest.mkdir(exist_ok=True)
        q.write_h5ad(dest/'EC_embeddings_metadata_only.h5ad',compression='gzip')
        q.obs.to_csv(dest/'EC_annotations.tsv.gz',sep='\t',na_rep='NA')
        counts=ad.AnnData(x[ii],obs=q.obs.copy(),var=var.copy());counts.layers['counts_raw']=raw[ii]
        counts.obsm['X_umap']=q.obsm['X_umap'];counts.uns['count_definition']='X: inherited DecontX-corrected fractional counts; counts_raw: original integer-valued counts; gene symbols collapsed by sum.'
        counts.write_h5ad(dest/'EC_corrected_and_raw_counts.h5ad',compression='gzip')
        saves.append(q.obs);summaries+=rows;markers+=mr;panels+=ev
    ann=pd.concat(saves)
    ann.to_csv(folder/'EC_annotations.tsv.gz',sep='\t',na_rep='NA')
    save(pd.DataFrame(summaries),folder/'EC_cluster_summary.tsv')
    save(pd.DataFrame(markers),folder/'EC_cluster_markers.tsv')
    save(pd.DataFrame(panels),folder/'EC_family_marker_evidence.tsv')
    save(pd.DataFrame(mixing),folder/'EC_library_mixing.tsv')
    save(pd.DataFrame(sens),folder/'EC_clustering_sensitivity.tsv')
    # RNA-library x cluster pseudobulks preserve technical/biological unit labels.
    indices=obs.index.get_indexer(ann.index)
    meta=ann[['age','sample','EC_cluster_id','EC_gene_name','EC_family']].copy()
    meta['group']=meta['sample'].astype(str)+'::'+meta.EC_cluster_id.astype(str)
    cat=pd.Categorical(meta.group)
    design=sparse.csr_matrix((np.ones(len(meta)),(cat.codes,np.arange(len(meta)))),shape=(len(cat.categories),len(meta)))
    for name, counts in [('corrected',x),('raw',raw)]:
        pb=pd.DataFrame((design@counts[indices]).toarray().T,index=genes,columns=cat.categories)
        pb.index.name='gene';pb.to_csv(folder/f'EC_{name}_pseudobulk_counts.tsv.gz',sep='\t',float_format='%.7g')
    pm=meta.groupby('group',observed=True).first().reindex(cat.categories)
    pm['nuclei']=meta.groupby('group',observed=True).size().reindex(cat.categories)
    pm.to_csv(folder/'EC_pseudobulk_metadata.tsv',sep='\t')
    save(obs.groupby(['age','sample','QC_disposition'],observed=True).size().reset_index(name='nuclei'),folder/'QC_counts_by_library.tsv')
    write_json(folder/'RNA_complete.json',{'source_EC_nuclei':len(obs),'primary_EC_nuclei':len(ann),'review_retained':int((~obs.primary_EC_eligible).sum()),
        'all_inherited_doublet_calls_singlet':True,'all_inherited_postambient_QC_pass':True,'counts_changed':False,
        'RNA_input_sha256':sha(RNA_INPUT),'ambient_method':'Inherited converged per-library DecontX, called nuclei only, no empty-droplet reference',
        'doublet_method':'Inherited per-library scDblFinder; residual incompatible-lineage program review; Scrublet retained diagnostic only',
        'two_pass_mural_review':'Provisional within-age clusters with >=3 mural diagnostic genes detected in>=15% of nuclei, both raw/corrected, in both libraries (>=10 nuclei/library) excluded as group review; final primary graph refitted. Not proof all group members are doublets.',
        'batch':'Harmony only between two libraries within each age; no age regression',
        'support':'n>=50, >=10 in both libraries, max-library fraction<=.9, >=3 markers positive in both, minimum seed/QC/marker-withheld cluster Jaccard>=.5',
        'mixed_program_limit':'Operational review excludes from primary reference, not proof of doublet or false biology; review nuclei and all original counts retained.'})
    print(pd.DataFrame(summaries)[['stage','gene_name','nuclei','family','grouping_support','top_markers']].to_string(index=False),flush=True)

if __name__=='__main__':main()
