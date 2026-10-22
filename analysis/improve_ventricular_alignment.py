"""Compare fixed weighted projections with native CytoSPACE on withheld genes.

The benchmark conditions on the previously frozen ventricular/QC cohort.
Neither published spatial labels nor held-out expression enter the fit.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import precision_recall_fscore_support, roc_auc_score

from run_cytospace_xenium import ROOT, SOURCE, CYTO, BASE, STAGES, HOLDOUT, write_matrix, digest

OUT = ROOT / 'results/ventricular_alignment_improvement'
CACHE = BASE / 'results/ventricular_improvement_cache'
METHODS = ['native_CytoSPACE', 'Pearson_nearest1', 'Pearson_soft30',
           'assay_scaled_soft30', 'assay_scaled_spatial', 'assay_scaled_shuffled_spatial']


def normalized(counts, training):
    # Exclude every withheld gene even from the size-factor denominator.
    total = counts[:, training].sum(1).astype(np.float64)
    assert (total > 0).all()
    return np.log2(1 + 1e6 * counts / total[:, None]).astype(np.float32)


def unit_features(ref, query, scaled):
    r, q = ref.copy(), query.copy()
    if scaled:
        r = np.clip((r-r.mean(0))/np.maximum(r.std(0), .25), -5, 5)
        q = np.clip((q-q.mean(0))/np.maximum(q.std(0), .25), -5, 5)
    else:
        r -= r.mean(1, keepdims=True)
        q -= q.mean(1, keepdims=True)
    r /= np.maximum(np.linalg.norm(r, axis=1, keepdims=True), 1e-8)
    q /= np.maximum(np.linalg.norm(q, axis=1, keepdims=True), 1e-8)
    return r.astype(np.float32), q.astype(np.float32)


def neighbors(ref, query, k=30):
    indices, weights, similarity = [], [], []
    k = min(k, len(ref))
    for start in range(0, len(query), 1024):
        score = query[start:start+1024] @ ref.T
        ix = np.argpartition(score, -k, axis=1)[:, -k:]
        value = np.take_along_axis(score, ix, axis=1)
        order = np.argsort(-value, axis=1, kind='stable')
        ix, value = np.take_along_axis(ix, order, axis=1), np.take_along_axis(value, order, axis=1)
        w = np.exp((value-value[:, :1])/.10)
        w /= w.sum(1, keepdims=True)
        indices.append(ix.astype(np.int32)); weights.append(w); similarity.append(value[:, 0])
    return np.concatenate(indices), np.concatenate(weights), np.concatenate(similarity)


def project(values, ix, weights):
    result = np.empty((len(ix), values.shape[1]), dtype=np.float32)
    for start in range(0, len(ix), 2048):
        sl = slice(start, start+2048)
        result[sl] = np.einsum('nk,nkg->ng', weights[sl], values[ix[sl]])
    return result


def posteriors(labels, ix, weights):
    classes = sorted(set(labels))
    encoded = np.array([classes.index(x) for x in labels])
    result = np.zeros((len(ix), len(classes)), dtype=np.float32)
    for j in range(len(classes)):
        result[:, j] = (weights*(encoded[ix] == j)).sum(1)
    return classes, result


def graph(coords, features):
    d, ix = cKDTree(coords).query(coords, k=13)
    # Self is removed; no neighbors across distances >75 µm.
    d, ix = d[:, 1:], ix[:, 1:]
    sim = np.einsum('ng,nkg->nk', features, features[ix])
    w = np.exp(-d/30) * np.maximum(sim, 0)**2 * (d <= 75)
    w /= np.maximum(w.sum(1, keepdims=True), 1e-12)
    return sparse.csr_matrix((w.ravel(), (np.repeat(np.arange(len(ix)), ix.shape[1]), ix.ravel())), shape=(len(ix), len(ix)))


def safe_r(x, y, rank=False):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float((spearmanr if rank else pearsonr)(x, y).statistic)


def native(stage, q, r, genes, qids, rids, coords):
    dest = OUT / stage / 'native_heldout'
    dest.mkdir(parents=True, exist_ok=True)
    result = dest / 'assigned_locations.csv.gz'
    if not result.exists():
        work = CACHE / stage
        rng = np.random.default_rng(42)
        qi, ri = rng.permutation(len(q)), rng.permutation(len(r))
        write_matrix(work/'reference', r[ri], genes, rids[ri])
        write_matrix(work/'query', q[qi], genes, qids[qi])
        for name, ids in [('reference',rids[ri]),('query',qids[qi])]:
            pd.DataFrame({'CellType':'Ventricular_CM'},index=ids).to_csv(work/f'{name}_labels.tsv',sep='\t',index_label='CellID' if name=='reference' else 'SpotID')
        pd.DataFrame(coords[qi],index=qids[qi],columns=['x_centroid','y_centroid']).to_csv(work/'coordinates.tsv',sep='\t',index_label='SpotID')
        command = [str(CYTO),'--single-cell','--scRNA-path',str(work/'reference/matrix.mtx.gz'),
            '--cell-type-path',str(work/'reference_labels.tsv'),'--st-path',str(work/'query/matrix.mtx.gz'),
            '--st-cell-type-path',str(work/'query_labels.tsv'),'--coordinates-path',str(work/'coordinates.tsv'),
            '--output-folder',str(work/'output'),'--number-of-selected-spots','3000','--number-of-processors','2',
            '--seed','42','--downsample-off','--plot-off']
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
            'NUMBA_CACHE_DIR':'/tmp/cytospace-numba','MPLCONFIGDIR':'/tmp/heart-mpl'}
        with open(dest/'upstream_stdout.txt','w') as handle:
            subprocess.run(command,env=env,stdout=handle,stderr=subprocess.STDOUT,check=True)
        assigned=pd.read_csv(work/'output/assigned_locations.csv',float_precision='round_trip')
        assigned.to_csv(result,index=False)
        shutil.copyfile(work/'output/log.txt',dest/'upstream_log.txt')
        (dest/'run_manifest.json').write_text(json.dumps({'command':command,'query_cells':len(q),'reference_cells':len(r),
            'training_genes':list(genes),'reference_values':'DecontX corrected unnormalized counts',
            'single_cell':True,'broad_label_only':True,'upstream_source_modified':False,
            'upstream_commit':'bb7fe12ab4466ef5e38dcac4c3ece07012060d60'},indent=2))
        shutil.rmtree(work)
    assigned=pd.read_csv(result,float_precision='round_trip').set_index('SpotID')
    assert assigned.index.is_unique and set(assigned.index)==set(qids)
    np.testing.assert_array_equal(assigned.loc[qids,['x_centroid','y_centroid']],coords)
    index=pd.Index(rids).get_indexer(assigned.loc[qids,'OriginalCID'])
    assert (index>=0).all()
    return index


def metrics(stage, qobs, truth, predicted, tests, counts, normalization_depth, training_detection):
    records=[]
    coords=qobs[['x_centroid','y_centroid']].to_numpy()
    grid=pd.DataFrame(np.floor(coords/250).astype(int),columns=['gx','gy'])
    cov=np.column_stack([np.ones(len(qobs)),np.log1p(normalization_depth),np.log1p(qobs.cell_area),np.log1p(training_detection)])
    basis=np.linalg.qr(cov)[0]
    for mi,method in enumerate(METHODS):
        for j,gene in enumerate(tests):
            x,y=predicted[mi,:,j],truth[:,j]
            binned=grid.assign(x=x,y=y).groupby(['gx','gy']).agg(x=('x','mean'),y=('y','mean'),n=('x','size'))
            binned=binned[binned.n>=20]
            detected=counts[:,j]>0
            records.append({'stage':stage,'method':method,'gene':gene,'cells':len(y),
                'cell_spearman':safe_r(x,y,True),'cell_pearson':safe_r(x,y),
                'depth_area_training_detection_adjusted_r':safe_r(x-basis@(basis.T@x),y-basis@(basis.T@y)),
                'detection_AUC':roc_auc_score(detected,x) if len(np.unique(detected))==2 else np.nan,
                'grid250um_spearman':safe_r(binned.x,binned.y,True),'grid_bins':len(binned),
                'positive_fraction':float(detected.mean())})
    return records


def library_benchmark(stage, counts, meta):
    records=[]
    libraries=sorted(meta['sample'].unique())
    for library in libraries:
        qi=np.flatnonzero(meta['sample'].eq(library));ri=np.flatnonzero(~meta['sample'].eq(library))
        rn,qn=normalized(counts[ri],np.ones(counts.shape[1],bool)),normalized(counts[qi],np.ones(counts.shape[1],bool))
        for scaled in [False,True]:
            rf,qf=unit_features(rn,qn,scaled)
            ix,w,_=neighbors(rf,qf)
            for field in ['mapping_family','ventricular_group']:
                classes,p=posteriors(meta.iloc[ri][field].astype(str).to_numpy(),ix,w)
                yh=np.array(classes)[p.argmax(1)];y=meta.iloc[qi][field].astype(str).to_numpy()
                precision,recall,f1,n=precision_recall_fscore_support(y,yh,labels=classes,zero_division=0)
                for label,a,b,c,d in zip(classes,precision,recall,f1,n):
                    records.append({'stage':stage,'query_library':library,'method':'assay_scaled_soft30' if scaled else 'Pearson_soft30',
                        'resolution':field,'label':label,'precision':a,'recall':b,'f1':c,'query_nuclei':int(d),
                        'training_nuclei':int(meta.iloc[ri][field].eq(label).sum())})
    return records


def final_mapping(stage, query, reference, qobs, rmeta, allobs):
    train=np.ones(query.shape[1],bool)
    qn,rn=normalized(query,train),normalized(reference,train)
    rf,qf=unit_features(rn,qn,True)
    ix,w,sim=neighbors(rf,qf)
    fam,fp=posteriors(rmeta.mapping_family.astype(str).to_numpy(),ix,w)
    groups,gp=posteriors(rmeta.ventricular_group.astype(str).to_numpy(),ix,w)
    result=qobs[['x_centroid','y_centroid','cell_area','transcript_counts']].copy()
    result['weighted_family']=np.array(fam)[fp.argmax(1)]
    result['family_neighbor_weight_fraction']=fp.max(1)
    result['weighted_fine_group']=np.array(groups)[gp.argmax(1)]
    result['fine_neighbor_weight_fraction']=gp.max(1)
    result['nearest_reference_similarity']=sim
    result['nearest_reference_id']=rmeta.index.to_numpy()[ix[:,0]]
    result['effective_reference_neighbors']=1/np.square(w).sum(1)
    for library in sorted(rmeta['sample'].unique()):
        sub=np.flatnonzero(rmeta['sample'].eq(library))
        rfeat,qfeat=unit_features(rn[sub],qn,True)
        si,sw,_=neighbors(rfeat,qfeat)
        sf,sp=posteriors(rmeta.iloc[sub].mapping_family.astype(str).to_numpy(),si,sw)
        sg,sgp=posteriors(rmeta.iloc[sub].ventricular_group.astype(str).to_numpy(),si,sw)
        result[f'family_from_{library}']=np.array(sf)[sp.argmax(1)]
        result[f'group_from_{library}']=np.array(sg)[sgp.argmax(1)]
    family_columns=[c for c in result if c.startswith('family_from_')]
    group_columns=[c for c in result if c.startswith('group_from_')]
    result['same_family_from_both_libraries']=result[family_columns].eq(result.weighted_family,axis=0).all(axis=1)
    result['same_group_from_both_libraries']=result[group_columns].eq(result.weighted_fine_group,axis=0).all(axis=1)
    # This is an evidence flag, not a calibrated cross-platform probability.
    result['family_concentrated_and_library_consistent']=(result.family_neighbor_weight_fraction>=.7)&result.same_family_from_both_libraries
    result['group_concentrated_and_library_consistent']=(result.fine_neighbor_weight_fraction>=.7)&result.same_group_from_both_libraries&result.family_concentrated_and_library_consistent
    result['family_candidate']=np.where(result.family_concentrated_and_library_consistent,result.weighted_family,'Ambiguous')
    result['fine_group_candidate']=np.where(result.group_concentrated_and_library_consistent,result.weighted_fine_group,'Ambiguous')
    for j,label in enumerate(fam):result['neighbor_weight_family__'+label]=fp[:,j]
    # Keep whole section geometry/QC flags without silently deleting exclusions.
    joined=allobs[['x_centroid','y_centroid','spatial_lineage_candidate','pass_expression_QC','pass_transfer_information']].copy()
    joined['included_in_weighted_mapping']=joined.index.isin(result.index)
    joined=joined.join(result.drop(columns=['x_centroid','y_centroid']),how='left')
    joined.to_csv(OUT/stage/'weighted_CM_annotations.tsv.gz',sep='\t',na_rep='NA')
    np.savez_compressed(OUT/stage/'weighted_reference_neighbors.npz',query_cell_ids=qobs.index.to_numpy(dtype=str),
        reference_cell_ids=rmeta.index.to_numpy(dtype=str),neighbor_indices=ix,neighbor_weights=w,family_classes=np.array(fam),
        fine_classes=np.array(groups),family_weight_fractions=fp,fine_weight_fractions=gp)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stages',nargs='+',default=STAGES,choices=STAGES);args=parser.parse_args()
    upstream=BASE/'tools/cytospace'
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=upstream,text=True).strip()=='bb7fe12ab4466ef5e38dcac4c3ece07012060d60'
    subprocess.run(['git','diff','--exit-code','HEAD','--','cytospace'],cwd=upstream,check=True,stdout=subprocess.DEVNULL)
    OUT.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
    with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
        raw,corrected,genes,ids=b['original'],b['corrected'],b['genes'].astype(str),b['cell_ids'].astype(str)
    np.testing.assert_array_equal(ids,meta.index)
    development=[g for g in HOLDOUT if g in genes]
    # Frozen new test genes: deterministic hash order, excluding prior marker
    # holdouts and direct gate markers. This list is fixed before new fits.
    excluded=set(development)|set('Tnnt2 Actc1 Myl2 Myl4 Myl7 Nppa Nppb Pecam1 Kdr Cdh5 Npr3 Upk3b Wt1'.split())
    test=sorted([g for g in genes if g not in excluded],key=lambda g:hashlib.sha256(('heart-alignment-audit-v1:'+g).encode()).hexdigest())[:20]
    held=development+test;training=~np.isin(genes,held)
    design={'development_genes':development,'new_test_genes':test,'training_genes':genes[training].tolist(),
        'methods':METHODS,'k':30,'temperature':.10,'gene_sd_floor':.25,'feature_clip':5,
        'spatial_neighbors':12,'spatial_radius_um':75,'distance_scale_um':30,'spatial_mix':.25,
        'selection':'No hyperparameters selected using new test genes. All fixed methods reported.',
        'normalization':'log2(1+CPM), denominator uses training genes only for held-out evaluation',
        'author_spatial_labels_used_in_fit':False,'native_reference_abundance_constraints':True,
        'weighted_projection_reference_abundance_constraints':False,
        'cohort':'Previously frozen independently marker-gated ventricular/QC cohort; not author final states. Conditional gene prediction, not blinded end-to-end annotation.',
        'source_reference_sha256':digest(SOURCE/'reference_shared_panel_counts.npz')}
    existing=OUT/'benchmark_design.json'
    if existing.exists():assert json.loads(existing.read_text())==design
    else:existing.write_text(json.dumps(design,indent=2)+'\n')
    valid=meta.broad_class.eq('Ventricular-enriched CM')&(corrected.sum(1)>=20)&((corrected>0).sum(1)>=8)&(raw.sum(1)>=20)&((raw>0).sum(1)>=8)
    for stage in args.stages:
        folder=OUT/stage;folder.mkdir(exist_ok=True)
        if (folder/'completion.json').exists():print('SKIP_COMPLETED',stage,flush=True);continue
        start=time.time();print('START',stage,flush=True)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        elig=pd.read_csv(ROOT/'results/cytospace_xenium'/stage/'all_cells_eligibility.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        np.testing.assert_array_equal(elig.index,data.obs_names)
        qi=np.flatnonzero(elig.cytospace_eligible&elig.spatial_lineage_candidate.eq('Ventricular-enriched CM'))
        ri=np.flatnonzero(valid&meta.age.eq(stage))
        q=data[qi,genes].X.toarray().astype(np.float32);r=corrected[ri].astype(np.float32)
        # A train-zero row would be excluded identically from all methods.
        qkeep=q[:,training].sum(1)>0;rkeep=r[:,training].sum(1)>0
        qi,ri,q,r=qi[qkeep],ri[rkeep],q[qkeep],r[rkeep]
        qobs=data.obs.iloc[qi];rmeta=meta.iloc[ri]
        qn,rn=normalized(q,training),normalized(r,training)
        hi=np.array([np.flatnonzero(genes==g)[0] for g in held]);truth=qn[:,hi]
        native_ix=native(stage,q[:,training],r[:,training],genes[training],qobs.index.to_numpy(),rmeta.index.to_numpy(),qobs[['x_centroid','y_centroid']].to_numpy())
        predictions=[rn[native_ix][:,hi]]
        rf,qf=unit_features(rn[:,training],qn[:,training],False)
        ix,w,_=neighbors(rf,qf)
        predictions.extend([rn[ix[:,0]][:,hi],project(rn[:,hi],ix,w)])
        rf,qf=unit_features(rn[:,training],qn[:,training],True)
        ix,w,_=neighbors(rf,qf)
        scaled=project(rn[:,hi],ix,w);predictions.append(scaled)
        coords=qobs[['x_centroid','y_centroid']].to_numpy()
        for positions in [coords,coords[np.random.default_rng(931).permutation(len(coords))]]:
            g=graph(positions,qf);mass=np.asarray(g.sum(1)).ravel()
            smoothed=np.asarray(g@scaled)
            predictions.append(np.where(mass[:,None]>0,.75*scaled+.25*smoothed,scaled))
        predictions=np.stack(predictions)
        results=pd.DataFrame(metrics(stage,qobs,truth,predictions,held,q[:,hi],q[:,training].sum(1),(q[:,training]>0).sum(1)))
        results['gene_set']=np.where(results.gene.isin(test),'new_test','development')
        results.to_csv(folder/'heldout_gene_metrics.tsv',sep='\t',index=False)
        np.savez_compressed(folder/'heldout_predictions.npz',query_cell_ids=qobs.index.to_numpy(dtype=str),
            withheld_genes=np.array(held),methods=np.array(METHODS),predictions=predictions,measured_log2CPM=truth)
        print('HELDOUT_DONE',stage,flush=True)
        cross=library_benchmark(stage,r,rmeta)
        pd.DataFrame(cross).to_csv(folder/'cross_library_metrics.tsv',sep='\t',index=False)
        print('LIBRARY_DONE',stage,flush=True)
        final=final_mapping(stage,q,r,qobs,rmeta,data.obs)
        final.groupby(['family_candidate','fine_group_candidate']).size().rename('cells').to_csv(folder/'candidate_counts.tsv',sep='\t')
        (folder/'completion.json').write_text(json.dumps({'stage':stage,'query_cells':len(q),'reference_nuclei':len(r),
            'all_original_cells':len(data),'training_genes':int(training.sum()),'withheld_genes':len(held),
            'excluded_zero_training_query':int((~qkeep).sum()),'excluded_zero_training_reference':int((~rkeep).sum()),
            'elapsed_seconds':time.time()-start,'source_query_sha256':digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad'),
            'author_labels_used_in_fit':False,'production_projection':'assay_scaled_soft30, all381 genes, unsmoothed; candidate interpretation contingent on benchmark'},indent=2)+'\n')
        print('DONE',stage,round(time.time()-start,1),flush=True)
        del data,q,r,qn,rn,predictions


if __name__=='__main__':main()
