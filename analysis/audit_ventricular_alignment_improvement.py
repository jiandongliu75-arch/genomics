"""Audit the frozen held-out comparison and annotate technical uncertainty."""
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from improve_ventricular_alignment import (OUT, SOURCE, ROOT, STAGES, METHODS,
    normalized, unit_features, neighbors, posteriors, metrics, digest)


def main():
    design=json.loads((OUT/'benchmark_design.json').read_text())
    meta=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
    with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
        counts,raw,genes,ids=b['corrected'],b['original'],b['genes'].astype(str),b['cell_ids'].astype(str)
    np.testing.assert_array_equal(ids,meta.index)
    training=np.isin(genes,design['training_genes'])
    test_genes=design['development_genes']+design['new_test_genes']
    hi=np.array([np.flatnonzero(genes==g)[0] for g in test_genes])
    assert set(design['training_genes']).isdisjoint(test_genes)
    assert len(set(test_genes))==39 and training.sum()==342
    allmetrics,validation,coverage,library_calls=[],[],[],[]
    for stage in STAGES:
        folder=OUT/stage
        while not (folder/'completion.json').exists():
            print('WAIT',stage,flush=True);time.sleep(5)
        completion=json.loads((folder/'completion.json').read_text())
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        eligibility=pd.read_csv(ROOT/'results/cytospace_xenium'/stage/'all_cells_eligibility.tsv.gz',sep='\t',index_col=0)
        qi=np.flatnonzero(eligibility.cytospace_eligible&eligibility.spatial_lineage_candidate.eq('Ventricular-enriched CM'))
        ri=np.flatnonzero(meta.age.eq(stage)&meta.broad_class.eq('Ventricular-enriched CM')&(counts.sum(1)>=20)&((counts>0).sum(1)>=8)&(raw.sum(1)>=20)&((raw>0).sum(1)>=8))
        q=data[qi,genes].X.toarray().astype(np.float32);r=counts[ri].astype(np.float32)
        qkeep=q[:,training].sum(1)>0;rkeep=r[:,training].sum(1)>0
        qi,ri,q,r=qi[qkeep],ri[rkeep],q[qkeep],r[rkeep]
        qobs=data.obs.iloc[qi];rmeta=meta.iloc[ri]
        with np.load(folder/'heldout_predictions.npz') as b:
            np.testing.assert_array_equal(b['query_cell_ids'],qobs.index)
            np.testing.assert_array_equal(b['withheld_genes'],test_genes)
            np.testing.assert_array_equal(b['methods'],METHODS)
            prediction,truth=b['predictions'],b['measured_log2CPM']
        expected=np.log2(1+1e6*q[:,hi]/q[:,training].sum(1)[:,None]).astype(np.float32)
        np.testing.assert_allclose(truth,expected,rtol=2e-7,atol=1e-6)
        m=pd.DataFrame(metrics(stage,qobs,truth,prediction,test_genes,q[:,hi],q[:,training].sum(1),(q[:,training]>0).sum(1)))
        m['gene_set']=np.where(m.gene.isin(design['new_test_genes']),'new_test','development')
        m.to_csv(folder/'heldout_gene_metrics.tsv',sep='\t',index=False)
        allmetrics.append(m)
        native=pd.read_csv(folder/'native_heldout/assigned_locations.csv.gz',float_precision='round_trip').set_index('SpotID').loc[qobs.index]
        assert native.index.is_unique and native.CellType.eq('Ventricular_CM').all()
        np.testing.assert_array_equal(native[['x_centroid','y_centroid']],qobs[['x_centroid','y_centroid']])
        native_ids=rmeta.index.get_indexer(native.OriginalCID);assert (native_ids>=0).all()
        reuse=native.OriginalCID.value_counts().reindex(rmeta.index,fill_value=0)
        assert reuse.ge(1).all() and reuse.sum()==len(q)
        rn=normalized(r,training);qn=normalized(q,training)
        np.testing.assert_allclose(prediction[0],rn[native_ids][:,hi],rtol=2e-7,atol=1e-6)
        sampled=np.random.default_rng(617).choice(len(q),40,replace=False)
        for scaled,mi in [(False,2),(True,3)]:
            rf,qf=unit_features(rn[:,training],qn[:,training],scaled)
            # Independent full sorting rather than production argpartition.
            score=qf[sampled]@rf.T
            ix=np.argsort(-score,axis=1,kind='stable')[:,:30]
            s=np.take_along_axis(score,ix,axis=1)
            weights=np.exp((s-s[:,:1])/.1);weights/=weights.sum(1,keepdims=True)
            expected=np.sum(weights[:,:,None]*rn[ix][:,:,hi],axis=1)
            np.testing.assert_allclose(prediction[mi,sampled],expected,rtol=2e-5,atol=2e-5)
            if not scaled:
                np.testing.assert_allclose(prediction[1,sampled],rn[ix[:,0]][:,hi],rtol=2e-7,atol=1e-6)
        for mi in [4,5]:
            assert prediction[mi].shape==truth.shape and np.isfinite(prediction[mi]).all()
        out=pd.read_csv(folder/'weighted_CM_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        np.testing.assert_array_equal(out.index,data.obs_names)
        np.testing.assert_array_equal(out[['x_centroid','y_centroid']],data.obs[['x_centroid','y_centroid']])
        assert set(out.index[out.included_in_weighted_mapping])==set(qobs.index)
        with np.load(folder/'weighted_reference_neighbors.npz') as b:
            np.testing.assert_array_equal(b['query_cell_ids'],qobs.index)
            np.testing.assert_array_equal(b['reference_cell_ids'],rmeta.index)
            ix,w=b['neighbor_indices'],b['neighbor_weights']
            production_pair_count=int(ix.size)
            np.testing.assert_allclose(w.sum(1),1,atol=3e-7)
            assert (ix>=0).all() and (ix<len(rmeta)).all() and (w>=0).all()
            fam,p=posteriors(rmeta.mapping_family.astype(str).to_numpy(),ix,w)
            group,gp=posteriors(rmeta.ventricular_group.astype(str).to_numpy(),ix,w)
            selected=out.loc[qobs.index]
            np.testing.assert_array_equal(np.array(fam)[p.argmax(1)],selected.weighted_family)
            np.testing.assert_allclose(p.max(1),selected.family_neighbor_weight_fraction,rtol=1e-6)
            np.testing.assert_array_equal(np.array(group)[gp.argmax(1)],selected.weighted_fine_group)
        fcols=[c for c in out if c.startswith('family_from_')];gcols=[c for c in out if c.startswith('group_from_')]
        assert len(fcols)==len(gcols)==2
        family_flag=out.included_in_weighted_mapping & out.family_neighbor_weight_fraction.ge(.7)&out[fcols].eq(out.weighted_family,axis=0).all(axis=1)
        group_flag=family_flag & out.fine_neighbor_weight_fraction.ge(.7)&out[gcols].eq(out.weighted_fine_group,axis=0).all(axis=1)
        out['family_concentrated_and_library_consistent']=family_flag
        out['group_concentrated_and_library_consistent']=group_flag
        out.loc[out.included_in_weighted_mapping,'family_candidate']=np.where(family_flag[out.included_in_weighted_mapping],out.loc[out.included_in_weighted_mapping,'weighted_family'],'Ambiguous')
        out.loc[out.included_in_weighted_mapping,'fine_group_candidate']=np.where(group_flag[out.included_in_weighted_mapping],out.loc[out.included_in_weighted_mapping,'weighted_fine_group'],'Ambiguous')
        out.to_csv(folder/'weighted_CM_annotations.tsv.gz',sep='\t',na_rep='NA')
        out.loc[out.included_in_weighted_mapping].groupby(['family_candidate','fine_group_candidate']).size().rename('cells').to_csv(folder/'candidate_counts.tsv',sep='\t')
        coverage.append({'stage':stage,'original_cells':len(out),'included_cells':len(q),
            'family_evidence_flag_cells':int(family_flag.sum()),'fine_evidence_flag_cells':int(group_flag.sum()),
            'family_flag_fraction':float(family_flag.sum()/len(q)),'fine_flag_fraction':float(group_flag.sum()/len(q))})
        for library in sorted(rmeta['sample'].unique()):
            query_index=np.flatnonzero(rmeta['sample'].eq(library));ref_index=np.flatnonzero(~rmeta['sample'].eq(library))
            train=np.ones(r.shape[1],bool)
            rf,qf=unit_features(normalized(r[ref_index],train),normalized(r[query_index],train),True)
            ix,w,_=neighbors(rf,qf)
            for field in ['mapping_family','ventricular_group']:
                classes,p=posteriors(rmeta.iloc[ref_index][field].astype(str).to_numpy(),ix,w)
                library_calls.append(pd.DataFrame({'stage':stage,'query_library':library,'resolution':field,
                    'query_cell_id':rmeta.index[query_index],'known_RNA_label':rmeta.iloc[query_index][field].to_numpy(),
                    'predicted_RNA_label':np.array(classes)[p.argmax(1)],'neighbor_weight_fraction':p.max(1)}))
        assert digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')==completion['source_query_sha256']
        validation.append({'stage':stage,'original_rows_conserved':len(out),'native_pairs_checked':len(q),
            'production_neighbor_pairs_checked':production_pair_count,
            'sampled_neighbor_projection_rows_independently_recalculated':len(sampled)*2,
            'heldout_metrics_recomputed':len(m),'native_withheld_predictions_recomputed':int(len(q)*39),
            'source_counts_and_coordinates_unchanged':True,'author_spatial_labels_used_in_fit':False})
        print('AUDITED',stage,coverage[-1],flush=True)
        del data,q,r,prediction,truth
    m=pd.concat(allmetrics,ignore_index=True);m.to_csv(OUT/'heldout_gene_metrics.tsv',sep='\t',index=False)
    fields=['cell_spearman','cell_pearson','depth_area_training_detection_adjusted_r','detection_AUC','grid250um_spearman']
    m.groupby(['stage','gene_set','method'])[fields].median().reset_index().to_csv(OUT/'heldout_summary.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(OUT/'mapping_coverage.tsv',sep='\t',index=False)
    calls=pd.concat(library_calls,ignore_index=True);calls.to_csv(OUT/'cross_library_predictions.tsv.gz',sep='\t',index=False)
    confidence=[]
    for (stage,resolution),part in calls.groupby(['stage','resolution']):
        for cutoff in [0,.5,.7,.9]:
            selected=part.neighbor_weight_fraction.ge(cutoff)
            confidence.append({'stage':stage,'resolution':resolution,'minimum_neighbor_weight_fraction':cutoff,
                'all_query_nuclei':len(part),'selected_nuclei':int(selected.sum()),
                'selected_fraction':float(selected.mean()),'label_agreement':float(part.loc[selected,'known_RNA_label'].eq(part.loc[selected,'predicted_RNA_label']).mean())})
    pd.DataFrame(confidence).to_csv(OUT/'cross_library_selective_agreement.tsv',sep='\t',index=False)
    assert digest(SOURCE/'reference_shared_panel_counts.npz')==design['source_reference_sha256']
    (OUT/'validation.json').write_text(json.dumps({'status':'PASS','stages':validation,
        'scope':'Numerical conservation and internal held-out predictions; not independent anatomical truth or animal-level validation.'},indent=2)+'\n')


if __name__=='__main__':main()
