"""Independent accounting, molecule-score and physical-neighborhood checks."""
import hashlib
import re
import json
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from anndata.io import read_elem
from scipy.spatial import cKDTree
from fibroblast_state_common import *


def main():
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    clusters=pd.read_csv(OUT/'RNA_cluster_summary.tsv',sep='\t')
    assert ann.index.is_unique and len(ann)==14985
    assert len(clusters)==31 and clusters.nuclei.sum()==len(ann)
    assert ann.lineage_mixing_review_flag.sum()==657
    with h5py.File(BASE/'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad') as f:
        original=read_elem(f['obs'])
    assert set(ann.index)==set(original.index[original.reviewed_identity.eq('Fibroblast')])
    np.testing.assert_allclose(ann.total_counts,original.loc[ann.index,'total_counts'])
    meta=pd.read_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t',index_col=0)
    assert meta.nuclei.sum()==14985
    for layer,col in [('corrected','total_counts'),('raw','raw_total_counts')]:
        pb=pd.read_csv(OUT/f'RNA_{layer}_pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
        key=ann['sample'].astype(str)+'::'+ann.fibroblast_cluster_id.astype(str)
        expected=original.loc[ann.index,col].groupby(key).sum().reindex(pb.columns)
        np.testing.assert_allclose(pb.sum(axis=0),expected,rtol=3e-6,atol=.1)
        if layer=='raw':assert np.allclose(pb.to_numpy(),np.rint(pb.to_numpy()))
    defs=pd.read_csv(OUT/'RNA_pseudobulk_program_definitions.tsv',sep='\t')
    for r in defs[defs.version.eq('annotation_markers_withheld')].itertuples():
        assert not set(r.genes.split(';'))&ANNOTATION_GENES
    stats=pd.read_csv(OUT/'spatial_adjusted_coefficients.tsv',sep='\t')
    assert np.isfinite(stats.standardized_coefficient).all()
    assert stats.standardized_coefficient.abs().max()<2,'unstable coefficient/design'
    assert not ((stats.stage.isin(['P0','P7']))&stats.predictor.eq('family__Postn-associated')).any()
    rng=np.random.default_rng(617)
    summary=[]
    controls=pd.read_csv(OUT/'CM_maturation_fixed_controls.tsv',sep='\t')
    assert controls.maturation_gene.sum()==21
    assert not set(controls.loc[controls.matched_control,'gene'])&set(WITHHOLD)
    for stage in STAGES:
        folder=OUT/stage
        atlas=pd.read_csv(ROOT/'results/whole_heart_spatial_atlas'/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        fb=pd.read_csv(folder/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        assert set(fb.index)==set(atlas.index[atlas.atlas_identity.eq('Fibroblast')])
        assert set(cm.index)==set(atlas.index[atlas.ventricular_CM_strict])
        np.testing.assert_allclose(fb[['x_centroid','y_centroid']],atlas.loc[fb.index,['x_centroid','y_centroid']])
        np.testing.assert_allclose(cm[['x_centroid','y_centroid']],atlas.loc[cm.index,['x_centroid','y_centroid']])
        p=fb[['weight__'+f for f in FAMILIES+[UNRESOLVED]]].to_numpy()
        assert np.allclose(p.sum(1),1,atol=1e-5)
        assert fb.loc[fb.supported_family_assignment,'measured_state_markers'].ge(2).all()
        assert fb.loc[fb.supported_family_assignment,'both_library_agreement'].all()
        with np.load(folder/'fibroblast_projection_calibration.npz') as b:
            assert not set(b['training_genes'].astype(str))&set(WITHHOLD)
            assert set(b['reference_ids'].astype(str)).issubset(set(ann.index[ann.fine_reference_eligible]))
        source=ad.read_h5ad(ROOT/'results/xenium_ventricular_mapping'/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad',backed='r')
        genes=source.var_names.to_numpy(dtype=str)
        selected=rng.choice(len(cm),min(64,len(cm)),replace=False)
        ids=cm.index[selected]
        x=source.X[source.obs_names.get_indexer(ids)].toarray().astype(float)
        logs=np.log1p(1e4*x/np.maximum(x.sum(1),1)[:,None])
        expected=logs[:,np.isin(genes,MATURATION)].mean(1)
        np.testing.assert_allclose(expected,cm.loc[ids,'maturation_21gene_mean'],rtol=1e-6)
        expected=expected-logs[:,np.isin(genes,controls.loc[controls.matched_control,'gene'])].mean(1)
        np.testing.assert_allclose(expected,cm.loc[ids,'maturation_matched_control'],rtol=1e-6,atol=1e-7)
        ank=logs[:,list(genes).index('Ankrd1')]
        np.testing.assert_allclose(ank,cm.loc[ids,'Ankrd1_full_denominator'],rtol=1e-6)
        tree=cKDTree(fb[['x_centroid','y_centroid']].to_numpy())
        xy=cm.loc[ids,['x_centroid','y_centroid']].to_numpy()
        for radius in [30,50,100]:
            env=pd.read_csv(folder/f'CM_neighborhoods_{radius}um.tsv.gz',sep='\t',index_col=0)
            neighbors=tree.query_ball_point(xy,radius)
            for cid,ii in zip(ids,neighbors):
                assert len(ii)==env.loc[cid,'fibroblast_neighbors']
                if len(ii):
                    for fam in FAMILIES:
                        np.testing.assert_allclose(fb.iloc[ii]['weight__'+fam].mean(),env.loc[cid,'family__'+fam],rtol=1e-5,atol=1e-8)
                    np.testing.assert_allclose(fb.iloc[ii].matrix__Basement_membrane.mean(),env.loc[cid,'matrix__Basement_membrane'],rtol=1e-5,atol=1e-8)
        # Independent reconstruction of the physical train/test separation counts.
        env=pd.read_csv(folder/'CM_neighborhoods_50um.tsv.gz',sep='\t',index_col=0)
        d=cm.join(env);d=d[d.fibroblast_neighbors.ge(3)&d.resolved_neighbors.ge(5)&d.arterial_distance_um.notna()]
        xy=d[['x_centroid','y_centroid']].to_numpy()
        bx=np.floor(xy[:,0]/400).astype(int);by=np.floor(xy[:,1]/400).astype(int)
        fold=((bx*73856093)^(by*19349663)^271)%5
        saved=pd.read_csv(OUT/'spatial_blocked_prediction_folds.tsv',sep='\t')
        saved=saved[saved.stage.eq(stage)&saved.outcome.eq('maturation_21gene_mean')&saved.feature_set.eq('matrix')&saved.model.eq('baseline')].set_index('fold')
        for k in range(5):
            test=fold==k;distance=cKDTree(xy[test]).query(xy[~test])[0]
            assert int(test.sum())==saved.loc[k,'test_cells']
            assert int((distance>=100).sum())==saved.loc[k,'training_cells']
        source.file.close()
        summary.append({'stage':stage,'broad_FB_profiles':len(fb),'strict_ventricular_CM_profiles':len(cm),
            'independent_CM_score_checks':len(ids),'independent_neighborhood_checks':len(ids)*3,
            'primary_spatial_CM_profiles':len(d),'primary_profiles_without_nuclear_area':int(d.nuclear_area_missing.sum())})
    # Exact same subsets must be used within each boundary comparison.
    bd=pd.read_csv(OUT/'boundary_spatial_association_sensitivity.tsv',sep='\t')
    assert bd.groupby(['stage','scope','outcome','predictor']).same_CM_subset.nunique().le(1).all()
    methods=json.loads((OUT/'association_methods.json').read_text())
    methods['composition_reference']='Postn omitted from early simplex models; unresolved-program remainder is the late reference; constant features dropped; Gram pseudoinverse cutoff1e-9'
    write_json(OUT/'association_methods.json',methods)
    write_json(OUT/'validation.json',{'status':'PASS','all_14985_source_FB_nuclei_accounted_for':True,
        '31_clusters_accounted_for':True,'657_lineage_review_profiles_preserved':True,
        'original_and_corrected_pseudobulk_count_accounting':True,'matrix_annotation_genes_withheld':True,
        'outcome_genes_excluded_from_fine_transfer':True,'coordinates_preserved':True,
        'independent_native_CM_score_and_physical_neighbor_checks':summary,
        'spatial_block_fold_sizes_and_100um_training_buffer_reconstructed':True,
        'same_boundary_comparison_subsets_verified':True,'numerically_stable_compositional_design':True})
    # File-manifest hashes are evidence of artifact identity, not biological validity.
    records=[]
    for p in sorted(OUT.rglob('*')):
        if not p.is_file() or p.name=='file_manifest.tsv':continue
        records.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    for p in sorted((ROOT/'analysis').glob('*fibroblast*.py')):
        records.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    assert max(r['bytes'] for r in records)<95*1024**2
    save_tsv(pd.DataFrame(records),OUT/'file_manifest.tsv')
    print('VALIDATION_PASS',summary,flush=True)


if __name__=='__main__':main()
