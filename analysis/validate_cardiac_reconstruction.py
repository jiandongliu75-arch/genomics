"""Independent source, score, surface-distance, model and inference checks."""
from reconstruct_cardiac_spatial import ROOT,OUT,UP,CFG,STAGES,MAT,read,save,dump,polys,sha,native_path
from test_cardiac_anatomy_mechanisms import exposures,covariates,demean
from cardiac_developmental_ordering import MODULES
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.interpolate import PchipInterpolator
from scipy.sparse.csgraph import breadth_first_order
from shapely import line_interpolate_point,distance
from PIL import Image
import gc,json,subprocess,re

def validate():
    rng=np.random.default_rng(691);checks=[];cal=read(OUT/'score_reference_calibration.tsv',index_col=None)
    expected={'P0':6298,'P7':11086,'P14':22532,'P21':24588};counts=read(OUT/'cohort_reconstruction_counts.tsv',index_col=None).set_index('stage')
    for stage in STAGES:
        print('VALIDATE_RECONSTRUCTION',stage,flush=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz');q=read(OUT/stage/'all_profile_observations.tsv.gz');c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');e=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz')
        assert q.index.equals(f.index) and q.index.is_unique
        np.testing.assert_allclose(q[['x_centroid','y_centroid']],f[['x_centroid','y_centroid']],atol=1e-6)
        assert q.author_Figure4d_CM.sum()==expected[stage]==counts.loc[stage,'published_vCM_Figure4d']
        assert q.author_Figure4e_label_eligible.equals(q.author_atlas_cell_state.isin(CFG['CM_author_states'][stage]))
        assert q.reconstructed_Figure4e_CM.equals(q.author_Figure4e_label_eligible&q.ventricular_morphology_selector)
        assert not q.loc[~q.ventricular_morphology_selector,'morphology_region'].isin(['Compact-compatible dense myocardium','Cavity-associated porous myocardium']).any()
        assert c.index.equals(e.index)
        for label in ['author_capillary','author_capillary_plus_HighMyo','QC_capillary','author_pericyte']:
            conflict=e['source_role_conflict__'+label]
            assert e.loc[conflict,label+'__apposition_fraction_2um'].isna().all()
        a=ad.read_h5ad(native_path(stage),backed='r');b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        assert a.var_names.equals(b.var_names);genes=a.var_names.to_numpy(dtype=str);ids=rng.choice(c.index,32,replace=False);ix=a.obs_names.get_indexer(ids)
        for view,x in [('native',a.X[ix]),('core',b.layers['high_confidence'][ix])]:
            dense=x.toarray();l=np.log1p(dense*1e4/np.maximum(dense.sum(1),1)[:,None]);v=cal[cal.view.eq(view)&cal.scope.eq('all')].set_index('gene').loc[genes]
            z=(l-v.reference_mean.to_numpy())/v.reference_SD.to_numpy()
            for name,panel in {'maturation21':MAT.gene.tolist(),**CFG['CM_components']}.items():
                np.testing.assert_allclose(z[:,np.isin(genes,panel)].mean(1),c.loc[ids,view+'_Z__'+name],rtol=3e-6,atol=2e-7)
            np.testing.assert_allclose(dense.sum(1),q.loc[ids,view+'_panel_depth'])
        a.file.close();b.file.close()
        with np.load(OUT/stage/'morphology_masks.npz') as masks:
            grid=float(masks['resolution_um']);xy=q[['x_centroid','y_centroid']].to_numpy();ij=np.rint(xy[:,::-1]/grid).astype(int);ij=np.clip(ij,[0,0],np.array(masks['envelope'].shape)-1)
            np.testing.assert_array_equal(q.ventricular_morphology_selector,masks['envelope'][tuple(ij.T)]&masks['tissue'][tuple(ij.T)])
        pp,path=polys(stage,q.index);source=np.flatnonzero(q.author_atlas_cell_state.eq('cap.EC').to_numpy());sourcepolys=np.array([pp[i] for i in source],dtype=object)
        # Brute-force every source polygon for two sampled CM boundaries. This
        # independently checks the STRtree distances and point-order handling.
        chosen=rng.choice(c.index[c.reconstructed_Figure4e_CM],2,replace=False)
        for cid in chosen:
            poly=pp[q.index.get_loc(cid)];pts=line_interpolate_point(poly.boundary,(np.arange(64)+.5)/64,normalized=True)
            nearest=np.array([np.min(distance(p,sourcepolys)) for p in pts])
            for cutoff in [1,2,5]:np.testing.assert_allclose(np.mean(nearest<=cutoff),e.loc[cid,f'author_capillary__apposition_fraction_{cutoff}um'],atol=1e-8)
            np.testing.assert_allclose(nearest.min(),e.loc[cid,'author_capillary__minimum_boundary_distance_um'],rtol=3e-6,atol=1e-7)
            np.testing.assert_allclose(poly.length,e.loc[cid,'CM_native_perimeter_um'],rtol=3e-6)
        g=sparse.load_npz(UP/stage/'proximity50_graph.npz').tocsr();ci=q.index.get_indexer(ids);n=np.asarray(g[ci]@q.author_atlas_cell_state.eq('cap.EC').to_numpy(float)).ravel()
        np.testing.assert_allclose(n,e.loc[ids,'author_capillary__nearby_profiles_50um'])
        ring=read(OUT/stage/'reconstructed_neighbor_levels.tsv.gz');g=sparse.load_npz(OUT/stage/'paper_like_symmetric6NN_graph.npz').tocsr()
        assert (g-g.T).nnz==0 and g.diagonal().sum()==0
        # Independent bounded BFS from a focal CM to capillary labels.
        seeds=set(source)
        for cid in ids[:4]:
            frontier={q.index.get_loc(cid)};seen=set(frontier);answer=np.nan
            for level in range(11):
                if frontier&seeds:answer=level;break
                nxt=set()
                for node in frontier:nxt.update(g.indices[g.indptr[node]:g.indptr[node+1]])
                frontier=nxt-seen;seen.update(frontier)
                if not frontier:break
            stored=ring.loc[cid,'paper_like6NN__cap.EC']
            if np.isfinite(answer):np.testing.assert_allclose(min(answer,10),stored)
            else:assert stored==10 or pd.isna(stored)
        checks.append({'stage':stage,'profiles':len(q),'Figure4d_count_exact':expected[stage],'score_rows_independently_recomputed':32,
            'CM_surface_samples_brute_force_checked':128,'independent_graph_BFS_checks':4,'original_ids_coordinates_preserved':True,'status':'PASS'})
        del pp,sourcepolys,q,c,e,f;gc.collect()
    # Recalculate a primary P7 coefficient by independent residualization.
    c=read(OUT/'P7/CM_outcomes_and_anatomy.tsv.gz');en=read(OUT/'P7/vascular_interface_and_matrix_context.tsv.gz');e=exposures(c,en);b=covariates(c,en,True)
    valid=c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20)&e.notna().all(axis=1)&b.notna().all(axis=1)
    c=c.loc[valid];b=b.loc[valid];e=e.loc[valid];name='capillary_apposition';outcome='core_Z__metabolic'
    v=np.column_stack([b,e]);sd=v.std(0);use=sd>1e-9;vv=(v[:,use]-v[:,use].mean(0))/sd[use]
    names=np.array(list(b)+list(e))[use];k=int(np.flatnonzero(names==name)[0]);vv=demean(vv,c.anatomy_block_200um.to_numpy())
    yy=c[outcome].to_numpy();yy=(yy-yy.mean())/yy.std();yy=demean(yy[:,None],c.anatomy_block_200um.to_numpy()).ravel()
    xx=vv[:,k];bb=np.column_stack([np.ones(len(c)),np.delete(vv,k,axis=1)]);xr=xx-bb@np.linalg.lstsq(bb,xx,rcond=1e-9)[0];yr=yy-bb@np.linalg.lstsq(bb,yy,rcond=1e-9)[0]
    estimate=xr@yr/(xr@xr);coeff=read(OUT/'anatomy_conditioned_coefficients.tsv.gz',index_col=None)
    row=coeff[coeff.stage.eq('P7')&coeff.cohort.eq('paper_label_mask')&coeff.capillary_source.eq('author_capillary')&coeff.adjustment.eq('within_region_capture')&coeff.exposure.eq(name)&coeff.outcome.eq(outcome)].iloc[0]
    np.testing.assert_allclose(estimate,row.standardized_coefficient,rtol=3e-6,atol=1e-8)
    cov=read(OUT/'shared_anatomy_covariance_partition.tsv',index_col=None)
    np.testing.assert_allclose(cov.total_covariance,cov.within_200um_covariance+cov.between_200um_covariance,rtol=1e-6,atol=1e-9)
    folds=read(OUT/'spatial_prediction_fold_checks.tsv',index_col=None);assert folds.minimum_train_test_separation_um.ge(CFG['prediction_buffer_um']).all()
    pairs=read(OUT/'matched_apposition_CM_pairs.tsv.gz',index_col=None)
    for stage,z in pairs.groupby('stage'):
        assert z.high_apposition_CM.is_unique and z.zero_apposition_CM.is_unique
        assert z.centroid_separation_um.le(100).all() and z.matching_distance.le(4).all()
        c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz')
        np.testing.assert_array_equal(c.loc[z.high_apposition_CM,'anatomy_block_200um'],c.loc[z.zero_apposition_CM,'anatomy_block_200um'])
        np.testing.assert_array_equal(c.loc[z.high_apposition_CM,'author_atlas_cell_state'],c.loc[z.zero_apposition_CM,'author_atlas_cell_state'])
        env=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz');dd=c.loc[c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20)]
        variables=pd.DataFrame({'depth':np.log1p(dd.native_panel_depth),'area':np.log1p(dd.cell_area),
            'perimeter':np.log1p(env.loc[dd.index,'CM_native_perimeter_um']),'capillary_area':np.log1p(env.loc[dd.index,'author_capillary__nearby_cell_area_fraction_50um'])})
        std=variables.to_numpy().std(0)
        diff=abs(variables.loc[z.high_apposition_CM].to_numpy()-variables.loc[z.zero_apposition_CM].to_numpy())/std
        assert diff[:,-1].max()<=CFG['matched_maximum_capillary_area_SD']+1e-8
        assert diff[:,:3].max()<=CFG['matched_maximum_capture_shape_SD']+1e-8
    balance=read(OUT/'matched_interface_covariate_balance.tsv',index_col=None)
    assert balance.matched_mean_difference_SD.abs().max()<.05
    # Ordering was fitted without the validation/endpoint genes, in two
    # disjoint biological RNA replicate series as well as the pooled series.
    trajectory=read(OUT/'RNA_cell_developmental_ordering.tsv.gz');observed=read(OUT/'RNA_ventricular_observed_programs.tsv.gz')
    assert trajectory.age.isin(STAGES).all() and trajectory.index.isin(observed.index).all()
    rep1=trajectory.index[trajectory.ordering.eq('RNA_replicate1')];rep2=trajectory.index[trajectory.ordering.eq('RNA_replicate2')]
    assert len(rep1.intersection(rep2))==0 and len(rep1)+len(rep2)==len(observed)
    excluded=set(sum(MODULES.values(),[]))
    for name in ['pooled','RNA_replicate1','RNA_replicate2']:
        with np.load(OUT/(name+'_RNA_ordering_model.npz')) as m:
            assert not set(m['HVG_symbols'])&excluded
            root=str(m['root_id']);r=trajectory[trajectory.ordering.eq(name)].loc[root]
            assert r.age=='P0' and abs(r.graph_geodesic_order)<1e-8
    estimates=read(OUT/'inferred_P3_P5_P10_programs.tsv',index_col=None);assert estimates.new_biological_replicates.eq(0).all() and not estimates.observed_at_this_day.any()
    lib=read(OUT/'RNA_observed_library_programs.tsv',index_col=None)
    for (scope,program),z in estimates.groupby(['scope','program']):
        v=lib[lib.scope.eq(scope)&lib.program.eq(program)].groupby('day').pseudobulk_logCPM_mean.mean().reindex([0,7,14,21])
        np.testing.assert_allclose(np.interp(z.inferred_day,[0,7,14,21],v),z.linear_estimate,rtol=3e-6)
        np.testing.assert_allclose(PchipInterpolator([0,7,14,21],v)(z.inferred_day),z.PCHIP_estimate,rtol=3e-6)
    definition=json.loads((OUT/'author_definition.json').read_text());assert not definition['exact_Figure4e_replication_claimed']
    save(pd.DataFrame(checks),'validation_checks.tsv')
    # Link to previous immutable input hashes, rather than silently regenerating
    # a new source fingerprint after a possible original-count change.
    previous=read(ROOT/'results/cardiac_niche_dynamics/source_checksums.tsv',index_col=None)
    verified=0
    for stage in STAGES:
        for path in [native_path(stage),ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad']:
            hit=previous[previous.path.isin([str(path),str(path.relative_to(ROOT))])];assert len(hit)==1
            assert sha(path)==hit.iloc[0].sha256;verified+=1
    dump('numerical_validation.json',{'status':'PASS','stages':checks,'unchanged_count_files_compared_to_prior_release':verified,
        'coefficient_independently_residualized':True,'shared_anatomy_covariance_identity_checked':True,'spatial_holdout_buffer_checked':True,
        'RNA_replicate_series_disjoint':True,'maturation_validation_genes_withheld_from_ordering':True,'interpolation_recomputed':True,
        'interpolated_days_are_new_samples':False,'exact_author_Figure4e_inputs_recovered':False,
        'scope':'Numerical/source/geometry checks, not anatomical ground truth, new spatial replication, perfusion or causal validation'})
    print('RECONSTRUCTION VALIDATION PASS',flush=True)

if __name__=='__main__':validate()
