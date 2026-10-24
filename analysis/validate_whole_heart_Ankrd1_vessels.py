"""Independent numerical checks of lineage, vascular geometry and local contrasts."""
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import shapely
from shapely.geometry import shape

from reassess_whole_heart_lineages import OUT,SOURCE,ROOT,STAGES,WITHHOLD
from improve_ventricular_alignment import normalized,unit_features,digest
from define_heart_arterial_structures import DEFINITIONS,read_boundaries
from test_Ankrd1_arterial_neighborhoods import masks


def main():
    metadata=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
    with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
        raw,corrected,genes,reference_ids=b['original'],b['corrected'],b['genes'].astype(str),b['cell_ids'].astype(str)
    np.testing.assert_array_equal(reference_ids,metadata.index)
    training=~np.isin(genes,WITHHOLD);assert training.sum()==377
    design=json.loads((OUT/'lineage_design.json').read_text())
    assert set(design['training_genes']).isdisjoint(WITHHOLD)
    assert set(design['training_genes'])==set(genes[training])
    assert digest(SOURCE/'reference_shared_panel_counts.npz')==design['reference_sha256']
    summaries,comparisons,contrasts,profiles,effects,previous_rows,validation,balances=[],[],[],[],[],[],[],[]
    for stage in STAGES:
        folder=OUT/stage
        obs=pd.read_csv(folder/'whole_heart_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        distance=pd.read_csv(folder/'arterial_boundary_distances.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        outcomes=pd.read_csv(folder/'association_outcomes.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        np.testing.assert_array_equal(obs.index,data.obs_names);np.testing.assert_array_equal(distance.index,obs.index)
        np.testing.assert_array_equal(obs[['x_centroid','y_centroid']],data.obs[['x_centroid','y_centroid']])
        np.testing.assert_array_equal(distance[['x_centroid','y_centroid']],obs[['x_centroid','y_centroid']])
        query=data[:,genes].X.toarray().astype(np.float32)
        expected=data.obs.pass_expression_QC.to_numpy()&(query[:,training].sum(1)>=20)&((query[:,training]>0).sum(1)>=8)
        np.testing.assert_array_equal(expected,obs.eligible_whole_reference)
        target=data[:,['Ankrd1']].X.toarray().ravel();np.testing.assert_array_equal(target,obs.measured_Ankrd1)
        denominator=np.asarray(data.X.sum(1)).ravel()-np.asarray(data[:,WITHHOLD].X.sum(1)).ravel()
        np.testing.assert_array_equal(denominator,obs.outcome_excluded_transcripts)
        np.testing.assert_allclose(np.log1p(1e4*target/np.maximum(denominator,1)),obs.Ankrd1_log1pCP10k,atol=1e-12)
        with np.load(folder/'whole_reference_weights.npz') as b:
            qids,rids=b['query_cell_ids'].astype(str),b['reference_cell_ids'].astype(str)
            ix,w,classes,p,states,sp=b['neighbor_indices'],b['neighbor_weights'],b['lineage_classes'].astype(str),b['lineage_weights'],b['state_classes'].astype(str),b['state_weights']
        np.testing.assert_array_equal(qids,obs.index[expected])
        assert pd.Index(rids).is_unique and pd.Index(rids).isin(metadata.index).all()
        rm=metadata.loc[rids];assert rm.age.eq(stage).all() and rm.reviewed_identity.ne('Unresolved/mixed').all()
        np.testing.assert_allclose(w.sum(1),1,atol=3e-7);assert (ix>=0).all() and (ix<len(rids)).all()
        y=rm.reviewed_identity.astype(str).to_numpy()
        inferred=np.column_stack([(w*(y[ix]==label)).sum(1) for label in classes])
        np.testing.assert_array_equal(inferred,p)
        np.testing.assert_allclose(p.sum(1),1,atol=3e-7);np.testing.assert_allclose(sp.sum(1),1,atol=3e-7)
        np.testing.assert_array_equal(obs.loc[qids,'whole_reference_lineage_proposal'],classes[p.argmax(1)])
        # Brute-force sorting checks the primary correspondence on sampled cells.
        ri=metadata.index.get_indexer(rids)
        rn=normalized(corrected[ri],training)[:,training];qn=normalized(query[expected],training)[:,training]
        rf,qf=unit_features(rn,qn,True)
        sampled=np.random.default_rng(407).choice(len(qids),30,replace=False)
        scores=qf[sampled]@rf.T;nearest=np.argsort(-scores,axis=1,kind='stable')[:,:30]
        value=np.take_along_axis(scores,nearest,axis=1);weights=np.exp((value-value[:,:1])/.1);weights/=weights.sum(1,keepdims=True)
        independent=np.column_stack([(weights*(y[nearest]==label)).sum(1) for label in classes])
        np.testing.assert_allclose(p[sampled],independent,rtol=3e-5,atol=3e-5)
        cm=obs.CM_strict
        assert obs.loc[cm,'whole_reference_lineage_proposal'].eq('Cardiomyocyte').all()
        assert obs.loc[cm,'both_RNA_libraries_agree'].eq(True).all() and obs.loc[cm,'raw_reference_agrees'].eq(True).all()
        assert obs.loc[cm,'lineage_weight'].ge(.7).all() and obs.loc[cm,'lineage_margin'].ge(.15).all()
        assert obs.loc[cm,'marker_Cardiomyocyte_genes'].ge(3).all() and obs.loc[cm,'marker_Cardiomyocyte_molecules'].ge(10).all()
        assert not obs.loc[cm,'mixed_CM_nonCM_marker_review'].any()
        assert obs.ventricular_CM_strict.eq(cm&obs.ventricular_chamber_marker_support).all()
        authors=obs.author_atlas_cell_state.isin(['v.CM','a.CM','av.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM'])
        for cohort,take in [('previous_CM_gate',obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')),('new_strict_CM',cm),('new_strict_ventricular',obs.ventricular_CM_strict)]:
            comparisons.append({'stage':stage,'cohort':cohort,'cells':int(take.sum()),'published_CM_cells':int((take&authors).sum()),
                'published_CM_fraction':float(authors[take].mean()),'published_nonCM_or_absent':int((take&~authors).sum())})
        features=json.loads((folder/'arterial_segments.geojson').read_text())['features']
        gcounts={}
        for definition in DEFINITIONS:
            selected=[f for f in features if f['properties']['definition']==definition]
            ec_used=[];sm_used=[]
            for f in selected:
                props=f['properties'];a=props['arterial_EC_ids'];s=props['SMC_ids'];ec_used+=a;sm_used+=s
                assert len(a)>=2 and len(s)>=2
                assert set(a+s).issubset(obs.index)
                ec=obs.loc[a,['x_centroid','y_centroid']].to_numpy();sm=obs.loc[s,['x_centroid','y_centroid']].to_numpy()
                coverage=float((cKDTree(sm).query(ec)[0]<=35).mean());assert abs(coverage-props['arterial_EC_mural_coverage'])<1e-12 and coverage>=.5
                assert (cKDTree(ec).query(sm)[0]<=35+1e-8).all()
                if definition=='independent_strict':assert obs.loc[a,'arterial_EC_seed'].all() and obs.loc[s,'SMC_seed'].all()
                if definition=='independent_marker':assert obs.loc[a,'arterial_EC_marker_seed'].all() and obs.loc[s,'SMC_marker_seed'].all()
                if definition=='published_comparator':assert obs.loc[a,'author_atlas_cell_state'].eq('arteriole.EC').all() and obs.loc[s,'author_atlas_cell_state'].eq('SMC').all()
                geom=shape(f['geometry']);assert geom.is_valid and geom.area>0
            assert len(set(ec_used))==len(ec_used) and len(set(sm_used))==len(sm_used)
            selected=[f for f in selected if f['properties']['primary_small_segment']]
            gcounts[definition]=len(selected)
            points=shapely.points(obs[['x_centroid','y_centroid']].iloc[sampled[:20]].to_numpy())
            for metric,key in [('boundary','geometry'),('unbuffered','unbuffered_geometry'),('hull','hull_geometry')]:
                saved=distance[definition+'__'+metric+'_distance_um'].iloc[sampled[:20]].to_numpy()
                if selected:
                    computed=np.min(np.column_stack([shapely.distance(points,shape(f[key])) for f in selected]),axis=1)
                    np.testing.assert_allclose(saved,computed,atol=1e-9)
                else:assert np.isnan(saved).all()
        # Source polygon reconstruction checks the µm geometry independently.
        samplefeatures=[]
        for definition in DEFINITIONS:
            group=[f for f in features if f['properties']['definition']==definition]
            if group:samplefeatures.append(group[0])
        wanted=set(c for f in samplefeatures for key in ['arterial_EC_ids','SMC_ids'] for c in f['properties'][key])
        polygons,_=read_boundaries(stage,wanted,obs)
        for f in samplefeatures:
            prop=f['properties'];union=shapely.union_all([polygons[c] for c in prop['arterial_EC_ids']+prop['SMC_ids']])
            assert union.buffer(3).symmetric_difference(shape(f['geometry'])).area<1e-6
        local=pd.read_csv(folder/'local_Ankrd1_contrasts.tsv',sep='\t')
        pairs=pd.read_csv(folder/'default_contrast_matched_pairs.tsv.gz',sep='\t')
        obs=obs.join(outcomes)
        qa_pairs=0
        for (definition,cohort),part in pairs.groupby(['definition','cohort']):
            near=obs.loc[part.near_cell_id];far=obs.loc[part.far_cell_id]
            assert part.near_cell_id.is_unique and part.far_cell_id.is_unique
            allowed=masks(obs)[cohort];assert allowed.loc[part.near_cell_id].all() and allowed.loc[part.far_cell_id].all()
            d=distance[definition+'__boundary_distance_um'];owner=distance[definition+'__boundary_segment_id']
            assert d.loc[part.near_cell_id].gt(0).all() and d.loc[part.near_cell_id].le(30).all()
            assert d.loc[part.far_cell_id].ge(60).all() and d.loc[part.far_cell_id].le(150).all()
            np.testing.assert_array_equal(owner.loc[part.near_cell_id],part.segment_id)
            np.testing.assert_array_equal(owner.loc[part.far_cell_id],part.segment_id)
            separation=np.linalg.norm(near[['x_centroid','y_centroid']].to_numpy()-far[['x_centroid','y_centroid']].to_numpy(),axis=1)
            assert (separation<=100+1e-8).all()
            for key in ['outcome_excluded_transcripts','cell_area']:
                ratio=near[key].to_numpy()/far[key].to_numpy();assert ((ratio>=.5)&(ratio<=2)).all()
            assert (np.abs(near.nearest_published_endocardial_EC_um.to_numpy()-far.nearest_published_endocardial_EC_um.to_numpy())<=50+1e-8).all()
            np.testing.assert_allclose(part.near_Ankrd1,near.Ankrd1_log1pCP10k,atol=1e-12)
            np.testing.assert_allclose(part.far_Ankrd1,far.Ankrd1_log1pCP10k,atol=1e-12)
            record=local[local.definition.eq(definition)&local.cohort.eq(cohort)&local.near_limit_um.eq(30)&local.distance_metric.eq('boundary')].iloc[0]
            assert record.matched_pairs==len(part)
            expected=np.expm1(part.near_Ankrd1).mean()/np.expm1(part.far_Ankrd1).mean()
            assert abs(expected-record.near_far_Ankrd1_CP10k_ratio)<1e-12
            mean=(part.near_Ankrd1-part.far_Ankrd1).groupby(part.segment_id).mean().mean()
            assert abs(mean-record.equal_segment_mean_Ankrd1_log_difference)<1e-12
            qa_pairs+=len(part)
            if definition=='independent_strict' and cohort=='whole_reference_strict':
                balances.append({'stage':stage,'pairs':len(part),'mean_pair_separation_um':part.pair_separation_um.mean(),
                    'median_depth_ratio':(part.near_depth/part.far_depth).median(),
                    'mean_log_depth_ratio':np.log(part.near_depth/part.far_depth).mean(),
                    'median_area_ratio':(part.near_area/part.far_area).median(),
                    'mean_log_area_ratio':np.log(part.near_area/part.far_area).mean(),
                    'mean_endocardial_difference_um':(part.near_endocardial_distance_um-part.far_endocardial_distance_um).mean()})
        lc=json.loads((folder/'lineage_completion.json').read_text());vc=json.loads((folder/'vessel_completion.json').read_text())
        assert digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')==lc['source_query_sha256']
        assert digest(folder/'whole_heart_cell_annotations.tsv.gz')==vc['lineage_annotations_sha256']
        previous=pd.read_csv(ROOT/'results/ventricular_spatial_resolution'/stage/'fine_CM_spatial_annotations.tsv.gz',sep='\t',index_col=0,usecols=['cell_id','fine_group_proposal','prioritized_fine_candidate'])
        for group,part in previous[previous.prioritized_fine_candidate].groupby('fine_group_proposal'):
            current=obs.loc[part.index]
            previous_rows.append({'stage':stage,'previous_fine_group':group,'previous_prioritized_cells':len(part),
                'new_strict_CM_cells':int(current.CM_strict.sum()),'new_inclusive_CM_cells':int(current.CM_inclusive.sum()),
                'published_CM_cells':int(authors.loc[part.index].sum()),'new_strict_and_published_CM':int((current.CM_strict&authors.loc[part.index]).sum()),
                'new_nonCM_resolved_cells':int(current.reassessed_lineage.isin(['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial','Glial']).sum()),
                'new_ambiguous_or_mixed_cells':int(current.reassessed_lineage.isin(['Ambiguous/review','Mixed CM/non-CM; review']).sum())})
        summaries.append({'stage':stage,'original_cells':len(obs),'eligible_query_cells':int(obs.eligible_whole_reference.sum()),
            'resolved_CMs':int(cm.sum()),'strict_ventricular_CMs':int(obs.ventricular_CM_strict.sum()),
            **{k+'_segments':v for k,v in gcounts.items()}})
        contrasts.append(local);profiles.append(pd.read_csv(folder/'distance_expression_profiles.tsv',sep='\t'))
        effects.append(pd.read_csv(folder/'segment_effects.tsv',sep='\t'))
        validation.append({'stage':stage,'original_cells_checked':len(obs),'weighted_reference_pairs_checked':int(ix.size),
            'sampled_primary_projections_independently_recomputed':len(sampled),'arterial_segments_checked':len(features),
            'source_vascular_envelopes_rebuilt':len(samplefeatures),'default_matched_pairs_checked':qa_pairs,
            'outcome_genes_excluded_from_identity_fits':True,'source_counts_and_centroids_unchanged':True})
        print('VALIDATED',stage,validation[-1],flush=True)
        del data,query,rn,qn,rf,qf
    pd.DataFrame(summaries).to_csv(OUT/'stage_summary.tsv',sep='\t',index=False)
    pd.DataFrame(comparisons).to_csv(OUT/'CM_lineage_concordance.tsv',sep='\t',index=False)
    pd.DataFrame(previous_rows).to_csv(OUT/'previous_fine_candidates_reassessment.tsv',sep='\t',index=False)
    pd.DataFrame(balances).to_csv(OUT/'primary_matching_balance.tsv',sep='\t',index=False)
    contrast=pd.concat(contrasts,ignore_index=True);contrast.to_csv(OUT/'local_Ankrd1_contrasts.tsv',sep='\t',index=False,na_rep='NA')
    contrast[contrast.definition.eq('independent_strict')&contrast.cohort.eq('whole_reference_strict')&contrast.distance_metric.eq('boundary')&contrast.near_limit_um.eq(30)].to_csv(OUT/'primary_local_Ankrd1_results.tsv',sep='\t',index=False,na_rep='NA')
    pd.concat(profiles,ignore_index=True).to_csv(OUT/'distance_expression_profiles.tsv',sep='\t',index=False,na_rep='NA')
    pd.concat(effects,ignore_index=True).to_csv(OUT/'segment_effects.tsv',sep='\t',index=False,na_rep='NA')
    (OUT/'validation.json').write_text(json.dumps({'status':'PASS','stages':validation,
        'scope':'Numerical/geometry conservation, lineage constraints and local paired-expression calculations; not independent anatomical or animal-level validation.'},indent=2)+'\n')


if __name__=='__main__':main()
