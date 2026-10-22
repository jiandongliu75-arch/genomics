"""Local, vessel-matched Ankrd1 contrasts across identity/geometry sensitivities."""
import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from reassess_whole_heart_lineages import OUT,ROOT,STAGES,digest
from define_heart_arterial_structures import DEFINITIONS

COHORTS=['whole_reference_strict','whole_reference_inclusive','annotation_consensus','published_ventricular','strict_Pearson_agreement','previous_CM_gate']
CONTRASTS=[(20,40,100),(30,60,150),(50,100,200)]


def masks(obs):
    return {'whole_reference_strict':obs.ventricular_CM_strict,
        'whole_reference_inclusive':obs.ventricular_CM_inclusive,
        'annotation_consensus':obs.CM_consensus_ventricular,
        'published_ventricular':obs.author_ventricular,
        'strict_Pearson_agreement':obs.ventricular_CM_strict&obs.Pearson_method_agrees.eq(True),
        'previous_CM_gate':obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')}


def match_pairs(obs,selected,distance,segments,near_limit,far_low,far_high):
    valid=selected&np.isfinite(distance)&(distance>0)
    near=valid&(distance<=near_limit);far=valid&(distance>=far_low)&(distance<=far_high)
    xy=obs[['x_centroid','y_centroid']].to_numpy()
    depth=obs.outcome_excluded_transcripts.to_numpy();area=obs.cell_area.to_numpy()
    endo=obs.nearest_published_endocardial_EC_um.to_numpy()
    rows=[]
    rng=np.random.default_rng(184)
    for segment in sorted(set(segments[near])):
        ni=np.flatnonzero(near&(segments==segment));fi=np.flatnonzero(far&(segments==segment))
        if not len(fi):continue
        available=np.ones(len(fi),bool)
        for i in rng.permutation(ni):
            candidates=np.flatnonzero(available)
            if not len(candidates):break
            j=fi[candidates]
            separation=np.linalg.norm(xy[j]-xy[i],axis=1)
            logdepth=np.log(np.maximum(depth[j],1)/max(depth[i],1))
            logarea=np.log(np.maximum(area[j],1)/max(area[i],1))
            endodiff=np.abs(endo[j]-endo[i])
            keep=(separation<=100)&(np.abs(logdepth)<=np.log(2))&(np.abs(logarea)<=np.log(2))&(endodiff<=50)
            if not keep.any():continue
            score=(separation[keep]/100)**2+(logdepth[keep]/.5)**2+(logarea[keep]/.5)**2+(endodiff[keep]/50)**2
            chosen=candidates[keep][np.argmin(score)];j=fi[chosen];available[chosen]=False
            rows.append({'near_cell_id':obs.index[i],'far_cell_id':obs.index[j],'segment_id':segment,
                'near_distance_um':float(distance[i]),'far_distance_um':float(distance[j]),
                'pair_separation_um':float(np.linalg.norm(xy[j]-xy[i])),
                'near_depth':float(depth[i]),'far_depth':float(depth[j]),'near_area':float(area[i]),'far_area':float(area[j]),
                'near_endocardial_distance_um':float(endo[i]),'far_endocardial_distance_um':float(endo[j]),
                'near_Ankrd1':float(obs.Ankrd1_log1pCP10k.iloc[i]),'far_Ankrd1':float(obs.Ankrd1_log1pCP10k.iloc[j]),
                'near_module':float(obs.response_module_log1pCP10k_mean.iloc[i]),'far_module':float(obs.response_module_log1pCP10k_mean.iloc[j]),
                'near_Ankrd1_high':bool(obs.Ankrd1_high.iloc[i]),'far_Ankrd1_high':bool(obs.Ankrd1_high.iloc[j]),
                'near_module_high':bool(obs.response_module_high.iloc[i]),'far_module_high':bool(obs.response_module_high.iloc[j])})
    return pd.DataFrame(rows),int(near.sum()),int(far.sum()),int((selected&np.isfinite(distance)&(distance==0)).sum())


def summarize_pairs(pairs):
    if pairs.empty:return {'matched_pairs':0,'segments_with_pairs':0}
    difference=pairs.near_Ankrd1-pairs.far_Ankrd1
    module=pairs.near_module-pairs.far_module
    group=pairs.assign(Ankrd1_difference=difference,module_difference=module).groupby('segment_id').agg(
        pair_count=('near_cell_id','size'),Ankrd1_difference=('Ankrd1_difference','mean'),module_difference=('module_difference','mean'))
    effects=group.Ankrd1_difference.to_numpy()
    p=np.nan
    if len(effects)>=6:
        rng=np.random.default_rng(1387)
        null=(rng.choice([-1,1],size=(4096,len(effects)))*effects).mean(1)
        p=float((1+(np.abs(null)>=abs(effects.mean())).sum())/4097)
    return {'matched_pairs':len(pairs),'segments_with_pairs':len(group),
        'pair_mean_Ankrd1_log_difference':float(difference.mean()),'pair_median_Ankrd1_log_difference':float(difference.median()),
        'equal_segment_mean_Ankrd1_log_difference':float(effects.mean()),
        'equal_segment_mean_module_log_difference':float(group.module_difference.mean()),
        'fraction_segments_positive_Ankrd1_difference':float((effects>0).mean()),
        'near_mean_Ankrd1_CP10k':float(np.expm1(pairs.near_Ankrd1).mean()),
        'far_mean_Ankrd1_CP10k':float(np.expm1(pairs.far_Ankrd1).mean()),
        'near_far_Ankrd1_CP10k_ratio':float(np.expm1(pairs.near_Ankrd1).mean()/max(np.expm1(pairs.far_Ankrd1).mean(),1e-12)),
        'near_Ankrd1_high_fraction':float(pairs.near_Ankrd1_high.mean()),'far_Ankrd1_high_fraction':float(pairs.far_Ankrd1_high.mean()),
        'near_module_high_fraction':float(pairs.near_module_high.mean()),'far_module_high_fraction':float(pairs.far_module_high.mean()),
        'exploratory_segment_signflip_p':p}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stages',nargs='+',default=STAGES,choices=STAGES);args=parser.parse_args()
    design={'outcomes':['Measured Ankrd1 log1pCP10k','Mean Creb5/Tlr4/Atf3 log1pCP10k'],
        'normalization':'403-target total excluding all4 outcome genes; unchanged original integer counts',
        'primary_cohort':'whole_reference_strict ventricular marker-supported CMs',
        'cohorts':COHORTS,'definitions':DEFINITIONS,'primary_definition':'independent_strict',
        'near_far_thresholds_um':[list(t) for t in CONTRASTS],'primary_contrast':'near0<d<=30µm, controls60<=d<=150µm',
        'pair_matching':'Same vascular segment,1:1 without replacement; physical separation<=100µm, depth/area ratios0.5–2,endocardial-distance difference<=50µm; best covariate-distance score; seed184',
        'regional_control':'Same-vessel local neighborhood and endocardial distance; no validated LV/RV or wall-depth anatomical masks',
        'high_expression':'Above the75th percentile in strict ventricular CMs at the same stage, before vessel association tests',
        'permutation':'Optional4096 segment-level sign flips when>=6 segments; exploratory within-section diagnostic, no animal-level inference; correlated vessel fragments can remain',
        'primary_summary':'Near/far expression ratio and equal-segment mean paired log-expression difference; every sensitivity reported without selection'}
    if (OUT/'association_design.json').exists():assert json.loads((OUT/'association_design.json').read_text())==design
    else:(OUT/'association_design.json').write_text(json.dumps(design,indent=2)+'\n')
    for stage in args.stages:
        folder=OUT/stage
        if (folder/'association_completion.json').exists():print('SKIP_ASSOCIATION',stage,flush=True);continue
        start=time.time()
        obs=pd.read_csv(folder/'whole_heart_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        distance=pd.read_csv(folder/'arterial_boundary_distances.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        np.testing.assert_array_equal(obs.index,distance.index)
        np.testing.assert_array_equal(obs[['x_centroid','y_centroid']],distance[['x_centroid','y_centroid']])
        endo=obs.author_atlas_cell_state.eq('endocardial.EC')
        assert endo.any()
        xy=obs[['x_centroid','y_centroid']].to_numpy()
        obs['nearest_published_endocardial_EC_um']=cKDTree(xy[endo]).query(xy)[0]
        primary=obs.ventricular_CM_strict
        aq=float(obs.loc[primary,'Ankrd1_log1pCP10k'].quantile(.75));mq=float(obs.loc[primary,'response_module_log1pCP10k_mean'].quantile(.75))
        obs['Ankrd1_high']=obs.Ankrd1_log1pCP10k.gt(aq);obs['response_module_high']=obs.response_module_log1pCP10k_mean.gt(mq)
        outcomes=obs[['nearest_published_endocardial_EC_um','Ankrd1_high','response_module_high']]
        outcomes.to_csv(folder/'association_outcomes.tsv.gz',sep='\t')
        masks_by_name=masks(obs);records=[];allpairs=[];curves=[];segment_effects=[]
        for definition in DEFINITIONS:
            d=distance[definition+'__boundary_distance_um'].to_numpy()
            segment=distance[definition+'__boundary_segment_id'].to_numpy()
            for cohort,mask in masks_by_name.items():
                choose=mask.to_numpy();finite=choose&np.isfinite(d)&(d>0)
                bins=np.array([0,15,30,60,100,150,250,500,np.inf])
                for low,high in zip(bins[:-1],bins[1:]):
                    selected=finite&(d>low)&(d<=high)
                    if selected.any():curves.append({'stage':stage,'definition':definition,'cohort':cohort,'distance_low_um':low,'distance_high_um':high,
                        'cells':int(selected.sum()),'segments':int(pd.Series(segment[selected]).nunique()),
                        'mean_Ankrd1_log1pCP10k':float(obs.loc[selected,'Ankrd1_log1pCP10k'].mean()),
                        'median_Ankrd1_log1pCP10k':float(obs.loc[selected,'Ankrd1_log1pCP10k'].median()),
                        'mean_module_log1pCP10k':float(obs.loc[selected,'response_module_log1pCP10k_mean'].mean()),
                        'Ankrd1_high_fraction':float(obs.loc[selected,'Ankrd1_high'].mean()),'module_high_fraction':float(obs.loc[selected,'response_module_high'].mean())})
                for near,far_low,far_high in CONTRASTS:
                    pairs,nn,nf,overlap=match_pairs(obs,choose,d,segment,near,far_low,far_high)
                    records.append({'stage':stage,'definition':definition,'cohort':cohort,'near_limit_um':near,'far_low_um':far_low,'far_high_um':far_high,
                        'cohort_cells':int(choose.sum()),'near_cells_before_matching':nn,'far_cells_before_matching':nf,'centroids_inside_envelope_review':overlap,**summarize_pairs(pairs)})
                    if near==30 and not pairs.empty:
                        pairs=pairs.assign(stage=stage,definition=definition,cohort=cohort)
                        allpairs.append(pairs)
                        se=pairs.assign(Ankrd1_difference=pairs.near_Ankrd1-pairs.far_Ankrd1,
                            module_difference=pairs.near_module-pairs.far_module).groupby('segment_id').agg(
                                pairs=('near_cell_id','size'),Ankrd1_difference=('Ankrd1_difference','mean'),module_difference=('module_difference','mean')).reset_index()
                        segment_effects.append(se.assign(stage=stage,definition=definition,cohort=cohort))
                # Boundary-versus-centroid sensitivity on the identical cohort
                # and vascular segments. Keep the primary contrast thresholds.
                for alternate in ['unbuffered','hull','EC_centroid']:
                    ad=distance[definition+'__'+alternate+'_distance_um'].to_numpy()
                    ass=distance[definition+'__'+alternate+'_segment_id'].to_numpy()
                    pairs,nn,nf,overlap=match_pairs(obs,choose,ad,ass,30,60,150)
                    records.append({'stage':stage,'definition':definition,'distance_metric':alternate,'cohort':cohort,
                        'near_limit_um':30,'far_low_um':60,'far_high_um':150,'cohort_cells':int(choose.sum()),
                        'near_cells_before_matching':nn,'far_cells_before_matching':nf,'centroids_inside_envelope_review':overlap,**summarize_pairs(pairs)})
        contrasts=pd.DataFrame(records);contrasts['distance_metric']=contrasts.distance_metric.fillna('boundary')
        contrasts.to_csv(folder/'local_Ankrd1_contrasts.tsv',sep='\t',index=False,na_rep='NA')
        pd.DataFrame(curves).to_csv(folder/'distance_expression_profiles.tsv',sep='\t',index=False,na_rep='NA')
        pd.concat(allpairs,ignore_index=True).to_csv(folder/'default_contrast_matched_pairs.tsv.gz',sep='\t',index=False,na_rep='NA') if allpairs else pd.DataFrame().to_csv(folder/'default_contrast_matched_pairs.tsv.gz',index=False)
        pd.concat(segment_effects,ignore_index=True).to_csv(folder/'segment_effects.tsv',sep='\t',index=False,na_rep='NA') if segment_effects else pd.DataFrame().to_csv(folder/'segment_effects.tsv',index=False)
        (folder/'association_completion.json').write_text(json.dumps({'stage':stage,'contrasts':len(contrasts),
            'stage_Ankrd1_high_cutoff_log1pCP10k':aq,'stage_module_high_cutoff_log1pCP10k':mq,
            'vessel_completion_sha256':digest(folder/'vessel_completion.json'),'lineage_completion_sha256':digest(folder/'lineage_completion.json'),
            'elapsed_seconds':time.time()-start,'no_animal_level_inference':True},indent=2)+'\n')
        show=contrasts[contrasts.definition.eq('independent_strict')&contrasts.cohort.eq('whole_reference_strict')&contrasts.near_limit_um.eq(30)&contrasts.distance_metric.eq('boundary')]
        print('ASSOCIATION_DONE',stage,show.to_dict('records'),flush=True)


if __name__=='__main__':main()
