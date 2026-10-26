"""Exact neighborhood correspondence conditioned on known ventricular CM state.

The saved Wormhole encoder is reused only on its original molecular/lineage
features. CM identity restricts comparisons; maturation and fine source labels
are held-out descriptors. Uniform correspondence is not lineage or expansion.
"""
import json,itertools
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from cardiac_subtype_niches import ROOT,OUT,UP,STAGES,read,source,save,SOURCES,SEED,FB,EC
from cardiac_story_analysis import cloud_cost,solve

def main():
    prior=ROOT/'results/cardiac_niche_dynamics';meta=read(prior/'Wormhole_heldout_clouds.tsv')
    with np.load(source(prior/'Wormhole_sample_embeddings.npz')) as b:embedding=b['heldout_embeddings']
    clouds={};observations=[];rng=np.random.default_rng(SEED)
    for stage in STAGES:
        print('SUBTYPE_TRANSPORT_PREPARE',stage,flush=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz',index_col=0)
        cm=read(OUT/(stage+'_proximity50_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id')
        candidates=meta[meta.stage.eq(stage)&meta.cell_id.isin(cm.index)].copy()
        with np.load(source(prior/stage/'nonCM_state_points.npz')) as b:points=b['points'];source_ids=b['source_ids'].astype(str)
        member=sparse.load_npz(source(prior/stage/'state_cloud_membership.npz')).tocsr()
        for state in ['v.CM','Ankrd1.CM']:
            z=candidates[candidates.cell_id.map(cm.CM_state).eq(state)].copy()
            z['near_endocardium']=z.cell_id.map(cm.endocardial_proxy_distance_um).le(100)
            for near,maximum in [(True,4),(False,8)]:
                group=z[z.near_endocardium.eq(near)];chosen=[]
                for j in rng.permutation(len(group)):
                    row=group.iloc[j];cid=row.cell_id;xy=cm.loc[cid,['x_centroid','y_centroid']].to_numpy(float)
                    if chosen and min(np.linalg.norm(xy-v) for v in chosen)<75:continue
                    chosen.append(xy);ci=int(row.name);focal=int(row.focal_index)
                    ix=member.getrow(focal).indices
                    if len(ix)>32:ix=np.random.default_rng(947+focal).choice(ix,32,replace=False)
                    clouds[ci]=points[ix];sources=f.loc[source_ids[ix]];nfb=sources.atlas_identity.eq('Fibroblast').sum();nec=sources.atlas_identity.eq('Endothelial').sum()
                    observation={'stage':stage,'CM_state':state,'distance_stratum':'≤100µm' if near else '>100µm','cell_id':cid,'cloud_index':ci,
                        'x_um':xy[0],'y_um':xy[1],'maturation_percentile_within_CM_state':cm.loc[cid,'state_maturation_percentile'],
                        'core_maturation21':cm.loc[cid,'core_Z__maturation21'],'core_metabolic':cm.loc[cid,'core_Z__metabolic'],
                        'sampled_source_profiles':len(ix),'source_FB_profiles':nfb,'source_EC_profiles':nec,'broad_context_required':False}
                    for name in FB:observation['FB__'+name]=float((sources.atlas_identity.eq('Fibroblast')&sources.fine_identity.eq(name)).sum()/nfb) if nfb else np.nan
                    for name in EC:observation['EC__'+name]=float((sources.atlas_identity.eq('Endothelial')&sources.EC_nuclear_type.eq(name)).sum()/nec) if nec else np.nan
                    observations.append(observation)
                    if len(chosen)==maximum:break
                assert len(chosen)==maximum,(stage,state,near,len(chosen))
    obs=pd.DataFrame(observations);save(obs,'CM_identity_restricted_OT_neighborhoods.tsv')
    pairs=[];plans=[];summary=[];within=[]
    for state in ['v.CM','Ankrd1.CM']:
        for stratum in ['≤100µm','>100µm']:
            groups={s:obs[obs.stage.eq(s)&obs.CM_state.eq(state)&obs.distance_stratum.eq(stratum)].reset_index(drop=True) for s in STAGES}
            for stage,g in groups.items():
                for i,j in itertools.combinations(range(len(g)),2):
                    a,b=int(g.iloc[i].cloud_index),int(g.iloc[j].cloud_index)
                    within.append({'stage':stage,'CM_state':state,'distance_stratum':stratum,'pair_kind':'Within same CM state','exact_W2_squared':cloud_cost(clouds[a],clouds[b])[0]})
            for sa,sb in zip(STAGES[:-1],STAGES[1:]):
                a,b=groups[sa],groups[sb];cost=np.zeros((len(a),len(b)));worm=np.zeros_like(cost)
                for i,x in a.iterrows():
                    for j,y in b.iterrows():
                        ai,bi=int(x.cloud_index),int(y.cloud_index);val,cs,cc=cloud_cost(clouds[ai],clouds[bi]);cost[i,j]=val;worm[i,j]=np.mean((embedding[ai]-embedding[bi])**2)
                        pairs.append({'stage_a':sa,'stage_b':sb,'CM_state':state,'distance_stratum':stratum,'cloud_a':ai,'cloud_b':bi,'exact_W2_squared':val,
                            'state_cost':cs,'lineage_cost':cc,'Wormhole_squared_distance':worm[i,j],
                            'maturation_percentile_difference':float(y.maturation_percentile_within_CM_state-x.maturation_percentile_within_CM_state)})
                value,plan=solve(cost,True);_,wp=solve(worm,True)
                for i,x in a.iterrows():
                    for j,y in b.iterrows():plans.append({'stage_a':sa,'stage_b':sb,'CM_state':state,'distance_stratum':stratum,'cloud_a':int(x.cloud_index),'cloud_b':int(y.cloud_index),
                        'exact_correspondence_mass':plan[i,j],'Wormhole_correspondence_mass':wp[i,j],'exact_cost':cost[i,j],
                        'fate_probability_or_expansion_estimate':False})
                summary.append({'stage_a':sa,'stage_b':sb,'CM_state':state,'distance_stratum':stratum,'pairs':cost.size,'Wormhole_exact_rank_correlation':spearmanr(cost.ravel(),worm.ravel()).statistic,
                    'nearest_exact_match_agreement':np.mean(np.argmin(cost,axis=1)==np.argmin(worm,axis=1)),
                    'correspondence_mass_reassigned':float(.5*abs(plan-wp).sum()),'extra_exact_cost_fraction':float((wp*cost).sum()/value-1),
                    'scope':'CM-identity-restricted, coarse-anatomical-proxy correspondence; outcome not used for distance, no cell tracking'})
                print('SUBTYPE_TRANSPORT_COMPLETE',state,stratum,sa,sb,flush=True)
    # Same-age, different-CM-state distances provide a source-neighborhood contrast.
    for stage in STAGES:
        for stratum in ['≤100µm','>100µm']:
            a=obs[obs.stage.eq(stage)&obs.CM_state.eq('v.CM')&obs.distance_stratum.eq(stratum)]
            b=obs[obs.stage.eq(stage)&obs.CM_state.eq('Ankrd1.CM')&obs.distance_stratum.eq(stratum)]
            for _,x in a.iterrows():
                for _,y in b.iterrows():within.append({'stage':stage,'CM_state':'Ankrd1.CM versus v.CM','distance_stratum':stratum,'pair_kind':'Between CM states',
                    'exact_W2_squared':cloud_cost(clouds[int(x.cloud_index)],clouds[int(y.cloud_index)])[0]})
    save(pd.DataFrame(pairs),'CM_identity_restricted_OT_pair_costs.tsv');save(pd.DataFrame(plans),'CM_identity_restricted_OT_correspondence.tsv')
    save(pd.DataFrame(summary),'CM_identity_restricted_Wormhole_assessment.tsv');save(pd.DataFrame(within),'CM_state_within_section_neighborhood_distances.tsv')
    descriptors=pd.DataFrame(pairs)[['stage_a','stage_b','CM_state','distance_stratum','cloud_a','cloud_b','maturation_percentile_difference']]
    joined=pd.DataFrame(plans).merge(descriptors,on=['stage_a','stage_b','CM_state','distance_stratum','cloud_a','cloud_b'])
    rank_rows=[]
    for key,g in joined.groupby(['stage_a','stage_b','CM_state','distance_stratum']):
        gap=g.maturation_percentile_difference.abs()
        rank_rows.append({'stage_a':key[0],'stage_b':key[1],'CM_state':key[2],'distance_stratum':key[3],
            'exact_mean_absolute_maturation_rank_gap':float((gap*g.exact_correspondence_mass).sum()),
            'Wormhole_mean_absolute_maturation_rank_gap':float((gap*g.Wormhole_correspondence_mass).sum()),
            'all_candidate_pairs_mean_absolute_rank_gap':float(gap.mean()),
            'scope':'Held-out age-relative, within-CM-state rank descriptor; no calendar-time, fate, rate or biological-replication interpretation'})
    save(pd.DataFrame(rank_rows),'CM_maturation_rank_correspondence_audit.tsv')
    (OUT/'subtype_OT_design.json').write_text(json.dumps({'seed':SEED,'observed_neighborhoods':len(obs),'focal_CM_states':['v.CM','Ankrd1.CM'],
        'same_CM_state_required_across_stages':True,'no_two_context_cluster_selection':True,'within_stage_maturation':'Held-out rank within stage×CM state, not distance or coupling input',
        'known_FB_EC_subtypes':'Source descriptors retained; uncertainty/unknown identity not filled','source_features':'Original16-PC +6-lineage Wormhole representation; no claim of subtype-label-trained encoder',
        'anatomical_proxy':'Endocardial≤100µm versus>100µm,4/8 focal profiles perCMstate/age','minimum_spacing_um_within_state_stratum':75,
        'new_training':False,'uniform_mass_is_fate_or_abundance':False},indent=2)+'\n')
    p=OUT/'source_files.json';SOURCES.update(json.loads(p.read_text()));p.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')

if __name__=='__main__':main()
