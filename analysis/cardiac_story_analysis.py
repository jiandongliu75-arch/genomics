"""Evidence tables for a cardiac-development manuscript and calibrated niche OT.

Reuses frozen QC/segmentation. Neighborhood correspondence is not cell fate.
"""
from pathlib import Path
import os,json,argparse,hashlib,itertools,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.optimize import linprog
from scipy.special import logsumexp
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/cardiac_developmental_story'
STAGES=['P0','P7','P14','P21'];SEED=26026
NICHES=ROOT/'results/cardiac_niche_dynamics'
RECON=ROOT/'results/cardiac_spatial_reconstruction'
UP=ROOT/'results/cardiac_spatial_pipeline'
SOURCES=set()

def source(path):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    SOURCES.add(str(path));return path

def read(path,**kw):return pd.read_csv(source(path),sep='\t',**kw)
def save(data,name):
    OUT.mkdir(parents=True,exist_ok=True)
    data.to_csv(OUT/name,sep='\t',index=False,float_format='%.9g')
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def RNA_fibroblasts():
    import anndata as ad
    print('STORY_RNA_FIBROBLASTS',flush=True)
    ann=read('results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz',index_col=0)
    names=read('results/fibroblast_names_and_neighborhoods/P21_fibroblast_name_key.tsv').set_index('original_cluster_id').display_name.to_dict()
    genesets={
      'Fibrillar matrix':'Col1a1 Col1a2 Col3a1 Col6a1 Col6a2 Col6a3'.split(),
      'Basement matrix':'Col4a1 Col4a2 Col4a5 Col4a6 Lama4 Hspg2'.split(),
      'Elastic matrix':'Eln Fbln1 Fbln2 Thsd4'.split(),
      'Matrix remodeling':'Fn1 Tnc Thbs1 Vcan Mmp2'.split(),
      'Cycling':'Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split(),
      'Early developmental':'Postn Ptn Tnc Col24a1'.split()}
    markers='Postn Ptn Tnc Col24a1 Col8a1 Col15a1 Hmcn2 Svep1 Fap Gfpt2 Adgrd1 Cmah Sema3c Lsamp Slit3 Eln Col4a5 Col4a6 Kirrel3 Erbb4 Unc5c'.split()
    genes=sorted(set(sum(genesets.values(),[])+markers))
    path=source('/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad')
    a=ad.read_h5ad(path,backed='r');symbols=a.var.gene_symbol.astype(str).to_numpy()
    positions=np.flatnonzero(np.isin(symbols,genes));ix=a.obs.index.get_indexer(ann.index);assert (ix>=0).all()
    # Column selection first avoids duplicating a full whole-transcriptome matrix.
    corrected=a.X[:,positions].toarray()[ix];raw=a.layers['counts_raw'][:,positions].toarray()[ix]
    collapsed=[]
    for x in [corrected,raw]:
        out=np.zeros((len(ann),len(genes)),np.float64)
        for j,pos in enumerate(positions):out[:,genes.index(symbols[pos])]+=x[:,j]
        collapsed.append(out)
    obs=a.obs.loc[ann.index];a.file.close()
    gene_rows=[];module_rows=[]
    ann['named_state']=ann.fibroblast_cluster_id.map(names).fillna(ann.fibroblast_cluster_id)
    for level in ['fibroblast_family','named_state']:
        for (age,lib,group),z in ann.groupby(['age','sample',level],observed=True):
            ii=ann.index.get_indexer(z.index)
            for view,x,total in [('corrected',collapsed[0],obs.total_counts),('raw',collapsed[1],obs.raw_total_counts)]:
                cpm=x[ii].sum(0)/max(total.iloc[ii].sum(),1)*1e6
                for g in markers:
                    j=genes.index(g);gene_rows.append({'level':level,'stage':age,'library':lib,'group':group,'nuclei':len(ii),'view':view,'gene':g,
                        'log1p_CPM':np.log1p(cpm[j]),'fraction_detected':np.mean(x[ii,j]>=.001)})
                for label,panel in genesets.items():
                    js=[genes.index(g) for g in panel]
                    module_rows.append({'level':level,'stage':age,'library':lib,'group':group,'nuclei':len(ii),'view':view,'program':label,
                        'mean_log1p_CPM':np.log1p(cpm[js]).mean(),'genes':';'.join(panel),'new_biological_replicates':0})
    save(pd.DataFrame(gene_rows),'fibroblast_library_marker_profiles.tsv')
    save(pd.DataFrame(module_rows),'fibroblast_library_matrix_programs.tsv')
    previous=read('results/fibroblast_states/RNA_adjacent_stage_candidate_matches.tsv')
    key=ann.groupby('fibroblast_cluster_id').agg(family=('fibroblast_family','first'),nuclei=('age','size')).to_dict('index')
    previous=previous[previous.from_cluster.isin(key)&previous.to_cluster.isin(key)].copy()
    for role in ['from','to']:
        previous[role+'_family']=previous[role+'_cluster'].map(lambda k:key[k]['family'])
        previous[role+'_name']=previous[role+'_cluster'].map(names).fillna(previous[role+'_cluster'])
        previous[role+'_nuclei']=previous[role+'_cluster'].map(lambda k:key[k]['nuclei'])
    previous['interpretation']='Frozen gene-withheld RNA signature resemblance, with current reviewed family names; candidate state continuity, not lineage ancestry, transport probability or conversion rate'
    save(previous,'reviewed_fibroblast_candidate_state_continuities.tsv')
    save(pd.DataFrame([{'program':k,'genes':';'.join(v),'count_denominator':'Full RNA pseudobulk molecules','function_measured':False} for k,v in genesets.items()]),'fibroblast_RNA_program_definitions.tsv')

def solve(cost,return_plan=False):
    n,m=cost.shape
    A=sparse.vstack([sparse.kron(sparse.eye(n),np.ones((1,m))),sparse.kron(np.ones((1,n)),sparse.eye(m))]).tocsr()
    r=linprog(cost.ravel(),A_eq=A,b_eq=np.r_[np.full(n,1/n),np.full(m,1/m)],bounds=(0,None),method='highs')
    assert r.success,r.message
    plan=r.x.reshape(n,m)
    assert max(abs(plan.sum(1)-1/n).max(),abs(plan.sum(0)-1/m).max())<1e-7
    return (r.fun,plan) if return_plan else r.fun

def cloud_cost(x,y,state_weight=.5,maximum=32):
    # Original features have half state / half identity squared-norm weighting.
    # Deterministic truncation tests sensitivity without inventing molecules.
    x=x[:maximum];y=y[:maximum]
    cs=cdist(x[:,:16],y[:,:16],'sqeuclidean')*(state_weight/.5)
    cc=cdist(x[:,16:],y[:,16:],'sqeuclidean')*((1-state_weight)/.5)
    value,plan=solve(cs+cc,True)
    return value,float((plan*cs).sum()),float((plan*cc).sum())

def entropic_plan(cost):
    eps=max(np.median(cost)*.05,1e-5);n,m=cost.shape
    kernel=-cost/eps;u=np.zeros(n);v=np.zeros(m)
    for _ in range(5000):
        u=-np.log(n)-logsumexp(kernel+v[None,:],axis=1)
        v=-np.log(m)-logsumexp(kernel+u[:,None],axis=0)
        plan=np.exp(kernel+u[:,None]+v[None,:])
        if max(abs(plan.sum(1)-1/n).max(),abs(plan.sum(0)-1/m).max())<1e-9:break
    assert max(abs(plan.sum(1)-1/n).max(),abs(plan.sum(0)-1/m).max())<1e-6
    return plan

def neighborhood_OT():
    print('STORY_ANATOMY_STRATIFIED_NEIGHBORHOOD_OT',flush=True)
    meta=read(NICHES/'Wormhole_heldout_clouds.tsv')
    with np.load(source(NICHES/'Wormhole_sample_embeddings.npz')) as b:emb=b['heldout_embeddings']
    clouds={};observations=[];rng=np.random.default_rng(SEED)
    for stage in STAGES:
        f=read(UP/stage/'atlas_annotations.tsv.gz',index_col=0)
        c=read(RECON/stage/'CM_outcomes_and_anatomy.tsv.gz',index_col=0)
        niche=read(UP/stage/'niche_annotations.tsv.gz',index_col=0)
        candidates=meta[meta.stage.eq(stage)].copy()
        ids=candidates.cell_id
        z=f.loc[ids].join(c[['reconstructed_Figure4e_CM']],how='left').join(niche[['niche']])
        z['cloud_index']=candidates.index.to_numpy()
        valid=z.reconstructed_Figure4e_CM.fillna(False)&z.niche.str.startswith('N2:')&z.author_atlas_cell_state.str.contains('CM',na=False)
        z=z[valid].copy();z['distance_stratum']=np.where(z.endocardial_proxy_distance_um.le(100),'Endocardial-proximal ≤100µm','More distant >100µm')
        with np.load(source(NICHES/stage/'nonCM_state_points.npz')) as b:points=b['points'];source_ids=b['source_ids'].astype(str)
        membership=sparse.load_npz(source(NICHES/stage/'state_cloud_membership.npz')).tocsr()
        selected=[]
        for stratum,group in z.groupby('distance_stratum'):
            chosen=[]
            for j in rng.permutation(len(group)):
                row=group.iloc[j];xy=row[['x_centroid','y_centroid']].to_numpy(float)
                if chosen and min(np.linalg.norm(xy-k) for k in chosen)<100:continue
                chosen.append(xy);selected.append((stratum,group.index[j]))
                if len(chosen)==12:break
            assert len(chosen)==12,(stage,stratum,len(chosen))
        for stratum,cid in selected:
            row=z.loc[cid];ci=int(row.cloud_index);focal=int(meta.iloc[ci].focal_index)
            ix=membership.indices[membership.indptr[focal]:membership.indptr[focal+1]]
            if len(ix)>32:ix=np.random.default_rng(947+focal).choice(ix,32,replace=False)
            x=points[ix];clouds[ci]=x
            obs={'stage':stage,'cell_id':cid,'cloud_index':ci,'distance_stratum':stratum,'x_um':row.x_centroid,'y_um':row.y_centroid,
                'endocardial_proxy_distance_um':row.endocardial_proxy_distance_um,'sampled_source_profiles':len(x),'eligible_source_profiles':np.diff(membership.indptr)[focal],
                'niche':'CM/EC-enriched','spatial_block_400um':f'{int(row.x_centroid//400)}_{int(row.y_centroid//400)}',
                'CM_core_logmean_maturation21':c.loc[cid,'core_logmean__maturation21'],'CM_core_metabolic':c.loc[cid,'core_Z__metabolic']}
            for key in ['FB_basement_matrix','FB_remodeling','EC_cycling','source_log_capture']:
                source_values=f.loc[source_ids[ix],key].to_numpy(float)
                obs['source_mean_'+key]=float(np.nanmean(source_values)) if np.isfinite(source_values).any() else np.nan
            for j,label in enumerate(['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']):obs['source_fraction_'+label]=float(x[:,16+j].mean()*np.sqrt(2))
            observations.append(obs)
    obs=pd.DataFrame(observations);save(obs,'OT_selected_observed_neighborhoods.tsv')
    rows=[];plans=[];assess=[];baseline=[]
    for stratum in obs.distance_stratum.unique():
        groups={stage:obs[obs.stage.eq(stage)&obs.distance_stratum.eq(stratum)].reset_index(drop=True) for stage in STAGES}
        for stage,g in groups.items():
            for i,j in itertools.combinations(range(len(g)),2):
                ai,bi=int(g.iloc[i].cloud_index),int(g.iloc[j].cloud_index)
                value,cs,cc=cloud_cost(clouds[ai],clouds[bi]);baseline.append({'stage':stage,'distance_stratum':stratum,'cloud_a':ai,'cloud_b':bi,'exact_W2_squared':value,'state_cost':cs,'identity_cost':cc})
        for sa,sb in zip(STAGES[:-1],STAGES[1:]):
            ga,gb=groups[sa],groups[sb];matrices={k:np.empty((len(ga),len(gb))) for k in ['primary','state25','state75','subsample16','Wormhole','mean_features']}
            for i,a in ga.iterrows():
                for j,b in gb.iterrows():
                    ai,bi=int(a.cloud_index),int(b.cloud_index);x,y=clouds[ai],clouds[bi]
                    value,cs,cc=cloud_cost(x,y);matrices['primary'][i,j]=value
                    matrices['state25'][i,j]=cloud_cost(x,y,.25)[0];matrices['state75'][i,j]=cloud_cost(x,y,.75)[0]
                    matrices['subsample16'][i,j]=cloud_cost(x,y,maximum=16)[0]
                    matrices['Wormhole'][i,j]=np.mean((emb[ai]-emb[bi])**2)
                    matrices['mean_features'][i,j]=np.sum((x.mean(0)-y.mean(0))**2)
                    rows.append({'stage_a':sa,'stage_b':sb,'distance_stratum':stratum,'cloud_a':ai,'cloud_b':bi,'exact_W2_squared':value,
                        'state_cost_on_optimal_plan':cs,'identity_cost_on_optimal_plan':cc,
                        'lineage_L1_difference':np.abs(x[:,16:].mean(0)-y[:,16:].mean(0)).sum()*np.sqrt(2),
                        **{key:matrices[key][i,j] for key in matrices if key!='primary'}})
            value,plan=solve(matrices['primary'],True);soft=entropic_plan(matrices['primary'])
            for key,cost in matrices.items():
                _,p=solve(cost,True)
                assess.append({'stage_a':sa,'stage_b':sb,'distance_stratum':stratum,'metric':key,
                    'exact_cost_rank_correlation':spearmanr(matrices['primary'].ravel(),cost.ravel()).statistic,
                    'nearest_exact_neighbor_agreement':np.mean(np.argmin(matrices['primary'],axis=1)==np.argmin(cost,axis=1)),
                    'transport_plan_total_variation_from_primary':.5*np.abs(p-plan).sum(),
                    'primary_exact_cost_inflation':float((p*matrices['primary']).sum()/value-1),
                    'entropy_soft_plan_total_variation_from_primary':float(.5*abs(soft-plan).sum())})
            for i,a in ga.iterrows():
                for j,b in gb.iterrows():
                    plans.append({'stage_a':sa,'stage_b':sb,'distance_stratum':stratum,'cloud_a':int(a.cloud_index),'cloud_b':int(b.cloud_index),
                        'balanced_correspondence_mass':plan[i,j],'entropic_correspondence_mass':soft[i,j],'exact_cloud_distance':matrices['primary'][i,j],
                        'uniform_mass_is_an_abundance_estimate':False,'interpretation':'Restricted cross-sectional neighborhood correspondence; not tracked cells, lineage, actual temporal flow or transition probability'})
            print('OT_COMPLETED',sa,sb,stratum,flush=True)
    pairs=pd.DataFrame(rows);save(pairs,'anatomy_stratified_neighborhood_pair_costs.tsv')
    save(pd.DataFrame(plans),'adjacent_stage_neighborhood_correspondence_plans.tsv')
    save(pd.DataFrame(assess),'neighborhood_transport_metric_sensitivity.tsv')
    save(pd.DataFrame(baseline),'within_section_neighborhood_distance_baseline.tsv')
    summaries=[]
    for (sa,sb,stratum),g in pairs.groupby(['stage_a','stage_b','distance_stratum']):
        ss=g[g.lineage_L1_difference.le(.3)]
        if len(ss):summaries.append({'stage_a':sa,'stage_b':sb,'distance_stratum':stratum,'comparison':'Composition-similar pairs L1≤0.3',
            'pairs':len(ss),'median_exact_W2_squared':ss.exact_W2_squared.median(),'median_identity_cost':ss.identity_cost_on_optimal_plan.median(),
            'median_state_cost':ss.state_cost_on_optimal_plan.median(),'median_state_fraction_of_cost':np.median(ss.state_cost_on_optimal_plan/ss.exact_W2_squared),
            'scope':'Representation-dependent descriptive distance decomposition; not biological effect size or proof of state conversion'})
    save(pd.DataFrame(summaries),'composition_similar_neighborhood_cost_decomposition.tsv')
    (OUT/'OT_design.json').write_text(json.dumps({'seed':SEED,'observed_neighborhoods':len(obs),'per_stage':24,'per_distance_stratum_per_stage':12,
        'reserved_neighborhoods_only':True,'focal_cohort':'Strict ventricular and public-CM-compatible, independent ventricular selector, existing CM/EC-enriched context',
        'anatomical_control':'Same endocardial proxy distance category (≤100 versus >100 µm); no anatomical registration or verified LV/RV/compact/trabecular match',
        'selection_spacing':'At least100 µm within each stage×distance category; overlapping/adjacent sampling regions are not independent hearts',
        'state_features':'Saved16-PC molecular +6-lineage representation; CM maturation outcomes excluded from inputs and normalization',
        'two_transport_levels':'Exact non-CM cloud-to-cloud W2², then uniform-marginal neighborhood correspondence within strata',
        'sensitivities':['State/identity weights25/75,50/50,75/25','At most16 versus32 source profiles/cloud','Saved Wormhole distance','Mean features','Entropic versus unregularized correspondence'],
        'new_encoder_training':False,'no_fates_or_rates':True,'new_biological_replication':0,'stage_and_section_confounded':True},indent=2)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['RNA','OT','all'],default='all');a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    start=time.time()
    if a.step in ['RNA','all']:RNA_fibroblasts()
    if a.step in ['OT','all']:neighborhood_OT()
    prev=OUT/'analysis_source_files.json'
    if prev.exists():SOURCES.update(json.loads(prev.read_text()))
    prev.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')
    (OUT/f'execution_{a.step}.json').write_text(json.dumps({'completed':True,'seconds':time.time()-start,'step':a.step,'seed':SEED},indent=2)+'\n')

if __name__=='__main__':main()
