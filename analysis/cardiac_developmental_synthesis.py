"""Cross-stage diagnostic for the already trained neighborhood metric.

The biological synthesis is primary. This optional diagnostic does not reconstruct
lineages, vessel perfusion, actual contact or causal developmental transitions.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-synthesis-mpl')
from pathlib import Path
import argparse,json,hashlib,itertools
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.optimize import linprog
from scipy.stats import spearmanr
from fit_wormhole_cardiac_niches import cloud_subset
from revisit_vascular_coordination import ROOT,CFG

SOURCE=ROOT/'results/cardiac_niche_dynamics'
OUT=ROOT/'results/cardiac_developmental_synthesis'

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def benchmark():
    import jax
    import jax.numpy as jnp
    from wassersteinwormhole import utils_OT
    calc=jax.jit(lambda x,a,y,b:utils_OT.S2([x,a],[y,b],eps=.1,lse_mode=True,num_iter=100,ot_scale=1.))
    meta=pd.read_csv(SOURCE/'Wormhole_heldout_clouds.tsv',sep='\t')
    with np.load(SOURCE/'Wormhole_sample_embeddings.npz') as b:emb=b['heldout_embeddings']
    stores={};features={};rng=np.random.default_rng(3841)
    for stage in CFG['stages']:
        with np.load(SOURCE/stage/'nonCM_state_points.npz') as b:points=b['points']
        member=sparse.load_npz(SOURCE/stage/'state_cloud_membership.npz').tocsr()
        for i in np.flatnonzero(meta.stage.eq(stage)):
            pc,_=cloud_subset(points,member,int(meta.iloc[i].focal_index),CFG['seed']);stores[i]=pc;features[i]=pc.mean(0)
    records=[]
    for sa,sb in itertools.combinations(CFG['stages'],2):
        print('CROSS_STAGE',sa,sb,flush=True)
        ia=np.flatnonzero(meta.stage.eq(sa));ib=np.flatnonzero(meta.stage.eq(sb))
        mean_b=np.array([features[i] for i in ib]);composition_b=mean_b[:,-6:]*np.sqrt(2)
        for mode in ['random_cross_stage','composition_matched_cross_stage']:
            aa=rng.choice(ia,64,replace=False)
            for ai in aa:
                if mode=='random_cross_stage':bi=int(rng.choice(ib))
                else:
                    l1=np.abs(composition_b-features[ai][-6:]*np.sqrt(2)).sum(1)
                    bi=int(ib[rng.choice(np.argsort(l1)[:5])])
                x,y=stores[ai],stores[bi];n=CFG['wormhole_maximum_cloud_points'];d=x.shape[1]
                xx=np.zeros((n,d),np.float32);yy=np.zeros((n,d),np.float32);wa=np.zeros(n,np.float32);wb=np.zeros(n,np.float32)
                xx[:len(x)]=x;yy[:len(y)]=y;wa[:len(x)]=1/len(x);wb[:len(y)]=1/len(y)
                sd=float(calc(jnp.array(xx),jnp.array(wa),jnp.array(yy),jnp.array(wb)))
                na,nb=len(x),len(y);A=sparse.vstack([sparse.kron(sparse.eye(na),np.ones((1,nb))),sparse.kron(np.ones((1,na)),sparse.eye(nb))]).tocsr()
                solved=linprog(cdist(x,y,'sqeuclidean').ravel(),A_eq=A,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs');assert solved.success
                records.append({'stage_a':sa,'stage_b':sb,'pair_selection':mode,'CM_a':meta.iloc[ai].cell_id,'CM_b':meta.iloc[bi].cell_id,
                    'cloud_index_a':ai,'cloud_index_b':bi,'Sinkhorn_divergence':sd,'exact_W2_squared':solved.fun,
                    'Wormhole_mean_squared_distance':float(np.mean((emb[ai]-emb[bi])**2)),
                    'mean_state_squared_distance':float(np.sum((features[ai]-features[bi])**2)),
                    'lineage_composition_L1_difference':float(np.abs(features[ai][-6:]-features[bi][-6:]).sum()*np.sqrt(2)),
                    'scope':'Reserved neighborhoods from separate sections; composition matching has no anatomical registration, biological replication or causal interpretation'})
    result=pd.DataFrame(records);result.to_csv(OUT/'cross_stage_distance_pairs.tsv',sep='\t',index=False,float_format='%.9g')
    rows=[]
    for (sa,sb,mode),z in result.groupby(['stage_a','stage_b','pair_selection'],sort=False):
        rows.append({'stage_a':sa,'stage_b':sb,'pair_selection':mode,'pairs':len(z),
            'Wormhole_vs_Sinkhorn_Spearman':spearmanr(z.Sinkhorn_divergence,z.Wormhole_mean_squared_distance).statistic,
            'mean_vs_Sinkhorn_Spearman':spearmanr(z.Sinkhorn_divergence,z.mean_state_squared_distance).statistic,
            'Wormhole_vs_exact_W2_Spearman':spearmanr(z.exact_W2_squared,z.Wormhole_mean_squared_distance).statistic,
            'median_composition_L1_difference':z.lineage_composition_L1_difference.median(),
            'interpretation':'Metric calibration only; not cross-stage niche identity or developmental mechanism validation'})
    pd.DataFrame(rows).to_csv(OUT/'cross_stage_metric_assessment.tsv',sep='\t',index=False,float_format='%.9g')
    source=[SOURCE/'Wormhole_sample_embeddings.npz',SOURCE/'Wormhole_heldout_clouds.tsv',ROOT/'analysis/cardiac_niche_dynamics.json']
    for stage in CFG['stages']:source.extend([SOURCE/stage/'nonCM_state_points.npz',SOURCE/stage/'state_cloud_membership.npz'])
    pd.DataFrame([{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in source]).to_csv(OUT/'cross_stage_source_checksums.tsv',sep='\t',index=False)
    (OUT/'cross_stage_execution.json').write_text(json.dumps({'trained_encoder_reused':True,'additional_training':False,'reserved_clouds':len(meta),
        'age_pairs':6,'pairs_per_age_pair_per_selection':64,'selections':['random','top5_lineage_composition_match'],
        'composition_matching_uses_CM_maturation':False,'molecular_ground_metric':'Original outcome-withheld non-CM state PCA + lineage features',
        'absolute_cell_abundance_in_metric':False,'within_cloud_spatial_arrangement_in_metric':False,
        'anatomical_compartments_matched':False,'metric_evaluation_is_biological_validation':False,'seed':3841},indent=2)+'\n')

def validate():
    z=pd.read_csv(OUT/'cross_stage_distance_pairs.tsv',sep='\t');m=pd.read_csv(OUT/'cross_stage_metric_assessment.tsv',sep='\t')
    assert len(z)==768 and len(m)==12 and np.isfinite(z.select_dtypes('number')).all().all()
    meta=pd.read_csv(SOURCE/'Wormhole_heldout_clouds.tsv',sep='\t')
    with np.load(SOURCE/'Wormhole_sample_embeddings.npz') as b:emb=b['heldout_embeddings']
    for r in z.sample(12,random_state=5).itertuples():
        assert meta.iloc[r.cloud_index_a].stage==r.stage_a and meta.iloc[r.cloud_index_b].stage==r.stage_b
        np.testing.assert_allclose(np.mean((emb[r.cloud_index_a]-emb[r.cloud_index_b])**2),r.Wormhole_mean_squared_distance,rtol=3e-6)
    for r in m.itertuples():
        d=z[z.stage_a.eq(r.stage_a)&z.stage_b.eq(r.stage_b)&z.pair_selection.eq(r.pair_selection)]
        assert len(d)==64
        np.testing.assert_allclose(spearmanr(d.Sinkhorn_divergence,d.Wormhole_mean_squared_distance).statistic,r.Wormhole_vs_Sinkhorn_Spearman,rtol=3e-6)
    # Check one exact cost per stage contrast using independently constructed
    # dense marginal constraints instead of the production sparse Kronecker form.
    for (sa,sb),d in z.groupby(['stage_a','stage_b']):
        r=d.iloc[0];clouds=[]
        for stage,index in [(sa,int(r.cloud_index_a)),(sb,int(r.cloud_index_b))]:
            with np.load(SOURCE/stage/'nonCM_state_points.npz') as b:points=b['points']
            member=sparse.load_npz(SOURCE/stage/'state_cloud_membership.npz').tocsr()
            pc,_=cloud_subset(points,member,int(meta.iloc[index].focal_index),CFG['seed']);clouds.append(pc)
        x,y=clouds;na,nb=len(x),len(y);A=np.zeros((na+nb,na*nb))
        for i in range(na):A[i,i*nb:(i+1)*nb]=1
        for j in range(nb):A[na+j,j::nb]=1
        solution=linprog(cdist(x,y,'sqeuclidean').ravel(),A_eq=A,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs')
        assert solution.success;np.testing.assert_allclose(solution.fun,r.exact_W2_squared,rtol=3e-6)
        selected=d[d.pair_selection.eq('composition_matched_cross_stage')];random=d[d.pair_selection.eq('random_cross_stage')]
        assert selected.lineage_composition_L1_difference.median()<random.lineage_composition_L1_difference.median()
    for r in pd.read_csv(OUT/'cross_stage_source_checksums.tsv',sep='\t').itertuples():assert sha(ROOT/r.path)==r.sha256
    (OUT/'cross_stage_validation.json').write_text(json.dumps({'status':'PASS','cross_age_pairs':len(z),'stage_pair_selection_groups':len(m),
        'heldout_metadata_checked':True,'saved_embedding_distance_independently_recomputed':True,'rank_metrics_recomputed':True,'source_hashes_verified':True,
        'exact_OT_independent_dense_LP_checks':6,'composition_matching_reduces_lineage_difference':True,
        'scope':'Computational consistency; no biological validation of developmental correspondence'},indent=2)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['benchmark','validate','all'],default='all');a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if a.step in ['benchmark','all']:benchmark()
    if a.step in ['validate','all']:validate()

if __name__=='__main__':main()
