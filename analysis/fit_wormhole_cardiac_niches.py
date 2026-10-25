"""Richer, outcome-withheld non-CM niche distributions and a Wormhole benchmark.

This models cross-sectional neighborhood states. It cannot track cells, infer
endocardial-to-coronary ancestry, or restore missing developmental timepoints.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-dynamics-mpl')
from revisit_vascular_coordination import ROOT,CFG,OUT,STAGES,read,save,dump
from cardiac_spatial_pipeline import MAT,NC,MARKERS
from pathlib import Path
import argparse,json,gc,time
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.optimize import linprog
from scipy.stats import spearmanr,pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.cluster import KMeans

UP=ROOT/CFG['upstream']

def cloud_subset(points,membership,row,seed):
    ix=membership.indices[membership.indptr[row]:membership.indptr[row+1]]
    if len(ix)>CFG['wormhole_maximum_cloud_points']:
        ix=np.random.default_rng(seed+int(row)).choice(ix,CFG['wormhole_maximum_cloud_points'],replace=False)
    return points[ix],ix

def prepare():
    fit=[];stage_info=[];coverage=[]
    excluded=set(MAT.gene)|set(MARKERS['Cardiomyocyte'])|{'Ryr2_iso2','Ryr2_iso3','Myom1_E_1','Myom1_E_2','Myom1_M','Ankrd1','Atf3','Creb5','Tlr4'}
    for stage in STAGES:
        print('STATE_PCA_PREPARE',stage,flush=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz');cm=f.loc[f.ventricular_CM_strict]
        mask=f.atlas_identity.isin(NC)&f.program_source_eligible
        ids=f.index[mask];positions=np.flatnonzero(mask.to_numpy())
        b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        genes=b.var_names.to_numpy(dtype=str);use=~np.isin(genes,list(excluded));used=genes[use]
        assert not set(used)&set(MAT.gene)
        x=b.layers['high_confidence'][positions][:,use].toarray().astype(np.float32)
        depth=x.sum(1);x=np.log1p(x*1e4/np.maximum(depth,1)[:,None])
        membership=sparse.load_npz(UP/stage/'proximity50_graph.npz').tocsr()[f.index.get_indexer(cm.index)][:,positions].tocsr()
        eligible=(np.diff(membership.indptr)>=3)&cm.CM_outcome_excluded_depth.ge(20).to_numpy()
        xy=cm[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/CFG['block_um']).astype(int)
        folds=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
        test=np.flatnonzero(eligible&(folds==0));train=np.flatnonzero(eligible&(folds!=0))
        train=train[cKDTree(xy[test]).query(xy[train])[0]>=CFG['holdout_buffer_um']]
        # No PCA fitting source may occur in any reserved test CM cloud.
        test_source=np.unique(membership[test].indices);source_training=np.ones(len(ids),bool);source_training[test_source]=False
        rng=np.random.default_rng(CFG['seed']);available=np.flatnonzero(source_training)
        take=rng.choice(available,min(CFG['pca_training_nonCM_per_stage'],len(available)),replace=False)
        fit.append(x[take]);stage_info.append((stage,x,used,ids,positions,cm,train,test,take))
        sparse.save_npz(OUT/stage/'state_cloud_membership.npz',membership)
        coverage.append({'stage':stage,'source_nonCM_profiles':len(ids),'focal_CM_profiles':len(cm),'eligible_clouds':int(eligible.sum()),
            'training_region_CMs':len(train),'heldout_region_CMs':len(test),'PCA_training_sources':len(take),'PCA_sources_in_heldout_clouds':int(np.isin(take,test_source).sum()),
            'minimum_focal_train_test_distance_um':cKDTree(xy[test]).query(xy[train])[0].min()})
        b.file.close();del b;gc.collect()
    train_matrix=np.vstack(fit);mu=train_matrix.mean(0);sd=np.maximum(train_matrix.std(0),.25)
    scaled=np.clip((train_matrix-mu)/sd,-5,5)
    pca=PCA(n_components=CFG['pca_components'],random_state=CFG['seed']).fit(scaled)
    pcsd=np.maximum(pca.transform(scaled).std(0),.25)
    np.savez_compressed(OUT/'nonCM_state_PCA_model.npz',genes=stage_info[0][2],excluded_genes=np.array(sorted(excluded)),gene_mean=mu,gene_SD=sd,
        pca_components=pca.components_,pca_mean=pca.mean_,PC_SD=pcsd,explained_variance_ratio=pca.explained_variance_ratio_)
    for stage,x,genes,ids,positions,cm,train,test,take in stage_info:
        assert np.array_equal(genes,stage_info[0][2])
        pc=np.clip(pca.transform(np.clip((x-mu)/sd,-5,5))/pcsd,-4,4)
        f=read(UP/stage/'atlas_annotations.tsv.gz')
        one=np.column_stack([f.loc[ids,'atlas_identity'].eq(l).to_numpy(float) for l in NC])
        # Equal state/composition blocks define the empirical OT ground metric.
        points=np.column_stack([pc/np.sqrt(2*pc.shape[1]),one/np.sqrt(2)]).astype(np.float32)
        rng=np.random.default_rng(CFG['seed']);train_sample=rng.choice(train,min(len(train),CFG['wormhole_training_clouds_per_stage']),replace=False)
        test_sample=rng.choice(test,min(len(test),CFG['wormhole_assessment_clouds_per_stage']),replace=False)
        np.savez_compressed(OUT/stage/'nonCM_state_points.npz',points=points,source_ids=ids.to_numpy(dtype=str),
            focal_ids=cm.index.to_numpy(dtype=str),PCA_training_source_indices=take,source_lineages=f.loc[ids,'atlas_identity'].to_numpy(dtype=str),
            training_focal_indices=train,heldout_focal_indices=test,training_sample=train_sample,heldout_sample=test_sample)
    save(pd.DataFrame(coverage),'state_embedding_cohort_coverage.tsv')
    dump('state_feature_definition.json',{'genes':stage_info[0][2].tolist(),'gene_count':len(stage_info[0][2]),'PCs':CFG['pca_components'],
        'focal_CM_RNA_used':False,'maturation_genes_used':False,'composition_and_state_blocks':'Equal squared-norm weighting',
        'PCA_fit':'Stage-balanced sources excluding every source in reserved heldout CM clouds',
        'technical_limits':'Core RNA normalization; residual depth/segmentation effects possible; one section per age confounds age with section',
        'time_interpretation':CFG['interpretation']})

def training_data():
    clouds=[];testclouds=[];trainmeta=[];testmeta=[]
    for stage in STAGES:
        with np.load(OUT/stage/'nonCM_state_points.npz') as b:
            points=b['points'];ids=b['focal_ids'].astype(str);train=b['training_sample'];test=b['heldout_sample']
        membership=sparse.load_npz(OUT/stage/'state_cloud_membership.npz').tocsr()
        for indices,target,meta in [(train,clouds,trainmeta),(test,testclouds,testmeta)]:
            for i in indices:
                pc,_=cloud_subset(points,membership,int(i),CFG['seed']);target.append(pc)
                meta.append({'stage':stage,'cell_id':ids[i],'focal_index':int(i)})
    return clouds,testclouds,pd.DataFrame(trainmeta),pd.DataFrame(testmeta)

def fit_wormhole():
    import jax
    from wassersteinwormhole import Wormhole
    from wassersteinwormhole.DefaultConfig import DefaultConfig
    from flax.serialization import to_bytes
    train,test,tm,qm=training_data()
    config=DefaultConfig().replace(emb_dim=CFG['wormhole_embedding_dim'],mlp_dim=CFG['wormhole_hidden_dim'],num_layers=CFG['wormhole_layers'],
        num_heads=CFG['wormhole_heads'],eps_enc=.1,eps_dec=.1,num_sinkhorn_iter=100,num_particles_output=16,scale_ot=False)
    model=Wormhole(point_clouds=train,point_clouds_test=test,config=config)
    start=time.time();model.train(training_steps=CFG['wormhole_training_steps'],batch_size=CFG['wormhole_batch_size'],verbose=100,
        key=jax.random.key(CFG['seed']))
    train_embedding=np.asarray(model.encode(model.point_clouds,model.weights,max_batch=128))
    test_embedding=np.asarray(model.encode(model.point_clouds_test,model.weights_test,max_batch=128))
    assert np.isfinite(train_embedding).all() and np.isfinite(test_embedding).all()
    np.savez_compressed(OUT/'Wormhole_sample_embeddings.npz',training_embeddings=train_embedding,heldout_embeddings=test_embedding)
    save(tm,'Wormhole_training_clouds.tsv');save(qm,'Wormhole_heldout_clouds.tsv')
    (OUT/'Wormhole_parameters.msgpack').write_bytes(to_bytes(model.params))
    save(pd.DataFrame({'encoder_loss':np.asarray(model.enc_loss_curve),'decoder_loss':np.asarray(model.dec_loss_curve)}),'Wormhole_training_loss.tsv')
    # Encode all eligible held-out/training regional niches for shared-cohort
    # comparisons, not just the small clouds used by the transformer optimizer.
    for stage in STAGES:
        print('WORMHOLE_ENCODE',stage,flush=True)
        with np.load(OUT/stage/'nonCM_state_points.npz') as b:points=b['points'];ids=b['focal_ids'].astype(str);trainids=b['training_focal_indices'];testids=b['heldout_focal_indices']
        membership=sparse.load_npz(OUT/stage/'state_cloud_membership.npz').tocsr()
        selected=np.r_[trainids,testids];output=[];means=[];cloudsizes=[]
        for start in range(0,len(selected),128):
            batch=[]
            for i in selected[start:start+128]:
                pc,_=cloud_subset(points,membership,int(i),CFG['seed']);batch.append(pc);means.append(pc.mean(0));cloudsizes.append(len(pc))
            output.append(np.asarray(model.encode(batch,max_batch=128)))
        np.savez_compressed(OUT/stage/'state_niche_embeddings.npz',cell_ids=ids[selected],focal_indices=selected,
            heldout=np.r_[np.zeros(len(trainids),bool),np.ones(len(testids),bool)],embedding=np.vstack(output),mean_features=np.array(means),cloud_points=np.array(cloudsizes))
    # Holdout metric benchmark is performed before using any CM phenotype.
    dump('Wormhole_execution.json',{'executed':True,'package':'wassersteinwormhole0.3.7','JAX':jax.__version__,'devices':[str(x) for x in jax.devices()],
        'training_steps':CFG['wormhole_training_steps'],'training_clouds':len(train),'heldout_clouds':len(test),'seconds':time.time()-start,
        'embedding_dim':CFG['wormhole_embedding_dim'],'architecture':'2 transformer blocks,2 heads,hidden64; encoder/decoder Sinkhorn divergence epsilon0.1,100 iterations;16 decoder points',
        'paper_default_reproduction':False,'purpose':'CPU-constrained, region-held-out metric benchmark; adoption depends on reserved distance and prediction tests'})
    benchmark_distances(model,tm,qm,train_embedding,test_embedding)

def benchmark_distances(model,tm,qm,te,qe):
    from wassersteinwormhole import utils_OT
    import jax
    import jax.numpy as jnp
    # S2 takes two [points, weights] pairs, not a points pair and a weights
    # pair. Fixed padding avoids compiling a solver for each cloud size.
    sinkhorn=jax.jit(lambda xa,wa,xb,wb: utils_OT.S2([xa,wa],[xb,wb],eps=.1,lse_mode=True,num_iter=100,ot_scale=1.))
    rows=[];rng=np.random.default_rng(CFG['seed'])
    for stage in STAGES:
        m=np.flatnonzero(qm.stage.eq(stage));pairs=np.array([rng.choice(m,2,replace=False) for _ in range(64)])
        with np.load(OUT/stage/'nonCM_state_points.npz') as b:points=b['points']
        membership=sparse.load_npz(OUT/stage/'state_cloud_membership.npz').tocsr()
        for ai,bi in pairs:
            arow=int(qm.iloc[ai].focal_index);brow=int(qm.iloc[bi].focal_index)
            xa,_=cloud_subset(points,membership,arow,CFG['seed']);xb,_=cloud_subset(points,membership,brow,CFG['seed'])
            n=CFG['wormhole_maximum_cloud_points'];padded=np.zeros((2,n,xa.shape[1]),np.float32);weights=np.zeros((2,n),np.float32)
            padded[0,:len(xa)]=xa;padded[1,:len(xb)]=xb;weights[0,:len(xa)]=1/len(xa);weights[1,:len(xb)]=1/len(xb)
            sink=float(sinkhorn(jnp.asarray(padded[0]),jnp.asarray(weights[0]),jnp.asarray(padded[1]),jnp.asarray(weights[1])))
            exact=np.nan
            # Full small LP gives an entropic-regularization-independent check.
            na,nb=len(xa),len(xb);constraints=sparse.vstack([sparse.kron(sparse.eye(na),np.ones((1,nb))),sparse.kron(np.ones((1,na)),sparse.eye(nb))]).tocsr()
            solved=linprog(cdist(xa,xb,'sqeuclidean').ravel(),A_eq=constraints,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs')
            assert solved.success;exact=solved.fun
            rows.append({'stage':stage,'cell_id_a':qm.iloc[ai].cell_id,'cell_id_b':qm.iloc[bi].cell_id,'Sinkhorn_divergence':sink,'exact_W2_squared':exact,
                # The package trains mean squared coordinate distance.
                'Wormhole_mean_squared_embedding_distance':float(np.mean((qe[ai]-qe[bi])**2)),
                'mean_squared_feature_distance':float(np.sum((xa.mean(0)-xb.mean(0))**2))})
    save(pd.DataFrame(rows),'Wormhole_reserved_distance_benchmark.tsv')

def resume_benchmark():
    """Assess an already trained/saved model without repeating optimization."""
    import jax
    tm=read(OUT/'Wormhole_training_clouds.tsv',index_col=None);qm=read(OUT/'Wormhole_heldout_clouds.tsv',index_col=None)
    with np.load(OUT/'Wormhole_sample_embeddings.npz') as b:te=b['training_embeddings'];qe=b['heldout_embeddings']
    assert (OUT/'Wormhole_parameters.msgpack').stat().st_size>0
    for stage in STAGES:assert (OUT/stage/'state_niche_embeddings.npz').exists()
    benchmark_distances(None,tm,qm,te,qe)
    dump('Wormhole_execution.json',{'executed':True,'package':'wassersteinwormhole0.3.7','JAX':jax.__version__,'devices':[str(x) for x in jax.devices()],
        'training_steps':CFG['wormhole_training_steps'],'training_clouds':len(tm),'heldout_clouds':len(qm),
        'embedding_dim':CFG['wormhole_embedding_dim'],'architecture':'2 transformer blocks,2 heads,hidden64; encoder/decoder Sinkhorn divergence epsilon0.1,100 iterations;16 decoder points',
        'paper_default_reproduction':False,'benchmark_resumed_from_saved_embeddings':True,
        'purpose':'CPU-constrained, region-held-out metric benchmark; adoption depends on reserved distance and prediction tests'})

def prediction_and_correspondence():
    from sklearn.metrics import adjusted_rand_score
    reports=[];profiles=[];archetypes=[];correspondence=[]
    bench=read(OUT/'Wormhole_reserved_distance_benchmark.tsv',index_col=None)
    metric=[]
    for s,z in bench.groupby('stage'):
        w=spearmanr(z.Sinkhorn_divergence,z.Wormhole_mean_squared_embedding_distance).statistic
        mean=spearmanr(z.Sinkhorn_divergence,z.mean_squared_feature_distance).statistic
        metric.append({'stage':s,'heldout_pairs':len(z),'Wormhole_vs_Sinkhorn_Spearman':w,'mean_vs_Sinkhorn_Spearman':mean,
            'Wormhole_vs_exact_W2_Spearman':spearmanr(z.exact_W2_squared,z.Wormhole_mean_squared_embedding_distance).statistic,
            'Wormhole_distance_MAE':np.mean(abs(z.Sinkhorn_divergence-z.Wormhole_mean_squared_embedding_distance)),
            'metric_gate':w>=.8 and w>=mean-.05})
    metric=pd.DataFrame(metric);save(metric,'Wormhole_metric_assessment.tsv')
    metric_ok=bool(metric.metric_gate.all())
    for stage in STAGES:
        with np.load(OUT/stage/'state_niche_embeddings.npz') as b:ids=b['cell_ids'].astype(str);emb=b['embedding'];means=b['mean_features'];test=b['heldout'];rows=b['focal_indices']
        c=read(UP/stage/'atlas_annotations.tsv.gz').loc[ids]
        old=read(UP/stage/'proximity50_CM_neighborhoods.tsv.gz').reindex(ids)
        valid=old.geometry_eligible.fillna(False).to_numpy()&c.CM_maturation_core.notna().to_numpy()
        y=c.CM_maturation_core.to_numpy();base=pd.DataFrame(index=c.index)
        for k in ['cell_area','CM_outcome_excluded_depth']:base['log_'+k]=np.log1p(c[k])
        base['nucleus_missing']=c.nucleus_area.isna().astype(float)
        for k in ['Endothelial_source_capture','Fibroblast_source_capture','Endothelial_source_coverage','Fibroblast_source_coverage','review_fraction']:
            base[k]=old[k]
        base['log_resolved_neighbors']=np.log1p(old.resolved_neighbors)
        for lineage in NC:base['fraction__'+lineage]=old['fraction__'+lineage]
        for k in ['arterial_distance_um','endocardial_proxy_distance_um']:base['log_'+k]=np.log1p(c[k])
        xy=c[['x_centroid','y_centroid']].to_numpy();from sklearn.preprocessing import SplineTransformer
        spatial=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(xy)
        for j in range(spatial.shape[1]):base['XY_'+str(j)]=spatial[:,j]
        arrays={'baseline':base.to_numpy(),'mean_state':np.column_stack([base,means]),'Wormhole':np.column_stack([base,emb])}
        train=(~test)&valid;assessment=test&valid
        for name,x in arrays.items():
            scale=StandardScaler().fit(x[train]);model=Ridge(alpha=10).fit(scale.transform(x[train]),y[train]);pred=model.predict(scale.transform(x[assessment]))
            reports.append({'stage':stage,'predictors':name,'training_CMs':int(train.sum()),'heldout_CMs':int(assessment.sum()),
                'heldout_R2':r2_score(y[assessment],pred),'heldout_MAE':np.mean(abs(y[assessment]-pred)),
                'scope':'One reserved spatial-block fold, fixed Ridge; not independent animals or causal evidence'})
        # Archetypes remain descriptive if learned metric validation fails.
        representation=emb if metric_ok else means
        scaling=StandardScaler().fit(representation[~test]);z=scaling.transform(representation)
        k=CFG['state_archetypes_per_stage'];km=KMeans(n_clusters=k,n_init=10,random_state=CFG['seed']).fit(z[~test]);labels=km.predict(z)
        aliases=[]
        q=read(OUT/stage/'CM_vascular_exposures.tsv.gz').reindex(ids)
        for group in range(k):
            take=labels==group;medoid=np.flatnonzero(take)[np.argmin(((z[take]-km.cluster_centers_[group])**2).sum(1))]
            archetypes.append({'stage':stage,'state':group,'niches':int(take.sum()),'medoid_CM':ids[medoid],
                'mean_maturation':c.loc[take,'CM_maturation_core'].mean(),
                'mean_growth_proxy':q.loc[take,'proximity50__Plexus_growth_proxy'].mean(),
                'mean_Notch_response':q.loc[take,'proximity50__Notch_response_proxy'].mean(),
                'mean_EC_Ptprm':q.loc[take,'proximity50__EC_sender_Ptprm'].mean(),
                'representation':'Wormhole' if metric_ok else 'Mean gene-state features; Wormhole metric gate failed',
                'status':'Exploratory state archetype; fixed K6, not a stable biological niche taxonomy'})
        profiles.append(pd.DataFrame({'stage':stage,'cell_id':ids,'state_archetype':labels,'heldout':test}))
        np.savez_compressed(OUT/stage/'state_archetype_representation.npz',cell_ids=ids,representation=representation,labels=labels,cluster_centers=km.cluster_centers_,scaled_representation=z)
    save(pd.DataFrame(reports),'Wormhole_maturation_prediction.tsv');save(pd.concat(profiles,ignore_index=True),'cross_sectional_state_assignments.tsv.gz')
    archetypes=pd.DataFrame(archetypes);save(archetypes,'cross_sectional_state_archetypes.tsv')
    # Compare actual medoid distributions across ages using exact empirical OT,
    # not an assumed cell lineage or a forced probability-conserving flow.
    for a,b in zip(STAGES[:-1],STAGES[1:]):
        sa=archetypes[archetypes.stage.eq(a)];sb=archetypes[archetypes.stage.eq(b)]
        stores={}
        for stage in [a,b]:
            with np.load(OUT/stage/'nonCM_state_points.npz') as d:points=d['points'];focal=d['focal_ids'].astype(str)
            member=sparse.load_npz(OUT/stage/'state_cloud_membership.npz').tocsr();stores[stage]=(points,focal,member)
        costs=np.zeros((len(sa),len(sb)))
        for i,ra in enumerate(sa.itertuples()):
            for j,rb in enumerate(sb.itertuples()):
                clouds=[]
                for stage,cid in [(a,ra.medoid_CM),(b,rb.medoid_CM)]:
                    points,focal,member=stores[stage];row=int(np.flatnonzero(focal==cid)[0]);pc,_=cloud_subset(points,member,row,CFG['seed']);clouds.append(pc)
                xa,xb=clouds;na,nb=len(xa),len(xb)
                constraints=sparse.vstack([sparse.kron(sparse.eye(na),np.ones((1,nb))),sparse.kron(np.ones((1,na)),sparse.eye(nb))]).tocsr()
                solved=linprog(cdist(xa,xb,'sqeuclidean').ravel(),A_eq=constraints,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs');assert solved.success
                costs[i,j]=solved.fun
        for i,ra in enumerate(sa.itertuples()):
            for j,rb in enumerate(sb.itertuples()):
                correspondence.append({'from_stage':a,'to_stage':b,'from_state':ra.state,'to_state':rb.state,'exact_medoid_W2_squared':costs[i,j],
                    'nearest_state_in_next_section':j==int(costs[i].argmin()),'nearest_state_in_previous_section':i==int(costs[:,j].argmin()),
                    'interpretation':'Phenotypic correspondence between separate sections, not tracked cells, lineage conversion, chronological rate or true dynamical reconstruction'})
    save(pd.DataFrame(correspondence),'cross_age_state_correspondence.tsv')
    dump('dynamic_interpretation_gate.json',{'Wormhole_metric_accepted_all_ages':metric_ok,'spatial_prediction_compared_to_mean_features':True,
        'archetypes_fixed_K_exploratory':True,'correspondences_use_exact_OT_on_observed_state_medoid_clouds':True,
        'tracked_cells':False,'endocardial_to_coronary_lineage_inference':False,'causal_influence':False,
        'limitations':'Missing fetal/P1-P6 intervals, one section per age, approximate state cost and source selection prevent reconstructing genuine cell trajectories or ancestry'})

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['prepare','fit','benchmark','assess','all'],default='all');args=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.step in ['prepare','all']:prepare()
    if args.step in ['fit','all']:fit_wormhole()
    if args.step=='benchmark':resume_benchmark()
    if args.step in ['assess','all']:prediction_and_correspondence()

if __name__=='__main__':main()
