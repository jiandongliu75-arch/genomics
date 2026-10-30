"""Outcome-blinded graph niches, local-null enrichment and optional OT pilot."""
from cardiac_spatial_pipeline import *
from itertools import combinations
from scipy.optimize import linprog
from scipy.spatial.distance import cdist

def features(stage):
    f=frame(stage);g=graph(stage,'contact');p,n=normalize_rows(g)
    codes=pd.Categorical(f.atlas_identity,categories=TYPES).codes;assert (codes>=0).all()
    x=np.eye(len(TYPES),dtype=np.float32)[codes]
    first=np.asarray(p@x);second=np.asarray(p@first)
    diagonal=np.asarray(p.multiply(p.T).sum(1)).ravel()
    second=(second-diagonal[:,None]*x)/np.maximum(1-diagonal,1e-8)[:,None]
    second=np.clip(second,0,1)
    assert np.allclose(first[n>0].sum(1),1)
    first_mass=first[:,:len(LINEAGES)].sum(1);second_mass=second[:,:len(LINEAGES)].sum(1)
    # Review burden is diagnostic, never a biological clustering feature.
    first_res=first[:,:len(LINEAGES)]/np.maximum(first_mass,1e-8)[:,None]
    second_res=second[:,:len(LINEAGES)]/np.maximum(second_mass,1e-8)[:,None]
    values=np.column_stack([first_res,second_res]).astype(np.float32)
    names=[str(h)+'hop__'+t for h in [1,2] for t in LINEAGES]
    eligible=(n*first_mass>=3-1e-6)&(first_mass>=CFG['niche_minimum_resolved_neighbor_fraction'])&(second_mass>=CFG['niche_minimum_resolved_neighbor_fraction'])
    return f,values,eligible,names

def niches(selection_only=False):
    samples=[];tests=[];sizes=[];counts=[];columns=None
    for stage in STAGES:
        print('NICHE_FEATURES',stage,flush=True)
        f,x,accepted,names=features(stage);columns=names
        rng=np.random.default_rng(CFG['seed']);eligible=np.flatnonzero(accepted)
        xy=f[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/CFG['spatial_block_um']).astype(int)
        fold=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
        train=eligible[fold[eligible]!=0];test=eligible[fold[eligible]==0]
        take=rng.choice(train,min(len(train),CFG['niche_training_profiles_per_stage']),replace=False)
        test=rng.choice(test,min(len(test),3000),replace=False)
        samples.append(x[take]);tests.append(x[test]);sizes.append((stage,len(f)))
        np.savez_compressed(OUT/stage/'outcome_blind_niche_features.npz',features=x,columns=np.array(names),eligible=accepted,
            training_profile_indices=take,assessment_profile_indices=test)
    train=np.vstack(samples);assessment=np.vstack(tests)
    scaler=StandardScaler().fit(np.sqrt(train))
    scaled=scaler.transform(np.sqrt(train));qt=scaler.transform(np.sqrt(assessment))
    # Composition-only features use no focal RNA, maturity values or global XY.
    rows=[];models={};alllabels={}
    wanted=CFG['niche_k']
    if selection_only:
        previous=pd.read_csv(OUT/'niche_model_selection.tsv',sep='\t')
        accepted=previous.loc[previous.minimum_seed_ARI.ge(.8)&previous.all_fits_converged]
        selected=int(accepted.sort_values('train_BIC').iloc[0].k) if len(accepted) else int(previous.sort_values('median_seed_ARI',ascending=False).iloc[0].k)
        wanted=[selected]
    for k in wanted:
        labels=[];fits=[]
        for seed in CFG['niche_seeds']:
            print('NICHE_FIT',k,seed,flush=True)
            model=GaussianMixture(n_components=k,covariance_type='full',reg_covar=.01,max_iter=300,n_init=CFG['niche_initializations'],random_state=seed).fit(scaled)
            pred=model.predict(qt);labels.append(pred);fits.append(model)
        ari=[adjusted_rand_score(labels[i],labels[j]) for i,j in combinations(range(len(labels)),2)]
        best=int(np.argmin([m.bic(scaled) for m in fits]));models[k]=fits[best];alllabels[k]=labels[best]
        rows.append({'k':k,'median_seed_ARI':np.median(ari),'minimum_seed_ARI':min(ari),'train_BIC':fits[best].bic(scaled),
            'heldout_mean_log_likelihood':fits[best].score(qt),'all_fits_converged':all(m.converged_ for m in fits),
            'training_profiles':len(train),'assessment_profiles':len(assessment),'feature_set':'two-step contact-neighbor composition; focal return removed; no CM RNA'})
    comparison=previous if selection_only else pd.DataFrame(rows)
    stable=comparison.loc[comparison.minimum_seed_ARI.ge(.8)&comparison.all_fits_converged]
    if stable.empty:
        chosen=int(comparison.sort_values(['median_seed_ARI','train_BIC'],ascending=[False,True]).iloc[0].k)
        status='Exploratory: no tested K met seed-stability gate'
    else:
        chosen=int(stable.sort_values('train_BIC').iloc[0].k);status='Seed-stable computational niches; not independently validated anatomy'
    model=models[chosen];save(comparison,'niche_model_selection.tsv')
    # Name components by their one-hop resolved composition, without invoking anatomy.
    centers=scaler.inverse_transform(model.means_)**2
    aliases={}
    reference=train[:,:len(LINEAGES)].mean(0)
    for k in range(chosen):
        component=centers[k,:len(LINEAGES)];component=component/max(component.sum(),1e-8)
        enrichment=(component+.01)/(reference+.01)
        order=np.argsort(-enrichment)[:2]
        aliases[k]=f'N{k+1}: '+LINEAGES[order[0]]+' / '+LINEAGES[order[1]]+' enriched context'
    model_arrays={'means':model.means_,'covariances':model.covariances_,'weights':model.weights_,
                  'scaler_mean':scaler.mean_,'scaler_scale':scaler.scale_,'feature_columns':np.array(columns)}
    np.savez_compressed(OUT/'niche_model_parameters.npz',**model_arrays)
    for stage,_ in sizes:
        f=frame(stage)
        with np.load(OUT/stage/'outcome_blind_niche_features.npz') as b:x=b['features'];eligible=b['eligible']
        z=scaler.transform(np.sqrt(x));posterior=model.predict_proba(z);label=posterior.argmax(1)
        good=eligible&(posterior.max(1)>=CFG['niche_minimum_posterior'])
        f['niche_candidate']=[aliases[k] for k in label]
        f['niche']=np.where(good,f.niche_candidate,'Niche unresolved')
        f['niche_maximum_posterior']=posterior.max(1)
        f['niche_entropy']=-np.sum(posterior*np.log(np.maximum(posterior,1e-12)),axis=1)/np.log(chosen)
        f['niche_assignment_status']=status
        save(f[['stage_cell_id','x_centroid','y_centroid','atlas_identity','fine_identity','niche','niche_candidate',
                'niche_maximum_posterior','niche_entropy']],f'{stage}/niche_annotations.tsv.gz',index=True)
        exported=ad.AnnData(obs=f[['stage_cell_id','atlas_identity','fine_identity','niche']].copy())
        exported.obsm['spatial']=f[['x_centroid','y_centroid']].to_numpy()
        exported.obsm['X_neighborhood']=x
        exported.obsp['spatial_connectivities']=graph(stage,'contact')
        exported.uns['feature_names']=columns;exported.uns['measured_RNA_matrix']='Not duplicated; see original nuclear_anchored_counts.h5ad and source manifest'
        exported.uns['CellCharter_run']=False;exported.uns['Wasserstein_Wormhole_run']=False
        exported.write_h5ad(OUT/stage/'niche_analysis_ready.h5ad',compression='gzip')
        for niche,zg in f.groupby('niche',observed=True):
            counts.append({'stage':stage,'niche':niche,'profiles':len(zg),**{l:int(zg.atlas_identity.eq(l).sum()) for l in LINEAGES}})
    save(pd.DataFrame(counts),'niche_stage_counts.tsv')
    dump('niche_methods.json',{'chosen_k':chosen,'status':status,'aliases':aliases,'backend':'sklearn GaussianMixture, not CellCharter',
        'feature_names':columns,'outcome_blind':True,'focal_two_step_return_removed':True,
        'selection':'Stage-balanced block-split training; minimum pairwise seed ARI>=0.8 on reserved blocks; eligible K by BIC',
        'review_burden_policy':'Review categories excluded from features; require >=50% resolved neighbor mass in both steps and >=3 resolved direct neighbors',
        'limits':'Seed stability is not replicate/anatomical validation; no age or outcome used as clustering predictors; local compositions can still vary by age/batch'})

def enrichments():
    rows=[]
    for stage in STAGES:
        print('ENRICHMENT',stage,flush=True)
        f=frame(stage);g=sparse.triu(graph(stage,'contact'),1).tocoo()
        resolved=f.atlas_resolved.to_numpy(bool);keep=resolved[g.row]&resolved[g.col]
        ii,jj=g.row[keep],g.col[keep]
        types=f.atlas_identity.map({l:i for i,l in enumerate(LINEAGES)}).fillna(-1).to_numpy(int)
        def count(t):
            a=np.minimum(t[ii],t[jj]);b=np.maximum(t[ii],t[jj])
            return np.bincount(a*len(LINEAGES)+b,minlength=len(LINEAGES)**2)
        observed=count(types)
        xy=f[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/CFG['enrichment_shuffle_block_um']).astype(int)
        _,bins=np.unique(blocks,axis=0,return_inverse=True)
        groups=[np.flatnonzero((bins==b)&resolved) for b in np.unique(bins)]
        rng=np.random.default_rng(CFG['seed']);null=[]
        for _ in range(CFG['enrichment_permutations']):
            shuffled=types.copy()
            for take in groups:shuffled[take]=rng.permutation(shuffled[take])
            null.append(count(shuffled))
        null=np.asarray(null)
        for a in range(len(LINEAGES)):
            for b in range(a,len(LINEAGES)):
                k=a*len(LINEAGES)+b;expected=null[:,k].mean()
                p=(1+np.sum(np.abs(null[:,k]-expected)>=abs(observed[k]-expected)))/(len(null)+1)
                rows.append({'stage':stage,'cell_type_1':LINEAGES[a],'cell_type_2':LINEAGES[b],'observed_contacts':observed[k],
                    'local_shuffle_mean':expected,'log2_observed_over_local_null':np.log2((observed[k]+.5)/(expected+.5)),
                    'within_section_permutation_p':p,'permutations':len(null),'null':'labels permuted within200um blocks; review positions fixed',
                    'scope':'conditional spatial co-occurrence, not animal-level significance or communication'})
    table=pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    table['within_section_BH_q']=np.nan
    for s,z in table.groupby('stage'):table.loc[z.index,'within_section_BH_q']=multipletests(z.within_section_permutation_p,method='fdr_bh')[1]
    save(table,'contact_enrichment_local_null.tsv')

def transport_pilot():
    """Small exact OT/SW comparison before considering a learned accelerator.

    These are distributions of non-CM observed programs/types, not RNA-to-space
    alignments. No CM maturation or focal CM expression enters their cost.
    """
    rows=[];definitions=[]
    for stage in STAGES:
        print('TRANSPORT_PILOT',stage,flush=True)
        f=frame(stage);source=f.atlas_identity.isin(NC).to_numpy();ids=np.flatnonzero(source)
        points=np.column_stack([f.atlas_identity.eq(l).to_numpy(float) for l in NC]+[f[p].fillna(0).to_numpy() for p in PROGS]+
                               [f[p].notna().to_numpy(float) for p in PROGS])[source].astype(np.float32)
        names=['lineage__'+l for l in NC]+PROGS+['observed__'+p for p in PROGS]
        # Program scaling uses source state only; other-lineage zeros stay zero.
        for j,p in enumerate(PROGS,start=len(NC)):
            valid=f[p].notna().to_numpy()[source]
            if valid.sum():points[valid,j]/=max(points[valid,j].std(),.25)
        points=np.clip(points,0,5)/np.sqrt(points.shape[1])
        membership=graph(stage,'proximity50')[:,ids].tocsr()
        # All cells are available as focal neighborhoods; the pilot evaluates CMs.
        sparse.save_npz(OUT/stage/'nonCM_pointcloud_membership.npz',membership)
        np.savez_compressed(OUT/stage/'nonCM_pointcloud_features.npz',features=points,feature_names=np.array(names),
                           source_profile_ids=f.index[ids].to_numpy(dtype=str),focal_profile_ids=f.index.to_numpy(dtype=str))
        eligible=np.flatnonzero(f.ventricular_CM_strict.to_numpy()&(np.diff(membership.indptr)>=3))
        rng=np.random.default_rng(CFG['seed']);projections=rng.normal(size=(points.shape[1],CFG['sliced_wasserstein_projections']))
        projections/=np.linalg.norm(projections,axis=0)
        xy=f[['x_centroid','y_centroid']].to_numpy();quant=np.linspace(0,1,32)
        audit_a=np.full((CFG['sliced_wasserstein_pairs'],32),-1,np.int32);audit_b=audit_a.copy()
        for repeat in range(CFG['sliced_wasserstein_pairs']):
            a,b=rng.choice(eligible,2,replace=False)
            for _ in range(20):
                if np.linalg.norm(xy[a]-xy[b])>=400:break
                b=int(rng.choice(eligible))
            ia=membership.indices[membership.indptr[a]:membership.indptr[a+1]]
            ib=membership.indices[membership.indptr[b]:membership.indptr[b+1]]
            if len(ia)>32:ia=rng.choice(ia,32,replace=False)
            if len(ib)>32:ib=rng.choice(ib,32,replace=False)
            audit_a[repeat,:len(ia)]=ia;audit_b[repeat,:len(ib)]=ib
            xa,xb=points[ia],points[ib]
            qa=np.quantile(xa@projections,quant,axis=0);qb=np.quantile(xb@projections,quant,axis=0)
            sw=np.sqrt(np.mean((qa-qb)**2));mean=np.linalg.norm(xa.mean(0)-xb.mean(0))
            exact=np.nan
            if repeat<64:
                na,nb=len(xa),len(xb)
                constraints=sparse.vstack([sparse.kron(sparse.eye(na),np.ones((1,nb))),sparse.kron(np.ones((1,na)),sparse.eye(nb))]).tocsr()
                result=linprog(cdist(xa,xb,'sqeuclidean').ravel(),A_eq=constraints,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs')
                assert result.success,result.message
                exact=np.sqrt(max(result.fun,0))
            rows.append({'stage':stage,'CM_id_a':f.index[a],'CM_id_b':f.index[b],'points_a':len(xa),'points_b':len(xb),
                'centroid_separation_um':np.linalg.norm(xy[a]-xy[b]),'sliced_W2':sw,'mean_feature_distance':mean,'exact_empirical_W2':exact})
        definitions.append({'stage':stage,'point_features':len(names),'source_nonCM_profiles':len(ids),'pilot_pairs':CFG['sliced_wasserstein_pairs'],
                            'exact_pairs':64,'Wormhole_trained':False,'meaning':'Within-stage observed non-CM type/program distributions, uniform mass; no CM outcomes or inferred sprouting'})
        np.savez_compressed(OUT/stage/'transport_pilot_point_indices.npz',point_indices_a=audit_a,point_indices_b=audit_b,projections=projections)
    save(pd.DataFrame(rows),'transport_distance_pilot.tsv.gz');save(pd.DataFrame(definitions),'transport_pilot_scope.tsv')
    summaries=[];d=pd.DataFrame(rows)
    for s,z in d.groupby('stage'):
        exact=z.loc[z.exact_empirical_W2.notna()]
        summaries.append({'stage':s,'pairs':len(z),'exact_pairs':len(exact),
            'SW_vs_exact_Spearman':spearmanr(exact.sliced_W2,exact.exact_empirical_W2).statistic,
            'mean_vs_exact_Spearman':spearmanr(exact.mean_feature_distance,exact.exact_empirical_W2).statistic,
            'SW_vs_mean_Spearman':spearmanr(z.sliced_W2,z.mean_feature_distance).statistic,
            'scope':'Distance pilot only; no evidence of improved CM prediction or causal mechanism'})
    save(pd.DataFrame(summaries),'transport_pilot_summary.tsv')

def optional_cellcharter():
    """Explicit opt-in backend on the exported, outcome-blinded niche features."""
    try:from cellcharter.tl import Cluster
    except ImportError as e:raise RuntimeError('CellCharter is optional and is not installed. Primary results use sklearn GMM; install CellCharter in a separate environment before this benchmark.') from e
    with np.load(OUT/'niche_model_parameters.npz') as b:mu=b['scaler_mean'];sd=b['scaler_scale']
    training=[];arrays=[]
    for stage in STAGES:
        with np.load(OUT/stage/'outcome_blind_niche_features.npz') as b:x=(np.sqrt(b['features'])-mu)/sd;take=b['training_profile_indices']
        training.append(x[take]);arrays.append((stage,x))
    k=json.loads((OUT/'niche_methods.json').read_text())['chosen_k']
    data=ad.AnnData(X=np.vstack(training))
    model=Cluster(n_clusters=k,random_state=CFG['seed']);model.fit(data,use_rep=None)
    for stage,x in arrays:
        q=ad.AnnData(X=x);labels=np.asarray(model.predict(q,use_rep=None)).astype(str)
        f=frame(stage);save(pd.DataFrame({'cell_id':f.index,'CellCharter_benchmark_cluster':labels}),f'{stage}/optional_CellCharter_benchmark.tsv.gz')
    dump('optional_CellCharter_status.json',{'executed':True,'purpose':'Backend sensitivity on fixed, outcome-blinded neighbor features; not the default published taxonomy'})

def optional_wormhole():
    """Verified0.3.7 API, explicit opt-in training in a separate JAX environment."""
    try:
        from wassersteinwormhole import Wormhole
        import jax
    except ImportError as e:raise RuntimeError('Wasserstein Wormhole is optional and is not installed. Observed non-CM point clouds and the exact/sliced OT pilot are already exported; install wassersteinwormhole==0.3.7 in a separate environment before training.') from e
    trainclouds=[];testclouds=[];trainids=[];testids=[];rng=np.random.default_rng(CFG['seed'])
    for stage in STAGES:
        f=frame(stage)
        with np.load(OUT/stage/'nonCM_pointcloud_features.npz') as b:points=b['features']
        membership=sparse.load_npz(OUT/stage/'nonCM_pointcloud_membership.npz').tocsr()
        eligible=np.flatnonzero(f.ventricular_CM_strict.to_numpy()&(np.diff(membership.indptr)>=3))
        xy=f[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/400).astype(int)
        fold=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
        test=eligible[fold[eligible]==0];train=eligible[fold[eligible]!=0]
        train=train[cKDTree(xy[test]).query(xy[train])[0]>=150]
        for indices,target,names,limit in [(train,trainclouds,trainids,2000),(test,testclouds,testids,500)]:
            for i in rng.choice(indices,min(limit,len(indices)),replace=False):
                ix=membership.indices[membership.indptr[i]:membership.indptr[i+1]]
                if len(ix)>32:ix=rng.choice(ix,32,replace=False)
                target.append(points[ix]);names.append(stage+'::'+str(f.index[i]))
    model=Wormhole(point_clouds=trainclouds,point_clouds_test=testclouds)
    model.train(training_steps=10000,key=jax.random.key(CFG['seed']))
    train_embeddings=np.asarray(model.encode(model.point_clouds,model.weights))
    test_embeddings=np.asarray(model.encode(model.point_clouds_test,model.weights_test))
    np.savez_compressed(OUT/'optional_Wormhole_embeddings.npz',training_ids=np.array(trainids),assessment_ids=np.array(testids),
        training_embeddings=train_embeddings,assessment_embeddings=test_embeddings)
    dump('optional_Wormhole_status.json',{'executed':True,'training_steps':10000,'source':'Official wassersteinwormhole0.3.7 API',
        'scope':'Distribution embedding benchmark; requires held-out OT and CM-prediction assessment before biological interpretation',
        'no_CM_maturation_or_focal_CM_RNA_in_cost':True})

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['niches','niches_select','enrichment','transport','cellcharter','wormhole','all'],default='all');args=p.parse_args()
    for name,fn in [('niches',niches),('enrichment',enrichments),('transport',transport_pilot)]:
        if args.step in [name,'all']:fn()
    if args.step=='cellcharter':optional_cellcharter()
    if args.step=='wormhole':optional_wormhole()
    if args.step=='niches_select':niches(selection_only=True)

if __name__=='__main__':main()
