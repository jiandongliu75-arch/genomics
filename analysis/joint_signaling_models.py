"""Receiver-controlled prediction with source-disjoint held-out regions."""
import os,time,json,argparse,gc
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from joint_signaling_common import *

SETS=['A','B','C','C_fine_context','D','E']
def feature_sets(columns,definitions):
    groups=definitions.set_index('name').feature_group.to_dict()
    included={'A':{'A'},'B':{'A','B'},'C':{'A','B','C'},'C_fine_context':{'A','B','C','C_fine_context'},
              'D':{'A','B','C','C_fine_context','D'},'E':{'A','B','C','C_fine_context','D','E'}}
    return {name:np.array([i for i,c in enumerate(columns) if ('A' if c.startswith('A__') else groups.get(c)) in allowed]) for name,allowed in included.items()}

def source_feature_groups(columns,definitions):
    meta=definitions.set_index('name');result={}
    for name in BROAD:
        keys=[]
        for i,col in enumerate(columns):
            if col not in meta.index or meta.loc[col,'feature_group'] not in ['C','D','E']:continue
            source=str(meta.loc[col,'source']);matches=(source=='Broad:'+name or broad(source)==name)
            if matches:keys.append(i)
        if keys:result['Broad:'+name]=np.array(keys)
    for source in ['EC:'+x for x in EC]+['FB:'+x for x in FB]:
        keys=[i for i,col in enumerate(columns) if col in meta.index and meta.loc[col,'feature_group'] in ['C','D','E'] and meta.loc[col,'source']==source]
        if keys:result[source]=np.array(keys)
    return result

def permutation_indices(labels,rng):
    perm=np.arange(len(labels))
    for label in np.unique(labels):
        ix=np.flatnonzero(labels==label)
        if len(ix)>1:perm[ix]=rng.permutation(ix)
    return perm

def fit_view(stage,view,refresh_ridge=False):
    from joint_signaling_capture import augment
    augment(stage,view)
    start=time.monotonic();cfg=config();folder=OUT/stage
    with np.load(folder/(view+'_predictor_features.npz')) as b:
        ids=b['cell_ids'];columns=b['columns'];X=b['features'].astype(float);Y=b['outcomes'].astype(float)
        eligible=b['eligible'];domain=b['domain'];states=b['CM_states'];xy=b['xy'].astype(float);outcomes=b['outcome_names']
    definitions=read(folder/'feature_definitions.tsv');sets=feature_sets(columns,definitions);sourcegroups=source_feature_groups(columns,definitions)
    assert all(len(sets[a])<=len(sets[b]) for a,b in zip(SETS,SETS[1:]))
    pred=np.full((2,len(SETS),len(ids),len(outcomes)),np.nan,np.float32)
    if refresh_ridge:
        with np.load(folder/(view+'_OOF_predictions.npz')) as z:pred[1]=z['predictions'][1]
    fingerprint=sha(folder/(view+'_predictor_features.npz'))+sha(__file__)+json.dumps(cfg,sort_keys=True)+str(refresh_ridge)
    fold_rows=[];metrics=[];importance=[];coefficients=[]
    counts=read(folder/'source_domains.tsv.gz',index_col=0)
    for fold in range(4):
        test=eligible&(domain==fold);train=eligible&(domain!=fold)
        if train.sum()<300 or test.sum()<100:raise RuntimeError(('Unsupported spatial fold',stage,view,fold,train.sum(),test.sum()))
        gap=float(cKDTree(xy[test]).query(xy[train])[0].min());assert gap>=400-1e-3,(stage,fold,gap)
        rng=np.random.default_rng(SEED+int(stage[1:])*100+fold)
        train_ids=np.flatnonzero(train)
        if len(train_ids)>cfg['maximum_training_CMs']:train_ids=np.sort(rng.choice(train_ids,cfg['maximum_training_CMs'],replace=False))
        test_ids=np.flatnonzero(test);yscale=StandardScaler().fit(Y[train_ids]);target=yscale.transform(Y[train_ids])
        checkpoint=folder/f'{view}_model_fold{fold}_checkpoint.npz'
        if checkpoint.exists():
            with np.load(checkpoint) as z:
                if str(z['fingerprint'])==fingerprint:
                    pred[:,:,test_ids]=z['predictions']
                    for dest,key in [(fold_rows,'fold_rows'),(metrics,'metrics'),(importance,'importance'),(coefficients,'coefficients')]:dest.extend(json.loads(str(z[key])))
                    print('JOINT_MODEL_CHECKPOINT_RESTORED',stage,view,fold,flush=True);continue
        starts=[len(v) for v in [fold_rows,metrics,importance,coefficients]]
        fold_rows.append({'stage':stage,'view':view,'fold':fold,'available_training_CMs':int(train.sum()),'fitted_training_CMs':len(train_ids),
                          'test_CMs':len(test_ids),'minimum_receiver_gap_um':gap,'train_test_source_ID_overlap':0,
                          'interpretation':'Each source is assigned to one domain; receivers have 200-um interior margin; regions are not independent hearts'})
        for learner_index,learner in enumerate(['Ridge','RandomForest']):
            if refresh_ridge and learner=='RandomForest':continue
            for set_index,name in enumerate(SETS):
                use=sets[name];a=X[train_ids][:,use];b=X[test_ids][:,use]
                scale=StandardScaler().fit(a)
                if learner=='Ridge':
                    molecular=np.array([not columns[j].startswith('A__') for j in use])
                    scale.scale_[molecular]=np.maximum(scale.scale_[molecular],cfg['ridge_minimum_molecular_scale'])
                a=scale.transform(a);b=scale.transform(b)
                model=(Ridge(alpha=cfg['ridge_alpha']) if learner=='Ridge' else RandomForestRegressor(
                    n_estimators=cfg['forest_trees'],max_depth=cfg['forest_max_depth'],min_samples_leaf=cfg['forest_min_leaf'],
                    max_features=cfg['forest_max_features'],n_jobs=4,random_state=SEED+fold))
                model.fit(a,target);prediction=yscale.inverse_transform(model.predict(b));pred[learner_index,set_index,test_ids]=prediction
                for j,outcome in enumerate(outcomes):
                    metrics.append({'stage':stage,'view':view,'fold':fold,'learner':learner,'feature_set':name,'outcome':outcome,'test_CMs':len(test_ids),
                                    'R2':r2_score(Y[test_ids,j],prediction[:,j]),'mean_squared_error':np.mean((Y[test_ids,j]-prediction[:,j])**2),
                                    'scope':'Unseen source-disjoint region in one observed section; fixed model complexity'})
                if learner=='Ridge' and name=='E':
                    for k,col in enumerate(columns[use]):
                        if col.startswith('A__'):continue
                        for j,outcome in enumerate(outcomes):coefficients.append({'stage':stage,'view':view,'fold':fold,'feature':col,'outcome':outcome,'training_scaled_ridge_coefficient':model.coef_[j,k]})
                if learner=='RandomForest' and name=='E':
                    baseline_error=np.mean((Y[test_ids]-prediction)**2,axis=0)
                    local_xy=(xy[test_ids]//200).astype(int)
                    for nullname,labels in [('within_CM_state',states[test_ids]),('within_CM_state_region200',np.array([s+'::'+str(x)+'_'+str(y) for s,(x,y) in zip(states[test_ids],local_xy)]))]:
                        perm=permutation_indices(labels,rng);moved=float(np.mean(perm!=np.arange(len(perm))))
                        for source,cols in sourcegroups.items():
                            local=np.flatnonzero(np.isin(use,cols))
                            if not len(local):continue
                            changed=b.copy();changed[:,local]=b[perm][:,local]
                            pp=yscale.inverse_transform(model.predict(changed));err=np.mean((Y[test_ids]-pp)**2,axis=0)
                            for j,outcome in enumerate(outcomes):importance.append({'stage':stage,'view':view,'fold':fold,'source_group':source,'perturbation':nullname,
                                'outcome':outcome,'test_CMs':len(test_ids),'perturbed_features':len(local),'fraction_rows_moved':moved,
                                'increase_MSE':err[j]-baseline_error[j],'increase_MSE_fraction_of_test_variance':(err[j]-baseline_error[j])/max(np.var(Y[test_ids,j]),1e-12),
                                'scope':'Predictive source-feature perturbation; correlated pooled/fine inputs remain; not a causal source effect'})
                print('JOINT_MODEL_FIT',stage,view,fold,learner,name,len(use),flush=True)
                del model,a,b;gc.collect()
        payload={key:np.array(json.dumps(v[lo:])) for key,v,lo in zip(['fold_rows','metrics','importance','coefficients'],[fold_rows,metrics,importance,coefficients],starts)}
        np.savez_compressed(checkpoint,fingerprint=np.array(fingerprint),predictions=pred[:,:,test_ids],**payload)
    np.savez_compressed(folder/(view+'_OOF_predictions.npz'),cell_ids=ids,domains=domain,CM_states=states,xy=xy.astype(np.float32),eligible=eligible,
        outcomes=Y.astype(np.float32),outcome_names=outcomes,predictions=pred,learners=np.array(['Ridge','RandomForest']),feature_sets=np.array(SETS),
        capture_controls_version=np.array(1),predictor_columns=columns,ridge_minimum_molecular_scale=np.array(cfg['ridge_minimum_molecular_scale']))
    if refresh_ridge:
        existing=read(folder/(view+'_fold_performance.tsv'))
        metrics=pd.concat([existing[existing.learner.eq('RandomForest')],pd.DataFrame(metrics)],ignore_index=True)
    else:metrics=pd.DataFrame(metrics)
    save(pd.DataFrame(fold_rows),stage+'/'+view+'_prediction_fold_checks.tsv');save(metrics,stage+'/'+view+'_fold_performance.tsv')
    if not refresh_ridge:save(pd.DataFrame(importance),stage+'/'+view+'_source_feature_perturbations.tsv')
    save(pd.DataFrame(coefficients),stage+'/'+view+'_ridge_feature_coefficients.tsv.gz')
    summary=[]
    for li,learner in enumerate(['Ridge','RandomForest']):
        for si,name in enumerate(SETS):
            for state in ['All',*sorted(np.unique(states[eligible]))]:
                take=eligible if state=='All' else eligible&(states==state)
                if take.sum()<100:continue
                for j,outcome in enumerate(outcomes):summary.append({'stage':stage,'view':view,'learner':learner,'feature_set':name,'CM_state':state,
                    'outcome':outcome,'OOF_CMs':int(take.sum()),'OOF_R2':r2_score(Y[take,j],pred[li,si,take,j]),
                    'OOF_MSE':np.mean((Y[take,j]-pred[li,si,take,j])**2),'scope':'Held-out contiguous domains within one section'})
    save(pd.DataFrame(summary),stage+'/'+view+'_OOF_performance.tsv')
    for checkpoint in folder.glob(view+'_model_fold*_checkpoint.npz'):checkpoint.unlink()
    print('JOINT_MODEL_VIEW_COMPLETE',stage,view,round(time.monotonic()-start,1),flush=True)

def summarize():
    performance=[];folds=[];importance=[];gain_rows=[]
    rng=np.random.default_rng(SEED);draws=512
    comparisons=[('B','A'),('C','B'),('C_fine_context','C'),('D','C_fine_context'),('E','C_fine_context'),('E','D')]
    for stage in STAGES:
        for view in ['core','native']:
            folder=OUT/stage
            performance.append(read(folder/(view+'_OOF_performance.tsv')));folds.append(read(folder/(view+'_prediction_fold_checks.tsv')))
            importance.append(read(folder/(view+'_source_feature_perturbations.tsv')))
            with np.load(folder/(view+'_OOF_predictions.npz')) as z:
                take=z['eligible'];xy=z['xy'][take];Y=z['outcomes'][take];P=z['predictions'][:,:,take];outcomes=z['outcome_names'];states=z['CM_states'][take]
                domain=z['domains'][take];sets=list(z['feature_sets']);learners=list(z['learners'])
            # Stratified spatial-block uncertainty conditional on already fitted OOF models.
            blocks=np.floor(xy/800).astype(int);key=np.array([str(d)+':'+str(x)+'_'+str(y) for d,(x,y) in zip(domain,blocks)])
            unique,code=np.unique(key,return_inverse=True);nb=len(unique);w=np.zeros((draws,nb),int)
            for d in np.unique(domain):
                indices=np.unique(code[domain==d])
                for j in range(draws):w[j]+=np.bincount(rng.choice(indices,len(indices),replace=True),minlength=nb)
            for li,learner in enumerate(learners):
                for larger,smaller in comparisons:
                    hi=sets.index(larger);lo=sets.index(smaller)
                    improvement=(Y-P[li,lo])**2-(Y-P[li,hi])**2
                    for state in ['All',*sorted(np.unique(states))]:
                        mask=np.ones(len(Y),bool) if state=='All' else states==state
                        if mask.sum()<200:continue
                        for j,outcome in enumerate(outcomes):
                            v=improvement[mask,j];cc=code[mask];var=np.var(Y[mask,j]);counts=np.bincount(cc,minlength=nb)
                            sums=np.bincount(cc,weights=v,minlength=nb);boot=(w@sums)/np.maximum(w@counts,1)/max(var,1e-12)
                            point=v.mean()/max(var,1e-12);a,b=np.quantile(boot,[.025,.975]);se=boot.std(ddof=1)
                            gain_rows.append({'stage':stage,'view':view,'learner':learner,'comparison':larger+' minus '+smaller,'CM_state':state,
                                'outcome':outcome,'CM_profiles':int(mask.sum()),'blocks800':len(np.unique(cc)),'delta_R2':point,'CI_low':a,'CI_high':b,
                                'within_section_p':2*norm.sf(abs(point)/max(se,1e-12)),
                                'primary_analysis':view=='core' and learner=='RandomForest' and state=='All' and larger in ['D','E'] and smaller=='C_fine_context',
                                'scope':'Fixed-model paired OOF error contrast; spatial uncertainty within one section; no training/annotation uncertainty propagated'})
    gain=pd.DataFrame(gain_rows);gain['primary_grid_BH_q']=np.nan
    ix=gain.primary_analysis&gain.blocks800.ge(8)
    p=gain.loc[ix,'within_section_p'].to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    gain.loc[ix,'primary_grid_BH_q']=q
    save(pd.concat(performance),'OOF_performance.tsv');save(pd.concat(folds),'prediction_fold_checks.tsv')
    save(pd.concat(importance),'source_feature_perturbations.tsv');save(gain,'incremental_prediction_gains.tsv')
    print('JOINT_MODEL_SUMMARY_COMPLETE',len(gain),'primary tests',int(ix.sum()),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=STAGES);p.add_argument('--view',choices=['core','native']);p.add_argument('--summarize',action='store_true');p.add_argument('--resume',action='store_true');p.add_argument('--wait',action='store_true');p.add_argument('--refresh-ridge',action='store_true');a=p.parse_args()
    if a.summarize:summarize();return
    for stage in ([a.stage] if a.stage else STAGES):
        for view in ([a.view] if a.view else ['core','native']):
            refresh=a.refresh_ridge
            if a.resume and (OUT/stage/(view+'_OOF_performance.tsv')).exists():
                with np.load(OUT/stage/(view+'_OOF_predictions.npz')) as z:
                    if 'capture_controls_version' in z and int(z['capture_controls_version'])==1:
                        if 'ridge_minimum_molecular_scale' in z and float(z['ridge_minimum_molecular_scale'])==config()['ridge_minimum_molecular_scale'] and not refresh:continue
                        refresh=True
            while a.wait and not all((OUT/stage/(view+'_'+suffix)).exists() for suffix in ['predictor_features.npz','transport_diagnostics.tsv','source_pair_communication.tsv']):time.sleep(5)
            fit_view(stage,view,refresh_ridge=refresh)

if __name__=='__main__':main()
