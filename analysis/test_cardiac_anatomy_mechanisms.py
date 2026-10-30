"""Descriptive reconstruction, shared-anatomy tests and mechanism comparison.

Within-section block inference does not add independent animals. Components,
mediator conditioning and source/geometry sensitivities remain explicit.
"""
from reconstruct_cardiac_spatial import ROOT,OUT,UP,CFG,STAGES,read,save,dump
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr,t
from scipy.spatial import cKDTree
from sklearn.preprocessing import StandardScaler,SplineTransformer
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import multipletests
import argparse,gc

COMPONENTS=['maturation21','metabolic','sarcomere','excitation_contraction','withdrawal']

def exposures(c,e,source='author_capillary'):
    x=pd.DataFrame(index=c.index)
    x['capillary_area']=np.log1p(e[source+'__nearby_cell_area_fraction_50um'])
    x['capillary_apposition']=e[source+'__apposition_fraction_2um']
    x['pericyte_context']=np.log1p(e.author_pericyte__nearby_profiles_50um)
    x['FB_basement_program']=e.FB__core__FB_basement_matrix
    x['FB_fibrillar_program']=e.FB__core__FB_fibrillar_matrix
    x['EC_cycling']=e.EC__core__EC_cycling
    x['EC_Ptprm']=e.EC__nuclear__Ptprm
    return x

def covariates(c,e,technical=False,state=False):
    b=pd.DataFrame(index=c.index)
    for key in ['distance_to_large_void_um','distance_to_outer_um']:b['log_'+key]=np.log1p(c[key])
    b['local_void_fraction']=c.local_void_fraction;b['relative_wall_geometry']=c.relative_wall_geometry
    if technical:
        for key in ['native_panel_depth','cell_area']:b['log_'+key]=np.log1p(c[key])
        b['nucleus_area_missing']=c.nucleus_area_missing.astype(float)
        b['EC_source_capture']=e.EC__log_core_capture;b['FB_source_capture']=e.FB__log_core_capture
        b['log_neighbor_count']=np.log1p(e.resolved_neighbors_50um)
    if state:b=pd.concat([b,pd.get_dummies(c.author_atlas_cell_state,prefix='CMstate',dtype=float)],axis=1)
    return b

def demean(v,labels):
    _,inv=np.unique(labels,return_inverse=True);counts=np.bincount(inv)
    sums=np.zeros((len(counts),v.shape[1]));np.add.at(sums,inv,v)
    return v-sums[inv]/counts[inv,None]

def fit(y,b,e,c,within=True):
    values=np.column_stack([b.to_numpy(float),e.to_numpy(float)]);names=list(b)+list(e)
    sd=values.std(0);use=sd>1e-9;xx=(values[:,use]-values[:,use].mean(0))/sd[use]
    yy=y.to_numpy(float);ys=np.maximum(yy.std(0),1e-9);yy=(yy-yy.mean(0))/ys
    if within:
        xx=demean(xx,c.anatomy_block_200um.to_numpy());yy=demean(yy,c.anatomy_block_200um.to_numpy())
    xx=np.column_stack([np.ones(len(xx)),xx]);beta=np.linalg.lstsq(xx,yy,rcond=1e-9)[0]
    xy=c[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/CFG['bootstrap_block_um']).astype(int);_,inv=np.unique(blocks,axis=0,return_inverse=True);nb=inv.max()+1
    H=np.zeros((nb,xx.shape[1],xx.shape[1]));R=np.zeros((nb,xx.shape[1],yy.shape[1]))
    for k in range(nb):
        z=inv==k;H[k]=xx[z].T@xx[z];R[k]=xx[z].T@yy[z]
    draws=[];rng=np.random.default_rng(CFG['seed'])
    for _ in range(CFG['bootstrap_draws']):
        w=np.bincount(rng.integers(nb,size=nb),minlength=nb);draws.append(np.linalg.pinv(np.einsum('b,bij->ij',w,H),rcond=1e-9)@np.einsum('b,bik->ik',w,R))
    draws=np.asarray(draws);lo,hi=np.quantile(draws,[.025,.975],axis=0);se=draws.std(0);selected=np.flatnonzero(use);rows=[]
    for j,name in enumerate(names):
        if j<len(b.columns) or j not in selected:continue
        k=int(np.flatnonzero(selected==j)[0])+1;keep=np.arange(xx.shape[1])!=k
        er=xx[:,k]-xx[:,keep]@np.linalg.lstsq(xx[:,keep],xx[:,k],rcond=1e-9)[0]
        for out,outcome in enumerate(y):rows.append({'exposure':name,'outcome':outcome,'standardized_coefficient':beta[k,out],
            'block_CI_low':lo[k,out],'block_CI_high':hi[k,out],'CMs':len(c),'bootstrap_blocks':nb,
            'approximate_block_p':2*t.sf(abs(beta[k,out])/max(se[k,out],1e-10),max(nb-1,1)),
            'partial_exposure_variance':np.var(er),'collinearity_review':np.var(er)<.05,
            'design_rank':np.linalg.matrix_rank(xx),'design_columns':xx.shape[1],
            'independent_spatial_animals_inferred':False})
    return rows

def descriptive():
    curves=[];associations=[];decomposition=[]
    for stage in STAGES:
        c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');r=read(OUT/stage/'reconstructed_neighbor_levels.tsv.gz');e=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz')
        q=read(OUT/stage/'all_profile_observations.tsv.gz')
        print('RECONSTRUCT_DESCRIPTIVE',stage,flush=True)
        cohorts={'author_Figure4d_vCM':c.author_Figure4d_CM,'Figure4e_independent_mask':c.reconstructed_Figure4e_CM,'our_ventricular_mask':c.broader_ventricular_CM}
        outcomes=['native_Z__maturation21','native_logmean__maturation21','native_within_stage_Z__maturation21','native_matched3_100__maturation21','core_Z__maturation21']
        for cohort,take in cohorts.items():
            for graph in ['paper_like6NN','tissue_Delaunay']:
                for source in ['cap.EC','p.EC','Peri','cap.EC_plus_HighMyo']:
                    level=r[graph+'__'+source]
                    labels=['cap.EC','cap.EC.High.Myo.'] if source=='cap.EC_plus_HighMyo' else [source]
                    role_conflict=c.author_atlas_cell_state.isin(labels)
                    for outcome in outcomes:
                        valid=take&~role_conflict&level.notna()&c[outcome].notna();xx=level.loc[valid];yy=c.loc[valid,outcome]
                        summary=pd.DataFrame({'level':xx,'score':yy}).groupby('level').score.agg(['mean','std','size'])
                        for lev,z in summary.iterrows():curves.append({'stage':stage,'cohort':cohort,'graph':graph,'source':source,'outcome':outcome,'level':lev,
                            'mean':z['mean'],'SD':z['std'],'CMs':z['size'],'naive_cell_SE_descriptive_only':z['std']/np.sqrt(z['size'])})
                        associations.append({'stage':stage,'cohort':cohort,'graph':graph,'source':source,'outcome':outcome,'CMs':len(xx),
                            'source_role_conflicted_CMs_excluded':int((take&role_conflict).sum()),
                            'individual_Spearman':spearmanr(xx,yy).statistic,'binned_mean_Pearson':np.corrcoef(summary.index,summary['mean'])[0,1] if len(summary)>2 else np.nan,
                            'scope':'Descriptive within-section association, not independently replicated causal effect'})
        take=c.reconstructed_Figure4e_CM
        for outcome in ['native_Z__'+k for k in COMPONENTS]:
            x=e.author_capillary__nearby_cell_area_fraction_50um.loc[take];y=c.loc[take,outcome]
            d=pd.DataFrame({'x':x,'y':y,'region':c.loc[take,'anatomy_block_200um']}).dropna();means=d.groupby('region')[['x','y']].transform('mean')
            total=np.mean((d.x-d.x.mean())*(d.y-d.y.mean()));within=np.mean((d.x-means.x)*(d.y-means.y));between=np.mean((means.x-d.x.mean())*(means.y-d.y.mean()))
            decomposition.append({'stage':stage,'outcome':outcome,'CMs':len(d),'total_covariance':total,'within_200um_covariance':within,'between_200um_covariance':between,
                'between_fraction_of_total':between/total if abs(total)>1e-9 else np.nan,
                'scope':'Exact descriptive covariance partition; not fraction of biology caused by anatomy, may be negative or exceed1'})
    save(pd.DataFrame(curves),'paper_style_neighbor_level_means.tsv');save(pd.DataFrame(associations),'descriptive_proximity_associations.tsv');save(pd.DataFrame(decomposition),'shared_anatomy_covariance_partition.tsv')

def models():
    rows=[];coverage=[];prediction=[];folds=[]
    for stage in STAGES:
        print('ANATOMY_MECHANISM_MODELS',stage,flush=True)
        c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');env=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz')
        for cohort,mask in [('paper_label_mask',c.reconstructed_Figure4e_CM),('our_ventricular_mask',c.broader_ventricular_CM)]:
            for source in ['author_capillary','author_capillary_plus_HighMyo','QC_capillary']:
                if source+'__apposition_fraction_2um' not in env:continue
                ee=exposures(c,env,source);bb=covariates(c,env,True)
                valid=mask&c.core_panel_depth.ge(20)&ee.notna().all(axis=1)&bb.notna().all(axis=1)
                d=c.loc[valid];e=ee.loc[valid];en=env.loc[valid]
                coverage.append({'stage':stage,'cohort':cohort,'capillary_source':source,'CMs':len(d),'nuclear_area_missing_retained':int(d.nucleus_area_missing.sum())})
                outcomes=[v+'_Z__'+k for v in ['native','core'] for k in COMPONENTS];y=d[outcomes]
                levels=['global_joint','within_region','within_region_capture','within_region_capture_CMstate'] if source=='author_capillary' else ['within_region_capture']
                for level in levels:
                    b=pd.DataFrame(index=d.index) if level=='global_joint' else covariates(d,en,'capture' in level,'CMstate' in level)
                    for r in fit(y,b,e,d,within=level!='global_joint'):rows.append({'stage':stage,'cohort':cohort,'capillary_source':source,'adjustment':level,**r})
                if cohort=='paper_label_mask' and source=='author_capillary':
                    # Interaction models expose candidate receiver/matrix gating;
                    # receiver expression is not silently conditioned out of the
                    # primary total association.
                    z=e.copy();z['CM_receiver_Ptprm']=d.core__Ptprm
                    centered=(z-z.mean())/z.std().replace(0,1)
                    z['apposition_x_CM_receiver']=centered.capillary_apposition*centered.CM_receiver_Ptprm
                    z['apposition_x_FB_basement']=centered.capillary_apposition*centered.FB_basement_program
                    for r in fit(y,covariates(d,en,True),z,d,within=True):rows.append({'stage':stage,'cohort':cohort,'capillary_source':source,'adjustment':'within_region_capture_interactions',**r})
                    ps,fs=blocked_prediction(d,en,e);prediction.extend({'stage':stage,**r} for r in ps);folds.extend({'stage':stage,**r} for r in fs)
        # Inclusive vascular-only cohort does not require an eligible fibroblast.
        e=exposures(c,env)[['capillary_area','capillary_apposition','pericyte_context','EC_cycling','EC_Ptprm']]
        b=covariates(c,env,False)
        valid=c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20)&e.notna().all(axis=1)&b.notna().all(axis=1)
        for r in fit(c.loc[valid,['native_Z__'+k for k in COMPONENTS]],b.loc[valid],e.loc[valid],c.loc[valid],within=True):rows.append({'stage':stage,'cohort':'inclusive_vascular_no_FB_requirement','capillary_source':'author_capillary','adjustment':'within_region',**r})
    rows=pd.DataFrame(rows)
    rows['within_section_FDR_BH']=rows.groupby(['adjustment','cohort','capillary_source'],dropna=False).approximate_block_p.transform(lambda v:multipletests(v,method='fdr_bh')[1])
    save(rows,'anatomy_conditioned_coefficients.tsv.gz');save(pd.DataFrame(coverage),'mechanism_model_cohort_coverage.tsv')
    save(pd.DataFrame(prediction),'spatially_heldout_mechanism_prediction.tsv');save(pd.DataFrame(folds),'spatial_prediction_fold_checks.tsv')

def blocked_prediction(c,en,e):
    base=covariates(c,en,True);xy=c[['x_centroid','y_centroid']].to_numpy()
    spl=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(xy)
    for j in range(spl.shape[1]):base['xy_'+str(j)]=spl[:,j]
    arrays={'anatomy_capture':base.to_numpy(),'vascular_access':np.column_stack([base,e[['capillary_area','pericyte_context']]]),
        'access_plus_apposition':np.column_stack([base,e[['capillary_area','pericyte_context','capillary_apposition','EC_Ptprm']]]),
        'all_context':np.column_stack([base,e])}
    blocks=np.floor(xy/CFG['bootstrap_block_um']).astype(int);fold=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
    predictions={k:np.full((len(c),4),np.nan) for k in arrays};checks=[];outcomes=['native_Z__'+k for k in COMPONENTS[1:]];yy=c[outcomes].to_numpy()
    for j in range(5):
        test=fold==j;train=~test
        if not test.any():continue
        ti=np.flatnonzero(train);train[ti[cKDTree(xy[test]).query(xy[train])[0]<CFG['prediction_buffer_um']]]=False
        if train.sum()<200:continue
        for name,x in arrays.items():
            scale=StandardScaler().fit(x[train]);model=Ridge(alpha=10).fit(scale.transform(x[train]),yy[train]);predictions[name][test]=model.predict(scale.transform(x[test]))
        checks.append({'fold':j,'training_CMs':int(train.sum()),'heldout_CMs':int(test.sum()),'minimum_train_test_separation_um':cKDTree(xy[test]).query(xy[train])[0].min()})
    rows=[]
    for name,pred in predictions.items():
        valid=np.isfinite(pred).all(axis=1)
        for j,outcome in enumerate(outcomes):rows.append({'predictors':name,'outcome':outcome,'OOF_CMs':int(valid.sum()),'OOF_R2':r2_score(yy[valid,j],pred[valid,j]),
            'scope':'Within-section buffered spatial generalization; not held-out animals or causal influence'})
    return rows,checks

def matched_interfaces():
    pairs=[];effects=[];rng=np.random.default_rng(CFG['seed'])
    for stage in STAGES:
        c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');e=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz');valid=c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20)
        d=c.loc[valid].copy();en=e.loc[d.index];d['apposition']=en.author_capillary__apposition_fraction_2um
        d['capillary_area']=en.author_capillary__nearby_cell_area_fraction_50um
        variables=np.column_stack([np.log1p(d.native_panel_depth),np.log1p(d.cell_area),np.log1p(en.CM_native_perimeter_um),np.log1p(d.capillary_area)])
        zz=(variables-np.nanmean(variables,axis=0))/np.nanstd(variables,axis=0)
        positive=d.apposition.gt(0);cut=max(.02,d.apposition[positive].quantile(.5));high=d.apposition.ge(cut);low=d.apposition.eq(0)
        records=[]
        for _,group in d.groupby(['anatomy_block_200um','author_atlas_cell_state'],observed=True):
            hi=group.index[high.loc[group.index]];lo=list(group.index[low.loc[group.index]])
            if len(hi)==0 or len(lo)==0:continue
            for cid in rng.permutation(hi):
                if not lo:break
                i=d.index.get_loc(cid);jj=d.index.get_indexer(lo);differences=abs(zz[jj]-zz[i]);dist=(differences**2).sum(1)
                physical=np.sqrt(((d.loc[lo,['x_centroid','y_centroid']].to_numpy()-d.loc[cid,['x_centroid','y_centroid']].to_numpy())**2).sum(1))
                fail=(physical>CFG['matched_maximum_distance_um'])|(differences[:,-1]>CFG['matched_maximum_capillary_area_SD'])|(differences[:,:3].max(1)>CFG['matched_maximum_capture_shape_SD'])
                dist[fail]=np.inf
                j=int(np.argmin(dist))
                if not np.isfinite(dist[j]) or dist[j]>4:continue
                control=lo.pop(j);records.append({'stage':stage,'high_apposition_CM':cid,'zero_apposition_CM':control,'anatomy_block':d.loc[cid,'anatomy_block_200um'],
                    'author_CM_state':d.loc[cid,'author_atlas_cell_state'],'matching_distance':dist[j],'centroid_separation_um':physical[j],
                    'high_apposition_fraction':d.loc[cid,'apposition'],'control_apposition_fraction':0,
                    'high_capillary_area':d.loc[cid,'capillary_area'],'control_capillary_area':d.loc[control,'capillary_area']})
        pairs.extend(records)
        if not records:continue
        rec=pd.DataFrame(records);xy=d.loc[rec.high_apposition_CM,['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/CFG['bootstrap_block_um']).astype(int);_,inv=np.unique(blocks,axis=0,return_inverse=True);nb=inv.max()+1
        for outcome in [v+'_Z__'+k for v in ['native','core'] for k in COMPONENTS]:
            delta=d.loc[rec.high_apposition_CM,outcome].to_numpy()-d.loc[rec.zero_apposition_CM,outcome].to_numpy();sd=d[outcome].std()
            sums=np.bincount(inv,weights=delta);counts=np.bincount(inv);boot=[]
            for _ in range(CFG['bootstrap_draws']):
                sel=rng.integers(nb,size=nb);boot.append(sums[sel].sum()/counts[sel].sum()/sd)
            effects.append({'stage':stage,'outcome':outcome,'matched_pairs':len(rec),'mean_standardized_difference':delta.mean()/sd,'block_CI_low':np.quantile(boot,.025),'block_CI_high':np.quantile(boot,.975),
                'scope':'Same200um region and author CM state; capture/area/perimeter/capillary-area matching, zero controls not reused; not a causal intervention'})
    save(pd.DataFrame(pairs),'matched_apposition_CM_pairs.tsv.gz');save(pd.DataFrame(effects),'matched_apposition_component_differences.tsv')

def contact_sensitivity(only_balance=False):
    rows=[];balance=[]
    for stage in STAGES:
        print('CONTACT_DEFINITION_MODELS',stage,flush=True)
        c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');en=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz');alt=read(OUT/stage/'alternative_contact_counts.tsv.gz')
        e=exposures(c,en);b=covariates(c,en,True)
        valid=c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20)&e.notna().all(axis=1)&b.notna().all(axis=1)
        d=c.loc[valid];y=d[[v+'_Z__'+k for v in ['native','core'] for k in COMPONENTS]]
        for name in ([] if only_balance else ['native_sampled_surface_fraction','refined_boundary__capillary_contacts','native_polygon_1um__capillary_contacts']):
            ee=e.loc[valid].copy()
            if name!='native_sampled_surface_fraction':ee['capillary_apposition']=np.log1p(alt.loc[valid,name])
            for r in fit(y,b.loc[valid],ee,d,within=True):rows.append({'stage':stage,'contact_definition':name,'cohort':'Identical CMs within stage',**r})
        pairs=read(OUT/'matched_apposition_CM_pairs.tsv.gz',index_col=None);pairs=pairs[pairs.stage.eq(stage)]
        for key,values in [('log_CM_depth',np.log1p(c.native_panel_depth)),('log_CM_area',np.log1p(c.cell_area)),
                           ('log_CM_perimeter',np.log1p(en.CM_native_perimeter_um)),('log_capillary_area',np.log1p(en.author_capillary__nearby_cell_area_fraction_50um))]:
            high=values.loc[pairs.high_apposition_CM].to_numpy();low=values.loc[pairs.zero_apposition_CM].to_numpy()
            balance.append({'stage':stage,'matching_covariate':key,'pairs':len(pairs),'matched_mean_difference_SD':(high.mean()-low.mean())/values.loc[c.reconstructed_Figure4e_CM].std(),
                'scope':'Balance diagnostic, no randomized treatment or new biological replication'})
    if rows:save(pd.DataFrame(rows),'paired_contact_definition_coefficients.tsv.gz')
    save(pd.DataFrame(balance),'matched_interface_covariate_balance.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['descriptive','models','matched','contacts','balance','all'],default='all');a=p.parse_args()
    if a.step in ['descriptive','all']:descriptive()
    if a.step in ['models','all']:models()
    if a.step in ['matched','all']:matched_interfaces()
    if a.step in ['contacts','all']:contact_sensitivity()
    if a.step=='balance':contact_sensitivity(only_balance=True)

if __name__=='__main__':main()
