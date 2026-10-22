"""Library-aware matrix contrasts and spatially blocked CM association tests."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-state-mpl')
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.preprocessing import SplineTransformer, StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from fibroblast_state_common import *


def design_matrix(data, terms):
    parts=[np.ones((len(data),1))]
    for t in terms:
        if t in ['age','family']:
            parts.append(pd.get_dummies(data[t],drop_first=True,dtype=float).to_numpy())
        else:parts.append(data[[t]].to_numpy(dtype=float))
    return np.column_stack(parts)


def explained(y,x):
    fit=x@np.linalg.lstsq(x,y,rcond=None)[0]
    return float(1-np.sum((y-fit)**2)/max(np.sum((y-y.mean())**2),1e-15))


def RNA_comparisons():
    meta=pd.read_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t',index_col=0)
    meta=meta.loc[meta.fine_reference_eligible].copy()
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    nuisance=ann.groupby(['sample','fibroblast_family'],observed=True).agg(
        mean_ambient=('decontx_contamination','mean'),median_depth=('total_counts','median'))
    tables={};units=None
    tested_panels={**MATRIX_PANELS,
        'Hyaluronan_synthesis':['Has1','Has2','Has3'],
        'Hyaluronan_turnover':['Hyal1','Hyal2','Cemip'],
        'Has2':['Has2'],
        'Crosslinking_enzymes':'Lox Loxl1 Loxl2 Loxl3 Loxl4 Tgm2'.split(),
        'Collagen_processing':'Plod1 Plod2 Plod3 Serpinh1'.split(),
        'Basement_laminins':'Lama2 Lama4 Lama5 Lamb1 Lamb2 Lamc1'.split(),
        'Basement_collagen_IV':['Col4a1','Col4a2'],
        'Basement_proteoglycan_Hspg2':['Hspg2']}
    definitions=[]
    for layer in ['corrected','raw']:
        counts=pd.read_csv(OUT/f'RNA_{layer}_pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
        counts=counts.loc[:,meta.index]
        groupkey=meta['sample'].astype(str)+'::'+meta.fibroblast_family.astype(str)
        counts=counts.T.groupby(groupkey).sum().T
        info=meta.assign(unit=groupkey).groupby('unit',observed=True).agg(
            sample=('sample','first'),age=('age','first'),family=('fibroblast_family','first'),nuclei=('nuclei','sum')).loc[counts.columns]
        counts=counts.loc[:,info.nuclei.ge(30)]
        info=info.loc[counts.columns]
        cpm=counts.divide(counts.sum(axis=0),axis=1)*1e6
        tables[layer]=(counts,cpm,info)
        rows=[]
        for program, panel in tested_panels.items():
            for version in ['full','annotation_markers_withheld']:
                gl=[g for g in panel if g in cpm.index and (version=='full' or g not in ANNOTATION_GENES)]
                if not gl:continue
                if layer=='corrected':definitions.append({'program':program,'version':version,'genes':';'.join(gl),'unit':'mean log1p pseudobulk CPM'})
                scores=np.log1p(cpm.loc[gl]).mean(0)
                for unit,value in scores.items():
                    r=info.loc[unit]
                    rows.append({'unit':unit,'sample':r['sample'],'age':r.age,'family':r.family,'nuclei':r.nuclei,
                        'layer':layer,'program':program,'version':version,'score_mean_log1p_CPM':value,
                        'mean_ambient':nuisance.loc[(r['sample'],r.family),'mean_ambient'],
                        'log_median_depth':np.log1p(nuisance.loc[(r['sample'],r.family),'median_depth'])})
        frame=pd.DataFrame(rows)
        save_tsv(frame,OUT/f'RNA_{layer}_matrix_pseudobulk_scores.tsv')
    save_tsv(pd.DataFrame(definitions),OUT/'RNA_pseudobulk_program_definitions.tsv')
    allscores=pd.concat([pd.read_csv(OUT/f'RNA_{l}_matrix_pseudobulk_scores.tsv',sep='\t') for l in tables],ignore_index=True)
    fits=[];contrasts=[];changes=[]
    for (layer,program,version),part in allscores.groupby(['layer','program','version'],observed=True):
        # All units receive equal weight; results describe this observed design.
        y=part.score_mean_log1p_CPM.to_numpy()
        full=design_matrix(part,['age','family'])
        stage=design_matrix(part,['age']);state=design_matrix(part,['family'])
        covars=design_matrix(part,['age','family','mean_ambient','log_median_depth'])
        rfull=explained(y,full)
        fits.append({'layer':layer,'program':program,'version':version,'library_family_units':len(part),
            'full_design_columns':full.shape[1],'full_design_rank':np.linalg.matrix_rank(full),
            'additive_R2':rfull,'delta_R2_stage':rfull-explained(y,state),
            'delta_R2_family':rfull-explained(y,stage),'with_QC_covariates_R2':explained(y,covars),
            'interpretation':'descriptive conditional variation; not proof of state independence or causal timing'})
        for scope,selected in [('interstitial_programs_only',part.family.isin(FAMILIES[:4])),
                ('late_P14_P21_interstitial',part.family.isin(FAMILIES[:4])&part.age.isin(['P14','P21']))]:
            restricted=part.loc[selected]
            if len(restricted)<8:continue
            ry=restricted.score_mean_log1p_CPM.to_numpy()
            rx=design_matrix(restricted,['age','family'])
            rr=explained(ry,rx)
            fits.append({'layer':layer,'program':program,'version':version,'scope':scope,'library_family_units':len(restricted),
                'full_design_columns':rx.shape[1],'full_design_rank':np.linalg.matrix_rank(rx),'additive_R2':rr,
                'delta_R2_stage':rr-explained(ry,design_matrix(restricted,['family'])),
                'delta_R2_family':rr-explained(ry,design_matrix(restricted,['age'])),
                'with_QC_covariates_R2':explained(ry,design_matrix(restricted,['age','family','mean_ambient','log_median_depth'])),
                'interpretation':'restricted shared-state support; descriptive and unbalanced, no state-independence claim'})
        for stage_name,s in part.groupby('age',observed=True):
            for fam,f in s.groupby('family',observed=True):
                if len(f)<2:continue
                deltas=[]
                for lib in f['sample']:
                    own=f.loc[f['sample'].eq(lib),'score_mean_log1p_CPM'].iloc[0]
                    other=s.loc[s['sample'].eq(lib)&s.family.ne(fam),'score_mean_log1p_CPM']
                    if len(other):deltas.append(own-other.mean())
                if len(deltas)==2:
                    contrasts.append({'layer':layer,'program':program,'version':version,'stage':stage_name,'family':fam,
                        'library1_difference':deltas[0],'library2_difference':deltas[1],'mean_difference':np.mean(deltas),
                        'both_library_directions_agree':np.sign(deltas[0])==np.sign(deltas[1])})
        for fam,f in part.groupby('family',observed=True):
            for a,b in zip(STAGES[:-1],STAGES[1:]):
                av=f.loc[f.age.eq(a),'score_mean_log1p_CPM'];bv=f.loc[f.age.eq(b),'score_mean_log1p_CPM']
                if len(av)!=2 or len(bv)!=2:continue
                cross=bv.to_numpy()[:,None]-av.to_numpy()[None,:]
                changes.append({'layer':layer,'program':program,'version':version,'family':fam,'from_stage':a,'to_stage':b,
                    'mean_difference':bv.mean()-av.mean(),'all_four_library_cross_comparisons_same_sign':bool((cross>0).all() or (cross<0).all()),
                    'minimum_cross_difference':cross.min(),'maximum_cross_difference':cross.max()})
    fitframe=pd.DataFrame(fits)
    fitframe['scope']=fitframe['scope'].fillna('all_families')
    save_tsv(fitframe,OUT/'RNA_matrix_stage_family_variation.tsv')
    save_tsv(pd.DataFrame(contrasts),OUT/'RNA_matrix_within_stage_family_contrasts.tsv')
    save_tsv(pd.DataFrame(changes),OUT/'RNA_matrix_within_family_stage_changes.tsv')
    # Unbiased early Postn-family stage comparisons, with a conservative CM-enrichment exclusion.
    a=pd.read_csv(ROOT/'results/cell_identities_subtypes/Cardiomyocyte/subcluster_pseudobulk_corrected_counts.tsv.gz',sep='\t',index_col=0)
    b=pd.read_csv(ROOT/'results/cell_identities_subtypes/Fibroblast/subcluster_pseudobulk_corrected_counts.tsv.gz',sep='\t',index_col=0)
    ca=a.drop(columns='gene_symbol').sum(1);fb=b.drop(columns='gene_symbol').sum(1)
    ca=pd.Series(ca.to_numpy(),index=a.gene_symbol).groupby(level=0).sum();fb=pd.Series(fb.to_numpy(),index=b.gene_symbol).groupby(level=0).sum()
    ca=ca/ca.sum()*1e6;fb=fb/fb.sum()*1e6
    excluded=set(ca.index[(ca>=10)&((ca+1)/(fb.reindex(ca.index).fillna(0)+1)>=4)])|set(CM_CHECK)|set(SEX)
    save_tsv(pd.DataFrame({'excluded_CM_enriched_gene':sorted(excluded)}),OUT/'RNA_stage_comparison_CM_exclusions.tsv')
    counts,cpm,info=tables['corrected'];rawcounts,rawcpm,_=tables['raw']
    results=[]
    for fam in FAMILIES:
        for early,late in zip(STAGES[:-1],STAGES[1:]):
            ea=info.index[info.family.eq(fam)&info.age.eq(early)];la=info.index[info.family.eq(fam)&info.age.eq(late)]
            if len(ea)!=2 or len(la)!=2:continue
            logs=np.log2(cpm+1);rl=np.log2(rawcpm+1)
            effect=logs[la].mean(1)-logs[ea].mean(1);raw_effect=rl[la].mean(1)-rl[ea].mean(1)
            crosses=np.stack([(logs[l]-logs[e]).to_numpy() for l in la for e in ea])
            keep=(cpm[ea.tolist()+la.tolist()].max(1)>=5)&(effect.abs()>=.5)
            keep&=((crosses>.25).all(0)|(crosses<-.25).all(0))&(np.sign(effect)==np.sign(raw_effect))
            keep&=~cpm.index.isin(excluded)&~cpm.index.str.match(r'^(mt-|Rpl|Rps|Gm\d)')
            for gene in cpm.index[keep]:
                results.append({'family':fam,'from_stage':early,'to_stage':late,'gene':gene,'corrected_log2CPM_difference':effect[gene],
                    'raw_log2CPM_difference':raw_effect[gene],'minimum_library_cross_difference':crosses[:,cpm.index.get_loc(gene)].min(),
                    'maximum_library_cross_difference':crosses[:,cpm.index.get_loc(gene)].max(),
                    'interpretation':'effect-size screen across two libraries/age; no nucleus-level inferential p-value'})
    save_tsv(pd.DataFrame(results),OUT/'RNA_within_family_stage_gene_changes.tsv')


def baseline_features(d):
    direct=pd.DataFrame(index=d.index)
    for col in ['outcome_excluded_depth','cell_area','arterial_distance_um','endocardial_proxy_distance_um','resolved_neighbors']:
        direct['log_'+col]=np.log1p(d[col].to_numpy(dtype=float))
    for col in ['nuclear_area_missing','fibroblast_fraction','endothelial_fraction','mural_fraction','immune_fraction','Slit2_log','cycling_log']:
        direct[col]=d[col].to_numpy(dtype=float)
    xy=d[['x_centroid','y_centroid']].to_numpy()
    sp=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(xy)
    for j in range(sp.shape[1]):direct['spatial_spline_'+str(j)]=sp[:,j]
    scaled=(xy-xy.mean(0))/np.maximum(xy.std(0),1)
    direct['spatial_xy_interaction']=scaled[:,0]*scaled[:,1]
    return direct


def coefficient_bootstrap(y,base,exposure,blocks,seed=45):
    names=exposure.columns.tolist()
    values=np.column_stack([base.to_numpy(),exposure.to_numpy()])
    sd=values.std(0);keep=sd>1e-8
    positions=np.flatnonzero(keep)
    xx=(values[:,keep]-values[:,keep].mean(0))/sd[keep]
    yy=(y-y.mean())/max(y.std(),1e-8)
    x=np.column_stack([np.ones(len(xx)),xx])
    beta=np.linalg.pinv(x.T@x,rcond=1e-9)@(x.T@yy)
    unique,inverse=np.unique(blocks,return_inverse=True)
    h=np.zeros((len(unique),x.shape[1],x.shape[1]));g=np.zeros((len(unique),x.shape[1]))
    for k in range(len(unique)):
        part=inverse==k;h[k]=x[part].T@x[part];g[k]=x[part].T@yy[part]
    rng=np.random.default_rng(seed);draws=[]
    for _ in range(96):
        sample=np.bincount(rng.integers(0,len(unique),len(unique)),minlength=len(unique))
        gram=np.einsum('b,bij->ij',sample,h);rhs=np.einsum('b,bi->i',sample,g)
        draws.append(np.linalg.pinv(gram,rcond=1e-9)@rhs)
    draws=np.asarray(draws)
    result=[]
    for j,name in enumerate(names):
        original=base.shape[1]+j
        if original not in positions:continue
        k=list(positions).index(original)+1
        result.append({'predictor':name,'standardized_coefficient':beta[k],
            'spatial_block_bootstrap_low':np.quantile(draws[:,k],.025),'spatial_block_bootstrap_high':np.quantile(draws[:,k],.975),
            'spatial_blocks':len(unique)})
    return result


def blocked_cv(y,base,extra,d,buffer=100):
    xy=d[['x_centroid','y_centroid']].to_numpy()
    bx=np.floor(xy[:,0]/400).astype(int);by=np.floor(xy[:,1]/400).astype(int)
    folds=((bx*73856093)^(by*19349663)^271)%5
    pred={name:np.full(len(y),np.nan) for name in ['baseline','augmented']};rows=[]
    values={'baseline':base.to_numpy(),'augmented':np.column_stack([base.to_numpy(),extra.to_numpy()])}
    for fold in range(5):
        test=folds==fold;train=~test
        if not test.any():continue
        distance=cKDTree(xy[test]).query(xy[train])[0]
        ti=np.flatnonzero(train);train[ti[distance<buffer]]=False
        if train.sum()<100:continue
        for name,x in values.items():
            scaler=StandardScaler().fit(x[train]);model=Ridge(alpha=10).fit(scaler.transform(x[train]),y[train])
            p=model.predict(scaler.transform(x[test]));pred[name][test]=p
            rows.append({'fold':fold,'model':name,'training_cells':int(train.sum()),'test_cells':int(test.sum()),
                         'R2':r2_score(y[test],p),'buffer_um':buffer})
    valid=np.isfinite(pred['baseline'])&np.isfinite(pred['augmented'])
    summary={'OOF_cells':int(valid.sum()),'baseline_OOF_R2':r2_score(y[valid],pred['baseline'][valid]),
             'augmented_OOF_R2':r2_score(y[valid],pred['augmented'][valid])}
    summary['incremental_OOF_R2']=summary['augmented_OOF_R2']-summary['baseline_OOF_R2']
    return summary,rows


def spatial_comparisons():
    coefficients=[];cv=[];foldrows=[];coverage=[];pairsummaries=[];correlations=[]
    for stage in STAGES:
        print('SPATIAL_MODELS',stage,flush=True)
        folder=OUT/stage
        cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        for radius in [30,50,100]:
            env=pd.read_csv(folder/f'CM_neighborhoods_{radius}um.tsv.gz',sep='\t',index_col=0)
            d=cm.join(env)
            eligible=d.fibroblast_neighbors.ge(3)&d.resolved_neighbors.ge(5)&d.arterial_distance_um.notna()
            d=d.loc[eligible].copy()
            coverage.append({'stage':stage,'radius_um':radius,'all_strict_ventricular_CMs':len(cm),'analyzed_CMs':len(d),
                             'missing_nuclear_area_retained':int(d.nuclear_area_missing.sum())})
            base=baseline_features(d)
            families=d[[c for c in d if c.startswith('family__')]]
            # Early reference families exhaust the simplex. Omit Postn as the
            # reference component instead of fitting numerically collinear fractions.
            if families['family__Postn-associated'].std()>1e-8:
                families=families.drop(columns='family__Postn-associated')
            matrix=d[[c for c in d if c.startswith('matrix__')]]
            sets={'families':families,'matrix':matrix,'both':pd.concat([families,matrix],axis=1)}
            xy=d[['x_centroid','y_centroid']].to_numpy()
            blocks=(np.floor(xy[:,0]/400).astype(int).astype(str)+'_'+np.floor(xy[:,1]/400).astype(int).astype(str))
            for outcome in ['maturation_21gene_mean','maturation_matched_control','Ankrd1_log1pCP10k']:
                y=d[outcome].to_numpy()
                for model,exposure in sets.items():
                    if model=='both' and radius!=50:continue
                    for row in coefficient_bootstrap(y,base,exposure,blocks):
                        coefficients.append({'stage':stage,'radius_um':radius,'outcome':outcome,'model':model,'CM_profiles':len(d),**row})
                    if radius==50:
                        summary,folds=blocked_cv(y,base,exposure,d)
                        cv.append({'stage':stage,'radius_um':radius,'outcome':outcome,'model':model,**summary})
                        foldrows.extend({'stage':stage,'outcome':outcome,'feature_set':model,**r} for r in folds)
                if radius==50 and outcome=='Ankrd1_log1pCP10k':
                    # Does any association survive conditioning on panel-defined CM maturation?
                    b=base.assign(CM_maturation=d.maturation_21gene_mean.to_numpy())
                    for row in coefficient_bootstrap(y,b,sets['both'],blocks):
                        coefficients.append({'stage':stage,'radius_um':radius,'outcome':outcome,'model':'both_plus_CM_maturation','CM_profiles':len(d),**row})
            if radius==50:
                for name in ['maturation_matched_control','maturation_excluded_denominator','Ankrd1_full_denominator']:
                    reference='Ankrd1_log1pCP10k' if name.startswith('Ankrd1') else 'maturation_21gene_mean'
                    correlations.append({'stage':stage,'variant':name,'reference':reference,'correlation':d[[name,reference]].corr().iloc[0,1]})
                pairs=pd.read_csv(ROOT/'results/whole_heart_Ankrd1_vessels'/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
                pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
                pairs=pairs[pairs.near_cell_id.isin(d.index)&pairs.far_cell_id.isin(d.index)]
                for col in list(families)+list(matrix):
                    a=d.loc[pairs.near_cell_id,col].to_numpy();b=d.loc[pairs.far_cell_id,col].to_numpy()
                    pairsummaries.append({'stage':stage,'predictor':col,'matched_pairs':len(pairs),'near_mean':a.mean(),'far_mean':b.mean(),
                        'near_minus_far':np.mean(a-b),'interpretation':'frozen same-vessel CM pairs; overlapping neighborhoods possible'})
        print('SPATIAL_MODELS_DONE',stage,flush=True)
    save_tsv(pd.DataFrame(coefficients),OUT/'spatial_adjusted_coefficients.tsv')
    save_tsv(pd.DataFrame(cv),OUT/'spatial_blocked_prediction.tsv')
    save_tsv(pd.DataFrame(foldrows),OUT/'spatial_blocked_prediction_folds.tsv')
    save_tsv(pd.DataFrame(coverage),OUT/'spatial_model_coverage.tsv')
    save_tsv(pd.DataFrame(pairsummaries),OUT/'vascular_matched_CM_fibroblast_environments.tsv')
    save_tsv(pd.DataFrame(correlations),OUT/'CM_score_definition_sensitivity.tsv')
    write_json(OUT/'association_methods.json',{'RNA_unit':'Library x reassessed family,>=30 nuclei; corrected and original pseudobulk scores',
        'RNA_variance':'Equal-unit OLS, stage+family; conditional DeltaR2; rank reported; no cell-level p-values or state-independence claim',
        'RNA_gene_stage_screen':'Adjacent ages with2libraries each,family>=30/library; abs corrected log2CPM effect>=0.5,all4library crosses same direction>=0.25,raw agreement; CM-enriched/sex/mt/ribosomal/Gm genes excluded',
        'spatial_primary':'Within-section 50um FB neighborhoods; >=3FB,>=5resolved neighbors; all original broad labels frozen',
        'controls':'RNA depth excluding all25outcome genes,cell area,nuclear-area missing,arterial distance,endocardial proxy,density,broad FB/EC/mural/immune composition,Slit2/cycling RNA,XY splines+interaction',
        'coefficient_units':'Standardized outcome per SD predictor,conditional joint-family or joint-matrix model; both fit at50um',
        'uncertainty':'96resamples of400um spatial blocks; descriptive within-section intervals,not independent animals',
        'prediction':'Five deterministic400um-block folds; training CMs within100um of any held-out CM removed; Ridge alpha10; train-only feature scaling; identical folds for baseline/augmented',
        'limitations':'Unreviewed anatomical regions,cell segmentation/reference mapping residuals and spatial confounding remain; no causal signaling,stiffness,ancestry or independent replicated developmental spatial claim'})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--part',choices=['RNA','spatial','all'],default='all');args=parser.parse_args()
    if args.part in ['RNA','all']:RNA_comparisons()
    if args.part in ['spatial','all']:spatial_comparisons()
