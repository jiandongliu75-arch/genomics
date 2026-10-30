"""Common-cohort joint models with progressively explicit conditioning targets."""
from cardiac_spatial_pipeline import *
from sklearn.ensemble import HistGradientBoostingRegressor

OUTCOMES=['CM_maturation_core','CM_maturation_native','CM_category_balanced','CM_Sarcomere','CM_Metabolism','CM_Electrophysiology','CM_TFs','CM_without_Myl2']

def design(f,d,level):
    b=pd.DataFrame(index=f.index)
    for key in ['CM_outcome_excluded_depth','cell_area']:
        b['log_'+key]=np.log1p(f[key])
    b['nucleus_missing']=f.nucleus_area.isna().astype(float)
    b['log_resolved_neighbors']=np.log1p(d.resolved_neighbors)
    b['review_fraction']=d.review_fraction
    for l in ['Endothelial','Fibroblast']:
        b[l+'_coverage']=d[l+'_source_coverage'];b[l+'_capture']=d[l+'_source_capture']
    for p in PROGS:
        if 'source_present__'+p in d:b['source_present__'+p]=d['source_present__'+p]
    if level in ['anatomy','CM_state']:
        for key in ['arterial_distance_um','endocardial_proxy_distance_um']:b['log_'+key]=np.log1p(f[key])
        xy=f[['x_centroid','y_centroid']].to_numpy()
        smooth=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(xy)
        for j in range(smooth.shape[1]):b['XY_spline_'+str(j)]=smooth[:,j]
        z=(xy-xy.mean(0))/np.maximum(xy.std(0),1);b['XY_interaction']=z[:,0]*z[:,1]
    if level=='CM_state':b['CM_cycling']=f.CM_cycling;b['CM_Slit2']=f.CM_Slit2
    return b

def exposures(d,families=False):
    cols=['fraction__'+l for l in NC]+['program__'+p for p in PROGS]
    if families:cols += [c for c in d if c.startswith(('FBfamily__','ECfamily__'))]
    return d[cols].copy()

def joint_fit(y,b,e,xy,draws):
    names=list(b)+list(e);values=np.column_stack([b.to_numpy(float),e.to_numpy(float)])
    means=values.mean(0);sd=values.std(0);use=sd>1e-8
    x=np.column_stack([np.ones(len(values)),(values[:,use]-means[use])/sd[use]])
    yy=y.to_numpy(float);ym=yy.mean(0);ys=np.maximum(yy.std(0),1e-8);yy=(yy-ym)/ys
    beta=np.linalg.lstsq(x,yy,rcond=1e-9)[0]
    blocks=np.floor(xy/CFG['spatial_block_um']).astype(int);_,inv=np.unique(blocks,axis=0,return_inverse=True);nb=inv.max()+1
    h=np.zeros((nb,x.shape[1],x.shape[1]));r=np.zeros((nb,x.shape[1],yy.shape[1]))
    for k in range(nb):
        take=inv==k;h[k]=x[take].T@x[take];r[k]=x[take].T@yy[take]
    rng=np.random.default_rng(CFG['seed']);samples=[]
    for _ in range(draws):
        w=np.bincount(rng.integers(nb,size=nb),minlength=nb)
        samples.append(np.linalg.pinv(np.einsum('b,bij->ij',w,h),rcond=1e-9)@np.einsum('b,big->ig',w,r))
    samples=np.asarray(samples);low,high=np.quantile(samples,[.025,.975],axis=0)
    selected=np.flatnonzero(use);rows=[]
    for j,name in enumerate(names):
        if j<len(b.columns) or j not in selected:continue
        k=int(np.flatnonzero(selected==j)[0])+1
        # A partial variance diagnostic marks nearly redundant exposure terms.
        keep=np.arange(x.shape[1])!=k
        er=x[:,k]-x[:,keep]@np.linalg.lstsq(x[:,keep],x[:,k],rcond=1e-9)[0]
        partial_var=np.var(er)
        for out,outcome in enumerate(y):
            rows.append({'exposure':name,'outcome':outcome,'standardized_coefficient':beta[k,out],
                'block_CI_low':low[k,out],'block_CI_high':high[k,out],'CMs':len(x),'spatial_blocks':nb,'bootstrap_draws':draws,
                'nonzero_exposure_CMs':int((e[name].abs()>1e-8).sum()),'partial_exposure_variance':partial_var,
                'collinearity_review':partial_var<.05,'design_rank':np.linalg.matrix_rank(x),'design_columns':x.shape[1],
                'conditional_exposure_SD':sd[j], 'raw_unit_coefficient':beta[k,out]*ys[out]/sd[j],
                'unadjusted_Spearman':spearmanr(y[outcome],e[name]).statistic})
    return rows

def blocked_prediction(y,b,e,xy):
    blocks=np.floor(xy/CFG['spatial_block_um']).astype(int)
    folds=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
    pred={name:np.full(len(y),np.nan) for name in ['baseline','composition','programs','joint']}
    comp=e[[c for c in e if c.startswith('fraction__')]].to_numpy()
    prog=e[[c for c in e if c.startswith('program__')]].to_numpy()
    arrays={'baseline':b.to_numpy(),'composition':np.column_stack([b,comp]),'programs':np.column_stack([b,prog]),'joint':np.column_stack([b,e])}
    rows=[]
    for fold in range(5):
        test=folds==fold;train=~test
        if not test.any():continue
        pos=np.flatnonzero(train);train[pos[cKDTree(xy[test]).query(xy[train])[0]<CFG['prediction_buffer_um']]]=False
        if train.sum()<100:continue
        for name,v in arrays.items():
            scaler=StandardScaler().fit(v[train]);model=Ridge(alpha=10).fit(scaler.transform(v[train]),y[train])
            pred[name][test]=model.predict(scaler.transform(v[test]))
        rows.append({'fold':fold,'train_CMs':int(train.sum()),'test_CMs':int(test.sum()),'minimum_train_test_separation_um':float(cKDTree(xy[test]).query(xy[train])[0].min())})
    valid=np.logical_and.reduce([np.isfinite(p) for p in pred.values()])
    summaries=[]
    for name,p in pred.items():
        summaries.append({'predictors':name,'OOF_CMs':int(valid.sum()),'OOF_R2':r2_score(y[valid],p[valid]) if valid.any() else np.nan,
            'OOF_mean_absolute_error':np.mean(np.abs(y[valid]-p[valid])) if valid.any() else np.nan})
    return summaries,rows

def models():
    coefficients=[];predictions=[];folds=[];coverage=[];dropped=[]
    for stage in STAGES:
        f=frame(stage);cm=f.loc[f.ventricular_CM_strict]
        tables={g:read(f"{CFG['output']}/{stage}/{g}_CM_neighborhoods.tsv.gz") for g in ['proximity50','contact','proximity30','proximity100']}
        for geometry in ['proximity50','contact','proximity30','proximity100','paired_proximity50','paired_contact']:
            g=geometry.replace('paired_','');d=tables[g]
            mask=d.geometry_eligible.copy()
            if geometry.startswith('paired_'):mask &= d.common_contact_50um_cohort
            data=d.loc[mask];c=cm.loc[data.index]
            print('JOINT_MODELS',stage,geometry,len(data),flush=True)
            coverage.append({'stage':stage,'geometry':geometry,'CMs':len(data),'nuclear_area_missing_retained':int(c.nucleus_area.isna().sum())})
            if len(data)<200:continue
            ee=exposures(data)
            # Exposure-only criteria; no maturation values select predictors.
            bad=[key for key in ee if ee[key].std()<1e-8 or (ee[key].abs()>1e-8).sum()<50]
            for key in bad:dropped.append({'stage':stage,'geometry':geometry,'exposure':key,'reason':'constant or fewer than50 nonzero CM neighborhoods'})
            ee=ee.drop(columns=bad)
            xy=c[['x_centroid','y_centroid']].to_numpy()
            for level in ['technical','anatomy','CM_state']:
                b=design(c,data,level)
                draws=CFG['bootstrap_primary'] if geometry=='proximity50' else CFG['bootstrap_sensitivity']
                for r in joint_fit(c[OUTCOMES],b,ee,xy,draws):coefficients.append({'stage':stage,'geometry':geometry,'adjustment':level,**r})
                if geometry=='proximity50':
                    ps,fs=blocked_prediction(c.CM_maturation_core.to_numpy(),b,ee,xy)
                    predictions.extend({'stage':stage,'adjustment':level,**r} for r in ps)
                    folds.extend({'stage':stage,'adjustment':level,**r} for r in fs)
            if geometry=='proximity50':
                other=exposures(data,families=True)
                extra=[c for c in other if c not in ee and (other[c].abs()>1e-8).sum()>=50 and other[c].std()>1e-8]
                # Only family variables are added; previously dropped core exposures stay dropped.
                extra=[key for key in extra if 'family__' in key]
                if extra:
                    e2=pd.concat([ee,other[extra]],axis=1)
                    for r in joint_fit(c[OUTCOMES],design(c,data,'anatomy'),e2,xy,CFG['bootstrap_sensitivity']):
                        coefficients.append({'stage':stage,'geometry':'proximity50','adjustment':'anatomy_plus_supported_families',**r})
        save(pd.DataFrame(coefficients),'joint_association_coefficients.tsv.gz')
        save(pd.DataFrame(predictions),'spatial_prediction.tsv');save(pd.DataFrame(folds),'spatial_prediction_folds.tsv')
        save(pd.DataFrame(coverage),'model_cohort_coverage.tsv');save(pd.DataFrame(dropped,columns=['stage','geometry','exposure','reason']),'unsupported_exposures.tsv')
    shortlist(pd.DataFrame(coefficients))

def shortlist(a):
    main=a.loc[a.geometry.eq('proximity50')&a.adjustment.eq('anatomy')&a.outcome.eq('CM_maturation_core')].copy()
    rows=[]
    for _,r in main.iterrows():
        sign=np.sign(r.standardized_coefficient)
        scorecheck=a.loc[a.stage.eq(r.stage)&a.geometry.eq('proximity50')&a.adjustment.eq('anatomy')&a.exposure.eq(r.exposure)&
                         a.outcome.isin(['CM_maturation_native','CM_category_balanced','CM_without_Myl2'])]
        boundary_agreement=bool(len(scorecheck)==3 and (np.sign(scorecheck.standardized_coefficient)==sign).all())
        paired=a.loc[a.stage.eq(r.stage)&a.geometry.isin(['paired_proximity50','paired_contact'])&a.adjustment.eq('anatomy')&
                     a.exposure.eq(r.exposure)&a.outcome.eq('CM_maturation_core')]
        contact_agreement=bool(len(paired)==2 and (np.sign(paired.standardized_coefficient)==sign).all())
        technical=a.loc[a.stage.eq(r.stage)&a.geometry.eq('proximity50')&a.adjustment.eq('technical')&a.exposure.eq(r.exposure)&a.outcome.eq('CM_maturation_core')]
        basic_agrees=bool(len(technical)==1 and np.sign(technical.iloc[0].standardized_coefficient)==sign)
        interval_direction=bool(r.block_CI_low>0 or r.block_CI_high<0)
        priority=boundary_agreement and contact_agreement and basic_agrees and interval_direction and not r.collinearity_review and abs(r.standardized_coefficient)>=.05
        rows.append({**r.to_dict(),'native_balanced_Myl2withheld_sign_agreement':boundary_agreement,
            'paired_contact_sign_agreement':contact_agreement,'technical_adjustment_sign_agreement':basic_agrees,
            'operational_effect_size_threshold_SD':.05,'priority_within_section_candidate':priority,
            'scope':'Exploratory robustness screen, not causal evidence, multiplicity-corrected significance or replicate validation'})
    save(pd.DataFrame(rows),'candidate_relationships.tsv')

def mechanisms():
    """RNA pseudobulk support for curated candidate pathways, not CellChat output."""
    pairs=[('Nrg1','Erbb4','neuregulin receptor component; co-receptors not fully tested'),
           ('Igf1','Igf1r','growth/trophic candidate'),('Igf2','Igf1r','growth/trophic candidate'),
           ('Tgfb1','Tgfbr1','TGF-beta receptor component; requires receptor complex'),
           ('Tgfb2','Tgfbr2','TGF-beta receptor component; requires receptor complex'),
           ('Fgf1','Fgfr1','FGF receptor component; co-factors not fully tested'),('Fgf2','Fgfr1','FGF receptor component; co-factors not fully tested'),
           ('Fgf9','Fgfr2','FGF receptor candidate'),('Bmp6','Bmpr1a','BMP receptor component; requires receptor complex'),
           ('Hbegf','Egfr','EGF-family receptor candidate'),('Jag1','Notch1','contact-associated candidate; not proof of receptor activation'),
           ('Dll4','Notch1','contact-associated candidate; not proof of receptor activation'),('Gas6','Axl','TAM receptor candidate')]
    meta=read('results/xenium_ventricular_mapping/reference_metadata.tsv.gz')
    with np.load(ROOT/'results/xenium_ventricular_mapping/reference_shared_panel_counts.npz') as b:
        genes=b['genes'].astype(str);ids=b['cell_ids'].astype(str);x=b['corrected']
    assert np.array_equal(ids,meta.index)
    needed=sorted(set(sum(([l,r] for l,r,_ in pairs),[])));positions=pd.Index(genes).get_indexer(needed)
    assert (positions>=0).all()
    rows=[]
    for stage in STAGES:
        libraries=sorted(meta.loc[meta.age.eq(stage),'sample'].unique())
        cm=meta.reviewed_identity.eq('Cardiomyocyte')&meta.proposed_subtype.isin(['Ventricular-like','Conduction-associated CM candidate'])
        for sender in NC:
            vals=[]
            for lib in libraries:
                s=meta['sample'].eq(lib)&meta.reviewed_identity.eq(sender);c=meta['sample'].eq(lib)&cm
                st=x[s.to_numpy()];ct=x[c.to_numpy()]
                if s.sum()<10 or c.sum()<10:continue
                # Shared381-target pseudobulk CPM, explicitly not whole-transcriptome CPM.
                sv=np.log1p(st[:,positions].sum(0)/max(st.sum(),1)*1e6)
                cv=np.log1p(ct[:,positions].sum(0)/max(ct.sum(),1)*1e6)
                vals.append((lib,s.sum(),c.sum(),np.asarray(sv).ravel(),np.asarray(cv).ravel()))
            for ligand,receptor,note in pairs:
                li=needed.index(ligand);ri=needed.index(receptor)
                enough=len(vals)==2
                support=enough and all(v[3][li]>np.log1p(5) and v[4][ri]>np.log1p(5) for v in vals)
                rows.append({'stage':stage,'sender':sender,'receiver':'Ventricular/conduction-associated CM RNA reference',
                    'ligand':ligand,'receptor_component':receptor,'both_measured_in_Xenium':True,
                    'RNA_libraries_evaluable':len(vals),'both_RNA_libraries_expression_support':support,
                    'minimum_sender_log1p_sharedpanelCPM':min([v[3][li] for v in vals],default=np.nan),
                    'minimum_receiver_log1p_sharedpanelCPM':min([v[4][ri] for v in vals],default=np.nan),
                    'interpretation':note,'scope':'Curated expression-supported hypothesis, not interaction probability, causal influence or validated receptor complex'})
    save(pd.DataFrame(rows),'candidate_signaling_expression.tsv')

def nonlinear_prediction():
    """Fixed-hyperparameter check that linear prediction misses interactions."""
    rows=[];audit=[]
    for stage in STAGES:
        print('NONLINEAR_PREDICTION',stage,flush=True)
        f=frame(stage);d=read(f"{CFG['output']}/{stage}/proximity50_CM_neighborhoods.tsv.gz")
        d=d.loc[d.geometry_eligible];c=f.loc[d.index];b=design(c,d,'anatomy');e=exposures(d)
        e=e[[key for key in e if e[key].std()>1e-8 and (e[key].abs()>1e-8).sum()>=50]]
        y=c.CM_maturation_core.to_numpy();xy=c[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/400).astype(int)
        folds=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
        arrays={'baseline':b.to_numpy(),'joint':np.column_stack([b,e])};pred={k:np.full(len(c),np.nan) for k in arrays}
        for fold in range(5):
            test=folds==fold;train=~test
            if not test.any():continue
            pos=np.flatnonzero(train);dist=cKDTree(xy[test]).query(xy[train])[0]
            train[pos[dist<CFG['prediction_buffer_um']]]=False
            if train.sum()<100:continue
            for name,x in arrays.items():
                model=HistGradientBoostingRegressor(max_iter=150,max_depth=3,max_leaf_nodes=8,learning_rate=.05,
                    min_samples_leaf=100,l2_regularization=10,early_stopping=False,random_state=CFG['seed']).fit(x[train],y[train])
                pred[name][test]=model.predict(x[test])
            audit.append({'stage':stage,'fold':fold,'training_CMs':int(train.sum()),'test_CMs':int(test.sum()),
                'minimum_train_test_separation_um':cKDTree(xy[test]).query(xy[train])[0].min()})
        valid=np.isfinite(pred['baseline'])&np.isfinite(pred['joint'])
        for name,p in pred.items():rows.append({'stage':stage,'adjustment':'anatomy','predictors':name,'OOF_CMs':int(valid.sum()),
            'OOF_R2':r2_score(y[valid],p[valid]),'OOF_MAE':np.mean(abs(y[valid]-p[valid])),
            'algorithm':'HistGradientBoosting; fixed150 iterations,depth3,leaf100,L2=10; no target-based tuning',
            'scope':'Nonlinear predictive sensitivity, not effect estimates or causal feature importance'})
    save(pd.DataFrame(rows),'nonlinear_prediction.tsv');save(pd.DataFrame(audit),'nonlinear_prediction_folds.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['models','mechanisms','nonlinear','all'],default='all');args=p.parse_args()
    if args.step in ['models','all']:models()
    if args.step in ['mechanisms','all']:mechanisms()
    if args.step in ['nonlinear','all']:nonlinear_prediction()

if __name__=='__main__':main()
