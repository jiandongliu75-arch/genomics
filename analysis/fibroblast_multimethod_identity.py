"""Compare genuinely different FB identity-transfer models; no lineage claims.

All models use only shared measured genes and training-reference transforms.
No query-wise gene centering, spatial smoothing or target label proportions.
Reference-library tests and naming-marker withholding precede spatial inference.
"""
from pathlib import Path
import os,argparse,json
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.special import softmax
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_fscore_support,balanced_accuracy_score
from fibroblast_state_common import FAMILIES,WITHHOLD,CM_CHECK,SEX,CYCLING,ANCHORS

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/fibroblast_multimethod_identity'
STAGES=['P0','P7','P14','P21'];METHODS=['balanced_logistic','rank_centroid','negative_binomial']
SEED=261026;INPUTS=set()
EXCLUDE=set(WITHHOLD+CM_CHECK+SEX+CYCLING)
NAMING=set(sum(ANCHORS.values(),[]))
def read(path,**kw):
    p=Path(path);INPUTS.add(str(p.relative_to(ROOT)))
    return pd.read_csv(p,sep='\t',**kw)
def save(d,name):
    OUT.mkdir(exist_ok=True,parents=True);d.to_csv(OUT/name,sep='\t',index=False,na_rep='NA',float_format='%.8g')
def refs():
    ann=read(ROOT/'results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz',index_col=0)
    p=ROOT/'results/xenium_ventricular_mapping/reference_shared_panel_counts.npz';INPUTS.add(str(p.relative_to(ROOT)))
    with np.load(p) as b:
        genes=b['genes'].astype(str);ids=b['cell_ids'].astype(str);ix=pd.Index(ids).get_indexer(ann.index);assert (ix>=0).all()
        x=b['corrected'][ix].astype(float);raw=b['original'][ix].astype(float)
    keep=ann.fibroblast_family.isin(FAMILIES)
    return ann[keep],x[keep],raw[keep],genes
def lognorm(x):return np.log1p(1e4*x/np.maximum(x.sum(1),1)[:,None])
def ranks(x):return rankdata(x,axis=1,method='average')/x.shape[1]
def normalize(x):return x/np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-10)

class Classifiers:
    def __init__(self,x,y,libraries,genes,withhold=False):
        self.classes=np.array([c for c in FAMILIES if np.any(y==c)])
        excluded=EXCLUDE|(NAMING if withhold else set())
        self.use=np.flatnonzero(~np.isin(genes,list(excluded))&((x>0).sum(0)>=max(5,int(.005*len(x)))))
        x=x[:,self.use];self.depth=np.median(x.sum(1));self.genes=genes[self.use]
        z=lognorm(x);self.mu=z.mean(0);self.sd=np.maximum(z.std(0),.25)
        self.logistic=LogisticRegression(C=.3,class_weight='balanced',solver='lbfgs',max_iter=500,random_state=SEED)
        self.logistic.fit(np.clip((z-self.mu)/self.sd,-5,5),y)
        r=ranks(x);templates=[];rankmeans=[]
        for c in self.classes:
            # Each available RNA library contributes equally within a family.
            p=[];q=[]
            for lib in np.unique(libraries[y==c]):
                m=(y==c)&(libraries==lib)
                p.append((x[m]/np.maximum(x[m].sum(1),1)[:,None]).mean(0));q.append(r[m].mean(0))
            templates.append(np.mean(p,axis=0));rankmeans.append(np.mean(q,axis=0))
        self.rank_templates=normalize(np.array(rankmeans)-.5)
        self.templates=np.maximum(np.array(templates),1e-6);self.templates/=self.templates.sum(1)[:,None]
        normalized_counts=x/np.maximum(x.sum(1),1)[:,None]*self.depth
        means=normalized_counts.mean(0);variances=normalized_counts.var(0)
        self.theta=np.clip(means**2/np.maximum(variances-means,1e-5),.5,100)
        self.empirical_priors=np.array([(y==c).mean() for c in self.classes])
    def predict(self,x):
        x=np.asarray(x[:,self.use],float);n=len(x);out={};scores={}
        z=np.clip((lognorm(x)-self.mu)/self.sd,-5,5)
        p=self.logistic.predict_proba(z)
        out['balanced_logistic']=p[:,pd.Index(self.logistic.classes_).get_indexer(self.classes)]
        r=normalize(ranks(x)-.5);cos=r@self.rank_templates.T
        out['rank_centroid']=softmax(cos/.02,axis=1)
        depth=np.maximum(x.sum(1),1);ll=np.zeros((n,len(self.classes)))
        for j,template in enumerate(self.templates):
            for lo in range(0,n,1500):
                hi=min(n,lo+1500);mu=np.maximum(depth[lo:hi,None]*template[None,:],1e-10)
                ll[lo:hi,j]=(x[lo:hi]*np.log(mu)-(x[lo:hi]+self.theta)*np.log(mu+self.theta)).sum(1)
        out['negative_binomial']=softmax(ll/np.sqrt(depth[:,None]),axis=1)
        for method,p in out.items():
            order=np.argsort(p,axis=1);top=order[:,-1]
            scores[method]={'label':self.classes[top],'confidence':p[np.arange(n),top],
                'margin':p[np.arange(n),top]-p[np.arange(n),order[:,-2]],'probabilities':p,
                'empirical_prior_label':self.classes[(np.log(np.maximum(p,1e-15))+np.log(self.empirical_priors)).argmax(1)]}
        return scores

def metrics(true,pred,stage,variant,method,library,condition):
    classes=np.array([c for c in FAMILIES if np.any(true==c)])
    pr,re,f1,n=precision_recall_fscore_support(true,pred,labels=classes,zero_division=0)
    rows=[{'stage':stage,'variant':variant,'method':method,'heldout_library':library,'condition':condition,'family':c,
           'precision':a,'recall':b,'F1':d,'reference_nuclei':int(k),
           'scope':'Opposite-library RNA reference recovery; not independent Xenium truth'} for c,a,b,d,k in zip(classes,pr,re,f1,n)]
    return rows
def calibrate(predictions):
    rows=[]
    d=pd.DataFrame(predictions)
    for (stage,variant,method,family),g in d.groupby(['stage','variant','method','prediction']):
        # Reference-only operational threshold; no Xenium labels used.
        q=g[g.condition.eq('observed')];candidates=np.unique(np.r_[0,np.quantile(q.confidence,np.linspace(0,1,41))])
        threshold=np.inf;support=0;precision=np.nan
        for t in candidates:
            z=q[q.confidence.ge(t)];p=z.correct.mean()
            if len(z)>=20 and p>=.8 and z.heldout_library.nunique()==2:
                threshold=t;support=len(z);precision=p;break
        rows.append({'stage':stage,'variant':variant,'method':method,'family':family,'reference_confidence_threshold':threshold,
            'reference_calibration_predictions':support,'reference_calibration_precision':precision,
            'scope':'Selected on reference CV predictions; not prospective accuracy or calibrated Xenium probability'})
    return pd.DataFrame(rows)

def benchmark():
    ann,x,raw,genes=refs();results=[];predictions=[];definitions=[]
    for stage in STAGES:
        take=ann.age.eq(stage).to_numpy();a=ann.loc[take];xx=x[take];yy=a.fibroblast_family.to_numpy();libs=a['sample'].to_numpy()
        for variant in ['shared_panel','naming_markers_withheld']:
            for i,lib in enumerate(np.unique(libs)):
                print('FB_MODEL_BENCHMARK',stage,variant,lib,flush=True)
                train=libs!=lib;test=~train;model=Classifiers(xx[train],yy[train],libs[train],genes,variant!='shared_panel')
                rng=np.random.default_rng(SEED+17*i+len(xx));q=xx[test];depth=q.sum(1)
                thinned=rng.poisson(q*np.minimum(1,60/np.maximum(depth,1))[:,None]).astype(float)
                dropout=q.copy();dropout[rng.random(dropout.shape)<.15]=0
                for condition,values in [('observed',q),('Poisson_thinned_to60',thinned),('15percent_gene_dropout',dropout)]:
                    outputs=model.predict(values)
                    for method,d in outputs.items():
                        results.extend(metrics(yy[test],d['label'],stage,variant,method,lib,condition))
                        if condition=='observed':
                            predictions.extend({'stage':stage,'variant':variant,'method':method,'heldout_library':lib,'condition':condition,
                                'cell_id':cid,'true_family':true,'prediction':p,'confidence':float(c),'margin':float(m),'correct':bool(p==true)}
                                for cid,true,p,c,m in zip(a.index[test],yy[test],d['label'],d['confidence'],d['margin']))
                definitions.extend({'stage':stage,'variant':variant,'heldout_library':lib,'gene':g} for g in model.genes)
    # Late pooled references test dependence on the age-specific family template.
    # Entire replicate-number groups are held out, never nuclei randomly split.
    late=ann.age.isin(['P14','P21']).to_numpy();a=ann.loc[late];xx=x[late];yy=a.fibroblast_family.to_numpy();libs=a['sample'].to_numpy()
    replicates=a.replicate.to_numpy()
    for replicate in np.unique(replicates):
        train=replicates!=replicate;test=~train
        model=Classifiers(xx[train],yy[train],libs[train],genes,False)
        for stage in ['P14','P21']:
            select=test&a.age.eq(stage).to_numpy();lib=str(a.loc[select,'sample'].iloc[0]);q=xx[select]
            print('FB_POOLED_REFERENCE_BENCHMARK',stage,lib,flush=True)
            rng=np.random.default_rng(SEED+int(replicate));thin=rng.poisson(q*np.minimum(1,60/np.maximum(q.sum(1),1))[:,None]).astype(float)
            drop=q.copy();drop[rng.random(drop.shape)<.15]=0
            for condition,values in [('observed',q),('Poisson_thinned_to60',thin),('15percent_gene_dropout',drop)]:
                for method,d in model.predict(values).items():
                    results.extend(metrics(yy[select],d['label'],stage,'late_pooled_reference',method,lib,condition))
                    if condition=='observed':predictions.extend({'stage':stage,'variant':'late_pooled_reference','method':method,'heldout_library':lib,'condition':condition,
                        'cell_id':cid,'true_family':true,'prediction':p,'confidence':float(c),'margin':float(m),'correct':bool(p==true)}
                        for cid,true,p,c,m in zip(a.index[select],yy[select],d['label'],d['confidence'],d['margin']))
            definitions.extend({'stage':stage,'variant':'late_pooled_reference','heldout_library':lib,'gene':g} for g in model.genes)
    save(pd.DataFrame(results),'RNA_opposite_library_model_benchmark.tsv')
    save(pd.DataFrame(predictions),'RNA_opposite_library_predictions.tsv.gz')
    save(calibrate(predictions),'RNA_reference_confidence_screens.tsv')
    save(pd.DataFrame(definitions),'RNA_fold_training_genes.tsv')
    (OUT/'benchmark_input_paths.json').write_text(json.dumps(sorted(INPUTS),indent=2)+'\n')

def queries(stage,genes):
    atlas=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',index_col=0)
    ids=atlas.index[atlas.atlas_identity.eq('Fibroblast')]
    native=ROOT/'results/xenium_ventricular_mapping'/stage/(stage+'_xenium_counts_and_candidate_mapping.h5ad');INPUTS.add(str(native.relative_to(ROOT)))
    b=ad.read_h5ad(native,backed='r');ix=b.obs_names.get_indexer(ids);gi=pd.Index(b.var_names).get_indexer(genes);assert (ix>=0).all() and (gi>=0).all()
    values={'native':b.X[ix][:,gi].toarray().astype(float)};b.file.close()
    path=ROOT/'results/xenium_nuclear_anchored_segmentation'/stage/'nuclear_anchored_counts.h5ad';INPUTS.add(str(path.relative_to(ROOT)))
    b=ad.read_h5ad(path,backed='r');ix=b.obs_names.get_indexer(ids);gi=pd.Index(b.var_names).get_indexer(genes)
    values['refined']=b.X[ix][:,gi].toarray().astype(float)
    values['core']=b.layers['high_confidence'][ix][:,gi].toarray().astype(float)
    values['nuclear']=b.layers['nuclear_only'][ix][:,gi].toarray().astype(float)
    b.file.close()
    return atlas.loc[ids],values

def consensus(labels,supported=None):
    labels=np.column_stack(labels);n=len(labels);out=np.full(n,'Unresolved',object);votes=np.zeros(n,int)
    for family in FAMILIES:
        v=(labels==family)
        if supported is not None:v&=np.column_stack(supported)
        counts=v.sum(1);keep=counts>=2;out[keep]=family;votes[keep]=counts[keep]
    return out,votes

def mapping():
    ann,x,raw,genes=refs();screens=read(OUT/'RNA_reference_confidence_screens.tsv')
    summary=[];depthchecks=[]
    for stage in STAGES:
        print('FB_MODEL_MAP',stage,flush=True);take=ann.age.eq(stage).to_numpy();a=ann.loc[take];xx=x[take];rr=raw[take];y=a.fibroblast_family.to_numpy();libs=a['sample'].to_numpy()
        atlas,q=queries(stage,genes);outputs={};primary_model=None
        variants=['shared_panel','naming_markers_withheld']+(['late_pooled_reference'] if stage in ['P14','P21'] else [])
        shape=(len(variants),len(q),len(METHODS),len(atlas))
        compact={k:np.zeros(shape,dtype=dtype) for k,dtype in [('label_code',np.uint8),('confidence',np.float32),('margin',np.float32),
                    ('reference_screen_pass',bool),('empirical_prior_label_code',np.uint8),('informative',bool)]}
        for variant in variants:
            if variant=='late_pooled_reference':
                late=ann.age.isin(['P14','P21']).to_numpy();model=Classifiers(x[late],ann.loc[late,'fibroblast_family'].to_numpy(),ann.loc[late,'sample'].to_numpy(),genes,False)
            else:model=Classifiers(xx,y,libs,genes,variant=='naming_markers_withheld')
            if variant=='shared_panel':primary_model=model
            for view,counts in q.items():
                pred=model.predict(counts);outputs[(variant,view)]=pred
                for method,r in pred.items():
                    s=screens[screens.stage.eq(stage)&screens.variant.eq(variant)&screens.method.eq(method)].set_index('family')
                    threshold=np.array([s.loc[t,'reference_confidence_threshold'] if t in s.index else np.inf for t in r['label']])
                    informative=(counts[:,model.use].sum(1)>=20)&((counts[:,model.use]>0).sum(1)>=8)
                    r['supported']=(r['confidence']>=threshold)&informative
                    at=(variants.index(variant),list(q).index(view),METHODS.index(method))
                    compact['label_code'][at]=pd.Index(FAMILIES).get_indexer(r['label']).astype(np.uint8)
                    compact['confidence'][at]=r['confidence'];compact['margin'][at]=r['margin']
                    compact['reference_screen_pass'][at]=r['supported'];compact['informative'][at]=informative
                    compact['empirical_prior_label_code'][at]=pd.Index(FAMILIES).get_indexer(r['empirical_prior_label']).astype(np.uint8)
        np.savez_compressed(OUT/(stage+'_model_predictions.npz'),**compact,cell_ids=atlas.index.to_numpy(str),
                            variants=np.array(variants),views=np.array(list(q)),methods=np.array(METHODS),families=np.array(FAMILIES))
        # Raw-reference sensitivity is reported, not imposed as an extra veto.
        rawmodel=Classifiers(rr,y,libs,genes,False);rawpred=rawmodel.predict(q['refined'])
        d=atlas[['stage','x_centroid','y_centroid','program_source_eligible','author_atlas_cell_state','fine_identity']].copy()
        for variant in variants:
            for view in q:
                pred=outputs[(variant,view)];lab,votes=consensus([pred[m]['label'] for m in METHODS])
                support,sv=consensus([pred[m]['label'] for m in METHODS],[pred[m]['supported'] for m in METHODS])
                d[variant+'__'+view+'__majority']=lab;d[variant+'__'+view+'__reference_supported']=support
                for family in FAMILIES:
                    eligible=d.program_source_eligible
                    summary.append({'stage':stage,'variant':variant,'view':view,'method':'majority_vote','family':family,
                        'eligible_FB_profiles':int(eligible.sum()),'proposal_profiles':int((eligible&pd.Series(lab,index=d.index).eq(family)).sum()),
                        'reference_screened_profiles':int((eligible&pd.Series(support,index=d.index).eq(family)).sum())})
                for method,r in pred.items():
                    d[variant+'__'+view+'__'+method]=r['label']
                    if variant=='shared_panel' and view=='refined':d[method+'__raw_reference_agrees']=rawpred[method]['label']==r['label']
                    for family in FAMILIES:
                        eligible=d.program_source_eligible.to_numpy();match=r['label']==family
                        summary.append({'stage':stage,'variant':variant,'view':view,'method':method,'family':family,
                            'eligible_FB_profiles':int(eligible.sum()),'proposal_profiles':int((eligible&match).sum()),
                            'reference_screened_profiles':int((eligible&match&r['supported']).sum())})
        refined=d.shared_panel__refined__reference_supported
        core=d.shared_panel__core__reference_supported
        d['cell_supported_candidate']=refined.where(d.program_source_eligible,'Unresolved')
        d['core_counterevidence']=core.ne('Unresolved')&refined.ne('Unresolved')&core.ne(refined)
        d['naming_marker_withheld_agrees']=d.shared_panel__refined__majority.eq(d.naming_markers_withheld__refined__majority)&d.shared_panel__refined__majority.ne('Unresolved')
        nuclear=d.shared_panel__nuclear__reference_supported
        d['nuclear_supports_cell_candidate']=nuclear.eq(refined)&refined.ne('Unresolved')
        d['nuclear_supported_alternative']=nuclear.ne('Unresolved')&refined.ne('Unresolved')&nuclear.ne(refined)
        d['evidence_tier']='Unresolved or method disagreement'
        valid=d.cell_supported_candidate.ne('Unresolved')
        d.loc[valid,'evidence_tier']='Cell-supported; core inconclusive'
        d.loc[valid&refined.eq(core),'evidence_tier']='Cell and core supported'
        d.loc[valid&d.core_counterevidence,'evidence_tier']='Cell/core disagreement; review'
        d['recommended_candidate']=d.cell_supported_candidate.where(~d.core_counterevidence,'Unresolved')
        # Same-profile depth-matched thinning distinguishes RNA loss from geometry.
        rng=np.random.default_rng(SEED+int(stage[1:]));ratios=np.minimum(1,q['core'].sum(1)/np.maximum(q['refined'].sum(1),1))
        full=d.shared_panel__refined__majority.to_numpy();coremajor=d.shared_panel__core__majority.to_numpy()
        for repeat in range(3):
            thin=rng.binomial(q['refined'].astype(np.int64),ratios[:,None]).astype(float)
            r=primary_model.predict(thin);thinmajor,_=consensus([r[m]['label'] for m in METHODS])
            for family in FAMILIES:
                use=d.program_source_eligible.to_numpy()&(full==family)
                depthchecks.append({'stage':stage,'family':family,'repeat':repeat,'refined_majority_profiles':int(use.sum()),
                    'core_majority_agreement':float((coremajor[use]==full[use]).mean()) if use.any() else np.nan,
                    'depth_matched_majority_agreement':float((thinmajor[use]==full[use]).mean()) if use.any() else np.nan,
                    'scope':'Same cells and retained-depth expectation; random thinning is not a replicate or a new segmentation'})
        save(d.reset_index(),stage+'_multimethod_FB_annotations.tsv.gz')
    save(pd.DataFrame(summary),'assignment_counts_by_model_and_view.tsv')
    save(pd.DataFrame(depthchecks),'core_vs_depth_matched_assignment_stability.tsv')
    path=OUT/'benchmark_input_paths.json';inputs=set(json.loads(path.read_text()))|INPUTS
    path.write_text(json.dumps(sorted(inputs),indent=2)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['benchmark','map','all'],default='all');a=p.parse_args();OUT.mkdir(exist_ok=True,parents=True)
    if a.step in ['benchmark','all']:benchmark()
    if a.step in ['map','all']:mapping()
if __name__=='__main__':main()
