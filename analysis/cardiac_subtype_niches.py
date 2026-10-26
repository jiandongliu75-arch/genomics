"""Known-subtype neighborhoods and within-stage, within-CM-state maturation.

No de novo niche clustering is used. Unknown subtypes are retained in lineage
denominators. Frozen native coordinates, labels and count views are unchanged.
"""
from pathlib import Path
import os,json,argparse,hashlib,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_subtype_niches'
UP=ROOT/'results/cardiac_spatial_pipeline';RECON=ROOT/'results/cardiac_spatial_reconstruction'
STAGES=['P0','P7','P14','P21'];CM=['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM']
FB=['Postn-associated','Gfpt2-associated','Col8a1-associated','Lsamp-associated','Kirrel3-associated']
EC=['Capillary-like','Endocardial-like','Arterial-like','Venous-like','Lymphatic-like']
SEED=271026;DRAWS=512;SOURCES=set()

def source(path):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;SOURCES.add(str(path));return path
def read(path,**kw):return pd.read_csv(source(path),sep='\t',**kw)
def save(data,name):
    OUT.mkdir(exist_ok=True,parents=True);data.to_csv(OUT/name,sep='\t',index=False,float_format='%.9g',na_rep='NA')
def stage_data(stage):
    f=read(UP/stage/'atlas_annotations.tsv.gz',index_col=0)
    c=read(RECON/stage/'CM_outcomes_and_anatomy.tsv.gz',index_col=0)
    take=f.ventricular_CM_strict&f.author_atlas_cell_state.isin(CM)
    d=f.loc[take].join(c[['core_panel_depth','native_panel_depth','reconstructed_Figure4e_CM','morphology_region','distance_to_large_void_um','distance_to_outer_um']+[k for k in c if k.startswith(('core_Z__','native_Z__'))]])
    d=d[d.reconstructed_Figure4e_CM.fillna(False)&d.core_panel_depth.ge(20)&d.CM_outcome_excluded_depth.ge(20)].copy()
    d['CM_state']=d.author_atlas_cell_state
    d['stage_maturation_percentile']=d.core_Z__maturation21.rank(pct=True)
    d['state_maturation_percentile']=d.groupby('CM_state').core_Z__maturation21.rank(pct=True)
    d['native_state_maturation_percentile']=d.groupby('CM_state').native_Z__maturation21.rank(pct=True)
    d['maturation_group']=np.where(d.state_maturation_percentile.le(.25),'Low',np.where(d.state_maturation_percentile.ge(.75),'High','Middle'))
    d['block200']=(d.x_centroid//200).astype(int).astype(str)+'_'+(d.y_centroid//200).astype(int).astype(str)
    d['block400']=(d.x_centroid//400).astype(int).astype(str)+'_'+(d.y_centroid//400).astype(int).astype(str)
    return f,d

def prepare():
    coverage=[];means=[];cohorts=[];defs=[]
    for name in FB:defs.append({'exposure':'FB__'+name,'source':'Cell+nuclear-supported fibroblast family','denominator':'All eligible broad fibroblast neighbors, including unknown family'})
    for name in EC:defs.append({'exposure':'EC__'+name,'source':'Nuclear-supported EC type, with broad EC identity/QC','denominator':'All eligible broad EC neighbors, including unknown type'})
    for stage in STAGES:
        print('SUBTYPE_PREPARE',stage,flush=True);f,d=stage_data(stage);allids=f.index;positions=allids.get_indexer(d.index)
        eligible=f.program_source_eligible.fillna(False)
        fb=f.atlas_identity.eq('Fibroblast')&eligible;ec=f.atlas_identity.eq('Endothelial')&eligible
        masks={**{'FB__'+t:fb&f.fine_identity.eq(t) for t in FB},**{'EC__'+t:ec&f.EC_nuclear_type.eq(t) for t in EC}}
        strict={**{'FB__'+t:fb&f.fine_identity.eq(t) for t in FB},**{'EC__'+t:ec&f.fine_identity.eq(t) for t in EC}}
        aliases={'FB__Postn-associated':'Postn.Fib','FB__Gfpt2-associated':'Gfpt2.Fib','FB__Col8a1-associated':'Col8a1.Fib','FB__Lsamp-associated':'Lsamp.Fib',
            'EC__Capillary-like':'cap.EC','EC__Endocardial-like':'endocardial.EC','EC__Arterial-like':'arteriole.EC','EC__Venous-like':'venous.EC','EC__Lymphatic-like':'lym.EC'}
        author={key:(fb if key.startswith('FB') else ec)&f.author_atlas_cell_state.eq(aliases.get(key,'__no_public_equivalent__')) for key in masks}
        for key,mask in masks.items():coverage.append({'stage':stage,'subtype':key,'primary_source_profiles':int(mask.sum()),'strict_cell_nuclear_profiles':int(strict[key].sum()),
            'published_comparator_profiles':int(author[key].sum()),
            'all_eligible_lineage_profiles':int((fb if key.startswith('FB') else ec).sum()),'supported_for_primary_model':int(mask.sum())>=30,
            'zero_or_sparse_is_biological_absence':False})
        for geometry in ['proximity50','proximity30','contact']:
            graph=sparse.load_npz(source(UP/stage/(geometry+'_graph.npz'))).tocsr()[positions]
            nfb=np.asarray(graph@fb.to_numpy(float)).ravel();nec=np.asarray(graph@ec.to_numpy(float)).ravel()
            resolved=f.atlas_resolved.to_numpy(float);nr=np.asarray(graph@resolved).ravel()
            out=d[['stage','CM_state','x_centroid','y_centroid','cell_area','nucleus_area','CM_outcome_excluded_depth','core_panel_depth','native_panel_depth',
                'core_Z__maturation21','native_Z__maturation21','core_Z__metabolic','core_Z__sarcomere','core_Z__excitation_contraction','core_Z__withdrawal',
                'stage_maturation_percentile','state_maturation_percentile','native_state_maturation_percentile','maturation_group','block200','block400',
                'arterial_distance_um','endocardial_proxy_distance_um','morphology_region','distance_to_large_void_um','distance_to_outer_um']].copy()
            out['cell_id']=out.index;out['neighbors_resolved']=nr;out['eligible_FB_neighbors']=nfb;out['eligible_EC_neighbors']=nec
            for lineage in ['Fibroblast','Endothelial','Pericyte','Smooth muscle','Immune','Epicardial']:
                out['broad__'+lineage]=np.asarray(graph@f.atlas_identity.eq(lineage).to_numpy(float)).ravel()/np.maximum(nr,1)
            for key,mask in masks.items():
                denom=nfb if key.startswith('FB') else nec
                out[key]=np.asarray(graph@mask.to_numpy(float)).ravel()/np.maximum(denom,1)
                out['strict__'+key]=np.asarray(graph@strict[key].to_numpy(float)).ravel()/np.maximum(denom,1)
                out['author__'+key]=np.asarray(graph@author[key].to_numpy(float)).ravel()/np.maximum(denom,1)
            for label,mask,denom in [('FB',fb,nfb),('EC',ec,nec)]:
                out[label+'_mean_log_capture']=np.asarray(graph@(mask.to_numpy(float)*f.source_log_capture.fillna(0).to_numpy(float))).ravel()/np.maximum(denom,1)
            out['FB_typed_fraction']=out[['FB__'+t for t in FB]].sum(axis=1);out['EC_typed_fraction']=out[['EC__'+t for t in EC]].sum(axis=1)
            out['geometry']=geometry;out['analysis_eligible']=(nr>=5)&(nfb>=1)&(nec>=1)
            save(out,stage+'_'+geometry+'_CM_subtype_neighborhoods.tsv.gz')
            if geometry!='proximity50':continue
            for state,g in out.groupby('CM_state',observed=True):
                cohorts.append({'stage':stage,'CM_state':state,'spatial_CM_profiles':len(g),'analysis_eligible_CM_profiles':int(g.analysis_eligible.sum()),
                    'source_typing_note':'Primary EC is nuclear supported; FB is cell+nuclear supported; unknowns retained in lineage denominator'})
                for group,h in g[g.analysis_eligible].groupby('maturation_group'):
                    for ex in masks:
                        means.append({'stage':stage,'CM_state':state,'maturation_group':group,'source_subtype':ex,'CM_profiles':len(h),
                            'mean_within_lineage_fraction':h[ex].mean(),'fraction_CMs_with_subtype_neighbor':h[ex].gt(0).mean(),
                            'mean_FB_typed_fraction':h.FB_typed_fraction.mean(),'mean_EC_typed_fraction':h.EC_typed_fraction.mean(),
                            'scope':'Descriptive source fractions, not anatomical contact, unbiased abundance or a causal effect'})
    save(pd.DataFrame(coverage),'source_subtype_coverage.tsv');save(pd.DataFrame(cohorts),'CM_state_cohort_coverage.tsv')
    save(pd.DataFrame(means),'CM_state_maturation_subtype_neighborhood_means.tsv');save(pd.DataFrame(defs),'exposure_definitions.tsv')

def center(a,groups):
    if a.ndim==1:a=a[:,None]
    _,ix=np.unique(groups,return_inverse=True);counts=np.bincount(ix);means=np.zeros((len(counts),a.shape[1]))
    np.add.at(means,ix,a);means/=counts[:,None]
    return a-means[ix]

def model_arrays(d,exposures,geometry='proximity50',view='primary',adjustment='within_region'):
    ykeys=['core_Z__maturation21','native_Z__maturation21','core_Z__metabolic','core_Z__sarcomere','core_Z__excitation_contraction','core_Z__withdrawal']
    cov={'log_CM_capture':np.log1p(d.CM_outcome_excluded_depth),'log_CM_area':np.log1p(d.cell_area),
        'nucleus_missing':d.nucleus_area.isna().astype(float),'log_neighbors':np.log1p(d.neighbors_resolved),
        'log_FB_count':np.log1p(d.eligible_FB_neighbors),'log_EC_count':np.log1p(d.eligible_EC_neighbors),
        'log_arterial_distance':np.log1p(d.arterial_distance_um),'log_endocardial_distance':np.log1p(d.endocardial_proxy_distance_um),
        'log_void_distance':np.log1p(d.distance_to_large_void_um),'log_outer_distance':np.log1p(d.distance_to_outer_um),
        'FB_capture':d.FB_mean_log_capture,'EC_capture':d.EC_mean_log_capture,
        **{k:d[k] for k in d if k.startswith('broad__')}}
    if adjustment=='technical_composition':
        for k in ['log_arterial_distance','log_endocardial_distance','log_void_distance','log_outer_distance']:cov.pop(k)
    X=pd.DataFrame(cov,index=d.index)
    prefix={'primary':'','strict_cell_nuclear':'strict__','published_reference_comparator':'author__'}[view]
    for ex in exposures:X[ex]=d[prefix+ex]
    y=d[ykeys].to_numpy(float);x=X.to_numpy(float)
    # Common stage-level SDs preserve coefficient comparability across CM states.
    x=(x-np.nanmean(x,0))/np.maximum(np.nanstd(x,0),1e-6)
    y=(y-np.nanmean(y,0))/np.maximum(np.nanstd(y,0),1e-6)
    return X.columns.tolist(),x,y,ykeys

def fit_states(d,exposures,geometry,view='primary',adjustment='within_region'):
    names,x,y,ykeys=model_arrays(d,exposures,geometry,view,adjustment)
    finite=np.isfinite(x).all(1)&np.isfinite(y).all(1)
    allblocks=np.unique(d.loc[finite,'block400']);blockmap={b:i for i,b in enumerate(allblocks)}
    rng=np.random.default_rng(SEED);weights=np.zeros((DRAWS,len(allblocks)))
    for b in range(DRAWS):weights[b]=np.bincount(rng.integers(len(allblocks),size=len(allblocks)),minlength=len(allblocks))
    rows=[];state_results={}
    for state in CM:
        take=finite&d.CM_state.eq(state).to_numpy()
        if take.sum()<200:continue
        z=d.loc[take]
        groups=z.block200.to_numpy() if adjustment=='within_region' else np.full(len(z),'state_intercept')
        xx=center(x[take],groups);yy=center(y[take],groups)
        keep=xx.std(0)>.015
        if keep.sum()<3:continue
        xx=xx[:,keep];kept=np.array(names)[keep];nb=len(allblocks);p=xx.shape[1];q=yy.shape[1]
        xtx=np.zeros((nb,p,p));xty=np.zeros((nb,p,q))
        for key,g in z.groupby('block400'):
            ii=z.index.get_indexer(g.index);j=blockmap[key];xtx[j]=xx[ii].T@xx[ii];xty[j]=xx[ii].T@yy[ii]
        beta=np.linalg.pinv(xtx.sum(0),rcond=1e-9)@xty.sum(0)
        bxx=np.einsum('bn,nij->bij',weights,xtx);bxy=np.einsum('bn,nij->bij',weights,xty)
        boot=np.stack([np.linalg.pinv(A,rcond=1e-8)@B for A,B in zip(bxx,bxy)])
        cond=float(np.linalg.cond(xtx.sum(0)))
        for j,name in enumerate(kept):
            if name not in exposures:continue
            for k,outcome in enumerate(ykeys):
                low,high=np.quantile(boot[:,j,k],[.025,.975]);se=np.std(boot[:,j,k],ddof=1)
                rows.append({'stage':d.stage.iloc[0],'CM_state':state,'geometry':geometry,'source_view':view,'adjustment':adjustment,
                    'source_subtype':name,'outcome':outcome,'standardized_coefficient':beta[j,k],'block_CI_low':low,'block_CI_high':high,
                    'within_section_approx_p':2*norm.sf(abs(beta[j,k])/max(se,1e-9)),'CM_profiles':len(z),'blocks400':z.block400.nunique(),
                    'design_condition_number':cond,'source_exposure_partial_SD':xx[:,j].std(),
                    'scope':'Known subtype predictors within CM identity, with explicit technical/regional adjustment levels; one section per age'})
        state_results[state]=(kept,beta,boot)
    differences=[]
    if all(s in state_results for s in ['v.CM','Ankrd1.CM']):
        a,A,Ab=state_results['Ankrd1.CM'];v,V,Vb=state_results['v.CM']
        for name in exposures:
            if name not in a or name not in v:continue
            ai=np.flatnonzero(a==name)[0];vi=np.flatnonzero(v==name)[0]
            for k,outcome in enumerate(ykeys):
                delta=A[ai,k]-V[vi,k];bs=Ab[:,ai,k]-Vb[:,vi,k];lo,hi=np.quantile(bs,[.025,.975])
                differences.append({'stage':d.stage.iloc[0],'geometry':geometry,'source_view':view,'adjustment':adjustment,
                    'source_subtype':name,'outcome':outcome,'Ankrd1_minus_working_slope':delta,'block_CI_low':lo,'block_CI_high':hi,
                    'scope':'Shared spatial-block bootstrap, common stage SD units; not independent-animal interaction validation'})
    return rows,differences

def models():
    coverage=read(OUT/'source_subtype_coverage.tsv');rows=[];differences=[]
    for stage in STAGES:
        for geometry in ['proximity50','proximity30','contact']:
            print('SUBTYPE_MODELS',stage,geometry,flush=True)
            d=read(OUT/(stage+'_'+geometry+'_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id');d=d[d.analysis_eligible].copy()
            for view,countcol,prefix in [('primary','primary_source_profiles',''),('strict_cell_nuclear','strict_cell_nuclear_profiles','strict__'),('published_reference_comparator','published_comparator_profiles','author__')]:
                if geometry!='proximity50' and view!='primary':continue
                cov=coverage[coverage.stage.eq(stage)&coverage[countcol].ge(30)]
                ex=[s for s in cov.subtype if d[prefix+s].gt(0).sum()>=50]
                if not ex:continue
                for adjustment in (['technical_composition','anatomy_proxies','within_region'] if geometry=='proximity50' else ['within_region']):
                    r,delta=fit_states(d,ex,geometry,view,adjustment);rows.extend(r);differences.extend(delta)
    out=pd.DataFrame(rows)
    # BH is a within-section exploratory annotation; it is not animal-level FDR.
    out['within_section_BH_q']=np.nan
    for _,g in out.groupby(['stage','geometry','source_view','adjustment']):
        p=g.within_section_approx_p.to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
        out.loc[g.index[order],'within_section_BH_q']=np.minimum(q,1)
    save(out,'subtype_conditioned_maturation_coefficients.tsv')
    save(pd.DataFrame(differences),'Ankrd1_versus_working_CM_slope_differences.tsv')

def matched():
    effects=[];records=[];rng=np.random.default_rng(SEED)
    for stage in STAGES:
        print('SUBTYPE_MATCHED',stage,flush=True)
        d=read(OUT/(stage+'_proximity50_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id');d=d[d.analysis_eligible].copy()
        x=np.column_stack([np.log1p(d.CM_outcome_excluded_depth),np.log1p(d.cell_area),np.log1p(d.eligible_FB_neighbors),np.log1p(d.eligible_EC_neighbors),d.FB_mean_log_capture,d.EC_mean_log_capture])
        x=(x-x.mean(0))/np.maximum(x.std(0),1e-6)
        for state in ['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM']:
            z=d[d.CM_state.eq(state)];pairs=[]
            for block,g in z.groupby('block200'):
                hi=g.index[g.maturation_group.eq('High')];lo=list(g.index[g.maturation_group.eq('Low')])
                if len(hi)==0 or len(lo)==0:continue
                for cid in rng.permutation(hi):
                    if not lo:break
                    i=d.index.get_loc(cid);jj=d.index.get_indexer(lo);diff=np.abs(x[jj]-x[i]);distance=(diff**2).sum(1)
                    xy=d.loc[lo,['x_centroid','y_centroid']].to_numpy()-d.loc[cid,['x_centroid','y_centroid']].to_numpy()
                    sep=np.sqrt((xy**2).sum(1));endo=np.abs(d.loc[lo,'endocardial_proxy_distance_um'].to_numpy()-d.loc[cid,'endocardial_proxy_distance_um'])
                    fail=(sep>100)|(diff[:,:4].max(1)>.5)|(diff[:,4:].max(1)>.75)|(endo>50)
                    distance[fail]=np.inf;j=int(np.argmin(distance))
                    if not np.isfinite(distance[j]):continue
                    control=lo.pop(j);pairs.append((cid,control));records.append({'stage':stage,'CM_state':state,'high_maturation_CM':cid,'low_maturation_CM':control,
                        'block200':block,'block400':d.loc[cid,'block400'],'distance_um':sep[j],'covariate_distance':distance[j]})
            if len(pairs)<20:continue
            high,low=zip(*pairs);blocks=d.loc[list(high),'block400'];codes=pd.factorize(blocks)[0];nb=codes.max()+1
            for ex in ['FB__'+f for f in FB]+['EC__'+t for t in EC]:
                delta=d.loc[list(high),ex].to_numpy()-d.loc[list(low),ex].to_numpy();sums=np.bincount(codes,weights=delta);counts=np.bincount(codes);boot=[]
                for _ in range(DRAWS):
                    ii=rng.integers(nb,size=nb);boot.append(sums[ii].sum()/counts[ii].sum())
                lo95,hi95=np.quantile(boot,[.025,.975]);effects.append({'stage':stage,'CM_state':state,'source_subtype':ex,'matched_pairs':len(pairs),
                    'high_minus_low_fraction':delta.mean(),'block_CI_low':lo95,'block_CI_high':hi95,
                    'scope':'Same CM state and local200µm region; capture/area/lineage-count/source-capture/endocardial calipers; within-section exploratory contrast'})
    save(pd.DataFrame(records),'within_CM_state_maturation_matched_pairs.tsv.gz');save(pd.DataFrame(effects),'matched_maturation_subtype_neighborhood_differences.tsv')

def paired_geometry():
    coverage=read(OUT/'source_subtype_coverage.tsv');rows=[]
    for stage in ['P7','P14']:
        a=read(OUT/(stage+'_proximity50_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id')
        for comparison in ['proximity30','contact']:
            b=read(OUT/(stage+'_'+comparison+'_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id')
            ids=a.index[a.analysis_eligible&b.analysis_eligible]
            ex=coverage.loc[coverage.stage.eq(stage)&coverage.primary_source_profiles.ge(30),'subtype'].tolist()
            for geometry,z in [('proximity50',a.loc[ids]),(comparison,b.loc[ids])]:
                print('SUBTYPE_PAIRED_GEOMETRY',stage,comparison,geometry,len(z),flush=True)
                result,_=fit_states(z,ex,geometry,'primary','within_region')
                rows.extend([{'paired_comparison':comparison,'same_CM_cohort':True,**r} for r in result])
    save(pd.DataFrame(rows),'same_CM_geometry_sensitivity.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['prepare','models','matched','paired','all'],default='all');a=p.parse_args();OUT.mkdir(exist_ok=True,parents=True);start=time.time()
    if a.step in ['prepare','all']:prepare()
    if a.step in ['models','all']:models()
    if a.step in ['matched','all']:matched()
    if a.step in ['paired','all']:paired_geometry()
    path=OUT/'source_files.json'
    if path.exists():SOURCES.update(json.loads(path.read_text()))
    path.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')
    (OUT/('execution_'+a.step+'.json')).write_text(json.dumps({'completed':True,'seconds':time.time()-start,'step':a.step,'seed':SEED,'new_biological_replicates':0},indent=2)+'\n')

if __name__=='__main__':main()
