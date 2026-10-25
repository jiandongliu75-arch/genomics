"""RNA-reference sprouting annotations, explicitly separate from measured RNA.

Same-age soft-neighbor transfer is evaluated across opposite RNA libraries at
single-nucleus and disjoint metacell resolution before exploratory spatial use.
No Xenium gene counts, nuclei, identities or physical coordinates are modified.
"""
from vascular_detection_common import *
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from improve_ventricular_alignment import neighbors
import argparse

PROJ=ROOT/'results/projected_sprouting_Xenium'
SPROUT=PROGRAMS['Original_sprouting_core']
UNMEASURED=[g for g in SPROUT if g!='Angpt2']
OUTCOMES=SPROUT+CYCLING

def fitting_genes():
    with np.load(ROOT/'results/xenium_ventricular_mapping/reference_shared_panel_counts.npz') as b:g=b['genes'].astype(str)
    # Angpt2 is observed in Xenium and allowed in the primary alignment. All
    # genuinely unmeasured outcomes and cycling positive-control targets are excluded.
    return g[~np.isin(g,UNMEASURED+CYCLING+LINEAGES['Cardiomyocyte']+SEX)]

def norm_features(x):
    return np.log1p(1e4*x/np.maximum(x.sum(axis=1),1)[:,None]).astype(np.float32)

def unit(r,q,scaled=True):
    if scaled:
        r=np.clip((r-r.mean(axis=0))/np.maximum(r.std(axis=0),.25),-5,5)
        q=np.clip((q-q.mean(axis=0))/np.maximum(q.std(axis=0),.25),-5,5)
    else:
        r=r-r.mean(axis=1,keepdims=True);q=q-q.mean(axis=1,keepdims=True)
    r=r/np.maximum(np.linalg.norm(r,axis=1,keepdims=True),1e-8)
    q=q/np.maximum(np.linalg.norm(q,axis=1,keepdims=True),1e-8)
    return r.astype(np.float32),q.astype(np.float32)

def target_expression(counts,positions):
    total=np.asarray(counts.sum(axis=1)).ravel()
    x=counts[:,positions].toarray() if sparse.issparse(counts) else counts[:,positions]
    return np.log1p(1e4*x/np.maximum(total,1)[:,None]).astype(np.float32)

def metrics(y):
    return {'sprouting_6g':y[:,:6].mean(axis=1),'unmeasured_5g':y[:,[SPROUT.index(g) for g in UNMEASURED]].mean(axis=1),
            'cycling_heldout_control':y[:,6:].mean(axis=1)}

def transfer(r,q,y,k=30,scaled=True):
    rf,qf=unit(r,q,scaled);ix,w,sim=neighbors(rf,qf,k)
    mean=np.einsum('nk,nkg->ng',w,y[ix])
    second=np.einsum('nk,nkg->ng',w,y[ix]**2)
    sd=np.sqrt(np.maximum(second-mean**2,0))
    return mean,sd,ix,w,sim

def metacells(x,obs,genes,size,layer='corrected'):
    """Disjoint near-neighbor groups built only from measured fitting genes."""
    fit=fitting_genes();pos=pd.Index(genes).get_indexer(fit);assert (pos>=0).all()
    f=norm_features(x[:,pos].toarray())
    f=np.clip((f-f.mean(axis=0))/np.maximum(f.std(axis=0),.25),-5,5)
    pc=PCA(n_components=min(20,len(obs)-2),random_state=42).fit_transform(f)
    rng=np.random.default_rng(42);order=rng.permutation(len(obs));remaining=np.ones(len(obs),bool);groups=[]
    for anchor in order:
        if not remaining[anchor]:continue
        candidates=np.flatnonzero(remaining)
        distance=((pc[candidates]-pc[anchor])**2).sum(axis=1)
        take=candidates[np.argsort(distance,kind='stable')[:size]]
        remaining[take]=False;groups.append(take)
    if len(groups)>1 and len(groups[-1])<max(2,size//2):
        leftover=groups.pop()
        centers=np.array([pc[g].mean(axis=0) for g in groups])
        nearest=int(np.argmin(((centers-pc[leftover].mean(axis=0))**2).sum(axis=1)))
        groups[nearest]=np.concatenate([groups[nearest],leftover])
    assignment=np.empty(len(obs),int)
    for j,take in enumerate(groups):assignment[take]=j
    design=sparse.csr_matrix((np.ones(len(obs)),(assignment,np.arange(len(obs)))),shape=(len(groups),len(obs)))
    pooled=(design@x).astype(np.float32).tocsr()
    rows=[]
    for i,take in enumerate(groups):
        families=obs.iloc[take].EC_family.astype(str).value_counts()
        rows.append({'metacell':f'{obs["sample"].iloc[0]}__{layer}__m{size}_{i}',
            'sample':str(obs['sample'].iloc[0]),'age':str(obs.age.iloc[0]),'source_nuclei':len(take),
            'dominant_family':families.index[0],'dominant_family_fraction':families.iloc[0]/len(take),
            'source_cell_ids':';'.join(obs.index[take])})
    meta=pd.DataFrame(rows).set_index('metacell')
    return pooled,meta,design

def family_residual(v,family):
    out=np.asarray(v).copy()
    for f in set(family):
        m=family==f
        out[m]-=out[m].mean()
    return out

def build_reference(stage):
    folder=PROJ/stage;folder.mkdir(parents=True,exist_ok=True)
    data=ad.read_h5ad(OUT/'snRNA'/stage/'EC_corrected_and_raw_counts.h5ad')
    genes=data.var_names.to_numpy(dtype=str);fit=fitting_genes()
    fitpos=pd.Index(genes).get_indexer(fit);targetpos=pd.Index(genes).get_indexer(OUTCOMES)
    assert (fitpos>=0).all() and (targetpos>=0).all()
    libs=sorted(data.obs['sample'].unique());refs={}
    for layer,x in [('corrected',data.X),('raw',data.layers['counts_raw'])]:
        for lib in libs:
            m=data.obs['sample'].eq(lib).to_numpy();o=data.obs.loc[m].copy();c=x[m].tocsr()
            refs[(layer,1,lib)]={'counts':c,'meta':o,'fitting_counts':c[:,fitpos].tocsr(),
                'features':norm_features(c[:,fitpos].toarray()),'target':target_expression(c,targetpos)}
            for size in [5,20]:
                pooled,meta,design=metacells(c,o,genes,size,layer)
                refs[(layer,size,lib)]={'counts':pooled,'meta':meta,'fitting_counts':pooled[:,fitpos].tocsr(),
                    'features':norm_features(pooled[:,fitpos].toarray()),'target':target_expression(pooled,targetpos)}
                meta.to_csv(folder/f'RNA_{lib}_{layer}_metacells{size}.tsv.gz',sep='\t')
    return data,genes,fit,refs

def benchmark(stage,refs,fit,folder):
    libs=sorted({key[2] for key in refs});rows=[];predictions=[]
    for size in [1,5,20]:
        for querylib in libs:
            trainlib=next(x for x in libs if x!=querylib)
            r,q=refs[('corrected',size,trainlib)],refs[('corrected',size,querylib)]
            family=(q['meta'].EC_family if size==1 else q['meta'].dominant_family).to_numpy(dtype=str)
            truth=metrics(q['target'])
            for method in ['assay_scaled_soft','Pearson_soft']:
                for k in [10,30,60]:
                    y,sd,ix,w,similarity=transfer(r['features'],q['features'],r['target'],k,method=='assay_scaled_soft')
                    scores=metrics(y)
                    for outcome,actual in truth.items():
                        pred=scores[outcome]
                        rho=corr(actual,pred);within=corr(family_residual(actual,family),family_residual(pred,family))
                        r2=float(r2_score(actual,pred)) if np.std(actual)>1e-9 else np.nan
                        rows.append({'stage':stage,'query_library':querylib,'RNA_resolution_nuclei':size,'method':method,'neighbors':k,'effective_neighbors':min(k,len(r['features'])),
                            'outcome':outcome,'query_units':len(q['features']),'reference_units':len(r['features']),
                            'Spearman':rho,'within_family_Spearman':within,'R2_vs_query_mean':r2,
                            'source_gene_targets':';'.join(SPROUT if outcome=='sprouting_6g' else UNMEASURED if outcome=='unmeasured_5g' else CYCLING),
                            'target_leakage_check':'Unmeasured5 and cycling targets absent; measured Angpt2 allowed for primary6g score',
                            'technical_support_gate':bool(np.isfinite(rho) and rho>=.2 and np.isfinite(within) and within>=.15 and np.isfinite(r2) and r2>0)})
                    if method=='assay_scaled_soft' and k==30:
                        for i,unit_id in enumerate(q['meta'].index):
                            predictions.append({'stage':stage,'query_library':querylib,'RNA_resolution_nuclei':size,'unit_id':unit_id,
                                'family':family[i],**{'observed_'+name:float(v[i]) for name,v in truth.items()},
                                **{'projected_'+name:float(v[i]) for name,v in scores.items()},'nearest_reference_similarity':float(similarity[i])})
            # Explicit measured-Angpt2 withholding sensitivity at the primary setting.
            use=fit!='Angpt2'
            y,*_=transfer(norm_features(r['fitting_counts'][:,use].toarray()),
                norm_features(q['fitting_counts'][:,use].toarray()),r['target'],30,True)
            for outcome,actual in truth.items():
                pred=metrics(y)[outcome];rho=corr(actual,pred);within=corr(family_residual(actual,family),family_residual(pred,family))
                r2=float(r2_score(actual,pred)) if np.std(actual)>1e-9 else np.nan
                rows.append({'stage':stage,'query_library':querylib,'RNA_resolution_nuclei':size,'method':'assay_scaled_Angpt2_withheld','neighbors':30,'effective_neighbors':min(30,len(r['features'])),
                    'outcome':outcome,'query_units':len(q['features']),'reference_units':len(r['features']),'Spearman':rho,
                    'within_family_Spearman':within,'R2_vs_query_mean':r2,'source_gene_targets':';'.join(SPROUT if outcome=='sprouting_6g' else UNMEASURED if outcome=='unmeasured_5g' else CYCLING),
                    'target_leakage_check':'All original sprouting and cycling targets absent',
                    'technical_support_gate':bool(np.isfinite(rho) and rho>=.2 and np.isfinite(within) and within>=.15 and np.isfinite(r2) and r2>0)})
    result=pd.DataFrame(rows);save(result,folder/'opposite_library_projection_benchmark.tsv')
    save(pd.DataFrame(predictions),folder/'opposite_library_primary_predictions.tsv.gz')
    return result

def project_stage(stage,data,genes,fit,refs,bench,folder):
    source=ad.read_h5ad(OUT/'Xenium'/stage/'EC_candidate_count_views.h5ad')
    pos=source.var_names.get_indexer(fit);assert (pos>=0).all()
    eligible=(source.obs.native_EC_broad&source.obs.nuclear_broad_EC_clean&~source.obs.multiple_projected_nuclear_components_review&~source.obs.geometry_unresolved_review).to_numpy()
    obs=source.obs.loc[eligible].copy();ids=obs.index
    values={'native':source.layers['native'][eligible].tocsr(),'nuclear':source.layers['nuclear_only'][eligible].tocsr(),
            'boundary_core':source.layers['high_confidence'][eligible].tocsr()}
    libs=sorted({key[2] for key in refs});frames={};gene_predictions={};neighbor_files={}
    for view,counts in values.items():
        q=norm_features(counts[:,pos].toarray())
        f=obs[['x_centroid','y_centroid','atlas_identity','nuclear_supported_family']].copy()
        for layer in ['corrected','raw']:
            for resolution in [1,5,20]:
                reference=np.vstack([refs[(layer,resolution,l)]['features'] for l in libs])
                targets=np.vstack([refs[(layer,resolution,l)]['target'] for l in libs])
                y,sd,ix,w,sim=transfer(reference,q,targets,30,True)
                for name,v in metrics(y).items():f[f'{layer}__RNA_resolution{resolution}__inferred_{name}']=v
                if layer=='corrected' and resolution==1:
                    gene_predictions[view]=y[:,:6]
                    f['reference_neighbor_effective_n']=1/np.maximum((w*w).sum(axis=1),1e-8)
                    f['nearest_reference_similarity']=sim
                    f['neighbor_sprouting_score_SD']=np.sqrt(np.maximum(np.sum(w*metrics(targets)['sprouting_6g'][ix]**2,axis=1)-metrics(y)['sprouting_6g']**2,0))
                    refs_ids=np.concatenate([refs[(layer,resolution,l)]['meta'].index.to_numpy(dtype=str) for l in libs])
                    np.savez_compressed(folder/f'{view}_RNA_reference_neighbors.npz',query_cell_ids=ids.to_numpy(dtype=str),reference_cell_ids=refs_ids,
                        neighbor_indices=ix,neighbor_weights=w)
                    for lib in libs:
                        r=refs[(layer,resolution,lib)]
                        yl,*_=transfer(r['features'],q,r['target'],30,True)
                        for name,v in metrics(yl).items():f[f'RNA_library_{lib}__inferred_{name}']=v
                    # Single-cell reference neighborhood-size sensitivity.
                    for k in [10,60]:
                        ys,*_=transfer(reference,q,targets,k,True)
                        for name,v in metrics(ys).items():f[f'neighbors{k}__inferred_{name}']=v
                    use=fit!='Angpt2'
                    refcounts=sparse.vstack([refs[(layer,resolution,l)]['fitting_counts'] for l in libs],format='csr')
                    ya,*_=transfer(norm_features(refcounts[:,use].toarray()),
                        norm_features(counts[:,pos][:,use].toarray()),targets,30,True)
                    for name,v in metrics(ya).items():f[f'Angpt2_withheld__inferred_{name}']=v
        f.to_csv(folder/f'{view}_inferred_sprouting_annotations.tsv.gz',sep='\t',na_rep='NA');frames[view]=f
    # Spatially pooled query projection: true observed counts of5 ECs, not
    # smoothing a categorical state label or pretending each pool is a new cell.
    xy=obs[['x_centroid','y_centroid']].to_numpy();pools_by_size={}
    for poolsize,radius in [(5,75),(20,150)]:
        distance,index=cKDTree(xy).query(xy,k=poolsize,distance_upper_bound=radius)
        keep=np.isfinite(distance).all(axis=1);ix=index[keep]
        design=sparse.csr_matrix((np.ones(ix.size),(np.repeat(np.arange(len(ix)),poolsize),ix.ravel())),shape=(len(ix),len(obs)))
        pools=obs.loc[keep,['x_centroid','y_centroid']].copy();pools['source_EC_ids']=[';'.join(ids[x]) for x in ix]
        pools['farthest_source_distance_um']=distance[keep,-1];pools['source_EC_profiles']=poolsize
        r=np.vstack([refs[('corrected',poolsize,l)]['features'] for l in libs]);yref=np.vstack([refs[('corrected',poolsize,l)]['target'] for l in libs])
        for view,counts in values.items():
            q=norm_features((design@counts[:,pos]).toarray())
            y,*_=transfer(r,q,yref,30,True)
            for name,v in metrics(y).items():pools[f'{view}__inferred_neighborhood{poolsize}_{name}']=v
        pools.to_csv(folder/f'neighborhood{poolsize}_inferred_sprouting.tsv.gz',sep='\t',na_rep='NA')
        pools_by_size[poolsize]=pools
    primary=frames['nuclear'];annotation=source.obs.copy()
    annotation['stage_cell_id']=[stage+'::'+c for c in annotation.index]
    annotation['sprouting_inference_eligible']=eligible
    annotation['sprouting_annotation_status']='Not projected: source EC/nuclear/geometry QC criteria unmet'
    annotation.loc[ids,'sprouting_annotation_status']='Exploratory RNA-reference projection; unmeasured genes inferred'
    annotation['inferred_sprouting_score_6g']=primary['corrected__RNA_resolution1__inferred_sprouting_6g'].reindex(annotation.index)
    annotation['inferred_unmeasured_sprouting_score_5g']=primary['corrected__RNA_resolution1__inferred_unmeasured_5g'].reindex(annotation.index)
    annotation['inferred_cycling_control_score']=primary['corrected__RNA_resolution1__inferred_cycling_heldout_control'].reindex(annotation.index)
    for poolsize in [5,20]:
        pools=pools_by_size[poolsize]
        annotation[f'inferred_neighborhood{poolsize}_sprouting_score_6g']=pools[f'nuclear__inferred_neighborhood{poolsize}_sprouting_6g'].reindex(annotation.index)
        annotation[f'inferred_neighborhood{poolsize}_unmeasured_score_5g']=pools[f'nuclear__inferred_neighborhood{poolsize}_unmeasured_5g'].reindex(annotation.index)
    annotation['measured_Angpt2_native_molecules']=np.asarray(source.layers['native'][:,source.var_names.get_loc('Angpt2')].toarray()).ravel()
    annotation['inferred_sprouting_boundary_difference']=np.abs(primary.corrected__RNA_resolution1__inferred_sprouting_6g-frames['boundary_core'].corrected__RNA_resolution1__inferred_sprouting_6g).reindex(annotation.index)
    library_columns=[f'RNA_library_{l}__inferred_sprouting_6g' for l in libs]
    annotation['inferred_sprouting_RNA_library_difference']=np.abs(primary[library_columns[0]]-primary[library_columns[1]]).reindex(annotation.index)
    reference_predictions=pd.read_csv(folder/'opposite_library_primary_predictions.tsv.gz',sep='\t')
    similarity_floor=float(reference_predictions.loc[reference_predictions.RNA_resolution_nuclei.eq(1),'nearest_reference_similarity'].quantile(.05))
    annotation['inferred_sprouting_nearest_reference_similarity']=primary.nearest_reference_similarity.reindex(annotation.index)
    annotation['alignment_reference_similarity_status']=np.select([~annotation.sprouting_inference_eligible,
        annotation.inferred_sprouting_nearest_reference_similarity.lt(similarity_floor)],
        ['Not projected','Below opposite-library RNA reference 5th percentile'],default='Within opposite-library RNA reference similarity range')
    gate=[]
    for size in [1,5,20]:
        d=bench[(bench.method=='assay_scaled_soft')&(bench.RNA_resolution_nuclei==size)&bench.outcome.eq('unmeasured_5g')&bench.neighbors.eq(30)]
        passed=len(d)==2 and d.technical_support_gate.all()
        gate.append({'stage':stage,'RNA_resolution_nuclei':size,'unmeasured_5g_both_library_gate_pass':bool(passed)})
    annotation['unmeasured_component_individual_technical_support']=gate[0]['unmeasured_5g_both_library_gate_pass']
    annotation['unmeasured_component_metacell5_technical_support']=gate[1]['unmeasured_5g_both_library_gate_pass']
    annotation['unmeasured_component_metacell20_technical_support']=gate[2]['unmeasured_5g_both_library_gate_pass']
    annotation.to_csv(folder/'EC_annotations_with_inferred_sprouting.tsv.gz',sep='\t',na_rep='NA')
    # Measured genes remain the only403 var features. Reference predictions are
    # separate normalized-expression annotations, never fabricated count layers.
    augmented=ad.AnnData(values['native'].copy(),obs=annotation.loc[ids].copy(),var=source.var.copy())
    augmented.layers['nuclear_only']=values['nuclear'];augmented.layers['boundary_core']=values['boundary_core']
    augmented.obsm['spatial']=obs[['x_centroid','y_centroid']].to_numpy()
    augmented.obsm['inferred_sprouting_genes_RNA_log1pCP10k']=pd.DataFrame(gene_predictions['nuclear'],index=ids,columns=SPROUT)
    augmented.uns['sprouting_projection']={'status':'Exploratory reference-derived expected RNA, not measured Xenium counts',
        'measured_var_features':403,'inferred_gene_units':'Weighted mean of same-age corrected snRNA log1pCP10k,30 soft cosine neighbors; not transcript counts',
        'target_genes':SPROUT,'unmeasured_targets':UNMEASURED,'measured_Angpt2_preserved':True,
        'score_6g':'Mean inferred expression over six genes','score_5g':'Mean over the five unmeasured genes; reported separately to expose limited missing-gene recovery',
        'feature_genes':fit.tolist(),'alignment':'Separate assay gene scaling,SDfloor0.25,clip5,L2 cosine,temp0.1; matched age only; no spatial coordinates or program labels in fit',
        'neighborhood':'Sum observed RNA of5 nuclear-clean ECs within75um or20 within150um; project to matched-resolution disjoint RNA metacells. Pools overlap and are not independent cells.'}
    augmented.uns['spatial_coordinate_units']='micrometers'
    augmented.write_h5ad(folder/'Xenium_ECs_with_inferred_sprouting.h5ad',compression='gzip')
    atlas=pd.read_csv(ROOT/'results/whole_heart_spatial_atlas'/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    whole=atlas[['x_centroid','y_centroid','atlas_identity']].copy();whole['stage']=stage
    whole['stage_cell_id']=[stage+'::'+c for c in whole.index]
    for c in ['sprouting_inference_eligible','sprouting_annotation_status','inferred_sprouting_score_6g','inferred_unmeasured_sprouting_score_5g',
              'inferred_neighborhood5_sprouting_score_6g','inferred_neighborhood5_unmeasured_score_5g','inferred_sprouting_RNA_library_difference','inferred_sprouting_boundary_difference']:
        whole[c]=annotation[c].reindex(whole.index)
    for c in ['inferred_neighborhood20_sprouting_score_6g','inferred_neighborhood20_unmeasured_score_5g']:
        whole[c]=annotation[c].reindex(whole.index)
    whole['sprouting_inference_eligible']=whole.sprouting_inference_eligible.fillna(False)
    whole['sprouting_annotation_status']=whole.sprouting_annotation_status.fillna('Not projected: non-EC or unavailable EC evidence')
    whole.to_csv(folder/'all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz',sep='\t',na_rep='NA')
    save(pd.DataFrame(gate),folder/'projection_support_gates.tsv')
    write_json(folder/'complete.json',{'stage':stage,'source_EC_candidates':len(annotation),'projected_ECs':len(ids),
        'projected_spatial_neighborhoods5':len(pools_by_size[5]),'projected_spatial_neighborhoods20':len(pools_by_size[20]),
        'whole_heart_profiles':len(whole),'observed_gene_counts_unchanged':True,
        'reference_nuclei':len(data),'inferred_gene_targets':SPROUT,'unmeasured_gene_targets':UNMEASURED,
        'opposite_library_reference_similarity_5pct':similarity_floor,
        'RNA_libraries':libs,'gate_results':gate,'label':'Exploratory inferred sprouting propensity, not measured angiogenesis',
        'primary_reference_resolution':1,'primary_neighbors':30,'all_spatial_outputs_are_exploratory':True})
    print('SPROUTING_PROJECTION_COMPLETE',stage,len(ids),gate,flush=True)

def main(stage):
    folder=PROJ/stage;folder.mkdir(parents=True,exist_ok=True)
    data,genes,fit,refs=build_reference(stage)
    rows=[]
    for lib in sorted(data.obs['sample'].unique()):
        m=data.obs['sample'].eq(lib).to_numpy()
        for g in SPROUT:
            j=data.var_names.get_loc(g)
            rows.append({'stage':stage,'sample':lib,'gene':g,'nuclei':int(m.sum()),'raw_molecules':float(data.layers['counts_raw'][m,j].sum()),
                'corrected_expected_mass':float(data.X[m,j].sum()),'raw_observed_and_corrected_ge0p5_nuclei':int(((data.layers['counts_raw'][m,j]>=1).multiply(data.X[m,j]>=.5)).sum())})
    save(pd.DataFrame(rows),folder/'RNA_reference_sprouting_gene_support.tsv')
    save(pd.DataFrame({'gene':fit,'measured_in_Xenium':True,'is_measured_sprouting_target':fit=='Angpt2'}),folder/'alignment_feature_genes.tsv')
    bench=benchmark(stage,refs,fit,folder)
    project_stage(stage,data,genes,fit,refs,bench,folder)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=STAGES);a=p.parse_args();main(a.stage)
