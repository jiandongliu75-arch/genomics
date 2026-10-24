"""Re-evaluate fine FB programs and measured CM outcomes on audited boundaries."""
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from fibroblast_state_common import *
from map_fibroblast_states_Xenium import FibroblastProjection,annotate,direct_programs,cm_scores,SOURCE,ATLAS
from compare_fibroblast_matrix_and_spatial import baseline_features


def partial_r(y,x,base):
    b=np.column_stack([np.ones(len(base)),base.to_numpy()])
    b=(b-b.mean(0))/np.maximum(b.std(0),1e-8);b[:,0]=1
    yr=y-b@np.linalg.lstsq(b,y,rcond=None)[0]
    xr=x-b@np.linalg.lstsq(b,x,rcond=None)[0]
    return np.corrcoef(yr,xr)[0,1] if min(yr.std(),xr.std())>1e-10 else np.nan


def main(save_profile_assignments=False):
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    retention=[];outcomes=[];associations=[];profile_assignments=[]
    for stage in STAGES:
        folder=OUT/stage;qc=ROOT/'results/xenium_boundary_qc'/stage
        print('BOUNDARY_STAGE',stage,flush=True)
        atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        genes=data.var_names.to_numpy(dtype=str)
        fb=pd.read_csv(folder/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        env=pd.read_csv(folder/'CM_neighborhoods_50um.tsv.gz',sep='\t',index_col=0)
        native_fc=data.X[data.obs_names.get_indexer(fb.index)].toarray()
        model=FibroblastProjection(stage,native_fc,genes,ann)
        ftree=cKDTree(fb[['x_centroid','y_centroid']].to_numpy())
        cxy=cm[['x_centroid','y_centroid']].to_numpy()
        edge=cKDTree(cxy).sparse_distance_matrix(ftree,50,output_type='coo_matrix')
        neighborhood=sparse.csr_matrix((np.ones(len(edge.data)),(edge.row,edge.col)),shape=(len(cm),len(fb)))
        nn=np.asarray(neighborhood.sum(1)).ravel()
        with np.load(qc/'geometry_sensitivity_counts.npz') as b:
            ids=b['cell_ids'].astype(str);np.testing.assert_array_equal(b['genes'].astype(str),genes)
            geo={v:b[v] for v in ['native_Q20','core_1um']}
        with np.load(qc/'depth_matched_control_counts.npz') as b:
            np.testing.assert_array_equal(b['cell_ids'].astype(str),ids);geo['thin_to_core1']=b['thin_to_core1']
        cohorts=[('geometry',ids,geo)]
        with np.load(qc/'morphology_comparison_counts.npz') as b:
            np.testing.assert_array_equal(b['genes'].astype(str),genes)
            cohorts.append(('image',b['cell_ids'].astype(str),{v:b[v] for v in ['native_Q20','membrane_watershed']}))
        for scope,ids,variants in cohorts:
            fi=pd.Index(ids).get_indexer(fb.index);cf=fi>=0
            ci=pd.Index(ids).get_indexer(cm.index);cc=ci>=0
            audited_neighbors=np.asarray(neighborhood@cf.astype(float)).ravel()/np.maximum(nn,1)
            common=cc&(audited_neighbors>=.9)&(nn>=3)&env.resolved_neighbors.ge(5).to_numpy()
            for values in variants.values():
                selected=np.maximum(ci,0)
                common&=(values[selected].sum(1)>=30)&((values[selected]>0).sum(1)>=10)
            baseframe=cm.loc[common].join(env.loc[common]);base=baseline_features(baseframe) if common.sum()>100 else None
            for variant,values in variants.items():
                counts=values[fi[cf]]
                pred=model.predict(counts)
                classified=annotate(counts,genes,model,pred);classified.index=fb.index[cf]
                original=fb.loc[classified.index]
                if save_profile_assignments:
                    audit=classified.copy()
                    audit.index.name='cell_id';audit=audit.reset_index()
                    audit['stage']=stage;audit['scope']=scope;audit['variant']=variant
                    audit['native_mapped_family']=original.mapped_family.to_numpy()
                    audit['native_supported_family_assignment']=original.supported_family_assignment.to_numpy()
                    profile_assignments.append(audit)
                for family in FAMILIES+[UNRESOLVED]:
                    mask=original.mapped_family.eq(family)
                    if not mask.any():continue
                    retention.append({'stage':stage,'scope':scope,'variant':variant,'original_supported_family':family,
                        'profiles':int(mask.sum()),'information_QC_fraction':classified.loc[mask,'state_information_QC'].mean(),
                        'same_proposal_fraction':classified.loc[mask,'family_proposal'].eq(original.loc[mask,'family_proposal']).mean(),
                        'same_supported_label_fraction':classified.loc[mask,'mapped_family'].eq(original.loc[mask,'mapped_family']).mean()})
                stats,_=cm_scores(values[ci[cc]],genes)
                original_cm=cm.loc[cc]
                for name in ['maturation_21gene_mean','Ankrd1_log1pCP10k']:
                    y=stats[name].to_numpy();old=original_cm[name].to_numpy()
                    outcomes.append({'stage':stage,'scope':scope,'variant':variant,'outcome':name,'CM_profiles':len(y),
                        'native_variant_correlation':np.corrcoef(old,y)[0,1],'mean_variant_minus_native':np.mean(y-old),
                        'profiles_without_nuclear_area':int(original_cm.nuclear_area_missing.sum())})
                if base is None:continue
                patched_weights=fb[['weight__'+f for f in FAMILIES]].to_numpy().copy()
                patched_weights[cf]=pred['primary'][0][:,:len(FAMILIES)]
                # Information-limited profiles remain counted in the parent FB census with uncertain fine state.
                patched_weights[np.flatnonzero(cf)[~classified.state_information_QC.to_numpy()]]=0
                patched_programs=fb[[c for c in fb if c.startswith('matrix__')]].copy()
                newprog=direct_programs(counts,genes);newprog.index=classified.index
                patched_programs.loc[classified.index]=newprog
                exposures={}
                for j,family in enumerate(FAMILIES):exposures['family__'+family]=np.asarray(neighborhood@patched_weights[:,j]).ravel()/np.maximum(nn,1)
                for name in patched_programs:exposures[name]=np.asarray(neighborhood@patched_programs[name].to_numpy()).ravel()/np.maximum(nn,1)
                ys,_=cm_scores(values[ci[common]],genes)
                for outcome in ['maturation_21gene_mean','Ankrd1_log1pCP10k']:
                    y=ys[outcome].to_numpy()
                    for predictor,x in exposures.items():
                        associations.append({'stage':stage,'scope':scope,'variant':variant,'outcome':outcome,'predictor':predictor,
                            'same_CM_subset':int(common.sum()),'partial_correlation':partial_r(y,x[common],base),
                            'coverage_rule':'CM audited and>=90% of its 50um FB neighbors audited;>=3FB;all variants CM information-QC;fixed original covariates'})
        print('BOUNDARY_DONE',stage,flush=True)
    save_tsv(pd.DataFrame(retention),OUT/'boundary_fibroblast_state_retention.tsv')
    save_tsv(pd.DataFrame(outcomes),OUT/'boundary_CM_outcome_sensitivity.tsv')
    save_tsv(pd.DataFrame(associations),OUT/'boundary_spatial_association_sensitivity.tsv')
    if save_profile_assignments:
        pd.concat(profile_assignments,ignore_index=True).to_csv(OUT/'boundary_profile_assignments.tsv.gz',sep='\t',index=False,na_rep='NA')
    write_json(OUT/'boundary_methods.json',{'scope':'Existing transcript-level boundary-audited profiles; no new molecule acquisition',
        'comparisons':['Native vs1um core vs count-matched random thinning','Same image-covered native vs membrane-guided alternative'],
        'fine_reference':'Frozen original stage reference and query calibration; broad parent labels retained',
        'same_CM_spatial_subset':'>=90% of 50um FB neighbors audited,>=3FB,>=5resolved neighbors;CM information-QC in every variant',
        'covariates':'Original frozen covariates; changed measured outcomes and FB exposures evaluated on identical subsets',
        'limitations':'Alternative boundaries are sensitivity models,not anatomical truth; sampled ROIs not representative independent sections; partial correlations not directly comparable to full-stage joint coefficients'})


if __name__=='__main__':main()
