"""Freeze matched CM pairs and condition their Ankrd1 contrast on FB environments."""
import numpy as np
import pandas as pd
from fibroblast_state_common import *


def main():
    records=[]
    for stage in STAGES:
        folder=OUT/stage
        cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        env=pd.read_csv(folder/'CM_neighborhoods_50um.tsv.gz',sep='\t',index_col=0)
        d=cm.join(env)
        ok=d.fibroblast_neighbors.ge(3)&d.resolved_neighbors.ge(5)
        pairs=pd.read_csv(ROOT/'results/whole_heart_Ankrd1_vessels'/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
        pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
        pairs=pairs[pairs.near_cell_id.isin(d.index[ok])&pairs.far_cell_id.isin(d.index[ok])]
        n=d.loc[pairs.near_cell_id].reset_index(drop=True);f=d.loc[pairs.far_cell_id].reset_index(drop=True)
        y=n.Ankrd1_log1pCP10k.to_numpy()-f.Ankrd1_log1pCP10k.to_numpy()
        nuisance=pd.DataFrame()
        for col in ['outcome_excluded_depth','cell_area','endocardial_proxy_distance_um','resolved_neighbors']:
            nuisance[col]=np.log1p(n[col])-np.log1p(f[col])
        for col in ['fibroblast_fraction','endothelial_fraction','mural_fraction','immune_fraction','Slit2_log','cycling_log','maturation_21gene_mean']:
            nuisance[col]=n[col]-f[col]
        family_cols=list(identifiable_family_exposures(d[[c for c in d if c.startswith('family__')]]))
        matrix_cols=[c for c in d if c.startswith('matrix__')]
        sets={'frozen_pairs_unadjusted':pd.DataFrame(index=range(len(pairs))),
              'matched_covariates_plus_CM_maturation':nuisance,
              'plus_FB_families':pd.concat([nuisance,n[family_cols]-f[family_cols]],axis=1),
              'plus_FB_matrix':pd.concat([nuisance,n[matrix_cols]-f[matrix_cols]],axis=1),
              'plus_FB_families_and_matrix':pd.concat([nuisance,n[family_cols+matrix_cols]-f[family_cols+matrix_cols]],axis=1)}
        groups=pairs.segment_id.to_numpy();unique,inverse=np.unique(groups,return_inverse=True)
        for name,frame in sets.items():
            vals=frame.to_numpy(dtype=float);sd=vals.std(0);vals=vals[:,sd>1e-8]/sd[sd>1e-8]
            # No centering: intercept is the contrast at zero near-far covariate difference.
            x=np.column_stack([np.ones(len(y)),vals])
            beta=np.linalg.pinv(x.T@x,rcond=1e-9)@(x.T@y)
            gram=np.zeros((len(unique),x.shape[1],x.shape[1]));rhs=np.zeros((len(unique),x.shape[1]))
            for j in range(len(unique)):
                m=inverse==j;gram[j]=x[m].T@x[m];rhs[j]=x[m].T@y[m]
            rng=np.random.default_rng(718);est=[]
            for _ in range(192):
                counts=np.bincount(rng.integers(0,len(unique),len(unique)),minlength=len(unique))
                est.append((np.linalg.pinv(np.einsum('b,bij->ij',counts,gram),rcond=1e-9)@np.einsum('b,bi->i',counts,rhs))[0])
            records.append({'stage':stage,'model':name,'same_matched_pairs':len(pairs),'vascular_segments':len(unique),
                'native_near_far_CP10k_ratio_same_subset':np.expm1(n.Ankrd1_log1pCP10k).mean()/np.expm1(f.Ankrd1_log1pCP10k).mean(),
                'conditional_near_minus_far_log_expression':beta[0],
                'segment_bootstrap_low':np.quantile(est,.025),'segment_bootstrap_high':np.quantile(est,.975),
                'interpretation':'conditional observational contrast in the same pairs; not causal or independent-animal evidence'})
    save_tsv(pd.DataFrame(records),OUT/'Ankrd1_vascular_pairs_after_FB_adjustment.tsv')


if __name__=='__main__':main()
