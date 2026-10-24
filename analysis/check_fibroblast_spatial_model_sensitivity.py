"""Test soft-transfer, supported-label, outcome-definition and vascular subsets."""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from fibroblast_state_common import *
from compare_fibroblast_matrix_and_spatial import baseline_features,coefficient_bootstrap,blocked_cv


def main():
    coefficients=[];predictions=[]
    for stage in STAGES:
        folder=OUT/stage
        cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        env=pd.read_csv(folder/'CM_neighborhoods_50um.tsv.gz',sep='\t',index_col=0)
        fb=pd.read_csv(folder/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        xy=cm[['x_centroid','y_centroid']].to_numpy()
        edge=cKDTree(xy).sparse_distance_matrix(cKDTree(fb[['x_centroid','y_centroid']].to_numpy()),50,output_type='coo_matrix')
        m=sparse.csr_matrix((np.ones(len(edge.data)),(edge.row,edge.col)),shape=(len(cm),len(fb)))
        nf=np.asarray(m.sum(1)).ravel()
        other=pd.DataFrame(index=cm.index)
        for fam in FAMILIES:
            other['family__'+fam]=np.asarray(m@fb['matrix_withheld_weight__'+fam].to_numpy()).ravel()/np.maximum(nf,1)
        d=cm.join(env)
        valid=d.fibroblast_neighbors.ge(3)&d.resolved_neighbors.ge(5)&d.arterial_distance_um.notna()
        d=d.loc[valid];base=baseline_features(d)
        families=d[[c for c in d if c.startswith('family__')]].copy()
        withheld=other.loc[d.index].copy()
        hard=d[[c for c in d if c.startswith('hardfamily__')]].copy()
        hard.columns=[c.replace('hardfamily__','family__') for c in hard]
        families,withheld,hard=[identifiable_family_exposures(frame) for frame in [families,withheld,hard]]
        matrix=d[[c for c in d if c.startswith('matrix__')]]
        bx=np.floor(d.x_centroid.to_numpy()/400).astype(int);by=np.floor(d.y_centroid.to_numpy()/400).astype(int)
        blocks=bx*100000+by
        jobs=[('mapping_ECM_genes_withheld',withheld,'maturation_21gene_mean',np.ones(len(d),bool)),
              ('supported_labels_only',hard,'maturation_21gene_mean',np.ones(len(d),bool)),
              ('mapping_ECM_genes_withheld',withheld,'Ankrd1_log1pCP10k',np.ones(len(d),bool)),
              ('supported_labels_only',hard,'Ankrd1_log1pCP10k',np.ones(len(d),bool)),
              ('outcome_denominator_excludes_25_genes',matrix,'maturation_excluded_denominator',np.ones(len(d),bool)),
              ('outcome_full_denominator',matrix,'Ankrd1_full_denominator',np.ones(len(d),bool)),
              ('within_150um_arterial_envelope',pd.concat([families,matrix],axis=1),'Ankrd1_log1pCP10k',d.arterial_distance_um.le(150).to_numpy())]
        for name,extra,outcome,mask in jobs:
            if mask.sum()<300:continue
            y=d.loc[mask,outcome].to_numpy()
            for row in coefficient_bootstrap(y,base.loc[mask],extra.loc[mask],blocks[mask]):
                coefficients.append({'stage':stage,'sensitivity':name,'outcome':outcome,'CM_profiles':int(mask.sum()),**row})
            if name!='within_150um_arterial_envelope':
                summary,_=blocked_cv(y,base.loc[mask],extra.loc[mask],d.loc[mask])
                predictions.append({'stage':stage,'sensitivity':name,'outcome':outcome,**summary})
        print('MODEL_SENSITIVITY_DONE',stage,flush=True)
    save_tsv(pd.DataFrame(coefficients),OUT/'spatial_model_sensitivity_coefficients.tsv')
    save_tsv(pd.DataFrame(predictions),OUT/'spatial_model_sensitivity_prediction.tsv')


if __name__=='__main__':main()
