"""Expose disagreement between focal-CM and source-EC annotation systems."""
from reconstruct_cardiac_spatial import OUT,UP,STAGES,read,save
import numpy as np
import pandas as pd

def run():
    rows=[]
    for stage in STAGES:
        q=read(OUT/stage/'all_profile_observations.tsv.gz');c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');e=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz');f=read(UP/stage/'atlas_annotations.tsv.gz')
        sources={'author_capillary':q.author_atlas_cell_state.eq('cap.EC'),
            'author_capillary_plus_HighMyo':q.author_atlas_cell_state.isin(['cap.EC','cap.EC.High.Myo.']),
            'QC_capillary':q.atlas_identity.eq('Endothelial')&q.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like'),
            'author_pericyte':q.author_atlas_cell_state.eq('Peri')}
        for name,source in sources.items():
            conflict=source.reindex(c.index).to_numpy(bool);e['source_role_conflict__'+name]=conflict
            for key in e:
                if key.startswith(name+'__') and ('apposition_fraction' in key or 'minimum_boundary_distance' in key):e.loc[conflict,key]=np.nan
            for cohort,mask in [('paper_label_mask',c.reconstructed_Figure4e_CM),('our_ventricular_mask',c.broader_ventricular_CM)]:
                rows.append({'stage':stage,'cohort':cohort,'source_definition':name,'focal_CMs':int(mask.sum()),'source_role_conflicts':int((mask.to_numpy()&conflict).sum()),
                    'action':'Retained in atlas; excluded from this source-specific CM boundary comparison to prevent self-apposition. No relabeling or proof of which identity is correct.'})
        save(e,f'{stage}/vascular_interface_and_matrix_context.tsv.gz',index=True)
    save(pd.DataFrame(rows),'focal_source_identity_conflicts.tsv')

if __name__=='__main__':run()
