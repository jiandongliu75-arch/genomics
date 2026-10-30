"""Describe broad-cell neighbors of FB states and direct within-stage matrix pairs."""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from fibroblast_state_common import *


def main():
    pairrows=[]
    for layer in ['corrected','raw']:
        frame=pd.read_csv(OUT/f'RNA_{layer}_matrix_pseudobulk_scores.tsv',sep='\t')
        frame=frame[frame.version.eq('annotation_markers_withheld')]
        for (stage,program),part in frame.groupby(['age','program'],observed=True):
            for j,a in enumerate(FAMILIES[:4]):
                for b in FAMILIES[:4][j+1:]:
                    aa=part[part.family.eq(a)].set_index('sample');bb=part[part.family.eq(b)].set_index('sample')
                    libs=aa.index.intersection(bb.index)
                    if len(libs)!=2:continue
                    delta=aa.loc[libs,'score_mean_log1p_CPM']-bb.loc[libs,'score_mean_log1p_CPM']
                    pairrows.append({'layer':layer,'stage':stage,'program':program,'family_A':a,'family_B':b,
                        'library1_A_minus_B':delta.iloc[0],'library2_A_minus_B':delta.iloc[1],
                        'mean_A_minus_B':delta.mean(),'both_library_directions_agree':np.sign(delta.iloc[0])==np.sign(delta.iloc[1])})
    save_tsv(pd.DataFrame(pairrows),OUT/'RNA_matrix_direct_family_pairs.tsv')
    rows=[]
    for stage in STAGES:
        folder=OUT/stage
        atlas=pd.read_csv(ROOT/'results/whole_heart_spatial_atlas'/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        fb=pd.read_csv(folder/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        dist=pd.read_csv(ROOT/'results/whole_heart_Ankrd1_vessels'/stage/'arterial_boundary_distances.tsv.gz',sep='\t',index_col=0)
        xy=fb[['x_centroid','y_centroid']].to_numpy()
        resolved=atlas.atlas_resolved
        total=cKDTree(atlas.loc[resolved,['x_centroid','y_centroid']].to_numpy()).query_ball_point(xy,30,return_length=True)-1
        d=pd.DataFrame(index=fb.index);d['resolved_other_neighbors']=total
        d['arterial_boundary_distance_um']=dist.loc[fb.index,'independent_strict__boundary_distance_um']
        for lineage in ['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']:
            m=atlas.atlas_identity.eq(lineage)
            num=cKDTree(atlas.loc[m,['x_centroid','y_centroid']].to_numpy()).query_ball_point(xy,30,return_length=True)
            if lineage=='Fibroblast':num-=1
            d['neighbor_fraction__'+lineage]=num/np.maximum(total,1)
        d.to_csv(folder/'FB_broad_neighbor_context_30um.tsv.gz',sep='\t',na_rep='NA')
        for family in FAMILIES:
            for kind in ['soft_reference_weight','supported_labels']:
                w=fb['weight__'+family].to_numpy() if kind.startswith('soft') else fb.mapped_family.eq(family).to_numpy(dtype=float)
                valid=(total>=5)&np.isfinite(d.arterial_boundary_distance_um)
                w=np.where(valid,w,0)
                if w.sum()<1:continue
                for col in d.columns:
                    if col=='resolved_other_neighbors':continue
                    rows.append({'stage':stage,'family':family,'summary':kind,'measure':col,
                        'sum_weights_or_supported_profiles':w.sum(),'weighted_mean':np.sum(d[col].fillna(0).to_numpy()*w)/w.sum(),
                        'all_FB_mean':d.loc[valid,col].mean(),
                        'interpretation':'descriptive native-coordinate neighborhood; soft weights are not calibrated identities; no physical-contact claim'})
    save_tsv(pd.DataFrame(rows),OUT/'FB_state_broad_spatial_context.tsv')


if __name__=='__main__':main()
