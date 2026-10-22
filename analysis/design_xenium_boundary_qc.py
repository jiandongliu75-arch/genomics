"""Freeze outcome-blind regions and boundary sensitivity settings."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from build_whole_heart_spatial_atlas import ROOT,STAGES
from reassess_whole_heart_lineages import OUT as FIT,digest

OUT=ROOT/'results/xenium_boundary_qc'
ATLAS=ROOT/'results/whole_heart_spatial_atlas'
DESIGN={'scope':'All4 stages; molecule-level boundary audit and segmentation-assignment sensitivity, not an animal-replicated validation',
 'nucleus_policy':'Do not reject or require a nucleus for CMs; preserve all cell seeds including profiles without measured nuclear area',
 'selection':'Per stage:2 arteriole-like segments with many pre-existing matched pairs and separated centers; ordinary/flagged CM–EC and CM–fibroblast interfaces;2 additional nuclear-area-missing CM interfaces when available. Selection excludes Ankrd1 outcomes.',
 'molecule_QV_threshold':20,'strict_QV_sensitivity':30,
 'geometric_variants':['native_Q20','core_1um','core_2um','dilation_2um','native_Q30'],
 'core_rule':'Keep originally assigned Q20 gene transcripts inside the deposited polygon and farther than1/2µm from its boundary. Nucleus not required.',
 'dilation_rule':'Preserve native Q20 assignments; add previously unassigned Q20 genes only within2µm of a unique nearest deposited cell polygon. Competing cell polygons are considered; do not steal assigned transcripts.',
 'morphology_variants':['centroid_Voronoi','membrane_watershed'],
 'morphology_grid_um':0.5,'morphology_mask':'Union of all deposited ROI cell polygons buffered2µm; constrain boundary testing to the observed cellular territory',
 'morphology_seeds':'One deposited cell-centroid seed per profile in the contextual ROI, including anucleate CM profiles; no independent cell-number inference',
 'membrane_cost':'Gaussian-smoothed ATP1A1/CD45/E-Cadherin intensity, robust1st–99th percentile scaling; seeded watershed compactness0.001, no RNA/lineage/outcome input',
 'context_margin_um':40,'interface_halfwidth_um':100,'arterial_margin_um':175,
 'boundary_shell_um':1,'CM_facing_boundary_distance_um':2,
 'expression_models':'Frozen reference/query gene scaling from the original full-stage dataset; primary, raw-reference and both RNA libraries; outcome genes excluded as before',
 'Ankrd1_analysis':'Recalculate outcomes in the same original3104 near/far CM pairs across all eligible vascular segments; report retained-informative-pair sensitivity and nuclear-area strata. Image-based resegmentation is restricted to selected ROIs.',
 'forbidden_inference':'No statement that coexpression proves contamination, that absence of a nucleus is an error, or that transcript purity validates new biological boundaries'}


def main():
    OUT.mkdir(exist_ok=True)
    path=OUT/'boundary_design.json'
    if path.exists():assert json.loads(path.read_text())==DESIGN
    else:path.write_text(json.dumps(DESIGN,indent=2)+'\n')
    rows=[];selection_sources=[]
    for stage in STAGES:
        obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        pairs=pd.read_csv(FIT/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
        pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
        features=json.loads((FIT/stage/'arterial_segments.geojson').read_text())['features']
        lookup={f['properties']['segment_id']:f for f in features if f['properties']['definition']=='independent_strict' and f['properties']['primary_small_segment']}
        chosen=[];cm=obs[obs.atlas_identity.eq('Cardiomyocyte')];tree=cKDTree(cm[['x_centroid','y_centroid']].to_numpy())
        for sid,n in pairs.groupby('segment_id').size().sort_values(ascending=False).items():
            f=lookup[sid];p=f['properties'];center=np.array([p['center_x'],p['center_y']])
            if chosen and np.min(np.linalg.norm(np.array(chosen)-center,axis=1))<300:continue
            from shapely.geometry import shape
            x0,y0,x1,y1=shape(f['geometry']).bounds;b=[x0-175,x1+175,y0-175,y1+175]
            rows.append({'stage':stage,'roi_id':stage+'_arterial_'+str(len(chosen)+1),'category':'Peri-arterial',
                'selection_cell_id':'NA','partner_CM_id':'NA','segment_id':sid,'original_matched_pairs':int(n),
                'x_min':max(0,b[0]),'x_max':b[1],'y_min':max(0,b[2]),'y_max':b[3],
                'selection_reason':'Many frozen matched pairs; source Ankrd1 not consulted'})
            chosen.append(center)
            if len(chosen)==2:break
        for label in ['Endothelial','Fibroblast']:
            for flagged in [False,True]:
                take=(obs.previous_reassessed_lineage.eq(label)&obs.atlas_identity.eq('Mixed-marker review')) if flagged else obs.atlas_identity.eq(label)
                part=obs[take];d,ix=tree.query(part[['x_centroid','y_centroid']].to_numpy())
                candidates=np.flatnonzero(d<=20);order=np.random.default_rng(20261022).permutation(candidates)
                selected=None
                for i in order:
                    xy=part[['x_centroid','y_centroid']].iloc[i].to_numpy()
                    if chosen and np.min(np.linalg.norm(np.array(chosen)-xy,axis=1))<200:continue
                    selected=(i,xy);break
                if selected is None:
                    assert len(order)>0,(stage,label,flagged);i=order[0];selected=(i,part[['x_centroid','y_centroid']].iloc[i].to_numpy())
                i,xy=selected;chosen.append(xy)
                rows.append({'stage':stage,'roi_id':stage+'_'+label+'_'+('flagged' if flagged else 'ordinary'),
                    'category':'CM–'+label+' interface','selection_cell_id':part.index[i],'partner_CM_id':cm.index[ix[i]],
                    'segment_id':'NA','original_matched_pairs':0,'x_min':max(0,xy[0]-100),'x_max':xy[0]+100,
                    'y_min':max(0,xy[1]-100),'y_max':xy[1]+100,
                    'selection_reason':'Mixed-marker profile' if flagged else 'Ordinary supported non-CM profile'})
        non=obs[obs.atlas_identity.isin(['Endothelial','Fibroblast'])]
        no=cm[cm.nucleus_area.isna()];ntree=cKDTree(non[['x_centroid','y_centroid']].to_numpy())
        d,ix=ntree.query(no[['x_centroid','y_centroid']].to_numpy());candidates=np.flatnonzero(d<=20)
        added=0
        for i in np.random.default_rng(741).permutation(candidates):
            xy=no[['x_centroid','y_centroid']].iloc[i].to_numpy()
            if chosen and np.min(np.linalg.norm(np.array(chosen)-xy,axis=1))<200:continue
            added+=1;chosen.append(xy)
            rows.append({'stage':stage,'roi_id':stage+'_anucleate_CM_'+str(added),'category':'CM without measured nuclear area',
                'selection_cell_id':no.index[i],'partner_CM_id':non.index[ix[i]],'segment_id':'NA','original_matched_pairs':0,
                'x_min':max(0,xy[0]-100),'x_max':xy[0]+100,'y_min':max(0,xy[1]-100),'y_max':xy[1]+100,
                'selection_reason':'Nuclear-area-missing CM adjacent to a supported non-CM; retained as a legitimate candidate CM profile'})
            if added==2:break
        selection_sources.append({'stage':stage,'atlas_annotations_sha256':digest(ATLAS/stage/'all_cell_annotations.tsv.gz'),
            'matched_pairs_sha256':digest(FIT/stage/'default_contrast_matched_pairs.tsv.gz'),
            'vascular_polygons_sha256':digest(FIT/stage/'arterial_segments.geojson')})
    pd.DataFrame(rows).to_csv(OUT/'selected_ROIs.tsv',sep='\t',index=False,na_rep='NA')
    (OUT/'selection_provenance.json').write_text(json.dumps(selection_sources,indent=2)+'\n')
    print('FROZEN_ROIS',len(rows),flush=True)


if __name__=='__main__':main()
