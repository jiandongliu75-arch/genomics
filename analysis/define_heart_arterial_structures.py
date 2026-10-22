"""Construct EC-plus-SMC vascular segments without using Ankrd1 outcomes."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN
import shapely
from shapely.geometry import Polygon,shape,mapping
from shapely.strtree import STRtree

from reassess_whole_heart_lineages import OUT,SOURCE,ROOT,STAGES,digest

RAW=Path('/workspace/heart-atlas/data/xenium')
DEFINITIONS=['independent_strict','independent_marker','published_comparator']


def clusters(obs,definition):
    if definition=='independent_strict':ec,smc=obs.arterial_EC_seed,obs.SMC_seed
    elif definition=='independent_marker':ec,smc=obs.arterial_EC_marker_seed,obs.SMC_marker_seed
    else:ec,smc=obs.author_atlas_cell_state.eq('arteriole.EC'),obs.author_atlas_cell_state.eq('SMC')
    arterial=obs.loc[ec,['x_centroid','y_centroid']];muscle=obs.loc[smc,['x_centroid','y_centroid']]
    if len(arterial)<2 or len(muscle)<2:return [],{'arterial_seed_cells':len(arterial),'SMC_seed_cells':len(muscle),'accepted_segments':0}
    labels=DBSCAN(eps=35,min_samples=2).fit_predict(arterial.to_numpy())
    valid=labels>=0;arterial=arterial.iloc[np.flatnonzero(valid)];labels=labels[valid]
    if not len(arterial):return [],{'arterial_seed_cells':int(ec.sum()),'SMC_seed_cells':len(muscle),'accepted_segments':0}
    # Assign a mural seed to its nearest clustered arterial EC only, preventing
    # the same vascular cell from being used in two objects of one definition.
    distance,nearest=cKDTree(arterial.to_numpy()).query(muscle.to_numpy())
    assignment=labels[nearest]
    accepted=[];rejected=0
    for cluster in sorted(set(labels)):
        a=arterial.iloc[np.flatnonzero(labels==cluster)]
        s=muscle.iloc[np.flatnonzero((assignment==cluster)&(distance<=35))]
        if len(s)<2:rejected+=1;continue
        adjacent,_=cKDTree(s.to_numpy()).query(a.to_numpy())
        coverage=float((adjacent<=35).mean())
        center=a.to_numpy().mean(0);angles=np.arctan2(s.y_centroid-center[1],s.x_centroid-center[0])
        sectors=np.unique(np.floor((angles+np.pi)/(np.pi/2)).astype(int)%4)
        if coverage<.5 or len(sectors)<2:rejected+=1;continue
        accepted.append({'segment_id':f'{definition}_{len(accepted):03d}','definition':definition,
            'arterial_EC_ids':a.index.tolist(),'SMC_ids':s.index.tolist(),'arterial_ECs':len(a),'SMC_cells':len(s),
            'arterial_EC_mural_coverage':coverage,'mural_angular_sectors':len(sectors),
            'center_x':float(center[0]),'center_y':float(center[1])})
    return accepted,{'arterial_seed_cells':int(ec.sum()),'SMC_seed_cells':len(muscle),'accepted_segments':len(accepted),'rejected_clusters_mural_support':rejected}


def read_boundaries(stage,wanted,obs):
    raw_stage='P07' if stage=='P7' else stage
    candidates=list(RAW.glob(f'*{raw_stage}__*/cell_boundaries.csv.gz'))
    if len(candidates)!=1:
        candidates=[p for p in RAW.glob('*/cell_boundaries.csv.gz') if p.parent.name.split('__')[0]==f'WT_{raw_stage}']
    assert len(candidates)==1,(stage,candidates)
    path=candidates[0];frames=[]
    for chunk in pd.read_csv(path,usecols=['cell_id','vertex_x','vertex_y'],chunksize=200000):
        keep=chunk.cell_id.isin(wanted)
        if keep.any():frames.append(chunk.loc[keep])
    polygons={}
    if frames:
        vertices=pd.concat(frames,ignore_index=True)
        for cell,part in vertices.groupby('cell_id',sort=False):
            coords=part[['vertex_x','vertex_y']].to_numpy()
            if len(coords)>=3:
                polygon=shapely.make_valid(Polygon(coords))
                if not polygon.is_empty and polygon.area>0:polygons[cell]=polygon
    assert set(wanted).issubset(polygons),f'Missing measured boundaries: {set(wanted)-set(polygons)}'
    # Native polygons and cell centroids must refer to the same µm coordinate
    # system; modest centroid/boundary discrepancies are explicitly retained.
    distances=[polygons[c].distance(shapely.Point(*obs.loc[c,['x_centroid','y_centroid']])) for c in wanted]
    return polygons,{'source_boundary_path':str(path),'source_boundary_sha256':digest(path),
        'vascular_cells_with_native_boundaries':len(polygons),'median_centroid_outside_own_polygon_um':float(np.median(distances)) if distances else 0,
        'maximum_centroid_outside_own_polygon_um':float(max(distances)) if distances else 0}


def nearest(points,geometries):
    distances=np.full(len(points),np.nan);ids=np.full(len(points),-1,int)
    if geometries:
        indices,values=STRtree(geometries).query_nearest(points,all_matches=False,return_distance=True)
        distances[indices[0]]=values;ids[indices[0]]=indices[1]
    return distances,ids


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stages',nargs='+',default=STAGES,choices=STAGES);args=parser.parse_args()
    design={'arterial_seed_cluster_eps_um':35,'minimum_arterial_ECs':2,'minimum_SMCs':2,
        'SMC_distance_to_clustered_arterial_EC_um':35,'minimum_EC_mural_coverage':.5,'minimum_mural_angular_quadrants':2,
        'primary_boundary':'Union of deposited arterial EC and assigned SMC cell polygons, buffered3µm; vascular-cell envelope, not histologically traced vessel wall',
        'primary_small_segment':'Convex-hull equivalent diameter<=150µm and maximum boundary span<=300µm',
        'sensitivity_boundaries':['unbuffered vascular-cell union','convex hull','clustered arterial EC centroid'],
        'definitions':DEFINITIONS,'outcome_expression_used_to_define_vessels':False,
        'interpretation':'Arteriole-like candidate vascular segments; a fragmented vessel can form multiple segments; no independently validated anatomical segmentation.'}
    if (OUT/'vessel_design.json').exists():assert json.loads((OUT/'vessel_design.json').read_text())==design
    else:(OUT/'vessel_design.json').write_text(json.dumps(design,indent=2)+'\n')
    for stage in args.stages:
        folder=OUT/stage
        if (folder/'vessel_completion.json').exists():print('SKIP_VESSELS',stage,flush=True);continue
        if not (folder/'lineage_completion.json').exists():raise FileNotFoundError(f'Complete lineage mapping before {stage}')
        start=time.time();obs=pd.read_csv(folder/'whole_heart_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        objects=[];counts={}
        for definition in DEFINITIONS:
            current,summary=clusters(obs,definition);objects.extend(current);counts[definition]=summary
        wanted=set(c for o in objects for key in ['arterial_EC_ids','SMC_ids'] for c in o[key])
        polygons,provenance=read_boundaries(stage,wanted,obs) if wanted else ({},{'vascular_cells_with_native_boundaries':0})
        geo=[];records=[]
        for o in objects:
            cellids=o['arterial_EC_ids']+o['SMC_ids']
            union=shapely.union_all([polygons[c] for c in cellids]);envelope=union.buffer(3);hull=union.convex_hull
            vertices=np.asarray(hull.exterior.coords)
            span=float(np.max(np.linalg.norm(vertices[:,None,:]-vertices[None,:,:],axis=2)))
            equivalent=float(2*np.sqrt(hull.area/np.pi))
            primary=equivalent<=150 and span<=300
            props={**o,'stage':stage,'vascular_cell_union_area_um2':float(union.area),'convex_hull_area_um2':float(hull.area),
                'hull_equivalent_diameter_um':equivalent,'maximum_span_um':span,'primary_small_segment':bool(primary)}
            geo.append({'type':'Feature','properties':props,'geometry':mapping(envelope),
                'unbuffered_geometry':mapping(union),'hull_geometry':mapping(hull)})
            records.append({k:v for k,v in props.items() if k not in ['arterial_EC_ids','SMC_ids']})
        (folder/'arterial_segments.geojson').write_text(json.dumps({'type':'FeatureCollection','features':geo})+'\n')
        pd.DataFrame(records).to_csv(folder/'arterial_segment_summary.tsv',sep='\t',index=False)
        pts=shapely.points(obs[['x_centroid','y_centroid']].to_numpy())
        output=obs[['x_centroid','y_centroid']].copy()
        for definition in DEFINITIONS:
            features=[f for f in geo if f['properties']['definition']==definition and f['properties']['primary_small_segment']]
            for geomkey,suffix in [('geometry','boundary'),('unbuffered_geometry','unbuffered'),('hull_geometry','hull')]:
                d,ix=nearest(pts,[shape(f[geomkey]) for f in features])
                names=np.array([f['properties']['segment_id'] for f in features]+['No vascular segment'])
                output[definition+'__'+suffix+'_distance_um']=d
                output[definition+'__'+suffix+'_segment_id']=names[np.where(ix>=0,ix,len(features))]
            ecids=[c for f in features for c in f['properties']['arterial_EC_ids']]
            ecowners=np.array([f['properties']['segment_id'] for f in features for c in f['properties']['arterial_EC_ids']])
            if ecids:
                d,idx=cKDTree(obs.loc[ecids,['x_centroid','y_centroid']].to_numpy()).query(obs[['x_centroid','y_centroid']].to_numpy())
                owner=ecowners[idx]
            else:d=np.full(len(obs),np.nan);owner=np.full(len(obs),'No vascular segment')
            output[definition+'__EC_centroid_distance_um']=d
            output[definition+'__EC_centroid_segment_id']=owner
            counts[definition]['primary_small_segments']=len(features)
        output.to_csv(folder/'arterial_boundary_distances.tsv.gz',sep='\t',na_rep='NA')
        (folder/'vessel_completion.json').write_text(json.dumps({'stage':stage,'definition_counts':counts,'segments_before_size_filter':len(objects),
            'boundary_provenance':provenance,'lineage_annotations_sha256':digest(folder/'whole_heart_cell_annotations.tsv.gz'),
            'elapsed_seconds':time.time()-start,'outcome_used_in_vessel_definition':False},indent=2)+'\n')
        print('VESSELS_DONE',stage,counts,flush=True)


if __name__=='__main__':main()
