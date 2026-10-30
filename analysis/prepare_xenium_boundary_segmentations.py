"""Morphology-informed fixed-seed boundary sensitivity, preserving anucleate CMs."""
import argparse
import gzip
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
import shapely
from shapely.geometry import mapping,Polygon
from skimage.draw import polygon as raster_polygon
from skimage.measure import find_contours
from skimage.segmentation import watershed
import tifffile

from design_xenium_boundary_qc import OUT,ATLAS,STAGES,DESIGN
from acquire_xenium_boundary_inputs import CACHE
from define_heart_arterial_structures import read_boundaries,RAW
from reassess_whole_heart_lineages import digest


def pieces(geometry):
    if geometry.geom_type=='Polygon':return [geometry]
    return [g for g in geometry.geoms if g.geom_type=='Polygon']


def rasterize(geometries,shape,x0,y0,res):
    result=np.zeros(shape,np.int32);conflict=0
    for label,geometry in enumerate(geometries,1):
        for part in pieces(geometry):
            xy=np.asarray(part.exterior.coords);r,c=raster_polygon((xy[:,1]-y0)/res,(xy[:,0]-x0)/res,shape)
            conflict+=int((result[r,c]>0).sum());result[r,c]=label
            for hole in part.interiors:
                xy=np.asarray(hole.coords);r,c=raster_polygon((xy[:,1]-y0)/res,(xy[:,0]-x0)/res,shape);result[r,c]=0
    return result,conflict


def vectorize(labels,cell_ids,x0,y0,res):
    result=[]
    for label,slices in enumerate(ndi.find_objects(labels),1):
        if slices is None:continue
        patch=np.pad(labels[slices]==label,1)
        contours=find_contours(patch,.5)
        polygons=[]
        for contour in contours:
            if len(contour)<4:continue
            y=(contour[:,0]+slices[0].start-1)*res+y0
            x=(contour[:,1]+slices[1].start-1)*res+x0
            geometry=shapely.make_valid(Polygon(np.column_stack([x,y])))
            if geometry.area>0:polygons.append(geometry)
        if polygons:
            # XOR of nested contours retains holes; union would incorrectly fill them.
            geometry=shapely.symmetric_difference_all(polygons)
            result.append({'type':'Feature','properties':{'cell_id':str(cell_ids[label-1]),'grid_label':label},'geometry':mapping(geometry)})
    return {'type':'FeatureCollection','features':result}


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=STAGES,required=True);args=p.parse_args();stage=args.stage
    folder=OUT/stage;folder.mkdir(exist_ok=True)
    if (folder/'morphology_completion.json').exists():print('SKIP_MORPHOLOGY',stage,flush=True);return
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t');rois=rois[rois.stage.eq(stage)]
    member=CACHE/stage/'morphology_focus_0001.ome.tif'
    alias='P07' if stage=='P7' else stage
    raw=next(f for f in RAW.iterdir() if f.name.split('__')[0]=='WT_'+alias)
    dapi_path=raw/'morphology_focus_0000.ome.tif'
    px=float(json.loads((raw/'experiment.xenium').read_text())['pixel_size'])
    # Read only page0 of each physical file. OME series would combine missing files.
    with tifffile.TiffFile(member,is_ome=False) as image:
        membrane=image.pages[0].asarray(maxworkers=2)
        xml=ET.fromstring(image.pages[0].description)
        channels=[c.attrib['Name'] for c in xml.iter() if c.tag.endswith('Channel')]
    assert channels[1]=='ATP1A1/CD45/E-Cadherin',channels
    with tifffile.TiffFile(dapi_path,is_ome=False) as image:dapi=image.pages[0].asarray(maxworkers=2)
    assert membrane.shape==dapi.shape and membrane.ndim==2
    needed=set()
    for roi in rois.itertuples():
        take=obs.x_centroid.between(max(0,roi.x_min-41),roi.x_max+41)&obs.y_centroid.between(max(0,roi.y_min-41),roi.y_max+41)
        needed.update(obs.index[take])
    polygons,provenance=read_boundaries(stage,needed,obs)
    nuclei=pd.read_parquet(CACHE/stage/'nucleus_boundaries.parquet')
    nuclei=nuclei[nuclei.cell_id.isin(needed)]
    summaries=[];roi_meta=[];seed_records=[]
    for roi in rois.itertuples():
        x0=max(0,roi.x_min-40);x1=min(membrane.shape[1]*px,roi.x_max+40)
        y0=max(0,roi.y_min-40);y1=min(membrane.shape[0]*px,roi.y_max+40);res=.5
        x0=np.floor(x0/res)*res;y0=np.floor(y0/res)*res
        shape=(int(np.ceil((y1-y0)/res))+1,int(np.ceil((x1-x0)/res))+1)
        yy,xx=np.meshgrid(y0+np.arange(shape[0])*res,x0+np.arange(shape[1])*res,indexing='ij')
        rows=np.clip(np.rint(yy/px).astype(int),0,membrane.shape[0]-1);cols=np.clip(np.rint(xx/px).astype(int),0,membrane.shape[1]-1)
        mem=membrane[rows,cols];dna=dapi[rows,cols]
        chosen=obs[obs.x_centroid.between(x0,x1)&obs.y_centroid.between(y0,y1)]
        ids=chosen.index.to_numpy(dtype=str);geometries=[polygons[c] for c in ids]
        native,overlap_pixels=rasterize(geometries,shape,x0,y0,res)
        mask=ndi.binary_dilation(native>0,iterations=4)
        markers=np.zeros(shape,np.int32)
        for label,(cell,o) in enumerate(chosen.iterrows(),1):
            r=int(round((o.y_centroid-y0)/res));c=int(round((o.x_centroid-x0)/res))
            original=(r,c)
            if markers[r,c]>0:
                candidates=[]
                for dr in range(-4,5):
                    for dc in range(-4,5):
                        nr,nc=r+dr,c+dc
                        if 0<=nr<shape[0] and 0<=nc<shape[1] and markers[nr,nc]==0 and native[nr,nc]==label:
                            candidates.append((dr*dr+dc*dc,nr,nc))
                assert candidates,(stage,roi.roi_id,cell)
                _,r,c=min(candidates)
            markers[r,c]=label;mask[r,c]=True
            seed_records.append({'stage':stage,'roi_id':roi.roi_id,'cell_id':cell,'grid_label':label,
                'seed_x_um':x0+c*res,'seed_y_um':y0+r*res,
                'rounding_or_collision_shift_um':float(np.hypot(x0+c*res-o.x_centroid,y0+r*res-o.y_centroid)),
                'nuclear_area_missing':bool(pd.isna(o.nucleus_area)),'nucleus_required':False})
        nearest=ndi.distance_transform_edt(markers==0,return_distances=False,return_indices=True)
        voronoi=markers[tuple(nearest)];voronoi[~mask]=0
        positive=mem[mask];lo,hi=np.quantile(positive,[.01,.99])
        cost=ndi.gaussian_filter(np.clip((mem.astype(np.float32)-lo)/max(hi-lo,1),0,1),sigma=.8)
        segmented=watershed(cost,markers=markers,mask=mask,compactness=.001).astype(np.int32)
        assert np.unique(markers[markers>0]).size==len(ids)
        assert set(np.unique(segmented))-set([0])==set(range(1,len(ids)+1))
        np.savez_compressed(folder/(roi.roi_id+'_segmentation_maps.npz'),cell_ids=ids,
            native_grid=native,centroid_Voronoi=voronoi,membrane_watershed=segmented,
            membrane=mem,DAPI=dna,x0_um=x0,y0_um=y0,grid_um=res,
            inner_bounds_um=np.array([roi.x_min,roi.x_max,roi.y_min,roi.y_max]))
        for name,labels in [('centroid_Voronoi',voronoi),('membrane_watershed',segmented)]:
            features=vectorize(labels,ids,x0,y0,res)
            with gzip.open(folder/(roi.roi_id+'_'+name+'_polygons.geojson.gz'),'wt') as handle:json.dump(features,handle)
            # Sparse integer contingency table yields source-raster overlap/IoU.
            stride=len(ids)+1;key=native.astype(np.int64)*stride+labels
            unique,n=np.unique(key,return_counts=True)
            old_area=np.bincount(native.ravel(),minlength=stride);new_area=np.bincount(labels.ravel(),minlength=stride)
            diagonal=np.zeros(stride,np.int64);a=unique//stride;b=unique%stride
            same=a==b;diagonal[a[same]]=n[same]
            for i,(cell,o) in enumerate(chosen.iterrows(),1):
                core=roi.x_min<=o.x_centroid<=roi.x_max and roi.y_min<=o.y_centroid<=roi.y_max
                summaries.append({'stage':stage,'roi_id':roi.roi_id,'variant':name,'cell_id':cell,
                    'atlas_identity':o.atlas_identity,'inner_evaluation_cell':core,'nuclear_area_missing':bool(pd.isna(o.nucleus_area)),
                    'source_raster_area_um2':float(old_area[i]*res*res),'variant_area_um2':float(new_area[i]*res*res),
                    'intersection_um2':float(diagonal[i]*res*res),'IoU':float(diagonal[i]/max(old_area[i]+new_area[i]-diagonal[i],1)),
                    'source_area_retained_fraction':float(diagonal[i]/max(old_area[i],1))})
        original={'type':'FeatureCollection','features':[{'type':'Feature','properties':{'cell_id':c},'geometry':mapping(polygons[c])} for c in ids]}
        with gzip.open(folder/(roi.roi_id+'_native_polygons.geojson.gz'),'wt') as handle:json.dump(original,handle)
        nuclei[nuclei.cell_id.isin(ids)].to_parquet(folder/(roi.roi_id+'_nuclear_boundaries.parquet'),index=False,compression='zstd')
        roi_meta.append({'roi_id':roi.roi_id,'context_bounds_um':[x0,x1,y0,y1],'grid_shape_yx':list(shape),'grid_um':res,
            'cell_seeds':len(ids),'native_raster_overlap_pixels':overlap_pixels,'membrane_scaling_1_99':[float(lo),float(hi)]})
        print('MORPH_ROI',stage,roi.roi_id,len(ids),shape,flush=True)
    pd.DataFrame(summaries).to_csv(folder/'segmentation_geometry_comparison.tsv.gz',sep='\t',index=False)
    pd.DataFrame(seed_records).to_csv(folder/'segmentation_seed_audit.tsv.gz',sep='\t',index=False)
    (folder/'morphology_completion.json').write_text(json.dumps({'stage':stage,'ROIs':roi_meta,'membrane_source_sha256':digest(member),
        'DAPI_source_sha256':digest(dapi_path),'nuclear_boundaries_sha256':digest(CACHE/stage/'nucleus_boundaries.parquet'),
        'membrane_acquisition':json.loads((CACHE/stage/'morphology_focus_0001.ome.tif.acquisition.json').read_text()),
        'channels':channels,'used_channels':['DAPI','ATP1A1/CD45/E-Cadherin'],'native_pixel_size_um':px,
        'same_centroid_seed_count_in_all_variants':True,'nucleus_required':False,'no_RNA_or_lineage_input_to_watershed':True,
        'scope':'Fixed-seed morphology-informed boundary sensitivity; does not independently infer merges/splits or cell number',
        'boundary_provenance':provenance,'ROI_selection_sha256':digest(OUT/'selected_ROIs.tsv')},indent=2)+'\n')
    print('MORPHOLOGY_SEGMENTATIONS_COMPLETE',stage,flush=True)


if __name__=='__main__':main()
