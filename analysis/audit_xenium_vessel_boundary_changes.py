"""Selected arterial envelopes rebuilt from alternative cell boundaries."""
import gzip
import json

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape,mapping

from design_xenium_boundary_qc import OUT,ATLAS,FIT,STAGES
from reassess_whole_heart_lineages import WITHHOLD


def main():
    rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t');rows=[];pair_rows=[];features=[]
    for stage in STAGES:
        folder=OUT/stage
        assert (folder/'ROI_reassignment_completion.json').exists()
        obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        original=json.loads((FIT/stage/'arterial_segments.geojson').read_text())['features']
        original={f['properties']['segment_id']:f for f in original}
        pairs=pd.read_csv(FIT/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
        pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
        for roi in rois[rois.stage.eq(stage)&rois.category.eq('Peri-arterial')].itertuples():
            f=original[roi.segment_id];seed=f['properties']['arterial_EC_ids']+f['properties']['SMC_ids'];old=shape(f['geometry'])
            with np.load(folder/(roi.roi_id+'_ROI_count_variants.npz')) as b:
                ids=b['cell_ids'].astype(str);genes=b['genes'].astype(str);eligible=b['comparison_eligible']
                matrices={k:b[k] for k in ['native_ROI_Q20','native_raster_Q20','centroid_Voronoi','membrane_watershed']}
            index=pd.Index(ids);local=pairs[pairs.segment_id.eq(roi.segment_id)].copy()
            ni=index.get_indexer(local.near_cell_id);fi=index.get_indexer(local.far_cell_id);valid=(ni>=0)&(fi>=0)
            k=valid.copy();valid[k]&=eligible[ni[k]]&eligible[fi[k]];local=local[valid];ni=ni[valid];fi=fi[valid]
            xy_near=shapely.points(obs.loc[local.near_cell_id,['x_centroid','y_centroid']].to_numpy())
            xy_far=shapely.points(obs.loc[local.far_cell_id,['x_centroid','y_centroid']].to_numpy())
            for variant in ['native_ROI_Q20','native_raster_Q20','centroid_Voronoi','membrane_watershed']:
                if variant in ['native_ROI_Q20','native_raster_Q20']:envelope=old
                else:
                    with gzip.open(folder/(roi.roi_id+'_'+variant+'_polygons.geojson.gz'),'rt') as h:polygons=json.load(h)['features']
                    lookup={g['properties']['cell_id']:shape(g['geometry']) for g in polygons}
                    assert set(seed).issubset(lookup)
                    envelope=shapely.union_all([lookup[c] for c in seed]).buffer(3)
                dn=shapely.distance(xy_near,envelope);df=shapely.distance(xy_far,envelope)
                same=(dn>0)&(dn<=30)&(df>=60)&(df<=150)
                counts=matrices[variant];den=counts.sum(1)-counts[:,np.isin(genes,WITHHOLD)].sum(1)
                value=1e4*counts[:,list(genes).index('Ankrd1')]/np.maximum(den,1)
                for subset,choose in [('same_original_pairs',np.ones(len(local),bool)),('distance_categories_preserved',same)]:
                    near,far=ni[choose],fi[choose]
                    ratio=float(value[near].mean()/max(value[far].mean(),1e-12)) if len(near) else np.nan
                    rows.append({'stage':stage,'roi_id':roi.roi_id,'segment_id':roi.segment_id,'variant':variant,'pair_subset':subset,
                        'pairs':int(choose.sum()),'Ankrd1_near_far_ratio':ratio,
                        'mean_paired_log_difference':float((np.log1p(value[near])-np.log1p(value[far])).mean()) if len(near) else np.nan,
                        'envelope_IoU':float(envelope.intersection(old).area/envelope.union(old).area),
                        'fraction_original_distance_pair_categories_preserved':float(same.mean()) if len(same) else np.nan,
                        'median_near_distance_um':float(np.median(dn)) if len(dn) else np.nan,
                        'median_far_distance_um':float(np.median(df)) if len(df) else np.nan})
                for i,row in enumerate(local.itertuples()):
                    pair_rows.append({'stage':stage,'roi_id':roi.roi_id,'variant':variant,'near_cell_id':row.near_cell_id,'far_cell_id':row.far_cell_id,
                        'near_distance_um':float(dn[i]),'far_distance_um':float(df[i]),'distance_categories_preserved':bool(same[i])})
                features.append({'type':'Feature','properties':{'stage':stage,'roi_id':roi.roi_id,'segment_id':roi.segment_id,'variant':variant},'geometry':mapping(envelope)})
    pd.DataFrame(rows).to_csv(OUT/'selected_vessel_envelope_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(pair_rows).to_csv(OUT/'selected_vessel_pair_distances.tsv.gz',sep='\t',index=False,na_rep='NA')
    with gzip.open(OUT/'selected_alternative_vessel_envelopes.geojson.gz','wt') as h:json.dump({'type':'FeatureCollection','features':features},h)
    print('SELECTED_VESSEL_BOUNDARY_SENSITIVITY_COMPLETE',len(rows),flush=True)


if __name__=='__main__':main()
