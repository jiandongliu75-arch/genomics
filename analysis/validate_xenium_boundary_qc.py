"""Independent source, molecule-assignment, geometry and outcome checks."""
import gzip
import json

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely
from shapely.geometry import shape

from design_xenium_boundary_qc import OUT,ATLAS,FIT,STAGES
from reassess_whole_heart_lineages import SOURCE,WITHHOLD,digest
from process_xenium_boundary_molecules import VARIANTS
from analyze_xenium_boundary_depth_controls import CONTROLS
from acquire_xenium_boundary_inputs import CACHE


def main():
    rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t');checks=[];rng=np.random.default_rng(610)
    primary_old=pd.read_csv(FIT/'primary_local_Ankrd1_results.tsv',sep='\t').set_index('stage')
    for stage in STAGES:
        print('VALIDATE',stage,flush=True);folder=OUT/stage
        obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        source=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        completion=json.loads((folder/'molecule_completion.json').read_text());reconstruction=json.loads((folder/'molecule_count_reconstruction.json').read_text())
        assert reconstruction['mismatched_entries']==0 and reconstruction['reconstructed_QV20_integer_counts_exact']
        assert reconstruction['source_total_gene_transcripts']==int(source.X.sum())==reconstruction['reconstructed_total_gene_transcripts']
        assert digest(OUT/'boundary_design.json')==completion['design_sha256']
        assert digest(OUT/'selected_ROIs.tsv')==completion['ROI_selection_sha256']
        with np.load(folder/'geometry_sensitivity_counts.npz') as b:
            ids=b['cell_ids'].astype(str);genes=b['genes'].astype(str);counts={k:b[k] for k in VARIANTS}
        np.testing.assert_array_equal(genes,source.var_names)
        np.testing.assert_array_equal(counts['native_Q20'],source[ids].X.toarray())
        assert (counts['core_2um']<=counts['core_1um']).all() and (counts['core_1um']<=counts['native_Q20']).all()
        assert (counts['native_Q30']<=counts['native_Q20']).all() and (counts['dilation_2um']>=counts['native_Q20']).all()
        with np.load(folder/'depth_matched_control_counts.npz') as b:
            np.testing.assert_array_equal(b['cell_ids'].astype(str),ids)
            for variant,(matched,_) in CONTROLS.items():
                assert (b[variant]>=0).all() and (b[variant]<=counts['native_Q20']).all()
                np.testing.assert_array_equal(b[variant].sum(1),counts[matched].sum(1))
        original=pd.read_csv(folder/'native_Q20_projection_assignments.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        np.testing.assert_array_equal(original.index,ids)
        np.testing.assert_array_equal(original.new_identity,obs.loc[ids,'atlas_identity'])
        # All target profiles and same nucleus-independent rules exist in each variant.
        for variant in VARIANTS+list(CONTROLS):
            assignments=pd.read_csv(folder/(variant+'_projection_assignments.tsv.gz'),sep='\t',index_col=0,float_precision='round_trip')
            np.testing.assert_array_equal(assignments.index,ids)
            np.testing.assert_array_equal(assignments.nuclear_area_missing,obs.loc[ids,'nucleus_area'].isna())
        q=pd.read_csv(folder/'molecule_gene_totals.tsv',sep='\t')
        run=json.loads((folder/'verified_gene_QV_metrics.json').read_text())
        assert abs(q.QV20.sum()/q.all_QV.sum()-run['vendor_Q20_fraction'])<1e-12
        assert abs(q.QV20_assigned.sum()/q.QV20.sum()-run['vendor_gene_assigned_fraction'])<1e-12
        quality=pd.read_csv(folder/'transcript_quality_distribution.tsv',sep='\t')
        assert quality.all_features.sum()==completion['rows_streamed']
        sampled_distances=0;polygons_checked=0;ROI_molecules=0;seed_count=0;anucleate_CM_seeds=0
        for roi in rois[rois.stage.eq(stage)].itertuples():
            name=roi.roi_id;m=pq.read_table(folder/(name+'_molecules.parquet')).to_pandas();ROI_molecules+=len(m)
            assert m.transcript_id.is_unique and np.isfinite(m[['x_location','y_location','z_location','qv']]).all().all()
            assert len(m)==completion['ROI_rows'][name]
            with np.load(folder/(name+'_segmentation_maps.npz')) as b:
                cells=b['cell_ids'].astype(str);x0,y0,res=float(b['x0_um']),float(b['y0_um']),float(b['grid_um'])
                label_maps={k:b[k] for k in ['native_grid','centroid_Voronoi','membrane_watershed']}
            seed_count+=len(cells);anucleate_CM_seeds+=int((obs.loc[cells,'atlas_identity'].eq('Cardiomyocyte')&obs.loc[cells,'nucleus_area'].isna()).sum())
            for variant in ['centroid_Voronoi','membrane_watershed']:
                assert set(np.unique(label_maps[variant]))-set([0])==set(range(1,len(cells)+1))
                with gzip.open(folder/(name+'_'+variant+'_polygons.geojson.gz'),'rt') as h:geometries=json.load(h)['features']
                lookup={g['properties']['cell_id']:shape(g['geometry']) for g in geometries}
                assert set(lookup)==set(cells)
                # Pixel-center membership independently checks vectorization, including holes.
                yy=rng.integers(0,label_maps[variant].shape[0],1000);xx=rng.integers(0,label_maps[variant].shape[1],1000)
                labs=label_maps[variant][yy,xx];take=labs>0
                pts=shapely.points(np.column_stack([x0+xx[take]*res,y0+yy[take]*res]))
                own=np.array([lookup[cells[i-1]] for i in labs[take]],dtype=object)
                assert shapely.covered_by(pts,own).all()
                # Test both positive and negative pixel-center membership in cell bounding boxes.
                for label in rng.choice(np.arange(1,len(cells)+1),min(20,len(cells)),replace=False):
                    g=lookup[cells[label-1]];xmin,ymin,xmax,ymax=g.bounds
                    xa=np.clip(np.rint((np.array([xmin,xmax])-x0)/res).astype(int),0,label_maps[variant].shape[1]-1)
                    ya=np.clip(np.rint((np.array([ymin,ymax])-y0)/res).astype(int),0,label_maps[variant].shape[0]-1)
                    xx=rng.integers(xa[0],xa[1]+1,50);yy=rng.integers(ya[0],ya[1]+1,50)
                    pts=shapely.points(np.column_stack([x0+xx*res,y0+yy*res]))
                    np.testing.assert_array_equal(shapely.covered_by(pts,g),label_maps[variant][yy,xx]==label)
                polygons_checked+=len(lookup)
            gi=pd.Index(genes).get_indexer(m.feature_name);ci=pd.Index(cells).get_indexer(m.cell_id)
            valid=m.qv.ge(20).to_numpy()&(gi>=0)
            rr=np.rint((m.y_location.to_numpy()-y0)/res).astype(int);cc=np.rint((m.x_location.to_numpy()-x0)/res).astype(int)
            inside=(rr>=0)&(cc>=0)&(rr<label_maps['native_grid'].shape[0])&(cc<label_maps['native_grid'].shape[1])
            with np.load(folder/(name+'_ROI_count_variants.npz')) as b:
                recalculated=np.zeros_like(b['native_ROI_Q20']);select=valid&(ci>=0)
                np.add.at(recalculated,(ci[select],gi[select]),1);np.testing.assert_array_equal(recalculated,b['native_ROI_Q20'])
                for variant,map_name in [('native_raster_Q20','native_grid'),('centroid_Voronoi','centroid_Voronoi'),('membrane_watershed','membrane_watershed')]:
                    assigned=np.full(len(m),-1,int);assigned[inside]=label_maps[map_name][rr[inside],cc[inside]]-1
                    keep=valid&(assigned>=0);rebuilt=np.zeros_like(recalculated)
                    np.add.at(rebuilt,(assigned[keep],gi[keep]),1);np.testing.assert_array_equal(rebuilt,b[variant])
                    assert rebuilt.sum()<=valid.sum()
            with gzip.open(folder/(name+'_native_polygons.geojson.gz'),'rt') as h:native={f['properties']['cell_id']:shape(f['geometry']) for f in json.load(h)['features']}
            known=m.cell_id.isin(native)&m.qv.ge(20)&m.feature_name.isin(genes)&m.native_boundary_distance_um_Q20_target.notna()
            select=np.flatnonzero(known);sample=rng.choice(select,min(100,len(select)),replace=False)
            p=shapely.points(m[['x_location','y_location']].iloc[sample].to_numpy());g=np.array([native[c] for c in m.cell_id.iloc[sample]],dtype=object)
            np.testing.assert_allclose(shapely.distance(p,shapely.boundary(g)),m.native_boundary_distance_um_Q20_target.iloc[sample],rtol=1e-6,atol=2e-6)
            np.testing.assert_array_equal(shapely.covered_by(p,g),m.native_inside_polygon_Q20_target.iloc[sample])
            sampled_distances+=len(sample)
        # Recompute every primary fixed-pair mean/ratio from the published alternative counts.
        pairs=pd.read_csv(FIT/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
        pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
        ix=pd.Index(ids);near=ix.get_indexer(pairs.near_cell_id);far=ix.get_indexer(pairs.far_cell_id);assert (near>=0).all() and (far>=0).all()
        table=pd.read_csv(folder/'Ankrd1_boundary_sensitivity.tsv',sep='\t')
        table=table[table.scope.eq('geometry_targets_and_all_frozen_pairs')&table.pair_filter.eq('all_frozen_pairs')&table.nuclear_area_stratum.eq('all')].set_index('variant')
        for variant,values in counts.items():
            denominator=values.sum(1)-values[:,np.isin(genes,WITHHOLD)].sum(1)
            expression=1e4*values[:,list(genes).index('Ankrd1')]/np.maximum(denominator,1)
            ratio=expression[near].mean()/expression[far].mean();difference=np.log1p(expression[near])-np.log1p(expression[far])
            equal=pd.Series(difference).groupby(pairs.segment_id.to_numpy()).mean().mean()
            assert table.loc[variant,'pairs']==len(pairs)
            assert abs(table.loc[variant,'near_far_mean_Ankrd1_CP10k_ratio']-ratio)<1e-12
            assert abs(table.loc[variant,'equal_segment_mean_log_difference']-equal)<1e-12
        assert abs(table.loc['native_Q20','near_far_mean_Ankrd1_CP10k_ratio']-primary_old.loc[stage,'near_far_Ankrd1_CP10k_ratio'])<1e-12
        checks.append({'stage':stage,'source_profiles_with_exact_Q20_reconstruction':source.n_obs,'whole_run_molecule_rows':completion['rows_streamed'],
            'target_profiles_with_exact_native_identity_reproduction':len(ids),'native_and_boundary_variant_count_checks':True,
            'ROI_molecule_records_checked_including_context_overlap':ROI_molecules,'alternative_cell_polygons_checked':polygons_checked,
            'source_boundary_distances_independently_recomputed':sampled_distances,'cell_seeds_in_context_ROIs':seed_count,
            'nuclear_area_missing_CM_seeds_preserved':anucleate_CM_seeds,'all_frozen_CM_pairs_checked':len(pairs)})
        print('BOUNDARY_VALIDATION_PASS',stage,checks[-1],flush=True)
    (OUT/'validation.json').write_text(json.dumps({'status':'PASS','stages':checks,
        'scope':'Source count reconstruction, native identity reproduction, count-depth controls, RNA assignment, polygon vectorization and paired outcomes; not independently adjudicated biological cell boundaries or animal replication'},indent=2)+'\n')


if __name__=='__main__':main()
