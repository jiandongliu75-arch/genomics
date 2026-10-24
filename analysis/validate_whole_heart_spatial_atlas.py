"""Source conservation and independent marker/spatial checks for the atlas."""
import gzip
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from shapely.geometry import shape,Point

from build_whole_heart_spatial_atlas import OUT,FIT,SOURCE,STAGES,LINEAGES,REVIEW,ANCHORS,PANELS,DISPLAY_GENES,digest


def main():
    counts=pd.read_csv(OUT/'lineage_counts.tsv',sep='\t')
    profiles=pd.read_csv(OUT/'measured_gene_profiles.tsv',sep='\t')
    neighbor=pd.read_csv(OUT/'neighborhood_composition.tsv',sep='\t')
    tiles=pd.read_csv(OUT/'spatial_100um_tile_counts.tsv.gz',sep='\t')
    audit=pd.read_csv(OUT/'nonCM_identity_audit.tsv',sep='\t')
    provenance=json.loads((OUT/'source_provenance.json').read_text())
    checks=[];allids=[]
    for stage in STAGES:
        folder=OUT/stage;print('VALIDATE',stage,flush=True)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        atlas=ad.read_h5ad(folder/f'{stage}_whole_heart_atlas.h5ad')
        obs=pd.read_csv(folder/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        previous=pd.read_csv(FIT/stage/'whole_heart_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        assert data.shape==atlas.shape and data.n_vars==403
        np.testing.assert_array_equal(atlas.obs_names,data.obs_names);np.testing.assert_array_equal(obs.index,data.obs_names)
        np.testing.assert_array_equal(data.var_names,atlas.var_names)
        assert (data.X!=atlas.X).nnz==0
        np.testing.assert_array_equal(atlas.obsm['spatial'],data.obs[['x_centroid','y_centroid']])
        np.testing.assert_array_equal(obs[['x_centroid','y_centroid']],data.obs[['x_centroid','y_centroid']])
        np.testing.assert_array_equal(obs.atlas_identity.astype(str),atlas.obs.atlas_identity.astype(str))
        assert set(obs.atlas_identity).issubset(set(LINEAGES+REVIEW))
        assert obs.stage.eq(stage).all() and obs.stage_cell_id.is_unique
        np.testing.assert_array_equal(obs.stage_cell_id,stage+'::'+obs.index)
        allids.extend(obs.stage_cell_id)
        assert obs.atlas_resolved.eq(obs.atlas_identity.isin(LINEAGES)).all()
        expected=previous.reassessed_lineage.replace({'Mixed CM/non-CM; review':'Mixed-marker review'}).copy()
        total=np.asarray(data.X.sum(axis=1)).ravel();panel_counts={}
        for label,genes in PANELS.items():
            subset=[g for g in genes if g in data.var_names];x=data[:,subset].X
            molecules=np.asarray(x.sum(axis=1)).ravel();detected=np.asarray((x>0).sum(axis=1)).ravel()
            panel_counts[label]=(molecules,detected)
            key=label.replace(' ','_')
            np.testing.assert_array_equal(obs['marker_'+key+'_molecules'],molecules)
            np.testing.assert_array_equal(obs['marker_'+key+'_genes'],detected)
            np.testing.assert_allclose(obs['marker_'+key+'_fraction'],molecules/np.maximum(total,1),atol=1e-12)
        cm=panel_counts['Cardiomyocyte'];cmstrong=(cm[1]>=3)&(cm[0]/np.maximum(total,1)>=.1)
        mixed=np.zeros(len(obs),bool)
        for label in ['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune']:
            mol,det=panel_counts[label];mixed|=(det>=3)&(mol/np.maximum(total,1)>=.1)
        mixed&=cmstrong
        np.testing.assert_array_equal(mixed,previous.mixed_CM_nonCM_marker_review)
        for label,genes in ANCHORS.items():
            subset=[g for g in genes if g in data.var_names];det=np.asarray((data[:,subset].X>0).sum(axis=1)).ravel()
            np.testing.assert_array_equal(obs['anchor_'+label.replace(' ','_')+'_genes'],det)
            take=previous.reassessed_lineage.eq(label).to_numpy()
            expected.iloc[np.flatnonzero(take&(det==0))]='Ambiguous/review'
            expected.iloc[np.flatnonzero(take&mixed)]='Mixed-marker review'
            row=audit[audit.stage.eq(stage)&audit.lineage.eq(label)].iloc[0]
            assert row.previous_resolved_cells==take.sum() and row.no_anchor_cells==(take&(det==0)).sum()
            assert row.mixed_CM_nonCM_cells==(take&mixed).sum()
        np.testing.assert_array_equal(expected,obs.atlas_identity)
        assert obs.atlas_identity.eq('Cardiomyocyte').equals(previous.CM_strict)
        resolved=obs.atlas_resolved
        assert obs.loc[resolved,'eligible_whole_reference'].all()
        assert obs.loc[resolved,'lineage_weight'].ge(.7).all() and obs.loc[resolved,'lineage_margin'].ge(.15).all()
        assert obs.loc[resolved,'both_RNA_libraries_agree'].eq(True).all() and obs.loc[resolved,'raw_reference_agrees'].eq(True).all()
        for label in LINEAGES+REVIEW:
            take=obs.atlas_identity.eq(label).to_numpy();n=take.sum()
            row=counts[counts.stage.eq(stage)&counts.atlas_identity.eq(label)].iloc[0];assert row.cells==n
            if n:
                part=data.X[take];mean=np.asarray(part.mean(axis=0)).ravel();detect=np.asarray((part>0).mean(axis=0)).ravel()
                summary=profiles[profiles.stage.eq(stage)&profiles.atlas_identity.eq(label)].set_index('gene').loc[data.var_names]
                np.testing.assert_allclose(summary.mean_raw_count,mean,rtol=1e-5,atol=1e-6)
                np.testing.assert_allclose(summary.detection_fraction,detect,rtol=1e-5,atol=1e-6)
                gi=data.var_names.get_indexer(DISPLAY_GENES);small=part[:,gi].toarray()
                transformed=np.log1p(1e4*small/np.maximum(total[take,None],1)).mean(axis=0)
                np.testing.assert_allclose(summary.loc[DISPLAY_GENES,'mean_log1pCP10k'],transformed,atol=1e-10)
        # Rebuild exact neighborhood summaries, then independently brute-force samples.
        xy=obs[['x_centroid','y_centroid']].to_numpy();code=pd.Index(LINEAGES).get_indexer(obs.atlas_identity)
        distance,ix=cKDTree(xy).query(xy,k=21,workers=2)
        valid=(ix!=np.arange(len(obs))[:,None])&(distance<=100)
        assert (valid.sum(1)<=20).all()
        edge=np.column_stack([((code[ix]==i)&valid).sum(1) for i in range(7)])
        np.testing.assert_array_equal(edge.sum(1),obs.resolved_neighbors_within_100um_among_20_nearest)
        assert obs.neighborhood_eligible.eq(resolved&(edge.sum(1)>=5)).all()
        for i,label in enumerate(LINEAGES):
            take=obs.neighborhood_eligible&obs.atlas_identity.eq(label)
            saved=neighbor[neighbor.stage.eq(stage)&neighbor.focal_lineage.eq(label)].set_index('neighbor_lineage').loc[LINEAGES]
            np.testing.assert_array_equal(saved.links_to_neighbor_lineage,edge[take].sum(0))
            np.testing.assert_allclose(saved.neighbor_fraction,edge[take].sum(0)/max(edge[take].sum(),1),atol=1e-12)
        sample=np.random.default_rng(917).choice(len(obs),50,replace=False)
        brute=cdist(xy[sample],xy);brute[np.arange(len(sample)),sample]=np.inf
        nearest=np.argsort(brute,axis=1,kind='stable')[:,:20]
        chosen=np.take_along_axis(brute,nearest,axis=1)<=100
        reference=np.column_stack([((code[nearest]==i)&chosen).sum(1) for i in range(7)])
        np.testing.assert_array_equal(reference,edge[sample])
        # Native-position tile totals independently reconstructed.
        t=tiles[tiles.stage.eq(stage)].set_index(['tile_x','tile_y'])
        positions=pd.MultiIndex.from_arrays([np.floor(xy[:,0]/100).astype(int),np.floor(xy[:,1]/100).astype(int)])
        grouped=pd.Series(1,index=positions).groupby(level=[0,1]).sum()
        np.testing.assert_array_equal(t.loc[grouped.index,'original_cells'],grouped)
        assert t.original_cells.sum()==len(obs) and t[LINEAGES].sum().sum()==resolved.sum()
        # Exported zoom geometry is the deposited segmentation, not inferred boundaries.
        zoom=json.loads((folder/'geometry_zoom_provenance.json').read_text());bounds=zoom['bounds_um']
        z=pd.read_csv(folder/'geometry_zoom_cells.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        selected=obs.x_centroid.between(bounds[0],bounds[1])&obs.y_centroid.between(bounds[2],bounds[3])
        np.testing.assert_array_equal(z.index,obs.index[selected])
        with gzip.open(folder/'geometry_zoom_polygons.geojson.gz','rt') as handle:g=json.load(handle)
        assert set(f['properties']['cell_id'] for f in g['features'])==set(z.index)
        for feature in g['features']:
            cell=feature['properties']['cell_id'];poly=shape(feature['geometry'])
            assert poly.is_valid and poly.area>0 and feature['properties']['atlas_identity']==obs.at[cell,'atlas_identity']
            assert poly.distance(Point(*obs.loc[cell,['x_centroid','y_centroid']]))<=zoom['maximum_centroid_outside_own_polygon_um']+1e-8
        record=next(p for p in provenance if p['stage']==stage)
        assert digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')==record['source_sha256']
        assert digest(FIT/stage/'whole_heart_cell_annotations.tsv.gz')==record['previous_annotations_sha256']
        checks.append({'stage':stage,'original_cells_checked':len(obs),'count_entries_checked':int(data.X.nnz),
            'resolved_cells':int(resolved.sum()),'native_zoom_polygons_checked':len(g['features']),
            'independent_brute_force_neighborhood_cells':len(sample),'count_and_coordinate_conservation':True,
            'nonCM_anchor_and_mixed_rules_independently_recomputed':True})
        print('PASS',stage,checks[-1],flush=True)
    assert len(allids)==408054 and len(set(allids))==len(allids)
    (OUT/'validation.json').write_text(json.dumps({'status':'PASS','stages':checks,
        'scope':'Source count/coordinate conservation, marker rules, exported identities and spatial calculations; not independent anatomical or animal-level validation'},indent=2)+'\n')


if __name__=='__main__':main()
