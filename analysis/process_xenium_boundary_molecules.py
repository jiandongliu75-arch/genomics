"""Stream original molecules, verify Q20 counts, audit boundary-localized RNA."""
import argparse
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import shapely
from shapely.strtree import STRtree

from design_xenium_boundary_qc import OUT,ATLAS,FIT,STAGES
from acquire_xenium_boundary_inputs import CACHE
from reassess_whole_heart_lineages import SOURCE,PANELS,digest
from define_heart_arterial_structures import read_boundaries

VARIANTS=['native_Q20','core_1um','core_2um','dilation_2um','native_Q30']


def nearest_distance(tree,points,maximum=5):
    d=np.full(len(points),np.inf)
    if len(tree.geometries) and len(points):
        pairs,value=tree.query_nearest(points,max_distance=maximum,all_matches=False,return_distance=True)
        d[pairs[0]]=value
    return d


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=STAGES,required=True);args=parser.parse_args();stage=args.stage
    folder=OUT/stage;folder.mkdir(exist_ok=True);start=time.time()
    if (folder/'molecule_completion.json').exists():print('SKIP_MOLECULES',stage,flush=True);return
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    source=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
    genes=source.var_names.to_numpy(dtype=str);gene_index=pd.Index(genes);ng=len(genes)
    rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t');rois=rois[rois.stage.eq(stage)]
    pairs=pd.read_csv(FIT/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
    pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
    keep=obs.index.isin(pd.unique(pairs[['near_cell_id','far_cell_id']].to_numpy().ravel()))
    in_core=np.zeros(len(obs),bool)
    for roi in rois.itertuples():
        in_core|=(obs.x_centroid.between(roi.x_min,roi.x_max)&obs.y_centroid.between(roi.y_min,roi.y_max)).to_numpy()
        keep|=(obs.x_centroid.between(roi.x_min-40,roi.x_max+40)&obs.y_centroid.between(roi.y_min-40,roi.y_max+40)).to_numpy()
    target_ids=obs.index[keep];target=obs.loc[target_ids].copy();target['ROI_inner_evaluation_cell']=in_core[keep]
    target['audit_lineage']=target.atlas_identity
    unresolved=~target.atlas_resolved
    target.loc[unresolved,'audit_lineage']=target.loc[unresolved,'previous_reassessed_lineage']
    target['nuclear_area_missing']=target.nucleus_area.isna()
    target.index.name='cell_id';target.to_csv(folder/'boundary_target_cells.tsv.gz',sep='\t',na_rep='NA')
    mapping=np.full(len(obs),-1,np.int32);mapping[np.flatnonzero(keep)]=np.arange(keep.sum())
    polygons,provenance=read_boundaries(stage,set(obs.index),obs)
    geom=np.array([polygons[c] for c in obs.index],dtype=object);bound=shapely.boundary(geom)
    tree=STRtree(geom)
    cm_tree=STRtree(geom[obs.atlas_identity.eq('Cardiomyocyte').to_numpy()])
    non_tree=STRtree(geom[obs.atlas_identity.isin(['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']).to_numpy()])
    counts={name:np.zeros((len(target),ng),np.int32) for name in VARIANTS}
    region_names=['deep_core','opposite_facing_shell','away_shell']
    region_counts=np.zeros((3,len(target),ng),np.int32)
    native=np.zeros(source.shape,np.int32)
    stats=np.zeros((len(target),7),np.int64)
    stat_names=['Q20_assigned_gene_transcripts','outside_native_polygon','inside_shell1','inside_core1','inside_core2','Q20_nucleus_overlap','dilation_added']
    global_counts={'rows':0,'Q20_all_features':0,'Q20_gene_transcripts':0,'Q20_assigned_genes':0,
        'Q20_unassigned_genes':0,'Q30_gene_transcripts':0,'unassigned_dilation_ties_rejected':0,'unknown_assigned_cell_ids':0}
    quality=np.zeros((6,3),np.int64);qedges=np.array([-np.inf,10,20,30,40,50,np.inf])
    gene_totals=np.zeros((ng,4),np.int64)
    writers={};roi_rows={r.roi_id:0 for r in rois.itertuples()}
    path=CACHE/stage/'transcripts.parquet';parquet=pq.ParquetFile(path)
    print('MOLECULE_SCHEMA',stage,str(parquet.schema),flush=True)
    expected_columns=['transcript_id','cell_id','overlaps_nucleus','feature_name','x_location','y_location','z_location','qv']
    assert set(expected_columns).issubset(parquet.schema_arrow.names)
    print('STREAM',stage,'all_rows',parquet.metadata.num_rows,'targets',len(target),'ROIs',len(rois),flush=True)
    try:
        for batch_number,batch in enumerate(parquet.iter_batches(batch_size=200000,columns=expected_columns)):
            table=pa.Table.from_batches([batch]);frame=table.to_pandas()
            feature=frame.feature_name.map(lambda x:x.decode() if isinstance(x,bytes) else x)
            names=frame.cell_id.map(lambda x:x.decode() if isinstance(x,bytes) else x)
            ci=obs.index.get_indexer(names);gi=gene_index.get_indexer(feature)
            ti=np.full(len(frame),-1,np.int32);assigned=ci>=0;ti[assigned]=mapping[ci[assigned]]
            gene=gi>=0;q=frame.qv.to_numpy();high=q>=20;strict=q>=30
            xy=frame[['x_location','y_location']].to_numpy();overlap=frame.overlaps_nucleus.fillna(False).astype(bool).to_numpy()
            hgene=high&gene;original=hgene&assigned;selected=original&(ti>=0)
            global_counts['rows']+=len(frame);global_counts['Q20_all_features']+=int(high.sum())
            global_counts['Q20_gene_transcripts']+=int(hgene.sum());global_counts['Q20_assigned_genes']+=int(original.sum())
            global_counts['Q20_unassigned_genes']+=int((hgene&~assigned).sum());global_counts['Q30_gene_transcripts']+=int((strict&gene).sum())
            unknown=~assigned&~names.isin(['UNASSIGNED','Unassigned','unassigned','-1']).to_numpy()
            global_counts['unknown_assigned_cell_ids']+=int(unknown.sum())
            bins=np.clip(np.digitize(q,qedges)-1,0,5)
            for j,mask in enumerate([np.ones(len(frame),bool),gene,gene&assigned]):quality[:,j]+=np.bincount(bins[mask],minlength=6)
            for j,mask in enumerate([gene,hgene,original,strict&gene]):np.add.at(gene_totals[:,j],gi[mask],1)
            np.add.at(native,(ci[original],gi[original]),1)
            for name,mask in [('native_Q20',selected),('native_Q30',strict&gene&assigned&(ti>=0))]:
                np.add.at(counts[name],(ti[mask],gi[mask]),1)
            inside=np.zeros(len(frame),bool);distance=np.full(len(frame),np.nan,np.float32)
            selected_indices=np.flatnonzero(selected)
            if len(selected_indices):
                points=shapely.points(xy[selected]);sg=geom[ci[selected]]
                contained=shapely.covered_by(points,sg);dd=shapely.distance(points,bound[ci[selected]])
                inside[selected]=contained;distance[selected]=dd
                for name,threshold in [('core_1um',1),('core_2um',2)]:
                    ok=contained&(dd>threshold)
                    np.add.at(counts[name],(ti[selected][ok],gi[selected][ok]),1)
                events=[np.ones(len(points),bool),~contained,contained&(dd<=1),contained&(dd>1),contained&(dd>2),overlap[selected]]
                for column,event in enumerate(events):np.add.at(stats[:,column],ti[selected][event],1)
                # Spatial shell direction uses only resolved opposite-lineage polygons.
                cm_mask=target.audit_lineage.eq('Cardiomyocyte').to_numpy()[ti[selected]]
                opposite=np.full(len(points),np.inf)
                opposite[~cm_mask]=nearest_distance(cm_tree,points[~cm_mask])
                opposite[cm_mask]=nearest_distance(non_tree,points[cm_mask])
                regions=[contained&(dd>2),contained&(dd<=1)&(opposite<=2),contained&(dd<=1)&(opposite>5)]
                for j,mask in enumerate(regions):np.add.at(region_counts[j],(ti[selected][mask],gi[selected][mask]),1)
            # Limited expansion: only currently unassigned high-quality gene molecules.
            unassigned=hgene&~assigned
            if unassigned.any():
                points=shapely.points(xy[unassigned])
                indices,dd=tree.query_nearest(points,max_distance=2,all_matches=True,return_distance=True)
                if indices.size:
                    query,cell=indices;hits=np.bincount(query,minlength=len(points));unique=hits[query]==1
                    global_counts['unassigned_dilation_ties_rejected']+=int((hits>1).sum())
                    target_row=mapping[cell];accept=unique&(target_row>=0)
                    source_rows=np.flatnonzero(unassigned)[query[accept]]
                    np.add.at(counts['dilation_2um'],(target_row[accept],gi[source_rows]),1)
                    np.add.at(stats[:,6],target_row[accept],1)
            for roi in rois.itertuples():
                mask=(xy[:,0]>=max(0,roi.x_min-40))&(xy[:,0]<=roi.x_max+40)&(xy[:,1]>=max(0,roi.y_min-40))&(xy[:,1]<=roi.y_max+40)
                if not mask.any():continue
                piece=table.filter(pa.array(mask))
                piece=piece.append_column('source_global_cell_index',pa.array(ci[mask],type=pa.int32()))
                piece=piece.append_column('native_inside_polygon_Q20_target',pa.array(inside[mask]))
                piece=piece.append_column('native_boundary_distance_um_Q20_target',pa.array(distance[mask]))
                if roi.roi_id not in writers:
                    dest=folder/(roi.roi_id+'_molecules.parquet');writers[roi.roi_id]=pq.ParquetWriter(dest,piece.schema,compression='zstd')
                writers[roi.roi_id].write_table(piece);roi_rows[roi.roi_id]+=len(piece)
            if batch_number%20==0:print('BATCH',stage,batch_number,global_counts['rows'],flush=True)
    finally:
        for writer in writers.values():writer.close()
    expected=source.X.toarray()
    discrepancy=native.astype(np.int64)-expected
    mismatch=int(np.count_nonzero(discrepancy));maximum=int(np.abs(discrepancy).max())
    recon={'stage':stage,'source_gene_count_entries':int(source.X.nnz),'reconstructed_QV20_integer_counts_exact':mismatch==0,
        'mismatched_entries':mismatch,'maximum_absolute_difference':maximum,'source_total_gene_transcripts':int(expected.sum()),
        'reconstructed_total_gene_transcripts':int(native.sum()),'QV_threshold':20}
    (folder/'molecule_count_reconstruction.json').write_text(json.dumps(recon,indent=2)+'\n')
    assert mismatch==0,recon
    counts['dilation_2um']+=counts['native_Q20']
    np.savez_compressed(folder/'geometry_sensitivity_counts.npz',cell_ids=target_ids.to_numpy(dtype=str),genes=genes,
        **counts,region_counts=region_counts,region_names=np.array(region_names))
    stat=pd.DataFrame(stats,index=target_ids,columns=stat_names);stat.index.name='cell_id'
    stat=target[['atlas_identity','audit_lineage','nuclear_area_missing','ROI_inner_evaluation_cell']].join(stat)
    stat.to_csv(folder/'cell_boundary_audit.tsv.gz',sep='\t',na_rep='NA')
    quality_table=pd.DataFrame(quality,columns=['all_features','gene_targets','assigned_gene_targets'])
    quality_table['QV_bin']=['<10','10–20','20–30','30–40','40–50','>=50'];quality_table.to_csv(folder/'transcript_quality_distribution.tsv',sep='\t',index=False)
    pd.DataFrame(gene_totals,index=genes,columns=['all_QV','QV20','QV20_assigned','QV30']).rename_axis('gene').to_csv(folder/'molecule_gene_totals.tsv',sep='\t')
    run=pd.read_csv(CACHE/stage/'metrics_summary.csv').iloc[0]
    global_counts['observed_Q20_fraction']=global_counts['Q20_all_features']/global_counts['rows']
    global_counts['observed_Q20_gene_assigned_fraction']=global_counts['Q20_assigned_genes']/global_counts['Q20_gene_transcripts']
    global_counts['vendor_Q20_fraction']=float(run.fraction_transcripts_decoded_q20)
    global_counts['vendor_gene_assigned_fraction']=float(run.fraction_transcripts_assigned)
    (folder/'molecule_run_metrics.json').write_text(json.dumps(global_counts,indent=2)+'\n')
    record={'stage':stage,'rows_streamed':parquet.metadata.num_rows,'target_cells':len(target),'ROI_rows':roi_rows,
        'source_transcripts_sha256':digest(path),'source_native_count_matrix_sha256':digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad'),
        'source_transcript_acquisition':json.loads((CACHE/stage/'transcripts.parquet.acquisition.json').read_text()),
        'native_boundary_provenance':provenance,'count_reconstruction':recon,'nucleus_required':False,
        'design_sha256':digest(OUT/'boundary_design.json'),'ROI_selection_sha256':digest(OUT/'selected_ROIs.tsv'),
        'elapsed_seconds':time.time()-start}
    (folder/'molecule_completion.json').write_text(json.dumps(record,indent=2)+'\n')
    print('MOLECULE_AUDIT_COMPLETE',stage,global_counts,flush=True)


if __name__=='__main__':main()
