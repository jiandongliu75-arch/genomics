"""Rebuild ROI expression from actual Q20 molecules under boundary alternatives."""
import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd

from design_xenium_boundary_qc import OUT,ATLAS,STAGES
from reassess_whole_heart_lineages import SOURCE,PANELS,digest


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=STAGES,required=True);args=p.parse_args();stage=args.stage
    folder=OUT/stage
    if (folder/'ROI_reassignment_completion.json').exists():print('SKIP_ROI_REASSIGNMENT',stage,flush=True);return
    assert (folder/'molecule_completion.json').exists() and (folder/'morphology_completion.json').exists()
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    source=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
    genes=source.var_names.to_numpy(dtype=str);rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t');rois=rois[rois.stage.eq(stage)]
    summaries=[];transitions=[];coverage=[]
    grouping=np.array(['Other_gene']*len(genes),dtype=object)
    for name in ['Cardiomyocyte','Endothelial','Fibroblast']:grouping[np.isin(genes,PANELS[name])]=name+'_panel'
    grouping[genes=='Ankrd1']='Ankrd1'
    for roi in rois.itertuples():
        name=roi.roi_id;molecules=pd.read_parquet(folder/(name+'_molecules.parquet'))
        with np.load(folder/(name+'_segmentation_maps.npz')) as b:
            ids=b['cell_ids'].astype(str);nativegrid=b['native_grid'];voronoi=b['centroid_Voronoi'];membrane=b['membrane_watershed']
            x0,y0,res=float(b['x0_um']),float(b['y0_um']),float(b['grid_um'])
        ci=pd.Index(ids).get_indexer(molecules.cell_id);gi=pd.Index(genes).get_indexer(molecules.feature_name)
        valid=molecules.qv.ge(20).to_numpy()&(gi>=0)
        r=np.rint((molecules.y_location.to_numpy()-y0)/res).astype(int);c=np.rint((molecules.x_location.to_numpy()-x0)/res).astype(int)
        inside=(r>=0)&(c>=0)&(r<membrane.shape[0])&(c<membrane.shape[1])
        matrices={'native_ROI_Q20':np.zeros((len(ids),len(genes)),np.int32)}
        take=valid&(ci>=0);np.add.at(matrices['native_ROI_Q20'],(ci[take],gi[take]),1)
        full=source.X[obs.index.get_indexer(ids)].toarray();fraction=matrices['native_ROI_Q20'].sum(1)/np.maximum(full.sum(1),1)
        cells=obs.loc[ids];inner=cells.x_centroid.between(roi.x_min,roi.x_max)&cells.y_centroid.between(roi.y_min,roi.y_max)
        touches=np.zeros((3,len(ids)),bool)
        original_global=molecules.source_global_cell_index.to_numpy();old=np.array(['Unassigned']*len(molecules),dtype=object)
        known=original_global>=0;old[known]=obs.atlas_identity.to_numpy()[original_global[known]]
        for i,(variant,labels) in enumerate([('native_raster_Q20',nativegrid),('centroid_Voronoi',voronoi),('membrane_watershed',membrane)]):
            assigned=np.full(len(molecules),-1,np.int32);assigned[inside]=labels[r[inside],c[inside]]-1
            counts=np.zeros((len(ids),len(genes)),np.int32);keep=valid&(assigned>=0)
            np.add.at(counts,(assigned[keep],gi[keep]),1);matrices[variant]=counts
            edge=np.unique(np.concatenate([labels[0],labels[-1],labels[:,0],labels[:,-1]]));edge=edge[edge>0]-1;touches[i,edge]=True
            new=np.array(['Unassigned']*len(molecules),dtype=object);known_new=assigned>=0
            new[known_new]=cells.atlas_identity.to_numpy()[assigned[known_new]]
            moved=valid&(ci>=0)&(assigned!=ci)
            summaries.append({'stage':stage,'roi_id':name,'variant':variant,'Q20_gene_molecules':int(valid.sum()),
                'assigned_Q20_genes':int(keep.sum()),'originally_ROI_assigned_Q20_genes':int((valid&(ci>=0)).sum()),
                'originally_ROI_assigned_changed_cell_or_unassigned':int(moved.sum()),
                'previously_unassigned_Q20_genes_now_assigned':int((valid&(original_global<0)&known_new).sum()),
                'outside_context_Q20_genes':int((valid&~inside).sum()),'total_seed_cells':len(ids)})
            frame=pd.DataFrame({'original_identity':old[valid],'new_identity':new[valid],'gene_group':grouping[gi[valid]]})
            for keys,n in frame.groupby(['original_identity','new_identity','gene_group']).size().items():
                transitions.append({'stage':stage,'roi_id':name,'variant':variant,'original_identity':keys[0],'new_identity':keys[1],
                    'gene_group':keys[2],'transcripts':int(n)})
        reliable=inner.to_numpy()&(fraction>=.95)&~touches[1:].any(0)
        cov=pd.DataFrame({'cell_id':ids,'stage':stage,'roi_id':name,'atlas_identity':cells.atlas_identity.to_numpy(),
            'inner_evaluation_cell':inner.to_numpy(),'native_molecule_coverage_fraction':fraction,
            'alternative_region_touches_context_edge':touches[1:].any(0),'comparison_eligible':reliable,
            'nuclear_area_missing':cells.nucleus_area.isna().to_numpy()})
        coverage.append(cov)
        np.savez_compressed(folder/(name+'_ROI_count_variants.npz'),cell_ids=ids,genes=genes,comparison_eligible=reliable,**matrices)
        print('REASSIGNED_ROI',stage,name,'inner',int(inner.sum()),'reliable',int(reliable.sum()),flush=True)
    pd.DataFrame(summaries).to_csv(folder/'ROI_assignment_summary.tsv',sep='\t',index=False)
    pd.DataFrame(transitions).to_csv(folder/'ROI_assignment_transitions.tsv.gz',sep='\t',index=False)
    pd.concat(coverage,ignore_index=True).to_csv(folder/'ROI_cell_comparison_coverage.tsv.gz',sep='\t',index=False)
    (folder/'ROI_reassignment_completion.json').write_text(json.dumps({'stage':stage,'ROIs':len(rois),
        'comparison_rule':'Inner centroid,>=95% native molecule coverage,and neither alternative region touching contextual crop edge; nucleus never required',
        'source_molecule_completion_sha256':digest(folder/'molecule_completion.json'),
        'source_morphology_completion_sha256':digest(folder/'morphology_completion.json'),
        'native_raster_control':'Additional0.5µm numerical control to distinguish grid-assignment effects from changed boundaries',
        'QV_threshold':20,'molecule_counts_regenerated':True},indent=2)+'\n')
    print('ROI_REASSIGNMENT_COMPLETE',stage,flush=True)


if __name__=='__main__':main()
