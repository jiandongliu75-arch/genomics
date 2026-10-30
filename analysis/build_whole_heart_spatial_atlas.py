"""Audit broad non-CM assignments and export all-cell Xenium atlases.

Uses the completed full-reference projection; does not refit or force ambiguous
cells into a lineage. All expression counts and native coordinates are retained.
"""
from pathlib import Path
import importlib.metadata
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from reassess_whole_heart_lineages import OUT as FIT,ROOT,SOURCE,STAGES,PANELS,digest

OUT=ROOT/'results/whole_heart_spatial_atlas'
LINEAGES=['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
REVIEW=['Mixed-marker review','Ambiguous/review','QC/information flag']
ANCHORS={'Endothelial':'Pecam1 Cdh5 Esam Tek Vwf'.split(),
    'Fibroblast':'Dcn Lum Pdgfra Tcf21'.split(),
    'Pericyte':'Rgs5 Cspg4 Pdgfrb'.split(),
    'Smooth muscle':'Myh11 Cnn1'.split(),
    'Immune':'Ptprc Csf1r Mertk Cd22 Cd28 Cd8b1'.split(),
    'Epicardial':'Upk3b Wt1'.split()}
DISPLAY_GENES='Tnnt2 Tnni3 Pecam1 Cdh5 Dcn Pdgfra Rgs5 Cspg4 Myh11 Cnn1 Ptprc Csf1r Upk3b Wt1'.split()


def array_counts(data,genes):
    return data[:,[g for g in genes if g in data.var_names]].X.toarray()


def author_broad(label):
    if label.endswith('.CM'):return 'Cardiomyocyte'
    if '.EC' in label:return 'Endothelial'
    if '.Fib' in label:return 'Fibroblast'
    return {'Peri':'Pericyte','SMC':'Smooth muscle','Epic':'Epicardial',
        'Blood':'Blood (published)','Neur':'Neural (published)'}.get(label,'Not classified/present')


def write(frame,path,index=False):
    frame.to_csv(path,sep='\t',index=index,na_rep='NA')


def main():
    OUT.mkdir(exist_ok=True)
    design={'source':'Completed 377-gene full-reference projection in whole_heart_Ankrd1_vessels; no refitting',
        'inherited_identity_thresholds':'weight>=0.70,margin>=0.15,canonical marker support,raw-reference and both RNA-library agreement',
        'nonCM_anchor_panels':ANCHORS,'anchor_rule':'At least one detected lineage anchor in an already resolved non-CM; no upward reassignment of ambiguous cells',
        'mixed_rule':'Previously resolved non-CM with>=3 CM genes and>=10% CM RNA plus>=3 non-CM panel genes and>=10% non-CM RNA moves to Mixed-marker review; coexpression is not proof of a doublet',
        'CM_rule':'Retain previous strict CM assignments; previous mixed-CM assignments are Mixed-marker review',
        'other_cross_lineage_flags':'Two non-CM panels with>=3 genes and>=10% RNA each flag coexpression but do not automatically demote; related vascular/stromal programs can overlap',
        'glial_limit':'No Plp1/Sox10/S100b/Mpz in the403-target panel; no supported glial assignments asserted',
        'coordinates':'Native Xenium x/y micrometers; no warping or cross-stage registration',
        'neighborhood':'20 nearest other original cells within100µm; summarize resolved-neighbor types when>=5 resolved neighbors; descriptive, not causal or animal-level inference',
        'tile_size_um':100,'tile_description':'Original-cell count and resolved-lineage composition per100µm square; not a histological tissue-area denominator',
        'figures':'Direct categorical assignments, consistent stage colors; no label smoothing or imputation',
        'column_names':'Replace slashes in column names with underscores for HDF5-safe keys; categorical identity values are unchanged'}
    (OUT/'atlas_design.json').write_text(json.dumps(design,indent=2)+'\n')
    counts=[];audits=[];markers=[];neighbors=[];tiles=[];comparisons=[];source_records=[];tech=[]
    for stage in STAGES:
        print('BUILD',stage,flush=True);folder=OUT/stage;folder.mkdir(exist_ok=True)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        base=pd.read_csv(FIT/stage/'whole_heart_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        np.testing.assert_array_equal(data.obs_names,base.index)
        np.testing.assert_array_equal(data.obs[['x_centroid','y_centroid']],base[['x_centroid','y_centroid']])
        obs=base.rename(columns={c:c.replace('/','_') for c in base.columns}).copy()
        obs['stage']=stage;obs['stage_cell_id']=stage+'::'+obs.index.astype(str)
        obs['previous_reassessed_lineage']=obs.reassessed_lineage
        obs['atlas_identity']=obs.reassessed_lineage.replace({'Mixed CM/non-CM; review':'Mixed-marker review'})
        obs['atlas_review_reason']=''
        for label,genes in ANCHORS.items():
            present=[g for g in genes if g in data.var_names]
            values=array_counts(data,present)
            detected=(values>0).sum(1)
            key=label.replace(' ','_');obs['anchor_'+key+'_genes']=detected
            previous=base.reassessed_lineage.eq(label)
            weak=previous&pd.Series(detected==0,index=obs.index)
            obs.loc[weak,'atlas_identity']='Ambiguous/review'
            obs.loc[weak,'atlas_review_reason']='No detected lineage-specific anchor'
            mixed=previous&base.mixed_CM_nonCM_marker_review
            obs.loc[mixed,'atlas_identity']='Mixed-marker review'
            obs.loc[mixed,'atlas_review_reason']='Strong CM and non-CM marker coexpression'
            audits.append({'stage':stage,'lineage':label,'previous_resolved_cells':int(previous.sum()),
                'no_anchor_cells':int(weak.sum()),'mixed_CM_nonCM_cells':int(mixed.sum()),
                'retained_after_both_checks':int(obs.atlas_identity.eq(label).sum()),'anchor_genes_assayed':','.join(present)})
        obs.loc[base.reassessed_lineage.eq('Mixed CM/non-CM; review'),'atlas_review_reason']='Strong CM and non-CM marker coexpression'
        obs.loc[base.reassessed_lineage.eq('Ambiguous/review'),'atlas_review_reason']='Inherited reference/marker/library uncertainty'
        obs.loc[base.reassessed_lineage.eq('QC/information flag'),'atlas_review_reason']='Inherited expression QC or training-information flag'
        panels=[]
        for label in ANCHORS:
            key=label.replace(' ','_')
            panels.append((obs['marker_'+key+'_genes'].ge(3)&obs['marker_'+key+'_fraction'].ge(.10)).to_numpy())
        obs['multiple_nonCM_marker_panels_flag']=np.column_stack(panels).sum(1)>=2
        obs['atlas_resolved']=obs.atlas_identity.isin(LINEAGES)
        obs['assay_Pearson_disagreement_flag']=obs.eligible_whole_reference&~obs.Pearson_method_agrees.eq(True)
        obs['author_broad_comparator']=obs.author_atlas_cell_state.astype(str).map(author_broad)
        # Neighbor identities use all original cell positions before selecting confident types.
        xy=obs[['x_centroid','y_centroid']].to_numpy();tree=cKDTree(xy)
        dist,ix=tree.query(xy,k=21,workers=2)
        self_index=np.arange(len(obs))[:,None]
        valid=(ix!=self_index)&(dist<=100)
        codes=pd.Index(LINEAGES).get_indexer(obs.atlas_identity)
        nc=np.column_stack([((codes[ix]==j)&valid).sum(1) for j in range(len(LINEAGES))])
        nvalid=nc.sum(1);obs['resolved_neighbors_within_100um_among_20_nearest']=nvalid
        obs['neighborhood_eligible']=obs.atlas_resolved&(nvalid>=5)
        global_fraction=np.bincount(codes[codes>=0],minlength=len(LINEAGES))/max(obs.atlas_resolved.sum(),1)
        for i,label in enumerate(LINEAGES):
            take=obs.neighborhood_eligible&obs.atlas_identity.eq(label);n=int(take.sum())
            total=int(nc[take].sum())
            for j,other in enumerate(LINEAGES):
                edge_count=int(nc[take,j].sum());fraction=edge_count/max(total,1)
                neighbors.append({'stage':stage,'focal_lineage':label,'neighbor_lineage':other,'eligible_focal_cells':n,
                    'resolved_neighbor_links':total,'links_to_neighbor_lineage':edge_count,'neighbor_fraction':fraction,
                    'global_resolved_fraction':global_fraction[j],
                    'observed_expected_ratio':fraction/global_fraction[j] if global_fraction[j]>0 else np.nan})
        local=pd.DataFrame({'tile_x':np.floor(xy[:,0]/100).astype(int),'tile_y':np.floor(xy[:,1]/100).astype(int),'identity':obs.atlas_identity.to_numpy()})
        tile=pd.crosstab([local.tile_x,local.tile_y],local.identity).reindex(columns=LINEAGES+REVIEW,fill_value=0)
        tile['original_cells']=tile.sum(axis=1);tile['resolved_cells']=tile[LINEAGES].sum(axis=1);tile=tile.reset_index();tile['stage']=stage
        tiles.append(tile)
        totals=np.asarray(data.X.sum(1)).ravel();values=data.X.toarray()
        for label in LINEAGES+REVIEW:
            take=obs.atlas_identity.eq(label).to_numpy();n=int(take.sum())
            counts.append({'stage':stage,'atlas_identity':label,'cells':n,'fraction_all_original':n/len(obs),
                'fraction_resolved':n/obs.atlas_resolved.sum() if label in LINEAGES else np.nan})
            if n:
                part=values[take];normal=np.log1p(1e4*part/np.maximum(totals[take,None],1))
                for j,gene in enumerate(data.var_names):
                    markers.append({'stage':stage,'atlas_identity':label,'cells':n,'gene':gene,
                        'detection_fraction':float((part[:,j]>0).mean()),'mean_raw_count':float(part[:,j].mean()),
                        'mean_log1pCP10k':float(normal[:,j].mean())})
            for author,nc_ in obs.loc[take,'author_broad_comparator'].value_counts().items():
                comparisons.append({'stage':stage,'atlas_identity':label,'published_comparator':author,'cells':int(nc_)})
        table=obs.copy();table.index.name='cell_id';write(table,folder/'all_cell_annotations.tsv.gz',True)
        # Keep the downloadable count atlas compact: all source genes/counts plus audited metadata.
        atlas=ad.AnnData(X=data.X.copy(),obs=obs.copy(),var=data.var.copy())
        atlas.obsm['spatial']=xy.copy();atlas.uns['atlas_design_json']=json.dumps(design)
        atlas.uns['spatial_units']='micrometers';atlas.uns['stage']=stage
        atlas.uns['source_count_matrix_sha256']=digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        atlas.uns['DAPI_overview_repository_path']=f'results/xenium_ventricular_mapping/{stage}/DAPI_overview.png'
        # Object booleans from joins must be nullable booleans or explicit strings for HDF5.
        for column in atlas.obs:
            if atlas.obs[column].dtype==object:
                atlas.obs[column]=atlas.obs[column].fillna('NA').astype(str)
        atlas.write_h5ad(folder/f'{stage}_whole_heart_atlas.h5ad',compression='gzip')
        source_records.append({'stage':stage,'cells':len(obs),'resolved_cells':int(obs.atlas_resolved.sum()),
            'genes':data.n_vars,'count_sum':int(data.X.sum()),'source_sha256':atlas.uns['source_count_matrix_sha256'],
            'previous_annotations_sha256':digest(FIT/stage/'whole_heart_cell_annotations.tsv.gz'),
            'DAPI_sha256':digest(SOURCE/stage/'DAPI_overview.png')})
        technical=pd.read_csv(FIT/stage/'cross_library_technical_metrics.tsv',sep='\t')
        tech.append(technical[technical.resolution.eq('lineage')])
        print('EXPORTED',stage,obs.atlas_identity.value_counts().to_dict(),flush=True)
        del data,atlas,values,nc,normal,obs,base
    for name,rows in [('lineage_counts.tsv',counts),('nonCM_identity_audit.tsv',audits),('measured_gene_profiles.tsv',markers),
        ('neighborhood_composition.tsv',neighbors),('published_broad_comparison.tsv',comparisons)]:write(pd.DataFrame(rows),OUT/name)
    write(pd.concat(tiles,ignore_index=True),OUT/'spatial_100um_tile_counts.tsv.gz')
    write(pd.concat(tech,ignore_index=True),OUT/'RNA_opposite_library_broad_metrics.tsv')
    (OUT/'source_provenance.json').write_text(json.dumps(source_records,indent=2)+'\n')
    (OUT/'runtime_versions.json').write_text(json.dumps({'python':sys.version,'packages':{k:importlib.metadata.version(k) for k in ['anndata','numpy','pandas','scipy','matplotlib','shapely']}},indent=2)+'\n')
    print('ATLAS_BUILD_COMPLETE',flush=True)


if __name__=='__main__':main()
