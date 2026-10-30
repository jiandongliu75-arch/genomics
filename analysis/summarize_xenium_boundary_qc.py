"""Aggregate all completed boundary-QC outputs without selecting best settings."""
import importlib.metadata
import json
import sys
import numpy as np
import pandas as pd

from design_xenium_boundary_qc import OUT,STAGES


def main():
    collections={name:[] for name in ['Ankrd1_boundary_sensitivity','identity_and_information_retention','boundary_identity_transitions',
        'boundary_marker_localization_summary','Ankrd1_depth_control_sensitivity','depth_control_identity_retention',
        'ROI_assignment_summary','segmentation_geometry_comparison']}
    run=[];profiles=[];facing=[]
    for stage in STAGES:
        folder=OUT/stage
        assert (folder/'sensitivity_completion.json').exists() and (folder/'depth_controls_completion.json').exists()
        for name,parts in collections.items():
            path=folder/(name+('.tsv.gz' if name=='segmentation_geometry_comparison' else '.tsv'))
            parts.append(pd.read_csv(path,sep='\t'))
        r=json.loads((folder/'verified_gene_QV_metrics.json').read_text())
        metrics=json.loads((folder/'molecule_completion.json').read_text())
        source_metrics=pd.read_csv('/workspace/heart-atlas/results/xenium_boundary_inputs/'+stage+'/metrics_summary.csv').iloc[0]
        run.append({'stage':stage,'all_decoded_molecule_rows':r['rows'],'gene_Q20_fraction':r['observed_gene_only_Q20_fraction'],
            'gene_Q20_transcripts':r['Q20_gene_transcripts'],'assigned_Q20_gene_transcripts':r['Q20_assigned_genes'],
            'assigned_Q20_gene_fraction':r['observed_Q20_gene_assigned_fraction'],
            'gene_Q30_transcripts':r['Q30_gene_transcripts'],'exact_source_counts_reconstructed':metrics['count_reconstruction']['reconstructed_QV20_integer_counts_exact'],
            'vendor_cell_stain_segmentation_fraction':float(source_metrics.segmented_cell_stain_frac),
            'vendor_nuclear_expansion_fraction':float(source_metrics.segmented_cell_nuc_expansion_frac),
            'original_target_cells':metrics['target_cells'],'cached_ROI_rows':sum(metrics['ROI_rows'].values())})
        cell=pd.read_csv(folder/'cell_shell_marker_profiles.tsv.gz',sep='\t')
        for lineage in ['Endothelial','Fibroblast']:
            p=cell[cell.atlas_identity.eq(lineage)]
            for region in ['deep_core','opposite_facing_shell','away_shell']:
                q=p[p.region.eq(region)];total=q.gene_transcripts.sum();foreign=q.counter_lineage_marker_transcripts.sum()
                facing.append({'stage':stage,'lineage':lineage,'region':region,'total_transcripts':int(total),
                    'CM_marker_transcripts':int(foreign),'CM_marker_fraction':float(foreign/max(total,1)),
                    'cells_with_transcripts':int(q.gene_transcripts.gt(0).sum())})
        target=pd.read_csv(folder/'boundary_target_cells.tsv.gz',sep='\t')
        for identity,p in target.groupby('atlas_identity'):
            profiles.append({'stage':stage,'atlas_identity':identity,'target_profiles':len(p),
                'without_nuclear_area':int(p.nucleus_area.isna().sum())})
    for name,parts in collections.items():pd.concat(parts,ignore_index=True).to_csv(OUT/(name+'.tsv'),sep='\t',index=False,na_rep='NA')
    pd.DataFrame(run).to_csv(OUT/'whole_run_molecule_QC.tsv',sep='\t',index=False)
    pd.DataFrame(facing).to_csv(OUT/'CM_marker_boundary_localization.tsv',sep='\t',index=False)
    pd.DataFrame(profiles).to_csv(OUT/'target_profile_nuclear_area_summary.tsv',sep='\t',index=False)
    a=pd.concat(collections['Ankrd1_boundary_sensitivity'],ignore_index=True)
    primary=a[a.scope.eq('geometry_targets_and_all_frozen_pairs')&a.pair_filter.eq('all_frozen_pairs')&a.nuclear_area_stratum.eq('all')]
    primary.to_csv(OUT/'primary_Ankrd1_assignment_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
    image=a[a.scope.eq('image_ROI_subset')&a.pair_filter.eq('all_frozen_pairs')&a.nuclear_area_stratum.eq('all')]
    image.to_csv(OUT/'image_ROI_Ankrd1_assignment_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
    versions={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','anndata','pyarrow','tifffile','scikit-image','imagecodecs','shapely','matplotlib','requests','remotezip']}
    (OUT/'runtime_versions.json').write_text(json.dumps({'python':sys.version,'packages':versions},indent=2)+'\n')
    print('ALL_BOUNDARY_RESULTS_AGGREGATED',flush=True)


if __name__=='__main__':main()
