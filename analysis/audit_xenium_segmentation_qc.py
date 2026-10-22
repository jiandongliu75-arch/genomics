"""Available Xenium cell-level QC diagnostics; does not infer new segmentation."""
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from build_whole_heart_spatial_atlas import OUT,SOURCE,STAGES,LINEAGES,REVIEW,digest
from plot_whole_heart_spatial_atlas import COLORS
from define_heart_arterial_structures import RAW


def main():
    summaries=[];samples=[];provenance=[]
    with PdfPages(OUT/'Xenium_available_QC_diagnostics.pdf') as pdf:
        for stage in STAGES:
            print('QC',stage,flush=True)
            folder=OUT/stage
            obs=pd.read_csv(folder/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
            data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
            alias='P07' if stage=='P7' else stage
            paths=[p for p in RAW.glob('*/cells.csv.gz') if p.parent.name.split('__')[0]=='WT_'+alias]
            assert len(paths)==1
            native=pd.read_csv(paths[0],index_col=0,float_precision='round_trip').loc[obs.index]
            np.testing.assert_array_equal(np.asarray(data.X.sum(axis=1)).ravel(),native.transcript_counts)
            np.testing.assert_allclose(native[['x_centroid','y_centroid']],obs[['x_centroid','y_centroid']],rtol=0,atol=1e-9)
            measured=np.asarray((data.X>0).sum(axis=1)).ravel()
            qc=obs[['stage','stage_cell_id','x_centroid','y_centroid','atlas_identity','pass_expression_QC','eligible_whole_reference']].copy()
            for c in native.columns:
                if c not in qc:qc[c]=native[c]
            qc['detected_target_genes']=measured
            qc['negative_probe_fraction']=qc.control_probe_counts/np.maximum(qc.transcript_counts+qc.control_probe_counts,1)
            qc['negative_codeword_fraction']=qc.control_codeword_counts/np.maximum(qc.total_counts,1)
            qc['unassigned_codeword_fraction']=qc.unassigned_codeword_counts/np.maximum(qc.total_counts,1)
            qc['nuclear_area_missing_flag']=qc.nucleus_area.isna()
            qc['nuclear_to_cell_area_ratio']=qc.nucleus_area/qc.cell_area
            qc['nuclear_area_exceeds_cell_flag']=qc.nucleus_area.gt(qc.cell_area*1.01)
            qc['transcripts_per_cell_area_um2']=qc.transcript_counts/np.maximum(qc.cell_area,1e-12)
            qc['CM_nonCM_marker_coexpression_flag']=obs.mixed_CM_nonCM_marker_review
            qc['extreme_area_review_flag']=False
            # Broad-lineage and stage-specific area review, not a global CM doublet cut.
            for label,part in qc.groupby('atlas_identity'):
                value=np.log1p(part.cell_area.to_numpy());median=np.median(value)
                mad=1.4826*np.median(np.abs(value-median))
                flag=np.abs(value-median)>5*max(mad,.05)
                if len(part)>=100:qc.loc[part.index,'extreme_area_review_flag']=flag
            expected=qc.transcript_counts.ge(30)&qc.detected_target_genes.ge(10)&qc.cell_area.gt(0)&np.isfinite(qc.x_centroid)&np.isfinite(qc.y_centroid)&qc.negative_probe_fraction.le(.1)
            np.testing.assert_array_equal(expected,qc.pass_expression_QC)
            qc['priority_segmentation_review_flag']=qc.CM_nonCM_marker_coexpression_flag|qc.extreme_area_review_flag|qc.nuclear_area_exceeds_cell_flag
            qc.to_csv(folder/'available_Xenium_QC_metrics.tsv.gz',sep='\t',na_rep='NA')
            for label in LINEAGES+REVIEW:
                part=qc[qc.atlas_identity.eq(label)]
                summaries.append({'stage':stage,'atlas_identity':label,'cells':len(part),
                    'median_transcripts':float(part.transcript_counts.median()),'median_detected_genes':float(part.detected_target_genes.median()),
                    'median_cell_area_um2':float(part.cell_area.median()),'median_nuclear_cell_ratio':float(part.nuclear_to_cell_area_ratio.median()),
                    'missing_nuclear_area_cells':int(part.nuclear_area_missing_flag.sum()),'nucleus_exceeds_cell_cells':int(part.nuclear_area_exceeds_cell_flag.sum()),
                    'extreme_area_review_cells':int(part.extreme_area_review_flag.sum()),
                    'marker_coexpression_review_cells':int(part.CM_nonCM_marker_coexpression_flag.sum()),
                    'source_QC_pass_cells':int(part.pass_expression_QC.sum()),
                    'control_probe_counts':int(part.control_probe_counts.sum()),'control_codeword_counts':int(part.control_codeword_counts.sum())})
                for flagged,group in part.groupby('priority_segmentation_review_flag'):
                    if len(group):samples.append(group.sample(n=min(5,len(group)),random_state=304).assign(sample_review_stratum='flagged' if flagged else 'comparison'))
            fig,axes=plt.subplots(2,2,figsize=(13,10))
            for label in LINEAGES:
                part=qc[qc.atlas_identity.eq(label)]
                axes[0,0].hist(np.log10(part.cell_area.clip(lower=.1)),bins=np.linspace(0,4,65),density=True,histtype='step',lw=1.4,color=COLORS[label],label=label)
                axes[0,1].hist(part.nuclear_to_cell_area_ratio.dropna(),bins=np.linspace(0,1.05,45),density=True,histtype='step',lw=1.4,color=COLORS[label])
            axes[0,0].set_xlabel('log10(cell area / µm²)');axes[0,0].set_ylabel('Within-lineage density');axes[0,0].legend(fontsize=7,frameon=False)
            axes[0,1].set_xlabel('Nuclear area / cell area');axes[0,1].set_ylabel('Within-lineage density')
            subset=qc.sample(n=min(12000,len(qc)),random_state=51)
            axes[1,0].scatter(subset.cell_area,subset.transcript_counts,c=[COLORS[v] for v in subset.atlas_identity],s=2,alpha=.35,rasterized=True)
            axes[1,0].set_xscale('log');axes[1,0].set_yscale('log');axes[1,0].set_xlabel('Cell area (µm²)');axes[1,0].set_ylabel('Measured transcripts per cell')
            axes[1,1].scatter(subset.x_centroid,subset.y_centroid,c='#dddddd',s=1,edgecolors='none',rasterized=True)
            missing=qc[qc.nuclear_area_missing_flag];axes[1,1].scatter(missing.x_centroid,missing.y_centroid,c='#b6487d',s=1,edgecolors='none',rasterized=True)
            axes[1,1].invert_yaxis();axes[1,1].set_aspect('equal');axes[1,1].set_xlabel('X (µm)');axes[1,1].set_ylabel('Y (µm)')
            axes[1,1].set_title(f'Missing nuclear area: {len(missing):,} ({len(missing)/len(qc):.1%})',fontsize=10)
            fig.suptitle(f'{stage}: available Xenium QC diagnostics',fontsize=18,weight='bold')
            fig.text(.5,.02,'Diagnostics preserve all original cells and counts. Area/coexpression flags prioritize review; they do not establish doublets.\nMissing nuclear area is not automatically rejected. No transcript-level QV audit, transcript reassignment or new segmentation was performed.',ha='center',fontsize=9)
            fig.subplots_adjust(left=.08,right=.96,top=.91,bottom=.12,hspace=.30,wspace=.30)
            fig.savefig(OUT/f'{stage}_available_Xenium_QC.png',dpi=160);pdf.savefig(fig);plt.close(fig)
            provenance.append({'stage':stage,'source_cell_table':str(paths[0]),'sha256':digest(paths[0]),'original_cells_checked':len(qc),
                'QC_filters_recomputed':True,'molecule_QV_audited':False,'cell_boundaries_resegmented':False})
    pd.DataFrame(summaries).to_csv(OUT/'available_Xenium_QC_summary.tsv',sep='\t',index=False,na_rep='NA')
    pd.concat(samples).to_csv(OUT/'candidate_segmentation_review_sample.tsv.gz',sep='\t',na_rep='NA')
    (OUT/'available_QC_provenance.json').write_text(json.dumps({'stages':provenance,
        'unassigned_codewords_are_not_unassigned_cell_transcripts':True,
        'area_review_rule':'Within stage and assigned broad category, abs(log1p area - median)>5*MADsigma;>=100 profiles required; review only',
        'sample':'Up to5 profiles per stage×category×flagged/comparison stratum,seed304; selected for future review, not manually adjudicated'},indent=2)+'\n')
    print('AVAILABLE_QC_DIAGNOSTICS_COMPLETE',flush=True)


if __name__=='__main__':main()
