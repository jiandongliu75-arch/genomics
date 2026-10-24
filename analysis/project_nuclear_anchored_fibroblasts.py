"""Compare whole-cell and nuclear-only fibroblast transfer after segmentation."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/nuclear-anchored-mpl')
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from fibroblast_state_common import FAMILIES,UNRESOLVED
from map_fibroblast_states_Xenium import FibroblastProjection,annotate
from reassess_whole_heart_lineages import PANELS as BROAD_PANELS
from nuclear_anchored_Xenium_segmentation import ROOT,OUT,SOURCE,ATLAS,STAGES,write_json

REVIEWED=ROOT/'results/fibroblast_Xenium_reviewed_refit'


def main(stage):
    folder=OUT/stage
    if (folder/'fibroblast_projection_complete.json').exists():print('SKIP_PROJECTION',stage,flush=True);return
    data=ad.read_h5ad(folder/'nuclear_anchored_counts.h5ad')
    atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
    np.testing.assert_array_equal(data.obs_names,atlas.index)
    native=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad',backed='r')
    mask=atlas.atlas_identity.eq('Fibroblast').to_numpy();ids=atlas.index[mask]
    genes=data.var_names.to_numpy(dtype=str)
    ann=pd.read_csv(REVIEWED/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    values={'native':native.X[mask].toarray().astype(np.float32),
            'nuclear_anchored':data.X[mask].toarray().astype(np.float32)}
    for name in ['nuclear_only','high_confidence','depth_matched','nuclear_depth_matched']:
        values[name]=data.layers[name][mask].toarray().astype(np.float32)
    values['cytoplasmic_only']=values['nuclear_anchored']-values['nuclear_only']
    assert (values['cytoplasmic_only']>=0).all()
    model=FibroblastProjection(stage,values['native'],genes,ann)
    frames={};rows=[]
    for variant,counts in values.items():
        pred=model.predict(counts)
        f=annotate(counts,genes,model,pred);f.index=ids
        f['total_gene_molecules']=counts.sum(1);f['detected_genes']=(counts>0).sum(1)
        cm=counts[:,np.isin(genes,BROAD_PANELS['Cardiomyocyte'])]
        fb=counts[:,np.isin(genes,BROAD_PANELS['Fibroblast'])]
        f['CM_marker_fraction']=cm.sum(1)/np.maximum(counts.sum(1),1)
        f['FB_marker_fraction']=fb.sum(1)/np.maximum(counts.sum(1),1)
        f['mixed_CM_FB_marker_review']=((cm>0).sum(1)>=3)&(cm.sum(1)>=10)&(f.CM_marker_fraction>=.1)&((fb>0).sum(1)>=3)&(f.FB_marker_fraction>=.1)
        f.to_csv(folder/f'fibroblast_{variant}_assignments.tsv.gz',sep='\t',na_rep='NA')
        frames[variant]=f
        for fam in FAMILIES+[UNRESOLVED]:
            rows.append({'stage':stage,'variant':variant,'family':fam,'broad_FB_profiles':len(f),
                'informative_profiles':int(f.state_information_QC.sum()),
                'supported_family_profiles':int(f.mapped_family.eq(fam).sum()),
                'proposals':int(f.family_proposal.eq(fam).sum()),
                'median_gene_molecules':f.total_gene_molecules.median(),
                'mixed_marker_review_profiles':int(f.mixed_CM_FB_marker_review.sum())})
        print('PROJECTED',stage,variant,int(f.supported_family_assignment.sum()),flush=True)
    old=pd.read_csv(REVIEWED/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
    assert frames['native'].mapped_family.equals(old.mapped_family)
    # Separate nuclear-modality calibration is reported as a sensitivity, not
    # quietly substituted for the frozen-native-calibration comparison above.
    nuclear=values['nuclear_only'];informative=frames['nuclear_only'].state_information_QC.to_numpy()
    if informative.sum()>=30:
        nuclear_model=FibroblastProjection(stage,nuclear[informative],genes,ann)
        f=annotate(nuclear,genes,nuclear_model,nuclear_model.predict(nuclear));f.index=ids
    else:
        f=frames['nuclear_only'].copy()
    f.to_csv(folder/'fibroblast_nuclear_modality_calibrated_assignments.tsv.gz',sep='\t',na_rep='NA')
    primary=frames['nuclear_anchored']
    output=atlas.loc[ids,['x_centroid','y_centroid','nucleus_area']].copy()
    output['native_family']=old.mapped_family
    output['nuclear_anchored_family']=primary.mapped_family
    output['nuclear_anchored_supported']=primary.supported_family_assignment
    output['nuclear_only_informative']=frames['nuclear_only'].state_information_QC
    output['nuclear_only_calibrated_family']=f.mapped_family
    output['nuclear_only_calibrated_supported']=f.supported_family_assignment
    output['nuclear_anchored_mixed_marker_review']=primary.mixed_CM_FB_marker_review
    output['evidence_tier']='Fine family unresolved'
    cell=primary.supported_family_assignment&~primary.mixed_CM_FB_marker_review
    nuc=f.supported_family_assignment
    same=f.mapped_family.eq(primary.mapped_family)
    output.loc[cell,'evidence_tier']='Cell-supported; nuclear confirmation unavailable'
    output.loc[cell&frames['nuclear_only'].state_information_QC,'evidence_tier']='Cell-supported; nuclear fine identity unresolved'
    output.loc[cell&nuc&same,'evidence_tier']='Cell and nuclear family concordant'
    conflict=cell&nuc&~same
    output.loc[conflict,'evidence_tier']='Nuclear-cell family disagreement; review'
    output.loc[~cell&nuc&primary.family_proposal.eq(f.mapped_family)&~primary.mixed_CM_FB_marker_review,
               'evidence_tier']='Nuclear-supported candidate; whole-cell strict criteria unmet'
    output.loc[primary.mixed_CM_FB_marker_review,'evidence_tier']='Mixed CM-FB marker review'
    output['recommended_cell_supported_family']=np.where(cell&~conflict,primary.mapped_family,UNRESOLVED)
    output.to_csv(folder/'fibroblast_nuclear_cell_evidence.tsv.gz',sep='\t',na_rep='NA')
    pd.DataFrame(rows).to_csv(folder/'fibroblast_variant_summary.tsv',sep='\t',index=False)
    tiers=output.evidence_tier.value_counts().rename_axis('evidence_tier').reset_index(name='profiles')
    tiers.to_csv(folder/'fibroblast_evidence_tiers.tsv',sep='\t',index=False)
    make_plots(stage,output,frames,f)
    write_json(folder/'fibroblast_projection_complete.json',{'stage':stage,'broad_fibroblasts':len(ids),
        'native_supported':int(old.supported_family_assignment.sum()),
        'nuclear_anchored_supported':int(primary.supported_family_assignment.sum()),
        'nuclear_only_informative':int(informative.sum()),
        'nuclear_modality_calibrated_supported':int(f.supported_family_assignment.sum()),
        'cell_nuclear_concordant':int((cell&nuc&same).sum()),
        'nuclear_cell_disagreement':int(conflict.sum()),
        'variants_use_frozen_native_query_calibration':list(values),
        'separate_nuclear_calibration':'Sensitivity uses only nuclear-information-eligible FB profiles to estimate query gene scaling; same RNA reference, same marker gates and library requirements',
        'nuclear_evidence_is_not_independent_ground_truth':True,
        'native_refit_assignments_reproduced_exactly':True})
    print('FIBROBLAST_PROJECTION_COMPLETE',stage,tiers.to_dict('records'),flush=True)


def make_plots(stage,output,frames,nuclear):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from fibroblast_state_common import COLORS
    fig,axes=plt.subplots(1,3,figsize=(16,7),constrained_layout=True)
    for ax,series,title in [(axes[0],output.native_family,'Native supported labels'),
            (axes[1],output.recommended_cell_supported_family,'Nuclear-anchored cell-supported labels'),
            (axes[2],nuclear.mapped_family,'Nuclear-only calibrated sensitivity')]:
        for fam,c in COLORS.items():
            d=output.loc[series.eq(fam)]
            ax.scatter(d.x_centroid,d.y_centroid,s=1.5,c=c,rasterized=True)
        ax.set_title(title);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
    fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,label=f) for f,c in COLORS.items()],
               loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle(stage+' · Whole-section segmentation; gray retains uncertain fine identity')
    fig.savefig(OUT/stage/f'{stage}_fibroblast_nuclear_anchored_maps.png',dpi=180)
    fig.savefig(OUT/stage/f'{stage}_fibroblast_nuclear_anchored_maps.pdf');plt.close(fig)
    labels=['native','nuclear_anchored','high_confidence','depth_matched','nuclear_only','nuclear_depth_matched','cytoplasmic_only']
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    axes[0].bar(range(len(labels)),[frames[x].supported_family_assignment.sum() for x in labels],color='#267ea8')
    axes[0].set_xticks(range(len(labels)),labels,rotation=35,ha='right');axes[0].set_ylabel('Profiles passing unchanged fine-label gates')
    axes[1].boxplot([frames[x].total_gene_molecules.to_numpy() for x in labels],tick_labels=labels,showfliers=False)
    axes[1].tick_params(axis='x',rotation=35);axes[1].set_ylabel('Gene molecules per broad FB profile')
    fig.suptitle(stage+' · More stringent ownership can reduce RNA information')
    fig.savefig(OUT/stage/f'{stage}_assignment_information_comparison.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=STAGES);a=p.parse_args();main(a.stage)
