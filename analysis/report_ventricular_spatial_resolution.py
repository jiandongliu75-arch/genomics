"""Render finest exploratory CM maps, evidence layers and resolution curves."""
import argparse
import importlib.metadata
import json
import os
import platform

os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from push_ventricular_spatial_resolution import OUT,PREVIOUS,SOURCE,ROOT,STAGES,digest
from plot_xenium_ventricular import background,frame

NAMES={
 'P0_M0':'Cycling','P0_M1':'Ddc / metabolic transition','P0_M2':'Ncam1 / neonatal working',
 'P0_M3':'Ankrd1 remodeling A','P0_M4':'Fatty-acid-use associated','P0_M5':'Conduction-program candidate',
 'P0_M6':'Ankrd1 remodeling B','P0_M7':'DNA-maintenance associated',
 'P7_M0':'Ppara / working CM','P7_M1':'Ankrd1 remodeling','P7_M2':'Regional electrical candidate',
 'P7_M3':'Cycling','P7_M4':'Ivns1abp / working CM','P7_M5':'DNA-maintenance associated',
 'P7_M6':'Conduction-program candidate','P7_M7':'Uncharacterized A','P7_M8':'Uncharacterized B',
 'P14_M0':'Slit2 / regional electrical','P14_M1':'Working / maturation','P14_M2':'Ankrd1 remodeling',
 'P14_M3':'Repolarization-associated RNA group','P14_M4':'Conduction-program candidate',
 'P21_M0':'Working / electrical candidate','P21_M1':'Working CM, weak positive markers',
 'P21_M2':'Ankrd1 remodeling A','P21_M3':'Slit2 / regional electrical',
 'P21_M4':'Ankrd1 remodeling B','P21_M5':'Conduction-program candidate'}


def mdtable(df):
    return '\n'.join(['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']+
        ['| '+' | '.join(map(str,r))+' |' for r in df.itertuples(index=False,name=None)])


def save(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=175,facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),facecolor='white');pdf.savefig(fig,facecolor='white');plt.close(fig)


def category(ax,obs,labels,mask,classes,colors):
    ax.scatter(obs.x_centroid,obs.y_centroid,s=.25,color='#e6e6e6',edgecolors='none',rasterized=True)
    handles=[]
    for i,label in enumerate(classes):
        select=mask&labels.eq(label)
        ax.scatter(obs.loc[select,'x_centroid'],obs.loc[select,'y_centroid'],s=.9,color=colors[i],edgecolors='none',rasterized=True)
        text=f'{label}: {NAMES.get(label,"merged union")} ({select.sum():,})' if len(label)<12 else f'Union {i+1}: {label} ({select.sum():,})'
        handles.append(Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[i],markersize=5,label=text))
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.01,1),fontsize=7.5,frameon=False)


def figures(summary,groups,frontier):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    with PdfPages(OUT/'ventricular_fine_resolution_report.pdf') as pdf:
        fig,axes=plt.subplots(2,4,figsize=(20,10))
        for col,stage in enumerate(STAGES):
            ax=axes[0,col]
            table=frontier[frontier.stage.eq(stage)]
            for split,color in [('validation','#a2b4c3'),('assessment','#237f85')]:
                part=table[table.phase.eq(split)].sort_values('classes_at_cut')
                ax.plot(part.classes_at_cut,part.macro_f1,'o-',label=split,color=color,lw=2)
            k=int(summary.set_index('stage').loc[stage,'validation_selected_stagewide_groups'])
            ax.axvline(k,color='#d17a28',ls='--',lw=1.5,label='Chosen on validation')
            ax.set_title(stage,weight='bold');ax.set_ylim(0,1.02);ax.set_xlabel('Groups retained at hierarchy cut');ax.grid(alpha=.15)
            if col==0:ax.set_ylabel('Mean per-class RNA F1');ax.legend(fontsize=8,loc='lower left')
            ax=axes[1,col];g=groups[groups.stage.eq(stage)].sort_values('group')
            positions=np.arange(len(g))
            ax.scatter(positions,g.min_precision,color='#8d9dad',label='Unrestricted RNA class',s=35)
            selected=g.selected_weight_threshold.notna()
            ax.scatter(positions[selected],g.loc[selected,'assessment_min_precision'],marker='D',color='#188f75',label='Frozen selective threshold',s=40)
            ax.axhline(.75,ls=':',color='#188f75',lw=1)
            ax.set_xticks(positions,[x.split('_')[1] for x in g.group],rotation=45)
            ax.set_ylim(-.02,1.02);ax.set_xlabel('Original fine RNA group');ax.grid(alpha=.15)
            if col==0:ax.set_ylabel('Minimum assessment precision across libraries');ax.legend(fontsize=7,loc='lower left')
        fig.suptitle('How far can the shared panel distinguish the existing ventricular RNA groups?',fontsize=18,weight='bold')
        fig.text(.5,.025,'Hierarchy merges use calibration nuclei; cut selection uses validation nuclei; the assessment quarter evaluates the frozen choices.\n'
            'Bottom: selective precision is shown only where a threshold was chosen on validation. These are internal RNA-label tests, not Xenium anatomical accuracy.',ha='center',fontsize=10)
        fig.subplots_adjust(left=.06,right=.985,top=.92,bottom=.13,wspace=.28,hspace=.40)
        save(fig,'resolution_and_selective_evidence',pdf)
        for stage in STAGES:
            folder=OUT/stage
            obs=pd.read_csv(folder/'fine_CM_spatial_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
            with np.load(folder/'fine_probabilities_and_sensitivity.npz') as b:
                ids,classes,p=b['query_cell_ids'].astype(str),b['classes'].astype(str).tolist(),b['full_panel']
            bg,extent,bounds=background(stage,obs)
            colors=[plt.get_cmap('tab10')(i) for i in range(len(classes))]
            mask=obs.included_in_weighted_mapping
            fig,axes=plt.subplots(2,2,figsize=(23,18))
            frame(axes[0,0],bg,extent,bounds,'Finest exploratory proposals: every existing RNA group')
            category(axes[0,0],obs,obs.fine_group_proposal,mask,classes,colors)
            frame(axes[0,1],bg,extent,bounds,'Fine candidates with either internal support check')
            category(axes[0,1],obs,obs.fine_group_proposal,obs.prioritized_fine_candidate,classes,colors)
            frame(axes[1,0],bg,extent,bounds,'Coarser cut selected on reserved RNA validation nuclei')
            unions=sorted(obs.loc[mask,'selected_resolution_union'].unique())
            category(axes[1,0],obs,obs.selected_resolution_union,mask,unions,[plt.get_cmap('tab10')(i) for i in range(len(unions))])
            frame(axes[1,1],bg,extent,bounds,'Fine-label weight concentration, without smoothing')
            im=axes[1,1].scatter(obs.loc[mask,'x_centroid'],obs.loc[mask,'y_centroid'],c=obs.loc[mask,'fine_weight'],s=.9,
                vmin=0,vmax=1,cmap='viridis',edgecolors='none',rasterized=True)
            fig.colorbar(im,ax=axes[1,1],fraction=.035,pad=.02,label='Maximum neighbor weight (not probability)')
            fig.suptitle(f'{stage}: pushing ventricular spatial mapping to {len(classes)} existing RNA groups',fontsize=21,weight='bold')
            fig.text(.5,.035,'Top left is a permissive hypothesis map; top right adds internal RNA support, library/raw consistency and measured-marker co-detection.\n'
                'The merged unions are computational resolution checks, not newly discovered biological cell types. No coordinates are moved or labels spatially smoothed.\n'
                'Conduction labels remain candidates. Kcnd2/Irx5 are unassayed; all working/electrical labels originate from the RNA reference.',ha='center',fontsize=11)
            fig.subplots_adjust(left=.055,right=.82,top=.94,bottom=.12,wspace=.85,hspace=.27)
            save(fig,f'{stage}_finest_CM_overview',pdf)
            rows=int(np.ceil(len(classes)/3));fig,axes=plt.subplots(rows,3,figsize=(19,rows*7.4),squeeze=False)
            evidence=groups[groups.stage.eq(stage)].set_index('group')
            for j,(ax,label) in enumerate(zip(axes.ravel(),classes)):
                frame(ax,bg,extent,bounds,f'{label} · {NAMES[label]}')
                order=np.argsort(p[:,j])
                coordinates=obs.loc[ids,['x_centroid','y_centroid']].to_numpy()
                ax.scatter(coordinates[order,0],coordinates[order,1],c=p[order,j],cmap='viridis',vmin=0,vmax=1,
                    s=.9,edgecolors='none',rasterized=True)
                chosen=obs.prioritized_fine_candidate&obs.fine_group_proposal.eq(label)
                ax.scatter(obs.loc[chosen,'x_centroid'],obs.loc[chosen,'y_centroid'],s=2.2,facecolors='none',edgecolors='black',linewidths=.23,rasterized=True)
                row=evidence.loc[label];markers=str(row.panel_supporting_genes).split(';')[:3]
                marker_text=', '.join(markers) if markers!=['nan'] else 'none meeting marker rule'
                ax.text(.01,-.13,f'Proposed: {int(row.proposed_cells):,}; prioritized: {int(row.prioritized_fine_candidate_cells):,}\nMeasured supporting genes: {marker_text}',transform=ax.transAxes,fontsize=8,va='top')
                ax.title.set_fontsize(10)
            for ax in axes.ravel()[len(classes):]:ax.axis('off')
            cax=fig.add_axes([.945,.23,.012,.54]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=cax)
            cb.set_label('RNA-neighbor weight for this group; not calibrated probability',fontsize=9)
            fig.suptitle(f'{stage}: individual spatial evidence maps for all {len(classes)} ventricular RNA groups',fontsize=19,weight='bold',y=.98)
            fig.text(.5,.023,'All original locations are retained. Colors show continuous group weight; black outlines mark prioritized candidate placements.\n'
                'Measured markers corroborate labels but also contributed to matching. Low weights and absent outlines do not establish biological absence.',ha='center',fontsize=10)
            fig.subplots_adjust(left=.055,right=.915,top=.94,bottom=.09,wspace=.29,hspace=.45)
            save(fig,f'{stage}_individual_fine_CM_maps',pdf)
        comparison=OUT/'comparison_with_published_states.png'
        if comparison.exists():
            fig,ax=plt.subplots(figsize=(13,7));ax.imshow(plt.imread(comparison));ax.axis('off')
            fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
            pdf.savefig(fig,facecolor='white');plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--skip-plots',action='store_true');args=parser.parse_args()
    validation=json.loads((OUT/'validation.json').read_text());assert validation['status']=='PASS'
    summary=pd.read_csv(OUT/'stage_resolution_summary.tsv',sep='\t')
    groups=pd.read_csv(OUT/'fine_group_evidence_summary.tsv',sep='\t')
    frontier=pd.read_csv(OUT/'resolution_frontier.tsv',sep='\t')
    if not args.skip_plots:figures(summary,groups,frontier)
    rows=[]
    for row in summary.itertuples():
        selected=groups[groups.stage.eq(row.stage)&groups.prioritized_fine_candidate_cells.gt(0)]
        rows.append({'Stage':row.stage,'Fine RNA groups mapped':row.fine_RNA_groups_tested,
            'Groups with prioritized local candidates':row.prioritized_fine_groups,'Prioritized cells':f'{row.prioritized_fine_candidate_cells:,}',
            'Stage-wide validation-selected groups':row.validation_selected_stagewide_groups})
    lines=['# Finest ventricular CM placements supported by the existing Xenium panel','',
        '**All 28 existing stage-specific ventricular RNA groups now have spatial maps**, including continuous weight maps for every individual group, permissive full-resolution proposals and evidence-prioritized subsets. No new RNA clusters were invented, original RNA memberships were unchanged, and no spatial coordinates were moved.','',
        mdtable(pd.DataFrame(rows)),'',
        '**Fine local placement and stage-wide separability are different questions.** Some fine groups have identifiable subsets even when a clean subdivision of every cell is unsupported. The rightmost column tests a whole-stage partition in which every union must meet minimum library-recovery criteria; it is not the number of biological cell types present or an upper bound on all local spatial detail.','',
        '[Ten-page PDF: resolution tests, four overview maps, four individual-group galleries and published-label comparison](ventricular_fine_resolution_report.pdf) · [resolution figure](resolution_and_selective_evidence.png) · [28-group evidence table](fine_group_evidence_summary.tsv)','',
        '**The cell-by-cell comparison reveals substantial lineage conflicts:** 10,842 of13,757 prioritized cells (78.8%) are also labeled CM in the published atlas; 2,899 (21.1%) have a published non-CM label and16 are absent. These conflicts are especially important for P7 cycling and P14 remodeling candidates. Internal RNA consistency did not establish their CM lineage on Xenium. Every conflict is flagged in the per-cell table; `prioritized_and_published_CM_compatible` identifies the annotation-compatible subset without silently removing the conflicting results.','',
        '## Finest placements and their criteria','',
        'A permissive proposal is supplied for every included cell using its largest fine-group RNA-neighbor weight. These are hypothesis maps at the maximum available reference resolution: 8 groups at P0, 9 at P7, 5 at P14 and 6 at P21. Numerical IDs are stage-specific; the same M number at different ages does not imply the same subtype.','',
        'The prioritized map retains either of two documented evidence routes. Both require agreement with separate-library projections and the original-count reference, plus at least two directly measured supporting genes. Route A requires weight≥0.70, margin≥0.20 and acceptable unrestricted fine-class recovery on the reserved RNA assessment nuclei. Route B uses a class-specific weight threshold chosen only on the validation quarter; the frozen threshold must then achieve≥75% selective precision with≥10 called nuclei in each library on the assessment quarter, with margin≥0.10. These complementary rules identify candidate placements; they are not calibrated probabilities or independent Xenium anatomical validation.','',
        'Other cells retain their proposals, weights, margins, entropy, original-reference result, separate-library results and two nonoverlapping gene-half results. The `placement_evidence` tiers additionally expose library/raw-consistent and permissive exploratory placements. We do not discard those weaker signals when showing the resolution limit.','']
    for stage in STAGES:
        part=groups[groups.stage.eq(stage)].sort_values('group')
        table=[]
        for r in part.itertuples():
            markers=str(r.panel_supporting_genes).split(';')[:3]
            table.append({'Group':r.group,'RNA interpretation':NAMES[r.group],'Proposed cells':f'{r.proposed_cells:,}',
                'Prioritized cells':f'{r.prioritized_fine_candidate_cells:,}',
                'Measured supporting genes':', '.join(markers) if markers!=['nan'] else 'None meeting rule',
                'RNA precision, unrestricted / selective':f'{r.min_precision:.2f} / '+(f'{r.assessment_min_precision:.2f}' if pd.notna(r.assessment_min_precision) else '—')})
        lines += [f'## {stage}','',mdtable(pd.DataFrame(table)),'',
            f'[Overview PNG]({stage}_finest_CM_overview.png) · [Every-group gallery PNG]({stage}_individual_fine_CM_maps.png) · [Gallery PDF]({stage}_individual_fine_CM_maps.pdf) · [All-cell annotations]({stage}/fine_CM_spatial_annotations.tsv.gz)','']
    lines += ['## Agreement with the published atlas','',
        '![Published final labels on our prioritized candidates](comparison_with_published_states.png)','',
        '[Detailed overlap by group](published_state_overlap_summary.tsv) · [broad-lineage compatibility](published_broad_lineage_compatibility.tsv) · [comparison PDF](comparison_with_published_states.pdf)','',
        'Agreement is partial. Among all prioritized cells, the P0/P7/P14 conduction-program groups carry the published Slit2.CM label in71.1% /79.6% /80.7% of cases. Cycling agrees with p.CM in86.8% at P0 and58.9% at P7. P14 Ankrd1-associated candidates carry Ankrd1.CM in45.3%. The regional electrical groups mostly fall within published working v.CM (54.0% at P14,61.9% at P21), with smaller Slit2.CM fractions; they do not correspond one-to-one to the authors’ Slit2 state. These percentages retain non-CM conflicts in the denominator and are annotation overlap, not accuracy.','',
        'The lineage discrepancies can reflect our ventricular gate, reference-transfer errors, segmentation/mixed RNA or differences in author annotation. This comparison cannot identify the correct label for each disputed cell. It is a reason to review the disputed placements before biological interpretation, rather than asserting that a fine-group evidence flag validates CM identity. The published labels were not used in fitting or evidence selection; filtering to the compatible subset is an explicit annotation-consensus restriction, not an independent validation.','',
        'The precision columns are **internal RNA-label recovery**, with the smaller value across the two library directions shown. A selective value appears only when a threshold was chosen using validation; values below0.75 do not pass Route B. Near-cutoff cases remain visible rather than being declared biologically absent. A zero or missing value can reflect inadequate rare-class support as well as poor recovery.','',
        'Positive markers are selected from the measured shared panel: group-minus-rest mean log1p(CP10k) difference≥0.15 in both libraries using both raw and corrected reference counts, and raw detection≥15% within the group in both libraries. Up to six genes are retained; tables show the first three. Generic chamber anchors are excluded. These marker checks use genes available to matching and are descriptive corroboration, not held-out validation. Full marker/effect tables are in each stage folder.','',
        'The electrical/repolarization names originate from full RNA profiles. **Kcnd2 and Irx5 are not measured by Xenium** and do not locate or orient these cells. Conduction-program labels are not verified Purkinje identities. Ankrd1 does not establish actual wall stress, and metabolic RNA states do not establish substrate flux. Uncharacterized P7 groups remain uncharacterized.','',
        '## Measuring the resolution frontier','',
        'Within each original RNA library×fine group, hashed barcode order defines disjoint calibration halves, validation quarters and assessment quarters. For each library direction, the reference contains only calibration nuclei from the opposite library; gene-scaling statistics also use calibration nuclei. Validation and assessment nuclei never enter those reference pools or fitted scaling statistics. All381 shared genes are used for this technical label-recovery test. Existing RNA labels were defined previously using the broader RNA data, so the labels themselves are not independently discovered truth.','',
        'Starting at the original fine resolution, the most mutually confused pair is merged using class-normalized calibration confusion only. The largest cut whose every union meets the fixed criteria on validation is selected: minimum precision0.60, recall0.35, F1 0.50, at least10 query and15 reference nuclei per library. The assessment quarter evaluates this frozen cut. All four selected cuts pass their assessment checks here; P21 selects only the unsplit union under this demanding whole-stage criterion. This does not erase its selectively supported regional group.','',
        'The selected unions are computational combinations of existing groups, not newly claimed biological types. Full memberships, library provenance, per-class metrics at every cut, RNA split probabilities and the resolution curve are downloadable. In particular, the original Ankrd1 A/B splits and many working/metabolic subdivisions are still chiefly exploratory on this panel.','',
        'Two deterministic, nonoverlapping halves of the shared genes provide an additional sensitivity check. A cell with zero RNA in a half is explicitly marked insufficient for that half. These halves are neither tuned to favor specific states nor used to impose spatial neighborhoods. Their agreement is exposed in the tables, not required as a universal filter.','',
        'The preceding [withheld-gene benchmark](../ventricular_alignment_improvement/README.md) establishes improved measured-expression prediction for the weighted projection. This new study addresses a different question: how much of the existing fine RNA labeling the panel can distinguish, and where its candidate placements lie. Prediction improvement alone is not treated as proof of fine-subtype identity.','',
        '## Files and reproducibility','',
        '- [Stage resolution summary](stage_resolution_summary.tsv), [all group-level evidence](fine_group_evidence_summary.tsv), [all resolution cuts](resolution_frontier.tsv), [numerical validation](validation.json), [design](design.json), [local threshold design](local_specificity_design.json), [methods](methods.json), [checksums](manifest.json).',
        '- Per-stage `fine_CM_spatial_annotations.tsv.gz` retains every original cell and centroid, the published author state as a separate comparator, fine proposals/weights, sensitivity results, markers detected, evidence tiers, prioritized-candidate flags and selected union. Exclusions remain in the table. Prefix stage before pooling section-specific barcodes.',
        '- `fine_probabilities_and_sensitivity.npz` contains query IDs, class order, full-panel weights and original-reference/gene-half weights. `RNA_split_probabilities.npz` contains RNA IDs, group truths, split assignments, library IDs and probabilities. These are weights, not measured expression or calibrated cell-identity probabilities.',
        '- `hierarchy.json` records every cut and the exact calibration reference/scaler IDs. `RNA_resolution_metrics.tsv` and `class_specific_RNA_thresholds.tsv` retain support counts and precision/recall diagnostics. [Published author-state comparison](../author_atlas_baseline/README.md) remains separate; author labels never enter fitting or threshold selection.','',
        '```bash',
        'OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python analysis/push_ventricular_spatial_resolution.py',
        'OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python analysis/calibrate_ventricular_spatial_resolution.py',
        'python analysis/compare_fine_CM_published_states.py',
        'python analysis/report_ventricular_spatial_resolution.py','```','',
        'Use the existing scientific environment and source matrices from the [weighted alignment report](../ventricular_alignment_improvement/README.md). Completed stage files skip fits; remove/rebuild the relevant derived outputs before deliberately changing inputs or design. Validation checks all original rows/coordinates, fine-weight conservation, RNA partitions/reference isolation, validation-only cut selection, measured-marker detection and sampled original/gene-half projections. No independently replicated animal-level spatial inference is claimed.','',
        'Sources: GSE290989 and Xenium Zenodo [P0/P7](https://zenodo.org/records/14961929), [P14/P21](https://zenodo.org/records/14961942); Wang, Dong, Song, Qian and Liu, CC BY4.0. [Article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/). DAPI and section plots are newly generated from deposited data; published final annotations are attributed to the article’s authors.','']
    prose='\n'.join(lines)
    for a,b in {'below0.75':'below 0.75','All381':'All 381','precision0.60':'precision 0.60','recall0.35':'recall 0.35',
        'least10':'least 10','and15':'and 15','CC BY4.0':'CC BY 4.0','of13,757':'of 13,757','and16':'and 16',
        'in71.1% /79.6% /80.7%':'in 71.1% / 79.6% / 80.7%','in86.8%':'in 86.8%','and58.9%':'and 58.9%',
        'in45.3%':'in 45.3%',',61.9%':', 61.9%'}.items():prose=prose.replace(a,b)
    (OUT/'README.md').write_text(prose)
    methods={'python':platform.python_version(),'versions':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','anndata','scikit-learn','matplotlib']},
        'fine_RNA_groups':28,'source_cohort':'168588 previously frozen ventricular-enriched candidates',
        'new_RNA_clustering':False,'source_counts_or_coordinates_modified':False,'spatial_smoothing':False,
        'author_labels_used_in_fit':False,'known_RNA_labels_independently_defined':False,
        'withheld_assessment_used_to_select_cut_or_threshold':False,'Xenium_probabilities_calibrated':False,
        'thresholds_transfer_from_RNA_to_Xenium_validated':False,'independent_anatomical_validation':False,
        'Kcnd2_Irx5_assayed':False,'design':json.loads((OUT/'design.json').read_text()),
        'local_design':json.loads((OUT/'local_specificity_design.json').read_text())}
    (OUT/'methods.json').write_text(json.dumps(methods,indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    scripts=[ROOT/'analysis'/n for n in ['push_ventricular_spatial_resolution.py','calibrate_ventricular_spatial_resolution.py','compare_fine_CM_published_states.py','report_ventricular_spatial_resolution.py']]
    (OUT/'manifest.json').write_text(json.dumps({'results':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(files)],
        'analysis':[{'path':str(p.relative_to(ROOT)),'sha256':digest(p)} for p in scripts],'validation':validation},indent=2)+'\n')
    print('REPORT_READY',len(files),flush=True)


if __name__=='__main__':main()
