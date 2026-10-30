"""Plot and report the fixed internal alignment comparison."""
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
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from improve_ventricular_alignment import OUT,ROOT,SOURCE,STAGES,METHODS,digest
from plot_xenium_ventricular import background,frame,COLORS


def mdtable(df):
    return '\n'.join(['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']+
        ['| '+' | '.join(map(str,r))+' |' for r in df.itertuples(index=False,name=None)])


def save(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=170,facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),facecolor='white')
    pdf.savefig(fig,facecolor='white');plt.close(fig)


def plot(metrics,technical):
    short=['CytoSPACE','Nearest 1','Soft 30','Assay-scaled','Spatial','Shuffled spatial']
    colors=['#6c757d','#a5a9ad','#5299ca','#098f81','#b78241','#b5a3bf']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    with PdfPages(OUT/'alignment_improvement_report.pdf') as pdf:
        fig,axes=plt.subplots(2,4,figsize=(21,10))
        for col,stage in enumerate(STAGES):
            sub=metrics[metrics.stage.eq(stage)&metrics.gene_set.eq('new_test')]
            for row,field in enumerate(['cell_spearman','grid250um_spearman']):
                ax=axes[row,col]
                for i,method in enumerate(METHODS):
                    values=sub.loc[sub.method.eq(method),field].dropna().to_numpy()
                    jitter=np.linspace(-.13,.13,len(values))
                    ax.scatter(i+jitter,values,s=15,color=colors[i],alpha=.5,edgecolors='none')
                    ax.plot([i-.20,i+.20],[np.median(values)]*2,lw=3,color=colors[i])
                ax.axhline(0,color='#aaaaaa',lw=.7);ax.grid(axis='y',alpha=.15)
                ax.set_xticks(range(6),short,rotation=45,ha='right',fontsize=8)
                ax.set_title(stage,weight='bold')
                if col==0:ax.set_ylabel('Single-cell gene correlation' if row==0 else '250 µm regional gene correlation')
                ax.set_ylim((-.1,.55) if row==0 else (-.35,1.0))
        fig.suptitle('Prediction of 20 newly withheld genes: same cells and training genes for every method',fontsize=18,weight='bold')
        fig.text(.5,.018,'Each point is one withheld gene; bars are medians. All 39 held-out genes were excluded from fitting and normalization denominators.\n'
            'These are internal gene-recovery tests conditional on the frozen ventricular gate, not independent anatomical subtype accuracy. Shuffled spatial is a control.',ha='center',fontsize=10)
        fig.subplots_adjust(left=.06,right=.98,top=.92,bottom=.15,wspace=.25,hspace=.6)
        save(fig,'withheld_gene_benchmark',pdf)
        for stage in STAGES:
            obs=pd.read_csv(OUT/stage/'weighted_CM_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
            bg,extent,bounds=background(stage,obs)
            fig,axes=plt.subplots(1,2,figsize=(21,10))
            for ax,field,title in zip(axes,['weighted_family','family_candidate'],
                ['Highest-weight RNA family: all included ventricular candidates','Concentrated and library-consistent family candidates']):
                frame(ax,bg,extent,bounds,title)
                ax.scatter(obs.x_centroid,obs.y_centroid,s=.25,color='#e6e6e6',edgecolors='none',rasterized=True)
                handles=[]
                mask=obs.included_in_weighted_mapping
                labels=obs.loc[mask,field]
                for label in sorted(labels.unique()):
                    select=mask&obs[field].eq(label)
                    color=COLORS.get(label,'#bfbfbf')
                    ax.scatter(obs.loc[select,'x_centroid'],obs.loc[select,'y_centroid'],s=.9,color=color,edgecolors='none',rasterized=True)
                    handles.append(Line2D([0],[0],marker='o',color='none',markerfacecolor=color,markersize=5,label=f'{label} ({select.sum():,})'))
                ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.01,1),fontsize=8,frameon=False)
            fig.suptitle(f'{stage}: unsmoothed weighted projection of ventricular RNA programs',fontsize=19,weight='bold')
            fig.text(.5,.04,'Original Xenium coordinates and counts are preserved. Left: strongest reference family; right: ≥0.70 neighbor weight and the same family from both RNA libraries.\n'
                'Gray cells are ambiguous or outside the frozen cohort. Weight fractions are not calibrated probabilities. These flags differ from the previous CytoSPACE acceptance rules.\n'
                'Conduction-program candidates are not verified Purkinje cells. Kcnd2/Irx5 are unassayed; neither is used to orient these sections.',ha='center',fontsize=10)
            fig.subplots_adjust(left=.055,right=.81,top=.9,bottom=.12,wspace=.85)
            save(fig,f'{stage}_weighted_family_maps',pdf)
        fig,axes=plt.subplots(1,2,figsize=(13,6))
        for ax,resolution,title in zip(axes,['mapping_family','ventricular_group'],['RNA program families','Stage-specific RNA groups']):
            for method,color in [('native_CytoSPACE','#6c757d'),('Pearson_soft30','#5299ca'),('assay_scaled_soft30','#098f81')]:
                part=technical[technical.resolution.eq(resolution)&technical.method.eq(method)]
                values=part.groupby('stage').f1.mean().reindex(STAGES)
                ax.plot(STAGES,values,'o-',label=method,color=color,lw=2)
            ax.set_ylim(0,1);ax.set_title(title,weight='bold');ax.set_ylabel('Mean per-label F1 across both library folds');ax.grid(alpha=.15)
        axes[1].legend(loc='lower left',fontsize=8)
        fig.suptitle('Technical transfer between the two RNA libraries at each age',fontsize=16,weight='bold')
        fig.text(.5,.035,'All381 shared genes; each query nucleus is excluded from its opposite-library reference. Existing RNA labels are not independent biological truth.\n'
            'Native CytoSPACE values come from the preceding same-panel library benchmark; input ordering/partition realization can differ.',ha='center',fontsize=9)
        fig.subplots_adjust(bottom=.20,top=.85,wspace=.3)
        save(fig,'cross_library_technical_comparison',pdf)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--skip-plots',action='store_true');args=parser.parse_args()
    validation=json.loads((OUT/'validation.json').read_text());assert validation['status']=='PASS'
    metrics=pd.read_csv(OUT/'heldout_gene_metrics.tsv',sep='\t')
    summary=pd.read_csv(OUT/'heldout_summary.tsv',sep='\t')
    coverage=pd.read_csv(OUT/'mapping_coverage.tsv',sep='\t')
    newtechnical=pd.concat([pd.read_csv(OUT/s/'cross_library_metrics.tsv',sep='\t') for s in STAGES])
    native=pd.read_csv(ROOT/'results/cytospace_xenium/reference_panel_benchmark/heldout_library_metrics.tsv',sep='\t')
    native['method']='native_CytoSPACE'
    technical=pd.concat([native,newtechnical],ignore_index=True)
    technical.to_csv(OUT/'cross_library_technical_comparison.tsv',sep='\t',index=False)
    if not args.skip_plots:
        plot(metrics,technical)
    chosen=summary[summary.gene_set.eq('new_test')].set_index(['stage','method'])
    rows=[]
    for s in STAGES:
        a,b=chosen.loc[(s,'native_CytoSPACE')],chosen.loc[(s,'assay_scaled_soft30')]
        rows.append({'Stage':s,'CytoSPACE cell ρ':f'{a.cell_spearman:.3f}','Weighted cell ρ':f'{b.cell_spearman:.3f}',
            'CytoSPACE regional ρ':f'{a.grid250um_spearman:.3f}','Weighted regional ρ':f'{b.grid250um_spearman:.3f}'})
    cov=coverage[['stage','included_cells','family_evidence_flag_cells','fine_evidence_flag_cells']].copy()
    cov.columns=['Stage','Included ventricular candidates','Family evidence flag','Fine-group evidence flag']
    for c in cov.columns[1:]:cov[c]=cov[c].map('{:,}'.format)
    design=json.loads((OUT/'benchmark_design.json').read_text())
    lines=['# Improving ventricular RNA–Xenium alignment with internal tests','',
        '**Weighted reference projection improves held-out measured-gene recovery at all four ages in this comparison.** It combines multiple reference nuclei, adjusts gene scales within each assay and allows projected subtype proportions to differ from the source RNA proportions. This improves expression prediction; it does not by itself prove the correct anatomical placement of every fine subtype.','',
        mdtable(pd.DataFrame(rows)),'',
        'Values are median Spearman correlations across **20 newly withheld genes**, evaluated on exactly the same cells and 342 training genes. Regional values average cells in fixed 250 µm bins with at least20 cells. No gene test p-values or animal-level confidence intervals are implied. The small single-cell correlations remain a material limit; regional recovery is stronger.','',
        '[Six-page PDF report: benchmark, four section maps and library comparison](alignment_improvement_report.pdf) · [benchmark PNG](withheld_gene_benchmark.png) · [all gene-level metrics](heldout_gene_metrics.tsv)','',
        '![Newly withheld-gene benchmark](withheld_gene_benchmark.png)','',
        '## What changed and why','',
        '1. **Allow subtype proportions to differ.** The weighted projection does not enforce the subtype mixture in the sampled RNA reference. Native CytoSPACE duplicate sampling inherits that mixture under the supplied broad label. The reference composition can still influence nearest neighbors, so absence of a quota does not remove all reference bias.',
        '2. **Combine evidence from30 nuclei.** A softmax-weighted set of similar RNA nuclei reduces the influence of sparse counts in any one nucleus. This is a weighted expression/program projection, not a one-to-one physical cell correspondence. The plain30-neighbor comparison separates this averaging benefit from assay scaling.',
        '3. **Adjust assay differences using training genes.** Log-normalized training genes are centered/scaled separately in the same-stage RNA and Xenium populations, with an SD floor0.25 and clipping at±5; cosine matching then uses30 neighbors with temperature0.10. All settings were fixed before reading the new test results. This is an unsupervised feature adjustment, not a fitted gene-wise sequencing-efficiency model. Differences in true population composition can also affect these statistics.',
        '4. **Report ambiguity.** Each included cell retains family/fine weight fractions and independent assignments using each RNA library alone. A candidate flag requires≥0.70 pooled neighbor weight and the same label from both libraries. Fine flags additionally require the family flag. These are descriptive evidence rules, not calibrated probabilities of correct Xenium identity.',
        '5. **Test spatial smoothing against a control.** A fixed light spatial blend was tested alongside a coordinate-shuffled graph. Their similar gains indicate that averaging, rather than verified anatomical neighborhood information, explains much of the change. The published primary maps therefore use the **unsmoothed** weighted projection.','',
        '## Internal evaluation without circular marker validation','',
        'The new test set contains: '+', '.join(design['new_test_genes'])+'. It was fixed by a deterministic gene-name hash, excluding the previously evaluated19 marker genes and a list of direct chamber/lineage anchors, before fitting. These20 genes were not used to choose hyperparameters. All six fixed methods are reported; the new test set was not used for a search over methods or thresholds.','',
        'The19 previously examined marker genes form a separate development/interpretation set. **All39 genes are removed together from the matching features and even the library-size normalization denominators.** Target values from RNA are used only after matching to form predictions. Native CytoSPACE, nearest-one, plain soft30, assay-scaled soft30 and its spatial/shuffled controls use the same held-out query cells and RNA pool. These comparisons distinguish predictive gain from agreement on genes supplied to the match. [Frozen design](benchmark_design.json).','',
        'Native baseline: unmodified upstream CytoSPACE1.1.0, seed42, single-cell mode, a single supplied `Ventricular_CM` label,3,000-cell partitions, two workers and molecule downsampling off, with the same corrected-reference values as the alternative methods. This is a controlled independent benchmark, not an exact author-pipeline reproduction. Commands, source hashes and training-gene lists are retained per stage.','',
        'The cohort is the preceding independently gated, QC/information-pass ventricular-enriched set:21,994 /38,580 /49,656 /58,358 cells, totaling168,588. Published spatial labels were not used to fit or select this mapping. The previously frozen gate/classifier and QC used the original full panel, however, so these are **conditional gene-prediction tests within an existing ventricular cohort**, not completely blinded end-to-end cell annotation. No claim is made that all possible gate dependence was removed.','',
        'Evaluation additionally reports detection AUC, Pearson correlation, training-depth/cell-area/training-detected-gene-adjusted correlation and fixed-grid regional correlation. The adjustment covariates use training genes only. Cell and regional correlations measure different resolutions. Regional averaging can improve correlations by reducing noise; it is not proof of cell-level precision.','',
        '## RNA-library effects and subtype confidence','',
        '[Technical library comparison](cross_library_technical_comparison.png) · [per-label metrics](cross_library_technical_comparison.tsv) · [per-nucleus predictions](cross_library_predictions.tsv.gz) · [agreement versus weight threshold](cross_library_selective_agreement.tsv)','',
        'Each age is analyzed separately. RNA queries are tested against the opposite library, with no shared query/reference nuclei. Production maps also compare separate-library projections. This addresses sensitivity to the two available RNA libraries; it cannot establish independent-animal replication or remove stage/section confounding. RNA labels were previously defined using the combined data, so technical label recovery is not independent biological identity validation. The weight fractions become more selective at higher cutoffs in this technical test, but they are not calibrated Xenium probabilities.','',
        mdtable(cov),'',
        '**These flags are different from the previous stringent CytoSPACE panel/seed/raw/library/marker acceptance rules.** Their counts must not be treated as an increase in independently validated subtypes or a biological abundance estimate. Both proposed labels and the numerical evidence are retained so that weak or rare distinctions can be evaluated separately. Conduction-program candidates remain RNA interpretations, not anatomically verified Purkinje cells.','',
        '## Maps and downloadable data','',
        '| Stage | Spatial maps | All-cell annotations | Weighted reference neighbors |',
        '| --- | --- | --- | --- |']
    for s in STAGES:
        lines.append(f'| {s} | [PNG]({s}_weighted_family_maps.png) · [PDF]({s}_weighted_family_maps.pdf) | [TSV.gz]({s}/weighted_CM_annotations.tsv.gz) | [NPZ]({s}/weighted_reference_neighbors.npz) |')
    lines += ['',
        'Maps retain all original section coordinates and include unassigned/excluded cells in gray. Left panels show the highest-weight family for all included cells; right panels show candidates meeting the stated evidence rules. Every original segment remains in the annotation tables, including exclusions. `weighted_family`/`weighted_fine_group` are proposals; `family_candidate`/`fine_group_candidate` add the evidence flags. `neighbor_weight_family__*` columns expose the full family weight distribution. Cell IDs are section-specific.','',
        'Each neighbor NPZ contains query and reference barcode arrays, query×30 neighbor indices/weights, label classes and family/fine weight fractions. Reference rows index the existing [RNA metadata](../xenium_ventricular_mapping/reference_metadata.tsv.gz). Prediction NPZs retain all six methods×cells×39 held-out expression predictions and measured normalized targets, with explicit gene/method/barcode arrays. No imputed values replace original Xenium counts.','',
        '## What this supports','',
        'The same datasets can support stronger **internal validation**: withhold measured genes, separate RNA libraries, compare simple baselines, test spatial controls and expose uncertainty. The improvement here is in reconstructing measured-gene variation and obtaining graded reference-program evidence. An exact RNA nucleus–Xenium cell pairing is generally not identifiable from these data. Stronger fine-subtype anatomy still benefits from independent histology, additional markers or sections.','',
        '**Kcnd2 and Irx5 remain unassayed.** Improved recovery of measured genes makes transfer more defensible but does not validate an unmeasured Kcnd2 gradient. No new Kcnd2 spatial map or guaranteed wall-depth coordinate is inferred here. The [actual published author-state maps](../author_atlas_baseline/README.md) remain a separate comparison layer.','',
        '## Reproduce and audit','',
        '```bash',
        'OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python analysis/improve_ventricular_alignment.py',
        'OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python analysis/audit_ventricular_alignment_improvement.py',
        'python analysis/report_ventricular_alignment_improvement.py','```','',
        'Use the scientific Python and separate pinned CytoSPACE runtime described in the [previous report](../cytospace_xenium/README.md). Completion files skip existing stage fits; remove/rebuild the relevant derived stage outputs if deliberately changing inputs or parameters. The audit recomputes metrics and evidence flags from saved predictions and checks source hashes, all cell/coordinate joins, native IDs/reuse, neighbor bounds/weights and independently sampled projections. [Numerical validation](validation.json), [methods](methods.json) and [manifest/checksums](manifest.json).','',
        'Sources: GSE290989; Xenium Zenodo [P0/P7](https://zenodo.org/records/14961929) and [P14/P21](https://zenodo.org/records/14961942), Wang, Dong, Song, Qian and Liu, CC BY4.0; [article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/). Figures are newly generated from deposited data. CytoSPACE: Vahid et al., Nature Biotechnology2023, [DOI](https://doi.org/10.1038/s41587-023-01697-9), [source/license](https://github.com/digitalcytometry/cytospace), pinned commit`bb7fe12ab4466ef5e38dcac4c3ece07012060d60`.','']
    prose='\n'.join(lines)
    for before,after in {'at least20':'at least 20','from30':'from 30','plain30-neighbor':'plain 30-neighbor',
        'floor0.25':'floor 0.25','at±5':'at ±5','uses30':'uses 30','temperature0.10':'temperature 0.10',
        'requires≥0.70':'requires ≥0.70','evaluated19':'evaluated 19','These20':'These 20','The19':'The 19',
        'All39':'All 39','CytoSPACE1.1.0':'CytoSPACE 1.1.0','seed42':'seed 42','label,3,000':'label, 3,000',
        'set:21,994 /38,580 /49,656 /58,358':'set: 21,994 / 38,580 / 49,656 / 58,358',
        'totaling168,588':'totaling 168,588','CC BY4.0':'CC BY 4.0','Biotechnology2023':'Biotechnology 2023',
        'commit`':'commit `','All381':'All 381'}.items():prose=prose.replace(before,after)
    (OUT/'README.md').write_text(prose)
    methods={'python':platform.python_version(),'versions':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','anndata','scikit-learn','matplotlib']},
        'design':design,'primary_projection':'assay_scaled_soft30, unsmoothed',
        'raw_native_baseline_is_exact_author_pipeline':False,'new_test_genes_used_to_tune':False,
        'all39_genes_excluded_from_normalization_denominators':True,
        'primary_production_gene_count':381,'source_expression_counts_modified':False,
        'source_coordinates_modified':False,'animal_replicated_validation':False,
        'weight_fractions_are_calibrated_probabilities':False,'Kcnd2_Irx5_measured':False}
    (OUT/'methods.json').write_text(json.dumps(methods,indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    scripts=[ROOT/'analysis'/n for n in ['improve_ventricular_alignment.py','audit_ventricular_alignment_improvement.py','report_ventricular_alignment_improvement.py']]
    (OUT/'manifest.json').write_text(json.dumps({'results':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(files)],
        'analysis':[{'path':str(p.relative_to(ROOT)),'sha256':digest(p)} for p in scripts],'validation':validation},indent=2)+'\n')
    print('REPORT_READY',len(files),'files',flush=True)


if __name__=='__main__':main()
