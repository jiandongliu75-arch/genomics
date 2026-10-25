"""Report legitimate detection gains, unresolved specificity and spatial proxies."""
from vascular_detection_common import *
from enhance_vascular_Xenium_detection import SPATIAL_PANELS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import textwrap

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})

def finish(fig,name,pdf):
    fig.savefig(DEST/(name+'.png'),dpi=170,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def textpage(pdf,title,paragraphs):
    fig=plt.figure(figsize=(12,8.5));fig.text(.05,.94,title,fontsize=19,weight='bold',va='top');y=.87
    for p in paragraphs:
        for line in textwrap.wrap(p,110):fig.text(.05,y,line,fontsize=11,va='top');y-=.027
        y-=.02
    assert y>=.025,('Page text overflow',title,y)
    pdf.savefig(fig);plt.close(fig)

def main():
    rna=pd.read_csv(DEST/'RNA/fractional_count_and_background_audit.tsv',sep='\t')
    support=pd.read_csv(DEST/'RNA/gene_support_across_two_libraries.tsv',sep='\t')
    cycle=pd.read_csv(DEST/'RNA/cycling_fractional_threshold_sensitivity.tsv',sep='\t')
    benchmark=pd.read_csv(DEST/'RNA/withheld_program_opposite_library_benchmark.tsv',sep='\t')
    rscore=pd.read_csv(DEST/'RNA/rank_score_library_summary.tsv',sep='\t')
    coherence=pd.read_csv(DEST/'RNA/within_library_gene_coherence.tsv',sep='\t')
    coverage=[];spatial_coherence=[];pool_coverage=[];consistency=[];cover=[]
    for stage in STAGES:
        folder=DEST/'Xenium'/stage
        coverage.append(pd.read_csv(folder/'observed_detection_coverage.tsv',sep='\t'))
        spatial_coherence.append(pd.read_csv(folder/'observed_gene_coherence.tsv',sep='\t'))
        pool_coverage.append(pd.read_csv(folder/'adaptive_pool_detection_coverage.tsv',sep='\t'))
        consistency.append(pd.read_csv(folder/'cross_view_proxy_consistency.tsv',sep='\t'))
        if stage=='P0':cover=pd.read_csv(folder/'panel_coverage.tsv',sep='\t')
    coverage=pd.concat(coverage);pool_coverage=pd.concat(pool_coverage);spatial_coherence=pd.concat(spatial_coherence);consistency=pd.concat(consistency)
    save(coverage,DEST/'all_stage_spatial_detection_coverage.tsv')
    save(pool_coverage,DEST/'all_stage_adaptive_pool_coverage.tsv')
    save(spatial_coherence,DEST/'all_stage_spatial_gene_coherence.tsv')
    save(consistency,DEST/'all_stage_cross_view_consistency.tsv')
    save(cover,DEST/'RNA_Xenium_program_coverage.tsv')
    cg=cycle.groupby(['stage','layer','minimum_expected_count'],observed=True)[['nuclei','cycling_support_nuclei']].sum().reset_index()
    cg['fraction']=cg.cycling_support_nuclei/cg.nuclei
    save(cg,DEST/'cycling_detection_rule_comparison.tsv')
    # The main spatial display deliberately uses measured genes only.
    with PdfPages(DEST/'vascular_detection_enhancement_report.pdf') as pdf:
        textpage(pdf,'Enhanced vascular-program detection · findings',[
            'An important correction: applying an integer threshold to fractional DecontX counts understated retained low-count RNA. P7 Apln has56 raw molecules and55.6 corrected expected molecules. The old >=1 rule counted5 nuclei; raw-supported fractional thresholds recover47. This is a quantification correction, not restoration of removed ambient RNA.',
            'RNA enhancement: preserve fractional counts; report expected-mass pseudobulks and threshold sensitivity; use expression/detection-matched rank enrichment; require age-matched library support; audit EC versus non-EC raw RNA and gene coherence. A cellular-background estimate is only a sensitivity proxy because no empty-droplet background is available.',
            'The earlier cycling percentages also used the strict integer rule. A fractional-aware >=0.5 expected-count threshold plus raw observation gives71/1249 atP0,71/833 atP7,20/1242 atP14 and3/1273 atP21. Thresholds0.1/0.25/0.5 agree here. These are marker-support fractions in the frozen QC cohort, not independently calibrated biological proliferation rates.',
            'Xenium directly measures onlyAngpt2 from the original six-gene sprouting panel. Apln/Esm1/Cxcr4/Pgf/Adm are unmeasured. Observed guidance, matrix and cycling proxies are kept separate. Neither lower RNA-quality cutoffs nor invented counts add these missing genes.',
            'Spatial enhancement: frozen nuclear-clean ECs, nuclear/core/native and matched-depth controls, non-overlapping25/50/100um tiles, and five-EC pools within75um. The latter cover89-99% of the selected EC positions but overlap, so they are displays at neighborhood resolution, not independent observations or cell-specific states.'])
        fig,axes=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
        for ax,g in zip(axes,['Apln','Angpt2']):
            a=rna[rna.gene.eq(g)].groupby('stage')[['raw_molecules','corrected_expected_mass']].sum().reindex(STAGES)
            ax.bar(np.arange(4)-.18,a.raw_molecules,width=.36,label='Raw molecules',color='#777777')
            ax.bar(np.arange(4)+.18,a.corrected_expected_mass,width=.36,label='Corrected expected mass',color='#0072b2')
            ax.set_xticks(range(4),STAGES);ax.set_title(g);ax.set_ylabel('Summed RNA counts / expected counts');ax.legend(frameon=False,fontsize=9)
        fig.suptitle('Most Apln and Angpt2 RNA mass was retained\nA strict >=1 rule missed near-one fractional values',fontsize=15)
        finish(fig,'RNA_raw_vs_retained_signal_mass',pdf)
        fig,axes=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
        for ax,g in zip(axes,['Apln','Angpt2']):
            a=rna[rna.gene.eq(g)].groupby('stage')[['nuclei','raw_positive_nuclei','corrected_positive_ge1','corrected_positive_ge0p5']].sum().reindex(STAGES)
            for j,(col,label,c) in enumerate([('raw_positive_nuclei','Raw >=1','#777777'),('corrected_positive_ge1','Legacy corrected >=1','#d55e00'),('corrected_positive_ge0p5','Corrected >=0.5','#0072b2')]):
                ax.bar(np.arange(4)+(j-1)*.25,100*a[col]/a.nuclei,width=.25,label=label,color=c)
            ax.set_xticks(range(4),STAGES);ax.set_title(g);ax.set_ylabel('Nuclei meeting the specified rule (%)');ax.legend(frameon=False,fontsize=8)
        fig.suptitle('Detection-rule sensitivity · raw observation required\nThese are observed/fractional RNA-support criteria, not calibrated cell-state probabilities',fontsize=14)
        finish(fig,'RNA_fractional_detection_rule_recovery',pdf)
        fig,axes=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
        for layer,threshold,label,c in [('DecontX',1,'Legacy corrected >=1','#d55e00'),('DecontX',.5,'Fractional-aware >=0.5 + raw evidence','#0072b2')]:
            a=cg[cg.layer.eq(layer)&cg.minimum_expected_count.eq(threshold)].set_index('stage').reindex(STAGES)
            axes[0].plot(STAGES,100*a.fraction,marker='o',color=c,label=label)
        axes[0].set_ylabel('Cycling-marker support (%)');axes[0].legend(fontsize=8,frameon=False);axes[0].set_title('Earlier cycling fractions need qualification')
        for i,p in enumerate(PROGRAMS):
            a=coherence[coherence.program.eq(p)&coherence.layer.eq('DecontX')]
            axes[1].scatter(np.full(len(a),i),a.median_leave_one_gene_out_partial_Spearman,s=30,alpha=.7)
        axes[1].set_xticks(range(4),['Original\ncore','Guidance\nsupport','Matrix\nremodeling','Cycling']);axes[1].axhline(0,color='gray',lw=.8)
        axes[1].set_ylabel('Within-library partial gene-coherence correlation');axes[1].set_title('Coverage does not establish coherent activation')
        fig.suptitle('Fractional-aware detection preserves the early-growth conclusion',fontsize=15)
        finish(fig,'cycling_correction_and_program_coherence',pdf)
        fig,axes=plt.subplots(1,4,figsize=(17,4.8),layout='constrained')
        for ax,program in zip(axes,PROGRAMS):
            for layer,c in [('raw','#999999'),('DecontX','#0072b2')]:
                a=rscore[rscore.layer.eq(layer)&rscore.maxrank.eq(1500)&rscore.program.eq(program)]
                means=[]
                for i,s in enumerate(STAGES):
                    v=a.loc[a.stage.eq(s),'mean_rank_excess'].to_numpy();means.append(v.mean());ax.scatter([i-.05,i+.05],v,s=22,color=c)
                ax.plot(range(4),means,color=c,label=layer)
            ax.set_xticks(range(4),STAGES);ax.set_title(program.replace('_',' '),fontsize=10);ax.set_ylabel('Mean matched-control rank excess');ax.legend(frameon=False,fontsize=8)
        fig.suptitle('Fractional-aware rank scoring · target-matched controls held fixed across stages',fontsize=14)
        finish(fig,'RNA_fractional_aware_rank_programs',pdf)
        fig,axes=plt.subplots(1,2,figsize=(14,5.5),layout='constrained')
        c=cover.groupby('program').agg(total=('gene','size'),observed=('Xenium_measured','sum')).reindex(PROGRAMS)
        axes[0].barh(range(4),c.total,color='#dddddd',label='Full RNA panel');axes[0].barh(range(4),c.observed,color='#0072b2',label='Measured on Xenium')
        axes[0].set_yticks(range(4),[p.replace('_',' ') for p in PROGRAMS]);axes[0].set_xlabel('Genes');axes[0].legend(frameon=False,fontsize=9)
        for i,p in enumerate(PROGRAMS):
            d=benchmark[benchmark.program.eq(p)];axes[1].scatter(np.full(len(d),i),d.Spearman,s=25)
        axes[1].axhline(.2,color='#d55e00',ls='--',label='One of three required gates')
        axes[1].set_xticks(range(4),['Original\ncore','Guidance\nsupport','Matrix','Cycling']);axes[1].set_ylabel('Opposite-library withheld-program Spearman');axes[1].legend(fontsize=8,frameon=False)
        fig.suptitle('Missing targets cannot be treated as measured RNA\nNo program passed all predefined inference gates in this benchmark',fontsize=14)
        finish(fig,'Xenium_coverage_and_withheld_inference_limits',pdf)
        titles=['Angpt2 · one gene','Guidance-marker composite','Matrix-marker composite','Cycling markers']
        for stage in STAGES:
            d=pd.read_csv(DEST/'Xenium'/stage/'adaptive_5EC_75um_nuclear_proxy_maps.tsv.gz',sep='\t',index_col=0)
            source=pd.read_csv(OUT/'Xenium'/stage/'EC_consensus_annotations.tsv.gz',sep='\t',index_col=0)
            fig,axes=plt.subplots(1,4,figsize=(20,7),layout='constrained')
            for ax,(program,panel),title in zip(axes,SPATIAL_PANELS.items(),titles):
                ax.scatter(source.x_centroid,source.y_centroid,s=.2,color='#dddddd',rasterized=True)
                supported=d[program+'__observed_genes'].ge(1 if len(panel)==1 else 2)
                score=d[program+'__matched_control_excess'];limit=max(.1,float(np.nanquantile(np.abs(score),.98)))
                p=ax.scatter(d.loc[supported,'x_centroid'],d.loc[supported,'y_centroid'],c=score[supported],s=2.2,cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
                ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(title,fontsize=11);ax.set(xlabel='X (µm)',ylabel='Y (µm)')
                fig.colorbar(p,ax=ax,shrink=.45,label='Matched-control RNA excess')
            fig.suptitle(f'{stage} · observed nuclear RNA in five-EC neighborhoods within75µm\nOverlapping display pools; guidance/matrix composites are exploratory proxies, not validated sprouting activity',fontsize=14)
            finish(fig,stage+'_observed_EC_neighborhood_proxy_maps',pdf)
        textpage(pdf,'Interpretation and remaining limits',[
            'Detection has improved without manufacturing RNA. Corrected expected mass and fractional-aware scoring recover low-count EC support that the earlier integer cutoff obscured. Raw and corrected matrices are unchanged; lowering a detection threshold does not identify true biology by itself.',
            'Apln and Angpt2 have support in both libraries at several stages. Esm1 remains genuinely sparse in this RNA dataset; it is absent in raw P7 EC counts. The full six-gene program does not show uniformly strong within-library coherence. Expanded guidance and remodeling panels are supportive functional composites, not replacements for the original biological definition.',
            'Spatial pooling increases the chance of observing multiple genes, at the cost of resolution. The five-EC displays overlap and can combine neighboring vessels; they must not be counted as independent active cells. Non-overlapping tiles, common source profiles and molecule-depth controls supply the more conservative sensitivity/coherence comparisons.',
            'OnlyAngpt2 from the original six-gene panel is directly assayed spatially. This specific cross-library Ridge benchmark does not justify imputing missing program genes; it does not prove every possible inference method must fail. Direct Apln/Esm1 validation would require additional measured targets, such as targeted RNA in situ with EC and membrane/nuclear markers.',
            'Non-EC raw cellular RNA is not measured ambient soup. Candidate overattenuation flags therefore cannot authorize restoring raw counts. The five vascular EC identity labels and the original reference cohort remain fixed; formerly reviewed ECs are evaluated only in a separate fractional-aware sensitivity table.',
            'Two RNA libraries per age are not assumed independent animals. One Xenium section per age does not support organ-wide developmental rates or causal vessel growth. Improved coverage is not independently validated program accuracy.'])
    # Reviewable interpretation, including an explicit correction to earlier rates.
    lines=['| Stage | Legacy ≥1 rule | Fractional-aware ≥0.5 + raw evidence |','|---|---:|---:|']
    for s in STAGES:
        a=cg[cg.stage.eq(s)&cg.layer.eq('DecontX')].set_index('minimum_expected_count')
        lines.append(f'| {s} | {100*a.loc[1,"fraction"]:.2f}% | {100*a.loc[.5,"fraction"]:.2f}% |')
    apln=support[support.gene.eq('Apln')].set_index('stage').reindex(STAGES)
    reads=['| Stage | Apln raw molecules | Corrected expected mass | Legacy positive nuclei | Fractional-aware support ≥0.25 | Both-library pseudobulk support |',
           '|---|---:|---:|---:|---:|---|']
    for s,d in apln.iterrows():reads.append(f'| {s} | {int(d.raw_molecules_total)} | {d.corrected_expected_mass_total:.2f} | {int(d.corrected_positive_ge1_total)} | {int(d.corrected_positive_ge0p25_total)} | {d.pseudobulk_support_both_libraries} |')
    (DEST/'README.md').write_text('''# Improved vascular-program detection: fractional RNA, pooling and validation

The main correctable problem was **using an integer detection threshold on fractional DecontX estimates**, rather than loss of most sprouting RNA. P7 Apln retains55.6 of56 raw molecules, but the legacy>=1 corrected-count rule counted only5 nuclei. Fractional-aware support recovers47 nuclei. The correction is to quantification and interpretation; no raw molecules are restored and original count matrices are unchanged.

[Illustrated report](vascular_detection_enhancement_report.pdf) · [Fractional/gene/background audit](RNA/fractional_count_and_background_audit.tsv) · [Both-library support](RNA/gene_support_across_two_libraries.tsv) · [Spatial panel coverage](RNA_Xenium_program_coverage.tsv) · [Adaptive pool coverage](all_stage_adaptive_pool_coverage.tsv) · [Validation](validation.json)

## Corrected interpretation of low-count RNA

'''+ '\n'.join(reads)+'''

Fractional-aware expected-count thresholds0.1,0.25,0.5 and1 are all reported. The alternative support rule always requires an original raw observation; positive numerical dust alone does not qualify. Both-library pseudobulk support means corrected expected mass>=5 plus>=3 nuclei with expected count>=0.25 in each same-age library; it is an exploratory criterion, not independently calibrated accuracy. The full corrected mass and original observations are more informative than one binary threshold.

**Earlier cycling percentages also need qualification.** The old>=1 corrected-count rule was overly conservative for this fractional matrix. Raw-observed support at corrected expected count>=0.5 gives:

'''+ '\n'.join(lines)+'''

For cycling, thresholds0.1/0.25/0.5 give identical counts here and agree with the original raw two-marker rule. This supports an early-to-late decline and a higher observed P7 fraction than P0 in this frozen cohort; it does not establish a population-level biological peak or a histological proliferation rate. The earlier3.04/2.88/1.13/0.08% values remain reproducible **strict-threshold counts**, not the best current prevalence estimates. Existing five-type identities and prior QC cohorts remain unchanged.

![RNA mass retained](RNA_raw_vs_retained_signal_mass.png)

## RNA enhancement and background audit

The 4,597 primary EC nuclei are rescored without overwriting counts. Library-level expected-mass pseudobulks reduce dropout; fractional thresholds expose cutoff effects. Top1500 rank enrichment subtracts20 expression/detection-matched control genes per target; top3000 is a sensitivity. Controls are selected using equal-weight library raw expression and detection, excluding programs, lineage diagnostics, mitochondrial/ribosomal and sex genes. All stages share the same controls. A0.001 numerical-dust floor avoids treating arbitrarily tiny residuals as expressed genes. These are diagnostic enrichment scores, not pathway activation probabilities.

The original six-gene panel is kept separate from an expanded sprouting/guidance-support panel, an EC matrix-remodeling panel and a cycling control. New panels include multifunctional and inhibitory genes and do not redefine the original tip-cell biology. `RNA/within_library_gene_coherence.tsv` reports leave-one-gene-out correlations after log depth, ambient fraction and source EC-family adjustment. Cycling is comparatively coherent. Guidance/remodeling composites have weak single-nucleus coherence; more covered genes must not be described as established coordinated activation.

The whole RNA reference is independently audited for target-gene abundance in each lineage/library. Non-EC cellular RNA supplies a background **sensitivity proxy**, not measured empty-droplet soup. The output compares raw EC enrichment, retained expected mass, expected cellular background and a2-times-background residual sensitivity. Such evidence can flag candidate overattenuation but cannot validate a new ambient model or justify restoring counts. For P7 Apln/Angpt2, the mass-retention audit supports correcting the binary rule rather than inferring major decontamination loss.

`RNA/former_QC_review_EC_sensitivity.tsv` examines previously reviewed ECs using raw-supported fractional lineage markers with residual foreign-program and coherent mural-group exclusions. This is a secondary sensitivity cohort; no nucleus or subtype label is silently added to the frozen primary reference. Original all-EC capture, QC and library biases remain.

## Spatial enhancement with measured RNA

The original six-gene sprouting panel has only **1/6 targets on Xenium: Angpt2**. Apln, Esm1, Cxcr4, Pgf and Adm are not measured. We report the following distinct observed views:

| Spatial view | Directly measured genes | Interpretation |
|---|---|---|
| Angpt2 | Angpt2 | Single-gene evidence, not a six-gene program |
| Guidance-associated composite | Angpt2, Unc5b, Plxnd1, Itgb3, Lama4 | Exploratory guidance/remodeling proxy; weak gene coherence |
| EC matrix-associated composite | Mmp2, Itga5, Itgb3, Fn1 | Exploratory remodeling proxy; lineage/boundary sensitive and shared with other cells |
| Cycling | Mki67, Top2a, Cenpf | Proliferation-associated RNA, separate from tip/stalk identity |

The primary spatial cohort is fixed before outcome scoring: original broad EC identity plus nuclear-broad-clean evidence, excluding multiple projected nuclear components and unresolved geometry. It contains29,363 profiles across four sections. Native, refined whole-cell, boundary-core, nuclear and exact molecule-depth controls all use identical source profiles. No QV threshold is lowered, previously unassigned RNA added or missing gene imputed.

Non-overlapping25/50/100µm tiles require>=5 selected EC profiles and preserve RNA counts. The50µm cohort can be spatially selective because sparse regions fail the minimum. Five-nearest-EC pools within75µm provide a higher-coverage display:89-99% of the **selected** EC positions, not all heart ECs. Pools overlap and can span neighboring vessels. They are neighborhood expression estimates, never independent observations or individual-cell functional labels. Source cell IDs and maximum distances are supplied for each pool. The stricter cell-and-nuclear-clean cohort is a separate50µm sensitivity.

Pooling increases multi-gene detection at the cost of resolution. Profile-versus-pool coverage uses different units and is not a validated sensitivity ratio. Coherence and cross-view comparisons use non-overlapping tiles, with depth, EC count and CM-RNA fraction diagnostics. Displayed guidance/matrix scores are observed marker composites with matched-control excess. They do **not** establish coherent sprouting programs, especially where held-out gene correlations are weak or inconsistent. P21 matrix coherence is particularly uncertain. Color scales are descriptive within-section scores, not comparable rates across animals/ages.

![P7 observed spatial proxies](P7_observed_EC_neighborhood_proxy_maps.png)

## Missing-gene inference gate

All program genes are withheld from the364 shared RNA fitting genes. Same-age opposite-library Ridge regression (alpha10; SD floor0.25; clipping±5) predicts corrected rank scores. Required gates are Spearman>=0.2, within-family residual Spearman>=0.15 and positive query R² in both directions. No program passes all gates in this benchmark, so missing-gene inferred program maps are **not released**. This is one conservative technical test, not proof that every possible inference method must fail. Mapping broad cell identity is different from recovering a rare unmeasured state program.

## Result and remaining validation

We recovered legitimate low-count RNA evidence and increased spatial neighborhood coverage. We did not establish a universally coherent six-gene sprouting program or manufacture missing spatial measurements. Esm1 remains genuinely sparse, including zero raw P7 EC molecules. Direct Apln/Esm1 spatial validation would need additional measured targets, for example targeted RNA in situ combined with EC, membrane and nuclear markers. Proliferation is still the strongest early-growth signal; guidance and ECM composites retain their exploratory status.

There are two RNA libraries per age and one Xenium section per age. These are not assumed independent animal replicates, and no organ-wide rate, causal vascular-growth or independently calibrated detection-accuracy claim is made. The five requested EC identity categories are preserved.

Reproduce in order with `analysis/enhance_vascular_RNA_detection.py`, `analysis/enhance_vascular_Xenium_detection.py --stage P0` (andP7/P14/P21), `analysis/report_vascular_detection.py`, then `analysis/validate_vascular_detection.py`. Shared definitions are in `analysis/vascular_detection_common.py`. Earlier source matrices and QC remain under `results/endothelial_QC_subtypes/`. `validation.json` and `file_manifest.tsv` record data-level checks and SHA256 hashes.
''')
    print('VASCULAR_REPORT_COMPLETE',flush=True)

if __name__=='__main__':main()
