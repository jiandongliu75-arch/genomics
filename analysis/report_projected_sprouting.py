"""Figures and documentation for an explicitly inferred sprouting layer."""
from project_RNA_sprouting_to_Xenium import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import textwrap

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})

def finish(fig,name,pdf):
    fig.savefig(PROJ/(name+'.png'),dpi=180,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def textpage(pdf,title,paragraphs):
    fig=plt.figure(figsize=(12,8.5));fig.text(.05,.94,title,fontsize=19,weight='bold',va='top');y=.87
    for p in paragraphs:
        for line in textwrap.wrap(p,110):fig.text(.05,y,line,fontsize=11,va='top');y-=.027
        y-=.019
    assert y>=.02,(title,y)
    pdf.savefig(fig);plt.close(fig)

def main():
    benches=[];supports=[];counts=[];sensitivity=[];whole=[]
    for stage in STAGES:
        folder=PROJ/stage
        benches.append(pd.read_csv(folder/'opposite_library_projection_benchmark.tsv',sep='\t'))
        supports.append(pd.read_csv(folder/'projection_support_gates.tsv',sep='\t'))
        counts.append(json.loads((folder/'complete.json').read_text()))
        whole.append(pd.read_csv(folder/'all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz',sep='\t',index_col=0))
        a=pd.read_csv(folder/'nuclear_inferred_sprouting_annotations.tsv.gz',sep='\t',index_col=0)
        for view in ['native','boundary_core']:
            b=pd.read_csv(folder/f'{view}_inferred_sprouting_annotations.tsv.gz',sep='\t',index_col=0)
            for name in ['sprouting_6g','unmeasured_5g']:
                col='corrected__RNA_resolution1__inferred_'+name
                sensitivity.append({'stage':stage,'outcome':name,'comparison':view+'_vs_nuclear','Spearman':corr(a[col],b[col]),
                    'median_absolute_difference':float(np.median(np.abs(a[col]-b[col])))})
        for name in ['sprouting_6g','unmeasured_5g']:
            base='corrected__RNA_resolution1__inferred_'+name
            for col in ['raw__RNA_resolution1__inferred_'+name,'neighbors10__inferred_'+name,'neighbors60__inferred_'+name,'Angpt2_withheld__inferred_'+name]:
                sensitivity.append({'stage':stage,'outcome':name,'comparison':col,'Spearman':corr(a[base],a[col]),
                    'median_absolute_difference':float(np.median(np.abs(a[base]-a[col])))})
    bench=pd.concat(benches);support=pd.concat(supports)
    save(bench,PROJ/'all_stage_opposite_library_benchmarks.tsv')
    save(support,PROJ/'unmeasured_component_support_by_resolution.tsv')
    save(pd.DataFrame(sensitivity),PROJ/'projection_sensitivity_summary.tsv')
    allheart=pd.concat(whole);assert allheart.stage_cell_id.is_unique
    allheart.to_csv(PROJ/'all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz',sep='\t',na_rep='NA')
    row_counts=pd.DataFrame([{k:c[k] for k in ['stage','source_EC_candidates','projected_ECs','projected_spatial_neighborhoods5','projected_spatial_neighborhoods20','whole_heart_profiles']} for c in counts])
    save(row_counts,PROJ/'projection_coverage.tsv')
    primary=bench[bench.method.eq('assay_scaled_soft')&bench.neighbors.eq(30)]
    summary=primary.groupby(['stage','RNA_resolution_nuclei','outcome'],observed=True).agg(
        libraries=('query_library','nunique'),min_Spearman=('Spearman','min'),max_Spearman=('Spearman','max'),
        min_within_family_Spearman=('within_family_Spearman','min'),max_within_family_Spearman=('within_family_Spearman','max'),
        min_R2=('R2_vs_query_mean','min'),max_R2=('R2_vs_query_mean','max'),both_libraries_gate=('technical_support_gate','all')).reset_index()
    save(summary,PROJ/'primary_benchmark_summary.tsv')
    with PdfPages(PROJ/'RNA_inferred_sprouting_Xenium_report.pdf') as pdf:
        textpage(pdf,'RNA-inferred sprouting added to the Xenium atlas',[
            'The new layer is a reference-derived estimate of sprouting-associated RNA, not an additional Xenium measurement. It annotates29,363 nuclear-clean previously labeled ECs across P0/P7/P14/P21. Every source cardiac profile remains present in the whole-heart annotation table; non-EC and excluded EC positions have no inferred score.',
            'The RNA reference contains4,597 QC-selected endothelial nuclei. Same-age reference neighbors are selected with the shared measured genes. Their corrected RNA expression is transferred as a weighted expectation:30 cosine neighbors,temperature0.1, separate-assay scaling,SDfloor0.25 and clipping5. No sprouting label, desired score or physical coordinate selects neighbors.',
            'Two scores are supplied. The six-gene score averages inferred Apln,Esm1,Angpt2,Cxcr4,Pgf andAdm. The five-gene score excludes Angpt2 and exposes the genuinely unmeasured component. Measured Angpt2 counts remain separately available. Missing genes are never inserted into the observed403-gene count matrix.',
            'Individual scores,5-EC neighborhoods within75um and20-EC neighborhoods within150um are supplied. Neighborhoods are projected to disjoint, shared-gene-defined RNA metacells at matched target resolution. Pooled scores use a different reference resolution and should not be compared numerically with single-cell scores as if they were the same estimator.',
            'All spatial output is explicitly exploratory. The opposite-library RNA test supports parts of the total score, but recovery of the five unmeasured genes is weak at individual-cell resolution. Coarse neighborhood support is limited and age-dependent; it is not independent spatial validation or proof of new-vessel growth.'])
        fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
        for ax,outcome,title in zip(axes,['sprouting_6g','unmeasured_5g','cycling_heldout_control'],
                                   ['Six-gene total (Angpt2 allowed)','Five genuinely unmeasured genes','Cycling: target genes withheld']):
            d=primary[primary.outcome.eq(outcome)]
            for size,c,offset in [(1,'#0072b2',-.18),(5,'#e69f00',0),(20,'#009e73',.18)]:
                for i,stage in enumerate(STAGES):
                    v=d.loc[d.stage.eq(stage)&d.RNA_resolution_nuclei.eq(size),'Spearman']
                    ax.scatter(np.repeat(i+offset,len(v)),v,s=30,color=c,label=f'RNA groups of{size}' if i==0 else None)
            ax.axhline(.2,c='gray',ls='--',lw=1);ax.set_xticks(range(4),STAGES);ax.set_ylabel('Opposite-library Spearman');ax.set_title(title,fontsize=11)
            ax.legend(fontsize=8,frameon=False)
        fig.suptitle('Reference-transfer benchmark · both RNA-library directions\nCorrelation alone is insufficient: family-residual correlation and positive R² are also required',fontsize=14)
        finish(fig,'RNA_projection_benchmark',pdf)
        # Consistent display limits across the four sections for each annotation.
        ec=[pd.read_csv(PROJ/s/'EC_annotations_with_inferred_sprouting.tsv.gz',sep='\t',index_col=0) for s in STAGES]
        plotcols=['inferred_sprouting_score_6g','inferred_unmeasured_sprouting_score_5g','measured_Angpt2_native_molecules','inferred_sprouting_RNA_library_difference']
        titles=['Inferred six-gene expectation','Inferred unmeasured-five component','Observed native Angpt2','RNA-library projection difference']
        limits={c:max(float(np.nanquantile(np.concatenate([x.loc[x.sprouting_inference_eligible,c].to_numpy() for x in ec]),.99)),.001) for c in plotcols}
        for stage,a in zip(STAGES,ec):
            w=allheart[allheart.stage.eq(stage)];d=a[a.sprouting_inference_eligible]
            fig,axes=plt.subplots(1,4,figsize=(20,5.5),layout='constrained')
            for ax,col,title in zip(axes,plotcols,titles):
                ax.scatter(w.x_centroid,w.y_centroid,s=.12,color='#e4e4e4',rasterized=True)
                p=ax.scatter(d.x_centroid,d.y_centroid,c=d[col],s=2.5,cmap='viridis',vmin=0,vmax=limits[col],rasterized=True)
                ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(title,fontsize=10);ax.set(xlabel='X (µm)',ylabel='Y (µm)')
                fig.colorbar(p,ax=ax,shrink=.6,label='Observed molecules' if col=='measured_Angpt2_native_molecules' else 'RNA-reference log1pCP10k units')
            fig.suptitle(f'{stage} · RNA-reference sprouting projection onto{len(d):,} Xenium ECs\nExploratory inferred expectations; stronger colors do not establish active sprouting or validate missing-gene recovery',fontsize=14)
            finish(fig,stage+'_inferred_sprouting_and_measured_Angpt2',pdf)
        for stage,a in zip(STAGES,ec):
            w=allheart[allheart.stage.eq(stage)]
            fig,axes=plt.subplots(2,2,figsize=(12,11),layout='constrained')
            for row,size in enumerate([5,20]):
                for col,(key,title) in enumerate([('sprouting_score_6g','Six-gene expectation'),('unmeasured_score_5g','Unmeasured-five component')]):
                    ax=axes[row,col];c=f'inferred_neighborhood{size}_'+key;d=a[a[c].notna()]
                    ax.scatter(w.x_centroid,w.y_centroid,s=.12,color='#e4e4e4',rasterized=True)
                    top=max(float(d[c].quantile(.99)),.001)
                    p=ax.scatter(d.x_centroid,d.y_centroid,c=d[c],s=2,cmap='viridis',vmin=0,vmax=top,rasterized=True)
                    ax.invert_yaxis();ax.set_aspect('equal');ax.set(xlabel='X (µm)',ylabel='Y (µm)');ax.set_title(f'{title} ·{size} ECs',fontsize=11)
                    fig.colorbar(p,ax=ax,shrink=.7,label=f'RNA-metacell{size} normalized-expression units')
            fig.suptitle(f'{stage} · exploratory RNA-inferred neighborhood sprouting\nFive ECs within75µm or20 within150µm; overlapping pools, not independent cells; scales differ by RNA reference resolution',fontsize=13)
            finish(fig,stage+'_inferred_sprouting_neighborhoods',pdf)
        textpage(pdf,'How to use the inferred layer',[
            'Each Xenium_ECs_with_inferred_sprouting.h5ad preserves observed native gene counts inX and403 observed gene features invar. Nuclear and boundary-core counts are retained as observed layers. Projected six-gene expression is stored separately inobsm, in RNA log1pCP10k units. Program scores and uncertainty diagnostics are inobs, with micrometer coordinates inobsm[spatial].',
            'The whole-heart TSV includes all original cardiac positions and identity labels, with inferred values only where EC/nuclear/geometry eligibility was satisfied. Join onstage pluscell_id, or the supplied unique stage_cell_id. Do not interpret missing scores as biological absence or add inferred expression to measured counts.',
            'Raw-reference,10/30/60-neighbor, native/nuclear/boundary-core and Angpt2-withheld sensitivities accompany the primary projection. Reference-neighbor heterogeneity is not a calibrated prediction interval. RNA-library differences and low reference similarity flag possible instability; they do not independently validate a score.',
            'Technical validation uses the opposite RNA library at the same age, with no unmeasured target or cycling-control target among fitting genes. Angpt2 withholding excludes it even from predictor normalization. RNA metacells are disjoint groups formed only from measured features; final small remainders are merged with a nearby group. Spatial neighborhoods overlap and can span neighboring vessels.',
            'Projection can add a useful reference-derived hypothesis layer to Xenium without adding actual molecules. The six-gene score is partly informed by directly measured Angpt2; confidence in the unmeasured five-gene component must be assessed separately. Sparse reference genes and age/library composition limit fine spatial inference. Independent targeted RNA in situ would be needed to verify unmeasured sprouting genes.'])
    lines=['| Stage | ECs with inferred layer | Five-EC positions | Twenty-EC positions |','|---|---:|---:|---:|']
    for d in counts:lines.append(f'| {d["stage"]} | {d["projected_ECs"]:,} | {d["projected_spatial_neighborhoods5"]:,} | {d["projected_spatial_neighborhoods20"]:,} |')
    gates=['| Stage | Individual RNA reference | Five-nucleus RNA groups | Twenty-nucleus RNA groups |','|---|---|---|---|']
    for stage in STAGES:
        s=support[support.stage.eq(stage)].set_index('RNA_resolution_nuclei')
        gates.append('| '+stage+' | '+' | '.join('Pass in both libraries' if bool(s.loc[n,'unmeasured_5g_both_library_gate_pass']) else 'Below gate' for n in [1,5,20])+' |')
    (PROJ/'README.md').write_text('''# RNA-inferred sprouting layer on the Xenium atlas

**Yes: a program absent from the Xenium panel can be added as an explicitly inferred RNA-reference annotation.** This is now implemented for all four stages. It provides an expected sprouting-associated RNA score at eligible Xenium EC positions, separate from observed transcripts. The observed403-gene matrices and five-type EC taxonomy remain unchanged.

[Illustrated report](RNA_inferred_sprouting_Xenium_report.pdf) · [Whole-heart inferred annotation layer](all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz) · [Benchmark summary](primary_benchmark_summary.tsv) · [Projection coverage](projection_coverage.tsv) · [Raw/boundary/reference sensitivity](projection_sensitivity_summary.tsv) · [Validation](validation.json)

## What was added

- `inferred_sprouting_score_6g`: reference-derived expected expression of **Apln, Esm1, Angpt2, Cxcr4, Pgf and Adm**.
- `inferred_unmeasured_sprouting_score_5g`: the five genes **excluding measured Angpt2**, reported separately so measured-gene information cannot masquerade as validated recovery of missing genes.
- `inferred_neighborhood5_*` and `inferred_neighborhood20_*`: measured Xenium RNA pooled over five ECs within75µm or20 ECs within150µm, projected to corresponding RNA metacells. These overlap and are not independent observations or single-cell activity labels.
- Per-gene RNA-reference predictions in `obsm['inferred_sprouting_genes_RNA_log1pCP10k']`, plus measured Angpt2 counts and source/uncertainty metadata. The five unmeasured genes are **not** inserted into observed `var`/`X` or count layers.

'''+ '\n'.join(lines)+'''

There are29,363 projected ECs and408,054 original cardiac positions in the whole-heart annotation layer. Unprojected profiles retain missing scores and their explicit QC/non-EC status. No subtype label or score quota forces a profile into an angiogenic state. Scores condition on the matched-age RNA reference, so between-stage differences can reflect that reference rather than independent spatial validation.

## Reference alignment and target definition

Each stage uses its own previously QC-selected EC RNA reference, totaling4,597 nuclei across eight RNA libraries. Same-age shared measured genes determine30 soft cosine neighbors (temperature0.1), with separate-assay scaling,SDfloor0.25,clipping5 and L2 normalization. Per-gene reference log1pCP10k expression is averaged with these weights; program scores average the relevant genes. These are expected normalized RNA-reference values, **not Xenium transcript counts or probabilities of active sprouting**. No spatial coordinate, inferred target label or desired program score enters the neighbor fit.

Measured Angpt2 is allowed in the primary alignment because it is actually available on Xenium. The five genuinely unmeasured targets and six cycling-control targets are absent from fitting genes. Angpt2 withholding is a sensitivity that also removes it from the predictor normalization denominator. RNA reference gene-support tables show where target genes were sparse or absent before any transfer. In particular, an unexpressed reference gene cannot become independent spatial evidence through projection.

This differs from the earlier conservative Ridge/rank-excess benchmark: it transfers reference gene-expression expectations, allows observed Angpt2 and evaluates disjoint RNA metacells as well as individual nuclei. The earlier benchmark remains intact. The present six-gene total can benefit from measured Angpt2, which is why the unmeasured-five score and its validation are central.

## Technical support for the unmeasured component

Both RNA-library directions at the same age must achieve Spearman>=0.2, within-family residual Spearman>=0.15 and positive out-of-library R². The primary setting is assay-scaled30-neighbor transfer; Pearson,10/60-neighbor and Angpt2-withheld settings are reported transparently rather than selecting a winner from the same held-out data.

'''+ '\n'.join(gates)+'''

**Fine recovery of the five unmeasured genes remains weak.** Coarse20-nucleus P21 RNA groups meet this technical gate, which is limited support at a coarser level, not independent spatial validation. The total six-gene score is more recoverable at some ages, but part of that is measured Angpt2. All released spatial scores therefore remain explicitly **exploratory reference-derived hypotheses**. A high score alone does not establish tip-cell identity, completed sprouting or increased vessel formation.

![Projection benchmark](RNA_projection_benchmark.png)

## Files and use

For each `P*/` directory:

- **`Xenium_ECs_with_inferred_sprouting.h5ad`** is ready to load. `X` contains original observed native Xenium counts for403 targets. `layers['nuclear_only']` and `layers['boundary_core']` retain observed ownership views. `obs` contains the inferred scores, QC status, RNA-library difference, boundary difference and reference-similarity flags. `obsm['spatial']` holds original micrometer coordinates. Predicted gene expression is a separate `obsm` table, never a count layer.
- `EC_annotations_with_inferred_sprouting.tsv.gz` preserves every source EC candidate, including unprojected QC-review rows. `all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz` covers every source cardiac profile. The root all-stage table combines all four sections with unique `stage_cell_id` keys.
- `neighborhood5_inferred_sprouting.tsv.gz` and `neighborhood20_inferred_sprouting.tsv.gz` give source cell IDs, physical radius and native/nuclear/core projected scores. Neighborhood values use a different RNA reference resolution and should not be compared numerically to single-nucleus scores as the same estimator.
- Nuclear/native/core detailed tables contain raw-reference and neighborhood-size sensitivities. Neighbor NPZ files store exact RNA source IDs and weights for reproducibility. Reference-neighbor score SD describes reference heterogeneity, not a calibrated uncertainty interval.
- Opposite-library benchmarks retain observed versus inferred RNA scores. RNA metacells are disjoint near-neighbor groups formed only from measured fitting genes, with small final remainders merged into a nearby group. Group sizes/source nuclei are exported. Sparse reference groups cap the effective number of neighbors; requested and effective counts are both recorded.

```python
import scanpy as sc
ec = sc.read_h5ad("P7/Xenium_ECs_with_inferred_sprouting.h5ad")
sc.pl.embedding(ec, basis="spatial", color=[
    "inferred_sprouting_score_6g",
    "inferred_unmeasured_sprouting_score_5g",
])
inferred_genes = ec.obsm["inferred_sprouting_genes_RNA_log1pCP10k"]
```

![P7 inferred layer alongside observed Angpt2](P7_inferred_sprouting_and_measured_Angpt2.png)

## Limits and provenance

This projection preserves source counts, coordinates, gene targets and EC labels. A low nearest-reference similarity, source-library difference or boundary sensitivity flags possible instability; none independently calibrates biological accuracy. Pooling can span adjacent vessels and increases spatial smoothing. Two RNA libraries per age are not assumed independent animals; one Xenium section per age does not establish reproducible developmental spatial rates. Inferred Apln/Esm1 require direct experimental validation before being treated as measured locations.

Reproduce with `python analysis/project_RNA_sprouting_to_Xenium.py --stage P0` (andP7/P14/P21), then `analysis/report_projected_sprouting.py` and `analysis/validate_projected_sprouting.py`. Original input matrices remain in `results/endothelial_QC_subtypes/`. `validation.json`, source count hashes and the file manifest record exact count preservation, weighted-reference reconstruction, withheld targets, RNA-metacell membership, spatial pool radii and complete atlas coverage.
''')
    print('PROJECTED_SPROUTING_REPORT_COMPLETE',flush=True)

if __name__=='__main__':main()
