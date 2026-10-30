"""Module-first score transfer, robust-score alternatives and explicit equivalence.

The existing mean-of-six projection is already a module. This analysis tests
whole-module definitions without selecting a winner on held-out query results.
"""
from project_RNA_sprouting_to_Xenium import *
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

MODULE_OUT=ROOT/'results/EC_gene_modules'
MODULES={'Sprouting_6g':SPROUT,'Unmeasured_sprouting_5g':UNMEASURED,'Proliferation_6g':CYCLING}
DEFINITIONS=['mean_log_expression','gene_balanced','rank_enrichment']

def ranks(counts,genes):
    positions=pd.Index(genes).get_indexer(OUTCOMES);assert (positions>=0).all()
    values=counts[:,positions].toarray();result=np.zeros_like(values,dtype=np.float32)
    for i in range(len(values)):
        a=counts.data[counts.indptr[i]:counts.indptr[i+1]];a=np.sort(a[a>=.001])
        keep=values[i]>=.001;v=values[i,keep]
        rank=len(a)-.5*(np.searchsorted(a,v,'left')+np.searchsorted(a,v,'right'))+.5
        result[i,keep]=np.maximum(0,1-rank/1500)
    return result

def make_scores(y,rank,mu,sd):
    # Calibration is frozen from the RNA training library in every benchmark.
    values={'mean_log_expression':y,'gene_balanced':np.clip((y-mu)/sd,-3,3),'rank_enrichment':rank}
    data={}
    for name,panel in MODULES.items():
        ix=[OUTCOMES.index(g) for g in panel]
        for definition,v in values.items():data[name+'__'+definition]=v[:,ix].mean(axis=1)
    return pd.DataFrame(data)

def library_units(stage,rna,genes,lib,size):
    m=rna.obs['sample'].eq(lib).to_numpy();obs=rna.obs.loc[m].copy();counts=rna.X[m].tocsr()
    if size>1:
        meta=pd.read_csv(PROJ/stage/f'RNA_{lib}_corrected_metacells{size}.tsv.gz',sep='\t',index_col=0)
        local=obs.index
        row=[];col=[]
        for j,s in enumerate(meta.source_cell_ids):
            ix=local.get_indexer(s.split(';'));assert (ix>=0).all();row.extend([j]*len(ix));col.extend(ix)
        design=sparse.csr_matrix((np.ones(len(row)),(row,col)),shape=(len(meta),len(obs)))
        counts=(design@counts).astype(np.float32).tocsr();obs=meta
    target=target_expression(counts,pd.Index(genes).get_indexer(OUTCOMES))
    return {'counts':counts,'obs':obs,'target':target,'rank':ranks(counts,genes)}

def main():
    MODULE_OUT.mkdir(exist_ok=True);rows=[];summaries=[];equivalence=[];models=[];sources=[]
    fit=fitting_genes()
    for stage in STAGES:
        print('MODULE_STAGE',stage,flush=True)
        folder=MODULE_OUT/stage;folder.mkdir(exist_ok=True)
        rna=ad.read_h5ad(OUT/'snRNA'/stage/'EC_corrected_and_raw_counts.h5ad')
        genes=rna.var_names.to_numpy(dtype=str);fitpos=pd.Index(genes).get_indexer(fit)
        libs=sorted(rna.obs['sample'].unique());refs={}
        for size in [1,5,20]:
            for lib in libs:refs[(size,lib)]=library_units(stage,rna,genes,lib,size)
        for size in [1,5,20]:
            for qlib in libs:
                rlib=next(l for l in libs if l!=qlib);r,q=refs[(size,rlib)],refs[(size,qlib)]
                mu=r['target'].mean(axis=0);sd=np.maximum(r['target'].std(axis=0),.25)
                yr=make_scores(r['target'],r['rank'],mu,sd);yq=make_scores(q['target'],q['rank'],mu,sd)
                family=(q['obs'].EC_family if size==1 else q['obs'].dominant_family).to_numpy(dtype=str)
                for mode,use in [('Angpt2_observed_allowed',np.ones(len(fit),bool)),('all_sprouting_targets_withheld',fit!='Angpt2')]:
                    rf=norm_features(r['counts'][:,fitpos[use]].toarray());qf=norm_features(q['counts'][:,fitpos[use]].toarray())
                    pred,*_=transfer(rf,qf,yr.to_numpy(),30,True)
                    for j,col in enumerate(yr):
                        module,definition=col.split('__');actual=yq[col].to_numpy();p=pred[:,j]
                        rho=corr(actual,p);within=corr(family_residual(actual,family),family_residual(p,family))
                        r2=float(r2_score(actual,p)) if np.std(actual)>1e-9 else np.nan
                        rows.append({'stage':stage,'query_library':qlib,'RNA_resolution_nuclei':size,'alignment':mode,
                            'module':module,'score_definition':definition,'query_units':len(actual),'Spearman':rho,'within_family_Spearman':within,
                            'R2_vs_query_mean':r2,'technical_gate':bool(np.isfinite(rho) and rho>=.2 and np.isfinite(within) and within>=.15 and np.isfinite(r2) and r2>0),
                            'calibration':'Training-library mean/SD only for balanced module; SDfloor0.25,clip3; rankmax1500,numerical floor0.001'})
        # New module fields use the already fixed RNA-reference weights. No
        # alignment or spatial data is tuned to make a module appear stronger.
        source=ad.read_h5ad(PROJ/stage/'Xenium_ECs_with_inferred_sprouting.h5ad')
        y=target_expression(rna.X,rna.var_names.get_indexer(OUTCOMES));rank=ranks(rna.X.tocsr(),genes)
        # Equal RNA-library moments for within-stage calibration, not across-age activity.
        mu=np.mean([y[rna.obs['sample'].eq(l).to_numpy()].mean(axis=0) for l in libs],axis=0)
        second=np.mean([(y[rna.obs['sample'].eq(l).to_numpy()]**2).mean(axis=0) for l in libs],axis=0)
        sd=np.maximum(np.sqrt(np.maximum(second-mu**2,0)),.25)
        score=make_scores(y,rank,mu,sd);score.index=rna.obs_names
        score.to_csv(folder/'RNA_reference_gene_module_scores.tsv.gz',sep='\t')
        with np.load(PROJ/stage/'nuclear_RNA_reference_neighbors.npz') as b:
            np.testing.assert_array_equal(b['query_cell_ids'].astype(str),source.obs_names)
            ri=score.index.get_indexer(b['reference_cell_ids'].astype(str));assert (ri>=0).all()
            projected=np.einsum('nk,nkg->ng',b['neighbor_weights'],score.iloc[ri].to_numpy()[b['neighbor_indices']])
        annotation=source.obs.copy()
        for j,col in enumerate(score):annotation['inferred_module__'+col]=projected[:,j]
        annotation['inferred_sprouting_module']=annotation['inferred_module__Sprouting_6g__mean_log_expression']
        annotation['inferred_sprouting_module_balanced']=annotation['inferred_module__Sprouting_6g__gene_balanced']
        annotation['inferred_sprouting_module_rank']=annotation['inferred_module__Sprouting_6g__rank_enrichment']
        annotation['inferred_unmeasured_sprouting_module']=annotation['inferred_module__Unmeasured_sprouting_5g__mean_log_expression']
        annotation['inferred_proliferation_module']=annotation['inferred_module__Proliferation_6g__mean_log_expression']
        a=annotation['inferred_module__Sprouting_6g__mean_log_expression'].to_numpy()
        b=source.obs.inferred_sprouting_score_6g.to_numpy();error=float(np.max(np.abs(a-b)))
        assert np.allclose(a,b,atol=1e-6,rtol=1e-5)
        equivalence.append({'stage':stage,'projected_ECs':len(source),'module_first_vs_gene_first_max_abs_error':error,
            'meaning':'Weighted mean of a mean equals mean of weighted gene values; aggregation order alone does not improve recovery.'})
        source.obs=annotation
        source.uns['gene_module_scores']={'modules':MODULES,'definitions':{
            'mean_log_expression':'Mean normalized reference gene expression; original six-gene score, directly projected as a module.',
            'gene_balanced':'Mean clipped gene Z values using equal-library same-age RNA moments; SDfloor0.25,clip3. Within-stage relative score.',
            'rank_enrichment':'Mean within-nucleus top1500 target-gene rank enrichment; dustfloor0.001. Not probability of activation.'},
            'status':'Exploratory inferred RNA modules; no new measured genes or biological activation validation',
            'default':'Existing mean-log-expression module preserved; alternatives are sensitivity layers, not a selected winner.'}
        source.write_h5ad(folder/'Xenium_ECs_with_gene_modules.h5ad',compression='gzip')
        keep=['stage_cell_id','x_centroid','y_centroid','nuclear_supported_family','sprouting_annotation_status']+[c for c in annotation if c.startswith('inferred_module__')]+[
            'inferred_sprouting_module','inferred_sprouting_module_balanced','inferred_sprouting_module_rank','inferred_unmeasured_sprouting_module','inferred_proliferation_module']
        annotation[keep].to_csv(folder/'inferred_gene_module_annotations.tsv.gz',sep='\t',na_rep='NA')
        models.append(annotation[keep].assign(stage=stage))
    bench=pd.DataFrame(rows);save(bench,MODULE_OUT/'module_definition_opposite_library_benchmarks.tsv')
    group=['stage','RNA_resolution_nuclei','alignment','module','score_definition']
    summary=bench.groupby(group,observed=True).agg(min_Spearman=('Spearman','min'),max_Spearman=('Spearman','max'),
        min_within_family_Spearman=('within_family_Spearman','min'),min_R2=('R2_vs_query_mean','min'),
        both_libraries_gate=('technical_gate','all')).reset_index()
    save(summary,MODULE_OUT/'module_definition_benchmark_summary.tsv')
    save(pd.DataFrame(equivalence),MODULE_OUT/'module_aggregation_order_equivalence.tsv')
    pd.concat(models).to_csv(MODULE_OUT/'all_stage_gene_module_annotations.tsv.gz',sep='\t',na_rep='NA')
    coherence=pd.read_csv(DEST/'RNA/within_library_gene_coherence.tsv',sep='\t')
    save(coherence[coherence.program.eq('Original_sprouting_core')&coherence.layer.eq('DecontX')],MODULE_OUT/'sprouting_module_gene_coherence.tsv')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    with PdfPages(MODULE_OUT/'EC_gene_module_comparison.pdf') as pdf:
        fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
        for ax,module in zip(axes,MODULES):
            d=bench[bench.module.eq(module)&bench.alignment.eq('Angpt2_observed_allowed')&bench.RNA_resolution_nuclei.eq(1)]
            for i,stage in enumerate(STAGES):
                for j,(definition,c) in enumerate(zip(DEFINITIONS,['#0072b2','#e69f00','#009e73'])):
                    v=d.loc[d.stage.eq(stage)&d.score_definition.eq(definition),'Spearman'];ax.scatter(np.full(len(v),i+(j-1)*.18),v,s=28,color=c,label=definition if i==0 else None)
            ax.axhline(.2,c='gray',ls='--');ax.set_xticks(range(4),STAGES);ax.set_ylabel('Opposite-library Spearman');ax.set_title(module.replace('_',' '));ax.legend(fontsize=7,frameon=False)
        fig.suptitle('Direct RNA gene-module transfer · individual nuclei\nPredeclared module definitions, both library directions; no winner selected on held-out results',fontsize=14)
        fig.savefig(MODULE_OUT/'module_definition_comparison.png',dpi=180,bbox_inches='tight');pdf.savefig(fig);plt.close(fig)
        for stage,annotation in zip(STAGES,models):
            fig,axes=plt.subplots(1,3,figsize=(17,6),layout='constrained')
            for ax,definition in zip(axes,DEFINITIONS):
                col='inferred_module__Sprouting_6g__'+definition;lo,hi=annotation[col].quantile([.01,.99])
                p=ax.scatter(annotation.x_centroid,annotation.y_centroid,c=annotation[col],s=2.5,cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
                ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(definition.replace('_',' '));ax.set(xlabel='X (µm)',ylabel='Y (µm)')
                fig.colorbar(p,ax=ax,shrink=.6,label='Inferred RNA module units; definition-specific')
            fig.suptitle(f'{stage} · same EC positions, three whole-module score definitions\nExploratory sprouting module; color stability alone is not biological validation',fontsize=14)
            fig.savefig(MODULE_OUT/(stage+'_sprouting_module_definitions.png'),dpi=180,bbox_inches='tight');pdf.savefig(fig);plt.close(fig)
    # Export actual benchmarks rather than claiming a robust score must improve them.
    gates=summary[summary.module.eq('Unmeasured_sprouting_5g')&summary.alignment.eq('Angpt2_observed_allowed')]
    counts=gates.groupby('score_definition').both_libraries_gate.sum().to_dict()
    write_json(MODULE_OUT/'results_summary.json',{'default_module_already_present':True,'projection_aggregation_order_equivalent':True,
        'module_definitions_tested':DEFINITIONS,'unmeasured_module_passed_stage_resolution_pairs':counts,
        'source_counts_and_default_sprouting_scores_unchanged':True,'all_inferred_spatial_alternatives_exploratory':True})
    key=[]
    for module,panel in MODULES.items():
        for gene in panel:key.append({'module':module,'gene':gene,'Xenium_measured':gene=='Angpt2' if module!='Proliferation_6g' else gene in ['Mki67','Top2a','Cenpf']})
    save(pd.DataFrame(key),MODULE_OUT/'gene_module_definitions.tsv')
    (MODULE_OUT/'README.md').write_text('''# EC gene modules: whole-program scores rather than individual-gene claims

**The existing sprouting score already is a six-gene module score.** This analysis explicitly scores the RNA module first and transfers that single score onto Xenium. With fixed reference weights, `projection(mean of genes) = mean of projected genes`. The differences are only floating-point rounding (less than0.000001), so changing aggregation order alone adds no predictive information.

[Module comparison report](EC_gene_module_comparison.pdf) · [All-stage module annotations](all_stage_gene_module_annotations.tsv.gz) · [Module definitions](gene_module_definitions.tsv) · [Opposite-library benchmark summary](module_definition_benchmark_summary.tsv) · [Aggregation-equivalence check](module_aggregation_order_equivalence.tsv) · [Validation](validation.json)

## Primary module-level readout

Each `P*/Xenium_ECs_with_gene_modules.h5ad` preserves the previous observed403-target count matrix, nuclear/boundary count layers, source identities and spatial coordinates. The main fields in `obs` are:

- **`inferred_sprouting_module`**: mean reference log1pCP10k of Apln,Esm1,Angpt2,Cxcr4,Pgf,Adm. This is numerically the existing six-gene score, now explicitly presented as one module.
- **`inferred_unmeasured_sprouting_module`**: diagnostic five-gene component excluding Angpt2, so measured Angpt2 cannot be mistaken for recovery of missing genes.
- **`inferred_proliferation_module`**: Mki67,Top2a,Cenpf,Ccnb1,Ube2c,Cdk1, kept separate from sprouting. Proliferation does not define tip-cell identity.
- **`inferred_sprouting_module_balanced`** and **`inferred_sprouting_module_rank`**: sensitivity versions of the same sprouting gene set, not automatically adopted replacements or newly validated activity measures.

The missing-five component is a validation diagnostic, not another biological EC type or independent program. All scores are inferred RNA annotations; the five requested EC identities and all actual Xenium measurements are preserved. Previous per-gene predictions remain available for auditing but are not the main module-level interpretation.

## Alternative whole-module definitions tested

1. **Mean normalized expression:** the existing mean of gene log1pCP10k values.
2. **Gene-balanced module:** standardize genes using training-library mean/SD, an SDfloor0.25 and clipping±3, then average. This reduces the influence of a single highly expressed or unusually variable gene. Spatial calibration uses equal-library moments within the matched-age reference; these scores are relative within age, not cross-age absolute activity.
3. **Rank enrichment:** mean target-gene enrichment among the top1500 expressed genes of each RNA nucleus/metacell, retaining fractional counts and excluding numerical dust below0.001. This is a rank statistic, not a probability of pathway activation.

Every module is calculated in RNA before projection. The spatial comparison keeps exactly the previously fixed reference IDs/weights; changing a score definition cannot tune the alignment toward a desired map. The three resulting maps are sensitivity comparisons, not independent experiments or a selected best-performing model.

## Does a module help here?

Aggregation can reduce independent gene-level noise and is the right presentation level for a multi-gene hypothesis. It cannot create information absent from the measured predictors, and a gene list is not automatically a coordinated biological program. The six-gene signature has weak/inconsistent within-library coherence at early stages (`sprouting_module_gene_coherence.tsv`), with sparse reference genes. Individual-gene plots remain useful for identifying single-gene domination or contamination.

We tested all three definitions in both RNA-library directions, separately atP0/P7/P14/P21, using individual nuclei and the existing disjoint5-/20-nucleus metacells. Calibration for gene-balanced scores uses only the training library. Shared-gene alignment is tested with measuredAngpt2 allowed and with all sprouting targets withheld, including from predictor normalization. Cycling targets are absent from predictors in all tests. The usual technical gates require Spearman>=0.2, family-residual Spearman>=0.15 and positive R² in both library directions.

The withholding sensitivity reuses the previously fixed RNA metacells, whose construction allowed measuredAngpt2. It therefore tests withholding from alignment predictors, not a fully independent ablation of metacell construction. The five unmeasured genes were absent from both alignment and original metacell-construction features throughout.

**Alternative scoring did not consistently improve the unmeasured-five component.** Mean-expression and balanced modules each pass only the coarseP21/20-nucleus condition; the rank module passes none of those stage/resolution pairs. This is a limited technical result, not proof that no possible module could work. The current mean module remains the default, and no alternative is promoted by choosing a winner on the same held-out query data. Some proliferation-module benchmarks benefit from a different definition, underscoring that score behavior is program-dependent.

![Module-definition comparison](module_definition_comparison.png)

## Interpretation and use

Use the sprouting module map as an **exploratory reference-derived propensity**, with the missing-five component and source-library/boundary uncertainty beside it. A high inferred module score does not establish actual sprouting, a tip-cell identity or new-vessel growth. Xenium directly measures onlyAngpt2 from this particular six-gene signature. Broader remodeling and cycling signatures remain separate from the sprouting signature; adding generic EC markers would change the biological question rather than validate sprouting.

```python
import scanpy as sc
ec = sc.read_h5ad("P7/Xenium_ECs_with_gene_modules.h5ad")
sc.pl.embedding(ec, basis="spatial", color=[
    "inferred_sprouting_module",
    "inferred_unmeasured_sprouting_module",
    "inferred_proliferation_module",
])
```

These files cover the same29,363 eligible EC positions as the [previous projection](../projected_sprouting_Xenium/README.md). `stage_cell_id` uniquely identifies positions across sections. Coarse RNA groups and spatial neighborhoods are not independent biological replicates; one Xenium section per age and two RNA libraries per age limit generalization. A genuinely coherent, independently validated larger sprouting module could help, but that has not been established by this score-definition comparison.

Reproduce with `python analysis/evaluate_EC_gene_modules.py`, then `python analysis/validate_EC_gene_modules.py` in the existing environment. No source matrix, reference weight, EC identity or coordinate is changed. The validation checks observed counts, module-to-reference arithmetic, exact correspondence to the original mean score, atlas membership, benchmark coverage, image/PDF integrity and source hashes.
''')
    print(pd.DataFrame(equivalence).to_string(index=False),flush=True)
    print('UNMEASURED_MODULE_GATES',counts,flush=True)

if __name__=='__main__':main()
