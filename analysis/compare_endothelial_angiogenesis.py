"""Descriptive stage comparison of endothelial cycling and vascular programs."""
from endothelial_common import *
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEST=ROOT/'results/endothelial_angiogenesis'
PROGRAMS={
    'Cycling':CYCLING,
    'DNA_replication':'Mcm2 Mcm3 Mcm4 Mcm5 Mcm6 Pcna Tyms Rrm2'.split(),
    'Sprouting_associated':'Apln Esm1 Angpt2 Cxcr4 Pgf Adm'.split(),
    'VEGF_Notch_machinery':'Kdr Flt1 Nrp1 Nrp2 Dll4 Hey1 Ets1 Elk3'.split(),
    'Stabilization_associated':'Klf2 Klf4 Tek Cdh5 S1pr1 Cldn5 Ocln'.split(),
}
CAP_MARKERS='Dach1 Cd36 Rgcc Car4 Emcn Gpihbp1'.split()

def main():
    DEST.mkdir(exist_ok=True)
    library_rows=[];gene_rows=[];coverage=[];cell_rows=[]
    for stage in STAGES:
        a=ad.read_h5ad(OUT/'snRNA'/stage/'EC_corrected_and_raw_counts.h5ad')
        b=ad.read_h5ad(OUT/'snRNA'/stage/'EC_review_corrected_and_raw_counts.h5ad')
        data=ad.concat([a,b],merge='same')
        primary=data.obs.primary_EC_eligible.to_numpy()
        cap=(np.asarray((data.X[:,data.var_names.isin(CAP_MARKERS)]>=1).sum(axis=1)).ravel()>=2)
        cohorts={'primary_QC':primary,'all_previously_labeled_EC':np.ones(len(data),bool),
                 'primary_ambient_le30pct':primary&data.obs.decontx_contamination.le(.3).to_numpy(),
                 'primary_capillary_marker_supported':primary&cap}
        for layer,counts in [('DecontX_corrected',data.X),('raw',data.layers['counts_raw'])]:
            z=lognorm(counts)
            cycle=np.asarray((counts[:,data.var_names.isin(CYCLING)]>=1).sum(axis=1)).ravel()>=2
            if layer=='DecontX_corrected':
                np.testing.assert_array_equal(cycle,data.obs.cycling_marker_count.ge(2))
            scores={}
            for name,panel in PROGRAMS.items():
                genes=[g for g in panel if g in data.var_names]
                scores[name]=np.asarray(z[:,data.var_names.isin(genes)].mean(axis=1)).ravel()
                if stage=='P0' and layer=='DecontX_corrected':
                    coverage.extend({'program':name,'gene':g,'RNA_measured':g in data.var_names} for g in panel)
            for cohort,m in cohorts.items():
                for lib in sorted(data.obs['sample'].unique()):
                    mask=m&data.obs['sample'].eq(lib).to_numpy()
                    if not mask.any():continue
                    row={'stage':stage,'sample':lib,'cohort':cohort,'layer':layer,'nuclei':int(mask.sum()),
                        'cycling_nuclei':int(cycle[mask].sum()),'cycling_fraction':float(cycle[mask].mean()),
                        'median_UMI':float(np.median(np.asarray(counts[mask].sum(axis=1)).ravel())),
                        'median_estimated_ambient':float(data.obs.loc[mask,'decontx_contamination'].median())}
                    for name,score in scores.items():row[name+'_mean_log1pCP10k']=float(score[mask].mean())
                    library_rows.append(row)
                    for gene in dict.fromkeys(sum(PROGRAMS.values(),[])):
                        if gene not in data.var_names:continue
                        j=data.var_names.get_loc(gene)
                        gene_rows.append({'stage':stage,'sample':lib,'cohort':cohort,'layer':layer,'gene':gene,'nuclei':int(mask.sum()),
                            'mean_log1pCP10k':float(z[mask,j].mean()),'fraction_count_ge1':float((counts[mask,j]>=1).mean())})
            if layer=='DecontX_corrected':
                d=data.obs[['age','sample','primary_EC_eligible','cycling_marker_count','decontx_contamination']].copy()
                d['capillary_marker_supported']=cap
                for name,score in scores.items():d[name+'_mean_log1pCP10k']=score
                cell_rows.append(d)
    libraries=pd.DataFrame(library_rows);save(libraries,DEST/'programs_and_cycling_by_library.tsv')
    save(pd.DataFrame(gene_rows),DEST/'individual_genes_by_library.tsv')
    save(pd.DataFrame(coverage),DEST/'RNA_program_definitions.tsv')
    pd.concat(cell_rows).to_csv(DEST/'EC_program_annotations.tsv.gz',sep='\t',na_rep='NA')
    primary=libraries[libraries.cohort.eq('primary_QC')&libraries.layer.eq('DecontX_corrected')]
    summary=primary.groupby('stage',observed=True).agg(EC_nuclei=('nuclei','sum'),cycling_nuclei=('cycling_nuclei','sum'),
        library_min_fraction=('cycling_fraction','min'),library_max_fraction=('cycling_fraction','max')).reindex(STAGES)
    summary['cycling_fraction']=summary.cycling_nuclei/summary.EC_nuclei
    save(summary.reset_index(),DEST/'primary_EC_cycling_by_stage.tsv')
    # Compare equal-weight library averages. No cells-as-independent-animal p-values.
    early_late=[]
    scorecols=['cycling_fraction']+[k+'_mean_log1pCP10k' for k in PROGRAMS]
    for (cohort,layer),d in libraries.groupby(['cohort','layer'],observed=True):
        for score in scorecols:
            early=d.loc[d.stage.isin(['P0','P7']),score].mean()
            late=d.loc[d.stage.isin(['P14','P21']),score].mean()
            early_late.append({'cohort':cohort,'layer':layer,'measure':score,'early_equal_library_mean':early,
                'late_equal_library_mean':late,'early_minus_late':early-late,'early_late_ratio':early/late if late>0 else np.nan})
    save(pd.DataFrame(early_late),DEST/'early_vs_late_descriptive_comparison.tsv')
    fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    for ax,col,title in [(axes[0],'cycling_fraction','Cycling-marker-positive EC nuclei'),
                        (axes[1],'DNA_replication_mean_log1pCP10k','DNA-replication-associated RNA'),
                        (axes[2],'Sprouting_associated_mean_log1pCP10k','Sprouting-associated RNA panel')]:
        values=[]
        for i,s in enumerate(STAGES):
            v=primary.loc[primary.stage.eq(s),col].to_numpy();v=100*v if col=='cycling_fraction' else v
            ax.scatter([i-.06,i+.06],v,s=45,color='#0072b2',zorder=3)
            ax.plot([i,i],[min(v),max(v)],color='#0072b2',lw=1)
            values.append(v.mean())
        ax.plot(range(4),values,color='#0072b2',lw=1.5);ax.set_xticks(range(4),STAGES);ax.set_title(title,fontsize=12)
        ax.set_ylabel('Percent of EC nuclei' if col=='cycling_fraction' else 'Mean log1p CP10k')
        ax.spines[['right','top']].set_visible(False)
    fig.suptitle('Endothelial programs during postnatal development\nDots: two RNA libraries per stage; lines: equal-library means and observed range, not confidence intervals',fontsize=14)
    fig.savefig(DEST/'EC_cycling_and_sprouting_by_stage.png',dpi=180,bbox_inches='tight')
    fig.savefig(DEST/'EC_cycling_and_sprouting_by_stage.pdf',bbox_inches='tight');plt.close(fig)
    table=['| Stage | Cycling-positive / EC nuclei | Percentage |','|---|---:|---:|']
    for s,d in summary.iterrows():table.append(f'| {s} | {int(d.cycling_nuclei)} / {int(d.EC_nuclei):,} | {100*d.cycling_fraction:.2f}% |')
    (DEST/'README.md').write_text('''# Endothelial proliferation and angiogenesis-associated RNA across ages

The primary QC endothelial cohort shows more cycling-marker-positive nuclei at P0/P7, a decline at P14 and very few at P21. This supports greater **early endothelial proliferative activity**, an important component of angiogenesis. The six cycling genes are Mki67, Top2a, Cenpf, Ccnb1, Ube2c and Cdk1; at least two need a corrected count of at least one.

'''+ '\n'.join(table)+'''

![Cycling, DNA replication and sprouting panels](EC_cycling_and_sprouting_by_stage.png)

Proliferation is only one component of new-vessel growth. The non-cycling sprouting-associated panel (Apln, Esm1, Angpt2, Cxcr4, Pgf, Adm), VEGF/Notch machinery and stabilization-associated programs are reported separately, without assuming all genes peak together or calling RNA receptor levels direct signaling activity. Scores are mean log1p CP10k over the stated genes, not a calibrated angiogenesis index. Individual genes and raw/corrected comparisons accompany the scores. DNA-replication genes provide complementary evidence but can also participate in DNA repair.

The strongest stage difference is proliferation: the equal-library mean cycling fraction is approximately five times higher in the combined P0/P7 group than P14/P21. DNA-replication-associated RNA also declines. The direction persists with raw counts, lower estimated ambient contamination, all previously labeled ECs and the capillary-marker-supported subset. P0 and P7 cycling fractions are similar; these data do not establish a distinct P7 proliferation peak.

Broader sprouting genes are heterogeneous. Apln expression is highest at P7, Angpt2 is higher at P0 and remains detectable at P14, while Esm1 is sparse. The multi-gene sprouting panel changes much less than the cycling program. VEGF/Notch machinery expression is higher later in this comparison; receptor/Notch RNA is not a direct measure of signaling activation or vessel-growth rate. The appropriate conclusion is **greater early EC proliferative growth, with continued and changing vascular programs later**, rather than every angiogenic pathway declining in parallel.

The paper reports an overall decline in Ki67-positive ECs from P0 to P21 and early coronary plexus expansion (PMC13561863, Supplementary Fig10b). Our RNA cycling fractions are not the paper's immunostaining measurements and do not prove completed cell division or directly measure new perfused vessels.

`programs_and_cycling_by_library.tsv` preserves both RNA library labels at every age. Sensitivities include all 5,573 previously labeled ECs, primary nuclei with estimated ambient<=30%, original raw counts and a primary capillary-marker-supported subset (>=2 of Dach1/Cd36/Rgcc/Car4/Emcn/Gpihbp1). This last subset is an expression-based enrichment check, not a newly validated capillary annotation. It includes subtype-unresolved cells: restricting only to the previous Capillary-like label would omit the P7 cycling group, producing an artificial zero. Changing subtype proportions, RNA depth, capture and QC selection still limit causal or quantitative organ-wide conclusions.

`early_vs_late_descriptive_comparison.tsv` compares equal-weight means of four early versus four late RNA libraries. Two libraries per age are not assumed independent animals; no cell-level or animal-level significance claim is made. The per-library gene table is the evidence for heterogeneous gene-program trends. These are descriptive cross-sectional comparisons, not direct measurements of angiogenesis flux.

Source RNA QC, identities, count matrices and coordinates remain unchanged. Reproduce with `python analysis/compare_endothelial_angiogenesis.py` in the existing environment. See `RNA_program_definitions.tsv`, `validation.json` and `file_manifest.tsv` for definitions, accounting and provenance.
''')
    assert len(primary)==8 and int(summary.EC_nuclei.sum())==4597 and int(summary.cycling_nuclei.sum())==77
    write_json(DEST/'validation.json',{'status':'PASS','primary_EC_nuclei':4597,'cycling_positive_primary':77,
        'RNA_libraries':8,'cycling_counts_reproduced_from_actual_corrected_RNA':True,'source_counts_unchanged':True,
        'programs':PROGRAMS,'no_independent_animal_or_direct_angiogenesis_claim':True})
    files=[p for p in DEST.iterdir() if p.is_file() and p.name!='file_manifest.tsv']+[Path(__file__)]
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]),DEST/'file_manifest.tsv')
    print(summary.to_string(),flush=True)
    print(pd.DataFrame(early_late).to_string(index=False),flush=True)

if __name__=='__main__':main()
