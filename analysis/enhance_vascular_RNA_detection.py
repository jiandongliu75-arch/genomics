"""Fractional-aware EC RNA scoring, pseudobulks, ambient audit and holdouts."""
from vascular_detection_common import *
import anndata as ad
import h5py
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
import gc

def load_primary():
    parts=[ad.read_h5ad(OUT/'snRNA'/s/'EC_corrected_and_raw_counts.h5ad') for s in STAGES]
    return ad.concat(parts,merge='same')

def background_audit(genes):
    """Non-EC cellular RNA is a sensitivity proxy, not an empty-droplet soup."""
    rows=[]
    with h5py.File(RNA_INPUT) as f:
        obs,var=read_elem(f['obs']),read_elem(f['var'])
        wanted=np.flatnonzero(var.gene_symbol.isin(TARGET_GENES))
        symbols=var.gene_symbol.iloc[wanted].to_numpy(dtype=str)
        matrix=sparse_dataset(f['layers/counts_raw'])
        totals={};numbers={}
        for start in range(0,len(obs),2048):
            o=obs.iloc[start:start+2048]
            x=matrix[start:start+2048][:,wanted].toarray()
            for (lib,lineage),indices in o.groupby(['sample','reviewed_identity'],observed=True).indices.items():
                key=(str(lib),str(lineage))
                v=x[indices].sum(axis=0)
                totals[key]=totals.get(key,np.zeros(len(wanted)))+v
                n,umi=numbers.get(key,(0,0))
                numbers[key]=(n+len(indices),umi+float(o.iloc[indices].raw_total_counts.sum()))
        for (lib,lineage),v in totals.items():
            n,umi=numbers[(lib,lineage)]
            # Explicit sum if an input symbol has multiple Ensembl features.
            values=pd.Series(v,index=symbols).groupby(level=0).sum()
            for g,value in values.items():
                rows.append({'sample':lib,'lineage':lineage,'gene':g,'nuclei':n,'raw_gene_molecules':value,'raw_total_UMI':umi,
                    'raw_CPM':value*1e6/max(umi,1)})
    table=pd.DataFrame(rows);save(table,DEST/'RNA/nonEC_background_lineage_audit.tsv')
    return table

def benchmark(data,score):
    """Can measured shared genes predict withheld program RNA across libraries?"""
    with np.load(ROOT/'results/xenium_ventricular_mapping/reference_shared_panel_counts.npz') as b:shared=b['genes'].astype(str)
    fitting=[g for g in shared if g not in set(TARGET_GENES+LINEAGES['Cardiomyocyte']+SEX)]
    save(pd.DataFrame({'gene':fitting}),DEST/'RNA/benchmark_fitting_genes.tsv')
    positions=data.var_names.get_indexer(fitting);assert (positions>=0).all()
    counts=data.X[:,positions].toarray().astype(np.float32)
    features=np.log1p(1e4*counts/np.maximum(counts.sum(axis=1),1)[:,None])
    rows=[]
    for stage in STAGES:
        stage_mask=data.obs.age.eq(stage).to_numpy()
        for lib in sorted(data.obs.loc[stage_mask,'sample'].unique()):
            train=stage_mask&~data.obs['sample'].eq(lib).to_numpy();test=stage_mask&data.obs['sample'].eq(lib).to_numpy()
            mu=features[train].mean(axis=0);sd=np.maximum(features[train].std(axis=0),.25)
            r=np.clip((features[train]-mu)/sd,-5,5)
            q=np.clip((features[test]-mu)/sd,-5,5)
            model=Ridge(alpha=10).fit(r,score.loc[train,list(PROGRAMS)].to_numpy())
            pred=model.predict(q)
            family=data.obs.loc[test,'EC_family'].astype(str).to_numpy()
            for j,program in enumerate(PROGRAMS):
                true=score.loc[test,program].to_numpy();estimate=pred[:,j]
                tr=true.copy();pr=estimate.copy()
                for fam in set(family):
                    m=family==fam;tr[m]-=tr[m].mean();pr[m]-=pr[m].mean()
                rho=corr(true,estimate);within=corr(tr,pr)
                r2=float(r2_score(true,estimate)) if np.std(true)>1e-10 else np.nan
                passed=np.isfinite(rho) and rho>=.2 and np.isfinite(within) and within>=.15 and r2>0
                rows.append({'stage':stage,'query_library':lib,'program':program,'training_nuclei':int(train.sum()),'query_nuclei':int(test.sum()),
                    'withheld_all_program_genes':True,'fitting_genes':len(fitting),'Spearman':rho,'within_family_residual_Spearman':within,
                    'R2_vs_query_mean':r2,'technical_gate_pass':passed})
    result=pd.DataFrame(rows)
    save(result,DEST/'RNA/withheld_program_opposite_library_benchmark.tsv')
    return result

def main():
    folder=DEST/'RNA';folder.mkdir(parents=True,exist_ok=True)
    data=load_primary();assert len(data)==4597
    genes=data.var_names.to_numpy(dtype=str)
    raw=data.layers['counts_raw'].tocsr();corrected=data.X.tocsr()
    logs=lognorm(raw)
    # Equal library weight for matching controls; no stage pattern enters selection.
    means=[];detect=[]
    for lib in sorted(data.obs['sample'].unique()):
        m=data.obs['sample'].eq(lib).to_numpy()
        means.append(np.asarray(logs[m].mean(axis=0)).ravel());detect.append(np.asarray((raw[m]>=1).mean(axis=0)).ravel())
    controls=select_controls(genes,np.mean(means,axis=0),np.mean(detect,axis=0))
    write_json(folder/'rank_matched_controls.json',controls)
    background=background_audit(genes)
    background=background[~background.lineage.isin(['Endothelial','Unresolved/mixed'])]
    bg=background.groupby(['sample','gene'],observed=True)[['raw_gene_molecules','raw_total_UMI']].sum()
    # Gene-level expected mass and threshold sensitivity; no raw count is restored.
    audit=[];summaries=[];coherence=[];scores={};cycle_sensitivity=[]
    cell=data.obs[['age','sample','EC_family','EC_grouping_support','decontx_contamination','total_counts','raw_total_counts']].copy()
    for layer,counts in [('raw',raw),('DecontX',corrected)]:
        for lib in sorted(data.obs['sample'].unique()):
            m=data.obs['sample'].eq(lib).to_numpy();o=data.obs.loc[m]
            c=counts[m][:,np.isin(genes,CYCLING)].toarray()
            observed_raw=raw[m][:,np.isin(genes,CYCLING)].toarray()>=1
            for threshold in ([1] if layer=='raw' else [.1,.25,.5,1]):
                positive=((c>=threshold)&observed_raw).sum(axis=1)>=2
                cycle_sensitivity.append({'stage':str(o.age.iloc[0]),'sample':lib,'layer':layer,'minimum_expected_count':threshold,
                    'nuclei':int(m.sum()),'cycling_support_nuclei':int(positive.sum()),'cycling_support_fraction':float(positive.mean()),
                    'raw_observation_required':True,'limit':'Detection-threshold sensitivity; not independently calibrated proliferation fraction'})
        print('RANK_SCORES',layer,flush=True)
        for maxrank in [1500,3000]:
            score=rank_excess(counts,genes,controls,maxrank)
            score.index=data.obs_names
            if maxrank==1500:scores[layer]=score
            for program in PROGRAMS:cell[f'{layer}__rank{maxrank}__{program}']=score[program]
            for lib in sorted(data.obs['sample'].unique()):
                m=data.obs['sample'].eq(lib).to_numpy()
                for program in PROGRAMS:
                    summaries.append({'stage':str(data.obs.loc[m,'age'].iloc[0]),'sample':lib,'layer':layer,'maxrank':maxrank,
                        'program':program,'nuclei':int(m.sum()),'mean_rank_excess':float(score.loc[m,program].mean()),
                        'positive_rank_excess_fraction':float(score.loc[m,program].gt(0).mean())})
        norm=lognorm(counts)
        for lib in sorted(data.obs['sample'].unique()):
            m=data.obs['sample'].eq(lib).to_numpy();o=data.obs.loc[m]
            family=pd.get_dummies(o.EC_family,drop_first=True,dtype=float).to_numpy()
            cov=np.column_stack([np.log1p(np.asarray(counts[m].sum(axis=1)).ravel()),o.decontx_contamination.to_numpy(),family])
            for name,panel in PROGRAMS.items():
                js=data.var_names.get_indexer([g for g in panel if g in data.var_names])
                v=norm[m][:,js].toarray()
                eligible=np.asarray((counts[m][:,js]>=.1).sum(axis=0)).ravel()>=5
                v=residual(v,cov)
                rs=[]
                for j in np.flatnonzero(eligible):
                    other=np.flatnonzero(eligible & (np.arange(len(js))!=j))
                    if len(other):rs.append(corr(v[:,j],v[:,other].mean(axis=1)))
                coherence.append({'stage':str(o.age.iloc[0]),'sample':lib,'layer':layer,'program':name,'genes_ge5_nuclei_at_count0p1':int(eligible.sum()),
                    'median_leave_one_gene_out_partial_Spearman':float(np.nanmedian(rs)) if len(rs) and np.isfinite(rs).any() else np.nan,
                    'covariates':'log counts, estimated ambient fraction, five-type/unresolved source family indicators',
                    'limit':'Descriptive gene coherence, not an independent biological significance test'})
    for lib in sorted(data.obs['sample'].unique()):
        m=data.obs['sample'].eq(lib).to_numpy();o=data.obs.loc[m]
        exposure=float((o.raw_total_counts*o.decontx_contamination).sum())
        total_raw=float(o.raw_total_counts.sum());total_corr=float(o.total_counts.sum())
        for g in TARGET_GENES:
            if g not in genes:continue
            j=data.var_names.get_loc(g);r=raw[m,j].toarray().ravel();c=corrected[m,j].toarray().ravel()
            b=bg.loc[(lib,g)];rate=b.raw_gene_molecules/max(b.raw_total_UMI,1)
            expected=exposure*rate
            ec_rate=(r.sum()+.5)/(total_raw+.5*len(genes));bg_rate=(b.raw_gene_molecules+.5)/(b.raw_total_UMI+.5*len(genes))
            row={'stage':str(o.age.iloc[0]),'sample':lib,'gene':g,'nuclei':int(m.sum()),'raw_molecules':float(r.sum()),
                'corrected_expected_mass':float(c.sum()),'fraction_raw_mass_retained':float(c.sum()/r.sum()) if r.sum()>0 else np.nan,
                'raw_positive_nuclei':int((r>=1).sum()),'corrected_positive_ge0p1':int((c>=.1).sum()),'corrected_positive_ge0p25':int((c>=.25).sum()),
                'corrected_positive_ge0p5':int((c>=.5).sum()),'corrected_positive_ge1':int((c>=1).sum()),
                'raw_EC_vs_nonEC_enrichment':ec_rate/bg_rate,'cellular_background_expected_molecules':expected,
                'raw_excess_under_2x_cellular_background':max(0,float(r.sum()-2*expected)),
                'corrected_pseudobulk_CPM':float(c.sum()*1e6/max(total_corr,1)),
                'primary_pseudobulk_support':bool(c.sum()>=5 and (c>=.25).sum()>=3)}
            row['correction_sensitive_EC_enriched_candidate']=bool(r.sum()>=10 and row['raw_EC_vs_nonEC_enrichment']>=2 and row['fraction_raw_mass_retained']<.5 and row['raw_excess_under_2x_cellular_background']>=2*c.sum())
            audit.append(row)
    audit=pd.DataFrame(audit);save(audit,folder/'fractional_count_and_background_audit.tsv')
    save(pd.DataFrame(summaries),folder/'rank_score_library_summary.tsv')
    save(pd.DataFrame(coherence),folder/'within_library_gene_coherence.tsv')
    save(pd.DataFrame(cycle_sensitivity),folder/'cycling_fractional_threshold_sensitivity.tsv')
    cell.to_csv(folder/'EC_fractional_aware_program_scores.tsv.gz',sep='\t',na_rep='NA')
    accepted=[]
    for (stage,gene),sub in audit.groupby(['stage','gene'],observed=True):
        accepted.append({'stage':stage,'gene':gene,'libraries':len(sub),'pseudobulk_support_both_libraries':bool(len(sub)==2 and sub.primary_pseudobulk_support.all()),
            'correction_sensitive_EC_enriched_both_libraries':bool(len(sub)==2 and sub.correction_sensitive_EC_enriched_candidate.all()),
            'corrected_expected_mass_total':sub.corrected_expected_mass.sum(),'raw_molecules_total':sub.raw_molecules.sum(),
            'corrected_positive_ge1_total':sub.corrected_positive_ge1.sum(),'corrected_positive_ge0p25_total':sub.corrected_positive_ge0p25.sum()})
    save(pd.DataFrame(accepted),folder/'gene_support_across_two_libraries.tsv')
    # Assess formerly reviewed ECs as a separate sensitivity; never silently add
    # them to the frozen primary lineage/reference or overwrite subtype labels.
    review_rows=[]
    for stage in STAGES:
        a=ad.read_h5ad(OUT/'snRNA'/stage/'EC_review_corrected_and_raw_counts.h5ad')
        gs=a.var_names.to_numpy(dtype=str);r=a.layers['counts_raw'];c=a.X
        ec=np.isin(gs,LINEAGES['Endothelial']);ep=(c[:,ec]>=.25).multiply(r[:,ec]>=1)
        ec_n=np.asarray(ep.sum(axis=1)).ravel();ec_fraction=np.asarray(c[:,ec].sum(axis=1)).ravel()/np.maximum(np.asarray(c.sum(axis=1)).ravel(),1)
        foreign=np.zeros(len(a),bool)
        for name,panel in LINEAGES.items():
            if name=='Endothelial':continue
            ix=np.isin(gs,panel);v=c[:,ix];present=(v>=.25).multiply(r[:,ix]>=1)
            number=np.asarray(present.sum(axis=1)).ravel();mass=np.asarray(v.sum(axis=1)).ravel()
            foreign|=(number>=3)&(mass>=10)&(mass/np.maximum(np.asarray(c.sum(axis=1)).ravel(),1)>=.02)
        eligible=(ec_n>=2)&(ec_fraction>=.001)&~foreign&~a.obs.coherent_mural_group_review.to_numpy()
        for lib in sorted(a.obs['sample'].unique()):
            m=a.obs['sample'].eq(lib).to_numpy()&eligible
            for g in TARGET_GENES:
                j=a.var_names.get_loc(g)
                review_rows.append({'stage':stage,'sample':lib,'gene':g,'secondary_fractional_aware_EC_review_candidates':int(m.sum()),
                    'raw_molecules':float(r[m,j].sum()),'corrected_expected_mass':float(c[m,j].sum()),
                    'primary_reference_changed':False})
    save(pd.DataFrame(review_rows),folder/'former_QC_review_EC_sensitivity.tsv')
    bench=benchmark(data,scores['DecontX'])
    write_json(folder/'complete.json',{'status':'PASS','nuclei':len(data),'libraries':8,'source_counts_unchanged':True,
        'rank_score':'Mean top1500 enrichment minus20 expression/detection-matched controls per gene; top3000 sensitivity; numerical dust floor0.001; fractional counts retained.',
        'pseudobulk_support':'Corrected expected mass>=5 and>=3 nuclei with corrected count>=0.25, required in both age-matched libraries; exploratory support, not calibrated accuracy.',
        'background':'QC-retained non-EC raw cellular RNA is a sensitivity proxy, not measured empty-droplet ambient. No per-cell count restoration or primary replacement correction.',
        'inference_gate':'Withhold every program gene from shared-panel features; opposite-library same-age Ridge alpha10. Both directions must achieve Spearman>=.2, within-family residual Spearman>=.15 and R2>0 before any inference can be considered.',
        'inference_feature_scaling':'Training-library normalization; SD floor0.25 and clipping[-5,5] prevent unstable amplification of almost-constant fractional features.',
        'programs':PROGRAMS,'program_limits':LIMITS,
        'technical_benchmark_passes':int(bench.technical_gate_pass.sum())})
    print(pd.DataFrame(accepted).query("gene in ['Apln','Esm1','Angpt2','Cxcr4']").to_string(index=False),flush=True)
    print('RNA_ENHANCEMENT_COMPLETE',flush=True)

if __name__=='__main__':main()
