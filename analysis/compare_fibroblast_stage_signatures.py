"""Compare whole-transcriptome cluster signatures without interpreted panels."""
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from anndata.io import read_elem,sparse_dataset
from fibroblast_state_common import *


def main():
    meta=pd.read_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t',index_col=0)
    meta=meta.loc[meta.fine_reference_eligible]
    counts=pd.read_csv(OUT/'RNA_corrected_pseudobulk_counts.tsv.gz',sep='\t',index_col=0).loc[:,meta.index]
    bycluster=counts.T.groupby(meta.fibroblast_cluster_id).sum().T
    cluster_meta=meta.groupby('fibroblast_cluster_id',observed=True).agg(stage=('age','first'),family=('fibroblast_family','first'),nuclei=('nuclei','sum'))
    cpm=bycluster.divide(bycluster.sum(0),axis=1)*1e6
    log=np.log2(cpm+1)
    signatures=log.copy()
    for stage in STAGES:
        ids=cluster_meta.index[cluster_meta.stage.eq(stage)]
        background=bycluster[ids].sum(1);background=background/background.sum()*1e6
        signatures.loc[:,ids]=log[ids].subtract(np.log2(background+1),axis=0)
    excluded=ANNOTATION_GENES|set(sum(MATRIX_PANELS.values(),[]))|set(WITHHOLD)|set(SEX)|set(CM_CHECK)
    eligible=~signatures.index.isin(excluded)&~signatures.index.str.match(r'^(mt-|Rpl|Rps|Gm\d)')
    eligible&=(cpm.max(axis=1)>=5)&(signatures.abs().max(axis=1)>=.4)
    similarity=signatures.loc[eligible].corr()
    similarity.to_csv(OUT/'RNA_cross_stage_signature_similarity.tsv',sep='\t')
    rows=[]
    for a,b in zip(STAGES[:-1],STAGES[1:]):
        aa=cluster_meta.index[cluster_meta.stage.eq(a)];bb=cluster_meta.index[cluster_meta.stage.eq(b)]
        for source in aa:
            target=similarity.loc[source,bb].idxmax()
            reciprocal=similarity.loc[target,aa].idxmax()==source
            rows.append({'from_cluster':source,'from_family':cluster_meta.loc[source,'family'],'to_cluster':target,
                'to_family':cluster_meta.loc[target,'family'],'correlation':similarity.loc[source,target],
                'reciprocal_nearest_match':reciprocal,'genes_compared':int(eligible.sum()),
                'interpretation':'relative RNA resemblance with interpreted genes withheld; not lineage ancestry'})
    save_tsv(pd.DataFrame(rows),OUT/'RNA_adjacent_stage_candidate_matches.tsv')
    # Reassess gene-stage effects after within-stage library, depth, ambient and cycle adjustment.
    effects=pd.read_csv(OUT/'RNA_within_family_stage_gene_changes.tsv',sep='\t')
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    requested=sorted(set(effects.gene))
    with h5py.File(BASE/'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad') as f:
        sourceobs=read_elem(f['obs']);var=read_elem(f['var'])
        ii=sourceobs.index.get_indexer(ann.index);assert (ii>=0).all()
        vi=np.flatnonzero(var.gene_symbol.isin(requested))
        # Read only FB rows from the source; collapse duplicate symbols exactly for these genes.
        x=sparse_dataset(f['X'])[ii].tocsr()[:,vi].toarray()
        gene=var.iloc[vi].gene_symbol.to_numpy()
    vals=pd.DataFrame(x.T,index=gene).groupby(level=0).sum().T
    vals.index=ann.index
    logs=np.log1p(vals.divide(ann.total_counts,axis=0)*1e4)
    adjusted=[]
    for (family,early,late),part in effects.groupby(['family','from_stage','to_stage'],observed=True):
        mask=ann.fine_reference_eligible&ann.fibroblast_family.eq(family)&ann.age.isin([early,late])
        o=ann.loc[mask];age=o.age.eq(late).to_numpy(dtype=float)
        nuisance=[]
        for stage in [early,late]:
            m=o.age.eq(stage);libs=sorted(o.loc[m,'sample'].unique())
            nuisance.append(np.where(m,np.where(o['sample'].eq(libs[-1]),1.,-1.),0.))
        nuisance.extend([o.decontx_contamination.to_numpy(),np.log1p(o.total_counts.to_numpy()),o.cycling_marker_count.to_numpy()])
        cov=np.column_stack(nuisance);cov=(cov-cov.mean(0))/np.maximum(cov.std(0),1e-8)
        design=np.column_stack([np.ones(len(o)),age,cov])
        if np.linalg.matrix_rank(design)<design.shape[1]:
            # Zero cell-cycle columns are removed, not silently included as a rank-deficient fit.
            design=design[:,np.r_[True,True,cov.std(0)>1e-8]]
        gl=part.gene.tolist()
        beta=np.linalg.lstsq(design,logs.loc[o.index,gl].to_numpy(),rcond=None)[0][1]
        for gene,b in zip(gl,beta):
            row=part.loc[part.gene.eq(gene)].iloc[0]
            adjusted.append({'family':family,'from_stage':early,'to_stage':late,'gene':gene,
                'adjusted_log1pCP10k_stage_coefficient':b,'adjusted_direction_agrees_with_pseudobulk':np.sign(b)==np.sign(row.corrected_log2CPM_difference),
                'adjusted_effect_ge_0_05':abs(b)>=.05,'nuclei':len(o),
                'interpretation':'nuisance-adjusted effect sensitivity; no independent-cell significance claim'})
    save_tsv(pd.DataFrame(adjusted),OUT/'RNA_stage_gene_QC_adjustment.tsv')
    write_json(OUT/'cross_stage_methods.json',{'signature_genes':int(eligible.sum()),
        'excluded':'Exact annotation panels,matrix panels,maturation/Ankrd1 panels,sex/CM-check/mt/ribosomal/Gm genes',
        'signature':'Cluster log2CPM minus same-stage pooled FB log2CPM; adjacent-stage Pearson similarity',
        'gene_effect_adjustment':'Original corrected per-nucleus log1pCP10k; stage plus nested within-age library contrasts,ambient,depth,cycling; effect sensitivity only',
        'ancestry_claim':False})
    print(pd.DataFrame(rows).to_string(index=False),flush=True)


if __name__=='__main__':main()
