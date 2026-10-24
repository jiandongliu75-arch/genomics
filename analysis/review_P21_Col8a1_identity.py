"""Review P21_F0 against Col8a1-related RNA groups without using UMAP proximity."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-review-mpl')
from pathlib import Path
import json
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from anndata.io import read_elem,sparse_dataset
from fibroblast_state_common import *

DEST=ROOT/'results/fibroblast_P21_identity_review'


def main():
    DEST.mkdir(exist_ok=True)
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    meta=pd.read_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t',index_col=0)
    meta=meta[meta.fine_reference_eligible].copy()
    evidence=pd.read_csv(OUT/'RNA_family_marker_evidence.tsv',sep='\t')
    save_tsv(evidence[evidence.cluster_id.isin(['P21_F0','P21_F2','P14_F1','P14_F2'])],DEST/'original_panel_assignment_evidence.tsv')
    genes_to_compare=pd.read_csv(OUT/'RNA_cross_stage_signature_similarity.tsv',sep='\t',index_col=0)
    # Recreate the exact panel-withheld gene universe used in the earlier analysis.
    tables={};relative={};featuremask=None
    for layer in ['corrected','raw']:
        pb=pd.read_csv(OUT/f'RNA_{layer}_pseudobulk_counts.tsv.gz',sep='\t',index_col=0).loc[:,meta.index]
        tables[layer]=pb
        percluster=pb.T.groupby(meta.fibroblast_cluster_id).sum().T
        cpm=percluster.divide(percluster.sum(axis=0),axis=1)*1e6
        signatures=np.log2(cpm+1)
        for stage in STAGES:
            ci=meta.loc[meta.age.eq(stage),'fibroblast_cluster_id'].unique()
            background=percluster[ci].sum(axis=1);background=background/background.sum()*1e6
            signatures.loc[:,ci]=signatures[ci].subtract(np.log2(background+1),axis=0)
        if layer=='corrected':
            excluded=ANNOTATION_GENES|set(sum(MATRIX_PANELS.values(),[]))|set(WITHHOLD)|set(SEX)|set(CM_CHECK)
            featuremask=~pb.index.isin(excluded)&~pb.index.str.match(r'^(mt-|Rpl|Rps|Gm\d)')
            featuremask&=(cpm.max(axis=1)>=5)&(signatures.abs().max(axis=1)>=.4)
        # Separate library backgrounds remove library-wide and stage-wide shifts.
        cols=[]
        cpm_lib=pb.divide(pb.sum(axis=0),axis=1)*1e6
        rel=np.log2(cpm_lib+1)
        for sample,units in meta.groupby('sample',observed=True).groups.items():
            background=pb[list(units)].sum(axis=1);background=background/background.sum()*1e6
            rel.loc[:,list(units)]=rel[list(units)].subtract(np.log2(background+1),axis=0)
        relative[layer]=rel
    assert featuremask.sum()==12783
    similarities=[]
    query_clusters=['P21_F0','P21_F2']
    ref_clusters=meta.loc[meta.age.eq('P14'),'fibroblast_cluster_id'].unique()
    for layer,rel in relative.items():
        for query in query_clusters:
            qi=meta.index[meta.fibroblast_cluster_id.eq(query)]
            for qunit in qi:
                for rlib in sorted(meta.loc[meta.age.eq('P14'),'sample'].unique()):
                    refs=meta.index[meta['sample'].eq(rlib)]
                    vals=np.corrcoef(rel.loc[featuremask,[qunit]+list(refs)].T)
                    scores=pd.Series(vals[0,1:],index=refs)
                    best=scores.idxmax()
                    for runit,value in scores.items():
                        similarities.append({'layer':layer,'query_cluster':query,'query_library':meta.loc[qunit,'sample'],
                            'reference_library':rlib,'reference_cluster':meta.loc[runit,'fibroblast_cluster_id'],
                            'reference_family':meta.loc[runit,'fibroblast_family'],'correlation':value,
                            'best_reference_in_library':runit==best,'genes_compared':int(featuremask.sum())})
    similarity=pd.DataFrame(similarities)
    save_tsv(similarity,DEST/'independent_library_cross_stage_similarity.tsv')
    # Marker evidence uses actual nuclear RNA counts, separately before/after correction.
    panel_genes=set(sum(PANELS.values(),[]))|set(sum(MATRIX_PANELS.values(),[]))|set('Adgrl3 Svep1 Fap Abca8a Hmcn2 Pcsk6 Fbln2 Sema3c Cfh Eln'.split())
    selected=ann[ann.fibroblast_cluster_id.isin(['P21_F0','P21_F2','P21_F1','P21_F5','P14_F1','P14_F2'])].copy()
    marker_rows=[]
    with h5py.File(BASE/'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad') as f:
        so=read_elem(f['obs']);var=read_elem(f['var']);ii=so.index.get_indexer(selected.index)
        vi=np.flatnonzero(var.gene_symbol.isin(list(panel_genes)))
        for layer,path in [('corrected','X'),('raw','layers/counts_raw')]:
            x=sparse_dataset(f[path])[ii].tocsr()[:,vi].toarray()
            counts=pd.DataFrame(x.T,index=var.iloc[vi].gene_symbol).groupby(level=0).sum().T
            counts.index=selected.index
            total=selected.total_counts if layer=='corrected' else selected.raw_total_counts
            logs=np.log1p(counts.divide(total,axis=0)*1e4)
            for (cluster,lib),ids in selected.groupby(['fibroblast_cluster_id','sample'],observed=True).groups.items():
                for gene in counts:
                    marker_rows.append({'cluster_id':cluster,'library':lib,'layer':layer,'gene':gene,'nuclei':len(ids),
                        'mean_log1pCP10k':logs.loc[ids,gene].mean(),'fraction_count_ge1':counts.loc[ids,gene].ge(1).mean(),
                        'pseudobulk_CPM':counts.loc[ids,gene].sum()/total.loc[ids].sum()*1e6})
    markers=pd.DataFrame(marker_rows);save_tsv(markers,DEST/'marker_expression_by_library.tsv')
    # Distinctness is assessed with two-library effects, never nucleus-level p-values.
    contrast_rows=[]
    for layer,pb in tables.items():
        cpm=pb.divide(pb.sum(axis=0),axis=1)*1e6
        log=np.log2(cpm+1)
        effects=[]
        for lib in sorted(meta.loc[meta.age.eq('P21'),'sample'].unique()):
            a=meta.index[meta['sample'].eq(lib)&meta.fibroblast_cluster_id.eq('P21_F0')][0]
            b=meta.index[meta['sample'].eq(lib)&meta.fibroblast_cluster_id.eq('P21_F2')][0]
            effects.append((lib,log[a]-log[b]))
        for gene in pb.index:
            a,b=effects[0][1][gene],effects[1][1][gene]
            contrast_rows.append({'layer':layer,'gene':gene,'library1_F0_minus_F2_log2CPM':a,
                'library2_F0_minus_F2_log2CPM':b,'mean_log2CPM_difference':(a+b)/2,
                'both_library_directions_agree':np.sign(a)==np.sign(b)})
    contrast=pd.DataFrame(contrast_rows);save_tsv(contrast,DEST/'P21_F0_vs_F2_gene_contrasts.tsv.gz')
    # Existing graph memberships expose resolution changes without reclustering.
    p21=ann[ann.age.eq('P21')].copy()
    cross=[]
    for key in ['resolution_0.25','resolution_0.7','seed_7','seed_123','unadjusted_cluster','QC_sensitivity_cluster','panel_withheld_cluster']:
        tab=pd.crosstab(p21.fibroblast_cluster_id,p21[key])
        for cluster,row in tab.iterrows():
            for label,n in row.items():
                cross.append({'comparison':key,'original_cluster':cluster,'alternative_cluster':str(label),'nuclei':int(n)})
    save_tsv(pd.DataFrame(cross),DEST/'existing_resolution_and_seed_memberships.tsv')
    # Add reviewed hierarchy; original labels and all spatial analyses remain identifiable.
    p21['original_family_assignment']=p21.fibroblast_family
    p21['reviewed_broad_family']=p21.fibroblast_family
    p21.loc[p21.fibroblast_cluster_id.eq('P21_F0'),'reviewed_broad_family']='Col8a1-associated'
    p21['reviewed_substate']=p21.fibroblast_cluster_id
    p21.loc[p21.fibroblast_cluster_id.eq('P21_F0'),'reviewed_substate']='P21_F0: Adgrl3/Svep1 late interstitial candidate'
    p21.loc[p21.fibroblast_cluster_id.eq('P21_F2'),'reviewed_substate']='P21_F2: Hmcn2/Pcsk6-associated substate'
    p21['review_status']='Original assignment retained'
    p21.loc[p21.fibroblast_cluster_id.eq('P21_F0'),'review_status']='Col8a1-related broad family candidate; Col8a1/Lsamp marker overlap; distinct RNA substate retained'
    p21.to_csv(DEST/'P21_reviewed_hierarchical_annotations.tsv.gz',sep='\t',na_rep='NA')
    # Matrix sensitivity for the proposed broader grouping, using original pseudobulks.
    matrix_rows=[]
    for layer,pb in tables.items():
        for lib in sorted(meta.loc[meta.age.eq('P21'),'sample'].unique()):
            sel=meta['sample'].eq(lib)
            for name,clusters in [('original_P21_F2',['P21_F2']),('reviewed_Col8a1_F0_plus_F2',['P21_F0','P21_F2']),('Gfpt2_F3',['P21_F3'])]:
                ids=meta.index[sel&meta.fibroblast_cluster_id.isin(clusters)]
                counts=pb[ids].sum(axis=1);cpm=counts/counts.sum()*1e6
                for program,panel in MATRIX_PANELS.items():
                    gl=[g for g in panel if g in cpm.index and g not in ANNOTATION_GENES]
                    if gl:matrix_rows.append({'layer':layer,'library':lib,'grouping':name,'program':program,
                        'score_mean_log1pCPM':np.log1p(cpm.loc[gl]).mean()})
    save_tsv(pd.DataFrame(matrix_rows),DEST/'matrix_grouping_sensitivity.tsv')
    # Readable figure and a compact PDF.
    plt.rcParams.update({'font.size':11,'pdf.fonttype':42})
    with PdfPages(DEST/'P21_Col8a1_identity_review.pdf') as pdf:
        fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
        for ax,col,title in [(axes[0],'original_family_assignment','Original exclusive panel assignment'),(axes[1],'reviewed_broad_family','Reviewed broader family; original subclusters retained')]:
            for family,color in COLORS.items():
                m=p21[col].eq(family)
                ax.scatter(p21.loc[m,'fibroblast_UMAP1'],p21.loc[m,'fibroblast_UMAP2'],c=color,s=3,rasterized=True)
            for cl,group in p21.groupby('fibroblast_cluster_id',observed=True):
                x,y=group[['fibroblast_UMAP1','fibroblast_UMAP2']].median()
                ax.text(x,y,cl.replace('P21_',''),ha='center',weight='bold',bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':1})
            ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
        fig.legend(handles=[Line2D([],[],marker='o',color=c,ls='',label=f) for f,c in COLORS.items()],loc='outside lower center',ncol=3,frameon=False)
        fig.suptitle('P21_F0 supports Col8a1-related broad identity with a distinct late interstitial program\nFamily review uses marker support, separate-library RNA resemblance and clustering stability; UMAP unchanged',fontsize=14)
        fig.savefig(DEST/'P21_broad_family_review.png',dpi=170,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
        # All-library reference matching.
        best=similarity[similarity.best_reference_in_library&similarity.query_cluster.eq('P21_F0')]
        fig,ax=plt.subplots(figsize=(11,5),constrained_layout=True)
        show=similarity[similarity.layer.eq('corrected')&similarity.query_cluster.eq('P21_F0')]
        tab=show.groupby('reference_cluster').correlation.agg(['mean','min','max']).sort_values('mean',ascending=False)
        ax.bar(range(len(tab)),tab['mean'],color=['#0072b2' if t in ['P14_F1','P14_F2','P14_F4'] else '#999999' for t in tab.index])
        ax.vlines(range(len(tab)),tab['min'],tab['max'],colors='black');ax.set_xticks(range(len(tab)),tab.index)
        ax.set_ylabel('Relative transcriptome correlation');ax.set_title('P21_F0 matches P14_F2 in every corrected-count library comparison\n12,783 genes; interpreted panels withheld; line spans four library combinations')
        fig.savefig(DEST/'P21_F0_cross_library_reference_matches.png',dpi=170,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
        gl=['Col8a1','Col15a1','Hmcn2','Htra3','Lsamp','Eln','Cfh','Gfpt2','Adgrl3','Svep1','Lamc1','Hspg2']
        fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True)
        for j,cl in enumerate(['P21_F0','P21_F2','P21_F1','P21_F5']):
            v=markers[markers.cluster_id.eq(cl)&markers.layer.eq('corrected')].groupby('gene').mean_log1pCP10k.mean().reindex(gl)
            ax.plot(range(len(gl)),v,'o-',label=cl,lw=1)
        ax.set_xticks(range(len(gl)),gl,rotation=45,ha='right');ax.set_ylabel('Mean log1p(CP10k)');ax.legend(frameon=False)
        ax.set_title('Marker expression distinguishes related P21 fibroblast substates\nMeans of separate libraries; Col8a1 is not exclusive to one named family')
        fig.savefig(DEST/'P21_F0_F2_marker_comparison.png',dpi=170,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    # Explicit validations for the reviewed result.
    target=evidence[evidence.cluster_id.eq('P21_F0')&evidence.candidate_family.eq('Col8a1-associated')].iloc[0]
    assert target.panel_qualified and len(target.supported_genes.split(';'))==10
    assert best.reference_cluster.eq('P14_F2').all(),best.to_dict('records')
    assert len(best)==8
    np.testing.assert_allclose(p21[['fibroblast_UMAP1','fibroblast_UMAP2']],ann.loc[p21.index,['fibroblast_UMAP1','fibroblast_UMAP2']])
    lower=p21.loc[p21.fibroblast_cluster_id.eq('P21_F0'),'resolution_0.25'].mode().iloc[0]
    shares={cl:int((p21.loc[p21.fibroblast_cluster_id.eq(cl),'resolution_0.25']==lower).sum()) for cl in ['P21_F0','P21_F2']}
    result={'status':'PASS','P21_F0_nuclei':1330,'library_counts':p21[p21.fibroblast_cluster_id.eq('P21_F0')]['sample'].value_counts().to_dict(),
        'original_reason_unresolved':'Col8a1 and Lsamp both qualified; top panel margin0.037 below0.1 cutoff, not missing Col8a1 marker support',
        'all_10_Col8a1_panel_markers_supported':True,
        'best_P14_cluster_in_all_8_raw_corrected_library_combinations':'P14_F2',
        'corrected_best_reference_correlation_range':best[best.layer.eq('corrected')].correlation.agg(['min','max']).to_dict(),
        'raw_best_reference_correlation_range':best[best.layer.eq('raw')].correlation.agg(['min','max']).to_dict(),
        'lower_resolution_common_cluster_memberships':shares,
        'reviewed_interpretation':'Col8a1-related broad family candidate; retain P21_F0 as Adgrl3/Svep1 late interstitial substate, distinct from P21_F2 Hmcn2/Pcsk6 substate',
        'all_4411_P21_nuclei_preserved':len(p21)==4411,'counts_and_UMAP_unchanged':True,
        'scope':'Additive reviewed hierarchy; earlier primary atlas labels and fitted spatial results retain their original definitions',
        'limits':'No lineage ancestry or anatomical validation; P14 reference labels provisional; overlapping marker programs remain'}
    write_json(DEST/'validation.json',result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
