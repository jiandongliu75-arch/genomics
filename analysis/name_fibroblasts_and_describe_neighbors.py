"""Name P21 RNA substates and describe confidence-qualified spatial families.

Names change presentation only. Spatial neighborhoods use existing family-level
calls and frozen broad cardiac identities; no finer RNA label is projected here.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-neighborhood-mpl')
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.preprocessing import StandardScaler
from fibroblast_state_common import ROOT,STAGES,FAMILIES,UNRESOLVED,COLORS

OUT=ROOT/'results/fibroblast_names_and_neighborhoods'
CONS=ROOT/'results/fibroblast_spatial_consensus'
ATLAS=ROOT/'results/whole_heart_spatial_atlas'
SEG=ROOT/'results/xenium_nuclear_anchored_segmentation'
LINEAGES=['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
RADII=[30,50,100]
ALIASES={
 'P21_F0':('Col8a1/Svep1 fibroblasts','Col8a1/Svep1','Col8a1-associated','Adgrl3;Svep1;Fap;Abca8a','Matrix organization/adhesion-associated RNA'),
 'P21_F2':('Col8a1/Hmcn2 fibroblasts','Col8a1/Hmcn2','Col8a1-associated','Hmcn2;Pcsk6;Adamtsl3;Fbln2;Lamc1','Matrix organization/processing-associated RNA'),
 'P21_F1':('Lsamp/Slit3 fibroblasts','Lsamp/Slit3','Lsamp-associated','Slit3;Rgs6;Chrm3;Cfh;Eln;Angpt1','Elastic-matrix and adhesion/signaling-associated RNA'),
 'P21_F5':('Lsamp/Col4a5 fibroblasts','Lsamp/Col4a5','Lsamp-associated','Lsamp;Col4a5;Col4a6;Itgbl1;Eln;Cfh','Elastic-matrix and collagen-IV-associated RNA'),
 'P21_F3':('Gfpt2-associated fibroblasts','Gfpt2','Gfpt2-associated','Gfpt2;Adgrd1;Cmah;Sema3c','Matrix/metabolic-associated RNA'),
 'P21_F4':('Kirrel3-associated fibroblasts','Kirrel3','Kirrel3-associated','Kirrel3;Erbb4;Pdzrn4;Dkk3','Adhesion/signaling-associated RNA; valve anatomy unverified'),
 'P21_F6':('Residual CM-program review','CM-program review',UNRESOLVED,'Fhl2;Rbm20;Fgf13;Cacna1c','Excluded from clean fibroblast reference; source retained')}


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')


def name_rna():
    a=pd.read_csv(ROOT/'results/fibroblast_Xenium_reviewed_refit/RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    a=a.loc[a.age.eq('P21')].copy();assert len(a)==4411
    rows=[]
    for cid,(name,short,family,markers,function) in ALIASES.items():
        m=a.fibroblast_cluster_id.eq(cid);assert m.any()
        a.loc[m,'display_name']=name;a.loc[m,'display_short_name']=short
        assert a.loc[m,'fibroblast_family'].eq(family).all()
        rows.append({'original_cluster_id':cid,'display_name':name,'short_name':short,'broad_family':family,
            'distinguishing_RNA_markers':markers,'functional_hypothesis':function,'nuclei':int(m.sum()),
            'spatial_resolution':'Broad family only; separate substate location unestablished' if cid not in ['P21_F6'] else 'Withheld from fine reference'})
    pd.DataFrame(rows).to_csv(OUT/'P21_fibroblast_name_key.tsv',sep='\t',index=False)
    a['original_clustering_support']=a.fibroblast_support
    a['identity_review_note']='Provisional RNA substate; separate anatomical location unestablished'
    a.loc[a.fibroblast_cluster_id.eq('P21_F0'),'identity_review_note']='Focused cross-library review supports Col8a1-related family; original pre-review support field retained separately'
    a.loc[a.fibroblast_cluster_id.eq('P21_F6'),'identity_review_note']='Residual CM-program review; excluded from clean fibroblast reference'
    columns=['age','sample','fibroblast_cluster_id','display_name','display_short_name','fibroblast_family',
        'original_clustering_support','identity_review_note','fine_reference_eligible','lineage_mixing_review_flag','fibroblast_UMAP1','fibroblast_UMAP2']
    a[columns].to_csv(OUT/'P21_RNA_named_annotations.tsv.gz',sep='\t',na_rep='NA')
    return a


def context(stage):
    c=pd.read_csv(CONS/stage/'fibroblast_consensus_annotations.tsv.gz',sep='\t',index_col=0)
    atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,
        usecols=['cell_id','x_centroid','y_centroid','atlas_identity','atlas_resolved','cell_area'])
    np.testing.assert_allclose(c[['x_centroid','y_centroid']],atlas.loc[c.index,['x_centroid','y_centroid']],rtol=0,atol=1e-7)
    assert atlas.loc[c.index,'atlas_identity'].eq('Fibroblast').all()
    resolved=atlas.loc[atlas.atlas_resolved]
    assert set(resolved.atlas_identity)<=set(LINEAGES)
    xy=c[['x_centroid','y_centroid']].to_numpy()
    whole=cKDTree(resolved[['x_centroid','y_centroid']].to_numpy())
    distances=pd.read_csv(ROOT/'results/whole_heart_Ankrd1_vessels'/stage/'arterial_boundary_distances.tsv.gz',sep='\t',index_col=0)
    primary=pd.read_csv(SEG/stage/'fibroblast_nuclear_anchored_assignments.tsv.gz',sep='\t',index_col=0)
    nuc=pd.read_csv(SEG/stage/'fibroblast_nuclear_only_assignments.tsv.gz',sep='\t',index_col=0)
    c['cell_area']=atlas.loc[c.index,'cell_area']
    c['gene_molecules']=primary.total_gene_molecules
    c['detected_genes']=primary.detected_genes
    c['nuclear_fraction']=nuc.total_gene_molecules/np.maximum(primary.total_gene_molecules,1)
    c['arterial_distance_um']=distances.loc[c.index,'independent_strict__boundary_distance_um']
    for radius in RADII:
        c[f'arterial_within_{radius}um']=np.where(c.arterial_distance_um.notna(),c.arterial_distance_um.le(radius).astype(float),np.nan)
    for radius in RADII:
        c[f'resolved_neighbors_{radius}um']=whole.query_ball_point(xy,radius,return_length=True)-1
    for lineage in LINEAGES:
        tree=cKDTree(resolved.loc[resolved.atlas_identity.eq(lineage),['x_centroid','y_centroid']].to_numpy())
        c['nearest_'+lineage+'_um']=tree.query(xy,k=2 if lineage=='Fibroblast' else 1)[0][:,1] if lineage=='Fibroblast' else tree.query(xy)[0]
        for radius in RADII:
            num=tree.query_ball_point(xy,radius,return_length=True)
            if lineage=='Fibroblast':num-=1
            c[f'n_{lineage}_{radius}um']=num
            c[f'fraction_{lineage}_{radius}um']=num/np.maximum(c[f'resolved_neighbors_{radius}um'],1)
    # Direct physical recomputation verifies the inherited 30um context, including
    # the resolved-neighbor denominator and self exclusion.
    old=pd.read_csv(ROOT/'results/fibroblast_states'/stage/'FB_broad_neighbor_context_30um.tsv.gz',sep='\t',index_col=0)
    np.testing.assert_array_equal(c.resolved_neighbors_30um,old.loc[c.index,'resolved_other_neighbors'])
    for lineage in LINEAGES:
        np.testing.assert_allclose(c[f'fraction_{lineage}_30um'],old.loc[c.index,'neighbor_fraction__'+lineage],atol=1e-7)
    for radius in RADII:
        assert (c[[f'n_{l}_{radius}um' for l in LINEAGES]].sum(1)==c[f'resolved_neighbors_{radius}um']).all()
    c.to_csv(OUT/stage/'FB_physical_neighbors.tsv.gz',sep='\t',na_rep='NA')
    return c


def summarize(stage,c):
    rows=[]
    for scope,label in [('A+B','high_confidence_family'),('A+B+C','assigned_family')]:
        for radius in RADII:
            eligible=c[f'resolved_neighbors_{radius}um'].ge(5)
            measures=[f'fraction_{l}_{radius}um' for l in LINEAGES]+['arterial_distance_um','arterial_within_30um','arterial_within_50um','arterial_within_100um']+[f'nearest_{l}_um' for l in LINEAGES]
            for family in FAMILIES:
                selected=eligible&c[label].eq(family)
                if not selected.any():continue
                for measure in measures:
                    values=c.loc[selected,measure].replace([np.inf,-np.inf],np.nan).dropna()
                    base=c.loc[eligible,measure].replace([np.inf,-np.inf],np.nan).dropna()
                    rows.append({'stage':stage,'scope':scope,'radius_um':radius,'family':family,
                        'assigned_profiles':int(c[label].eq(family).sum()),'eligible_profiles':int(selected.sum()),
                        'sparse_group_lt30':int(selected.sum())<30,'measure':measure,
                        'mean':values.mean(),'median':values.median(),'all_FB_mean':base.mean(),
                        'difference_from_all_FB':values.mean()-base.mean(),
                        'interpretation':'Descriptive same-section neighborhood, not cell contact, abundance, signaling or independent-animal evidence'})
    return rows


def matched_context(stage,c):
    # Predefined contrasts with adequate high-confidence groups in this dataset.
    comparisons={'P7':('Gfpt2-associated','Postn-associated'),'P14':('Col8a1-associated','Gfpt2-associated')}
    if stage not in comparisons:return [],[],[]
    target,reference=comparisons[stage]
    valid=c.resolved_neighbors_30um.ge(5)
    t=c.loc[valid&c.high_confidence_family.eq(target)].copy()
    r=c.loc[valid&c.high_confidence_family.eq(reference)].copy()
    cov=pd.DataFrame(index=c.index)
    for col in ['gene_molecules','cell_area','resolved_neighbors_50um']:cov['log_'+col]=np.log1p(c[col])
    cov['nuclear_fraction']=c.nuclear_fraction
    scaler=StandardScaler().fit(cov.loc[t.index.union(r.index)])
    tv=scaler.transform(cov.loc[t.index]);rv=scaler.transform(cov.loc[r.index])
    distance,ix=cKDTree(rv).query(tv)
    keep=distance<=2.;tt=t.iloc[np.flatnonzero(keep)];rr=r.iloc[ix[keep]]
    pairs=[{'stage':stage,'target_family':target,'reference_family':reference,'target_cell_id':i,'reference_cell_id':j,
            'standardized_covariate_distance':d} for i,j,d in zip(tt.index,rr.index,distance[keep])]
    balance=[]
    for k,name in enumerate(cov.columns):
        balance.append({'stage':stage,'covariate':name,'target_family':target,'reference_family':reference,
            'before_standardized_mean_difference':tv[:,k].mean()-rv[:,k].mean(),
            'after_standardized_mean_difference':tv[keep,k].mean()-rv[ix[keep],k].mean()})
    rows=[]
    for radius in RADII:
        for lineage in LINEAGES:
            col=f'fraction_{lineage}_{radius}um'
            rows.append({'stage':stage,'radius_um':radius,'measure':col,'target_family':target,'reference_family':reference,
                'target_profiles_before':len(t),'reference_profiles_before':len(r),'matched_pairs':len(tt),
                'unique_reference_profiles':rr.index.nunique(),'target_matched_mean':tt[col].mean(),
                'reference_matched_mean':rr[col].mean(),'matched_mean_difference':tt[col].mean()-rr[col].mean(),
                'unmatched_mean_difference':t[col].mean()-r[col].mean(),
                'interpretation':'Nearest matching with replacement and caliper2SD, on RNA depth, area, nuclear fraction and density; descriptive sensitivity, not causal adjustment or independent replication'})
    for col in ['arterial_within_30um','arterial_within_50um','arterial_within_100um']:
        rows.append({'stage':stage,'radius_um':30,'measure':col,'target_family':target,'reference_family':reference,
            'target_profiles_before':len(t),'reference_profiles_before':len(r),'matched_pairs':len(tt),
            'unique_reference_profiles':rr.index.nunique(),'target_matched_mean':tt[col].mean(),
            'reference_matched_mean':rr[col].mean(),'matched_mean_difference':tt[col].mean()-rr[col].mean(),
            'unmatched_mean_difference':t[col].mean()-r[col].mean(),
            'interpretation':'Same matched30um-eligible profiles; native independent vascular-envelope distance, not traced histological vessel wall'})
    return rows,pairs,balance


def plots(named,summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    colors=['#0072b2','#56b4e9','#cc79a7','#aa4499','#009e73','#8b5c38']
    selected=['P21_F0','P21_F2','P21_F1','P21_F5','P21_F3','P21_F4']
    with PdfPages(OUT/'named_fibroblasts_and_spatial_context.pdf') as pdf:
        fig,ax=plt.subplots(figsize=(10,7),constrained_layout=True)
        for cid,color in zip(selected,colors):
            d=named.loc[named.fibroblast_cluster_id.eq(cid)]
            ax.scatter(d.fibroblast_UMAP1,d.fibroblast_UMAP2,s=4,c=color,label=ALIASES[cid][1],rasterized=True)
            ax.text(d.fibroblast_UMAP1.median(),d.fibroblast_UMAP2.median(),ALIASES[cid][1].replace('/','/\n'),
                    fontsize=10,ha='center',va='center',bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1.5})
        ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
        ax.set_title('P21 RNA fibroblast substates · descriptive marker names')
        ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False,markerscale=3)
        fig.savefig(OUT/'P21_named_fibroblast_RNA_UMAP.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        s=summary.loc[summary.scope.eq('A+B')&summary.radius_um.eq(30)&~summary.sparse_group_lt30]
        keys=s[['stage','family','eligible_profiles']].drop_duplicates()
        labels=[];values=[]
        for r in keys.itertuples():
            p=s.loc[s.stage.eq(r.stage)&s.family.eq(r.family)].set_index('measure')
            labels.append(f'{r.stage} {r.family.replace("-associated", "")} (n={r.eligible_profiles})')
            values.append([p.loc[f'fraction_{l}_30um','mean'] for l in LINEAGES])
        values=np.array(values)
        fig,ax=plt.subplots(figsize=(12,6),constrained_layout=True)
        im=ax.imshow(values,cmap='YlGnBu',vmin=0,vmax=.65,aspect='auto')
        ax.set_yticks(range(len(labels)),labels)
        ax.set_xticks(range(len(LINEAGES)),['CM','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial'],rotation=30,ha='right')
        for i in range(len(labels)):
            for j in range(len(LINEAGES)):ax.text(j,i,f'{100*values[i,j]:.1f}%',ha='center',va='center',color='white' if values[i,j]>.4 else 'black')
        fig.colorbar(im,ax=ax,shrink=.6,label='Mean share of resolved neighbors')
        ax.set_title('30µm neighborhoods · high-confidence families with ≥30 eligible profiles\nDescriptive same-section context; source profile itself excluded')
        fig.savefig(OUT/'high_confidence_FB_neighbor_composition.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        for stage in STAGES:
            a=pd.read_csv(OUT/stage/'FB_physical_neighbors.tsv.gz',sep='\t',index_col=0)
            atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,
                usecols=['cell_id','x_centroid','y_centroid','atlas_identity','atlas_resolved'])
            fig,axes=plt.subplots(1,2,figsize=(15,7),constrained_layout=True)
            bg={'Cardiomyocyte':'#dddddd','Endothelial':'#a6cee3','Pericyte':'#cab2d6','Smooth muscle':'#6a3d9a','Immune':'#fdbf6f','Epicardial':'#b2df8a'}
            for lineage,color in bg.items():
                p=atlas.loc[atlas.atlas_resolved&atlas.atlas_identity.eq(lineage)]
                axes[0].scatter(p.x_centroid,p.y_centroid,s=.35,c=color,rasterized=True)
            for family in FAMILIES:
                p=a.loc[a.high_confidence_family.eq(family)]
                axes[0].scatter(p.x_centroid,p.y_centroid,s=4,c=COLORS[family],edgecolors='black',linewidths=.1,rasterized=True)
            axes[0].set_title('High-confidence fibroblast families over other cardiac profiles')
            im=axes[1].scatter(a.x_centroid,a.y_centroid,c=a.fraction_Cardiomyocyte_30um,s=2,cmap='viridis',vmin=0,vmax=1,rasterized=True)
            axes[1].set_title('CM share within30µm around all broad fibroblasts')
            fig.colorbar(im,ax=axes[1],shrink=.4,label='Resolved-neighbor CM fraction')
            for ax in axes:
                ax.set_xlim(atlas.x_centroid.min()-100,atlas.x_centroid.max()+100)
                ax.set_ylim(atlas.y_centroid.max()+100,atlas.y_centroid.min()-100)
                ax.set_aspect('equal');ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
            handles=[Line2D([],[],marker='o',ls='',color=c,label=l) for l,c in bg.items()]
            handles += [Line2D([],[],marker='o',ls='',color=COLORS[f],label=f.replace('-associated',' FB family')) for f in FAMILIES]
            fig.legend(handles=handles,loc='outside lower center',ncol=4,frameon=False,fontsize=8)
            fig.suptitle(stage+' · Frozen broad cardiac identities and original physical positions')
            fig.savefig(OUT/f'{stage}_fibroblast_other_cell_context.png',dpi=180);pdf.savefig(fig);plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    named=name_rna();rows=[];matched=[];pairs=[];balance=[];checks=[]
    for stage in STAGES:
        (OUT/stage).mkdir(exist_ok=True)
        c=context(stage);rows+=summarize(stage,c)
        a,b,d=matched_context(stage,c);matched+=a;pairs+=b;balance+=d
        checks.append({'stage':stage,'broad_FB_profiles':len(c),'physical_30um_context_matches_previous_exactly':True,
            'new_50_100um_contexts_recomputed':True,'original_coordinates_unchanged':True})
        print('CONTEXT_DONE',stage,len(c),flush=True)
    summary=pd.DataFrame(rows)
    summary.to_csv(OUT/'family_neighbor_context.tsv',sep='\t',index=False)
    pd.DataFrame(matched).to_csv(OUT/'matched_family_context_sensitivity.tsv',sep='\t',index=False)
    pd.DataFrame(pairs).to_csv(OUT/'RNA_coverage_matched_pairs.tsv.gz',sep='\t',index=False)
    pd.DataFrame(balance).to_csv(OUT/'matching_covariate_balance.tsv',sep='\t',index=False)
    assert sum(c['broad_FB_profiles'] for c in checks)==33790
    write_json(OUT/'validation.json',{'status':'PASS','stages':checks,'P21_RNA_nuclei':len(named),
        'P21_source_ids_and_labels_retained':True,'named_states_not_projected_to_unresolved_spatial_substates':True})
    write_json(OUT/'methods.json',{'names':'Presentation aliases for P21 RNA substates; not a new cell-type definition or universal double-positive requirement',
        'spatial_labels':'Existing consensus A+B primary, A+B+C sensitivity; no new label fit or reallocation',
        'neighbor_radii_um':RADII,'neighbor_denominator':'Resolved native broad cardiac profiles excluding the source fibroblast; require>=5neighbors',
        'nearest_distance':'Other-profile centroid distance; artery uses preceding native EC/SMC-envelope boundary',
        'minimum_group_for_main_interpretation':30,'low_n_groups':'Reported but marked sparse; no firm niche claim',
        'matching':'Predefined P7 Gfpt2 vsPostn and P14 Col8a1 vsGfpt2; nearest with replacement, standardized4covariates,Euclidean caliper2,30um-eligible cohorts',
        'covariates':'log1p gene molecules,log1p native area,log1p resolved50um neighbor count,nuclear/whole-cell molecule fraction',
        'limits':'One section/age; confidence selection and incomplete fine coverage; neighbors not physical contact or causal signaling; cell profiles not independent animal replicates; no substate-specific spatial attribution'})
    plots(named,summary)
    # The narrative is a reviewed static file alongside the reproducible tables.
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='file_manifest.tsv']+[Path(__file__)]
    manifest=[]
    for p in sorted(files):
        with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        manifest.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest})
    pd.DataFrame(manifest).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)


if __name__=='__main__':main()
