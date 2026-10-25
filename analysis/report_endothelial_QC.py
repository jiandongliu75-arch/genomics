"""Readable marker names, RNA/spatial figures, and actual EC-boundary RNA audit."""
from endothelial_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import h5py
import anndata as ad
import scanpy as sc
from anndata.io import read_elem
from scipy import sparse
import gzip
from shapely.geometry import shape,box
from plot_xenium_boundary_qc import outline

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
ATLAS=ROOT/'results/whole_heart_spatial_atlas'
SEG=ROOT/'results/xenium_nuclear_anchored_segmentation'
BOUND=ROOT/'results/xenium_boundary_qc'

def finish(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight')
    pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def points(ax,coords,labels,colors,s=3,legend=True):
    for name,c in colors.items():
        m=np.asarray(labels==name)
        if m.any():ax.scatter(coords[m,0],coords[m,1],s=s,color=c,alpha=.8,edgecolors='none',rasterized=True)
    if legend:
        handles=[Line2D([],[],ls='',marker='o',color=c,label=f'{n} (n={np.sum(labels==n):,})',markersize=4) for n,c in colors.items() if np.sum(labels==n)]
        ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=8)
    ax.set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2')

def textpage(pdf,title,lines):
    fig=plt.figure(figsize=(12,8.5));fig.text(.05,.94,title,fontsize=19,weight='bold',va='top')
    y=.87
    import textwrap
    for line in lines:
        for part in textwrap.wrap(line,width=110) if line else ['']:
            fig.text(.05,y,part,fontsize=11,va='top');y-=.028
        y-=.016
    pdf.savefig(fig);plt.close(fig)

def rnafigures(pdf):
    ann=pd.read_csv(OUT/'snRNA/EC_annotations.tsv.gz',sep='\t',index_col=0)
    summary=pd.read_csv(OUT/'snRNA/EC_cluster_summary.tsv',sep='\t')
    for stage in STAGES:
        a=ann[ann.age.eq(stage)];xy=a[['EC_UMAP1','EC_UMAP2']].to_numpy()
        fig,axes=plt.subplots(1,3,figsize=(20,6.8),layout='constrained')
        rows=summary[summary.stage.eq(stage)]
        palette=dict(zip(rows.cluster_id,plt.get_cmap('tab20').colors))
        points(axes[0],xy,a.EC_cluster_id.to_numpy(),palette,legend=False,s=7)
        for _,r in rows.iterrows():
            m=a.EC_cluster_id.eq(r.cluster_id);c=a.loc[m,['EC_UMAP1','EC_UMAP2']].median().to_numpy()
            # Labels are aliases; source computational IDs remain in tables.
            axes[0].annotate(r.gene_name.replace(' ECs',''),c,xytext=(7,5),textcoords='offset points',fontsize=8,
                bbox={'facecolor':'white','edgecolor':palette[r.cluster_id],'alpha':.88,'pad':2},arrowprops={'arrowstyle':'-','color':'#666666','lw':.5})
        axes[0].set_title('Discovered RNA groups · gene-enrichment names')
        points(axes[1],xy,a.EC_family.to_numpy(),COLORS,s=7);axes[1].set_title('Conservative subtype interpretation')
        libs=sorted(a['sample'].unique());points(axes[2],xy,a['sample'].to_numpy(),dict(zip(libs,['#0072b2','#d55e00'])),s=6)
        axes[2].set_title('RNA libraries · same coordinates')
        fig.suptitle(f'{stage} endothelial nuclei · {len(a):,} primary profiles\nWithin-stage library adjustment; different stages have independently fitted UMAP axes',fontsize=15)
        finish(fig,f'{stage}_snRNA_endothelial_named_UMAP',pdf)
    # One shared, unintegrated RNA embedding displays developmental separation.
    x=sparse.vstack([ad.read_h5ad(OUT/'snRNA'/s/'EC_corrected_and_raw_counts.h5ad').X for s in STAGES],format='csr')
    ordered=pd.concat([ann[ann.age.eq(s)] for s in STAGES])
    with h5py.File(OUT/'snRNA/P0/EC_corrected_and_raw_counts.h5ad') as h:var=read_elem(h['var'])
    a=ad.AnnData(lognorm(x),obs=ordered[['age','EC_family','sample']].copy(),var=var.copy())
    sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat')
    use=a.var.highly_variable.to_numpy() & ~var.mt.to_numpy() & ~var.ribo.to_numpy() & ~np.isin(var.index,WITHHOLD)
    a=a[:,use].copy();sc.pp.scale(a,zero_center=False,max_value=10);sc.pp.pca(a,n_comps=30,random_state=42)
    sc.pp.neighbors(a,n_neighbors=20,n_pcs=30,random_state=42);sc.tl.umap(a,random_state=42,min_dist=.4,maxiter=300)
    joint=ordered[['age','sample','EC_cluster_id','EC_gene_name','EC_family']].copy()
    joint['UMAP1'],joint['UMAP2']=a.obsm['X_umap'][:,0],a.obsm['X_umap'][:,1]
    joint.to_csv(OUT/'snRNA/all_stage_unintegrated_UMAP.tsv.gz',sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(17,7),layout='constrained')
    points(axes[0],a.obsm['X_umap'],joint.EC_family.to_numpy(),COLORS,s=5)
    points(axes[1],a.obsm['X_umap'],joint.age.to_numpy(),dict(zip(STAGES,['#0072b2','#e69f00','#009e73','#cc3311'])),s=5)
    axes[0].set_title('Subtype interpretation');axes[1].set_title('Postnatal stage')
    fig.suptitle('All-stage endothelial RNA · same UMAP, no batch/age integration\nStage and library effects can both contribute to separation',fontsize=15)
    finish(fig,'all_stage_snRNA_endothelial_UMAP',pdf)
    qc=pd.read_csv(OUT/'snRNA/all_endothelial_QC_decisions.tsv.gz',sep='\t')
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
    pc={'Primary endothelial reference':'#009e73','Insufficient canonical EC support; review':'#999999','Mixed-lineage program; excluded primary, retained review':'#d55e00','Coherent EC/mural program; excluded primary, retained review':'#8e5da3'}
    for name,c in pc.items():
        d=qc[qc.QC_disposition.eq(name)]
        axes[0].scatter(d.total_counts,d.n_genes_by_counts,c=c,s=3,label=f'{name.split(";")[0]} ({len(d):,})',rasterized=True)
        axes[1].scatter(d.decontx_contamination,d.Cardiomyocyte_fraction,c=c,s=3,rasterized=True)
    axes[0].set(xscale='log',xlabel='DecontX-corrected UMIs',ylabel='Detected genes')
    axes[0].legend(loc='upper left',fontsize=7)
    axes[1].axhline(.02,c='#444444',ls='--',lw=1);axes[1].set(xlabel='Estimated ambient fraction',ylabel='CM diagnostic fraction of full transcriptome')
    g=qc.groupby(['age','QC_disposition']).size().unstack(fill_value=0).reindex(STAGES)
    bottom=np.zeros(4)
    for c,color in pc.items():
        val=g[c].to_numpy() if c in g else np.zeros(4)
        axes[2].bar(STAGES,val,bottom=bottom,color=color);bottom+=val
    axes[2].set_ylabel('Previously labeled endothelial nuclei')
    fig.suptitle('RNA QC: inherited scDblFinder singlets and ambient correction, plus endothelial lineage review')
    finish(fig,'snRNA_endothelial_QC',pdf)
    # Module specificity is visible directly rather than inferred from UMAP position.
    panel=list(dict.fromkeys('Cdh5 Pecam1 Erg Dach1 Cd36 Rgcc Emcn Npr3 Nrxn3 Nfatc1 Gja5 Sox17 Sema3g Ephb4 Nr2f2 Prox1 Flt4 Ccl21a Top2a Mki67'.split()))
    values=[];detect=[];labels=[]
    for stage in STAGES:
        a=ad.read_h5ad(OUT/'snRNA'/stage/'EC_corrected_and_raw_counts.h5ad')
        z=lognorm(a.X);ix=a.var_names.get_indexer(panel)
        assert (ix>=0).all()
        for _,r in summary[summary.stage.eq(stage)].iterrows():
            m=a.obs.EC_cluster_id.eq(r.cluster_id).to_numpy();values.append(np.asarray(z[m][:,ix].mean(axis=0)).ravel())
            detect.append(np.asarray((a.X[m][:,ix]>=1).mean(axis=0)).ravel());labels.append(stage+' · '+r.gene_name)
    values=np.array(values);detect=np.array(detect)
    fig,ax=plt.subplots(figsize=(12,12),layout='constrained')
    y,xg=np.indices(values.shape);p=ax.scatter(xg.ravel(),y.ravel(),s=10+110*detect.ravel(),c=values.ravel(),cmap='viridis',edgecolors='none')
    ax.set_xticks(np.arange(len(panel)),panel,rotation=60,ha='right');ax.set_yticks(np.arange(len(labels)),labels);ax.invert_yaxis()
    ax.set_title('Endothelial marker evidence · color: mean log1p CP10k; dot area: detection fraction',pad=15)
    fig.colorbar(p,ax=ax,shrink=.45,label='Mean log1p CP10k')
    finish(fig,'RNA_endothelial_marker_evidence',pdf)

def spatialfigures(pdf):
    for stage in STAGES:
        folder=OUT/'Xenium'/stage
        a=pd.read_csv(folder/'EC_consensus_annotations.tsv.gz',sep='\t',index_col=0)
        atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        fig,axes=plt.subplots(1,3,figsize=(19,7),layout='constrained')
        for ax in axes:
            ax.scatter(atlas.x_centroid,atlas.y_centroid,s=.15,c='#e8e8e8',rasterized=True)
            ax.scatter(a.x_centroid,a.y_centroid,s=.5,c='#bbbbbb',alpha=.5,rasterized=True)
            ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='X (µm)',ylabel='Y (µm)')
        views=[(a.high_confidence_subtype,'fine_consensus_family','A/B: cell and nuclear concordance'),
               (a.assigned_including_moderate,'fine_consensus_family','A/B/C: adds boundary-stable cell calls'),
               (a.nuclear_supported_subtype,'nuclear_supported_family','Nuclear-supported subtypes · separate evidence')]
        for ax,(m,col,title) in zip(axes,views):
            for f,c in COLORS.items():
                if f==UNRESOLVED:continue
                d=a[m&a[col].eq(f)];ax.scatter(d.x_centroid,d.y_centroid,s=2.5,c=c,edgecolors='none',rasterized=True)
            ax.set_title(title+f'\nn={m.sum():,} / {len(a):,} EC candidates',fontsize=11)
        fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,label=f,markersize=5) for f,c in COLORS.items() if f!=UNRESOLVED],
            loc='outside lower center',ncol=3,frameon=False)
        fig.suptitle(f'{stage} endothelial spatial allocation\nGray EC candidates retain uncertain fine labels; nuclear-only support does not establish clean whole-cell RNA',fontsize=14)
        finish(fig,f'{stage}_Xenium_endothelial_confidence_maps',pdf)
        cl=pd.read_csv(folder/'EC_de_novo_clusters.tsv.gz',sep='\t',index_col=0)
        cs=pd.read_csv(folder/'EC_de_novo_cluster_summary.tsv',sep='\t')
        xy=cl[['Xenium_UMAP1','Xenium_UMAP2']].to_numpy()
        fig,axes=plt.subplots(1,2,figsize=(17,7),layout='constrained')
        cmap=dict(zip(cs.cluster_id,plt.get_cmap('tab20').colors))
        points(axes[0],xy,cl.Xenium_cluster_id.to_numpy(),cmap,s=3,legend=False)
        for _,r in cs.iterrows():
            m=cl.Xenium_cluster_id.eq(r.cluster_id);c=cl.loc[m,['Xenium_UMAP1','Xenium_UMAP2']].median().to_numpy()
            axes[0].annotate(r.gene_name.replace(' EC program',''),c,xytext=(7,4),textcoords='offset points',fontsize=8,
                bbox={'facecolor':'white','edgecolor':cmap[r.cluster_id],'alpha':.85,'pad':2},arrowprops={'arrowstyle':'-','color':'#777777','lw':.5})
        view=a.reindex(cl.index)
        points(axes[1],xy,view.fine_consensus_family.to_numpy(),COLORS,s=4)
        axes[0].set_title('Independent 403-panel RNA grouping');axes[1].set_title('Conservative transferred subtype · same UMAP')
        fig.suptitle(f'{stage} Xenium de novo endothelial subclustering · n={len(cl):,}\nGeometry and mixed-program QC; common RNA coverage across ownership views; no coordinates in clustering',fontsize=14)
        finish(fig,f'{stage}_Xenium_endothelial_UMAP',pdf)
    q=pd.read_csv(OUT/'Xenium/EC_RNA_ownership_depth_comparison.tsv',sep='\t')
    q=q[q.cohort.eq('native_broad_EC')]
    variants=['native','nuclear_anchored','high_confidence','depth_matched','nuclear_only','nuclear_depth_matched']
    names=['Native','New cell','Interior + nucleus','Core depth control','Nuclear','Nuclear depth control']
    fig,axes=plt.subplots(1,2,figsize=(16,6),layout='constrained')
    for s,c in zip(STAGES,['#0072b2','#e69f00','#009e73','#cc3311']):
        d=q[q.stage.eq(s)].set_index('variant').reindex(variants)
        axes[0].plot(names,100*d.mean_CM_fraction,marker='o',color=c,label=s)
        axes[1].plot(names,d.fine_supported,marker='o',color=c,label=s)
    for ax in axes:ax.tick_params(axis='x',rotation=35);ax.legend(frameon=False)
    axes[0].set_ylabel('Mean CM diagnostic RNA fraction (%)')
    axes[1].set_ylabel('Profiles passing each view’s strict subtype gates')
    axes[0].set_title('Measured program mixture · same native EC cohort')
    axes[1].set_title('Per-view support; consensus requires agreement')
    fig.suptitle('Ownership and molecule-depth sensitivity\nA reduced CM signal is diagnostic evidence; it does not measure the true contamination fraction',fontsize=14)
    finish(fig,'Xenium_endothelial_nuclear_and_boundary_sensitivity',pdf)

def boundaryfigures(pdf):
    from plot_xenium_boundary_qc import read_shapes
    rois=pd.read_csv(BOUND/'selected_ROIs.tsv',sep='\t');rois=rois[rois.category.eq('Endothelial')]
    # Older selection table uses non-CM interface categories; select exact fixed IDs.
    if len(rois)==0:
        rois=pd.read_csv(BOUND/'selected_ROIs.tsv',sep='\t');rois=rois[rois.roi_id.str.contains('_Endothelial_')]
    molecule_tables=[];audit=[]
    for stage in STAGES:
        atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        with np.load(SEG/stage/'whole_section_segmentation.npz') as b:
            grid=b['cell_labels'];res=float(b['grid_um']);cellids=b['cell_ids'].astype(str)
        ecids=set(atlas.index[atlas.atlas_identity.eq('Endothelial')])
        with h5py.File(SEG/stage/'nuclear_anchored_counts.h5ad') as h:genes=read_elem(h['var']).index.to_numpy()
        for roi in rois[rois.stage.eq(stage)].itertuples():
            folder=BOUND/stage;raw=pd.read_parquet(folder/(roi.roi_id+'_molecules.parquet'))
            native=read_shapes(folder/(roi.roi_id+'_native_polygons.geojson.gz'))
            frame=box(roi.x_min,roi.y_min,roi.x_max,roi.y_max)
            valid={c for c,g in native.items() if c in ecids and frame.covers(g)}
            m=raw.qv.ge(20)&raw.feature_name.isin(genes)&raw.cell_id.isin(valid)
            d=raw.loc[m].copy();d['stage']=stage;d['roi_id']=roi.roi_id
            d['region']=np.select([d.overlaps_nucleus.eq(1),~d.native_inside_polygon_Q20_target,
                d.native_boundary_distance_um_Q20_target.le(.85)],['Deposited nuclear','Outside native polygon','Cytoplasmic boundary <=0.85um'],default='Cytoplasmic interior >0.85um')
            d['CM_marker']=d.feature_name.isin(LINEAGES['Cardiomyocyte']);d['EC_marker']=d.feature_name.isin(LINEAGES['Endothelial'])
            molecule_tables.append(d[['stage','roi_id','transcript_id','cell_id','feature_name','region','CM_marker','EC_marker']])
            center=atlas.loc[roi.selection_cell_id,['x_centroid','y_centroid']].to_numpy()
            bounds=[center[0]-28,center[0]+28,center[1]-28,center[1]+28]
            with np.load(folder/(roi.roi_id+'_segmentation_maps.npz')) as b:
                mem=b['membrane'].astype(float);dna=b['DAPI'].astype(float);x0,y0,rr=float(b['x0_um']),float(b['y0_um']),float(b['grid_um'])
            extent=[x0,x0+mem.shape[1]*rr,y0+mem.shape[0]*rr,y0]
            lo,hi=np.quantile(mem,[.01,.995]);gray=np.sqrt(np.clip((mem-lo)/max(hi-lo,1),0,1))
            dlo,dhi=np.quantile(dna,[.02,.998]);dn=np.sqrt(np.clip((dna-dlo)/max(dhi-dlo,1),0,1))
            image=np.stack([1-.6*gray,1-.6*gray,1-.6*gray],axis=2);image[:,:,0]*=1-.3*dn;image[:,:,2]=np.maximum(image[:,:,2],.65*dn)
            fig,axes=plt.subplots(1,2,figsize=(14,7),layout='constrained')
            for ax in axes:
                ax.imshow(image,extent=extent,origin='upper',rasterized=True)
                ax.set(xlim=bounds[:2],ylim=[bounds[3],bounds[2]],xlabel='X (µm)',ylabel='Y (µm)');ax.set_aspect('equal')
            window=box(bounds[0],bounds[2],bounds[1],bounds[3])
            for c,g in native.items():
                if g.intersects(window):outline(axes[0],g,'#555555',.6)
            if roi.selection_cell_id in native:outline(axes[0],native[roi.selection_cell_id],'#007b88',1.6)
            xlo,xhi=int(bounds[0]/res),int(bounds[1]/res)+1;ylo,yhi=int(bounds[2]/res),int(bounds[3]/res)+1
            tile=grid[ylo:yhi,xlo:xhi];xx=np.arange(xlo,xlo+tile.shape[1])*res;yy=np.arange(ylo,ylo+tile.shape[0])*res
            for label in np.unique(tile):
                if label>0:
                    axes[1].contour(xx,yy,tile==label,levels=[.5],colors=['#007b88' if cellids[label-1]==roi.selection_cell_id else '#555555'],linewidths=[1.5 if cellids[label-1]==roi.selection_cell_id else .6])
            nuc=pd.read_parquet(folder/(roi.roi_id+'_nuclear_boundaries.parquet'))
            for _,part in nuc.groupby(['cell_id','label_id'],sort=False):
                xy=part[['vertex_x','vertex_y']].to_numpy()
                if len(xy)>2 and ((xy[:,0]>=bounds[0])&(xy[:,0]<=bounds[1])&(xy[:,1]>=bounds[2])&(xy[:,1]<=bounds[3])).any():
                    for ax in axes:ax.plot(xy[:,0],xy[:,1],c='#245bb0',lw=.8)
            p=raw[raw.qv.ge(20)&raw.x_location.between(bounds[0],bounds[1])&raw.y_location.between(bounds[2],bounds[3])]
            groups={'EC RNA':(LINEAGES['Endothelial'],'#008a95'),'CM RNA':(LINEAGES['Cardiomyocyte'],'#d97900'),'FB RNA':(LINEAGES['Fibroblast'],'#ac42ad')}
            for name,(panel,c) in groups.items():
                d=p[p.feature_name.isin(panel)]
                for ax in axes:
                    ax.scatter(d.x_location,d.y_location,c=c,s=8,alpha=.8,edgecolors='none',rasterized=True)
                    n=d[d.overlaps_nucleus.eq(1)]
                    ax.scatter(n.x_location,n.y_location,facecolors='none',edgecolors='#111111',s=20,lw=.5,rasterized=True)
            axes[0].set_title('Native cell boundaries + measured RNA');axes[1].set_title('Nuclear-anchored cell boundaries + same RNA')
            handles=[Line2D([],[],ls='',marker='o',color=c,label=n,markersize=5) for n,(_,c) in groups.items()]
            handles += [Line2D([],[],ls='',marker='o',markerfacecolor='none',markeredgecolor='black',label='Ring: deposited nuclear-overlap flag'),Line2D([],[],color='#245bb0',label='Deposited nuclear outline')]
            fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=9,frameon=False)
            fig.suptitle(f'{stage} · fixed {roi.roi_id.split("_")[-1]} EC–CM interface\n{roi.selection_cell_id} · boundaries are image-guided; RNA colors do not determine segmentation',fontsize=14)
            finish(fig,roi.roi_id+'_endothelial_nuclear_boundary_audit',pdf)
            audit.append({'stage':stage,'ROI':roi.roi_id,'fully_contained_native_EC_profiles':len(valid),'observed_Q20_gene_RNA':int(m.sum()),'source_selection_cell':roi.selection_cell_id})
        del grid
    allm=pd.concat(molecule_tables).drop_duplicates(['stage','transcript_id'])
    rows=[]
    for (stage,region),d in allm.groupby(['stage','region'],observed=True):
        rows.append({'stage':stage,'region':region,'molecules':len(d),'EC_profiles':d.cell_id.nunique(),
            'CM_diagnostic_RNA':int(d.CM_marker.sum()),'EC_diagnostic_RNA':int(d.EC_marker.sum()),
            'CM_diagnostic_fraction':d.CM_marker.mean(),'EC_diagnostic_fraction':d.EC_marker.mean()})
    save(pd.DataFrame(rows),OUT/'EC_interface_transcript_localization.tsv')
    save(pd.DataFrame(audit),OUT/'EC_interface_ROI_coverage.tsv')
    localization=pd.DataFrame(rows)
    regions=['Deposited nuclear','Cytoplasmic interior >0.85um','Cytoplasmic boundary <=0.85um']
    fig,axes=plt.subplots(1,2,figsize=(14,5.5),layout='constrained')
    for stage,c in zip(STAGES,['#0072b2','#e69f00','#009e73','#cc3311']):
        d=localization[localization.stage.eq(stage)].set_index('region').reindex(regions)
        axes[0].plot(['Nuclear','Cytoplasmic interior','Boundary band'],100*d.CM_diagnostic_fraction,marker='o',color=c,label=stage)
        axes[1].plot(['Nuclear','Cytoplasmic interior','Boundary band'],100*d.EC_diagnostic_fraction,marker='o',color=c,label=stage)
    for ax in axes:ax.set_ylabel('Diagnostic RNA / observed gene RNA (%)');ax.legend(frameon=False)
    axes[0].set_title('Cardiomyocyte-associated RNA');axes[1].set_title('Canonical endothelial RNA')
    fig.suptitle('Actual transcripts in the eight fixed EC–CM interfaces\nNative broad EC owners; fully contained polygons; selected ROIs; pooled molecules, not organ-wide contamination rates',fontsize=13)
    finish(fig,'EC_interface_nuclear_vs_cytoplasmic_RNA',pdf)

def catalog_and_readme():
    rna=pd.read_csv(OUT/'snRNA/EC_cluster_summary.tsv',sep='\t')
    qc=pd.read_csv(OUT/'snRNA/all_endothelial_QC_decisions.tsv.gz',sep='\t')
    q=pd.read_csv(OUT/'Xenium/EC_QC_summary.tsv',sep='\t')
    ownership=pd.read_csv(OUT/'Xenium/EC_RNA_ownership_depth_comparison.tsv',sep='\t')
    spatial=pd.read_csv(OUT/'Xenium/all_stage_EC_consensus.tsv.gz',sep='\t')
    with h5py.File(SEG/'P0/nuclear_anchored_counts.h5ad') as h:genes=read_elem(h['var']).index
    functions={FAMILIES[0]:'Microvascular exchange and transport; Cd36/Fabp4-associated lipid-handling programs are hypotheses from RNA, not measured transport.',
        FAMILIES[1]:'Endocardial/chamber-lining-associated signaling and interface programs; actual chamber-wall location requires spatial/histological validation.',
        FAMILIES[2]:'Arterial specification and endothelial coupling-associated program; arteriole anatomy is not proven from RNA alone.',
        FAMILIES[3]:'Candidate venous/venular specification; Nr2f2/Ephb4 can overlap other developing ECs and Vwf is not vein-specific.',
        FAMILIES[4]:'Candidate lymphatic specification and fluid-drainage biology; require endothelial identity plus Prox1/Flt4, not Lyve1 alone.'}
    keys=[]
    for family in FAMILIES:
        keys.append({'subtype_name':family,'RNA_panel':';'.join(PANELS[family]),'anchors':';'.join(ANCHORS[family]),
            'Xenium_measured_panel':';'.join(g for g in PANELS[family] if g in genes),
            'Xenium_measured_anchors':';'.join(g for g in ANCHORS[family] if g in genes),'potential_function_and_limit':functions[family]})
    save(pd.DataFrame(keys),OUT/'endothelial_subtype_name_key.tsv')
    # Keep RNA gene-program names and vascular interpretations as separate levels.
    rna['potential_function_and_limit']=rna.family.map(functions).fillna('Distinct RNA grouping with vascular subtype unresolved; avoid assigning a vessel compartment or novel biological type from separation alone.')
    rna['cycling_overlay']=np.where(rna.cycling_fraction.ge(.30),'Cycling-enriched transcriptional overlay','No dominant cycling overlay')
    save(rna,OUT/'RNA_named_subtype_catalog.tsv')
    family_rows=[]
    for stage in STAGES:
        d=spatial[spatial.stage.eq(stage)]
        for family in FAMILIES:
            family_rows.append({'stage':stage,'family':family,
                'RNA_nuclei':int(rna.loc[rna.stage.eq(stage)&rna.family.eq(family),'nuclei'].sum()),
                'Xenium_A_B_all_candidates':int((d.high_confidence_subtype&d.fine_consensus_family.eq(family)).sum()),
                'Xenium_A_B_C_all_candidates':int((d.assigned_including_moderate&d.fine_consensus_family.eq(family)).sum()),
                'Xenium_additional_nuclear_only_candidates':int((d.nuclear_only_subtype_candidate&d.nuclear_supported_family.eq(family)).sum())})
    save(pd.DataFrame(family_rows),OUT/'subtype_counts_by_stage_and_evidence.tsv')
    # Exact author marker list is reported separately from our curated classification.
    import zipfile,xml.etree.ElementTree as E
    z=zipfile.ZipFile(BASE/'metadata/paper_supplementary_tables.xlsx');ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    strings=[''.join(t.itertext()) for t in E.fromstring(z.read('xl/sharedStrings.xml'))]
    paper=[]
    for row in E.fromstring(z.read('xl/worksheets/sheet4.xml')).findall('s:sheetData/s:row',ns):
        v={}
        for c in row:
            t=c.find('s:v',ns)
            if t is not None:v[''.join(filter(str.isalpha,c.attrib['r']))]=strings[int(t.text)] if c.attrib.get('t')=='s' else t.text
        if 'EC' in v.get('A',''):
            paper.append({'author_state':v['A'],'gene':v['B'],'source_avg_log2FC':float(v['D']),'Xenium_measured':v['B'] in genes})
    save(pd.DataFrame(paper),OUT/'paper_endothelial_markers_and_coverage.tsv')
    qcstage=qc.groupby('age').agg(source=('age','size'),primary=('primary_EC_eligible','sum'),mural_review=('coherent_mural_group_review','sum')).reindex(STAGES)
    rows=['| Stage | Previously labeled EC nuclei | Primary RNA | Retained for QC review | RNA groups |',
          '|---|---:|---:|---:|---:|']
    for stage,d in qcstage.iterrows():rows.append(f'| {stage} | {d.source:,} | {d.primary:,} | {d.source-d.primary:,} | {len(rna[rna.stage.eq(stage)])} |')
    native=q[q.cohort.eq('native_broad_EC')].set_index('stage').reindex(STAGES)
    srows=['| Stage | Native broad EC profiles | Strict A/B fine labels | With moderate C | Extra nuclear-only candidates | Possible multi-nuclear merges |',
           '|---|---:|---:|---:|---:|---:|']
    for stage,d in native.iterrows():
        srows.append(f'| {stage} | {d.profiles:,} | {d.high_confidence_A_B:,} | {d.assigned_A_B_C:,} | {d.additional_nuclear_only_candidates:,} | {d.multiple_projected_nuclear_components_review:,} |')
    o=ownership[ownership.cohort.eq('native_broad_EC')]
    orows=['| Stage | Native whole cell | Refined whole cell | Boundary core | Core depth control | Nuclear | Nuclear depth control |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for stage in STAGES:
        d=o[o.stage.eq(stage)].set_index('variant')
        orows.append('| '+stage+' | '+' | '.join(f'{100*d.loc[v,"mean_CM_fraction"]:.1f}%' for v in ['native','nuclear_anchored','high_confidence','depth_matched','nuclear_only','nuclear_depth_matched'])+' |')
    text=f'''# Endothelial-specific QC, subclustering and Xenium nuclear/boundary reassessment

The completed reanalysis retains **{qc.primary_EC_eligible.sum():,} of 5,573** previously labeled endothelial nuclei for primary discovery. All 5,573 nuclei already passed whole-heart per-library doublet filtering and ambient correction. Additional endothelial-specific review retains the excluded nuclei and their original/corrected counts as an audit. The four Xenium sections contain **69,617 native broad EC profiles**; the new audit includes these plus uncertain EC proposals, totaling **95,316 candidates**. No original spatial profile or coordinate is deleted or moved.

[Full illustrated report, 27 pages](endothelial_QC_and_subtypes_report.pdf) · [Named RNA groups, markers and support](RNA_named_subtype_catalog.tsv) · [Vascular subtype definitions and possible functions](endothelial_subtype_name_key.tsv) · [Spatial QC summary](Xenium/EC_QC_summary.tsv) · [All spatial annotations](Xenium/all_stage_EC_consensus.tsv.gz) · [Validation](validation.json)

## RNA QC and discovery

{chr(10).join(rows)}

The inherited global QC uses scDblFinder **per entire library**, with Scrublet diagnostic only. DecontX fits converged separately in all eight libraries; called-nucleus clusters supply the contamination model because empty-droplet backgrounds are unavailable. The inherited filters require 750–7,500 detected features, at least 1,000 corrected UMIs, at most 5% mitochondrial RNA and estimated ambient fraction at most 50%. [Original QC summary](snRNA/inherited_QC_summary.tsv), [original validation/provenance](snRNA/inherited_QC_provenance.json) and unchanged upstream scripts under `analysis/endothelial_upstream_QC/` distinguish this completed upstream correction from the new analysis. We do not rerun droplet models on a purified lineage.

The new full-transcriptome screen requires at least two canonical EC genes with corrected count at least one and at least 0.1% EC-panel RNA. Strong incompatible programs (at least three genes, ten molecules and 2% of all RNA) enter review. A two-pass check also identifies **{int(qc.coherent_mural_group_review.sum()):,} nuclei** in provisional groups with at least three mural markers detected in at least 15% of nuclei in both raw and corrected counts in both libraries (at least ten nuclei/library). These groups retain coherent EC/mural programs; this is **not proof that every member is a doublet**. The final primary graph excludes them and is refitted. [Every QC decision](snRNA/all_endothelial_QC_decisions.tsv.gz), [provisional group evidence](snRNA/provisional_mural_group_review.tsv) and [threshold sensitivity](snRNA/QC_threshold_sensitivity.tsv) are supplied.

Discovery is independent at each age: 2,000 variable genes, 30 PCs, 20 neighbors, Leiden resolution0.45, seed42, with Harmony only between the two same-age libraries. Mitochondrial/ribosomal, sex and broad-lineage diagnostic genes do not drive the primary graph. Comparisons include unadjusted data, raw counts, nuisance covariates, explicit subtype-marker withholding, alternate seeds/resolutions, and a newly clustered ambient<=30% cohort. The [library mixing](snRNA/EC_library_mixing.tsv), [clustering sensitivity](snRNA/EC_clustering_sensitivity.tsv), [actual markers](snRNA/EC_cluster_markers.tsv) and [raw/corrected library-level family evidence](snRNA/family_markers_by_library.tsv) show the limits of finer grouping. Two RNA libraries per age are not automatically independent animals.

**Gene-enrichment names and vascular interpretations are separate levels.** Dach1/Cd36 and Ccdc85a/Ablim3 programs are examples of capillary-like groups; Tmem108/Npr3, Tmem108/Cdh11 and related programs support endocardial-like interpretation. Eln/Ptprj, Vegfc/Pcsk5 and Arl15/Prdm16 groups can support arterial-like interpretation when a broader arterial panel qualifies. Prox1/Adamtsl1 and Adamtsl1/Sema3d groups are small lymphatic-like candidates with Prox1/Flt4 support. A Vegfc/Hdac9 group has candidate venous-like panel support at P14. Cycling-associated groups are transcriptional overlays, not an additional vessel compartment. Gpc6-, Hdac9- and other unresolved programs are not asserted to be new vessel types. The catalog includes every group, its sample support, stability and proposed function. Some named groups are candidate/sensitivity-dependent; a stable graph alone does not prove anatomical identity.

The **all-stage RNA UMAP has no age or library integration**. Stage and technical effects can both contribute to its separation. The per-stage UMAPs use within-age library adjustment and independent axes; do not compare their geometric coordinates across stages. No layout is manually moved to exaggerate subtype separation.

![Same all-stage RNA coordinates colored by subtype and stage](all_stage_snRNA_endothelial_UMAP.png)

## Xenium RNA ownership and possible merged profiles

We reuse the completed **whole-section image-guided segmentation**, which already includes endothelial cells: deposited DAPI-derived nuclear outlines anchor watershed on the acquired ATP1A1/CD45/E-Cadherin boundary channel at0.425µm. Original profile IDs are retained; nuclear components of an existing profile share its parent ID. This is an alternative boundary refinement, **not independent new nuclear detection or demonstrated cell-number/merge-split correction**. Source nuclear transcript ownership is preserved; previously unassigned RNA is not forced into primary cells. Existing vendor segmentation also used interior stain information; the alternative is not assumed superior.

The endothelial-specific analysis compares native counts, refined whole-cell counts, nuclei plus cytoplasmic interior excluding a0.85µm boundary band, nuclear-only counts, cytoplasmic-only counts and exact per-profile molecule-depth controls. QV>=20 and the403 gene targets are retained. Nuclear and whole-cell fractions have a panel denominator, so RNA and Xenium percentage thresholds differ. Nuclear RNA is less exposed to a cell boundary but remains vulnerable to background, projected nuclear overlap and loss of genuine cytoplasmic state information. It is not independent ground truth.

Multiple deposited projected nuclear components can indicate a merged profile, a fragmented projection or true multinucleation. These profiles receive **merge review** and are excluded from the strict fine release and extra nuclear-only subtype candidate release. They are not deleted. Missing nuclear area or an absent outline alone does not remove a cell. Seed type, component count, native/refined area ratio and all original coordinates are exported. We make no claim that every homotypic merge has been detected or split.

## Fine assignment evidence

{chr(10).join(srows)}

The table uses the same **native broad EC cohort** in each stage. The all-candidate summary separately includes uncertain original EC proposals; added profiles remain explicit candidates. A new broad reference competition includes non-EC alternatives and excludes RNA ECs set aside by the new QC. Broad competition uses same-age RNA-library/subtype centroids with separate **whole-section** query calibration; fine transfer uses same-age eligible EC nuclei with separate **EC-candidate** query calibration. These calibrations intentionally use different populations.

Fine transfer requires a30-neighbor weighted vote>=0.6, margin>=0.15, agreement with raw reference and both separate RNA libraries, at least two actually measured subtype genes including an anchor, adequate RNA information (>=20 fitting molecules,>=8 genes), EC-lineage competition and no strong competing program. Lymphatic calls require **both Prox1 and Flt4**; Lyve1 alone is insufficient. RNA-reference unknown subtypes remain competitors. No spatial smoothing, subtype quota or published spatial label enters the fit. [Opposite-library technical benchmarks](Xenium/P21/opposite_RNA_library_benchmark.tsv) are supplied at all stages; reference labels are not independent biological truth.

- **A:** native, refined whole-cell and both nuclear calibrations support the same subtype.
- **B:** a newly supported whole-cell+nuclear call where the native fine call was unresolved.
- **C:** native, refined whole-cell and boundary-core subtype agreement, with nuclear confirmation unavailable.
- **Additional nuclear-only candidate:** both nuclear calibrations support the same subtype and no supported native/cell/core subtype conflicts, but the stringent whole-cell consensus is unmet. Many have mixed cell-body RNA. These are shown separately and never described as fully clean whole cells.

Strong mixed programs, incompatible supported subtype calls and unresolved/multiple-component geometry block the corresponding fine release. These tiers describe **operational robustness, not calibrated accuracy**. Shared molecules and references make the views correlated. The original parent atlas is not overwritten, and an unresolved fine label is not loss of EC lineage or evidence that a subtype is absent. All released assignments have explicit actual marker support and source views in their per-stage tables.

The independent Xenium subclusters use the measured panel after foreign-lineage/sex/outcome genes are withheld:20 PCs,20 neighbors,Leiden0.5,seed42, with alternate seeds/resolutions and core/nuclear/depth-control graphs. Final graphs exclude reviewed geometry and mixed cell-body programs and use **a common information-eligible cohort across all compared views**. Fine graph subdivisions are not automatically given vascular identities. The spatial maps show whole-section assignments; UMAPs show only the QC-eligible graph cohort.

## What the nuclear/boundary comparison establishes

Mean cardiomyocyte diagnostic RNA fraction among the **same native broad EC profiles**:

{chr(10).join(orows)}

Nuclear RNA has less CM-associated signal than molecule-depth-matched whole-cell controls at every age. Boundary-core counts also reduce that signal beyond the core depth control. However, the new whole-cell boundary refinement itself does **not** generally reduce CM-associated signal. Thus the useful improvement is a **nuclear-supported evidence layer and conservative boundary sensitivity**, not a claim that the alternative segmentation globally cleans every EC.

These fractions measure a diagnostic RNA program, **not the true ambient/contamination fraction**. The paper describes capillary ECs with myofiber-associated RNA. We do not erase that possible biology or establish it from mixed counts: strong mixed programs remain a review category, and their nuclear persistence, boundary sensitivity and depth controls are reported. Incomplete segmentation, ambient/background RNA and genuine expression remain alternatives.

The eight fixed, previously selected EC–CM interfaces are re-inspected using the **actual transcript coordinates, acquired stains, deposited nuclear outlines and both native/refined borders**. Rings denote deposited nuclear-overlap flags. [ROI coverage](EC_interface_ROI_coverage.tsv) and [transcript localization](EC_interface_transcript_localization.tsv) quantify ROI-observed QV>=20 molecules in fully contained native EC polygons, deduplicating transcript IDs across overlapping ROIs. These selected interfaces are descriptive; they are not a random, organ-wide estimate.

In these interfaces, nuclear RNA has lower CM-associated fractions and higher canonical EC fractions than cytoplasmic RNA. A CM-RNA spike specifically at the narrow boundary is **not consistent across stages**: cytoplasmic interiors can be similarly mixed. Thus simple edge erosion alone is insufficient, and projected nuclear/whole-cell RNA cannot establish whether all CM-associated signal is exogenous.

![Ownership and depth sensitivity](Xenium_endothelial_nuclear_and_boundary_sensitivity.png)

## Downloads and reproduction

Each `snRNA/P*/` directory contains primary and review H5ADs: **X is inherited DecontX-corrected fractional counts; `counts_raw` is original integer-valued RNA**. Together they reconstruct all5,573 original EC nuclei. The primary H5ADs include new annotations and UMAP coordinates. Duplicate full-transcriptome gene symbols are collapsed by sum; shared transfer targets have no duplicate symbols. Library-by-group raw and corrected pseudobulks are also supplied.

Each `Xenium/P*/EC_candidate_count_views.h5ad` contains **X=refined whole-cell counts**, with `native`, `high_confidence`, `nuclear_only`, `depth_matched` and `nuclear_depth_matched` layers, all403 targets, final annotations and native physical coordinates. Cytoplasmic counts equal X minus `nuclear_only`. These are alternative views of the same molecules; do not add them together or treat them as independent observations. All candidate rows, including review/unresolved rows, remain available. `EC_de_novo_embeddings_metadata_only.h5ad` contains the actual graph and independent Xenium UMAP. Per-view assignments, marker evidence, stage summaries, gene names, PDF/PNG figures and a SHA256 manifest accompany the counts.

Run the following in order with the existing input data and Python environment:

```bash
export MPLCONFIGDIR=/tmp/endothelial-mpl NUMBA_CACHE_DIR=/tmp/endothelial-numba
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
python analysis/reassess_endothelial_RNA.py
for stage in P0 P7 P14 P21; do python analysis/reassess_endothelial_Xenium.py --stage "$stage"; done
python analysis/consolidate_endothelial_QC.py
python analysis/report_endothelial_QC.py
python analysis/validate_endothelial_QC.py
```

`methods.json` and the per-stage completion files specify all gates and caveats. `runtime_versions.json`, original R session records and upstream QC provenance distinguish source and current environments. Validation checks exact downloaded count views against their sources, every retained/review RNA nucleus, actual spatial marker support, consensus logic, molecule-depth matching, original coordinates, nuclear-merge exclusions, all27 PDF pages and image integrity. No RNA count correction is silently repeated; prior CM maturation, Ankrd1 and fibroblast models are not refitted by this endothelial-specific task.

Sources: [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/), DOI10.1038/s44161-026-00861-z; [snRNA GSE290989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290989); [P0/P7 Xenium](https://zenodo.org/records/14961929), [P14/P21 Xenium](https://zenodo.org/records/14961942). The exact Supplementary Table3 EC markers and their coverage are provided separately from our curated identity panels. One section per stage and capture/selection biases preclude independent developmental abundance claims, causal cell signaling or proof of anatomical identity from RNA alone.
'''
    (OUT/'README.md').write_text(text)
    import importlib.metadata,sys
    packages=['numpy','pandas','scipy','scanpy','anndata','harmonypy','scikit-learn','h5py','matplotlib','shapely','pyarrow','Pillow']
    write_json(OUT/'runtime_versions.json',{'python':sys.version,'packages':{p:importlib.metadata.version(p) for p in packages}})

def main():
    q=pd.read_csv(OUT/'Xenium/EC_QC_summary.tsv',sep='\t');r=json.loads((OUT/'snRNA/RNA_complete.json').read_text())
    with PdfPages(OUT/'endothelial_QC_and_subtypes_report.pdf') as pdf:
        textpage(pdf,'Endothelial QC and subtypes · RNA + Xenium',[
            f"RNA: {r['source_EC_nuclei']:,} previously labeled endothelial nuclei; {r['primary_EC_nuclei']:,} in the primary reanalysis; {r['review_retained']:,} retained as weak-lineage or mixed-program review.",
            'All RNA input nuclei already passed per-library scDblFinder and converged DecontX ambient correction. This task verifies those calls and adds endothelial-specific residual-program QC; ambient correction is not repeated on an isolated lineage.',
            'Coherent EC/mural programs in both RNA libraries enter group review and are excluded from the final primary graph. Xenium parents with multiple projected nuclear components enter merge/multinucleation review; no successful cell splitting is asserted.',
            'RNA discovery: independent clustering at P0, P7, P14 and P21. Within-age library effects are adjusted; raw counts, unadjusted graphs, nuisance adjustment, marker withholding, seeds, resolutions and an ambient<=30% cohort are compared. Gene-based names identify enrichment programs, with internal IDs retained in tables.',
            'Xenium: all original broad ECs plus uncertain EC proposals are audited with the completed image-only whole-section nuclear-anchored segmentation. Native, refined cell, boundary core, nuclear-only, cytoplasmic-only and molecule-depth controls are compared. No droplet doublet or ambient algorithm is imposed on segmented profiles.',
            'Strict A/B labels require cell and nuclear agreement; C adds boundary-stable cell labels. Additional nuclear-supported subtype candidates are reported separately, including profiles whose cell bodies remain mixed. Unresolved profiles remain in the audit and gray map background.',
            'These are operational evidence levels, not calibrated biological accuracy. One Xenium section per age limits generalization. Original coordinates and parent atlas labels are preserved.'])
        lines=[]
        for _,d in q[q.cohort.eq('native_broad_EC')].iterrows():
            lines.append(f"{d.stage}: {d.profiles:,} native ECs; {d.high_confidence_A_B:,} strict A/B fine labels; {d.assigned_A_B_C:,} including C; {d.additional_nuclear_only_candidates:,} extra nuclear-only subtype candidates.")
        lines += ['Fine labels require RNA-reference agreement with both separate libraries and raw counts, multiple measured subtype markers including an anchor, broad-lineage competition, adequate RNA information and no conflicting strong evidence. Lymphatic candidates require both Prox1 and Flt4; Lyve1 alone is insufficient.',
            'A capillary-EC population with myofiber-associated RNA was described in the source paper. Our mixed-program flags do not prove those cells are doublets or that the RNA is biologically endogenous. Nuclear persistence and boundary/depth sensitivity are reported as evidence, not a binary biological verdict.',
            'Source: postnatal heart atlas, DOI10.1038/s44161-026-00861-z; snRNA GSE290989; Xenium Zenodo14961929 and14961942. Paper Supplementary Table3 markers supplement curated EC panels. Full methods, marker evidence, counts and validation accompany this report.']
        textpage(pdf,'Assignment coverage and interpretation limits',lines)
        rnafigures(pdf);spatialfigures(pdf);boundaryfigures(pdf)
    catalog_and_readme()
    print('REPORT_COMPLETE',OUT/'endothelial_QC_and_subtypes_report.pdf',flush=True)

if __name__=='__main__':main()
