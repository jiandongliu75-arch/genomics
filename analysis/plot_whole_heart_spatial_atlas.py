"""Consistent all-cell maps and measured spatial summaries of the Xenium atlas."""
import gzip
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from shapely.geometry import mapping
import numpy as np
import pandas as pd

from build_whole_heart_spatial_atlas import OUT,ROOT,SOURCE,STAGES,LINEAGES,REVIEW,DISPLAY_GENES
from define_heart_arterial_structures import read_boundaries
from plot_xenium_ventricular import background,frame

COLORS={'Cardiomyocyte':'#4878b7','Endothelial':'#00a087','Fibroblast':'#b65da3',
    'Pericyte':'#997344','Smooth muscle':'#e7a425','Immune':'#db4545','Epicardial':'#7850b4',
    'Mixed-marker review':'#6b6b6b','Ambiguous/review':'#b9b9b9','QC/information flag':'#e0e0e0'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})


def load(stage):
    return pd.read_csv(OUT/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')


def save(fig,name,pdf=None):
    fig.savefig(OUT/(name+'.png'),dpi=180,facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),facecolor='white')
    if pdf:pdf.savefig(fig,facecolor='white')
    plt.close(fig)


def footer(fig,text):fig.text(.5,.018,text,fontsize=9,ha='center',va='bottom',linespacing=1.4)


def categorical(ax,obs,size=.9):
    # Unresolved points below the resolved classes, without spatial smoothing.
    for label in REVIEW+LINEAGES:
        part=obs[obs.atlas_identity.eq(label)]
        ax.scatter(part.x_centroid,part.y_centroid,s=size,c=COLORS[label],edgecolors='none',rasterized=True)


def legend(ax,obs=None,fontsize=9,**kwargs):
    handles=[]
    for label in LINEAGES+REVIEW:
        n='' if obs is None else f' ({int(obs.atlas_identity.eq(label).sum()):,})'
        handles.append(Line2D([],[],marker='o',linestyle='',color=COLORS[label],markersize=5,label=label+n))
    return ax.legend(handles=handles,fontsize=fontsize,frameon=False,**kwargs)


def overview(pdf):
    fig,axes=plt.subplots(2,2,figsize=(14,14))
    for ax,stage in zip(axes.ravel(),STAGES):
        obs=load(stage);bg,extent,bounds=background(stage,obs)
        frame(ax,bg,extent,bounds,f'{stage}: {len(obs):,} original cells; {obs.atlas_resolved.mean():.1%} resolved')
        categorical(ax,obs,.6)
    handles=legend(axes[0,0],fontsize=9,loc='upper left',bbox_to_anchor=(0,-.1),ncol=5)
    handles.remove()
    fig.legend(handles=handles.legend_handles,labels=[x.get_text() for x in handles.texts],loc='lower center',bbox_to_anchor=(.5,.055),ncol=5,fontsize=10,frameon=False)
    fig.suptitle('All cardiac cells on the native P0–P21 Xenium sections',fontsize=20,weight='bold')
    footer(fig,'Same colors across stages. Direct broad-lineage assignments; uncertain and mixed profiles remain gray.\nEach panel is its own acquired section, shown in native µm coordinates; sections are not registered to one another.')
    fig.subplots_adjust(left=.06,right=.97,top=.93,bottom=.13,hspace=.23,wspace=.16)
    save(fig,'all_cardiac_cells_P0_P21',pdf)


def stage_maps(stage,pdf,tiles):
    obs=load(stage);bg,extent,bounds=background(stage,obs)
    fig,ax=plt.subplots(figsize=(12,10));frame(ax,bg,extent,bounds,f'{stage}: all cardiac cell identities')
    categorical(ax,obs,1);legend(ax,obs,fontsize=10,loc='upper left',bbox_to_anchor=(1.01,1))
    footer(fig,'Broad RNA-reference assignments with lineage-anchor and mixed-marker review.\nGray cells remain available in the downloadable atlas; colors do not imply fine subtype or anatomical validation.')
    fig.subplots_adjust(left=.08,right=.74,top=.92,bottom=.09)
    save(fig,f'{stage}_all_cell_atlas')
    fig,axes=plt.subplots(2,4,figsize=(21,11))
    for ax,label in zip(axes.ravel(),LINEAGES+['All review/QC flags']):
        frame(ax,bg,extent,bounds,label)
        ax.scatter(obs.x_centroid,obs.y_centroid,s=.15,c='#dddddd',edgecolors='none',rasterized=True)
        take=obs.atlas_identity.eq(label) if label in LINEAGES else ~obs.atlas_resolved
        part=obs[take];c=COLORS.get(label,'#777777')
        ax.scatter(part.x_centroid,part.y_centroid,s=1.3 if label!='Cardiomyocyte' else .6,c=c,edgecolors='none',rasterized=True)
        ax.text(.02,.98,f'{len(part):,} cells',transform=ax.transAxes,va='top',fontsize=10,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
    fig.suptitle(f'{stage}: spatial distributions of individual cardiac lineages',fontsize=20,weight='bold')
    footer(fig,'Each panel highlights one lineage on the same native section. Other original cells are light gray.\nNo assigned glial category is shown because the reference glial-marker panel is unassayed; this is not evidence of biological absence.')
    fig.subplots_adjust(left=.04,right=.98,top=.91,bottom=.10,hspace=.25,wspace=.25)
    save(fig,f'{stage}_lineage_spatial_panels',pdf)

    # Choose one region with a diverse resolved composition, without using anatomy/outcomes.
    tile=tiles[tiles.stage.eq(stage)].copy();p=tile[LINEAGES].to_numpy()/np.maximum(tile.resolved_cells.to_numpy()[:,None],1)
    entropy=-(p*np.log(np.maximum(p,1e-12))).sum(1)
    eligible=tile.resolved_cells.ge(40).to_numpy()&((p>0).sum(1)>=4)
    if not eligible.any():eligible=tile.resolved_cells.ge(20).to_numpy()
    score=np.where(eligible,entropy,-1);t=tile.iloc[np.argmax(score)]
    center=np.array([t.tile_x+.5,t.tile_y+.5])*100;b=[center[0]-250,center[0]+250,center[1]-250,center[1]+250]
    local=obs[obs.x_centroid.between(b[0],b[1])&obs.y_centroid.between(b[2],b[3])]
    polygons,provenance=read_boundaries(stage,set(local.index),obs)
    vertices=[];colors=[]
    for cell,g in polygons.items():
        parts=[g] if g.geom_type=='Polygon' else [part for part in g.geoms if part.geom_type=='Polygon']
        for part in parts:vertices.append(np.asarray(part.exterior.coords));colors.append(COLORS[local.at[cell,'atlas_identity']])
    fig,axes=plt.subplots(1,2,figsize=(15,8))
    frame(axes[0],bg,extent,b,f'{stage}: native cell polygons ({len(local):,} cells)')
    axes[0].add_collection(PolyCollection(vertices,facecolors=colors,edgecolors='white',linewidths=.15,alpha=.85,rasterized=True))
    frame(axes[1],bg,extent,b,'Same region: original centroids')
    categorical(axes[1],local,8);legend(axes[1],local,fontsize=8,loc='upper left',bbox_to_anchor=(1.01,1))
    footer(fig,'Example selected by high broad-lineage diversity in a 100 µm tile (≥40 resolved cells and ≥4 lineages when available).\nBoth panels show the same 500 × 500 µm region; no cell positions are moved. Cell segmentation is deposited, not newly inferred.')
    fig.subplots_adjust(left=.06,right=.82,top=.90,bottom=.13,wspace=.28)
    save(fig,f'{stage}_all_cell_geometry_zoom',pdf)
    local.to_csv(OUT/stage/'geometry_zoom_cells.tsv.gz',sep='\t',na_rep='NA')
    with gzip.open(OUT/stage/'geometry_zoom_polygons.geojson.gz','wt') as handle:
        json.dump({'type':'FeatureCollection','features':[{'type':'Feature','properties':{'cell_id':c,'atlas_identity':local.at[c,'atlas_identity']},'geometry':mapping(g)} for c,g in polygons.items()]},handle)
    (OUT/stage/'geometry_zoom_provenance.json').write_text(json.dumps({'stage':stage,'bounds_um':b,'selection_tile':[int(t.tile_x),int(t.tile_y)],'selection_entropy':float(entropy[np.argmax(score)]),'native_cells':len(local),**provenance},indent=2)+'\n')


def summaries(pdf):
    counts=pd.read_csv(OUT/'lineage_counts.tsv',sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(15,8))
    for ax,labels,title in zip(axes,[LINEAGES+REVIEW,LINEAGES[1:]],['All original cells: identities and review flags','Resolved non-CM composition']):
        base=np.zeros(4)
        den=counts[counts.atlas_identity.isin(labels)].groupby('stage').cells.sum().reindex(STAGES).to_numpy()
        for label in labels:
            n=counts[counts.atlas_identity.eq(label)].set_index('stage').cells.reindex(STAGES).to_numpy()
            height=n/den
            ax.bar(STAGES,height,bottom=base,color=COLORS[label],label=label);base+=height
        ax.set_ylim(0,1);ax.set_ylabel('Fraction of the indicated cells');ax.set_title(title,weight='bold')
    axes[1].legend(loc='upper left',bbox_to_anchor=(1.01,1),frameon=False,fontsize=9)
    handles=[Line2D([],[],marker='s',linestyle='',color=COLORS[k],label=k) for k in ['Cardiomyocyte']+REVIEW]
    axes[0].legend(handles=handles,loc='upper left',bbox_to_anchor=(0,-.08),frameon=False,fontsize=8,ncol=2)
    footer(fig,'Composition is conditional on the deposited sections, segmentation and conservative assignment criteria.\nOne section per age; proportions are descriptive and do not estimate replicated developmental changes in cell abundance.')
    fig.subplots_adjust(left=.06,right=.85,top=.90,bottom=.22,wspace=.28)
    save(fig,'cardiac_lineage_composition',pdf)

    profiles=pd.read_csv(OUT/'measured_gene_profiles.tsv',sep='\t')
    fig,axes=plt.subplots(2,2,figsize=(15,10))
    for ax,stage in zip(axes.ravel(),STAGES):
        part=profiles[profiles.stage.eq(stage)&profiles.atlas_identity.isin(LINEAGES)&profiles.gene.isin(DISPLAY_GENES)]
        for i,label in enumerate(LINEAGES):
            row=part[part.atlas_identity.eq(label)].set_index('gene').reindex(DISPLAY_GENES)
            ax.scatter(range(len(DISPLAY_GENES)),[i]*len(DISPLAY_GENES),s=160*row.detection_fraction,c=row.mean_log1pCP10k,cmap='viridis',vmin=0,vmax=6,edgecolors='#555555',linewidths=.3)
        ax.set_xticks(range(len(DISPLAY_GENES)),DISPLAY_GENES,rotation=60,ha='right',fontsize=9)
        ax.set_yticks(range(len(LINEAGES)),LINEAGES,fontsize=9);ax.invert_yaxis();ax.set_title(stage,weight='bold')
        ax.set_xlim(-.6,len(DISPLAY_GENES)-.4);ax.set_ylim(len(LINEAGES)-.5,-.5)
    m=plt.cm.ScalarMappable(cmap='viridis',norm=plt.Normalize(0,6))
    cax=fig.add_axes([.92,.30,.016,.40]);fig.colorbar(m,cax=cax,label='Mean measured log1p(CP10k)')
    fig.suptitle('Measured marker consistency of the audited broad lineages',fontsize=18,weight='bold')
    footer(fig,'Dot area: detected fraction (0–100%); color: mean measured expression, all403-target normalization.\nThese genes contribute to annotation and are not independent validation markers. Full gene summaries are downloadable.')
    fig.subplots_adjust(left=.12,right=.89,top=.91,bottom=.14,hspace=.48,wspace=.38)
    save(fig,'broad_lineage_measured_markers',pdf)

    d=pd.read_csv(OUT/'neighborhood_composition.tsv',sep='\t');fig,axes=plt.subplots(2,2,figsize=(14,13))
    for ax,stage in zip(axes.ravel(),STAGES):
        p=d[d.stage.eq(stage)].pivot(index='focal_lineage',columns='neighbor_lineage',values='observed_expected_ratio').reindex(index=LINEAGES,columns=LINEAGES)
        fraction=d[d.stage.eq(stage)].pivot(index='focal_lineage',columns='neighbor_lineage',values='neighbor_fraction').reindex(index=LINEAGES,columns=LINEAGES)
        img=ax.imshow(np.log2(p.clip(lower=.125)),vmin=-3,vmax=3,cmap='RdBu_r')
        ax.set_xticks(range(7),LINEAGES,rotation=55,ha='right',fontsize=8);ax.set_yticks(range(7),LINEAGES,fontsize=8)
        for i in range(7):
            for j in range(7):ax.text(j,i,f'{100*fraction.iloc[i,j]:.0f}%',ha='center',va='center',fontsize=8,color='white' if abs(np.log2(max(p.iloc[i,j],.125)))>1.8 else 'black')
        ax.set_title(stage,weight='bold');ax.set_ylabel('Focal cell lineage')
    cax=fig.add_axes([.92,.35,.016,.35]);fig.colorbar(img,cax=cax,label='log2(observed / global resolved fraction)')
    fig.suptitle('Broad-lineage spatial neighborhoods',fontsize=18,weight='bold')
    footer(fig,'Text: resolved-neighbor composition. Color: enrichment relative to each section’s overall resolved composition.\n20 nearest other original cells within100 µm; focal cells require≥5 resolved neighbors.\nDescriptive spatial co-occurrence; cell density, anatomy and ambiguity affect these values. No interaction or significance claim.')
    fig.subplots_adjust(left=.11,right=.88,top=.92,bottom=.15,hspace=.57,wspace=.28)
    save(fig,'cardiac_cell_neighborhoods',pdf)


def main():
    tiles=pd.read_csv(OUT/'spatial_100um_tile_counts.tsv.gz',sep='\t')
    with PdfPages(OUT/'whole_heart_spatial_atlas_report.pdf') as pdf:
        overview(pdf)
        for stage in STAGES:print('PLOT',stage,flush=True);stage_maps(stage,pdf,tiles)
        summaries(pdf)
    print('ATLAS_FIGURES_COMPLETE',flush=True)


if __name__=='__main__':main()
