"""Actual RNA locations and acquired stains under native/alternative borders."""
import argparse
import gzip
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from shapely.geometry import shape,box

from design_xenium_boundary_qc import OUT,ATLAS,FIT,STAGES
from reassess_whole_heart_lineages import PANELS

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
GCOLORS={'CM markers':'#dc7a00','EC markers':'#008fa5','Fibroblast markers':'#ad43b4','Ankrd1':'#d83948'}


def outline(ax,g,c='black',lw=.55,ls='-',alpha=1):
    if g.geom_type=='Polygon':
        xy=np.asarray(g.exterior.coords);ax.plot(xy[:,0],xy[:,1],c=c,lw=lw,ls=ls,alpha=alpha)
    elif hasattr(g,'geoms'):
        for part in g.geoms:outline(ax,part,c,lw,ls,alpha)


def read_shapes(path):
    with gzip.open(path,'rt') as h:features=json.load(h)['features']
    return {f['properties']['cell_id']:shape(f['geometry']) for f in features}


def render_roi(stage,roi,pdf=None):
    folder=OUT/stage
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    raw=pd.read_parquet(folder/(roi.roi_id+'_molecules.parquet'))
    native=read_shapes(folder/(roi.roi_id+'_native_polygons.geojson.gz'))
    alternative=read_shapes(folder/(roi.roi_id+'_membrane_watershed_polygons.geojson.gz'))
    with np.load(folder/(roi.roi_id+'_segmentation_maps.npz')) as b:
        mem=b['membrane'].astype(float);dna=b['DAPI'].astype(float)
        x0,y0,res=float(b['x0_um']),float(b['y0_um']),float(b['grid_um'])
    if roi.category=='Peri-arterial':
        fs=json.loads((FIT/stage/'arterial_segments.geojson').read_text())['features']
        f=next(f for f in fs if f['properties']['segment_id']==roi.segment_id)
        g=shape(f['geometry']);xlo,ylo,xhi,yhi=g.bounds
        bounds=[xlo-25,xhi+25,ylo-25,yhi+25]
        selected=f['properties']['arterial_EC_ids']+f['properties']['SMC_ids']
        info=f'{len(f["properties"]["arterial_EC_ids"])} arterial EC + {len(f["properties"]["SMC_ids"])} SMC seeds'
    else:
        center=obs.loc[roi.selection_cell_id,['x_centroid','y_centroid']].to_numpy()
        bounds=[center[0]-25,center[0]+25,center[1]-25,center[1]+25]
        selected=[roi.selection_cell_id,roi.partner_CM_id]
        o=obs.loc[roi.selection_cell_id]
        info=f'{roi.selection_cell_id}: {o.atlas_identity}; nuclear area '+('not measured' if pd.isna(o.nucleus_area) else f'{o.nucleus_area:.1f} µm²')
    extent=[x0,x0+mem.shape[1]*res,y0+mem.shape[0]*res,y0]
    lo,hi=np.quantile(mem,[.01,.995]);gray=np.sqrt(np.clip((mem-lo)/max(hi-lo,1),0,1))
    dlo,dhi=np.quantile(dna,[.02,.998]);d=np.sqrt(np.clip((dna-dlo)/max(dhi-dlo,1),0,1))
    image=np.stack([1-.55*gray,1-.55*gray,1-.55*gray],axis=2)
    image[:,:,2]=np.maximum(image[:,:,2],.65*d);image[:,:,0]*=1-.3*d
    window=box(bounds[0],bounds[2],bounds[1],bounds[3])
    shown={cell:g for cell,g in native.items() if g.intersects(window)}
    visible=raw.x_location.between(bounds[0],bounds[1])&raw.y_location.between(bounds[2],bounds[3])
    points=raw[visible&raw.qv.ge(20)].copy()
    fig,axes=plt.subplots(1,3,figsize=(18,6.5))
    for ax,title in zip(axes,['Acquired stains + native cell borders','Measured QV≥20 RNA + native borders','Same RNA + membrane-guided alternative']):
        ax.imshow(image,extent=extent,origin='upper',rasterized=True)
        ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[3],bounds[2]);ax.set_aspect('equal')
        ax.set_title(title,fontsize=12,weight='bold');ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
    for cell,g in shown.items():
        for ax in axes[:2]:outline(ax,g,'#333333',.6,alpha=.8)
        outline(axes[2],alternative.get(cell,g),'#278345',.65,alpha=.9)
    # Nuclear outlines are deposited projections; no nucleus criterion defines a seed.
    nuclei=pd.read_parquet(folder/(roi.roi_id+'_nuclear_boundaries.parquet'))
    for _,part in nuclei.groupby(['cell_id','label_id'],sort=False):
        xy=part[['vertex_x','vertex_y']].to_numpy()
        if len(xy)>2 and ((xy[:,0]>=bounds[0])&(xy[:,0]<=bounds[1])&(xy[:,1]>=bounds[2])&(xy[:,1]<=bounds[3])).any():
            axes[0].plot(xy[:,0],xy[:,1],c='#386ca2',lw=.7,alpha=.8)
    for cell in selected:
        if cell in native:
            outline(axes[0],native[cell],'#111111',1.3)
            if cell in alternative:outline(axes[2],alternative[cell],'#126b35',1.3)
    if roi.category!='Peri-arterial':
        o=obs.loc[roi.selection_cell_id]
        for ax in axes:ax.plot(o.x_centroid,o.y_centroid,marker='+',c='black',ms=9,mew=1.3)
    groups={'CM markers':PANELS['Cardiomyocyte'],'EC markers':PANELS['Endothelial'],'Fibroblast markers':PANELS['Fibroblast']}
    if roi.category=='Peri-arterial':groups['Ankrd1']=['Ankrd1']
    rng=np.random.default_rng(944);sampled=False
    for label,genes in groups.items():
        p=points[points.feature_name.isin(genes)]
        if len(p)>5000:p=p.iloc[rng.choice(len(p),5000,replace=False)];sampled=True
        for ax in axes[1:]:ax.scatter(p.x_location,p.y_location,s=3 if roi.category!='Peri-arterial' else 2,c=GCOLORS[label],alpha=.7,edgecolors='none',rasterized=True)
    handles=[Line2D([],[],marker='o',linestyle='',color=GCOLORS[k],markersize=4,label=k) for k in groups]
    axes[1].legend(handles=handles,loc='upper left',fontsize=7,framealpha=.9)
    axes[0].text(.02,.98,'Gray: boundary stain\nBlue: DAPI / nuclear outlines',transform=axes[0].transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    length=20 if bounds[1]-bounds[0]<100 else 50
    for ax in axes:
        x=bounds[0]+.06*(bounds[1]-bounds[0]);y=bounds[3]-.08*(bounds[3]-bounds[2])
        ax.plot([x,x+length],[y,y],c='black',lw=2);ax.text(x,y-.015*(bounds[3]-bounds[2]),f'{length} µm',fontsize=8)
    fig.suptitle(f'{roi.roi_id}: {roi.category}',fontsize=17,weight='bold')
    note='RNA points sampled to≤5,000 per marker group for display; all molecules retained in data.' if sampled else 'Every displayed marker molecule retains its original position; no expression smoothing.'
    fig.text(.5,.035,info+'\n'+note+'\nAlternative: fixed original cell seeds, including anucleate profiles; boundary stain only. This is a sensitivity test, not validated replacement segmentation.',ha='center',fontsize=9,linespacing=1.4)
    fig.subplots_adjust(left=.04,right=.98,top=.86,bottom=.17,wspace=.25)
    fig.savefig(folder/(roi.roi_id+'_transcript_boundary_inspection.png'),dpi=180,facecolor='white')
    if pdf:pdf.savefig(fig,facecolor='white')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--stages',nargs='+',choices=STAGES,default=STAGES);args=p.parse_args()
    rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t')
    for stage in args.stages:
        with PdfPages(OUT/(stage+'_transcript_boundary_gallery.pdf')) as pdf:
            for roi in rois[rois.stage.eq(stage)].itertuples():
                render_roi(stage,roi,pdf);print('RNA_INSPECTION_FIGURE',roi.roi_id,flush=True)
    print('TRANSCRIPT_BOUNDARY_GALLERIES_COMPLETE',flush=True)


if __name__=='__main__':main()
