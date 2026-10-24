"""Plot ventricular microenvironment hypotheses on unchanged stage UMAPs.

Consumes analysis outputs only. No scores, annotations, clusters or coordinates
are fitted or edited. Program colors are RNA-score hypotheses, not measured
oxygen, mechanical force or spatial position.
"""
from pathlib import Path
import argparse,json,os,re,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-microenvironment-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
ROOT=Path(__file__).resolve().parents[1]
STAGES=['P0','P7','P14','P21'];XY=['UMAP1','UMAP2']
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#6A3D9A','#17BECF','#B15928','#E31A1C','#526D35','#4D4D4D','#9467BD','#1B9E77','#A6761D','#E7298A']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'pdf.fonttype':42})
def sortkey(s):return tuple(int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)',str(s)))
def palette(ids):return {g:COLORS[i] if i<len(COLORS) else plt.get_cmap('hsv')(i/len(ids)) for i,g in enumerate(ids)}
def load(out):
 d=pd.read_csv(out/'cell_annotations_and_scores.csv.gz',index_col=0);s=pd.read_csv(out/'cluster_summary.tsv',sep='\t');cfg=json.loads((out/'program_definitions.json').read_text())
 need={'age','sample','gate_origin','microenvironment_cluster','microenvironment_label','support_status',*XY};assert need<=set(d),need-set(d)
 assert d.index.is_unique and d[list(need)].notna().all().all() and set(d.age)==set(STAGES) and np.isfinite(d[XY]).all().all()
 assert s.microenvironment_cluster.is_unique and s.nuclei.sum()==len(d)
 for cid,g in d.groupby('microenvironment_cluster'):
  r=s[s.microenvironment_cluster.eq(cid)].iloc[0];assert len(g)==r.nuclei and g.age.eq(r.stage).all()
  assert g.microenvironment_label.eq(r.microenvironment_label).all() and g.support_status.eq(r.support_status).all()
 assert set(d.microenvironment_cluster)==set(s.microenvironment_cluster)
 programs=cfg.get('display_programs',list(cfg['programs'])[:6]);assert programs
 for program in programs:
  assert program in cfg['programs'] and 'z_'+program in d
  col='coverage_'+program
  if col in d:assert d[col].dropna().between(-1e-9,1+1e-9).all()
 return d,s,cfg,programs

def base(ax,d,title):
 xy=d[XY].to_numpy();span=np.maximum(np.ptp(xy,axis=0),1e-9)
 ax.set(xlim=(xy[:,0].min()-.08*span[0],xy[:,0].max()+.08*span[0]),ylim=(xy[:,1].min()-.08*span[1],xy[:,1].max()+.08*span[1]),xlabel='Within-stage UMAP 1',ylabel='Within-stage UMAP 2',title=title)
 ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal',adjustable='box');ax.spines[['top','right']].set_visible(False)

def anchors(d):
 span=np.maximum(np.ptp(d[XY].to_numpy(),axis=0),1e-9);a={}
 for cid,g in d.groupby('microenvironment_cluster'):
  xy=g[XY].to_numpy();c=np.median(xy,axis=0);a[cid]=xy[np.argmin(np.sum(((xy-c)/span)**2,axis=1))]
 return a

def labels(ax,a,size=9):
 for cid,xy in a.items():ax.text(*xy,cid,ha='center',va='center',fontsize=size,weight='bold',zorder=8,bbox={'facecolor':'white','edgecolor':'#aaaaaa','linewidth':.35,'alpha':.9,'pad':1.5})

def categorical(ax,d,col,pal):
 z=d.sample(frac=1,random_state=42);colors=z[col].map(pal);assert colors.notna().all();ax.scatter(z.UMAP1,z.UMAP2,c=colors.tolist(),s=4,alpha=.82,linewidths=0,rasterized=True)

def save(fig,out,stem,pdf):
 fig.savefig(out/(stem+'.png'),dpi=180,facecolor='white');fig.savefig(out/(stem+'.pdf'),facecolor='white');pdf.savefig(fig);plt.close(fig)

def stage_clusters(d,s,stage,out,pdf):
 d=d[d.age.eq(stage)];ids=sorted(d.microenvironment_cluster.unique(),key=sortkey);pal=palette(ids);a=anchors(d);height=max(3,.58*(len(ids)+1))
 fig=plt.figure(figsize=(25,9+height));grid=fig.add_gridspec(2,3,height_ratios=[7,height],left=.035,right=.99,top=.85,bottom=.085,wspace=.14,hspace=.20)
 axes=[fig.add_subplot(grid[0,i]) for i in range(3)]
 categorical(axes[0],d,'microenvironment_cluster',pal);base(axes[0],d,'Data-driven subclusters');labels(axes[0],a)
 origin={'core_ventricular':'#bbbbbb','Slit2_overlap_reference':'#CC3311'};assert set(d.gate_origin)<=set(origin)
 # Draw reference cells last so their location remains visible.
 for originname in ['core_ventricular','Slit2_overlap_reference']:
  g=d[d.gate_origin.eq(originname)];axes[1].scatter(g.UMAP1,g.UMAP2,c=origin[originname],s=4.5 if originname.endswith('reference') else 2.7,alpha=.9 if originname.endswith('reference') else .5,linewidths=0,rasterized=True)
 base(axes[1],d,'Original ventricular gate + ambiguous reference');labels(axes[1],a,size=8)
 libraries=sorted(d['sample'].unique());assert len(libraries)==2;libs=dict(zip(libraries,['#0072B2','#D55E00']));categorical(axes[2],d,'sample',libs);base(axes[2],d,'Same UMAP colored by library');labels(axes[2],a,size=8)
 fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=f'{name.replace("_"," ")} (n={d.gate_origin.eq(name).sum():,})') for name,c in origin.items()],loc='upper center',bbox_to_anchor=(.51,.943),frameon=False,fontsize=11)
 fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=f'{lib} (n={d["sample"].eq(lib).sum():,})') for lib,c in libs.items()],loc='upper right',bbox_to_anchor=(.985,.943),frameon=False,fontsize=10)
 ax=fig.add_subplot(grid[1,:]);ax.axis('off');ss=s.set_index('microenvironment_cluster').loc[ids];rows=[]
 for cid,r in ss.iterrows():
  g=d[d.microenvironment_cluster.eq(cid)];nref=g.gate_origin.eq('Slit2_overlap_reference').sum()
  rows.append([cid,f'{int(r.nuclei):,}',textwrap.fill(str(r.microenvironment_label),52),textwrap.fill(str(r.support_status),50),f'{int(nref):,}'])
 table=ax.table(cellText=rows,colLabels=['Subcluster','Nuclei','Hypothesis-aligned label','Evidence / sensitivity status','Reference nuclei'],colWidths=[.11,.07,.35,.36,.11],cellLoc='left',colLoc='left',bbox=[0,0,1,1]);table.auto_set_font_size(False);table.set_fontsize(10.5)
 for (r,c),cell in table.get_celld().items():
  cell.set_edgecolor('#dddddd');cell.set_linewidth(.5)
  if r==0:cell.set_facecolor('#eeeeee');cell.get_text().set_weight('bold')
  elif c==0:cell.get_text().set_color(pal[ids[r-1]]);cell.get_text().set_weight('bold')
 fig.suptitle(f'{stage}: ventricular microenvironment-associated state candidates · n={len(d):,}',fontsize=22,weight='bold',y=.98)
 fig.text(.5,.035,'Same stage coordinates in all three panels · Ambiguous Slit2-overlap cells are a reference, not verified conduction or Purkinje cells\nLabels are RNA-based hypotheses, not measured wall location, oxygen exposure or mechanical force. Other ages have independent coordinates.',ha='center',fontsize=11)
 save(fig,out,stage+'_microenvironment_clusters',pdf)

def stage_programs(d,cfg,programs,stage,out,pdf,mincoverage):
 d=d[d.age.eq(stage)].sample(frac=1,random_state=42);cols=3;rows=int(np.ceil(len(programs)/cols));fig,axes=plt.subplots(rows,cols,figsize=(22,6.2*rows+1.7),squeeze=False);norm=Normalize(-2,2);cm=plt.get_cmap('RdBu_r')
 for ax,program in zip(axes.flat,programs):
  vals=d['z_'+program].to_numpy();valid=np.isfinite(vals);coverage='coverage_'+program
  if coverage in d:valid &=d[coverage].to_numpy()>=mincoverage
  ax.scatter(d.UMAP1,d.UMAP2,c='#cccccc',s=3,alpha=.35,linewidths=0,rasterized=True)
  points=ax.scatter(d.loc[valid,'UMAP1'],d.loc[valid,'UMAP2'],c=vals[valid],norm=norm,cmap=cm,s=4,alpha=.9,linewidths=0,rasterized=True)
  info=cfg['programs'][program];name=info.get('display_label',program.replace('_',' '));base(ax,d,textwrap.fill(name,40));ax.set_title(textwrap.fill(name,40)+f'\ncovered nuclei: {int(valid.sum()):,}/{len(d):,}',fontsize=13)
  cb=fig.colorbar(points,ax=ax,shrink=.72,pad=.02);cb.set_label('Within-library standardized RNA score',fontsize=9);cb.ax.tick_params(labelsize=9)
 for ax in list(axes.flat)[len(programs):]:ax.set_visible(False)
 fig.suptitle(f'{stage}: candidate spatial, oxygen-response and mechanical programs',fontsize=22,weight='bold',y=.985)
 fig.text(.5,.025,f'Same UMAP and axis limits throughout · Fixed color limits −2 to +2; scores are relative within libraries, not comparable oxygen measurements\nGray: nonfinite score or detected program-gene coverage <{mincoverage:.0%} where coverage is supplied. Programs can overlap and form gradients.\nRNA scores do not verify wall position, hypoxia, local force, or a discrete biological subtype.',ha='center',fontsize=11)
 fig.tight_layout(rect=(.01,.085,.99,.95));save(fig,out,stage+'_microenvironment_programs',pdf)

def overview(d,out,pdf):
 fig,axes=plt.subplots(2,2,figsize=(19,17));
 for ax,stage in zip(axes.flat,STAGES):
  s=d[d.age.eq(stage)];ids=sorted(s.microenvironment_cluster.unique(),key=sortkey);categorical(ax,s,'microenvironment_cluster',palette(ids));base(ax,s,f'{stage} · {len(s):,} nuclei · {len(ids)} groups');labels(ax,anchors(s),size=9)
 fig.suptitle('Ventricular microenvironment-associated candidates within each stage',fontsize=21,weight='bold',y=.985)
 fig.text(.5,.025,'Independent coordinate spaces between stages · Colors/IDs identify computational groups; see stage pages for labels and evidence status\nCandidate spatial or environmental interpretations require validation; the reference gate includes ambiguous Slit2-overlap cells.',ha='center',fontsize=11)
 fig.tight_layout(rect=(.015,.07,.985,.95));save(fig,out,'ventricular_microenvironment_overview',pdf)

def main():
 p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=ROOT/'results/ventricular_microenvironment');p.add_argument('--min-coverage',type=float,default=.2);a=p.parse_args();assert 0<=a.min_coverage<=1
 d,s,cfg,programs=load(a.outdir)
 with PdfPages(a.outdir/'ventricular_microenvironment_figures.pdf') as pdf:
  overview(d,a.outdir,pdf)
  for stage in STAGES:stage_clusters(d,s,stage,a.outdir,pdf);stage_programs(d,cfg,programs,stage,a.outdir,pdf,a.min_coverage)
 print(f'Saved stage clusters/origin/library figures, {len(programs)} program overlays per stage, overview and combined PDF; {len(d):,} nuclei passed alignment checks.')
if __name__=='__main__':main()
