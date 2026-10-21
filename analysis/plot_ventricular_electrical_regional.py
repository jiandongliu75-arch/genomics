"""Display replicate-aware electrical/regional gene contrasts without subtype claims.

Inputs are computed on original expression. Dots indicate concordant directions
across eligible libraries and raw-count checks, not statistical significance.
"""
from pathlib import Path
import argparse,os,re
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-electrical-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]
STAGES=['P0','P7','P14','P21']
GENES=['Slit2','Mpped2','Gja5','Etv1','Irx3','Cntn2','Scn10a','Pcp4','Corin','Mylk4','Tbx3','Irx5','Kcnd2','Kcnip2','Scn5a','Gja1','Hey2']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})

def sortkey(s):return tuple(int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)',str(s)))
def asbool(v):
 if pd.isna(v):return False
 if isinstance(v,(bool,np.bool_)):return bool(v)
 return str(v).lower() in {'true','1','1.0','yes'}

def load(path):
 t=pd.read_csv(path,sep='\t')
 need={'age','cluster','gene','nuclei','reference_nuclei','eligible_libraries','corrected_delta','raw_delta','both_libraries_corrected_positive','both_libraries_corrected_negative','both_libraries_raw_same_direction'}
 assert need<=set(t),need-set(t)
 assert not t.duplicated(['age','cluster','gene']).any()
 assert set(t.age)==set(STAGES) and set(GENES)<=set(t.gene)
 for (stage,cluster),g in t.groupby(['age','cluster']):
  assert g.nuclei.nunique()==1 and g.reference_nuclei.nunique()==1
  assert 0<=g.reference_nuclei.iloc[0]<=g.nuclei.iloc[0]
  assert set(GENES)<=set(g.gene),(stage,cluster)
 for c in ['both_libraries_corrected_positive','both_libraries_corrected_negative','both_libraries_raw_same_direction']:t[c]=t[c].map(asbool)
 return t

def stage_data(t,stage):
 z=t[t.age.eq(stage)].copy();ids=sorted(z.cluster.unique(),key=sortkey)
 effects=z.pivot(index='cluster',columns='gene',values='corrected_delta').reindex(index=ids,columns=GENES)
 eligible=z.pivot(index='cluster',columns='gene',values='eligible_libraries').reindex(index=ids,columns=GENES)
 valid=eligible.ge(2)&np.isfinite(effects)
 consistent=((z.corrected_delta.gt(0)&z.both_libraries_corrected_positive)|(z.corrected_delta.lt(0)&z.both_libraries_corrected_negative))&z.both_libraries_raw_same_direction
 z['direction_dot']=consistent
 dots=z.pivot(index='cluster',columns='gene',values='direction_dot').reindex(index=ids,columns=GENES).fillna(False)&valid
 info=z.drop_duplicates('cluster').set_index('cluster').loc[ids]
 labels=[f'{cid}  n={int(r.nuclei):,}  ref={r.reference_nuclei/r.nuclei:.0%}' for cid,r in info.iterrows()]
 return effects.where(valid),dots,info,labels

def block(ax,t,stage,norm,cmap,showx=True):
 effects,dots,info,labels=stage_data(t,stage)
 im=ax.imshow(np.ma.masked_invalid(effects.to_numpy(dtype=float)),cmap=cmap,norm=norm,aspect='auto',interpolation='none')
 ax.set_yticks(np.arange(len(effects)),labels,fontsize=10)
 for tick,(_,r) in zip(ax.get_yticklabels(),info.iterrows()):
  if r.reference_nuclei/r.nuclei>=.5:tick.set_color('#8E236B');tick.set_weight('bold')
 ax.set_xticks(np.arange(len(GENES)),GENES if showx else [],rotation=55,ha='right',fontsize=10)
 ax.set_title(f'{stage} · {len(effects)} computational groups',loc='left',weight='bold',fontsize=14,pad=8)
 ax.set_xticks(np.arange(-.5,len(GENES),1),minor=True);ax.set_yticks(np.arange(-.5,len(effects),1),minor=True)
 ax.grid(which='minor',color='white',linewidth=.65);ax.tick_params(which='minor',bottom=False,left=False);ax.tick_params(axis='both',length=0)
 y,x=np.where(dots.to_numpy());ax.scatter(x,y,s=13,c='black',marker='o',linewidths=0)
 for spine in ax.spines.values():spine.set_visible(False)
 return im

def notes(limit):return (f'Color: equal-library mean cluster-minus-other-CMs difference in mean log1p(counts / total × 10,000), clipped to ±{limit:g}; not log2 fold change.\n'
 'Gray: fewer than two eligible libraries or missing contrast. A library requires ≥30 cluster nuclei and ≥30 other CMs from the same stage.\n'
 'Black dot: same corrected-count effect direction in both libraries, with raw-count direction concordance; not a significance test.\n'
 'ref = ambiguous Slit2-overlap reference fraction; purple labels mark ≥50% reference cells. Slit2/Mpped2/Gja5/Etv1 partly overlap selection panels, not independent validation.\n'
 'Gene patterns do not verify wall location, mechanical load, oxygen exposure or Purkinje identity. Comparisons are within each stage.')

def main():
 p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=ROOT/'results/ventricular_electrical_regional');p.add_argument('--limit',type=float,default=1.5);a=p.parse_args();assert a.limit>0
 t=load(a.outdir/'cluster_gene_contrasts.tsv');norm=Normalize(-a.limit,a.limit);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d5d5d5')
 nrows=[t[t.age.eq(stage)].cluster.nunique() for stage in STAGES]
 fig=plt.figure(figsize=(18,21));grid=fig.add_gridspec(4,1,height_ratios=nrows,left=.255,right=.915,top=.93,bottom=.185,hspace=.39)
 for i,stage in enumerate(STAGES):im=block(fig.add_subplot(grid[i,0]),t,stage,norm,cmap)
 cb=fig.colorbar(im,cax=fig.add_axes([.935,.32,.018,.48]),extend='both');cb.set_label('Mean log-normalized expression difference',fontsize=12)
 fig.suptitle('Electrical and regional-associated gene contrasts within each stage',fontsize=21,weight='bold',y=.978)
 fig.text(.5,.025,notes(a.limit),ha='center',fontsize=10.5,linespacing=1.5)
 fig.savefig(a.outdir/'ventricular_electrical_regional_contrasts.png',dpi=180,facecolor='white');fig.savefig(a.outdir/'ventricular_electrical_regional_contrasts.pdf',facecolor='white');plt.close(fig)
 for stage in STAGES:
  nr=t[t.age.eq(stage)].cluster.nunique();fig,ax=plt.subplots(figsize=(18,max(7.5,3+.55*nr)))
  im=block(ax,t,stage,norm,cmap);fig.subplots_adjust(left=.255,right=.90,top=.83,bottom=.31)
  cb=fig.colorbar(im,cax=fig.add_axes([.925,.40,.018,.37]),extend='both');cb.set_label('Mean log-normalized expression difference',fontsize=11)
  fig.suptitle(f'{stage}: electrical and regional-associated gene contrasts',fontsize=20,weight='bold',y=.96)
  fig.text(.5,.02,notes(a.limit),ha='center',fontsize=9.5,linespacing=1.45)
  fig.savefig(a.outdir/(stage+'_electrical_regional_contrasts.png'),dpi=180,facecolor='white');fig.savefig(a.outdir/(stage+'_electrical_regional_contrasts.pdf'),facecolor='white');plt.close(fig)
 print(f'Saved four stage heatmaps and combined PNG/PDF; {len(t)} gene/group contrasts passed structural checks. Dots are direction consistency only.')
if __name__=='__main__':main()
