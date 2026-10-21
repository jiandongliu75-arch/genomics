"""Compare stage and proposed subtype on identical, existing cardiomyocyte UMAPs."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/cell_identities_subtypes/Cardiomyocyte/cell_subcluster_annotations.csv.gz'
OUT=ROOT/'results/cardiomyocyte_focus';OUT.mkdir(exist_ok=True)
d=pd.read_csv(SOURCE,index_col=0)
mask=d.reviewed_identity.str.contains('unresolved|mixed',case=False,na=False)|d.proposed_subtype.str.contains('unresolved|mixed',case=False,na=False)
d=d.loc[~mask].copy();assert len(d)==46485
STAGES={'P0':'#0072B2','P7':'#E69F00','P14':'#009E73','P21':'#CC3311'}
SUBTYPES={'Ventricular-like':'#0072B2','Atrial-like':'#8E44AD','Conduction-associated CM candidate':'#D55E00'}
assert set(d.age)==set(STAGES) and set(d.proposed_subtype)==set(SUBTYPES)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'pdf.fonttype':42})
coords=['subcluster_UMAP1','subcluster_UMAP2'];xy=d[coords].to_numpy();assert np.isfinite(xy).all()
span=np.ptp(xy,axis=0);lo=xy.min(axis=0)-span*.08;hi=xy.max(axis=0)+span*.08
anchors={}
for group,s in d.groupby('subcluster_id'):
 v=s[coords].to_numpy();median=np.median(v,axis=0);anchors[group]=v[np.argmin(np.sum(((v-median)/span)**2,axis=1))]
def base(ax,title):
 ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),title=title,xlabel='UMAP 1',ylabel='UMAP 2');ax.set_xticks([]);ax.set_yticks([]);ax.spines[['top','right']].set_visible(False)
 ax.set_aspect('equal',adjustable='box')
 for name,a in anchors.items():
  ax.text(*a,name,ha='center',va='center',fontsize=10,weight='bold',bbox=dict(facecolor='white',edgecolor='#aaaaaa',linewidth=.4,alpha=.92,pad=2),zorder=5)
def points(ax,s,column,palette):
 s=s.sample(frac=1,random_state=42)
 ax.scatter(s[coords[0]],s[coords[1]],s=3.2,c=s[column].map(palette),alpha=.82,linewidths=0,rasterized=True)
def handles(palette,count=False,column=None):
 return [Line2D([],[],marker='o',linestyle='',markersize=9,color=c,label=k+(f' (n={int(d[column].eq(k).sum()):,})' if count else '')) for k,c in palette.items()]
notes='Same 46,485 cardiomyocytes, coordinates and axis limits throughout · CAR IDs are computational subclusters\nSubtype labels are provisional; conduction-associated candidates are not verified Purkinje cells. No unresolved/mixed assignments displayed.'
with PdfPages(OUT/'cardiomyocyte_stage_and_subtypes.pdf') as pdf:
 fig,axes=plt.subplots(1,2,figsize=(21,12),sharex=True,sharey=True)
 points(axes[0],d,'age',STAGES);base(axes[0],'A  |  Developmental stage')
 points(axes[1],d,'proposed_subtype',SUBTYPES);base(axes[1],'B  |  Proposed cardiomyocyte subtype')
 axes[0].legend(handles=handles(STAGES,True,'age'),loc='upper center',bbox_to_anchor=(.5,1.19),ncol=2,frameon=False,fontsize=12)
 axes[1].legend(handles=handles(SUBTYPES,True,'proposed_subtype'),loc='upper center',bbox_to_anchor=(.5,1.19),frameon=False,fontsize=12)
 fig.suptitle('Cardiomyocytes: stage and subtype on the same UMAP',fontsize=23,weight='bold',y=.98)
 fig.text(.5,.035,notes,ha='center',fontsize=12)
 fig.subplots_adjust(left=.045,right=.98,top=.76,bottom=.10,wspace=.12)
 fig.savefig(OUT/'cardiomyocyte_stage_vs_subtype.png',dpi=200);fig.savefig(OUT/'cardiomyocyte_stage_vs_subtype.pdf');pdf.savefig(fig);plt.close(fig)
 fig,axes=plt.subplots(2,2,figsize=(18,17),sharex=True,sharey=True)
 for ax,stage in zip(axes.flat,STAGES):
  s=d[d.age.eq(stage)]
  ax.scatter(d[coords[0]],d[coords[1]],s=1.5,c='#dddddd',alpha=.35,linewidths=0,rasterized=True)
  points(ax,s,'proposed_subtype',SUBTYPES);base(ax,f'{stage} · n={len(s):,}')
 fig.suptitle('Cardiomyocyte subtypes across developmental stages',fontsize=23,weight='bold',y=.985)
 fig.legend(handles=handles(SUBTYPES),loc='upper center',bbox_to_anchor=(.5,.96),ncol=3,frameon=False,fontsize=12)
 fig.text(.5,.025,notes+'\nColor: proposed subtype within the named stage; gray: cardiomyocytes from other stages. Representation is not a differential-abundance test.',ha='center',fontsize=11)
 fig.subplots_adjust(left=.055,right=.98,top=.90,bottom=.09,hspace=.16,wspace=.12)
 fig.savefig(OUT/'cardiomyocyte_subtypes_by_stage.png',dpi=200);fig.savefig(OUT/'cardiomyocyte_subtypes_by_stage.pdf');pdf.savefig(fig);plt.close(fig)
pd.crosstab(d.proposed_subtype,d.age).reindex(columns=list(STAGES)).to_csv(OUT/'subtype_by_stage_counts.tsv',sep='\t')
pd.crosstab([d.subcluster_id,d.proposed_subtype],d.age).reindex(columns=list(STAGES)).to_csv(OUT/'subcluster_subtype_by_stage_counts.tsv',sep='\t')
manifest={'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'displayed_nuclei':len(d),'hidden_unresolved_nuclei':int(mask.sum()),'coordinates':coords,'coordinate_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'method':'Existing coordinates, labels and memberships preserved; same axis limits and aspect ratio in all panels. Stage is color in overview A, proposed subtype in overview B and stage facets. No expression reanalysis or UMAP refitting.','stage_colors':STAGES,'subtype_colors':SUBTYPES,'files':[]}
for p in sorted(OUT.iterdir()):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('Saved cardiomyocyte-only stage/subtype comparison and stage facets; all 46,485 cardiomyocytes preserved.')
