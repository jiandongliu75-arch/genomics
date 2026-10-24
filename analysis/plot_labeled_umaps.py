"""Redraw existing UMAP coordinates with readable external identity callouts."""
from pathlib import Path
import hashlib,json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'results/cell_identities_subtypes'; OUT=ROOT/'results/labeled_umaps';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'pdf.fonttype':42})
PALETTE=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#6A3D9A','#17BECF','#B15928','#E31A1C','#526D35','#4D4D4D','#9467BD','#1B9E77','#A6761D','#E7298A']
STAGES={'P0':'#0072B2','P7':'#E69F00','P14':'#009E73','P21':'#CC3311'}
records=[]
exclusions=[]
counts=[]
def filter_display(d,view):
 mask=d.reviewed_identity.str.contains('unresolved|mixed',case=False,na=False) | d.proposed_subtype.str.contains('unresolved|mixed',case=False,na=False)
 for idx,row in d.loc[mask].iterrows():
  exclusions.append({'view':view,'cell_id':idx,'reviewed_identity':row.reviewed_identity,'subcluster_id':row.subcluster_id,'proposed_subtype':row.proposed_subtype,'age':row.age})
 shown=d.loc[~mask].copy()
 assert not shown.proposed_subtype.str.contains('unresolved|mixed',case=False).any()
 assert set(shown.age)<=set(STAGES)
 for (identity,age),n in shown.groupby(['reviewed_identity','age']).size().items():counts.append({'view':view,'identity':identity,'stage':age,'displayed_nuclei':int(n)})
 return shown

def draw(d,title,group,coords,stem,subtype_colors=False,stage_colors=False):
 groups=sorted(d[group].unique(),key=lambda x:int(x.split('_')[-1]) if group=='subcluster_id' else x)
 cats=sorted(d.proposed_subtype.unique()) if subtype_colors else groups
 colors={c:PALETTE[i%len(PALETTE)] for i,c in enumerate(cats)}
 fig=plt.figure(figsize=(19,12));ax=fig.add_axes([.29,.14,.42,.72])
 xy=d[list(coords)].to_numpy(); assert np.isfinite(xy).all()
 span=np.ptp(xy,axis=0);lo=xy.min(axis=0)-span*.07;hi=xy.max(axis=0)+span*.07
 centers=[]
 if stage_colors:
  shuffled=d.sample(frac=1,random_state=42)
  ax.scatter(shuffled[coords[0]],shuffled[coords[1]],s=2.4 if len(d)>20000 else 6,c=shuffled.age.map(STAGES),alpha=.8,linewidths=0,rasterized=True)
  from matplotlib.lines import Line2D
  ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=color,label=f'{age} (n={int(d.age.eq(age).sum()):,})',markersize=8) for age,color in STAGES.items()],title='Postnatal stage',loc='upper left',fontsize=11,framealpha=.98)
 for g in groups:
  s=d[d[group]==g];v=s[list(coords)].to_numpy();center=np.median(v,axis=0)
  # Use an observed point near the median, so the connector touches the group.
  anchor=v[np.argmin(np.sum(((v-center)/span)**2,axis=1))]
  category=s.proposed_subtype.iloc[0] if subtype_colors else g;c='#333333' if stage_colors else colors[category]
  if not stage_colors:ax.scatter(v[:,0],v[:,1],s=2.4 if len(d)>20000 else 6,c=c,alpha=.8,linewidths=0,rasterized=True)
  if group=='subcluster_id':
   subtype=s.proposed_subtype.iloc[0];review=bool(s.subtype_review_flag.any())
   label=f'{g} · n={len(s):,}\n'+textwrap.fill(subtype,29)
   if review and 'candidate' not in subtype.lower() and 'unresolved' not in subtype.lower():label+=' [review]'
  else:
   label=textwrap.fill(g + (' candidate' if g=='Glial' else ''),27)+f'\nn={len(s):,}'
  centers.append((g,anchor,c,label))
 ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),xlabel='UMAP 1',ylabel='UMAP 2')
 ax.set_xticks([]);ax.set_yticks([]);ax.spines[['top','right']].set_visible(False)
 split=np.median([v[1][0] for v in centers])
 for left in [True,False]:
  side=sorted([v for v in centers if (v[1][0]<=split)==left],key=lambda v:v[1][1])
  for entry,y in zip(side,np.linspace(.06,.94,len(side))):
   g,anchor,c,label=entry
   ax.annotate(label,xy=anchor,xycoords='data',xytext=(-.13 if left else 1.13,y),textcoords='axes fraction',ha='right' if left else 'left',va='center',fontsize=11.5,color=c,weight='bold',annotation_clip=False,bbox=dict(boxstyle='round,pad=.35',facecolor='white',edgecolor=c,alpha=.97),arrowprops=dict(arrowstyle='-',color=c,lw=1.1,alpha=.8,connectionstyle='arc3,rad=0'))
   if group=='subcluster_id':
    ax.text(*anchor,g,fontsize=8,ha='center',va='center',weight='bold',bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1))
 fig.suptitle(title,fontsize=23,weight='bold',y=.96)
 fig.text(.5,.91,'Color: postnatal stage (P0, P7, P14, P21)' if stage_colors else 'Color: proposed subtype' if subtype_colors else 'Color: computational subcluster' if group=='subcluster_id' else 'Color: reviewed broad cell identity',ha='center',fontsize=14)
 fig.text(.5,.055,'Existing UMAP coordinates preserved · Labels are provisional RNA-based identities\nUnresolved/mixed assignments hidden from figures only · [review] = additional uncertainty · Candidate labels retained',ha='center',fontsize=11,color='#444444')
 fig.savefig(OUT/(stem+'.png'),dpi=200,facecolor='white');fig.savefig(OUT/(stem+'.pdf'),facecolor='white');pdf.savefig(fig);plt.close(fig)
 records.append({'plot':stem,'nuclei':len(d),'groups':len(groups),'coordinates':list(coords),'source':str((SRC/('all_cell_identity_and_subtype_annotations.csv.gz' if group!='subcluster_id' else title.split(' — ')[0]+'/cell_subcluster_annotations.csv.gz')).relative_to(ROOT))})
def stage_facets(d,title,group,coords,stem):
 fig,axes=plt.subplots(2,2,figsize=(18,14),sharex=True,sharey=True)
 xy=d[list(coords)].to_numpy();span=np.ptp(xy,axis=0);lo=xy.min(axis=0)-span*.08;hi=xy.max(axis=0)+span*.08
 anchors=[]
 for g,s in d.groupby(group):
  v=s[list(coords)].to_numpy();center=np.median(v,axis=0);anchor=v[np.argmin(np.sum(((v-center)/span)**2,axis=1))];anchors.append((g,anchor))
 for ax,(stage,color) in zip(axes.flat,STAGES.items()):
  s=d[d.age==stage]
  ax.scatter(d[coords[0]],d[coords[1]],s=1.5,c='#dddddd',alpha=.35,linewidths=0,rasterized=True)
  ax.scatter(s[coords[0]],s[coords[1]],s=3 if len(d)>20000 else 8,c=color,alpha=.85,linewidths=0,rasterized=True)
  for g,a in anchors:
   ax.annotate(g,xy=a,xytext=(3,3),textcoords='offset points',fontsize=8,weight='bold',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1))
  ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),title=f'{stage} · n={len(s):,}',xlabel='UMAP 1',ylabel='UMAP 2');ax.set_xticks([]);ax.set_yticks([]);ax.spines[['top','right']].set_visible(False)
 fig.suptitle(title+' — stage comparison',fontsize=22,weight='bold',y=.98)
 fig.text(.5,.025,'Identical coordinates and axes across stages · Gray: other stages · Unresolved/mixed assignments hidden from figures only\nCluster IDs refer to the accompanying labeled subtype views; stage plots show representation, not differential abundance tests.',ha='center',fontsize=11)
 fig.tight_layout(rect=(.02,.065,.98,.95));fig.savefig(OUT/(stem+'.png'),dpi=200);fig.savefig(OUT/(stem+'.pdf'));pdf.savefig(fig);plt.close(fig)
 records.append({'plot':stem,'nuclei':len(d),'groups':int(d[group].nunique()),'coordinates':list(coords),'stage_counts':{a:int(d.age.eq(a).sum()) for a in STAGES}})
with PdfPages(OUT/'all_labeled_umaps.pdf') as pdf:
 d=filter_display(pd.read_csv(SRC/'all_cell_identity_and_subtype_annotations.csv.gz',index_col=0),'All cell types')
 draw(d,'Postnatal heart — broad cell identities','reviewed_identity',('UMAP1','UMAP2'),'all_cell_types_labeled')
 draw(d,'Postnatal heart — cell identities and stage','reviewed_identity',('UMAP1','UMAP2'),'all_cell_types_stage_labeled',stage_colors=True)
 stage_facets(d,'Postnatal heart','reviewed_identity',('UMAP1','UMAP2'),'all_cell_types_stage_facets')
 for parent in ['Cardiomyocyte','Fibroblast','Endothelial','Mural','Immune']:
  d=filter_display(pd.read_csv(SRC/parent/'cell_subcluster_annotations.csv.gz',index_col=0),parent)
  for sub in [False,True]:draw(d,parent+' — subclusters and proposed identities','subcluster_id',('subcluster_UMAP1','subcluster_UMAP2'),parent+('_subtypes_labeled' if sub else '_subclusters_labeled'),sub)
  draw(d,parent+' — identities and stage','subcluster_id',('subcluster_UMAP1','subcluster_UMAP2'),parent+'_stage_labeled',stage_colors=True)
  stage_facets(d,parent,'subcluster_id',('subcluster_UMAP1','subcluster_UMAP2'),parent+'_stage_facets')
pd.DataFrame(exclusions).to_csv(OUT/'hidden_unresolved_assignments.csv.gz',index=False)
pd.DataFrame(counts).to_csv(OUT/'displayed_nuclei_by_stage.tsv',sep='\t',index=False)
manifest={'display_filter':'Exclude reviewed_identity or proposed_subtype containing unresolved or mixed (case-insensitive); candidate/review assignments retained. Exclusions apply to figures only.', 'method':'Display-only redraw; all existing coordinates, memberships and annotations unchanged. No geometric separation or UMAP refitting.','plots':records,'files':[]}
for p in sorted(OUT.glob('*')):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(f'Saved {len(records)} views, PNG/PDF, combined PDF and display-filter audit tables. Coordinates unchanged.')
