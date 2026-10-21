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
records=[]
def draw(d,title,group,coords,stem,subtype_colors=False):
 groups=sorted(d[group].unique(),key=lambda x:int(x.split('_')[-1]) if group=='subcluster_id' else x)
 cats=sorted(d.proposed_subtype.unique()) if subtype_colors else groups
 colors={c:PALETTE[i%len(PALETTE)] for i,c in enumerate(cats)}
 fig=plt.figure(figsize=(19,12));ax=fig.add_axes([.29,.14,.42,.72])
 xy=d[list(coords)].to_numpy(); assert np.isfinite(xy).all()
 span=np.ptp(xy,axis=0);lo=xy.min(axis=0)-span*.07;hi=xy.max(axis=0)+span*.07
 centers=[]
 for g in groups:
  s=d[d[group]==g];v=s[list(coords)].to_numpy();center=np.median(v,axis=0)
  # Use an observed point near the median, so the connector touches the group.
  anchor=v[np.argmin(np.sum(((v-center)/span)**2,axis=1))]
  category=s.proposed_subtype.iloc[0] if subtype_colors else g;c=colors[category]
  ax.scatter(v[:,0],v[:,1],s=2.4 if len(d)>20000 else 6,c=c,alpha=.8,linewidths=0,rasterized=True)
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
 fig.text(.5,.91,'Color: proposed subtype' if subtype_colors else 'Color: computational subcluster' if group=='subcluster_id' else 'Color: reviewed broad cell identity',ha='center',fontsize=14)
 fig.text(.5,.055,'Existing UMAP coordinates preserved · Labels are provisional RNA-based identities\n[review] = additional uncertainty; candidate/unresolved labels retained · Separate lineage UMAPs have independent coordinate spaces',ha='center',fontsize=11,color='#444444')
 fig.savefig(OUT/(stem+'.png'),dpi=200,facecolor='white');fig.savefig(OUT/(stem+'.pdf'),facecolor='white');pdf.savefig(fig);plt.close(fig)
 records.append({'plot':stem,'nuclei':len(d),'groups':len(groups),'coordinates':list(coords),'source':str((SRC/('all_cell_identity_and_subtype_annotations.csv.gz' if group!='subcluster_id' else title.split(' — ')[0]+'/cell_subcluster_annotations.csv.gz')).relative_to(ROOT))})
with PdfPages(OUT/'all_labeled_umaps.pdf') as pdf:
 d=pd.read_csv(SRC/'all_cell_identity_and_subtype_annotations.csv.gz',index_col=0)
 draw(d,'Postnatal heart — broad cell identities','reviewed_identity',('UMAP1','UMAP2'),'all_cell_types_labeled')
 for parent in ['Cardiomyocyte','Fibroblast','Endothelial','Mural','Immune']:
  d=pd.read_csv(SRC/parent/'cell_subcluster_annotations.csv.gz',index_col=0)
  for sub in [False,True]:draw(d,parent+' — subclusters and proposed identities','subcluster_id',('subcluster_UMAP1','subcluster_UMAP2'),parent+('_subtypes_labeled' if sub else '_subclusters_labeled'),sub)
manifest={'method':'Display-only redraw; all existing coordinates, memberships and annotations unchanged. No geometric separation or UMAP refitting.','plots':records,'files':[]}
for p in sorted(OUT.glob('*')):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('Saved 11 labeled views, PNG/PDF, and combined PDF. Coordinates unchanged.')
