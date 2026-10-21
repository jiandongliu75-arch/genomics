"""Apply conservative nomenclature and export criteria/decision tables, without recomputing expression."""
from pathlib import Path
import json,pandas as pd,numpy as np,h5py
from anndata.io import read_elem,write_elem
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/cell_identities_subtypes'
old='Purkinje-like candidate';new='Conduction-associated CM candidate'
for path in out.rglob('*'):
 if not path.is_file():continue
 if path.name.endswith('.csv.gz'):
  frame=pd.read_csv(path,index_col=0)
  for col in ['proposed_subtype','candidate_subtype']:
   if col in frame:frame[col]=frame[col].replace({old:new})
  frame.to_csv(path,compression={'method':'gzip','compresslevel':1})
 elif path.suffix=='.tsv':
  frame=pd.read_csv(path,sep='\t')
  for col in ['proposed_subtype','candidate_subtype']:
   if col in frame:frame[col]=frame[col].replace({old:new})
  if 'score_margin' in frame and 'runner_up_identity' in frame:
   frame.loc[frame.runner_up_identity.eq('no other marker-supported lineage'),'score_margin']=np.nan
  if old in frame.columns:frame=frame.rename(columns={old:new})
  frame.to_csv(path,sep='\t',index=False)
 elif path.suffix=='.h5ad':
  with h5py.File(path,'r+') as f:
   obs=read_elem(f['obs'])
   obs['proposed_subtype']=obs.proposed_subtype.astype(str).replace({old:new})
   write_elem(f,'obs',obs)
# Annotate exactly why a proposal remains under review.
s=pd.read_csv(out/'subcluster_identity_evidence.tsv',sep='\t')
why=[]
for _,r in s.iterrows():
 reasons=[]
 if r.confidence=='unresolved':reasons.append('no panel passed the required multi-marker/anchor gate')
 else:
  if r.nuclei<50:reasons.append('fewer than 50 nuclei')
  if r.libraries_with_ge10_nuclei<2:reasons.append('fewer than two libraries with >=10 nuclei')
  if r.largest_library_fraction>.8:reasons.append('>80% from one library')
  if r.panel_score<.3:reasons.append('low positive panel score')
  if r.score_margin<.25:reasons.append('insufficient score separation from alternative panel')
  if r.panel_enrichment<.1:reasons.append('panel not enriched >=0.1 mean log units vs parent rest')
  if r.positive_within_library_contrasts<2:reasons.append('fewer than two positive within-library contrasts')
 why.append('; '.join(reasons) if reasons else 'all prespecified marker/library support rules passed')
s['review_reason']=why;s['marker_identity_evidence']=np.where(s.supported_markers.fillna('').ne(''),'multi-marker profile present','insufficient marker profile')
s.to_csv(out/'subcluster_identity_evidence.tsv',sep='\t',index=False)
resolution=[]
for parent in ['Cardiomyocyte','Fibroblast','Endothelial','Mural','Immune']:
 path=out/parent/'subclusters_metadata_only.h5ad';m=ad.read_h5ad(path)
 for res in [.3,.6,1.0]:resolution.append({'parent':parent,'resolution':res,'clusters':m.obs[f'leiden_r{res:g}'].nunique(),'nuclei':m.n_obs})
 # Refresh subtype plot legends after the conservative nomenclature change.
 xy=m.obsm['X_umap'];primary=m.obs['leiden_r0.6'].astype(str);groups=sorted(primary.unique(),key=int);fig,axes=plt.subplots(1,3,figsize=(20,6));cmap=plt.get_cmap('tab20')
 for j,c in enumerate(groups):
  keep=primary.eq(c).to_numpy();axes[0].scatter(xy[keep,0],xy[keep,1],s=2,color=cmap(j%20),rasterized=True);center=np.median(xy[keep],axis=0);axes[0].text(*center,str(c),ha='center',fontsize=8,bbox={'facecolor':'white','alpha':.8,'pad':1,'edgecolor':'none'})
 for j,label in enumerate(m.obs.proposed_subtype.unique()):
  keep=m.obs.proposed_subtype.eq(label).to_numpy();axes[1].scatter(xy[keep,0],xy[keep,1],s=2,label=label,color=cmap(j%20),rasterized=True)
 axes[1].legend(fontsize=6,markerscale=4,loc='best')
 for age in ['P0','P7','P14','P21']:
  keep=m.obs.age.eq(age).to_numpy();axes[2].scatter(xy[keep,0],xy[keep,1],s=2,label=age,rasterized=True)
 axes[2].legend(markerscale=4,fontsize=7)
 for ax,title in zip(axes,['Within-lineage subclusters','Proposed subtype (including candidates)','Age, no integration']):ax.set(title=title,xlabel='UMAP1',ylabel='UMAP2')
 fig.suptitle(parent);fig.tight_layout();fig.savefig(out/parent/'subcluster_umap.png',dpi=160);fig.savefig(out/parent/'subcluster_umap.pdf');plt.close(fig)
pd.DataFrame(resolution).to_csv(out/'resolution_cluster_counts.tsv',sep='\t',index=False)
cfg=json.loads((ROOT/'analysis/subtype_criteria.json').read_text());criteria=[]
for parent,defs in cfg['subtypes'].items():
 for d in defs:
  criteria.append({'parent':parent,'proposed_identity':d['name'],'positive_panel':';'.join(d['markers']),'required_anchors':';'.join(d['anchors']),'min_anchors':d['anchor_min'],'exclusion_checks':';'.join(d['exclude']),'interpretation_limit':d['note']})
pd.DataFrame(criteria).to_csv(out/'subtype_identity_criteria.tsv',sep='\t',index=False)
ann=pd.read_csv(out/'all_cell_identity_and_subtype_annotations.csv.gz',index_col=0)
for col in ['reviewed_identity','proposed_subtype','subtype_confidence']:
 pd.crosstab(ann[col],ann['sample']).to_csv(out/f'{col}_counts_by_sample.tsv',sep='\t')
s.groupby(['parent','proposed_subtype','confidence'],observed=True).agg(subclusters=('subcluster_id','size'),nuclei=('nuclei','sum')).reset_index().to_csv(out/'subtype_counts_and_confidence.tsv',sep='\t',index=False)
assert len(ann)==73052 and ann.index.is_unique and int(s.nuclei.sum())==int(ann.subcluster_id.ne('not_subclustered').sum())
print('Finalized annotations:',len(ann),'nuclei;',len(s),'subclusters')
print(s.groupby(['parent','proposed_subtype','confidence']).nuclei.sum().to_string())
