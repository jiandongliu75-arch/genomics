"""New presentation plates assembled from frozen, validated atlas results."""
from pathlib import Path
import os,textwrap,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/unified-heart-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/heart_atlas_unified_report';FIG=OUT/'figures';STAGES=['P0','P7','P14','P21']
COL={'Cardiomyocyte':'#4477AA','Endothelial':'#228833','Fibroblast':'#AA3377','Pericyte':'#AA9944','Smooth muscle':'#EE7733','Immune':'#CC3311','Epicardial':'#55AABB'}
EC={'Capillary-like':'#4477AA','Endocardial-like':'#DDAA33','Arterial-like':'#CC3311','Venous-like':'#228833','Lymphatic-like':'#AA3377'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def read(p,**kw):return pd.read_csv(ROOT/p,sep='\t',**kw)
def panel(ax,letter,title,width=42):ax.set_title(letter+'  '+textwrap.fill(title,width),loc='left',weight='bold',pad=10)
def save(fig,stem):
 FIG.mkdir(exist_ok=True,parents=True)
 for ext in ['png','pdf','svg']:
  p=FIG/(stem+'.'+ext);fig.savefig(p,dpi=300,facecolor='white')
  if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text().splitlines())+'\n')
 plt.close(fig);print('UNIFIED_FIGURE',stem,flush=True)
def mapbase(ax,f):
 ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
 for s in ax.spines.values():s.set_visible(False)
 x0,x1=f.x_centroid.min(),f.x_centroid.max();y0,y1=f.y_centroid.min(),f.y_centroid.max();x=x0+.06*(x1-x0);y=y1-.05*(y1-y0)
 ax.plot([x,x+500],[y,y],color='#333',lw=1.8);ax.text(x+250,y-.02*(y1-y0),'500 µm',ha='center',fontsize=7)
def atlas():
 fig,axs=plt.subplots(2,4,figsize=(14.5,8.8),layout='constrained')
 for j,stage in enumerate(STAGES):
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  ax=axs[0,j]
  for name,color in COL.items():
   z=f[f.atlas_identity.eq(name)]
   # Display-only sampling reduces visual saturation; all FB/EC estimates remain unchanged.
   if len(z)>40000:z=z.sample(40000,random_state=947)
   ax.scatter(z.x_centroid,z.y_centroid,s=.5,color=color,linewidths=0,rasterized=True)
  mapbase(ax,f);panel(ax,chr(97+j),stage+' · broad cardiac identities',27)
  ax=axs[1,j];z=f[f.atlas_resolved]
  if len(z)>50000:z=z.sample(50000,random_state=947)
  ax.scatter(z.x_centroid,z.y_centroid,s=.2,color='#e8e8e8',linewidths=0,rasterized=True)
  eligible=f.atlas_identity.eq('Endothelial')&f.program_source_eligible
  un=f[eligible&~f.EC_nuclear_type.isin(EC)];ax.scatter(un.x_centroid,un.y_centroid,s=1.1,c='#bdbdbd',linewidths=0,rasterized=True)
  for name,color in EC.items():
   z=f[eligible&f.EC_nuclear_type.eq(name)];ax.scatter(z.x_centroid,z.y_centroid,s=2,color=color,linewidths=0,rasterized=True)
  mapbase(ax,f);panel(ax,chr(101+j),stage+' · supported EC subtypes',27)
 fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,label=k,ms=4) for k,c in COL.items()],loc='center',bbox_to_anchor=(.5,.505),ncol=4,fontsize=8,frameon=False).set_in_layout(False)
 fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,label=k.replace('-like',''),ms=4) for k,c in EC.items()],loc='lower center',bbox_to_anchor=(.5,.005),ncol=5,fontsize=8,frameon=False).set_in_layout(False)
 fig.get_layout_engine().set(rect=(0,.045,1,.955),hspace=.17,wspace=.04)
 save(fig,'Figure_1_cardiac_lineages_and_endothelial_subtypes')
def vascular():
 fig=plt.figure(figsize=(11,8.5),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.4]);a=fig.add_subplot(gs[0,:]);a.axis('off');panel(a,'a','Coronary growth provides an anatomical framework',85)
 cards=[('Late fetal–P7','Existing coronary vessels expand\nespecially into inner ventricular muscle'),('By P7','Inner myocardial vessel density\napproaches outer-wall density'),('P7–P14 and later','Early angiogenic activity declines;\nfunctional integration remains to be measured')]
 for i,(title,text) in enumerate(cards):
  x=.015+i*.335;a.add_patch(FancyBboxPatch((x,.18),.30,.65,boxstyle='round,pad=.01',facecolor=['#eef3fa','#edf7ef','#f8edf5'][i],edgecolor='#789'))
  a.text(x+.15,.72,title,ha='center',weight='bold',fontsize=10);a.text(x+.15,.43,'\n'.join(textwrap.fill(s,31) for s in text.split('\n')),ha='center',va='center',fontsize=9)
  if i<2:a.annotate('',xy=(x+.33,.5),xytext=(x+.306,.5),arrowprops={'arrowstyle':'->','color':'#555'})
 a.text(.5,.03,'Published vascular sequence from Lu et al. (2021); not new perfusion measurements',ha='center',fontsize=9)
 ax=fig.add_subplot(gs[1,0]);d=read('results/cardiac_spatial_reconstruction/nuclear_cycling_specificity_coefficients.tsv');d=d[d.geometry.eq('within50um')&d.outcome.eq('CM_nuclear_cycling')].set_index('stage').loc[STAGES]
 ax.errorbar(range(4),d.standardized_coefficient,yerr=[d.standardized_coefficient-d.block_CI_low,d.block_CI_high-d.standardized_coefficient],fmt='o-',color='#228833',capsize=3)
 ax.axhline(0,color='#999',lw=.8);ax.set_xticks(range(4),STAGES);ax.set_ylabel('EC-cycling / CM-cycling association');panel(ax,'b','Local muscle and vascular growth occur together early')
 ax=fig.add_subplot(gs[1,1]);d=read('results/cardiac_spatial_reconstruction/primary_matched_component_results.tsv')
 for k,(outcome,label,color) in enumerate([('core_Z__withdrawal','Lower cell-cycle expression','#aa3377'),('core_Z__metabolic','Higher metabolic-regulatory expression','#4477aa')]):
  q=d[d.outcome.eq(outcome)].set_index('stage').loc[STAGES];x=np.arange(4)+(k-.5)*.15
  ax.errorbar(x,q.mean_standardized_difference,yerr=[q.mean_standardized_difference-q.block_CI_low,q.block_CI_high-q.mean_standardized_difference],fmt='o-',color=color,capsize=3,label=label)
 ax.axhline(0,color='#999',lw=.8);ax.set_xticks(range(4),STAGES);ax.set_ylabel('Capillary-adjacent minus matched comparison');ax.legend(frameon=False,fontsize=8,loc='upper left');panel(ax,'c','Capillary proximity relates to different CM components')
 save(fig,'Figure_6_coronary_growth_and_CM_components')
def transport():
 fig,axs=plt.subplots(2,2,figsize=(12,8.7),layout='constrained');d=read('results/cardiac_subtype_niches/CM_identity_restricted_Wormhole_assessment.tsv')
 labels=[s+'→'+t+'\n'+c.replace('.CM','')+' '+r for s,t,c,r in d[['stage_a','stage_b','CM_state','distance_stratum']].itertuples(index=False,name=None)]
 axs[0,0].bar(range(len(d)),d.Wormhole_exact_rank_correlation,color='#4477aa');axs[0,0].set_xticks(range(len(d)),labels,rotation=70,fontsize=7);axs[0,0].set_ylim(0,1);axs[0,0].set_ylabel('Rank correlation');panel(axs[0,0],'a','Approximate and exact distances agree variably')
 axs[0,1].bar(range(len(d)),d.nearest_exact_match_agreement,color='#aa3377');axs[0,1].set_xticks(range(len(d)),labels,rotation=70,fontsize=7);axs[0,1].set_ylim(0,1);axs[0,1].set_ylabel('Nearest-match agreement');panel(axs[0,1],'b','The closest matches are less stable')
 q=read('results/cardiac_subtype_niches/CM_maturation_rank_correspondence_audit.tsv');q=q[q.stage_a.eq('P7')&q.stage_b.eq('P14')];ax=axs[1,0];labels=[s.replace('.CM','')+' '+r for s,r in q[['CM_state','distance_stratum']].itertuples(index=False,name=None)]
 for i,(col,label,color) in enumerate([('exact_mean_absolute_maturation_rank_gap','Exact matches','#4477aa'),('Wormhole_mean_absolute_maturation_rank_gap','Wormhole matches','#aa3377'),('all_candidate_pairs_mean_absolute_rank_gap','All possible pairs','#999999')]):ax.bar(np.arange(len(q))+(i-1)*.24,q[col],width=.23,label=label,color=color)
 ax.set_xticks(range(len(q)),labels,rotation=25,fontsize=8);ax.set_ylabel('Mean difference in relative CM maturity rank');ax.legend(frameon=False,fontsize=8);panel(ax,'c','Matching preserves some regional ordering')
 q=read('results/cardiac_fibroblast_CM_refit/updated_frozen_transport_neighborhood_descriptors.tsv');ax=axs[1,1];ax.hist(q.eligible_FB_sources,bins=np.arange(-.5,q.eligible_FB_sources.max()+1.5),color='#009e73',edgecolor='white');ax.set_xlabel('Eligible FBs among sampled source cells');ax.set_ylabel('Observed neighborhoods');panel(ax,'d','Sparse fibroblast sampling limits fine dynamics')
 save(fig,'Figure_S3_transport_calibration_and_sampling')
if __name__=='__main__':atlas();vascular();transport()
