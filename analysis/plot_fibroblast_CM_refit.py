"""Integrated plates for the updated FB–CM biological report."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-CM-mpl')
from pathlib import Path
import json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
from refit_fibroblast_CM_relationships import ROOT,OUT,STAGES,FAMILIES,MODULES
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':300})
COLORS={'Postn':'#e69f00','Gfpt2':'#009e73','Col8a1':'#0072b2','Lsamp':'#cc79a7','Kirrel3':'#8b5c38','Unresolved':'#bdbdbd'}
MNAME={'Postn.Fib':'Postn-associated','Gfpt2.Fib':'Gfpt2-associated','Col8a1.Fib':'Col8a1-associated','Lsamp.Fib':'Lsamp-associated','Fibrillar_matrix':'Fibrillar matrix','Basement_matrix':'Basement matrix','Elastic_matrix':'Elastic matrix'}
def rd(p):return pd.read_csv(ROOT/p,sep='\t')
R=rd('results/cardiac_fibroblast_CM_refit/CM_fibroblast_refitted_coefficients.tsv')
SAME=rd('results/cardiac_fibroblast_CM_refit/same_CM_radius_and_boundary_sensitivity.tsv')
VASC=rd('results/cardiac_fibroblast_CM_refit/same_CM_vascular_geometry_sensitivity.tsv')
NO_CM=rd('results/cardiac_fibroblast_CM_refit/CM_gene_excluded_fibroblast_program_sensitivity.tsv')
PROG=rd('results/cardiac_fibroblast_CM_refit/fibroblast_per_profile_programs.tsv.gz')
COVER=rd('results/cardiac_fibroblast_CM_refit/fibroblast_assignment_coverage.tsv')
PCOV=rd('results/cardiac_fibroblast_CM_refit/fibroblast_program_coverage.tsv')
INVENTORY=[]
def finish(fig,stem,title):
 for ext in ['png','pdf','svg']:
  fig.savefig(FIG/(stem+'.'+ext),dpi=300,facecolor='white')
 if (FIG/(stem+'.svg')).exists():
  p=FIG/(stem+'.svg');p.write_text('\n'.join(s.rstrip() for s in p.read_text().splitlines())+'\n')
 plt.close(fig);INVENTORY.append({'stem':stem,'title':title});print('REFIT_FIGURE',stem,flush=True)
def panel(ax,label,title):
 columns=ax.get_subplotspec().get_gridspec().ncols
 width={1:90,2:40,3:32,4:25}.get(columns,40)
 ax.set_title(label+'  '+textwrap.fill(title,width=width),loc='left',fontweight='bold',pad=9)
def frame(ax):
 ax.spines[['top','right']].set_visible(False)
def spatial(ax,stage,values=None,labels=None,vmin=None,vmax=None):
 f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',sep='\t',index_col=0)
 lab=pd.read_csv(ROOT/f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',sep='\t',index_col=0)
 source=lab[lab.program_source_eligible]
 bg=f[f.atlas_resolved]
 if len(bg)>60000:bg=bg.sample(60000,random_state=947)
 ax.scatter(bg.x_centroid,bg.y_centroid,s=.3,c='#e9e9e9',rasterized=True,linewidths=0)
 if labels is not None:
  for fam in FAMILIES+['Unresolved']:
   use=source.recommended_candidate.eq(fam+'-associated') if fam!='Unresolved' else source.recommended_candidate.eq('Unresolved')
   ax.scatter(source.loc[use,'x_centroid'],source.loc[use,'y_centroid'],s=2,c=COLORS[fam],rasterized=True,linewidths=0)
  im=None
 else:
  val=values.reindex(source.index)
  im=ax.scatter(source.x_centroid,source.y_centroid,c=val,s=2,cmap='viridis',vmin=vmin,vmax=vmax,rasterized=True,linewidths=0)
 ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
 for s in ax.spines.values():s.set_visible(False)
 xmin,xmax=ax.get_xlim();ymax,ymin=ax.get_ylim();x=xmin+.06*(xmax-xmin);y=ymax-.05*(ymax-ymin)
 ax.plot([x,x+500],[y,y],c='#333333',lw=2);ax.text(x+250,y-.02*(ymax-ymin),'500 µm',ha='center',va='bottom',fontsize=7)
 return im

def overview():
 fig,axs=plt.subplots(2,3,figsize=(13,8.6),layout='constrained',gridspec_kw={'width_ratios':[1,1,.95]})
 for ax,stage,letter in zip([axs[0,0],axs[0,1],axs[1,0],axs[1,1]],STAGES,'abcd'):
  spatial(ax,stage,labels=True);panel(ax,letter,stage+' fibroblast candidates')
 axs[1,0].legend(handles=[Patch(facecolor=COLORS[k],label=k) for k in FAMILIES+['Unresolved']],loc='lower left',bbox_to_anchor=(-.05,-.20),ncol=3,frameon=False,fontsize=8)
 ax=axs[0,2];a=COVER[COVER.label_layer.eq('candidate')].pivot(index='family',columns='stage',values='candidate_profiles').reindex(index=FAMILIES,columns=STAGES)
 m=ax.imshow(np.ma.masked_invalid(a.to_numpy()),cmap='Blues',aspect='auto',vmin=0,vmax=6000);ax.set_xticks(range(4),STAGES);ax.set_yticks(range(5),FAMILIES)
 for i in range(5):
  for j in range(4):ax.text(j,i,'N/A' if np.isnan(a.iloc[i,j]) else f'{int(a.iloc[i,j]):,}',ha='center',va='center',fontsize=8,color='white' if a.iloc[i,j]>3500 else 'black')
 panel(ax,'e','Assignments require an available reference');ax.text(.0,-.16,'N/A: family not assessable with that age’s reference\nCounts show coverage, not tissue abundance.',transform=ax.transAxes,fontsize=8,va='top')
 ax=axs[1,2];a=PCOV[PCOV.view.eq('refined')].pivot(index='program',columns='stage',values='mean_rank_excess').reindex(index=MODULES,columns=STAGES)
 im=ax.imshow(a.to_numpy(),cmap='coolwarm',vmin=-.3,vmax=.3,aspect='auto');ax.set_xticks(range(4),STAGES);ax.set_yticks(range(len(a)),[MNAME[m] for m in a.index],fontsize=8);panel(ax,'f','Gene patterns can be measured at every age')
 fig.colorbar(im,ax=ax,fraction=.04,pad=.02,label='Mean marker-rank excess')
 finish(fig,'Figure_1_fibroblast_states_and_shared_programs','Fibroblast identities and shared patterns give complementary maps')

def program_maps():
 fig,axs=plt.subplots(3,4,figsize=(13.5,10.8),layout='constrained')
 for i,mod in enumerate(['Postn.Fib','Gfpt2.Fib','Col8a1.Fib']):
  p=PROG[PROG.view.eq('refined')];vmin,vmax=np.quantile(p['program__'+mod],[.03,.97])
  for j,stage in enumerate(STAGES):
   val=p[p.stage.eq(stage)].set_index('cell_id')['program__'+mod];im=spatial(axs[i,j],stage,values=val,vmin=vmin,vmax=vmax)
   panel(axs[i,j],chr(97+i*4+j),stage+' · '+MNAME[mod])
  fig.colorbar(im,ax=axs[i,:],fraction=.014,pad=.01,label='Marker-pattern score')
 finish(fig,'Figure_2_measured_fibroblast_program_maps','Early and late fibroblast patterns overlap without requiring a subtype label')

def heatmaps():
 fig,axs=plt.subplots(3,1,figsize=(13.5,11.5),layout='constrained')
 primary=R[R.primary_analysis&R.source_kind.eq('program')]
 cols=[(s,c) for s in STAGES for c in ['v.CM','Ankrd1.CM','p.CM']]
 for ax,letter,outcome,title in zip(axs,'abc',['core_Z__metabolic','core_Z__sarcomere','core_Z__cycling'],['CM energy-metabolism expression','CM sarcomere expression','CM cell-cycle expression']):
  a=np.full((7,12),np.nan);ci=np.zeros_like(a,dtype=bool)
  for i,mod in enumerate(MODULES):
   for j,(stage,state) in enumerate(cols):
    r=primary[primary.stage.eq(stage)&primary.CM_state.eq(state)&primary.target.eq('refined__'+mod)&primary.outcome.eq(outcome)]
    if len(r):a[i,j]=r.coefficient.iloc[0];ci[i,j]=r.CI_low.iloc[0]*r.CI_high.iloc[0]>0
  cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e7e7e7');im=ax.imshow(np.ma.masked_invalid(a),aspect='auto',cmap=cmap,vmin=-.12,vmax=.12)
  ax.set_yticks(range(7),[MNAME[m] for m in MODULES]);ax.set_xticks(range(12),[s+'\n'+{'v.CM':'Working','Ankrd1.CM':'Ankrd1','p.CM':'Cycling'}[c] for s,c in cols],fontsize=8)
  for x in [2.5,5.5,8.5]:ax.axvline(x,color='white',lw=2)
  for y,x in zip(*np.where(ci)):ax.text(x,y,'•',ha='center',va='center',fontsize=10)
  panel(ax,letter,title);fig.colorbar(im,ax=ax,fraction=.015,pad=.01,label='Adjusted association')
 finish(fig,'Figure_3_local_program_relationships','Local matrix relationships differ by CM state and outcome')

def pick(stage,state,target,outcome,layer=None,adjust='within_region',geom='proximity50',data=R):
 mask=data.stage.eq(stage)&data.CM_state.eq(state)&data.target.eq(target)&data.outcome.eq(outcome)&data.geometry.eq(geom)&data.adjustment.eq(adjust)
 if layer is not None:mask&=data.source_layer.eq(layer)
 x=data[mask];assert len(x)<=1,(stage,state,target,outcome,layer,geom,len(x));return None if len(x)==0 else x.iloc[0]
def forest(ax,rows):
 valid=[(label,r,color) for label,r,color in rows if r is not None]
 for i,(label,r,c) in enumerate(valid):ax.errorbar(r.coefficient,i,xerr=[[max(0,r.coefficient-r.CI_low)],[max(0,r.CI_high-r.coefficient)]],fmt='o',color=c,capsize=2,ms=4,lw=1.3)
 ax.set_yticks(range(len(valid)),[x[0] for x in valid],fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#999999',lw=.8);ax.set_xlabel('Adjusted association (standardized units)');frame(ax)
def early_relationships():
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 rows=[(label,pick('P7','p.CM',layer+'__Postn','core_Z__metabolic',layer),'#e69f00') for layer,label in [('legacy','Earlier stringent labels'),('candidate','Broader candidates'),('crosschecked','Cell/core + hidden-gene support')]]
 forest(axs[0,0],rows);panel(axs[0,0],'a','Postn identity: the earlier association weakens')
 rows=[(label,pick('P7','p.CM',view+'__Postn.Fib','core_Z__metabolic',view),'#d79500') for view,label in [('refined','Whole-cell program'),('core','Inner-cell program'),('native','Deposited-cell program'),('depthmatched','Random RNA loss to core depth')]]
 rows.append(('CM genes excluded from FB score',pick('P7','p.CM','noCM_refined__Postn.Fib','core_Z__metabolic','noCM_refined',data=NO_CM),'#d79500'))
 forest(axs[0,1],rows);panel(axs[0,1],'b','Postn pattern: negative estimates remain tentative')
 rows=[(label,pick('P7','p.CM',view+'__Gfpt2.Fib','core_Z__metabolic',view),'#009e73') for view,label in [('refined','Whole-cell program'),('core','Inner-cell program'),('native','Deposited-cell program'),('depthmatched','Random RNA loss to core depth')]]
 rows.append(('CM genes excluded from FB score',pick('P7','p.CM','noCM_refined__Gfpt2.Fib','core_Z__metabolic','noCM_refined',data=NO_CM),'#009e73'))
 forest(axs[1,0],rows);panel(axs[1,0],'c','Gfpt2 pattern: a possible contrasting environment')
 rows=[]
 for mod,color in [('Postn.Fib','#d79500'),('Gfpt2.Fib','#009e73')]:
  for geom,label in [('proximity50','50 µm'),('proximity30','30 µm')]:
   data=SAME[SAME.paired_geometry.eq('proximity30')];r=pick('P7','p.CM','refined__'+mod,'core_Z__metabolic','refined',geom=geom,data=data)
   rows.append((MNAME[mod]+' · '+label,r,color))
 forest(axs[1,1],rows);panel(axs[1,1],'d','Smaller neighborhoods: same 2,753 cycling CMs')
 finish(fig,'Figure_4_early_fibroblast_CM_relationships','Early fibroblast environments suggest hypotheses without establishing metabolic instruction')

def structural_relationships():
 fig,axs=plt.subplots(2,2,figsize=(12.5,8.8),layout='constrained')
 rows=[]
 for state,label,color in [('Ankrd1.CM','Ankrd1 CM','#a04a9e'),('v.CM','Working CM','#2670a0')]:
  for outcome,desc in [('core_Z__sarcomere','Sarcomere'),('core_Z__sarcomere_noMyl2','Sarcomere without Myl2'),('core_Z__maturation21','Overall 21-gene score')]:rows.append((label+' · '+desc,pick('P14',state,'EC__Endocardial-like',outcome,'candidate'),color))
 forest(axs[0,0],rows);panel(axs[0,0],'a','P14 endocardial setting relates to a CM component')
 rows=[('Core CM RNA; FB identities',pick('P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere','candidate'),'#a04a9e'),
       ('Core CM RNA; FB programs',pick('P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere','programs'),'#a04a9e'),
       ('Deposited CM RNA',pick('P14','Ankrd1.CM','EC__Endocardial-like','native_Z__sarcomere','candidate'),'#777777')]
 for geom,label in [('proximity50','50 µm'),('proximity30','30 µm')]:
  d=VASC[VASC.paired_geometry.eq('proximity30')];r=pick('P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere','candidate',geom=geom,data=d);rows.append(('Same CMs · '+label,r,'#a04a9e'))
 forest(axs[0,1],rows);panel(axs[0,1],'b','The result depends on how CM RNA is assigned')
 rows=[(label,pick('P0','v.CM',view+'__Elastic_matrix','core_Z__sarcomere',view),'#009c9c') for view,label in [('refined','Whole-cell FB program'),('core','Inner-cell FB program'),('native','Deposited-cell FB program'),('depthmatched','Random RNA loss to core depth')]]
 for view,label in [('noCM_refined','CM genes excluded · whole FB'),('noCM_core','CM genes excluded · inner FB')]:rows.append((label,pick('P0','v.CM',view+'__Elastic_matrix','core_Z__sarcomere',view,data=NO_CM),'#009c9c'))
 forest(axs[1,0],rows);panel(axs[1,0],'c','P0 elastic-matrix setting: a candidate to test')
 rows=[(label,pick('P21','Ankrd1.CM',view+'__Basement_matrix','core_Z__sarcomere',view),'#0072b2') for view,label in [('refined','Whole-cell FB program'),('core','Inner-cell FB program'),('native','Deposited-cell FB program'),('depthmatched','Random RNA loss to core depth')]]
 forest(axs[1,1],rows);panel(axs[1,1],'d','P21 basement matrix: weaker inner-cell evidence')
 finish(fig,'Figure_5_structural_relationships_and_boundaries','Structural RNA relationships depend on CM identity and transcript ownership')

def spatial_and_time():
 fig,axs=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
 old=rd('results/fibroblast_names_and_neighborhoods/matched_family_context_sensitivity.tsv')
 new=rd('results/cardiac_fibroblast_CM_refit/updated_matched_fibroblast_spatial_context.tsv')
 rows=[]
 for layer,label in [('legacy','Earlier stringent'),('candidate','Broader candidates'),('crosschecked','Crosschecked candidates')]:
  d=old if layer=='legacy' else new[new.new_label_layer.eq(layer)]
  r=d[d.stage.eq('P14')&d.measure.eq('arterial_within_50um')].iloc[0]
  rows.append((label,100*r.target_matched_mean,100*r.reference_matched_mean))
 ax=axs[0,0];y=np.arange(3)
 ax.barh(y-.16,[r[1] for r in rows],height=.3,color=COLORS['Col8a1'],label='Col8a1');ax.barh(y+.16,[r[2] for r in rows],height=.3,color=COLORS['Gfpt2'],label='Gfpt2')
 ax.set_yticks(y,[r[0] for r in rows]);ax.invert_yaxis();ax.set_xlabel('Within 50 µm of an arterial structure (%)');ax.legend(frameon=False,fontsize=8);frame(ax);panel(ax,'a','P14 arterial association is smaller with broader labels')
 d=rd('results/cardiac_fibroblast_CM_refit/updated_fibroblast_spatial_and_program_summary.tsv');d=d[d.label_layer.eq('candidate')&d.measure.eq('arterial_within_50um')]
 a=d.pivot(index='family',columns='stage',values='mean').reindex(index=[f+'-associated' for f in FAMILIES],columns=STAGES)*100
 ax=axs[0,1];im=ax.imshow(np.ma.masked_invalid(a.to_numpy()),aspect='auto',cmap='Blues',vmin=0,vmax=25);ax.set_xticks(range(4),STAGES);ax.set_yticks(range(5),FAMILIES)
 for i in range(5):
  for j in range(4):ax.text(j,i,'N/A' if pd.isna(a.iloc[i,j]) else f'{a.iloc[i,j]:.0f}%',ha='center',va='center',fontsize=8)
 panel(ax,'b','Arterial proximity is one part of fibroblast location')
 u=rd('results/cardiac_fibroblast_CM_refit/updated_transport_program_comparisons.tsv');u=u[u.stage_a.eq('P7')&u.stage_b.eq('P14')&u.method.eq('exact')]
 ax=axs[1,0];contexts=[('v.CM','≤100µm'),('v.CM','>100µm'),('Ankrd1.CM','≤100µm'),('Ankrd1.CM','>100µm')]
 for k,(mod,color) in enumerate([('Postn.Fib',COLORS['Postn']),('Gfpt2.Fib',COLORS['Gfpt2']),('Col8a1.Fib',COLORS['Col8a1'])]):
  vals=[u[u.CM_state.eq(state)&u.distance_stratum.eq(dist)&u.program.eq(mod)].mean_target_minus_source.iloc[0] for state,dist in contexts]
  ax.plot(vals,np.arange(4)+(k-1)*.15,'o',color=color,label=MNAME[mod],ms=5)
 ax.axvline(0,color='#999',lw=.8);ax.set_yticks(range(4),['Working · near lining','Working · farther away','Ankrd1 · near lining','Ankrd1 · farther away'],fontsize=8);ax.invert_yaxis();ax.set_xlabel('P14 minus P7 pattern score in frozen exact matches');ax.legend(frameon=False,fontsize=7,loc='upper right');frame(ax);panel(ax,'c','Matched surroundings change from P7 to P14')
 t=rd('results/cardiac_fibroblast_CM_refit/updated_frozen_transport_neighborhood_descriptors.tsv');ax=axs[1,1]
 for i,stage in enumerate(STAGES):
  for state,color,offset in [('v.CM','#2670a0',-.12),('Ankrd1.CM','#a04a9e',.12)]:
   z=t[t.stage.eq(stage)&t.CM_state.eq(state)];jitter=np.linspace(-.07,.07,len(z));ax.scatter(i+offset+jitter,z.eligible_FB_sources,s=20,color=color,alpha=.7,label=state if i==0 else None)
 ax.set_xticks(range(4),STAGES);ax.set_ylabel('Eligible FBs per sampled neighborhood');ax.set_ylim(-.4,max(t.eligible_FB_sources)+1);ax.legend(frameon=False,fontsize=8);frame(ax);panel(ax,'d','Few FBs limit the resolution of temporal matches')
 finish(fig,'Figure_6_spatial_context_and_cross_age_comparisons','Fibroblast position and changing surroundings guide specific developmental questions')

def identity_supplement():
 fig,axs=plt.subplots(3,1,figsize=(13.5,10),layout='constrained');cols=[(s,c) for s in STAGES for c in ['v.CM','Ankrd1.CM','p.CM']]
 r=R[R.primary_analysis&R.source_kind.eq('identity')]
 for ax,letter,outcome,title in zip(axs,'abc',['core_Z__metabolic','core_Z__sarcomere','core_Z__cycling'],['Metabolic expression','Sarcomere expression','Cell-cycle expression']):
  a=np.full((5,12),np.nan)
  for i,fam in enumerate(FAMILIES):
   for j,(stage,state) in enumerate(cols):
    z=r[r.stage.eq(stage)&r.CM_state.eq(state)&r.target.eq('candidate__'+fam)&r.outcome.eq(outcome)]
    if len(z):a[i,j]=z.coefficient.iloc[0]
  cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e7e7e7');im=ax.imshow(np.ma.masked_invalid(a),cmap=cmap,vmin=-.12,vmax=.12,aspect='auto');ax.set_yticks(range(5),FAMILIES);ax.set_xticks(range(12),[s+'\n'+{'v.CM':'Working','Ankrd1.CM':'Ankrd1','p.CM':'Cycling'}[c] for s,c in cols],fontsize=8)
  for x in [2.5,5.5,8.5]:ax.axvline(x,c='white',lw=2)
  panel(ax,letter,title);fig.colorbar(im,ax=ax,fraction=.015,pad=.01,label='Adjusted association')
 finish(fig,'Figure_S1_fibroblast_identity_associations','Subtype labels do not establish one uniform effect on neighboring CMs')

if __name__=='__main__':
 overview();program_maps();heatmaps();early_relationships();structural_relationships();spatial_and_time();identity_supplement()
 (OUT/'figure_inventory.json').write_text(json.dumps(INVENTORY,indent=2)+'\n')
