"""Integrated plates for the observed-Xenium / replicated-RNA analysis."""
import os,json,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/joint-signaling-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from joint_signaling_common import *

FIG=OUT/'figures'
LABELS={'maturation21':'Overall 21-gene','metabolic':'Metabolic-regulatory','sarcomere':'Sarcomere',
        'excitation_contraction':'Excitation–contraction','cycling':'Cell-cycle RNA'}
COLORS={'A':'#777777','B':'#CCBB44','C':'#4477AA','C_fine_context':'#228833','D':'#EE7733','E':'#AA3377'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,'pdf.fonttype':42,
 'axes.spines.top':False,'axes.spines.right':False})
LEGENDS=[]
def title(ax,letter,text,width=25):ax.set_title(letter+'  '+textwrap.fill(text,width),loc='left',weight='bold',pad=10)
def finish(fig,stem,title_,legend):
 FIG.mkdir(parents=True,exist_ok=True)
 for ext in ['png','pdf']:fig.savefig(FIG/(stem+'.'+ext),dpi=300,facecolor='white')
 plt.close(fig);LEGENDS.append({'stem':stem,'title':title_,'legend':legend});print('JOINT_FIGURE',stem,flush=True)
def heat(ax,mat,labels,cmap='Blues',norm=None,fmt=None):
 im=ax.imshow(np.ma.masked_invalid(mat),aspect='auto',cmap=cmap,norm=norm);im.cmap.set_bad('#ececec')
 ax.set_xticks(range(4),STAGES);ax.set_yticks(range(len(labels)),labels,fontsize=10)
 if fmt:
  for y,x in np.ndindex(mat.shape):
   ax.text(x,y,fmt(mat[y,x]) if np.isfinite(mat[y,x]) else '—',ha='center',va='center',fontsize=9.5)
 return im
def overview():
 selected=read(OUT/'selected_RNA_supported_pairs.tsv');support=read(OUT/'RNA_pair_source_support.tsv.gz')
 support=support[support.pair_id.isin(selected.pair_id)]
 fig,ax=plt.subplots(2,3,figsize=(11,7),layout='constrained')
 a=ax[0,0];a.axis('off');title(a,'a','Independent roles of the assays')
 blocks=[('snRNA-seq · 2 libraries per age','All subunits expressed in both libraries;\nCM receptors and source identity'),
         ('Xenium · observed RNA and locations','Local source expression, CM receptors,\n21-gene and component outcomes'),
         ('Four buffered test domains','No shared source cells across folds;\n≥400 µm between training and test CMs')]
 for j,(heading,body) in enumerate(blocks):
  y=.85-j*.31;a.text(.02,y,heading,weight='bold',transform=a.transAxes,color='#203C54')
  a.text(.02,y-.07,body,va='top',transform=a.transAxes,linespacing=1.4)
 a=ax[0,1];s=read(OUT/'RNA_pair_stage_support.tsv');counts=[s[s.stage.eq(g)&s.RNA_supported&s.pair_id.isin(selected.pair_id)].shape[0] for g in STAGES]
 a.bar(STAGES,counts,color='#4477AA');a.set_ylim(0,64);a.set_ylabel('RNA-supported candidate pairs')
 for i,n in enumerate(counts):a.text(i,n+1,str(n),ha='center')
 title(a,'b','RNA-supported candidate pairs')
 a=ax[0,2];program=read(OUT/'RNA_CM_programs_by_library.tsv');program=program[program.CM_state.eq('All')]
 for name,color in [('Cell_cycle_full','#AA3377'),('Sarcomere_mechanics','#4477AA'),('Oxidative_metabolism','#228833'),('Fatty_acid_oxidation','#EE7733')]:
  z=program[program.program.eq(name)];m=z.groupby('stage').mean_log1p_CP10k.mean().reindex(STAGES)
  # Each panel has its own gene baseline; display proportional change from P0.
  display={'Cell_cycle_full':'Cell-cycle','Sarcomere_mechanics':'Sarcomere/cytoskeleton','Oxidative_metabolism':'Oxidative','Fatty_acid_oxidation':'Fatty-acid'}
  base=m.iloc[0];a.plot(range(4),m/base,'o-',label=display[name],color=color,ms=4)
  for i,g in enumerate(STAGES):a.scatter(np.full(2,i),z[z.stage.eq(g)].mean_log1p_CP10k/base,s=13,color=color,alpha=.6)
 a.set_xticks(range(4),STAGES);a.set_ylabel('Program / own P0 mean');a.legend(fontsize=9.5,frameon=False)
 title(a,'c','CM programs change differently')
 expression=read(OUT/'RNA_complex_expression_by_library.tsv.gz');expression=expression[expression.view.eq('corrected')]
 a=ax[1,0]
 for gene,color in [('Tnc','#AA3377'),('Fn1','#EE7733'),('Col1a1','#4477AA'),('Lamc1','#228833')]:
  z=expression[expression.source_group.eq('Broad:Fibroblast')&expression.species.eq(gene)]
  m=z.groupby('stage').mean_log1p_CP10k.mean().reindex(STAGES);baseline=m.iloc[0]
  a.plot(range(4),m/baseline,'o-',color=color,label=gene,ms=4)
  for i,stage in enumerate(STAGES):a.scatter(np.full(2,i),z[z.stage.eq(stage)].mean_log1p_CP10k/baseline,s=13,color=color,alpha=.6)
 a.set_xticks(range(4),STAGES);a.set_ylabel('Mean log-RNA / own P0 mean');a.legend(fontsize=10,frameon=False,ncol=2)
 title(a,'d','Fibroblast matrix expression changes')
 a=ax[1,1]
 for group,gene,color,label in [('Broad:Endothelial','Nrg1','#4477AA','ECs: Nrg1'),('CM','Erbb4','#AA3377','CMs: Erbb4')]:
  z=expression[expression.source_group.eq(group)&expression.species.eq(gene)]
  m=z.groupby('stage').fraction_complex_detected.mean().reindex(STAGES)
  a.plot(range(4),m,'o-',color=color,label=label,ms=4)
  for i,stage in enumerate(STAGES):a.scatter(np.full(2,i),z[z.stage.eq(stage)].fraction_complex_detected,s=13,color=color,alpha=.6)
 a.set_xticks(range(4),STAGES);a.set_ylim(0,1.08);a.set_ylabel('Nuclei with target RNA (fraction)');a.legend(fontsize=9.5,frameon=False)
 title(a,'e','A vascular candidate has CM support')
 a=ax[1,2];d=read(OUT/'P7/receiver_domains.tsv.gz',index_col=0)
 a.scatter(d.x_centroid,d.y_centroid,s=.3,c='#e5e5e5',linewidths=0,rasterized=True)
 z=d[d.analysis_eligible];a.scatter(z.x_centroid,z.y_centroid,c=z.domain,cmap='tab10',s=.7,linewidths=0,rasterized=True,vmin=0,vmax=9)
 a.set_aspect('equal');a.invert_yaxis();a.set_xticks([]);a.set_yticks([])
 a.text(.02,.02,f'{len(z):,} retained P7 CMs',transform=a.transAxes,fontsize=10);title(a,'f','P7 · buffered domains',30)
 fig.suptitle('Muscle development and matrix composition change on different schedules',fontsize=14,weight='bold')
 finish(fig,'Figure_1_joint_assay_design','CM development and fibroblast matrix programs follow different schedules',
  'a, Roles of the assays and spatial validation. b, Candidate pairs that pass the two-library RNA screen at each age, after outcome-independent pathway reduction. Counts reflect eligibility, not signaling activity. c, Canonical CM RNA programs, excluding candidate ligand/receptor and spatial outcome genes; lines show the two-library mean and points the libraries, scaled to each program’s own P0 mean. These are RNA programs, not flux measurements. d, Broad-fibroblast matrix transcripts on the same full-genome normalization; ratios describe changes in mean log-expression, not molecule fold changes or deposited protein. e, Endothelial Nrg1 and ventricular-CM Erbb4 detection in both RNA libraries; availability does not establish signaling or a particular ERBB receptor assembly. f, Eligible P7 ventricular CMs in the four source-disjoint test domains; gray marks excluded boundary-buffer CMs. All other ages use the same design.')

def spatial_maps():
 fig=plt.figure(figsize=(11,9),layout='constrained');gs=fig.add_gridspec(3,5,width_ratios=[1,1,1,1,.06])
 definitions=[('outcome','21-gene score'),('cot__Broad:Fibroblast__Lamc1--Dag1','FB laminin input'),('cot__Broad:Endothelial__Nrg1--Erbb4','EC neuregulin input')]
 for row,(feature,label) in enumerate(definitions):
  for j,stage in enumerate(STAGES):
   a=fig.add_subplot(gs[row,j]);f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
   bg=f.iloc[::max(1,len(f)//35000)];a.scatter(bg.x_centroid,bg.y_centroid,s=.15,c='#ededed',linewidths=0,rasterized=True)
   with np.load(OUT/stage/'core_predictor_features.npz') as z:
    take=z['eligible'];xy=z['xy'][take];names=list(z['columns'])
    values=z['outcomes'][take,0] if feature=='outcome' else z['features'][take,names.index(feature)] if feature in names else None
   if values is not None:
    # Ranks show local patterns, avoiding claims of calibrated inter-age signaling mass.
    positive=np.ones(len(values),bool) if feature=='outcome' else values>1e-12
    if feature!='outcome':a.scatter(xy[:,0],xy[:,1],c='#d1d1d1',s=1.1,linewidths=0,rasterized=True)
    ranks=pd.Series(values[positive]).rank(pct=True).to_numpy()
    a.scatter(xy[positive,0],xy[positive,1],c=ranks,cmap='viridis',norm=Normalize(0,1),s=1.1,linewidths=0,rasterized=True)
    if feature!='outcome':a.text(.02,.94,f'{positive.mean():.0%} positive',transform=a.transAxes,fontsize=9.5)
   else:a.text(.5,.5,'Not RNA-supported\nfor this source and age',ha='center',va='center',transform=a.transAxes)
   a.set_aspect('equal');a.invert_yaxis();a.set_xticks([]);a.set_yticks([])
   for spine in a.spines.values():spine.set_visible(False)
   x=f.x_centroid.min()+.06*np.ptp(f.x_centroid);y=f.y_centroid.max()-.04*np.ptp(f.y_centroid)
   a.plot([x,x+500],[y,y],color='#333',lw=1.7);a.text(x+250,y-.015*np.ptp(f.y_centroid),'500 µm',ha='center',fontsize=9.5)
   title(a,chr(97+row*4+j),stage+' · '+label,19)
  fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=fig.add_subplot(gs[row,4]),label='Within-age percentile' if row==0 else 'Percentile among positive inputs')
 fig.suptitle('Measured CM development alongside source-specific molecular candidates',fontsize=14,weight='bold')
 finish(fig,'Figure_2_observed_CM_and_source_maps','CM development alongside source-specific molecular candidates',
  'a–d, Directly observed core-RNA 21-gene score, ranked within age. e–h, FB-associated Lamc1–Dag1 collective-transport input. i–l, EC-associated Nrg1–Erbb4 input. Candidates were chosen on biological grounds before inspecting prediction gains. Candidate colors rank inputs above 10⁻¹² within each age; gray CMs have no allocation above that numerical threshold. Insets show the positive fraction, not active signaling. Ranks do not establish a common inter-age magnitude. Only eligible buffered ventricular CMs are shown as foreground, and all resolved cardiac source types compete within the same domain-specific solve. Light tissue and excluded buffers retain context. A 50-µm matrix range does not establish deposition or physical contact; a 200-µm secreted range does not establish diffusion or receptor activation. Scale bars, 500 µm.')

def prediction():
 p=read(OUT/'OOF_performance.tsv');g=read(OUT/'incremental_prediction_gains.tsv')
 p=p[p.view.eq('core')&p.learner.eq('RandomForest')&p.CM_state.eq('All')]
 fig=plt.figure(figsize=(11,8),layout='constrained')
 outer=fig.add_gridspec(3,1,height_ratios=[.13,1,1.25]);top=outer[1].subgridspec(1,3);bottom=outer[2].subgridspec(1,3)
 axs=np.array([[fig.add_subplot(top[j]) for j in range(3)],[fig.add_subplot(bottom[j]) for j in range(3)]])
 legend=fig.add_subplot(outer[0]);legend.axis('off')
 legend.legend(handles=[Line2D([],[],color=c,marker='o',label=('C + types' if n=='C_fine_context' else n),ms=4) for n,c in COLORS.items()],ncol=6,loc='center',frameon=False,fontsize=10)
 for ax,key,letter in [(axs[0,0],'maturation21','a'),(axs[0,1],'metabolic','b'),(axs[0,2],'cycling','c')]:
  for name,color in COLORS.items():
   z=p[p.feature_set.eq(name)&p.outcome.eq('core_Z__'+key)].set_index('stage').reindex(STAGES)
   ax.plot(range(4),z.OOF_R2,'o-',label=name,color=color,ms=4)
  ax.axhline(0,color='#bbb',lw=.7);ax.set_xticks(range(4),STAGES);ax.set_ylabel('Held-out-region R²');title(ax,letter,LABELS[key])
 q=g[g.primary_analysis];norm=TwoSlopeNorm(vmin=-max(.005,abs(q.delta_R2).max()),vcenter=0,vmax=max(.005,abs(q.delta_R2).max()))
 for ax,comparison,letter,label in [(axs[1,0],'D minus C_fine_context','d','Non-CM ligands · ΔR²'),(axs[1,1],'E minus C_fine_context','e','Ligands + transport · ΔR²')]:
  z=q[q.comparison.eq(comparison)];mat=z.pivot(index='outcome',columns='stage',values='delta_R2').reindex(index=['core_Z__'+k for k in LABELS],columns=STAGES).to_numpy()
  im=heat(ax,mat,['Overall 21','Metabolic','Sarcomere','Excitation','Cycling'],cmap='RdBu_r',norm=norm)
  for y,key in enumerate(LABELS):
   for x,stage in enumerate(STAGES):
    r=z[z.stage.eq(stage)&z.outcome.eq('core_Z__'+key)].iloc[0]
    ax.text(x,y,f'{r.delta_R2:+.3f}'.replace('+0.','+.').replace('-0.','-.')+('*' if r.primary_grid_BH_q<.05 else ''),ha='center',va='center',fontsize=9.5)
  title(ax,letter,label)
 fig.colorbar(ScalarMappable(norm=norm,cmap='RdBu_r'),ax=[axs[1,0],axs[1,1]],label='ΔR² over fine-context model',orientation='horizontal',shrink=.8,pad=.06,aspect=35)
 ax=axs[1,2];z=g[g.view.eq('core')&g.learner.eq('RandomForest')&g.CM_state.eq('All')&g.comparison.eq('E minus D')]
 for offset,(key,color) in enumerate(zip(['maturation21','metabolic','sarcomere','cycling'],['#AA3377','#228833','#4477AA','#EE7733'])):
  r=z[z.outcome.eq('core_Z__'+key)].set_index('stage').reindex(STAGES);xx=np.arange(4)+(offset-1.5)*.13
  short={'maturation21':'Overall 21','metabolic':'Metabolic','sarcomere':'Sarcomere','cycling':'Cycling'}
  ax.errorbar(xx,r.delta_R2,yerr=[r.delta_R2-r.CI_low,r.CI_high-r.delta_R2],fmt='o',color=color,capsize=2,ms=3,label=short[key])
 ax.axhline(0,color='#aaa',lw=.8);ax.set_xticks(range(4),STAGES);ax.set_ylabel('ΔR² over ligand-only')
 ax.legend(fontsize=9.5,frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.12))
 title(ax,'f','Transport beyond ligands',20)
 fig.suptitle('Testing additional information from neighboring molecular inputs',fontsize=14,weight='bold')
 finish(fig,'Figure_3_spatial_prediction','Testing additional predictive information in unseen regions',
  'a–c, Pooled out-of-fold R² across four source-disjoint spatial domains for the same CMs and fixed forest complexity. A: capture, anatomy, CM identity and broad context; B: add CM receptors; C: add CM-to-CM inputs; C + types (C_fine_context): add fine source counts and capture; D: add non-CM ligand exposure; E: add non-CM joint ligand–receptor transport. A negative R² means worse prediction than the held-out cohort mean. d,e, Incremental prediction over C_fine_context; asterisks mark BH q<0.05 in the supported primary 40-test family, including negative contrasts. f, Added transport over ligand-only features; 95% intervals use paired 800-µm block resampling conditional on fitted models, rather than independent hearts or causal effects. All components, learners, RNA views and folds are supplied in tables.')

def source_sensitivity():
 d=read(OUT/'source_feature_perturbations.tsv');d=d[d.view.eq('core')]
 sources=['Broad:'+g for g in BROAD]+['EC:'+g for g in EC]+['FB:'+g for g in FB]
 rows=[]
 for key,grp in d.groupby(['stage','perturbation','source_group','outcome']):
  rows.append(dict(zip(['stage','perturbation','source_group','outcome'],key),effect=np.average(grp.increase_MSE_fraction_of_test_variance,weights=grp.test_CMs)))
 avg=pd.DataFrame(rows);save(avg,'source_feature_sensitivity_summary.tsv')
 fig,axs=plt.subplots(2,3,figsize=(11,9),layout='constrained')
 vmax=max(.01,abs(avg.effect).quantile(.98));norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
 for row,null in enumerate(['within_CM_state','within_CM_state_region200']):
  for j,key in enumerate(['maturation21','metabolic','sarcomere']):
   ax=axs[row,j];q=avg[avg.perturbation.eq(null)&avg.outcome.eq('core_Z__'+key)]
   mat=q.pivot(index='source_group',columns='stage',values='effect').reindex(index=sources,columns=STAGES).to_numpy()
   labels=[g.replace('Broad:','').replace('EC:','EC: ').replace('FB:','FB: ').replace('Gfpt2','Gfpt2/Pi16') for g in sources]
   im=heat(ax,mat,labels,cmap='RdBu_r',norm=norm)
   ax.axhline(6.5,color='white',lw=2);ax.axhline(11.5,color='white',lw=2)
   short={'maturation21':'21-gene score','metabolic':'Metabolic','sarcomere':'Sarcomere'}
   title(ax,chr(97+row*3+j),short[key]+' · '+('regional' if row==0 else 'local · 200 µm'),20)
 fig.colorbar(ScalarMappable(norm=norm,cmap='RdBu_r'),ax=axs,label='Increase in prediction error / fold test variance',shrink=.6,pad=.02)
 fig.suptitle('How source information changes with spatial scale',fontsize=14,weight='bold')
 finish(fig,'Figure_4_source_feature_sensitivity','Source-feature sensitivity across spatial scales',
  'a–c, Held-out prediction change after jointly permuting the named source’s molecular features within CM state. d–f, Permutation additionally restricted to the same 200-µm spatial square. Values are test-CM-weighted averages of error changes normalized by each fold’s test variance; one fixed permutation is used per fold and scale. Positive values indicate that the fitted model uses these features, not that the source promotes maturation. Broad source perturbations include both broad and fine molecular descriptors; fine-source perturbations leave correlated pooled descriptors intact. Geometry, source counts and CM receptors remain fixed. Gray indicates no source-specific RNA-supported features, rather than absence of the cell type. This exploratory sensitivity analysis does not supply independent source-effect tests.')

def main():
 overview();spatial_maps();prediction();source_sensitivity()
 (OUT/'figure_legends.json').write_text(json.dumps(LEGENDS,indent=2,ensure_ascii=False)+'\n');export_inputs()
if __name__=='__main__':main()
