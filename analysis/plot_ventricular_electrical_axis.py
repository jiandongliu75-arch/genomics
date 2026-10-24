"""Plot literature-guided expression associations; no electrophysiological labels assigned."""
from pathlib import Path
import numpy as np,pandas as pd
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_additional_axes';d=pd.read_csv(OUT/'gene_profiles_by_stage_library_cluster.tsv',sep='\t')
genes=['Irx5','Kcnd2','Kcnip2','Smyd1'];sel=d[d.grouping.eq('age+microenvironment_cluster')&d.gene.isin(genes)]
wide=sel.pivot(index=['age','microenvironment_cluster','nuclei'],columns='gene',values=['mean_log1p_10k','fraction_corrected_ge1','fraction_raw_ge1']);wide.columns=['__'.join(c) for c in wide.columns];wide.to_csv(OUT/'electrical_axis_cluster_profiles.tsv',sep='\t')
fig,axes=plt.subplots(2,2,figsize=(16,13));plt.rcParams['pdf.fonttype']=42
for ax,stage in zip(axes.flat,['P0','P7','P14','P21']):
 pooled=sel[sel.age.eq(stage)].pivot(index='microenvironment_cluster',columns='gene',values='mean_log1p_10k');counts=sel[sel.age.eq(stage)].drop_duplicates('microenvironment_cluster').set_index('microenvironment_cluster').nuclei
 libs=d[d.grouping.eq('age+sample+microenvironment_cluster')&d.age.eq(stage)&d.gene.isin(genes)&d.nuclei.ge(30)]
 library_names=sorted(libs['sample'].unique());markers=['o','s'];palette={c:plt.get_cmap('tab10')(i) for i,c in enumerate(pooled.index)}
 for cluster,row in pooled.iterrows():
  ax.scatter(row.Irx5,row.Kcnd2,s=55,c=[palette[cluster]],edgecolors='black',linewidths=.5,zorder=4)
  offset={'P7_M1':(-45,10),'P7_M4':(10,-10),'P7_M3':(5,12),'P21_M1':(10,-13),'P21_M4':(5,12)}.get(cluster,(5,5))
  ax.annotate(cluster,(row.Irx5,row.Kcnd2),xytext=offset,textcoords='offset points',fontsize=9,color=palette[cluster])
  for j,lib in enumerate(library_names):
   values=libs[libs.microenvironment_cluster.eq(cluster)&libs['sample'].eq(lib)].set_index('gene').mean_log1p_10k
   if len(values):ax.scatter(values.Irx5,values.Kcnd2,s=30,marker=markers[j],facecolors='none',edgecolors=[palette[cluster]],alpha=.65,zorder=3)
 ax.set_title(stage,fontsize=17,weight='bold');ax.set_xlabel('Irx5 mean log1p RNA / 10,000');ax.set_ylabel('Kcnd2 mean log1p RNA / 10,000');ax.spines[['top','right']].set_visible(False);ax.margins(x=.2,y=.18)
 ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color='black',label='Pooled group mean')]+[Line2D([],[],marker=markers[j],linestyle='',markerfacecolor='none',color='gray',label=lib.replace('GSM','')) for j,lib in enumerate(library_names)],fontsize=8,loc='best',framealpha=.9)
fig.suptitle('Irx5–Kcnd2: candidate repolarization axis across ventricular RNA states',fontsize=19,weight='bold',y=.985)
fig.text(.5,.02,'Filled points: pooled expression groups; open points: original-library means (>=30 nuclei/group)\nOnly 5–9 groups per stage; individual-nucleus Irx5 detection is sparse. Association does not establish wall position, ion current or action-potential duration.\nStages use different axis ranges; raw-count and QC sensitivities are in the accompanying tables.',ha='center',fontsize=10)
fig.tight_layout(rect=(.02,.085,.98,.95));fig.savefig(OUT/'Irx5_Kcnd2_expression_axis.png',dpi=220);fig.savefig(OUT/'Irx5_Kcnd2_expression_axis.pdf');plt.close(fig)
print('Saved group/library expression plot and gene profiles; no new cell labels or clusters fitted.')
