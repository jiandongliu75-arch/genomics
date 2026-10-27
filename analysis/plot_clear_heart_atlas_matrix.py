"""Rearrange existing fibroblast results into four readable biological panels."""
from pathlib import Path
import os,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/clear-heart-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/heart_atlas_clear_language_revision'
STAGES=['P0','P7','P14','P21'];FAMILIES=['Postn-associated','Col8a1-associated','Gfpt2-associated','Lsamp-associated','Kirrel3-associated']
COL={'Col8a1-associated':'#0072B2','Gfpt2-associated':'#009E73'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.7,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'clear-heart-matrix-v1'})
def read(p):return pd.read_csv(ROOT/p,sep='\t')
def panel(ax,label,title):
    ax.set_title(textwrap.fill(title,43),loc='left',pad=10);ax.text(-.08,1.06,label,transform=ax.transAxes,fontsize=12,weight='bold')
def main():
    fig,axs=plt.subplots(2,2,figsize=(7.4,7.5),layout='constrained')
    d=read('results/cardiac_developmental_story/fibroblast_library_matrix_programs.tsv')
    d=d[d.view.eq('corrected')&d.level.eq('fibroblast_family')&~d.group.str.startswith('Unresolved')]
    order=[(s,f) for s in STAGES for f in FAMILIES if ((d.stage==s)&(d.group==f)).any()]
    programs=['Early developmental','Matrix remodeling','Fibrillar matrix','Basement matrix','Elastic matrix','Cycling']
    m=d.groupby(['stage','group','program']).mean_log1p_CPM.mean().unstack().reindex(index=pd.MultiIndex.from_tuples(order),columns=programs)
    z=(m-m.mean())/m.std(ddof=0);ax=axs[0,0];im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set(xticks=range(6),xticklabels=['Early growth','Remodeling','Fibrillar','Basement','Elastic','Cell cycle'],yticks=range(len(order)),yticklabels=[s+' '+f.split('-')[0] for s,f in order]);ax.tick_params(axis='x',rotation=50)
    fig.colorbar(im,ax=ax,shrink=.65,label='Relative expression within each program');panel(ax,'a','Fibroblast gene patterns change with age')
    d=read('results/cardiac_developmental_story/fibroblast_library_marker_profiles.tsv')
    names=['Col8a1/Svep1 fibroblasts','Col8a1/Hmcn2 fibroblasts','Lsamp/Slit3 fibroblasts','Lsamp/Col4a5 fibroblasts'];genes=['Svep1','Fap','Hmcn2','Col8a1','Lsamp','Slit3','Eln','Col4a5','Col4a6']
    d=d[d.view.eq('corrected')&d.stage.eq('P21')&d.level.eq('named_state')&d.group.isin(names)&d.gene.isin(genes)]
    v=d.groupby(['group','gene']).log1p_CPM.mean().unstack().reindex(index=names,columns=genes)
    det=d.groupby(['group','gene']).fraction_detected.mean().unstack().reindex(index=names,columns=genes);z=(v-v.mean())/v.std(ddof=0)
    ax=axs[0,1]
    for i in range(4):
        for j in range(9):ax.scatter(j,i,s=8+100*det.iloc[i,j],c=[z.iloc[i,j]],cmap='RdBu_r',vmin=-1.8,vmax=1.8,edgecolors='none')
    ax.set(xticks=range(9),xticklabels=genes,yticks=range(4),yticklabels=[n.replace(' fibroblasts','') for n in names],ylim=(3.6,-.6));ax.tick_params(axis='x',rotation=60)
    from matplotlib.cm import ScalarMappable
    fig.colorbar(ScalarMappable(norm=plt.Normalize(-1.8,1.8),cmap='RdBu_r'),ax=ax,shrink=.6,label='Relative gene expression')
    panel(ax,'b','P21 groups have different matrix-related features')
    d=read('results/fibroblast_names_and_neighborhoods/matched_family_context_sensitivity.tsv');d=d[d.stage.eq('P14')&d.radius_um.eq(30)]
    ax=axs[1,0]
    for j,key in enumerate(['arterial_within_50um','fraction_Smooth muscle_30um']):
        r=d[d.measure.eq(key)].iloc[0]
        for offset,value,family in [(-.16,r.target_matched_mean,'Col8a1-associated'),(.16,r.reference_matched_mean,'Gfpt2-associated')]:
            ax.bar(j+offset,100*value,width=.3,color=COL[family])
            ax.text(j+offset,100*value+.6,f'{100*value:.1f}',ha='center',fontsize=7)
    ax.set(xticks=[0,1],xticklabels=['Within 50 µm\nof an arterial structure','Smooth-muscle share\nof nearby cells'],ylabel='Earlier matched comparison (%)',ylim=(0,24))
    ax.legend(handles=[Rectangle((0,0),1,1,color=c,label=f.split('-')[0]) for f,c in COL.items()],frameon=False,fontsize=7)
    panel(ax,'c','Some Col8a1 fibroblasts lie closer to arteries')
    ax=axs[1,1];ax.axis('off');panel(ax,'d','RNA resemblance suggests comparisons across ages')
    points={'P7 Postn':(.13,.74),'P14 Col8a1':(.50,.74),'P21 Col8/Svep1':(.87,.74),'P7 Gfpt2':(.13,.29),'P14 Gfpt2':(.50,.29),'P21 Gfpt2':(.87,.29)}
    for name,(xx,yy) in points.items():ax.text(xx,yy,name.replace(' ','\n'),ha='center',va='center',fontsize=7.2,bbox={'boxstyle':'round,pad=.3','fc':'#f6f7f8','ec':'#999','lw':.6})
    d=read('results/cardiac_developmental_story/reviewed_fibroblast_candidate_state_continuities.tsv')
    pairs=[('P7 Postn','P14 Col8a1','P7_F2','P14_F2'),('P14 Col8a1','P21 Col8/Svep1','P14_F2','P21_F0'),('P7 Gfpt2','P14 Gfpt2','P7_F4','P14_F0'),('P14 Gfpt2','P21 Gfpt2','P14_F0','P21_F3')]
    for a,b,source,target in pairs:
        r=d[d.from_cluster.eq(source)&d.to_cluster.eq(target)].iloc[0].correlation
        xx,yy=points[a];uu,vv=points[b];ax.annotate('',xy=(uu-.10,vv),xytext=(xx+.10,yy),arrowprops={'arrowstyle':'->','ls':'--','color':'#555','lw':.9});ax.text((xx+uu)/2,yy+.15,f'r={r:.2f}',ha='center',fontsize=7)
    ax.text(.5,.05,'Comparisons to test; individual cells were not followed',ha='center',va='bottom',fontsize=6.8,transform=ax.transAxes)
    folder=OUT/'figures';folder.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']:
        p=folder/('Figure_2_fibroblast_matrix_patterns.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight',pad_inches=.09)
        if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text().splitlines())+'\n')
    plt.close(fig);print('CLEAR_MATRIX_FIGURE_CREATED',flush=True)
if __name__=='__main__':main()
