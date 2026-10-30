"""Integrated plates for alternate FB assignments and continuous programs."""
from pathlib import Path
import os,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/fb-multimethod-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fibroblast_multimethod_identity import ROOT,OUT,STAGES,FAMILIES,METHODS,read,save

COLORS={'Postn-associated':'#E69F00','Gfpt2-associated':'#009E73','Col8a1-associated':'#0072B2','Lsamp-associated':'#CC79A7','Kirrel3-associated':'#8B5C38','Unresolved':'#cfcfcf'}
MNAME={'balanced_logistic':'Balanced logistic','rank_centroid':'Rank centroid','negative_binomial':'Count model'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.5,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,
 'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'fb-multimethod-v1'})
def panel(ax,label,title):
    ax.set_title(textwrap.fill(title,45),loc='left',pad=9);ax.text(-.08,1.06,label,transform=ax.transAxes,weight='bold',fontsize=12)
def export(fig,stem):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    for ext in ['png','pdf','svg']:
        p=folder/(stem+'.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight',pad_inches=.09)
        if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)
def annotations(stage):return read(OUT/(stage+'_multimethod_FB_annotations.tsv.gz'),index_col='cell_id')
def mapbase(ax,stage,label,title):
    f=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',index_col=0)
    z=f[f.atlas_resolved];ax.scatter(z.x_centroid,z.y_centroid,s=.10,c='#e5e5e5',lw=0,rasterized=True)
    x=f.x_centroid.min()+100;y=f.y_centroid.max()+30;ax.plot([x,x+500],[y,y],c='#222',lw=1.1);ax.text(x+250,y+55,'500 µm',ha='center',va='top',fontsize=6)
    ax.set_aspect('equal');ax.invert_yaxis();ax.set_axis_off();panel(ax,label,title)
def scatter_labels(ax,d,labels):
    for family,color in COLORS.items():
        if family=='Unresolved':continue
        z=d[labels.eq(family)];ax.scatter(z.x_centroid,z.y_centroid,s=1.5,c=color,lw=0,rasterized=True)

def figure1():
    d=read(OUT/'RNA_opposite_library_model_benchmark.tsv');fig,axs=plt.subplots(2,2,figsize=(7.4,7.0),layout='constrained')
    z=d[d.stage.eq('P21')&d.variant.eq('shared_panel')&d.condition.eq('observed')&d.family.eq('Gfpt2-associated')]
    ax=axs[0,0]
    for i,key in enumerate(['precision','recall']):
        for j,method in enumerate(METHODS):
            values=z[z.method.eq(method)][key];ax.scatter(np.full(len(values),j+(i-.5)*.18),values,s=16,c=['#4477AA','#EE7733'][i],marker=['o','s'][i],label=key.title() if j==0 else None)
    ax.set(xticks=range(3),xticklabels=[MNAME[m] for m in METHODS],ylim=(0,1.05),ylabel='Reference recovery');ax.tick_params(axis='x',rotation=20);ax.legend(frameon=False);panel(ax,'a','P21 Gfpt2 recovery: sensitivity can cost specificity')
    ax=axs[0,1];z=d[d.stage.eq('P21')&d.variant.eq('shared_panel')&d.condition.eq('observed')]
    mat=z.pivot_table(index='method',columns='family',values='F1').reindex(index=METHODS,columns=FAMILIES[1:])
    im=ax.imshow(mat,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto');ax.set(xticks=range(4),xticklabels=[x.replace('-associated','') for x in mat.columns],yticks=range(3),yticklabels=[MNAME[m] for m in METHODS]);ax.tick_params(axis='x',rotation=30)
    for i in range(3):
        for j in range(4):ax.text(j,i,f'{mat.iloc[i,j]:.2f}',ha='center',va='center',color='white' if mat.iloc[i,j]>.6 else 'black',fontsize=7)
    fig.colorbar(im,ax=ax,shrink=.6,label='Mean opposite-library F1');panel(ax,'b','Recoverability differs by reference family')
    ax=axs[1,0]
    for j,variant in enumerate(['shared_panel','naming_markers_withheld','late_pooled_reference']):
        z=d[d.stage.eq('P21')&d.variant.eq(variant)&d.condition.eq('observed')&d.family.eq('Gfpt2-associated')]
        for i,method in enumerate(METHODS):
            v=z[z.method.eq(method)].F1;ax.scatter(np.full(len(v),i+(j-1)*.15),v,s=17,c=['#4477AA','#228833','#AA3377'][j],label=['Shared panel','Naming genes withheld','P14/P21 pooled'][j] if i==0 else None)
    ax.set(xticks=range(3),xticklabels=[MNAME[m] for m in METHODS],ylim=(0,1),ylabel='Gfpt2 reference F1');ax.tick_params(axis='x',rotation=20);ax.legend(frameon=False,fontsize=6.5);panel(ax,'c','Naming-marker and reference-age sensitivity')
    ax=axs[1,1]
    for j,condition in enumerate(['observed','Poisson_thinned_to60','15percent_gene_dropout']):
        z=d[d.stage.eq('P21')&d.variant.eq('shared_panel')&d.condition.eq(condition)&d.family.eq('Gfpt2-associated')]
        for i,method in enumerate(METHODS):
            values=z[z.method.eq(method)].F1
            ax.scatter(np.full(len(values),i+(j-1)*.15),values,s=17,c=['#4477AA','#EE7733','#AA3377'][j],label=['Observed RNA','Thinned to ~60','15% gene dropout'][j] if i==0 else None)
    ax.set(xticks=range(3),xticklabels=[MNAME[m] for m in METHODS],ylim=(0,1),ylabel='Gfpt2 reference F1');ax.tick_params(axis='x',rotation=20);ax.legend(frameon=False,fontsize=6.5);panel(ax,'d','Information loss limits even reference recovery')
    export(fig,'Figure_1_reference_recovery_and_sensitivity')

def figure2():
    stage='P21';d=annotations(stage);d=d[d.program_source_eligible];modules=read(OUT/(stage+'_published_module_scores.tsv.gz'),index_col='cell_id').loc[d.index]
    fig,axs=plt.subplots(3,2,figsize=(7.4,9.2),layout='constrained')
    for ax,label,title in zip(axs.flat,'abcdef',['Public named FB labels','Multi-method cell-supported candidates','Cell/core/naming-withheld candidates','Gfpt2-associated program intensity','Gfpt2 / Col8a1 program overlap','Public Gfpt2 labels: program evidence']):mapbase(ax,stage,label,'P21 · '+title)
    aliases={f.split('-')[0]+'.Fib':f for f in FAMILIES if f!='Kirrel3-associated'}
    scatter_labels(axs[0,0],d,d.author_atlas_cell_state.map(aliases).fillna('Unresolved'))
    scatter_labels(axs[0,1],d,d.recommended_candidate);scatter_labels(axs[1,0],d,d.crosschecked_candidate)
    im=axs[1,1].scatter(d.x_centroid,d.y_centroid,c=modules['Gfpt2.Fib__rank_excess'],s=1.3,cmap='RdBu_r',vmin=-.3,vmax=.3,lw=0,rasterized=True)
    fig.colorbar(im,ax=axs[1,1],shrink=.55,label='Rank score − matched controls')
    gf=modules['Gfpt2.Fib__enrichment_flag'];co=modules['Col8a1.Fib__enrichment_flag']
    for ax,mask in [(axs[2,0],pd.Series(True,index=d.index)),(axs[2,1],d.author_atlas_cell_state.eq('Gfpt2.Fib'))]:
        for label,flag,color in [('Gfpt2 program',gf&~co,'#009E73'),('Col8a1 program',co&~gf,'#0072B2'),('Both',co&gf,'#AA3377')]:
            z=d[mask&flag];ax.scatter(z.x_centroid,z.y_centroid,s=2,c=color,lw=0,rasterized=True)
    handles=[Line2D([],[],marker='o',ls='',c=c,label=f.replace('-associated','')) for f,c in COLORS.items() if f!='Unresolved']
    handles.append(Line2D([],[],marker='o',ls='',c='#AA3377',label='Both programs'))
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=7)
    export(fig,'Figure_2_P21_identity_and_program_maps')

def figure3():
    fig,axs=plt.subplots(3,2,figsize=(7.4,9.0),layout='constrained',gridspec_kw={'height_ratios':[1,1,.75]});rows=[]
    for i,stage in enumerate(STAGES):
        d=annotations(stage);take=d.program_source_eligible;ax=axs[i//2,i%2];mapbase(ax,stage,chr(97+i),stage+' · multi-method candidate identities');scatter_labels(ax,d[take],d.loc[take,'recommended_candidate'])
        for family in FAMILIES:
            selected=take&d.recommended_candidate.eq(family)
            rows.append({'stage':stage,'family':family,'eligible_FB_profiles':int(take.sum()),'candidate_profiles':int(selected.sum()),
                'cell_core_supported':int((take&d.boundary_supported_candidate.eq(family)).sum()),'crosschecked_profiles':int((take&d.crosschecked_candidate.eq(family)).sum()),
                'nuclear_corroborated':int((selected&d.nuclear_supports_cell_candidate).sum()),'supported_nuclear_alternative':int((selected&d.nuclear_supported_alternative).sum())})
    summary=pd.DataFrame(rows);save(summary,'final_candidate_coverage.tsv')
    ax=axs[2,0];mat=summary.pivot(index='stage',columns='family',values='candidate_profiles').reindex(index=STAGES,columns=FAMILIES).fillna(0)
    denom=summary.groupby('stage').eligible_FB_profiles.first().reindex(STAGES);bottom=np.zeros(4)
    for family in FAMILIES:
        v=100*mat[family]/denom;ax.bar(range(4),v,bottom=bottom,color=COLORS[family],label=family.replace('-associated',''));bottom+=v
    ax.set(xticks=range(4),xticklabels=STAGES,ylim=(0,100),ylabel='Eligible FBs (%)');panel(ax,'e','Candidate coverage, not unbiased abundance')
    ax=axs[2,1];p=read(OUT/'core_vs_depth_matched_assignment_stability.tsv');p=p[p.stage.eq('P21')];families=[f for f in FAMILIES if f!='Postn-associated']
    for j,(col,label,color) in enumerate([('core_majority_agreement','Boundary core','#4477AA'),('depth_matched_majority_agreement','Depth-matched thinning','#EE7733')]):
        z=p.groupby('family')[col].mean();ax.bar(np.arange(4)+(j-.5)*.3,[z.get(f,np.nan) for f in families],width=.29,color=color,label=label)
    ax.set(xticks=range(4),xticklabels=[x.replace('-associated','') for x in families],ylim=(0,1),ylabel='Same majority label');ax.tick_params(axis='x',rotation=25);ax.legend(frameon=False,fontsize=6.5);panel(ax,'f','P21 boundary sensitivity versus RNA loss')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',c=COLORS[f],label=f.replace('-associated','')) for f in FAMILIES],loc='outside lower center',ncol=5,frameon=False,fontsize=7)
    export(fig,'Figure_3_stage_maps_and_boundary_stability')

def figure4():
    fig,axs=plt.subplots(2,2,figsize=(7.4,7.2),layout='constrained')
    d=read(OUT/'heldout_naming_marker_profiles.tsv');d=d[d.stage.eq('P21')&d.view.eq('refined')&d.grouping.eq('withheld_majority')]
    families=['Gfpt2-associated','Col8a1-associated','Lsamp-associated','Kirrel3-associated','Unresolved'];genes=['Gfpt2','Pi16','Col8a1','Lsamp']
    ax=axs[0,0];mat=d.pivot(index='identity',columns='gene',values='mean_log1p_per10000_shared_panel').reindex(index=families,columns=genes)
    det=d.pivot(index='identity',columns='gene',values='detection_fraction').reindex(index=families,columns=genes)
    im=ax.imshow(mat,cmap='YlGnBu',aspect='auto',vmin=0,vmax=5)
    for i in range(5):
        for j in range(4):ax.text(j,i,f'{100*det.iloc[i,j]:.0f}%',ha='center',va='center',fontsize=7,c='white' if mat.iloc[i,j]>3 else '#222')
    ax.set(xticks=range(4),xticklabels=genes,yticks=range(5),yticklabels=[f.replace('-associated','') for f in families]);fig.colorbar(im,ax=ax,shrink=.6,label='Mean log1p shared-panel CP10k');panel(ax,'a','Genes withheld from classifier remain informative')
    ax=axs[0,1];d=read(OUT/'published_Gfpt2_Col8a1_program_overlap.tsv');d=d[d.stage.eq('P21')&d.scope.eq('public_Gfpt2')]
    labels=['Gfpt2_only','Col8a1_only','both','neither'];colors=['#009E73','#0072B2','#AA3377','#ccc'];s=d.set_index('program_combination')
    ax.bar(range(4),[100*s.loc[x,'profiles']/s.loc[x,'scope_profiles'] for x in labels],color=colors)
    ax.set(xticks=range(4),xticklabels=['Gfpt2 only','Col8a1 only','Both','Neither'],ylabel='Public P21 Gfpt2 profiles (%)');ax.tick_params(axis='x',rotation=30);panel(ax,'b','The public Gfpt2 population contains varied programs')
    ax=axs[1,0];d=read(OUT/'published_FB_module_coverage.tsv');d=d[d.author_family.isin(['Postn.Fib','Gfpt2.Fib','Col8a1.Fib','Lsamp.Fib','valve.Fib'])]
    ax.barh(d.author_family,d.measured_upregulated_markers,color=['#4477AA' if b else '#bbb' for b in d.module_scored]);ax.set_xlabel('Published markers in shared measured panel');panel(ax,'c','Published signatures are small and partly overlapping')
    ax=axs[1,1];ax.axis('off');panel(ax,'d','Identity and program evidence answer different questions')
    rows=[('Subtype candidate','Several reference models agree;\nRNA screens and boundary checks qualify it.'),('Program-positive profile','A measured multi-gene program is enriched;\nit need not form a separate subtype.'),('Developmental fate','Requires ancestry evidence. Shared programs\nor cross-stage resemblance do not establish it.')]
    for y,(h,b) in zip([.83,.53,.23],rows):ax.text(.01,y,h,weight='bold',fontsize=8);ax.text(.01,y-.075,b,va='top',fontsize=7.2)
    export(fig,'Figure_4_heldout_genes_and_published_programs')

def main():
    for fn in [figure1,figure2,figure3,figure4]:print('FB_MULTIMETHOD_FIGURE',fn.__name__,flush=True);fn()
if __name__=='__main__':main()
