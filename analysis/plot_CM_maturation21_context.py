"""Integrated plates for directly measured CM maturation and neighboring cells."""
from pathlib import Path
import os, textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/maturation21-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/CM_maturation21_nonCM_context';FIG=OUT/'figures'
STAGES=['P0','P7','P14','P21']
EC={'Capillary-like':'#0072B2','Endocardial-like':'#E69F00','Arterial-like':'#D55E00','Venous-like':'#009E73','Lymphatic-like':'#CC79A7'}
FB={'Postn':'#CC6677','Gfpt2':'#228833','Col8a1':'#4477AA','Lsamp':'#AA3377','Kirrel3':'#CCBB44'}
OTHERS={'Pericyte':'#44AA99','Smooth muscle':'#EE7733','Immune':'#882255','Epicardial':'#332288'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,
    'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})

def read(p,**kw):
    p=Path(p);return pd.read_csv(p if p.is_absolute() else ROOT/p,sep='\t',**kw)
def save(fig,stem):
    FIG.mkdir(exist_ok=True,parents=True)
    for ext in ['png','pdf']:fig.savefig(FIG/(stem+'.'+ext),dpi=300,facecolor='white')
    plt.close(fig);print('MATURATION21_FIGURE',stem,flush=True)
def title(ax,letter,text,width=42):ax.set_title(letter+'  '+textwrap.fill(text,width),loc='left',weight='bold',pad=9)
def base(ax,f):
    ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values():s.set_visible(False)
    x0,x1=f.x_centroid.min(),f.x_centroid.max();y0,y1=f.y_centroid.min(),f.y_centroid.max()
    x=x0+.06*(x1-x0);y=y1-.04*(y1-y0)
    ax.plot([x,x+500],[y,y],color='#222',lw=1.8);ax.text(x+250,y-.02*(y1-y0),'500 µm',ha='center',fontsize=7)
def get(stage):
    f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
    c=read(OUT/stage/'CM_maturation21_scores.tsv.gz',index_col='cell_id');c=c[c.strict_ventricular_map_cohort]
    return f,c

def score_maps():
    data={s:get(s) for s in STAGES};allvalues=np.concatenate([c.core_Z__maturation21.to_numpy() for _,c in data.values()])
    norm=Normalize(*np.quantile(allvalues,[.02,.98]));fig=plt.figure(figsize=(14.6,8.5),layout='constrained')
    gs=fig.add_gridspec(2,5,width_ratios=[1,1,1,1,.055])
    for j,stage in enumerate(STAGES):
        f,c=data[stage];bg=f[f.atlas_resolved]
        if len(bg)>45000:bg=bg.sample(45000,random_state=SEED)
        for row in [0,1]:
            ax=fig.add_subplot(gs[row,j]);ax.scatter(bg.x_centroid,bg.y_centroid,s=.15,c='#e9e9e9',linewidths=0,rasterized=True)
            values=c.core_Z__maturation21 if row==0 else c.core_state_percentile
            ax.scatter(c.x_centroid,c.y_centroid,c=values,cmap='viridis',norm=norm if row==0 else Normalize(0,1),s=1.2,linewidths=0,rasterized=True)
            base(ax,f);title(ax,chr(97+row*4+j),stage+(' · 21-gene score' if row==0 else ' · rank within CM identity'),31)
            ax.text(.02,.95,f'n = {len(c):,}',transform=ax.transAxes,fontsize=8)
    fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),cax=fig.add_subplot(gs[0,4]),label='Frozen-reference 21-gene score')
    fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=fig.add_subplot(gs[1,4]),label='Within-stage, within-identity percentile')
    fig.suptitle('Cardiomyocyte maturation varies across the same cardiac section',fontsize=15,weight='bold')
    save(fig,'Figure_1_direct_21gene_CM_maturation_maps')

SEED=271027
def neighbor_maps():
    data={s:get(s) for s in STAGES};norm=Normalize(*np.quantile(np.concatenate([c.core_Z__maturation21 for _,c in data.values()]),[.02,.98]))
    fig=plt.figure(figsize=(14.6,12.4),layout='constrained');gs=fig.add_gridspec(6,4,height_ratios=[1,.105,1,.105,1,.105])
    palettes=[EC,FB,OTHERS]
    for row,palette in enumerate(palettes):
        for j,stage in enumerate(STAGES):
            f,c=data[stage];ax=fig.add_subplot(gs[2*row,j]);bg=f[f.atlas_resolved]
            if len(bg)>35000:bg=bg.sample(35000,random_state=SEED)
            ax.scatter(bg.x_centroid,bg.y_centroid,s=.12,c='#eeeeee',linewidths=0,rasterized=True)
            ax.scatter(c.x_centroid,c.y_centroid,c=c.core_Z__maturation21,cmap='Greys',norm=norm,s=.65,linewidths=0,alpha=.65,rasterized=True)
            if row==1:
                labels=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
                for name,color in palette.items():
                    ids=labels.index[labels.program_source_eligible&labels.recommended_candidate.eq(name+'-associated')];z=f.loc[ids]
                    ax.scatter(z.x_centroid,z.y_centroid,s=3,c=color,linewidths=0,rasterized=True)
            else:
                for name,color in palette.items():
                    take=f.program_source_eligible&((f.atlas_identity.eq('Endothelial')&f.EC_nuclear_type.eq(name)) if row==0 else f.atlas_identity.eq(name))
                    z=f[take];ax.scatter(z.x_centroid,z.y_centroid,s=3,c=color,linewidths=0,rasterized=True)
            base(ax,f);title(ax,chr(97+row*4+j),stage+' · '+['endothelial subtypes','fibroblast families','other non-CMs'][row],30)
        legend=fig.add_subplot(gs[2*row+1,:]);legend.axis('off')
        handles=[Line2D([],[],ls='',marker='o',color=color,label=(name.replace('-like','') if row==0 else 'Gfpt2/Pi16' if name=='Gfpt2' else name),ms=5) for name,color in palette.items()]
        legend.legend(handles=handles,ncol=len(handles),loc='center',frameon=False,fontsize=9)
    fig.suptitle('Non-CM identities alongside directly measured CM maturation\nDarker cardiomyocytes have higher 21-gene scores; source colors show identity',fontsize=14,weight='bold')
    save(fig,'Figure_2_nonCM_subtypes_over_CM_maturation')

def association_heatmaps():
    d=read(OUT/'nonCM_maturation_associations.tsv');d=d[d.primary_analysis]
    sources=['Broad__'+s for s in ['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']]+['EC__'+s for s in EC]+['FB__'+s for s in FB]
    labels=[s.replace('Broad__','').replace('EC__','EC: ').replace('FB__','FB: ').replace('Gfpt2','Gfpt2/Pi16') for s in sources]
    outcomes=[('core_Z__maturation21','Overall 21-gene score'),('core_Z__metabolic','Metabolic-regulatory component'),('core_Z__sarcomere','Sarcomere component'),('core_Z__cycling','Cell-cycle expression · separate outcome')]
    vmax=max(.1,np.quantile(abs(d.coefficient),.98));norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    fig,axs=plt.subplots(2,2,figsize=(12.5,10.7),layout='constrained')
    for i,(key,label) in enumerate(outcomes):
        ax=axs.flat[i];q=d[d.outcome.eq(key)]
        mat=q.pivot(index='source',columns='stage',values='coefficient').reindex(index=sources,columns=STAGES)
        img=ax.imshow(np.ma.masked_invalid(mat.to_numpy()),cmap='RdBu_r',norm=norm,aspect='auto');img.cmap.set_bad('#eeeeee')
        ax.set_xticks(range(4),STAGES);ax.set_yticks(range(len(sources)),labels,fontsize=8);title(ax,chr(97+i),label,44)
        for yi,source in enumerate(sources):
            for xi,stage in enumerate(STAGES):
                z=q[q.source.eq(source)&q.stage.eq(stage)]
                if len(z)==0:ax.text(xi,yi,'—',ha='center',va='center',fontsize=8,color='#888');continue
                r=z.iloc[0];sig='*' if r.primary_grid_BH_q<.05 else ''
                ax.text(xi,yi,f'{r.coefficient:+.2f}{sig}',ha='center',va='center',fontsize=7,color='white' if abs(r.coefficient)>.65*vmax else '#222')
        ax.axhline(5.5,color='white',lw=2);ax.axhline(10.5,color='white',lw=2)
    fig.colorbar(ScalarMappable(norm=norm,cmap='RdBu_r'),ax=axs,label='Within-region standardized association',shrink=.7,pad=.02)
    fig.suptitle('Neighbor relationships differ across CM developmental components\n* Global primary BH q < 0.05; gray = unavailable or insufficient support',fontsize=14,weight='bold')
    save(fig,'Figure_3_nonCM_relationships_with_CM_components')

def forest(ax,d,x,y,lo,hi,label,color,offset=0):
    q=d.set_index(x).reindex(STAGES);xx=np.arange(4)+offset
    ax.errorbar(xx,q[y],yerr=[q[y]-q[lo],q[hi]-q[y]],fmt='o-',label=label,color=color,capsize=3,ms=4)
    ax.axhline(0,color='#999',lw=.8);ax.set_xticks(range(4),STAGES)

def capillary_audit():
    all_d=read(OUT/'matched_capillary_maturation_results.tsv');d=all_d[all_d.cohort.eq('strict_ventricular_pairs')]
    fig,axs=plt.subplots(2,3,figsize=(14,8.3),layout='constrained')
    ax=axs[0,0]
    for view,col,off in [('core','#0072B2',-.05),('native','#D55E00',.05)]:
        forest(ax,d[d.outcome.eq(view+'_Z__maturation21')],'stage','standardized_difference','CI_low','CI_high',view.capitalize()+' CM RNA',col,off)
    qc=all_d[all_d.cohort.eq('QC_capillary_strict_pairs')]
    forest(ax,qc[qc.outcome.eq('core_Z__maturation21')],'stage','standardized_difference','CI_low','CI_high','QC capillaries · core RNA','#009E73',.12)
    title(ax,'a','Full 21-gene score near capillaries');ax.set_ylabel('Apposed minus matched comparison (SD)');ax.legend(frameon=False,fontsize=8)
    for ax,key,label,letter in [(axs[0,1],'metabolic','Metabolic-regulatory expression','b'),(axs[0,2],'cycling','Cell-cycle expression','c')]:
        forest(ax,d[d.outcome.eq('core_Z__'+key)],'stage','standardized_difference','CI_low','CI_high','Core CM RNA','#0072B2')
        forest(ax,qc[qc.outcome.eq('core_Z__'+key)],'stage','standardized_difference','CI_low','CI_high','QC capillaries · core RNA','#009E73',.10)
        ax.legend(frameon=False,fontsize=7)
        title(ax,letter,label);ax.set_ylabel('Apposed minus comparison (component SD)')
    ax=axs[1,0]
    selected=[('core_Z__maturation21','Full 21'),('core_Z__without_Cebpb_Rxrg_Lpin1','Without Cebpb/Rxrg/Lpin1'),('core_Z__without_metabolic8','Without metabolic 8'),('core_Z__without_sarcomere3','Without sarcomere 3'),('core_Z__without_excitation5','Without excitation 5')]
    q=d[d.stage.eq('P7')].set_index('outcome')
    for j,(key,label) in enumerate(selected):
        r=q.loc[key];ax.errorbar(r.standardized_difference,j,xerr=[[r.standardized_difference-r.CI_low],[r.CI_high-r.standardized_difference]],fmt='o',color='#0072B2',capsize=3)
    ax.set_yticks(range(len(selected)),[x[1] for x in selected],fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#999',lw=.8)
    title(ax,'d','P7 · fixed-reference gene omissions');ax.set_xlabel('Difference in full-score SD units')
    ax=axs[1,1];q=all_d[all_d.cohort.eq('strict_ventricular_depth50_pairs')]
    for key,label,color,off in [('core_Z__maturation21','Original depth · same CMs','#0072B2',-.05),('core_Z__depth50_mean3','50 molecules · three draws','#009E73',.05)]:
        forest(ax,q[q.outcome.eq(key)],'stage','standardized_difference','CI_low','CI_high',label,color,off)
    title(ax,'e','Depth sensitivity uses identical pairs');ax.set_ylabel('Difference in full-score SD units');ax.legend(frameon=False,fontsize=8)
    ax=axs[1,2]
    for key,label,color,off in [('core_Z__sarcomere','Sarcomere','#CC79A7',-.05),('core_Z__excitation_contraction','Excitation–contraction','#E69F00',.05)]:
        forest(ax,d[d.outcome.eq(key)],'stage','standardized_difference','CI_low','CI_high',label,color,off)
    ax.set_ylabel('Apposed minus comparison (component SD)');ax.legend(frameon=False,fontsize=8)
    title(ax,'f','Other CM components show different patterns')
    save(fig,'Figure_4_capillary_score_and_robustness_tests')

def fibroblast_development():
    d=read(OUT/'fibroblast_library_matrix_scores.tsv');fig,axs=plt.subplots(2,3,figsize=(14,8.5),layout='constrained')
    for ax,program,label,letter in [(axs[0,0],'Fibronectin','Fibronectin expression','a'),(axs[0,1],'Fibrillar_collagen','Fibrillar-collagen expression','b'),(axs[0,2],'Basement_membrane','Basement-membrane expression','c')]:
        for family,color in FB.items():
            q=d[d.family.eq(family+'-associated')&d.program.eq(program)]
            q=q.assign(day=q.age.str[1:].astype(int));m=q.groupby('day').score_mean_log1p_CPM.mean()
            ax.plot(m.index,m.values,'o-',color=color,label=family,ms=4)
            ax.scatter(q.day,q.score_mean_log1p_CPM,s=14,color=color,alpha=.7)
        ax.set_xticks([0,7,14,21],STAGES);ax.set_ylabel('Library pseudobulk log1p CPM');title(ax,letter,label)
    axs[0,0].legend(frameon=False,fontsize=8,ncol=2)
    ax=axs[1,0];q=read(OUT/'fibroblast_matrix_gene_stage_changes.tsv')
    contrasts=[('Gfpt2-associated','P7','P14','Gfpt2/Pi16 · P7→P14'),('Col8a1-associated','P14','P21','Col8a1 · P14→P21')]
    genes=['Fn1','Col1a1','Col1a2','Col3a1','Col8a1']
    for j,(family,a,b,label) in enumerate(contrasts):
        z=q[q.family.eq(family)&q.from_stage.eq(a)&q.to_stage.eq(b)].set_index('gene').reindex(genes)
        ax.bar(np.arange(len(genes))+(j-.5)*.34,z.corrected_log2CPM_difference,width=.34,label=label,color=['#228833','#4477AA'][j])
    ax.axhline(0,color='#999',lw=.8);ax.set_xticks(range(len(genes)),genes,rotation=35);ax.set_ylabel('Later minus earlier log2 CPM');ax.legend(frameon=False,fontsize=7)
    title(ax,'d','Matrix changes are not a uniform collagen increase')
    ax=axs[1,1];ax.axis('off');title(ax,'e','Candidate developmental continuity · not ancestry')
    paths=[('Postn','Col8a1','Col8a1',[.410,.697]),('Gfpt2','Gfpt2','Gfpt2',[.554,.689])]
    for row,(a,b,c,corr) in enumerate(paths):
        y=.70-row*.38
        for j,(name,day) in enumerate(zip([a,b,c],['P7','P14','P21'])):
            x=.10+j*.39;ax.text(x,y,name,ha='center',va='center',fontsize=10,color=FB[name],weight='bold',transform=ax.transAxes)
            ax.text(x,y+.12,day,ha='center',fontsize=9,transform=ax.transAxes)
            if j<2:
                ax.annotate('',xy=(x+.29,y),xytext=(x+.10,y),xycoords='axes fraction',arrowprops={'arrowstyle':'->','linestyle':'--','color':'#777'})
                ax.text(x+.195,y-.09,f'r = {corr[j]:.3f}',ha='center',fontsize=8,transform=ax.transAxes)
    ax.text(.5,.02,'Selected marker-withheld expression matches;\nno measured fate or conversion rate',ha='center',fontsize=9,transform=ax.transAxes)
    ax=axs[1,2];q=read(OUT/'fibroblast_spatial_geometry_snapshots.tsv')
    for name,color in FB.items():
        z=q[q.family.eq(name)&q.reference_assessable&q.within_ventricular_selector_profiles.ge(30)].copy()
        ax.plot(z.stage.str[1:].astype(int),z.median_relative_wall_geometry,'o-',label=name,color=color)
    ax.set_xticks([0,7,14,21],STAGES);ax.set_ylim(0,1);ax.set_ylabel('Median operational wall-geometry coordinate')
    title(ax,'f','Spatial snapshots do not track migration');ax.text(.02,.02,'Geometry proxy, not validated transmural depth',transform=ax.transAxes,fontsize=7)
    save(fig,'Figure_6_fibroblast_matrix_development')

def score_audit():
    fig,axs=plt.subplots(2,3,figsize=(14.5,9),layout='constrained')
    d=read(OUT/'same_CM_depth_and_gene_omission_stage_means.tsv');d=d[d.CM_state.eq('All')&d.view.eq('core')]
    ax=axs[0,0]
    for key,label,color in [('maturation21','Original depth','#0072B2'),('depth50_mean3','50 molecules · three draws','#009E73')]:
        z=d[d.score.eq(key)].set_index('stage').reindex(STAGES);ax.plot(range(4),z['mean']-z['mean'].iloc[0],'o-',label=label,color=color)
    ax.set_xticks(range(4),STAGES);ax.set_ylabel('Mean score relative to P0');ax.legend(frameon=False,fontsize=8)
    title(ax,'a','Age trend persists but weakens at equal depth')
    ax=axs[0,1]
    for key,label,color in [('maturation21','Full 21','#0072B2'),('without_Cebpb_Rxrg_Lpin1','Without Cebpb/Rxrg/Lpin1','#CC79A7')]:
        z=d[d.score.eq(key)].set_index('stage').reindex(STAGES);ax.plot(range(4),z['mean']-z['mean'].iloc[0],'o-',label=label,color=color)
    ax.set_xticks(range(4),STAGES);ax.set_ylabel('Mean index relative to its P0 mean');ax.legend(frameon=False,fontsize=8)
    title(ax,'b','Joint omission changes the temporal signal')
    ax=axs[0,2];q=read(OUT/'gene_contributions_to_stage_score_change.tsv');q=q[q.view.eq('core')].sort_values('contribution_to_21gene_score_change')
    ax.barh(range(len(q)),q.contribution_to_21gene_score_change,color=np.where(q.prespecified_three_gene_omission,'#CC79A7','#4477AA'))
    ax.set_yticks(range(len(q)),q.gene,fontsize=7);ax.axvline(0,color='#999',lw=.8);ax.set_xlabel('Contribution to P21–P0 score difference')
    title(ax,'c','Temporal contributions are unequal')
    ax=axs[1,0]
    for stage,color in zip(STAGES,['#4477AA','#EE7733','#228833','#CC3311']):
        _,c=get(stage);z=c.sample(min(1500,len(c)),random_state=SEED)
        ax.scatter(np.log10(z.core_panel_depth),z.core_Z__maturation21,s=2,alpha=.25,color=color,label=stage,rasterized=True)
    ax.set_xlabel('log₁₀ core measured-gene molecules');ax.set_ylabel('21-gene score');ax.legend(frameon=False,fontsize=8,markerscale=3)
    title(ax,'d','RNA depth remains an important consideration')
    ax=axs[1,1];q=read(OUT/'leave_one_gene_out_correlations.tsv');q=q[q.scope.eq('strict_ventricular')]
    for view,color,offset in [('core','#0072B2',-.04),('native','#D55E00',.04)]:
        z=q[q.view.eq(view)].groupby('stage').Pearson_with_full.agg(['mean','min','max']).reindex(STAGES)
        ax.errorbar(np.arange(4)+offset,z['mean'],yerr=[z['mean']-z['min'],z['max']-z['mean']],fmt='o-',color=color,capsize=3,label=view.capitalize()+' CM RNA')
    ax.set_xticks(range(4),STAGES);ax.set_ylim(.95,1.002);ax.set_ylabel('Pearson r with full score');ax.legend(frameon=False,fontsize=8)
    title(ax,'e','Single-gene omission preserves score ordering')
    ax=axs[1,2];q=read(OUT/'incremental_spatial_prediction.tsv');q=q[q.outcome.eq('core_Z__maturation21')]
    baseline=q[q.model.eq('anatomy_technical_broad')].set_index('stage').reindex(STAGES).spatially_heldout_R2
    for j,(model,label,color) in enumerate([('plus_EC_FB_identities','Add EC/FB identities','#0072B2'),('plus_FB_measured_programs','Add identities + FB programs','#CC79A7')]):
        z=q[q.model.eq(model)].set_index('stage').reindex(STAGES);ax.bar(np.arange(4)+(j-.5)*.32,z.spatially_heldout_R2-baseline,width=.32,label=label,color=color)
    ax.axhline(0,color='#999',lw=.8);ax.set_xticks(range(4),STAGES);ax.set_ylabel('Change in spatially held-out R²');ax.legend(frameon=False,fontsize=7)
    title(ax,'f','Finer labels add no overall-score prediction here')
    save(fig,'Figure_5_maturation_score_depth_and_gene_audit')

if __name__=='__main__':
    score_maps();neighbor_maps();association_heatmaps();capillary_audit();score_audit();fibroblast_development()
