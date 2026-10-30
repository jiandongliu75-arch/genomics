"""Biological synthesis plates, with source-version-specific captions."""
from pathlib import Path
import os,json,shutil,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-biology-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Rectangle
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_biological_synthesis'
OLD=ROOT/'results/cardiac_subtype_niches';FIG=OUT/'figures'
STAGES=['P0','P7','P14','P21']
COLORS={'Postn-associated':'#E69F00','Gfpt2-associated':'#009E73','Col8a1-associated':'#0072B2','Lsamp-associated':'#CC79A7','Kirrel3-associated':'#8B5C38'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
                    'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def panel(ax,s,title):
    spec=ax.get_subplotspec()
    if spec.get_gridspec().ncols==2 and len(spec.colspan)==1:title=textwrap.fill(title,43)
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.08,1.05,s,transform=ax.transAxes,weight='bold',fontsize=12)
def export(fig,name):
    for ext in ['png','pdf','svg']:
        p=FIG/(name+'.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight',pad_inches=.09)
        if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)
def copy(old,new):
    for ext in ['png','pdf','svg']:shutil.copyfile(OLD/'figures'/(old+'.'+ext),FIG/(new+'.'+ext))

def fibroblasts():
    fig,axs=plt.subplots(3,2,figsize=(7.4,9.0),layout='constrained',gridspec_kw={'height_ratios':[1,1,.8]})
    for i,stage in enumerate(STAGES):
        ax=axs[i//2,i%2];d=read(OUT/(stage+'_FB_assignment_evidence.tsv.gz'));f=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',usecols=['x_centroid','y_centroid','atlas_identity','atlas_resolved'])
        z=f[f.atlas_resolved];ax.scatter(z.x_centroid,z.y_centroid,s=.12,c='#e2e2e2',lw=0,rasterized=True)
        selected=d[d.multigene_boundary_candidate.ne('Unresolved')]
        for fam,color in COLORS.items():
            z=selected[selected.multigene_boundary_candidate.eq(fam)];ax.scatter(z.x_centroid,z.y_centroid,s=1.5,c=color,lw=0,rasterized=True)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set_axis_off();panel(ax,chr(97+i),stage+' · boundary-supported FB candidates')
        x=f.x_centroid.min();y=f.y_centroid.max()+70;ax.plot([x,x+500],[y,y],c='black',lw=1.4);ax.text(x+250,y+55,'500 µm',ha='center',fontsize=6.5)
    ax=axs[2,0];d=read(OUT/'Postn_direct_Xenium_evidence.tsv')
    ax.plot(range(4),100*d.Postn_detection_fraction,'o-',c=COLORS['Postn-associated'],lw=1.7)
    ax.set(xticks=range(4),xticklabels=STAGES,ylim=(0,105),ylabel='FBs with ≥1 native Postn molecule (%)')
    panel(ax,'e','Postn detection declines after P7')
    ax=axs[2,1];d=read(OUT/'FB_assignment_coverage_comparison.tsv');s=d.groupby('stage')[['legacy_strict_family','boundary_supported_nuclear_optional','multigene_boundary_candidate']].sum().reindex(STAGES)
    den=d.groupby('stage').eligible_FB_profiles.first().reindex(STAGES)
    for col,label,color in [('legacy_strict_family','Nuclear-required legacy','#999999'),('boundary_supported_nuclear_optional','Nuclear optional + anchor','#4477AA'),('multigene_boundary_candidate','Multi-gene + boundary','#228833')]:
        ax.plot(range(4),100*s[col]/den,'o-',label=label,c=color,ms=3)
    ax.set(xticks=range(4),xticklabels=STAGES,ylim=(0,75),ylabel='Eligible FBs assigned (%)');ax.legend(frameon=False,fontsize=6.5)
    panel(ax,'f','Evidence layers change assignment coverage')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',c=c,label=k.replace('-associated','')+' FB') for k,c in COLORS.items()],loc='outside lower center',ncol=5,frameon=False,fontsize=7)
    export(fig,'Figure_2_fibroblast_development_and_assignment')

def forest(ax,old,new,stage,state,ex,outcome,letter,title):
    def choose(d,view,geom):
        return d[d.stage.eq(stage)&d.CM_state.eq(state)&d.source_subtype.eq(ex)&d.outcome.eq(outcome)&d.source_view.eq(view)&d.geometry.eq(geom)&d.adjustment.eq('within_region')]
    entries=[]
    for d,view,geom,label,c in [(old,'primary','proximity50','Legacy · 50 µm','#999999'),
                              (new,'boundary_supported_nuclear_optional','proximity50','Nuclear optional · 50 µm','#4477AA'),
                              (new,'multigene_boundary_candidate','proximity50','Multi-gene · 50 µm','#228833'),
                              (new,'multigene_boundary_candidate','proximity30','Multi-gene · 30 µm','#228833')]:
        z=choose(d,view,geom)
        if len(z):entries.append((z.iloc[0],label,c))
    for i,(r,label,c) in enumerate(entries):
        ax.errorbar(r.standardized_coefficient,i,xerr=[[r.standardized_coefficient-r.block_CI_low],[r.block_CI_high-r.standardized_coefficient]],fmt='o',c=c,ms=4,capsize=2)
    ax.axvline(0,c='#777',lw=.6);ax.set(yticks=range(len(entries)),yticklabels=[x[1] for x in entries],xlabel='Local association (standardized units)');ax.invert_yaxis();panel(ax,letter,title)

def vascular():
    fig=plt.figure(figsize=(7.4,9.1),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1.05,1,1]);a=fig.add_subplot(gs[0,:]);a.axis('off');panel(a,'a','Regional vascular expansion changes the setting in which CMs develop')
    cards=[('Late fetal–P7','Existing coronary plexus expands\npreferentially into inner myocardium', '#4477AA'),
           ('By P7','Inner vascular density approaches\nouter myocardial density', '#228833'),
           ('P7–P14 and later','Early angiogenic activity subsides;\nvascular integration remains a question', '#AA3377')]
    for i,(h,b,c) in enumerate(cards):
        x=.015+i*.335;a.add_patch(FancyBboxPatch((x,.22),.30,.59,boxstyle='round,pad=.01',fc=c+'10',ec=c))
        a.text(x+.15,.72,h,ha='center',weight='bold',fontsize=9);a.text(x+.15,.47,'\n'.join(textwrap.fill(line,28) for line in b.split('\n')),ha='center',va='center',fontsize=7.2)
        if i<2:a.annotate('',xy=(x+.33,.51),xytext=(x+.305,.51),arrowprops={'arrowstyle':'->','color':'#555'})
    a.text(.5,.08,'Published vascular findings (Lu et al., 2021); not newly measured perfusion in our sections',ha='center',fontsize=7,color='#555')
    old=read(OLD/'subtype_conditioned_maturation_coefficients.tsv');new=read(OUT/'nuclear_optional_FB_maturation_coefficients.tsv')
    forest(fig.add_subplot(gs[1,0]),old,new,'P7','p.CM','FB__Postn-associated','core_Z__metabolic','b','P7 cycling CMs · Postn FB / metabolic RNA')
    forest(fig.add_subplot(gs[1,1]),old,new,'P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere','c','P14 Ankrd1 CMs · endocardial EC / sarcomere RNA')
    ax=fig.add_subplot(gs[2,0]);d=read(ROOT/'results/cardiac_spatial_reconstruction/primary_matched_component_results.tsv')
    for k,(outcome,label,color) in enumerate([('core_Z__withdrawal','Lower cycling RNA','#AA3377'),('core_Z__metabolic','Metabolic-regulatory RNA','#4477AA')]):
        for i,s in enumerate(STAGES):
            r=d[d.stage.eq(s)&d.outcome.eq(outcome)].iloc[0];y=i+(k-.5)*.18
            ax.errorbar(r.mean_standardized_difference,y,xerr=[[r.mean_standardized_difference-r.block_CI_low],[r.block_CI_high-r.mean_standardized_difference]],fmt='o',ms=3,c=color,capsize=2,label=label if i==0 else None)
    ax.axvline(0,c='#777',lw=.6);ax.set(yticks=range(4),yticklabels=STAGES,xlabel='Apposed minus matched non-apposed CM');ax.invert_yaxis();ax.legend(frameon=False,fontsize=6.5);panel(ax,'d','Capillary apposition · different components over time')
    ax=fig.add_subplot(gs[2,1]);ax.axis('off');panel(ax,'e','Growth, access and functional integration are distinct')
    text=[('Growing vasculature','EC cycling / sprouting-associated RNA;\nnot measured oxygen delivery'),('Local vascular access','Capillary distance or boundary apposition;\nnot proof of a perfused vessel'),('CM developmental response','Metabolism, sarcomeres and cycling can differ;\nfunctional maturation requires direct measurements')]
    for y,(h,b) in zip([.82,.51,.20],text):
        ax.text(.02,y,h,weight='bold',fontsize=8);ax.text(.02,y-.075,b,fontsize=7.4,va='top')
    export(fig,'Figure_4_vascular_growth_and_CM_specialization')

def model():
    fig,axs=plt.subplots(3,1,figsize=(7.4,8.4),layout='constrained',gridspec_kw={'height_ratios':[1.0,1.1,1.05]})
    for ax in axs:ax.axis('off')
    ax=axs[0];panel(ax,'a','A developmental model: coordinated changes with different local timing')
    columns=[('Early growth · P0/P7','Postn-prominent remodeling matrix\nExpanding coronary microvasculature\nCycling and heterogeneous CM states'),('Transition · P7/P14','Inner myocardium gains vascular coverage\nFibroblast programs specialize\nCM metabolic and structural programs diverge'),('Later specialization · P14/P21','Col8a1 / Gfpt2 / Lsamp RNA states\nVascular and matrix contexts persist\nRegional CM maturation remains heterogeneous')]
    for i,(h,b) in enumerate(columns):
        x=.015+i*.335;ax.add_patch(FancyBboxPatch((x,.22),.30,.58,boxstyle='round,pad=.01',fc=['#fff5df','#eaf5f0','#eaf1fa'][i],ec='#9aa8b3'))
        ax.text(x+.15,.72,h,ha='center',weight='bold',fontsize=8.4);ax.text(x+.15,.44,'\n'.join(textwrap.fill(line,32) for line in b.split('\n')),ha='center',va='center',fontsize=7.2)
        if i<2:ax.annotate('',xy=(x+.33,.50),xytext=(x+.305,.50),arrowprops={'arrowstyle':'->'})
    ax.text(.5,.08,'Proposed synthesis; arrows denote developmental windows, not observed cell conversion or a causal chain',ha='center',fontsize=7)
    ax=axs[1];panel(ax,'b','Three competing explanations for regional CM development')
    items=[('Supply','Vascular expansion improves delivery,\nthen metabolic specialization follows.'),('Matrix / mechanics','Fibroblast specialization changes the\nsubstrate and load experienced by CMs.'),('Cellular integration','EC contact or local signals affect\nCM structure after vascular arrival.')]
    for i,(h,b) in enumerate(items):
        x=.015+i*.335;ax.add_patch(FancyBboxPatch((x,.20),.30,.59,boxstyle='round,pad=.01',fc='#f5f8fa',ec='#4477AA'))
        ax.text(x+.15,.70,h,ha='center',weight='bold',fontsize=9);ax.text(x+.15,.43,'\n'.join(textwrap.fill(line,34) for line in b.split('\n')),ha='center',va='center',fontsize=7.4)
    ax.text(.5,.07,'Ankrd1 and working CMs may respond differently within the same anatomical region',ha='center',fontsize=7.4)
    ax=axs[2];panel(ax,'c','Experiments that distinguish these explanations')
    rows=[('P3/P5/P10 in independent hearts','Co-measure perfusion, hypoxia, CM state and matrix organization; resolve temporal order.'),
          ('Match the starting CM identity','Separate perfusion/oxygen, EC contact and fibroblast-derived matrix; measure respiration, calcium and force.'),
          ('Trace proposed state continuity','Test early-to-late fibroblast ancestry; use transport to nominate comparisons, not to establish fate.')]
    for y,(h,b) in zip([.80,.51,.22],rows):
        ax.text(.015,y,h,weight='bold',fontsize=8);ax.text(.015,y-.075,textwrap.fill(b,112),fontsize=7.4,va='top')
    export(fig,'Figure_6_regional_developmental_model')

def audit():
    fig,axs=plt.subplots(2,2,figsize=(7.4,6.5),layout='constrained')
    d=read(OUT/'FB_assignment_coverage_comparison.tsv');q=d[d.family.eq('Gfpt2-associated')].set_index('stage').reindex(STAGES)
    ax=axs[0,0];x=np.arange(4)
    for k,(col,label,c) in enumerate([('public_comparator','Public labels','#999'),('whole_cell_supported','Whole cell','#4477AA'),('boundary_supported_nuclear_optional','Cell + boundary','#228833'),('legacy_strict_family','Nuclear-required','#AA3377')]):
        ax.bar(x+(k-1.5)*.18,q[col],width=.17,label=label,color=c)
    ax.set(xticks=x,xticklabels=STAGES,ylabel='Gfpt2-labelled / supported profiles');ax.legend(frameon=False,fontsize=6);panel(ax,'a','Nuclear confirmation is only one source of attrition')
    ax=axs[0,1];p=read(OUT/'FB_public_vs_whole_cell_proposals.tsv');p=p[p.stage.eq('P21')&p.public_label.eq('Gfpt2.Fib')].sort_values('profiles')
    ax.barh([x.replace('-associated','') for x in p.whole_cell_proposal],p.profiles,color='#4477AA');ax.set_xlabel('P21 public Gfpt2 FBs');panel(ax,'b','Whole-cell proposals often disagree before filtering')
    ax=axs[1,0];p=read(OUT/'Postn_direct_snRNA_evidence.tsv');p=p[p.view.eq('corrected')]
    for i,s in enumerate(STAGES):
        z=p[p.stage.eq(s)];ax.scatter(i+np.linspace(-.13,.13,len(z)),100*z.fraction_detected,c='#E69F00',s=18)
    ax.set(xticks=range(4),xticklabels=STAGES,ylim=(0,105),ylabel='Postn-detected nuclei (%)');panel(ax,'c','Residual Postn RNA across library × FB families')
    ax=axs[1,1];ax.axis('off');panel(ax,'d','Interpretation of the revised evidence layers')
    for y,h,b in [( .87,'Nuclear optional','Missing nuclear markers no longer reject a\nconcordant whole-cell / boundary-core call.'),(.58,'Multi-gene candidate','A missing naming marker alone no longer\nrejects a supported multi-gene signature.'),(.28,'Remaining uncertainty','Reference disagreement and incomplete panel\nseparation remain; labels are not forced.')]:
        ax.text(0,y,h,weight='bold',fontsize=8);ax.text(0,y-.08,b,va='top',fontsize=7.4)
    export(fig,'Figure_S1_fibroblast_assignment_audit')

def main():
    FIG.mkdir(exist_ok=True,parents=True)
    copy('Figure_1_CM_identity_and_maturation','Figure_1_CM_identity_and_maturation')
    fibroblasts()
    copy('Figure_4_fibroblast_specialization','Figure_3_fibroblast_matrix_specialization')
    vascular()
    copy('Figure_5_CM_identity_conditioned_transport','Figure_5_CM_identity_conditioned_transport')
    model();audit()
    print('BIOLOGICAL_FIGURES_COMPLETE',flush=True)
if __name__=='__main__':main()
