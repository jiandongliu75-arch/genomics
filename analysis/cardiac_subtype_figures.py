"""Subtype-first integrated figure plates; broad clustering is background only."""
from pathlib import Path
import os,json,shutil,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-subtype-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch,Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import spearmanr
from cardiac_subtype_niches import ROOT,OUT,STAGES,CM,FB,EC,read,save,source,SOURCES
from cardiac_story_figures import panel,mapbase,AGE

CP={'v.CM':'#4477AA','Ankrd1.CM':'#EE7733','Slit2.CM':'#228833','Ddc.CM':'#AA3377','p.CM':'#CC3311'}
CN={'v.CM':'Working','Ankrd1.CM':'Ankrd1','Slit2.CM':'Slit2','Ddc.CM':'Ddc','p.CM':'Cycling'}
EX=['FB__Postn-associated','FB__Gfpt2-associated','FB__Col8a1-associated','EC__Capillary-like','EC__Endocardial-like','EC__Arterial-like','EC__Venous-like','EC__Lymphatic-like']
XL=['Postn FB','Gfpt2 FB','Col8a1 FB','Capillary EC','Endocardial EC','Arterial EC','Venous EC','Lymphatic EC']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.8,'axes.titlesize':8.5,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,
    'legend.fontsize':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
PLATES=[];KEY=[]

def cmdata(stage):return read(OUT/(stage+'_proximity50_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id')
def base(stage):return read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
def export(fig,stem,message):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    for ext in ['png','pdf','svg']:
        p=folder/(stem+'.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight',pad_inches=.08)
        if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text().splitlines())+'\n')
    PLATES.append(stem);KEY.append({'figure':stem,'central_message':message,'raster_DPI':300,'formats':'PNG;PDF;SVG'})
    plt.close(fig)

def figure1():
    fig,axs=plt.subplots(3,2,figsize=(7.4,8.9),layout='constrained',gridspec_kw={'height_ratios':[1,1,.85]})
    summary=[]
    for j,stage in enumerate(['P0','P7']):
        d=cmdata(stage);f=base(stage);ax=axs[0,j];mapbase(ax,f,stage+' · ventricular CM states',chr(97+j))
        for state,color in CP.items():
            z=d[d.CM_state.eq(state)];ax.scatter(z.x_centroid,z.y_centroid,s=.5,c=color,lw=0,rasterized=True)
        ax=axs[1,j];mapbase(ax,f,stage+' · within-section maturation',chr(99+j))
        im=ax.scatter(d.x_centroid,d.y_centroid,c=d.stage_maturation_percentile,s=.55,cmap='viridis',vmin=0,vmax=1,lw=0,rasterized=True)
        fig.colorbar(im,ax=ax,shrink=.65,label='Within-stage score percentile')
        ax=axs[2,j];arrays=[d[d.CM_state.eq(s)].core_Z__maturation21.to_numpy() for s in CM]
        parts=ax.violinplot(arrays,positions=range(len(CM)),showmedians=True,showextrema=False)
        for state,b in zip(CM,parts['bodies']):b.set_facecolor(CP[state]);b.set_alpha(.65)
        parts['cmedians'].set_color('#222')
        ax.set(xticks=range(len(CM)),xticklabels=[CN[s] for s in CM],ylabel='Measured core21-gene index');ax.tick_params(axis='x',rotation=35)
        panel(ax,chr(101+j),stage+' · maturation varies within each identity')
        for state,g in d.groupby('CM_state'):
            q=g.core_Z__maturation21.quantile([.1,.5,.9]);summary.append({'stage':stage,'CM_state':state,'profiles':len(g),'P10':q.loc[.1],'median':q.loc[.5],'P90':q.loc[.9],
                'core_native_score_Spearman':spearmanr(g.core_Z__maturation21,g.native_Z__maturation21).statistic})
    fig.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,c=c,label=CN[s]) for s,c in CP.items()],loc='outside lower center',ncol=5,frameon=False)
    save(pd.DataFrame(summary),'P0_P7_CM_identity_maturation_distributions.tsv')
    export(fig,'Figure_1_CM_identity_and_maturation','CM identity and within-stage maturation are separate axes; no two-context cluster restriction')

def mean_matrix(ax,stage,label):
    d=read(OUT/'CM_state_maturation_subtype_neighborhood_means.tsv');coverage=read(OUT/'source_subtype_coverage.tsv').set_index(['stage','subtype'])
    rows=[(s,q) for s in CM for q in ['Low','High'] if ((d.stage==stage)&(d.CM_state==s)&(d.maturation_group==q)).any()]
    mat=np.full((len(rows),len(EX)),np.nan)
    for i,(state,q) in enumerate(rows):
        for j,ex in enumerate(EX):
            if coverage.loc[(stage,ex),'primary_source_profiles']<30:continue
            z=d[d.stage.eq(stage)&d.CM_state.eq(state)&d.maturation_group.eq(q)&d.source_subtype.eq(ex)]
            if len(z):mat[i,j]=100*z.mean_within_lineage_fraction.iloc[0]
    cmap=plt.colormaps['YlGnBu'].copy();cmap.set_bad('#e5e5e5');im=ax.imshow(np.ma.masked_invalid(mat),cmap=cmap,vmin=0,vmax=70,aspect='auto')
    for i in range(len(rows)):
        for j in range(len(EX)):
            if np.isfinite(mat[i,j]):ax.text(j,i,f'{mat[i,j]:.0f}',ha='center',va='center',fontsize=5.7,c='white' if mat[i,j]>40 else '#222')
    ax.set(xticks=range(len(EX)),xticklabels=XL,yticks=range(len(rows)),yticklabels=[CN[s]+' · '+q for s,q in rows]);ax.tick_params(axis='x',rotation=60)
    panel(ax,label,stage+' · subtype × maturation neighborhoods');return im

def figure2():
    fig,axs=plt.subplots(3,2,figsize=(7.4,9.7),layout='constrained')
    for i,stage in enumerate(STAGES):im=mean_matrix(axs[i//2,i%2],stage,chr(97+i))
    fig.colorbar(im,ax=list(axs[:2].ravel()),shrink=.45,pad=.01,label='Subtype share within source lineage (%)',extend='max')
    ax=axs[2,0];d=read(OUT/'source_subtype_coverage.tsv');mat=d.pivot(index='subtype',columns='stage',values='primary_source_profiles').reindex(index=EX,columns=STAGES)
    im=ax.imshow(np.log10(1+mat),cmap='Blues',aspect='auto',vmin=0,vmax=3.6)
    for i in range(len(EX)):
        for j in range(4):ax.text(j,i,str(int(mat.iloc[i,j])),ha='center',va='center',fontsize=6,c='white' if mat.iloc[i,j]>400 else 'black')
    ax.set(xticks=range(4),xticklabels=STAGES,yticks=range(len(EX)),yticklabels=XL);panel(ax,'e','Spatial source coverage · not abundance')
    ax=axs[2,1];d=read(OUT/'CM_state_maturation_subtype_neighborhood_means.tsv');d=d[d.CM_state.isin(['v.CM','Ankrd1.CM'])&d.maturation_group.isin(['High','Low'])]
    w=d.pivot(index=['stage','CM_state','source_subtype'],columns='maturation_group',values='mean_within_lineage_fraction');w['delta']=100*(w.High-w.Low)
    rows=[(s,t) for s in STAGES for t in ['v.CM','Ankrd1.CM']];mat=np.full((8,8),np.nan);cov=read(OUT/'source_subtype_coverage.tsv').set_index(['stage','subtype'])
    for i,(s,t) in enumerate(rows):
        for j,ex in enumerate(EX):
            if cov.loc[(s,ex),'primary_source_profiles']>=30:mat[i,j]=w.loc[(s,t,ex),'delta']
    cmap=plt.colormaps['RdBu_r'].copy();cmap.set_bad('#e5e5e5');im=ax.imshow(np.ma.masked_invalid(mat),cmap=cmap,vmin=-7,vmax=7,aspect='auto')
    ax.set(xticks=range(8),xticklabels=XL,yticks=range(8),yticklabels=[s+' '+CN[t] for s,t in rows]);ax.tick_params(axis='x',rotation=60)
    fig.colorbar(im,ax=ax,shrink=.6,pad=.01,label='High–low maturation difference (pp)');panel(ax,'f','Descriptive contrasts within CM identity')
    export(fig,'Figure_2_known_subtype_neighborhoods','Known FB and EC subtypes are analyzed around each CM identity and maturation group; unsupported cells are gray')

def coef_heat(ax,outcome,label,title):
    d=read(OUT/'subtype_conditioned_maturation_coefficients.tsv');d=d[d.geometry.eq('proximity50')&d.source_view.eq('primary')&d.adjustment.eq('within_region')&d.outcome.eq(outcome)]
    rows=[(s,t) for s in STAGES for t in CM if ((d.stage==s)&(d.CM_state==t)).any()];mat=np.full((len(rows),8),np.nan);qs=mat.copy()
    for i,(s,t) in enumerate(rows):
        for j,ex in enumerate(EX):
            z=d[d.stage.eq(s)&d.CM_state.eq(t)&d.source_subtype.eq(ex)]
            if len(z):mat[i,j]=z.standardized_coefficient.iloc[0];qs[i,j]=z.within_section_BH_q.iloc[0]
    cmap=plt.colormaps['RdBu_r'].copy();cmap.set_bad('#e5e5e5');im=ax.imshow(np.ma.masked_invalid(mat),cmap=cmap,vmin=-.08,vmax=.08,aspect='auto')
    for i in range(len(rows)):
        for j in range(8):
            if qs[i,j]<=.05:ax.text(j,i,'•',ha='center',va='center',fontsize=8)
    ax.set(xticks=range(8),xticklabels=XL,yticks=range(len(rows)),yticklabels=[s+' '+CN[t] for s,t in rows]);ax.tick_params(axis='x',rotation=60)
    panel(ax,label,title);return im

def sensitivity(ax,stage,state,ex,outcome,label,title):
    all=read(OUT/'subtype_conditioned_maturation_coefficients.tsv');p=read(OUT/'same_CM_geometry_sensitivity.tsv')
    z=all[all.stage.eq(stage)&all.CM_state.eq(state)&all.source_subtype.eq(ex)&all.outcome.eq(outcome)&all.adjustment.eq('within_region')&all.geometry.eq('proximity50')]
    r=[]
    for view,text in [('primary','50µm · nuclear EC'),('strict_cell_nuclear','50µm · strict EC'),('published_reference_comparator','50µm · public labels')]:
        v=z[z.source_view.eq(view)]
        if len(v):r.append((text,v.iloc[0]))
    z=p[p.stage.eq(stage)&p.CM_state.eq(state)&p.source_subtype.eq(ex)&p.outcome.eq(outcome)&p.paired_comparison.eq('proximity30')]
    for geom,text in [('proximity50','50µm · common CMs'),('proximity30','30µm · common CMs')]:
        v=z[z.geometry.eq(geom)]
        if len(v):r.append((text,v.iloc[0]))
    for i,(text,v) in enumerate(r):ax.errorbar(v.standardized_coefficient,i,xerr=[[v.standardized_coefficient-v.block_CI_low],[v.block_CI_high-v.standardized_coefficient]],fmt='o',c='#4477AA',ms=3,capsize=2)
    ax.axvline(0,c='gray',lw=.7);ax.set(yticks=range(len(r)),yticklabels=[t for t,_ in r],xlabel='Conditional coefficient (stage SD units)');ax.invert_yaxis();panel(ax,label,title)

def figure3():
    fig,axs=plt.subplots(3,2,figsize=(7.4,10.1),layout='constrained',gridspec_kw={'height_ratios':[1.4,1.4,.9]})
    im=coef_heat(axs[0,0],'core_Z__maturation21','a','Overall maturation · within CM state/region')
    coef_heat(axs[0,1],'core_Z__metabolic','b','Metabolic-regulatory RNA')
    coef_heat(axs[1,0],'core_Z__sarcomere','c','Sarcomere RNA')
    fig.colorbar(im,ax=[axs[0,0],axs[0,1],axs[1,0]],shrink=.35,pad=.01,label='Conditional coefficient (stage SD)',extend='both')
    sensitivity(axs[1,1],'P7','p.CM','FB__Postn-associated','core_Z__metabolic','d','P7 cycling CMs · Postn FB context')
    sensitivity(axs[2,0],'P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere','e','P14 Ankrd1-CMs · endocardial EC context')
    ax=axs[2,1];d=read(OUT/'Ankrd1_versus_working_CM_slope_differences.tsv');d=d[d.stage.eq('P14')&d.geometry.eq('proximity50')&d.source_view.eq('primary')&d.adjustment.eq('within_region')&d.outcome.eq('core_Z__sarcomere')]
    for i,r in enumerate(d.itertuples()):ax.errorbar(r.Ankrd1_minus_working_slope,i,xerr=[[r.Ankrd1_minus_working_slope-r.block_CI_low],[r.block_CI_high-r.Ankrd1_minus_working_slope]],fmt='o',c='#EE7733',ms=3,capsize=2)
    ax.axvline(0,c='gray',lw=.7);ax.set(yticks=range(len(d)),yticklabels=[XL[EX.index(x)] for x in d.source_subtype],xlabel='Ankrd1 minus working-CM slope');ax.invert_yaxis();panel(ax,'f','CM-state specificity · sarcomere at P14')
    export(fig,'Figure_3_CM_state_specific_associations','Subtype associations can be CM-state- and component-specific; local models and source/geometry sensitivity limit causal interpretation')

def figure4():
    old=ROOT/'results/cardiac_developmental_story/figures/Figure_3_fibroblast_state_remodeling';stem='Figure_4_fibroblast_specialization'
    for ext in ['png','pdf','svg']:shutil.copyfile(source(old.with_suffix('.'+ext)),OUT/'figures'/(stem+'.'+ext))
    PLATES.append(stem);KEY.append({'figure':stem,'central_message':'Reviewed fibroblast programs and candidate RNA continuity provide a temporal counterpart to subtype-resolved neighborhoods','raster_DPI':300,'formats':'PNG;PDF;SVG'})

def figure5():
    fig,axs=plt.subplots(3,2,figsize=(7.4,8.9),layout='constrained')
    ax=axs[0,0];ax.axis('off');panel(ax,'a','Condition transport on known CM identity')
    items=[('Focal unit','Stage + CM identity + measured maturation'),('Source environment','Known FB/EC subtype descriptors + molecular state'),('Allowed correspondence','Same CM identity and endocardial-distance category')]
    for y,(h,b) in zip([.79,.48,.17],items):
        ax.add_patch(FancyBboxPatch((.02,y-.12),.95,.25,boxstyle='round,pad=.01',fc='#f1f4f7',ec='#aab4bd',lw=.7))
        ax.text(.495,y+.04,h,ha='center',fontsize=8,weight='bold');ax.text(.495,y-.045,textwrap.fill(b,44),ha='center',va='center',fontsize=6.5)
    ax=axs[0,1];d=read(OUT/'CM_identity_restricted_OT_neighborhoods.tsv');mat=d.groupby(['stage','CM_state'])[EX].mean().reindex(pd.MultiIndex.from_product([STAGES,['v.CM','Ankrd1.CM']]))
    im=ax.imshow(100*mat,cmap='YlGnBu',vmin=0,vmax=70,aspect='auto');ax.set(xticks=range(8),xticklabels=XL,yticks=range(8),yticklabels=[s+' '+CN[t] for s,t in mat.index]);ax.tick_params(axis='x',rotation=60)
    panel(ax,'b','Known source descriptors in the transport subset');fig.colorbar(im,ax=ax,shrink=.6,pad=.01,label='Within-lineage source share (%)')
    d=read(OUT/'CM_identity_restricted_Wormhole_assessment.tsv')
    for ax,key,label,title in [(axs[1,0],'Wormhole_exact_rank_correlation','c','Distance rank calibration'),(axs[1,1],'nearest_exact_match_agreement','d','Fine nearest-neighborhood agreement')]:
        for i,state in enumerate(['v.CM','Ankrd1.CM']):
            for j,(a,b) in enumerate(zip(STAGES[:-1],STAGES[1:])):
                z=d[d.CM_state.eq(state)&d.stage_a.eq(a)&d.stage_b.eq(b)]
                ax.scatter(j+(i-.5)*.15+np.linspace(-.025,.025,len(z)),z[key],c=CP[state],s=18,label=CN[state] if j==0 else None)
        ax.set(xticks=range(3),xticklabels=['P0–P7','P7–P14','P14–P21'],ylim=(-.03,1.03),ylabel='Rank correlation' if label=='c' else 'Fraction agreeing with exact OT');ax.legend(frameon=False);panel(ax,label,title)
    ax=axs[2,0];d=read(OUT/'CM_state_within_section_neighborhood_distances.tsv');arrays=[];positions=[]
    for i,s in enumerate(STAGES):
        for k,kind in enumerate(['Within same CM state','Between CM states']):arrays.append(d[d.stage.eq(s)&d.pair_kind.eq(kind)].exact_W2_squared);positions.append(i+(k-.5)*.28)
    b=ax.boxplot(arrays,positions=positions,widths=.24,showfliers=False,patch_artist=True,medianprops={'color':'black'},whis=(10,90))
    for i,x in enumerate(b['boxes']):x.set_facecolor('#cccccc' if i%2==0 else '#EE7733')
    ax.set(xticks=range(4),xticklabels=STAGES,ylabel='Exact non-CM neighborhood W₂²');ax.legend(handles=[Rectangle((0,0),1,1,fc='#ccc',label='Within CM identity'),Rectangle((0,0),1,1,fc='#EE7733',label='Between identities')],frameon=False,fontsize=6);panel(ax,'e','Do CM states occupy different environments?')
    ax=axs[2,1];d=read(OUT/'CM_maturation_rank_correspondence_audit.tsv');d['difference']=d.exact_mean_absolute_maturation_rank_gap-d.all_candidate_pairs_mean_absolute_rank_gap
    for i,state in enumerate(['v.CM','Ankrd1.CM']):
        for j,(a,b) in enumerate(zip(STAGES[:-1],STAGES[1:])):
            z=d[d.CM_state.eq(state)&d.stage_a.eq(a)&d.stage_b.eq(b)]
            for r in z.itertuples():ax.scatter(j+(i-.5)*.15,r.difference,c=CP[state],marker='o' if r.distance_stratum=='≤100µm' else 's',s=24)
    ax.axhline(0,c='gray',lw=.7);ax.set(xticks=range(3),xticklabels=['P0–P7','P7–P14','P14–P21'],ylabel='Exact correspondence rank gap − all-pair gap');ax.text(.02,.98,'○ ≤100µm; □ >100µm\nBelow 0: closer held-out maturity ranks',transform=ax.transAxes,va='top',fontsize=6.2);panel(ax,'f','Maturation rank is held out from transport')
    export(fig,'Figure_5_CM_identity_conditioned_transport','Transport is conditioned on known CM identity; fine source descriptors and held-out maturation are audited rather than collapsed into broad contexts')

def figure6():
    fig,axs=plt.subplots(3,1,figsize=(7.4,7.6),layout='constrained')
    for ax in axs:ax.axis('off')
    ax=axs[0];panel(ax,'a','The analysis unit is a subtype-defined developmental neighborhood')
    columns=[('CM identity','Ankrd1; working; Slit2; Ddc; cycling'),('Maturation within stage','Continuous21-gene score and separate functional-RNA components'),('Non-CM subtypes','Postn / Gfpt2 / Col8a1 FBs; capillary / endocardial / arterial ECs')]
    for i,(title,body) in enumerate(columns):
        x=.02+i*.33;ax.add_patch(FancyBboxPatch((x,.18),.29,.66,boxstyle='round,pad=.012',fc='#f5f8fa',ec=['#EE7733','#4477AA','#228833'][i],lw=1));ax.text(x+.145,.76,textwrap.fill(title,24),ha='center',va='top',weight='bold',fontsize=8.5);ax.text(x+.145,.44,textwrap.fill(body,31),ha='center',va='center',fontsize=8)
    ax.text(.5,.04,'Broad lineage context is background; it is not a limit on subtype-defined analyses',ha='center',fontsize=7)
    ax=axs[1];panel(ax,'b','Component-specific candidates from the observed data')
    for i,(title,body,color) in enumerate([('P7 · cycling CM + Postn FB','Lower maturation/metabolic-regulatory RNA; annotation and radius sensitivity remain','#AA3377'),('P14 · Ankrd1-CM + endocardial EC','Higher sarcomere RNA in local models; not a uniform increase in all maturation components','#228833')]):
        x=.03+i*.49;ax.add_patch(FancyBboxPatch((x,.12),.45,.73,boxstyle='round,pad=.01',fc=color+'10',ec=color,lw=.9));ax.text(x+.225,.78,textwrap.fill(title,43),ha='center',va='top',weight='bold',fontsize=8.5,color=color);ax.text(x+.225,.43,textwrap.fill(body,48),ha='center',va='center',fontsize=8)
    ax=axs[2];panel(ax,'c','Choose subtype-matched functional and lineage experiments')
    for i,(title,body) in enumerate([('Fibroblast comparison','Postn versus Gfpt2 context around early cycling CMs; assess metabolism, matrix and completed division'),('EC × CM identity','Endocardial-like versus capillary ECs around Ankrd1 versus working CMs; sarcomeres, calcium and force'),('State continuity','Trace validated early FB labels into later matrix-specialized states; do not treat OT as ancestry')]):
        x=.02+i*.33;ax.add_patch(FancyBboxPatch((x,.10),.29,.75,boxstyle='round,pad=.01',fc='#f5f8fa',ec='#4477AA',lw=.8));ax.text(x+.145,.78,title,ha='center',va='top',weight='bold',fontsize=8);ax.text(x+.145,.44,textwrap.fill(body,31),ha='center',va='center',fontsize=7.5)
    export(fig,'Figure_6_subtype_developmental_framework','Known CM identities, within-stage maturation and specific FB/EC subtypes define the developmental questions and follow-up experiments')

def main():
    (OUT/'figures').mkdir(exist_ok=True)
    for fn in [figure1,figure2,figure3,figure4,figure5,figure6]:print('SUBTYPE_FIGURE',fn.__name__,flush=True);fn()
    # Merge the 300-dpi individual vector PDFs without degrading raster layers.
    save(pd.DataFrame(KEY),'figure_message_key.tsv')
    p=OUT/'source_files.json';SOURCES.update(json.loads(p.read_text()));p.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')

if __name__=='__main__':main()
