"""Integrated, publication-style panels for the developmental story.

All maps use original coordinates. Vector PDF/SVG and 300-dpi PNG are exported.
"""
from pathlib import Path
import os,json,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-story-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Rectangle,Ellipse
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr
from cardiac_story_analysis import ROOT,OUT,STAGES,read,save,source,SOURCES

COL={'Cardiomyocyte':'#4477AA','Endothelial':'#228833','Fibroblast':'#AA3377','Pericyte':'#AA9944','Smooth muscle':'#EE7733','Immune':'#CC3311','Epicardial':'#55AABB'}
FB={'Postn-associated':'#DDAA33','Col8a1-associated':'#4477AA','Gfpt2-associated':'#228833','Lsamp-associated':'#AA3377','Kirrel3-associated':'#66CCEE'}
EC={'Capillary-like':'#4477AA','Endocardial-like':'#DDAA33','Arterial-like':'#CC3311','Venous-like':'#228833','Lymphatic-like':'#AA3377'}
AGE=dict(zip(STAGES,['#0072B2','#E69F00','#009E73','#CC79A7']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.8,'axes.titlesize':8.5,'axes.labelsize':7.5,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,'lines.linewidth':1.2,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
FIGURES=[];MANIFEST=[]

def panel(ax,label,title):
    ax.text(-.12,1.07,label,transform=ax.transAxes,fontsize=12,weight='bold',va='bottom')
    ax.set_title(title,pad=5,loc='left',weight='normal')

def base(stage):return read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)

def mapbase(ax,f,title,label):
    ax.scatter(f.x_centroid,f.y_centroid,s=.12,c='#e6e6e6',lw=0,rasterized=True)
    ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
    for sp in ax.spines.values():sp.set_visible(False)
    xmin,xmax=f.x_centroid.min(),f.x_centroid.max();ymin,ymax=f.y_centroid.min(),f.y_centroid.max()
    x=xmin+.06*(xmax-xmin);y=ymax-.05*(ymax-ymin)
    ax.plot([x,x+500],[y,y],c='#202020',lw=1.8);ax.text(x+250,y-.02*(ymax-ymin),'500 µm',ha='center',fontsize=6)
    panel(ax,label,title)

def export(fig,stem,message,panels):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:
        path=folder/(stem+'.'+ext);fig.savefig(path,dpi=300,bbox_inches='tight',pad_inches=.08)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    FIGURES.append((fig,stem));MANIFEST.append({'figure':stem,'message':message,'panels':'; '.join(panels),'PNG_DPI':300,'vector_formats':'PDF;SVG'})

def figure1():
    fig=plt.figure(figsize=(7.4,8.3),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,1.05])
    for i,stage in enumerate(STAGES):
        ax=fig.add_subplot(gs[i//2,i%2]);f=base(stage);mapbase(ax,f,stage,chr(97+i))
        for name,color in COL.items():
            z=f[f.atlas_identity.eq(name)];ax.scatter(z.x_centroid,z.y_centroid,s=.22,c=color,lw=0,rasterized=True)
    ax=fig.add_subplot(gs[2,0]);d=read('results/cardiac_spatial_pipeline/contact_enrichment_local_null.tsv')
    mat=np.array([[d[d.stage.eq(s)&d.cell_type_1.eq(k)&d.cell_type_2.eq(k)].log2_observed_over_local_null.iloc[0] for s in STAGES] for k in COL])
    im=ax.imshow(mat,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=2.5),aspect='auto');ax.set(xticks=range(4),xticklabels=STAGES,yticks=range(7),yticklabels=list(COL))
    for i in range(7):
        for j in range(4):ax.text(j,i,f'{2**mat[i,j]:.1f}',ha='center',va='center',fontsize=7,color='white' if abs(mat[i,j])>.65 else 'black')
    panel(ax,'e','Local same-lineage contact enrichment');fig.colorbar(im,ax=ax,shrink=.65,label='log₂(observed / local null)')
    ax=fig.add_subplot(gs[2,1]);f=base('P14');n=read('results/cardiac_spatial_pipeline/P14/niche_annotations.tsv.gz',index_col=0);mapbase(ax,f,'Two stable broad contexts · P14','f')
    for prefix,color in [('N1:','#AA3377'),('N2:','#4477AA')]:
        z=f[n.niche.str.startswith(prefix)];ax.scatter(z.x_centroid,z.y_centroid,s=.3,c=color,lw=0,rasterized=True)
    ax.text(.02,.98,'Blue: CM/EC enriched\nPurple: fibroblast/mural enriched',transform=ax.transAxes,fontsize=6.5,va='top',bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    fig.legend(handles=[Line2D([],[],ls='',marker='o',ms=3,color=c,label=k) for k,c in COL.items()],loc='outside lower center',ncol=4,frameon=False)
    export(fig,'Figure_1_persistent_spatial_architecture','Persistent myocardial, vascular and surface/stromal organization anchors the developmental comparison',['a–d Native all-lineage section maps','e Homotypic contacts relative to200µm local shuffles; numbers are ratios','f Existing seed-stable broad context map; no new niche classes'])

def expression_map(ax,stage,gene,label,title,vmax):
    import anndata as ad
    f=base(stage);mapbase(ax,f,title,label)
    a=ad.read_h5ad(source(f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad'),backed='r')
    x=a.X[:];vals=np.asarray(x[:,a.var_names.get_loc(gene)].toarray()).ravel();depth=np.asarray(x.sum(1)).ravel()
    score=pd.Series(np.log1p(vals*1e4/np.maximum(depth,1)),index=a.obs_names).reindex(f.index)
    keep=f.ventricular_CM_strict&f.author_atlas_cell_state.str.contains('CM',na=False)
    z=f[keep];order=np.argsort(score.loc[z.index].to_numpy())
    im=ax.scatter(z.x_centroid.iloc[order],z.y_centroid.iloc[order],c=score.loc[z.index].iloc[order],s=.8,cmap='viridis',vmin=0,vmax=vmax,lw=0,rasterized=True)
    if gene=='Ankrd1':
        data=json.loads(source(f'results/whole_heart_Ankrd1_vessels/{stage}/arterial_segments.geojson').read_text())
        for feature in data['features']:
            p=feature['properties']
            if p['definition']!='independent_strict' or not p['primary_small_segment']:continue
            g=feature['geometry'];polys=g['coordinates'] if g['type']=='MultiPolygon' else [g['coordinates']]
            for poly in polys:
                xy=np.asarray(poly[0]);ax.plot(xy[:,0],xy[:,1],c='black',lw=.35)
    else:
        z=f[f.author_atlas_cell_state.eq('endocardial.EC')];ax.scatter(z.x_centroid,z.y_centroid,s=.3,c='#CC3311',lw=0,rasterized=True)
    a.file.close();return im

def figure2():
    fig=plt.figure(figsize=(7.4,8.1),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,.9])
    for i,(stage,gene) in enumerate([('P7','Slit2'),('P14','Slit2'),('P7','Ankrd1'),('P14','Ankrd1')]):
        ax=fig.add_subplot(gs[i//2,i%2]);im=expression_map(ax,stage,gene,chr(97+i),f'{stage} · measured {gene}',5 if gene=='Slit2' else 8)
        fig.colorbar(im,ax=ax,shrink=.6,pad=.01,label='logₑ(1 + panel CP10k)',extend='max')
    ax=fig.add_subplot(gs[2,0]);d=read('results/author_atlas_baseline/descriptive_distances_to_author_EC_states.tsv');d=d[d.landmark_author_state.eq('endocardial.EC')]
    for state,color in [('v.CM','#4477AA'),('Slit2.CM','#228833')]:
        z=d[d.author_CM_state.eq(state)].set_index('stage').loc[STAGES];ax.plot(range(4),z.median_nearest_distance_um,'o-',c=color,label=state)
    ax.set(xticks=range(4),xticklabels=STAGES,ylabel='Median endocardial-EC distance (µm)');ax.legend(frameon=False);panel(ax,'e','Endocardial-associated CM identity persists')
    ax=fig.add_subplot(gs[2,1]);d=read('results/whole_heart_Ankrd1_vessels/primary_local_Ankrd1_results.tsv').set_index('stage').loc[STAGES]
    ax.bar(range(4),d.near_far_Ankrd1_CP10k_ratio,color=[AGE[s] for s in STAGES],width=.6);ax.axhline(1,c='#444',lw=.7)
    for i,r in enumerate(d.itertuples()):ax.text(i,r.near_far_Ankrd1_CP10k_ratio+.03,f'{r.matched_pairs:,}\npairs',ha='center',fontsize=6.5)
    ax.set(xticks=range(4),xticklabels=STAGES,ylabel='Near / matched-control Ankrd1 RNA',ylim=(0,1.95));panel(ax,'f','A partial peri-arterial response program')
    export(fig,'Figure_2_regional_cardiomyocyte_programs','Persistent regional CM programs overlap with local vascular response states',['a,b Native measured Slit2; red=endocardial EC reference','c,d Native measured Ankrd1; black=candidate arterial envelopes','e Public-label distance comparison, not independent annotation validation','f Frozen local matched Ankrd1 contrasts; no animal-level error bars'])

def draw_FB_map(ax,stage,label):
    f=base(stage);mapbase(ax,f,stage+' · supported fibroblast families',label)
    for family,color in FB.items():
        z=f[f.atlas_identity.eq('Fibroblast')&f.fine_identity.eq(family)]
        ax.scatter(z.x_centroid,z.y_centroid,s=1.9,c=color,lw=0,rasterized=True)

def figure3():
    fig=plt.figure(figsize=(7.4,9.1),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1.1,1.1])
    for i,stage in enumerate(['P7','P14']):draw_FB_map(fig.add_subplot(gs[0,i]),stage,chr(97+i))
    ax=fig.add_subplot(gs[1,0]);d=read(OUT/'fibroblast_library_matrix_programs.tsv');d=d[d.view.eq('corrected')&d.level.eq('fibroblast_family')&~d.group.str.startswith('Unresolved')]
    order=[(s,f) for s in STAGES for f in FB if ((d.stage==s)&(d.group==f)).any()]
    programs=['Early developmental','Matrix remodeling','Fibrillar matrix','Basement matrix','Elastic matrix','Cycling']
    mat=d.groupby(['stage','group','program']).mean_log1p_CPM.mean().unstack().reindex(index=pd.MultiIndex.from_tuples(order),columns=programs)
    zz=(mat-mat.mean())/mat.std(ddof=0);im=ax.imshow(zz,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set(xticks=range(len(programs)),xticklabels=['Early','Remodel','Fibrillar','Basement','Elastic','Cycling'],yticks=range(len(order)),yticklabels=[s+' '+f.split('-')[0] for s,f in order]);ax.tick_params(axis='x',rotation=45)
    fig.colorbar(im,ax=ax,shrink=.7,label='Within-program z score');panel(ax,'c','RNA matrix programs change within families')
    ax=fig.add_subplot(gs[1,1]);d=read('results/fibroblast_names_and_neighborhoods/matched_family_context_sensitivity.tsv');d=d[d.stage.eq('P14')&d.radius_um.eq(30)]
    keys=['within_50um_arterial_envelope','fraction_Smooth muscle_30um'];available=d.measure.to_list()
    # Use exact frozen column names, never replace missing values by zeros.
    arterial=d[d.measure.eq('arterial_within_50um')].iloc[0];smooth=d[d.measure.eq('fraction_Smooth muscle_30um')].iloc[0]
    for j,(r,name) in enumerate([(arterial,'Near arterial\nenvelope'),(smooth,'Smooth-muscle\nneighbors')]):
        ax.bar(j-.16,100*r.target_matched_mean,width=.3,color=FB['Col8a1-associated']);ax.bar(j+.16,100*r.reference_matched_mean,width=.3,color=FB['Gfpt2-associated'])
    ax.set(xticks=[0,1],xticklabels=['Near arterial\nenvelope','Smooth-muscle\nneighbors'],ylabel='Matched neighborhood metric (%)')
    ax.legend(handles=[Rectangle((0,0),1,1,color=FB[f],label=f.split('-')[0]) for f in ['Col8a1-associated','Gfpt2-associated']],frameon=False)
    panel(ax,'d','Col8a1 has a local vascular-associated component')
    ax=fig.add_subplot(gs[2,0]);d=read(OUT/'fibroblast_library_marker_profiles.tsv');names=['Col8a1/Svep1 fibroblasts','Col8a1/Hmcn2 fibroblasts','Lsamp/Slit3 fibroblasts','Lsamp/Col4a5 fibroblasts'];genes=['Svep1','Fap','Hmcn2','Col8a1','Lsamp','Slit3','Eln','Col4a5','Col4a6']
    d=d[d.view.eq('corrected')&d.stage.eq('P21')&d.level.eq('named_state')&d.group.isin(names)&d.gene.isin(genes)]
    val=d.groupby(['group','gene']).log1p_CPM.mean().unstack().reindex(index=names,columns=genes);det=d.groupby(['group','gene']).fraction_detected.mean().unstack().reindex(index=names,columns=genes)
    zz=(val-val.mean())/val.std(ddof=0)
    for i in range(len(names)):
        for j in range(len(genes)):ax.scatter(j,i,s=8+100*det.iloc[i,j],c=[zz.iloc[i,j]],cmap='RdBu_r',vmin=-1.8,vmax=1.8,edgecolors='none')
    ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=range(4),yticklabels=[s.replace(' fibroblasts','') for s in names],ylim=(3.6,-.6));ax.tick_params(axis='x',rotation=60)
    from matplotlib.cm import ScalarMappable
    fig.colorbar(ScalarMappable(norm=plt.Normalize(-1.8,1.8),cmap='RdBu_r'),ax=ax,shrink=.55,pad=.01,label='Gene-wise z')
    panel(ax,'e','Late RNA substates exceed spatial resolution')
    ax=fig.add_subplot(gs[2,1]);ax.axis('off');panel(ax,'f','Candidate continuity, not lineage conversion')
    positions={'P7 Postn':(.12,.73),'P14 Col8a1':(.49,.73),'P21 Col8a1':(.87,.73),'P7 Gfpt2':(.12,.28),'P14 Gfpt2':(.49,.28),'P21 Gfpt2':(.87,.28)}
    for name,(x,y) in positions.items():ax.text(x,y,name.replace(' ','\n'),ha='center',va='center',fontsize=7.5,bbox={'boxstyle':'round,pad=.3','fc':'#f4f5f7','ec':'#999','lw':.6})
    for a,b,r in [('P7 Postn','P14 Col8a1',.410),('P14 Col8a1','P21 Col8a1',.697),('P7 Gfpt2','P14 Gfpt2',.554),('P14 Gfpt2','P21 Gfpt2',.689)]:
        x,y=positions[a];u,v=positions[b];ax.annotate('',xy=(u-.09,v),xytext=(x+.09,y),arrowprops={'arrowstyle':'->','ls':'--','lw':1,'color':'#555'});ax.text((x+u)/2,y+.13,f'r={r:.2f}',ha='center',fontsize=7)
    ax.text(.5,.03,'Selected gene-withheld RNA signature matches\nDashed arrows = hypotheses; widths encode no mass',ha='center',va='bottom',fontsize=6.5,transform=ax.transAxes)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,c=c,label=k.split('-')[0]) for k,c in FB.items()],loc='outside lower center',ncol=5,frameon=False)
    export(fig,'Figure_3_fibroblast_state_remodeling','Fibroblast programs specialize within overlapping interstitial environments; RNA continuity does not establish ancestry',['a,b Strict cell+nuclear family subsets; gray includes unresolved fibroblasts','c Corrected full-RNA library-family pseudobulks; column z scores','d Frozen coverage/morphology-matched P14 Col8a1 versus Gfpt2 metrics','e Named P21 RNA substates, colors relative expression, size detection','f Selected frozen gene-withheld signature correlations; reviewed P21_F0 family used'])

def EC_map(ax,stage,label):
    f=base(stage);mapbase(ax,f,stage+' · strict endothelial subtype subset',label)
    d=read('results/endothelial_five_types/Xenium_five_type_annotations.tsv.gz');d=d[d.stage.eq(stage)]
    for name,color in EC.items():
        z=d[d.EC_type.eq(name)];ax.scatter(z.x_centroid,z.y_centroid,s=3,c=color,lw=0,rasterized=True)

def figure4():
    fig=plt.figure(figsize=(7.4,8.7),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,1])
    for j,stage in enumerate(['P7','P14']):EC_map(fig.add_subplot(gs[0,j]),stage,chr(97+j))
    ax=fig.add_subplot(gs[1,0]);d=read('results/cardiac_spatial_reconstruction/RNA_observed_library_programs.tsv');d=d[d.scope.eq('QC_Endothelial')&d.program.eq('EC_cycling')]
    for i,stage in enumerate(STAGES):
        vals=d[d.stage.eq(stage)].pseudobulk_logCPM_mean.to_numpy();ax.scatter(np.full(len(vals),i)+[-.06,.06],vals,s=20,c=AGE[stage]);ax.plot([i,i],[vals.min(),vals.max()],c=AGE[stage])
    ax.set(xticks=range(4),xticklabels=STAGES,ylabel='Cycling module · mean logₑ(1 + CPM)');panel(ax,'c','EC cycling is concentrated in early development')
    ax=fig.add_subplot(gs[1,1]);d=read('results/cardiac_spatial_reconstruction/nuclear_cycling_specificity_coverage.tsv');d=d[d.geometry.eq('within50um')].set_index('stage').loc[STAGES]
    ax.plot(range(4),100*d.CM_fraction_two_cycling_genes_detected,'o-',c=COL['Cardiomyocyte'],label='CM');ax.plot(range(4),100*d.EC_fraction_two_cycling_genes_detected,'s-',c=COL['Endothelial'],label='EC')
    ax.set(xticks=range(4),xticklabels=STAGES,ylabel='Profiles detecting ≥2 cycling genes (%)');ax.legend(frameon=False);panel(ax,'d','Early cycling gives way to a late low-cycling floor')
    ax=fig.add_subplot(gs[2,0]);d=read('results/cardiac_spatial_reconstruction/primary_matched_component_results.tsv')
    for k,(comp,label,color) in enumerate([('metabolic','Metabolic-regulatory','#4477AA'),('withdrawal','Lower cycling','#CC6677')]):
        z=d[d.outcome.eq('core_Z__'+comp)].set_index('stage').loc[STAGES]
        ax.errorbar(np.arange(4)+(k-.5)*.15,z.mean_standardized_difference,yerr=[z.mean_standardized_difference-z.block_CI_low,z.block_CI_high-z.mean_standardized_difference],fmt='o-',ms=3.5,capsize=2,c=color,label=label)
    ax.axhline(0,c='gray',lw=.7);ax.set(xticks=range(4),xticklabels=STAGES,ylabel='Matched high–zero apposition difference (SD)');ax.legend(frameon=False,fontsize=6.5);panel(ax,'e','Apposition relates to distinct CM components')
    ax=fig.add_subplot(gs[2,1]);d=read('results/cardiac_spatial_reconstruction/nuclear_cycling_specificity_coefficients.tsv')
    for k,(key,label,color) in enumerate([('within50um','Within50 µm','#228833'),('annulus30to50um','30–50 µm shell','#999999')]):
        z=d[d.geometry.eq(key)&d.outcome.eq('CM_nuclear_cycling')].set_index('stage').loc[STAGES]
        ax.errorbar(np.arange(4)+(k-.5)*.15,z.standardized_coefficient,yerr=[z.standardized_coefficient-z.block_CI_low,z.block_CI_high-z.standardized_coefficient],fmt='o-',ms=3.5,capsize=2,c=color,label=label)
    ax.axhline(0,c='gray',lw=.7);ax.set(xticks=range(4),xticklabels=STAGES,ylabel='EC cycling–CM cycling association (SD)');ax.legend(frameon=False,fontsize=6.5);panel(ax,'f','A separate, very local early growth association')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,c=c,label=k.replace('-like','')) for k,c in EC.items()],loc='outside lower center',ncol=5,frameon=False)
    export(fig,'Figure_4_vascular_growth_and_CM_components','Vascular growth and integration relate to different dimensions of CM development',['a,b Strict nuclear+cell EC subtype subsets; sparse calls are not absence','c Two observed RNA libraries per age, no animal-level intervals','d Nuclear-QC subsets only; expression is not completed division','e Strict local matched core-RNA contrasts;400µm block95% intervals','f Nuclear EC/CM cycling models;400µm block95% intervals'])

def figure5():
    fig=plt.figure(figsize=(7.4,8.6),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[.95,1,1.05])
    ax=fig.add_subplot(gs[0,0]);ax.axis('off');panel(ax,'a','Compare tissue neighborhoods, not cell fates')
    for y,title,body in [(.83,'96 reserved CM-centered neighborhoods','24/age; CM/EC-rich; two endocardial-distance strata'),(.51,'Exact non-CM distribution distance','16 molecular PCs +6 lineage features; CM outcome withheld'),(.19,'Restricted cross-stage correspondence','Uniform sampled mass; no abundance or lineage inference')]:
        ax.add_patch(FancyBboxPatch((.02,y-.13),.96,.25,boxstyle='round,pad=.015',facecolor='#f1f4f7',edgecolor='#9da9b5',lw=.7,transform=ax.transAxes))
        ax.text(.5,y+.035,title,ha='center',va='center',weight='bold',fontsize=7.5,transform=ax.transAxes);ax.text(.5,y-.045,textwrap.fill(body,47),ha='center',va='center',fontsize=6.5,transform=ax.transAxes)
    ax=fig.add_subplot(gs[0,1]);d=read('results/cardiac_developmental_synthesis/cross_stage_metric_assessment.tsv')
    for i,(metric,label,color) in enumerate([('Wormhole_vs_Sinkhorn_Spearman','Wormhole','#4477AA'),('mean_vs_Sinkhorn_Spearman','Mean features','#999999')]):
        for j,mode in enumerate(['random_cross_stage','composition_matched_cross_stage']):
            v=d[d.pair_selection.eq(mode)][metric].to_numpy();ax.scatter(j+(i-.5)*.15+np.linspace(-.025,.025,len(v)),v,c=color,s=15,label=label if j==0 else None)
    ax.set(xticks=[0,1],xticklabels=['Broad mixtures','Composition-similar'],ylabel='Rank correlation with Sinkhorn distance',ylim=(.4,1));ax.legend(frameon=False);panel(ax,'b','Fine state comparison is more demanding')
    ax=fig.add_subplot(gs[1,0]);d=read(OUT/'anatomy_stratified_neighborhood_pair_costs.tsv');w=read(OUT/'within_section_neighborhood_distance_baseline.tsv')
    comparisons=list(zip(STAGES[:-1],STAGES[1:]));data=[];positions=[]
    for i,(sa,sb) in enumerate(comparisons):
        data.extend([w[w.stage.isin([sa,sb])].exact_W2_squared.to_numpy(),d[d.stage_a.eq(sa)&d.stage_b.eq(sb)].exact_W2_squared.to_numpy()]);positions.extend([i-.16,i+.16])
    bp=ax.boxplot(data,positions=positions,widths=.25,showfliers=False,patch_artist=True,whis=(10,90),medianprops={'color':'black','lw':.8})
    for i,b in enumerate(bp['boxes']):b.set_facecolor('#cccccc' if i%2==0 else '#4477AA');b.set_linewidth(.6)
    ax.set(xticks=range(3),xticklabels=['P0–P7','P7–P14','P14–P21'],ylabel='Exact neighborhood W₂²');ax.legend(handles=[Rectangle((0,0),1,1,fc='#ccc',label='Within-section'),Rectangle((0,0),1,1,fc='#4477AA',label='Across-stage')],frameon=False,fontsize=6.5);panel(ax,'c','Regional variability overlaps stage differences')
    ax=fig.add_subplot(gs[1,1]);d=read(OUT/'composition_similar_neighborhood_cost_decomposition.tsv')
    for j,stratum in enumerate(d.distance_stratum.unique()):
        z=d[d.distance_stratum.eq(stratum)].set_index(['stage_a','stage_b']).loc[comparisons];ax.bar(np.arange(3)+(j-.5)*.32,100*z.median_state_fraction_of_cost,width=.3,color=['#4477AA','#88CCEE'][j],label='≤100 µm' if j==0 else '>100 µm')
    ax.set(xticks=range(3),xticklabels=['P0–P7','P7–P14','P14–P21'],ylim=(0,100),ylabel='Molecular-state share of selected cost (%)');ax.legend(frameon=False,title='Endocardial proxy',title_fontsize=6.5);panel(ax,'d','Similar mixtures retain molecular-state differences')
    ax=fig.add_subplot(gs[2,0]);d=read(OUT/'neighborhood_transport_metric_sensitivity.tsv')
    for key,label,color in [('state25','25% state','#228833'),('state75','75% state','#AA3377'),('subsample16','16 source profiles','#DDAA33'),('Wormhole','Wormhole','#4477AA'),('mean_features','Mean features','#999999')]:
        z=d[d.metric.eq(key)];ax.scatter(100*z.transport_plan_total_variation_from_primary,100*z.primary_exact_cost_inflation,s=20,c=color,label=label)
    ax.set(xlabel='Correspondence mass reassigned (%)',ylabel='Increase in primary exact cost (%)');ax.legend(frameon=False,fontsize=6,loc='upper left');panel(ax,'e','Fine correspondence is not uniquely identified')
    ax=fig.add_subplot(gs[2,1]);z=d[d.metric.eq('Wormhole')]
    for i,(sa,sb) in enumerate(comparisons):
        q=z[z.stage_a.eq(sa)&z.stage_b.eq(sb)];ax.scatter(np.full(len(q),i)+[-.06,.06],100*q.nearest_exact_neighbor_agreement,c='#4477AA',s=28)
    ax.set(xticks=range(3),xticklabels=['P0–P7','P7–P14','P14–P21'],ylim=(0,100),ylabel='Exact nearest-neighborhood agreement (%)');panel(ax,'f','Wormhole does not establish a unique transition map')
    export(fig,'Figure_5_transport_constrained_correspondence','Optimal transport supports distributional comparison but not a unique developmental fate map',['a Explicit cohort and hierarchy of transport problems','b Earlier768-pair reserved cross-stage metric audit','c New within-section/adjacent-stage distance distributions; pairs are dependent','d Composition-similar pairs L1≤0.3; equal-weight representation makes this conditional, not variance explained','e Alternative correspondences can differ with little extra exact cost','f Saved Wormhole nearest matches against exact OT within coarse anatomical strata'])

def box(ax,xy,w,h,title,body,color='#4477AA'):
    x,y=xy;ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015',fc=color+'12',ec=color,lw=.9))
    ax.text(x+w/2,y+h-.055,textwrap.fill(title,23),ha='center',va='top',fontsize=8.2,weight='bold',color=color)
    ax.text(x+w/2,y+h/2-.045,textwrap.fill(body,28),ha='center',va='center',fontsize=7.5,linespacing=1.4)

def figure6():
    fig=plt.figure(figsize=(7.4,7.8),layout='constrained');gs=fig.add_gridspec(3,1,height_ratios=[1.1,1,.9])
    ax=fig.add_subplot(gs[0]);ax.axis('off');panel(ax,'a','Working model · persistent places, changing states')
    columns=[('P0–P7','Growth and matrix remodeling','Early EC/CM cycling; Postn-associated interstitium'),('P7–P14','Regional integration','Candidate cell-cycle withdrawal and matrix specialization'),('P14–P21','Continued specialization','Metabolic-regulatory CM programs; low cycling')]
    for i,(age,title,body) in enumerate(columns):
        left=.02+i*.33;color=list(AGE.values())[min(i,3)]
        ax.add_patch(FancyBboxPatch((left,.12),.29,.78,boxstyle='round,pad=.012',fc='#f8fafb',ec=color,lw=.9))
        ax.text(left+.145,.85,age,ha='center',va='top',fontsize=9,weight='bold',color=color)
        # Conceptual tissue unit: CMs, a capillary/EC track, stromal cell/matrix.
        for yy in [.51,.65]:
            width=.11+i*.025
            ax.add_patch(Ellipse((left+.145,yy),width,.085,fc='#4477AA25',ec=COL['Cardiomyocyte'],lw=.8))
            ax.add_patch(Ellipse((left+.145,yy),.018,.026,fc=COL['Cardiomyocyte'],ec='none'))
        xx=np.linspace(left+.035,left+.255,60);yy=.41+.018*np.sin((xx-left)*35)
        ax.plot(xx,yy,c=COL['Endothelial'],lw=3,solid_capstyle='round')
        for x in [left+.07,left+.145,left+.22]:ax.add_patch(Ellipse((x,.41+.018*np.sin((x-left)*35)),.026,.037,fc='white',ec=COL['Endothelial'],lw=.7))
        if i==0:
            ax.plot([left+.11,left+.11,left+.08],[.42,.48,.5],c=COL['Endothelial'],lw=1.5)
            ax.text(left+.035,.74,'Cycling / growth',fontsize=6.4,color=COL['Endothelial'])
        ax.add_patch(Ellipse((left+.07,.57),.055,.022,angle=30,fc='#AA337730',ec=COL['Fibroblast'],lw=.8))
        ax.plot([left+.035,left+.09,left+.2,left+.26],[.31,.34,.3,.33],c='#AA3377',ls='--',lw=.8)
        ax.text(left+.145,.23,textwrap.fill(title,26),ha='center',va='center',fontsize=7.3)
    ax.text(.5,.015,'Conceptual tissue units: blue CM; green EC/capillary; purple fibroblast/matrix. Shapes and stages are not morphometric measurements.',ha='center',va='bottom',fontsize=6.2)
    ax=fig.add_subplot(gs[1]);ax.axis('off');panel(ax,'b','Anatomy conditions developmental responses')
    for i,(title,body,color) in enumerate([('Endocardial interface','Persistent Slit2-associated CM context; regional identity is not simply immaturity','#228833'),('Arterial / mural interface','Partial Ankrd1 response and Col8a1-associated stromal component','#EE7733'),('Myocardial interstitium','CM–capillary–pericyte neighborhoods with changing fibroblast and vascular states','#4477AA')]):box(ax,(.02+i*.33,.09),.29,.78,title,body,color)
    ax=fig.add_subplot(gs[2]);ax.axis('off');panel(ax,'c','A discriminating experiment · oxygen × EC contact × CM state')
    for i,(title,body,color) in enumerate([('Delivery support','Measure oxygen and respiration; equal EC number and matrix','#4477AA'),('Direct contact','Direct versus separated ECs; verify interface and viability','#228833'),('CM competence','Early versus later starting states; calcium and contractile endpoints','#AA3377')]):box(ax,(.02+i*.33,.13),.29,.72,title,body,color)
    ax.text(.5,.01,'Independent donor preparations are inference units; perturb a candidate receptor only after a contact benefit is demonstrated',ha='center',fontsize=6.7)
    export(fig,'Figure_6_developmental_model_and_tests','Development is framed as regional, asynchronous remodeling of a multicellular tissue unit; causal alternatives remain testable',['a Proposed developmental coordination, not observed longitudinal transitions','b Anatomical contexts integrating observed maps and candidate functions','c Prespecified oxygen/contact/state experiment with functional endpoints'])

def main():
    for fn in [figure1,figure2,figure3,figure4,figure5,figure6]:
        print('STORY_FIGURE',fn.__name__,flush=True);fn()
    with PdfPages(OUT/'Cardiac_development_integrated_figures.pdf') as pdf:
        for fig,_ in FIGURES:pdf.savefig(fig,dpi=300,bbox_inches='tight',pad_inches=.08)
    for fig,_ in FIGURES:plt.close(fig)
    save(pd.DataFrame(MANIFEST),'figure_message_and_panel_key.tsv')
    p=OUT/'analysis_source_files.json';SOURCES.update(json.loads(p.read_text()));p.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')

if __name__=='__main__':main()
