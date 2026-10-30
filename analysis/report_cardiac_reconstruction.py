"""An evidence-led spatial reconstruction, component results and experiment plan."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-reconstruction-mpl')
from reconstruct_cardiac_spatial import ROOT,OUT,CFG,STAGES,read,save,dump,sha
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import textwrap,json

plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
COLORS=dict(zip(STAGES,['#0072b2','#e69f00','#009e73','#cc79a7']))
COMP={'metabolic':'Metabolic-regulatory RNA','sarcomere':'Sarcomere RNA','excitation_contraction':'Excitation/contraction RNA','withdrawal':'Lower cycling RNA'}

def write(fig,name,pdf=None,foot=None):
    if foot:fig.text(.03,-.055,textwrap.fill(foot,175),fontsize=8,va='top',color='#454545')
    fig.savefig(OUT/(name+'.png'),dpi=150,bbox_inches='tight')
    if pdf is not None:pdf.savefig(fig,bbox_inches='tight')
    plt.close(fig)

def textpage(title,sections,name,pdf):
    fig=plt.figure(figsize=(12,10));fig.text(.06,.94,title,fontsize=19,weight='bold');y=.86
    wrapped=[textwrap.fill(s,112,break_long_words=False,break_on_hyphens=False) for _,s in sections]
    total=sum(s.count('\n')+1 for s in wrapped);step=min(.021,(.79-.050*len(sections))/total);font=min(10.8,step*720/1.35)
    for (heading,_),body in zip(sections,wrapped):
        fig.text(.06,y,heading,fontsize=12,weight='bold',va='top');y-=.033
        fig.text(.06,y,body,fontsize=font,va='top',linespacing=1.35);y-=step*(body.count('\n')+1)+.017
    assert y>=.065
    write(fig,name,pdf)

def selections_and_decision():
    matched=read(OUT/'matched_apposition_component_differences.tsv',index_col=None)
    primary=matched[matched.outcome.isin(['core_Z__metabolic','core_Z__withdrawal'])]
    save(primary,'primary_matched_component_results.tsv')
    priorities=[
        {'priority':1,'hypothesis':'The CM phenotype associated with capillary apposition differs by developmental window: P7 cell-cycle withdrawal versus later metabolic-regulatory transcription',
         'existing_evidence':'Strict same-region/CM-state/capture/shape/capillary-area matched comparisons plus conditional native/core geometry analyses',
         'unresolved':'One spatial section per age; interface and oxygen delivery not experimentally separated; late cycling floor; strict capillary identity attenuates some signals',
         'experiment':'Oxygen condition x direct-versus-separated EC contact x initial CM developmental state, with matched matrix/EC density and donor-paired allocation',
         'falsification':'No reproducible contact or oxygen-dependent functional effect after comparable viability and baseline state; revise toward shared regional development or other cues'},
        {'priority':2,'hypothesis':'EC PTPRM is necessary for a contact-dependent CM effect',
         'existing_evidence':'Atlas protein/perturbation evidence; present mean sender RNA is weak after controls, receiver interactions are secondary and stage dependent',
         'unresolved':'Actual trans-binding, EC-specific necessity and MEIS1 mediation are not observed in this dataset',
         'experiment':'Conditional on detecting contact benefit: EC-specific Ptprm loss/rescue with WT CMs, reciprocal CM perturbation and MEIS1 rescue for the cell-cycle effect',
         'falsification':'Contact benefit persists after verified EC-specific loss with preserved EC viability/network support; or proposed MEIS1 rescue fails'},
        {'priority':3,'hypothesis':'Fibroblast matrix state gates the vascular-apposition effect',
         'existing_evidence':'Observed RNA matrix remodeling across libraries; local apposition-by-basement-program interactions are not consistently positive',
         'unresolved':'RNA does not measure deposited ECM, stiffness or architecture',
         'experiment':'Hold matrix constant in the first experiment; test controlled matrix composition/mechanics only if direct measurements indicate a modifier',
         'falsification':'No matrix-by-contact interaction under controlled EC density, oxygen and baseline CM state'}]
    save(pd.DataFrame(priorities),'experiment_priority_and_falsifiers.tsv')
    dump('selected_experiment.json',{'selected':'Functional discrimination of oxygen support versus EC contact across CM developmental states',
        'factors':['Oxygen availability: two measured, noninjurious conditions chosen by pilot','EC contact: direct coculture versus physically separated ECs with matched medium/EC number','CM initial developmental state: early versus later, baseline characterized'],
        'primary_functional_endpoints':['Mitochondrial respiration/substrate use','Calcium transient kinetics and contractile performance'],
        'secondary_endpoints':['Sarcomere organization','CM cycling/withdrawal and ploidy','RNA programs corresponding to the spatial readouts'],
        'essential_controls':['Donor-paired allocation; independent biological donors as inference units','Fixed matrix, EC density and culture duration','Measured oxygen; matched survival and baseline CM state','ROS/DNA-damage assessment to separate injury-associated arrest from maturation','CM-only and noncontact EC controls'],
        'pilot':'At least four independent donor preparations per initial CM state, with treatments paired within donor; use between-donor variance to power the confirmatory study. Wells are technical repeats.',
        'interpretation':'Culture oxygen manipulates delivery conditions, not in vivo capillary perfusion; verify the in vivo correlate separately with perfused-vessel labeling/hypoxia imaging',
        'molecular_followup':'EC-specific Ptprm loss/rescue only after a contact benefit is demonstrated at matched oxygen; reciprocal CM perturbation and downstream rescue distinguish mechanism',
        'new_animal_or_culture_experiments_performed':False})

def summary(pdf):
    m=read(OUT/'primary_matched_component_results.tsv',index_col=None)
    def result(s,o):
        r=m[m.stage.eq(s)&m.outcome.eq('core_Z__'+o)].iloc[0]
        return f"{s}: {r.mean_standardized_difference:+.3f} SD [{r.block_CI_low:+.3f}, {r.block_CI_high:+.3f}], {int(r.matched_pairs):,} pairs"
    textpage('Spatial reconstruction: what is now supported',[
        ('The paper’s descriptive direction is recovered','At P7, maturation-index values decrease with increasing graph distance from capillary ECs, proliferating ECs and pericytes under explicit reconstructed definitions. Strong correlations of distance-bin means coexist with modest cell-level correlations. This is a descriptive reconstruction, not recovery of the missing original graph, score or anatomical masks.'),
        ('Shared anatomy accounts for much of the broad pattern','At P7, approximately 71% of capillary-area/maturation-index covariance lies between 200µm regions. This is a covariance decomposition, not 71% of a causal effect or variance explained. Within-region and strict matched analyses then ask whether interface-level associations remain.'),
        ('Apposition relates to different components across stages','Strict matches balance local region, author CM state, capture, area, perimeter and capillary abundance. Core-RNA differences for lower cycling RNA: '+result('P7','withdrawal')+'. Metabolic-regulatory RNA: '+result('P14','metabolic')+'; '+result('P21','metabolic')+'. These small within-section contrasts motivate functional tests, not proof of a developmental switch.'),
        ('RNA ordering and inferred missing days','43,249 observed ventricular RNA nuclei have a reproducible developmental ordering in pooled and disjoint replicate series. P3/P5/P10 estimates are explicitly conditional interpolation of transcript-program expectations. They do not add observed hearts, sections, perfusion measurements or biological replication.'),
        ('Experiment selected','Prioritize a functional oxygen-by-contact comparison across CM developmental states, holding matrix and EC number constant. It distinguishes permissive supply from contact effects. PTPRM is a conditional molecular follow-up; the present data do not justify treating mean EC PTPRM RNA as the established driver.')
    ],'00_findings_and_scope',pdf)

def bg(ax,stage):
    from PIL import Image
    info=json.loads((ROOT/f'results/xenium_ventricular_mapping/{stage}/morphology_display.json').read_text());im=np.asarray(Image.open(ROOT/f'results/xenium_ventricular_mapping/{stage}/DAPI_overview.png'))
    ax.imshow(1-.4*im/255,cmap='gray',vmin=0,vmax=1,extent=info['extent_um'],origin='upper',rasterized=True)
    ax.set(xlim=(0,info['extent_um'][1]),ylim=(info['extent_um'][2],0),xlabel='X (µm)',ylabel='Y (µm)');ax.set_aspect('equal')

def anatomy_maps(pdf):
    fig,axs=plt.subplots(2,2,figsize=(15,14),layout='constrained')
    palette={'Compact-compatible dense myocardium':'#277c8e','Cavity-associated porous myocardium':'#d58b23','Transition/unclassified':'#aaa9a6'}
    for stage,ax in zip(STAGES,axs.ravel()):
        bg(ax,stage);c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');c=c[c.reconstructed_Figure4e_CM]
        for label,color in palette.items():
            z=c[c.morphology_region.eq(label)];ax.scatter(z.x_centroid,z.y_centroid,c=color,s=.6,lw=0,rasterized=True)
        with np.load(OUT/stage/'morphology_masks.npz') as b:
            res=float(b['resolution_um']);mask=b['envelope'];ax.contour(np.arange(mask.shape[1])*res,np.arange(mask.shape[0])*res,mask,levels=[.5],colors=['#742f72'],linewidths=.9)
        ax.set_title(f'{stage}: {len(c):,} author-label CMs in independent selector')
    fig.legend(handles=[Line2D([],[],ls='',marker='o',color=v,label=k) for k,v in palette.items()],loc='upper center',bbox_to_anchor=(.5,1.015),ncol=3,frameon=False,fontsize=9)
    write(fig,'01_independent_anatomy_masks',pdf,'Acquired DAPI and native cell-territory geometry. Purple: independent envelope. These are operational morphology classes, not recovered author LV/RV/trabecular source labels or anatomical ground truth.')

def reconstructed_curves(pdf):
    d=read(OUT/'paper_style_neighbor_level_means.tsv',index_col=None);smooth=read(OUT/'paper_style_R_LOESS_curves.tsv',index_col=None);a=read(OUT/'descriptive_proximity_associations.tsv',index_col=None)
    fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for j,source in enumerate(['cap.EC','p.EC','Peri']):
        for graph,color in [('paper_like6NN','#0072b2'),('tissue_Delaunay','#d58b23')]:
            take=d.stage.eq('P7')&d.cohort.eq('Figure4e_independent_mask')&d.source.eq(source)&d.outcome.eq('native_Z__maturation21')&d.graph.eq(graph)
            z=d[take];ss=smooth[smooth.stage.eq('P7')&smooth.cohort.eq('Figure4e_independent_mask')&smooth.source.eq(source)&smooth.outcome.eq('native_Z__maturation21')&smooth.graph.eq(graph)]
            axs[0,j].plot(ss.level,ss.LOESS_mean,c=color,label=graph);axs[0,j].scatter(z.level,z['mean'],c=color,s=22)
            axs[1,j].plot(z.level,z.CMs,'o-',c=color,ms=4)
        r=a[a.stage.eq('P7')&a.cohort.eq('Figure4e_independent_mask')&a.source.eq(source)&a.outcome.eq('native_Z__maturation21')&a.graph.eq('paper_like6NN')].iloc[0]
        axs[0,j].set(title=f'{source}: individual ρ={r.individual_Spearman:.3f}\nDistance-bin mean r={r.binned_mean_Pearson:.3f}',xlabel='Graph rings (≥10 pooled)',ylabel='Mean reconstructed 21-gene index')
        axs[1,j].set(yscale='log',xlabel='Graph rings (≥10 pooled)',ylabel='CMs contributing to each mean')
    axs[0,0].legend(fontsize=8,frameon=False);fig.suptitle('Reconstruct first: the P7 descriptive proximity pattern is positive',fontsize=16)
    write(fig,'02_P7_reconstructed_descriptive_relationship',pdf,'Native measured RNA with pooled per-gene standardization; R stats::loess reproduces the released plotting rule. Exact author score/mask/neighbor inputs are absent. Bin correlations are not independent-animal evidence.')

def component_maps():
    fig,axs=plt.subplots(2,4,figsize=(21,11),layout='constrained')
    for col,stage in enumerate(STAGES):
        c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');e=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz');c=c[c.reconstructed_Figure4e_CM]
        for ax in axs[:,col]:bg(ax,stage)
        im=axs[0,col].scatter(c.x_centroid,c.y_centroid,c=c.core_Z__metabolic,s=.7,cmap='viridis',vmin=0,vmax=4,lw=0,rasterized=True)
        fig.colorbar(im,ax=axs[0,col],shrink=.6,label='Core metabolic-regulatory RNA index')
        im=axs[1,col].scatter(c.x_centroid,c.y_centroid,c=e.loc[c.index,'author_capillary__apposition_fraction_2um'],s=.7,cmap='magma',vmin=0,vmax=.5,lw=0,rasterized=True)
        fig.colorbar(im,ax=axs[1,col],shrink=.6,label='CM boundary fraction within 2 µm of cap.EC')
        axs[0,col].set_title(stage+' · measured CM program');axs[1,col].set_title(stage+' · operational capillary apposition')
    write(fig,'03_all_stage_component_and_apposition_maps',None,'Observed stages only; common color scales. Metabolic-regulatory transcription is not measured respiration/perfusion. Sampled polygon apposition is not verified membrane protein contact.')

def anatomy_and_matches(pdf):
    fig,axs=plt.subplots(1,3,figsize=(18,6),layout='constrained')
    d=read(OUT/'shared_anatomy_covariance_partition.tsv',index_col=None);d=d[d.stage.eq('P7')&d.outcome.ne('native_Z__withdrawal')]
    labels={'native_Z__maturation21':'21-gene index',**{'native_Z__'+k:v for k,v in COMP.items()}}
    axs[0].barh(np.arange(len(d)),100*d.between_fraction_of_total,color='#799b9d');axs[0].set(yticks=np.arange(len(d)),yticklabels=[labels[k] for k in d.outcome],xlim=(0,100),xlabel='% total covariance between 200µm regions',title='P7: regional organization is substantial');axs[0].invert_yaxis()
    m=read(OUT/'matched_apposition_component_differences.tsv',index_col=None)
    for ax,component in zip(axs[1:],['withdrawal','metabolic']):
        for vi,view in enumerate(['native','core']):
            z=m[m.outcome.eq(view+'_Z__'+component)].set_index('stage').loc[STAGES]
            for j,stage in enumerate(STAGES):
                r=z.loc[stage];ax.errorbar(r.mean_standardized_difference,j+(vi-.5)*.14,xerr=[[r.mean_standardized_difference-r.block_CI_low],[r.block_CI_high-r.mean_standardized_difference]],fmt='o' if view=='core' else 's',mfc=COLORS[stage] if view=='core' else 'white',c=COLORS[stage],capsize=2)
        ax.axvline(0,c='gray',lw=.8);ax.set(yticks=np.arange(4),yticklabels=[s+f' ({int(m[m.stage.eq(s)].matched_pairs.iloc[0]):,} pairs)' for s in STAGES],xlabel='High versus zero apposition: standardized difference',title=COMP[component]);ax.invert_yaxis()
    fig.suptitle('After strict local matching, component-specific associations remain',fontsize=16)
    write(fig,'04_shared_anatomy_and_strict_matched_CM_results',pdf,'Squares: native RNA; circles: core RNA; 95% intervals from 400 µm blocks. Same 200µm region/CM state, ≤100µm apart; strict capture/shape and capillary-area calipers. Covariance partition is not causal attribution.')

def robustness(pdf):
    fig,axs=plt.subplots(2,2,figsize=(16,11),layout='constrained')
    d=read(OUT/'paired_contact_definition_coefficients.tsv.gz',index_col=None)
    for j,(key,label) in enumerate([('native_sampled_surface_fraction','Native surface fraction'),('refined_boundary__capillary_contacts','Refined contact count'),('native_polygon_1um__capillary_contacts','Native ≤1µm contact count')]):
        z=d[d.exposure.eq('capillary_apposition')&d.outcome.eq('core_Z__metabolic')&d.contact_definition.eq(key)].set_index('stage').loc[STAGES]
        axs[0,0].errorbar(np.arange(4)+(j-1)*.1,z.standardized_coefficient,yerr=[z.standardized_coefficient-z.block_CI_low,z.block_CI_high-z.standardized_coefficient],fmt='o-',ms=4,capsize=2,label=label)
    d=read(OUT/'anatomy_conditioned_coefficients.tsv.gz',index_col=None)
    for j,(key,label) in enumerate([('author_capillary','Public cap.EC'),('author_capillary_plus_HighMyo','Public cap.EC + High.Myo'),('QC_capillary','Strict nuclear-supported capillary')]):
        z=d[d.cohort.eq('paper_label_mask')&d.capillary_source.eq(key)&d.adjustment.eq('within_region_capture')&d.exposure.eq('capillary_apposition')&d.outcome.eq('core_Z__metabolic')].set_index('stage').loc[STAGES]
        axs[0,1].errorbar(np.arange(4)+(j-1)*.1,z.standardized_coefficient,yerr=[z.standardized_coefficient-z.block_CI_low,z.block_CI_high-z.standardized_coefficient],fmt='o-',ms=4,capsize=2,label=label)
    for ax in axs[0]:ax.axhline(0,c='gray',lw=.8);ax.set(xticks=np.arange(4),xticklabels=STAGES,ylabel='Adjusted apposition–metabolic coefficient (SD)');ax.legend(fontsize=8,frameon=False)
    axs[0,0].set_title('Geometry sensitivity: same CMs, different exposure definitions');axs[0,1].set_title('Identity sensitivity: strict capillary source attenuates effects')
    d=read(OUT/'nuclear_cycling_specificity_coefficients.tsv',index_col=None)
    for j,(key,label) in enumerate([('within50um','QC-EC nuclear cycling within 50µm'),('annulus30to50um','QC-EC nuclear cycling in 30–50µm shell')]):
        z=d[d.geometry.eq(key)&d.outcome.eq('CM_nuclear_cycling')].set_index('stage').loc[STAGES]
        axs[1,0].errorbar(np.arange(4)+(j-.5)*.15,z.standardized_coefficient,yerr=[z.standardized_coefficient-z.block_CI_low,z.block_CI_high-z.standardized_coefficient],fmt='o-',ms=4,capsize=2,label=label)
    axs[1,0].axhline(0,c='gray',lw=.8);axs[1,0].set(xticks=np.arange(4),xticklabels=STAGES,ylabel='EC cycling–CM cycling coefficient (SD)',title='Early shared cycling survives the nuclear/QC source check');axs[1,0].legend(fontsize=8,frameon=False)
    z=read(OUT/'nuclear_cycling_specificity_coverage.tsv',index_col=None);z=z[z.geometry.eq('within50um')].set_index('stage').loc[STAGES]
    axs[1,1].plot(np.arange(4),100*z.CM_fraction_two_cycling_genes_detected,'o-',label='CM nuclear profiles');axs[1,1].plot(np.arange(4),100*z.EC_fraction_two_cycling_genes_detected,'s-',label='QC EC nuclear profiles')
    axs[1,1].set(xticks=np.arange(4),xticklabels=STAGES,ylabel='% with ≥2 measured cycling genes',title='Late low cycling limits interpretation of withdrawal effects');axs[1,1].legend(fontsize=8,frameon=False)
    write(fig,'05_geometry_identity_and_nuclear_sensitivity',pdf,'All intervals are within-section 400µm block intervals, not animal replication. Native/refined metrics have different units before standardization. Nuclear-marker coherence is not proof of completed cell division; 30–50µm shell is a centroid-distance sensitivity.')

def coefficient_gallery():
    d=read(OUT/'anatomy_conditioned_coefficients.tsv.gz',index_col=None);d=d[d.cohort.eq('paper_label_mask')&d.capillary_source.eq('author_capillary')&d.adjustment.eq('within_region_capture')]
    exposures=['capillary_area','capillary_apposition','pericyte_context','FB_basement_program','FB_fibrillar_program','EC_cycling','EC_Ptprm']
    fig,axs=plt.subplots(2,2,figsize=(16,10),layout='constrained')
    for ax,(component,label) in zip(axs.ravel(),COMP.items()):
        z=d[d.outcome.eq('core_Z__'+component)];matrix=z.pivot(index='stage',columns='exposure',values='standardized_coefficient').reindex(index=STAGES,columns=exposures)
        im=ax.imshow(matrix,cmap='RdBu_r',vmin=-.15,vmax=.15,aspect='auto');ax.set(xticks=np.arange(len(exposures)),xticklabels=[s.replace('_',' ') for s in exposures],yticks=np.arange(4),yticklabels=STAGES,title=label);ax.tick_params(axis='x',rotation=40)
        for i,stage in enumerate(STAGES):
            for j,ex in enumerate(exposures):
                r=z[z.stage.eq(stage)&z.exposure.eq(ex)].iloc[0]
                ax.text(j,i,f'{r.standardized_coefficient:.2f}'+('•' if r.within_section_FDR_BH<=.05 and (r.block_CI_low>0 or r.block_CI_high<0) else ''),ha='center',va='center',fontsize=8)
        fig.colorbar(im,ax=ax,shrink=.7,label='Conditional standardized coefficient')
    write(fig,'06_joint_context_coefficients',None,'Core-RNA outcomes; within 200µm regions plus capture/geometry controls. Dot: block interval excludes zero and within-section BH ≤.05. Secondary nuclear/QC tests are required: late core-only cycling associations are not robust in those tests. No causal effect or biological replication is asserted.')

def RNA_ordering(pdf):
    d=read(OUT/'RNA_cell_developmental_ordering.tsv.gz');z=d[d.ordering.eq('pooled')]
    with np.load(OUT/'pooled_RNA_ordering_model.npz') as b:pc=b['PCs'];ids=b['cell_ids'].astype(str)
    z=z.loc[ids];fig,axs=plt.subplots(2,2,figsize=(15,11),layout='constrained')
    for stage,color in COLORS.items():
        m=z.age.eq(stage).to_numpy();axs[0,0].scatter(pc[m,0],pc[m,1],s=.4,c=color,lw=0,label=stage,rasterized=True)
    axs[0,0].set(xlabel='Outcome-withheld RNA PC1',ylabel='PC2',title='43,249 observed ventricular nuclei');axs[0,0].legend(markerscale=6,frameon=False)
    im=axs[0,1].scatter(pc[:,0],pc[:,1],c=z.graph_geodesic_order,s=.4,cmap='viridis',lw=0,rasterized=True);fig.colorbar(im,ax=axs[0,1],label='Graph distance ordering (not calendar age)');axs[0,1].set(xlabel='PC1',ylabel='PC2',title='P0-rooted ordering; no intermediate cells generated')
    for j,name in enumerate(['pooled','RNA_replicate1','RNA_replicate2']):
        v=d[d.ordering.eq(name)].groupby('age').graph_geodesic_order.median().reindex(STAGES);axs[1,0].plot(np.arange(4),v,'o-',label=name)
    axs[1,0].set(xticks=np.arange(4),xticklabels=STAGES,ylabel='Median ordering within each independently fitted series',title='Observed order reproduces in disjoint RNA replicate series');axs[1,0].legend(fontsize=8,frameon=False)
    edge=read(OUT/'RNA_cross_age_graph_edges.tsv',index_col=None);edge=edge[edge.ordering.eq('pooled')].pivot(index='from_stage',columns='to_stage',values='fraction_of_from_stage_edges').reindex(index=STAGES,columns=STAGES)
    show=edge.to_numpy(copy=True);np.fill_diagonal(show,np.nan);im=axs[1,1].imshow(100*show,cmap='Blues',vmin=0,vmax=15)
    for i in range(4):
        for j in range(4):axs[1,1].text(j,i,'—' if i==j else f'{100*edge.iloc[i,j]:.2f}%',ha='center',va='center',fontsize=9)
    axs[1,1].set(xticks=np.arange(4),xticklabels=STAGES,yticks=np.arange(4),yticklabels=STAGES,xlabel='Neighbor age',ylabel='Focal age',title='Early age gaps are sparsely bridged');fig.colorbar(im,ax=axs[1,1],shrink=.7,label='% directed neighbor edges')
    write(fig,'07_observed_RNA_ordering_and_gap_support',pdf,'Features exclude all validation modules and cycling genes; PCA training is library-balanced. Pooled ordering removes within-age library mean offsets without forcing age integration. Root/age association is not a clock, lineage proof or dense observation of missing intervals; depth remains associated with ordering.')

def inferred_programs(pdf):
    observed=read(OUT/'RNA_observed_library_programs.tsv',index_col=None);inferred=read(OUT/'inferred_P3_P5_P10_programs.tsv',index_col=None)
    from scipy.interpolate import PchipInterpolator
    programs=[('Ventricular_CM','CM_metabolic','CM metabolic-regulatory RNA'),('Ventricular_CM','CM_cycling','CM cycling RNA'),('Ventricular_CM','CM_respiratory_chain_RNA_extended','CM respiratory-chain RNA (not respiration)'),
        ('QC_Endothelial','EC_cycling','QC EC cycling RNA'),('QC_reassessed_Fibroblast','FB_basement_matrix','Reassessed FB basement-matrix RNA'),('QC_reassessed_Fibroblast','FB_remodeling','Reassessed FB remodeling RNA')]
    fig,axs=plt.subplots(2,3,figsize=(17,10),layout='constrained');days=np.array([0,7,14,21]);grid=np.linspace(0,21,150)
    for ax,(scope,program,label) in zip(axs.ravel(),programs):
        z=observed[observed.scope.eq(scope)&observed.program.eq(program)]
        for day in days:
            vals=z.loc[z.day.eq(day),'pseudobulk_logCPM_mean'].to_numpy();ax.scatter(day+np.linspace(-.25,.25,len(vals)),vals,c='black',s=25,zorder=4)
        mean=z.groupby('day').pseudobulk_logCPM_mean.mean().reindex(days)
        ax.plot(grid,np.interp(grid,days,mean),c='#377eb8',ls='--',lw=1,label='Linear endpoint interpolation')
        ax.plot(grid,PchipInterpolator(days,mean)(grid),c='#d58b23',lw=1.5,label='PCHIP endpoint interpolation')
        v=inferred[inferred.scope.eq(scope)&inferred.program.eq(program)]
        center=(v.PCHIP_estimate+v.linear_estimate)/2
        ax.errorbar(v.inferred_day,center,yerr=[center-v.conditional_library_and_model_low,v.conditional_library_and_model_high-center],fmt='D',mfc='white',mec='#ae4968',ecolor='#ae4968',ms=5,capsize=3,label='P3/P5/P10 conditional estimate')
        ax.set(xticks=[0,3,5,7,10,14,21],xlabel='Postnatal day (inferred points are open diamonds)',ylabel='Pseudobulk mean log1p(full-RNA CPM)',title=label)
    axs[0,0].legend(fontsize=7,frameon=False)
    write(fig,'08_conditional_P3_P5_P10_program_estimates',pdf,'Black points: two observed RNA libraries/age. Open diamonds: model-based expectations, no new samples. Intervals combine endpoint-library bootstrap and two interpolation choices but cannot measure unknown curvature. Respiratory-chain and other RNA panels need not increase together and do not establish function.')

def experiment_page(pdf):
    textpage('Roadmap completed to a testable experiment',[
        ('Hypothesis selected from these results','The phenotype associated with vascular apposition is developmentally contingent: P7 is a candidate cell-cycle-withdrawal window, whereas P14/P21 show metabolic-regulatory associations. This hypothesis follows component-specific, locally matched observations. It is not a demonstrated causal timing switch, and late low cycling limits detection.'),
        ('First experiment: oxygen support versus contact','Use a factorial comparison of two measured noninjurious oxygen conditions, direct versus physically separated EC contact, and early versus later initial CM developmental states. Match EC number, matrix, medium, culture duration, baseline CM state and viability. Allocate treatments within independent donor preparations; wells are technical repeats. Pilot at least four donors/state and power confirmation using donor-level variance.'),
        ('Readouts that decide among mechanisms','Primary: mitochondrial respiration/substrate use and calcium/contractile function. Secondary: sarcomere organization, cycling/withdrawal, ploidy and the corresponding RNA programs. Measure oxygen and assess ROS/DNA damage and survival so injury-associated arrest is not called maturation. An in vivo correlate requires perfused-vessel/hypoxia imaging; culture oxygen is not capillary perfusion.'),
        ('Decision rules','An oxygen effect without additional contact benefit prioritizes permissive support. A contact benefit at matched oxygen prioritizes an instructive/interface mechanism. Then test EC-specific Ptprm loss/rescue with WT CMs, reciprocal CM perturbation and MEIS1 rescue for a cell-cycle effect. If controlled contact/oxygen effects do not reproduce, reconsider shared regional drivers and other niche cues.'),
        ('Matrix and remaining reconstruction work','Keep matrix constant first: current apposition-by-FB basement-program interactions are not consistently positive, and RNA is not deposited matrix or stiffness. Exact original Figure 4e replication still requires the unpublished source masks, score-bearing RDS and neighbor table/build settings. Their absence is documented; these derived reconstructions and inferred days do not replace those observations.')
    ],'09_selected_functional_experiment',pdf)

def main():
    selections_and_decision()
    with PdfPages(OUT/'cardiac_spatial_reconstruction_report.pdf') as pdf:
        summary(pdf);anatomy_maps(pdf);reconstructed_curves(pdf);anatomy_and_matches(pdf);robustness(pdf);RNA_ordering(pdf);inferred_programs(pdf);experiment_page(pdf)
    component_maps();coefficient_gallery()

if __name__=='__main__':main()
