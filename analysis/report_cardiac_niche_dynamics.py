"""Figures and a paper-aware, evidence-bounded cardiac niche report."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-dynamics-mpl')
from revisit_vascular_coordination import ROOT,OUT,UP,STAGES,read,save,dump,CFG
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import textwrap

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
PAL={'Cardiomyocyte':'#d56c82','Endothelial':'#377eb8','Fibroblast':'#259c65','Pericyte':'#916bb1','Smooth muscle':'#dc9a27','Immune':'#cc3939','Epicardial':'#ab9e12'}
AGES=dict(zip(STAGES,['#0072b2','#e69f00','#009e73','#cc79a7']))

def write(fig,name,pdf,foot=None):
    if foot:fig.text(.03,-.08,textwrap.fill(foot,180),fontsize=8,color='#444444',va='top')
    fig.savefig(OUT/(name+'.png'),dpi=150,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def textpage(pdf,title,paragraphs,name):
    fig=plt.figure(figsize=(12,9));fig.text(.06,.94,title,fontsize=19,weight='bold')
    y=.86
    for heading,body in paragraphs:
        fig.text(.06,y,heading,fontsize=12,weight='bold',va='top');y-=.035
        body=textwrap.fill(body,118)
        fig.text(.06,y,body,fontsize=10.5,va='top',linespacing=1.5)
        y-=.026*(body.count('\n')+1)+.032
    write(fig,name,pdf)

def review_table():
    coverage=[]
    for stage in STAGES:
        q=read(OUT/stage/'measured_vascular_states.tsv.gz')
        for scope,mask in [('all_broad_EC',q.EC_broad),('primary_EC',q.EC_primary)]:
            for name,n in q.loc[mask,'EC_nuclear_type'].fillna('Unresolved').value_counts().items():
                coverage.append({'stage':stage,'source_population':scope,'nuclear_subtype':name,'profiles':n,'fraction':n/max(mask.sum(),1),
                    'scope':'Subtype nonassignment or no eligible profiles is not absence of that cell type in tissue'})
    save(pd.DataFrame(coverage),'EC_source_nuclear_subtype_coverage.tsv')
    rows=[
      ('Lu2021 Fig1-4, Methods','Late-fetal/neonatal pre-existing coronary plexus expands; late endocardium does not measurably supply coronary ECs in tested tracing windows','Separate coronary-compatible, endocardial and lymphatic states; no computational endocardial-to-coronary ancestry','Measured states and source audits; lineage inference prohibited'),
      ('Lu2021 Fig3','Flt4/Dll4-enriched inner plexus; regional Hypoxyprobe and EdU evidence; Aplnr/Emcn outer gradient','Growth, regional oxygenation, proliferation and venous identity are different axes','Flt4+Dll4 core co-detection proxy; Emcn missing; neither hypoxia nor lumen morphology claimed'),
      ('Lu2021 Fig6','Dll4 deletion increases plexus expansion; NICD activation suppresses growth','DLL4/NOTCH expression is not a monotonic positive sprouting measure','Growth proxy, Notch-response proxy and EC-EC expression compatibility stored separately'),
      ('Lu2021 Discussion','Rapid perinatal growth declines approximately P7-P14','P0/P7/P14/P21 cannot resolve P1-P6 sequence or fetal ancestry','Cross-sectional interpretation; per-library RNA cycling and sparse sprout-marker audit'),
      ('Wang2026 Fig2, CellCharter Methods','scVI expression latent, neighbor aggregation, ClusterAutoK, final K12','Prior composition-only K2 is a different feature model, not a reproduction or exhaustive anatomy','Richer outcome-withheld non-CM observed-gene state distribution added; no K12 replication claimed'),
      ('Wang2026 Fig4, Maturation Methods','21-gene index; body describes standardized mean, Methods specify AddModuleScore with 100 matched controls per gene','Exact pool/binning/assay is insufficiently specified for a measured403-only reproduction','Native/core means, prior CM-only controls, feasible all-cell8-bin20-control sensitivity; report ambiguity without alleging failure'),
      ('Wang2026 Fig4d, Proximity Methods','P7 capillary/proliferating EC and pericyte proximity associated with higher CM maturity; graph rings0-10','Test positive contact/integration hypotheses; do not assume proliferation implies delayed maturation','Pruned Delaunay and refined-contact rings, continuous distances, unreachable censor flags, author-label sensitivity'),
      ('Wang2026 Fig4, ExtendedData7','PTPRM boundary IF, MEIS1 association, CM-directed Ptprm perturbation affects cycling/maturity','CM receiver and CM-CM contact may contribute; EC-specific necessity is not established','EC nuclear sender, CM sender, receiver conditioning, native/refined geometry and source-label sensitivity'),
      ('Wang2026 Communication Methods','COMMOT0.0.3, CellPhoneDBv4, heteromeric complexes, paracrine200um and ECM/contact50um','Expression products alone do not reproduce transport-based signaling inference','Components and graph-specific compatibility reported; no COMMOT implementation or pathway activation claimed'),
      ('Wang2026 Segmentation Methods','Multimodal vendor segmentation includes boundary/interior stains;10FOV manual comparison; overlap>=90% and TP/(TP+FP+FN)','Refinement is not proven superior; CM/EC mixed programs may be biological or segmentation sensitive','Native polygon0.5/1/2um graphs, refined/core sensitivity and cap.EC.High.Myo source audit'),
      ('Wormhole benchmark, this release','Distribution transformer trained on3200 niches with800 held-out niches; empirical OT distance checked before phenotype use','Adopt only if preserved held-out distances and useful relative to simple means','Metric gate fails across ages; exact empirical OT and descriptive means retained; no genuine longitudinal reconstruction')]
    save(pd.DataFrame(rows,columns=['source','published_evidence_or_method','analysis_implication','implemented_or_boundary']),'paper_review_and_pipeline_changes.tsv')

def contexts(pdf):
    d=read(OUT/'coarse_context_composition_explanation.tsv',index_col=None)
    fig,axs=plt.subplots(1,2,figsize=(15,7),layout='constrained')
    p=d.pivot_table(index=['stage','computational_context'],columns='neighbor_type',values='mean_resolved_neighbor_fraction')
    p=p.reindex(sorted(p.index,key=lambda x:(STAGES.index(x[0]),x[1])))
    left=np.zeros(len(p));y=np.arange(len(p))
    for cell,color in PAL.items():axs[0].barh(y,p[cell],left=left,color=color,label=cell);left+=p[cell]
    axs[0].set(yticks=y,yticklabels=[a+' · '+(b.split(':')[0] if ':' in b else 'Unresolved') for a,b in p.index],xlabel='Mean direct-neighbor fraction',title='N1: FB/SMC enriched · N2: CM/EC enriched');axs[0].invert_yaxis()
    axs[0].legend(fontsize=8,bbox_to_anchor=(0,-.14),loc='upper left',ncol=3,frameon=False)
    a=read(OUT/'native_refined_architecture_enrichment.tsv',index_col=None)
    pairs=[('Smooth muscle','Smooth muscle'),('Epicardial','Epicardial'),('Immune','Immune'),('Fibroblast','Immune')]
    for j,(ca,cb) in enumerate(pairs):
        for k,geometry in enumerate(['native_boundary_1um','refined_contact']):
            z=a[a.type_a.eq(ca)&a.type_b.eq(cb)&a.geometry.eq(geometry)].set_index('stage').loc[STAGES]
            axs[1].plot(np.arange(4)+(j-1.5)*.12,z.log2_observed_over_local_null,marker='o' if k==0 else 'x',ls='-' if k==0 else '--',color=plt.get_cmap('tab10')(j),label=ca+'–'+cb if k==0 else None)
    axs[1].axhline(0,c='gray',lw=.8);axs[1].set(xticks=np.arange(4),xticklabels=STAGES,ylabel='log2(observed contacts / local-shuffle expectation)',title='Spatial concentration survives geometry changes')
    axs[1].legend(fontsize=8,frameon=False);fig.suptitle('Broad contexts and contact concentration',fontsize=17)
    write(fig,'01_contexts_and_concentration',pdf,'Solid circles: native polygons within1µm; dashed crosses: refined shared boundary. Null shuffles labels within200µm blocks. Association, not signaling.')

def RNA(pdf):
    d=read(OUT/'RNA_program_library_profiles.tsv',index_col=None);d=d[d.scope.eq('EC_primary')]
    fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    names=['Cycling_RNA','Sprouting_6g','Plexus_growth_2021','Notch_response_RNA','Hypoxia_compatible_RNA']
    for ax,name in zip(axs.ravel(),names):
        z=d[d.program.eq(name)]
        for j,stage in enumerate(STAGES):
            vals=z.loc[z.stage.eq(stage),'pseudobulk_mean_log1p_fullRNA_CPM'].to_numpy()
            ax.scatter(j+np.linspace(-.09,.09,len(vals)),vals,c=AGES[stage],s=45)
            if len(vals):ax.plot([j-.15,j+.15],[np.mean(vals)]*2,c='black',lw=1)
        ax.set(xticks=np.arange(4),xticklabels=STAGES,ylabel='Pseudobulk mean log1p(full-RNA CPM)',title=name.replace('_',' '))
    co=read(OUT/'RNA_vascular_module_coherence.tsv',index_col=None)
    z=co.pivot_table(index=['stage','library'],columns='program',values='median_leave_one_gene_out_partial_Spearman').reindex(columns=names)
    z=z.reindex(sorted(z.index,key=lambda a:(STAGES.index(a[0]),a[1])))
    im=axs[1,2].imshow(z.to_numpy(),vmin=-.5,vmax=.8,cmap='coolwarm',aspect='auto');axs[1,2].set(yticks=np.arange(len(z)),yticklabels=[s+' '+l.rsplit('_',1)[-1] for s,l in z.index],xticks=np.arange(len(names)),xticklabels=['Cycle','Sprout','Plexus','Notch','Hypoxia'],title='Depth/family-adjusted module coherence')
    for i in range(len(z)):
        for j in range(len(names)):
            val=z.iloc[i,j];axs[1,2].text(j,i,'NA' if pd.isna(val) else f'{val:.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=axs[1,2],shrink=.6,label='Median leave-one-gene-out partial Spearman')
    fig.suptitle('Whole-transcriptome snRNA-seq: cycling is better supported than one sprouting module',fontsize=16)
    write(fig,'02_RNA_vascular_programs',pdf,'Two libraries/age, corrected counts. Six-gene sprouting panel is sparse; scores do not establish sprout morphology. NA: <2 sufficiently detected genes.')

def maps(pdf):
    fig,axs=plt.subplots(2,2,figsize=(15,14),layout='constrained')
    colors={'Other eligible EC':'#a8bdd0','Endocardial-compatible':'#f0ae42','Growth-compatible':'#bd3734'}
    for stage,ax in zip(STAGES,axs.ravel()):
        f=read(UP/stage/'atlas_annotations.tsv.gz');q=read(OUT/stage/'measured_vascular_states.tsv.gz')
        ax.scatter(f.x_centroid,f.y_centroid,c='#e8e8e8',s=.1,linewidths=0,rasterized=True)
        masks=[q.EC_primary,q.EC_primary&q.EC_nuclear_type.eq('Endocardial-like'),q.EC_growth_compatible]
        for mask,(label,color) in zip(masks,colors.items()):
            z=q[mask];ax.scatter(z.x_centroid,z.y_centroid,c=color,s=2 if label=='Growth-compatible' else .8,linewidths=0,rasterized=True)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='X (µm)',ylabel='Y (µm)',title=f'{stage}: {int(q.EC_growth_compatible.sum()):,} growth-compatible / {int(q.EC_primary.sum()):,} eligible ECs')
    fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,label=l) for l,c in colors.items()],loc='upper center',bbox_to_anchor=(.5,1.02),ncol=3,frameon=False)
    write(fig,'03_measured_vascular_state_maps',pdf,'Growth-compatible = primary EC, not assigned endocardial/lymphatic identity, Flt4 AND Dll4 core co-detection. Unresolved subtypes remain possible; no verified sprouts or oxygenation.')
    fig,axs=plt.subplots(2,2,figsize=(15,14),layout='constrained')
    for row,stage in enumerate(['P0','P7']):
        f=read(UP/stage/'atlas_annotations.tsv.gz');q=read(OUT/stage/'measured_vascular_states.tsv.gz');cm=f[f.ventricular_CM_strict]
        for ax in axs[row]:ax.scatter(f.x_centroid,f.y_centroid,c='#eeeeee',s=.1,linewidths=0,rasterized=True)
        im=axs[row,0].scatter(cm.x_centroid,cm.y_centroid,c=cm.CM_maturation_core,cmap='viridis',vmin=1.5,vmax=3.7,s=.8,linewidths=0,rasterized=True)
        z=q[q.EC_growth_compatible];axs[row,0].scatter(z.x_centroid,z.y_centroid,c='#e63e23',s=3,linewidths=0,rasterized=True)
        fig.colorbar(im,ax=axs[row,0],shrink=.5,label='CM core21 mean log1p(panel CP10k)')
        z=q[q.EC_primary];im=axs[row,1].scatter(z.x_centroid,z.y_centroid,c=z.nuclear__Ptprm,cmap='magma',vmin=0,vmax=6,s=2,linewidths=0,rasterized=True)
        fig.colorbar(im,ax=axs[row,1],shrink=.5,label='EC nuclear Ptprm log1p(panel CP10k)')
        axs[row,0].set_title(stage+' · CM maturation and growth-compatible ECs (red)');axs[row,1].set_title(stage+' · PTPRM expression in eligible EC nuclei')
        for ax in axs[row]:ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='X (µm)',ylabel='Y (µm)')
    write(fig,'04_early_CM_vascular_maps',pdf,'Spatial overlap is descriptive. Adjacent sections/ages are different specimens; these maps do not measure EC-specific causal signals or local oxygen delivery.')

def rings(pdf):
    d=read(OUT/'vascular_ring_and_model_summary.tsv',index_col=None)
    fig,axs=plt.subplots(2,4,figsize=(18,9),layout='constrained')
    types={'Author_capillary':'#377eb8','Capillary_compatible':'#00522f','Endocardial_compatible':'#e69f00','Growth_compatible':'#b53336'}
    xnames=['1','2','3','4','5-6','7-10']
    for col,stage in enumerate(STAGES):
        for row,outcome in enumerate(['maturity_native_21','maturity_core_21']):
            z=d[d.stage.eq(stage)&d.graph.eq('Delaunay')&d.outcome.eq(outcome)]
            for name,color in types.items():
                v=z[z.EC_definition.eq(name)].copy();v['ring']=v.ring.astype(str)
                v=v.set_index('ring').reindex(xnames);axs[row,col].plot(np.arange(6),v.mean_maturity,c=color,marker='o',ms=3,label=name.replace('_',' '))
            axs[row,col].set(xticks=np.arange(6),xticklabels=xnames,xlabel='Nearest source: graph rings',title=stage+' · '+('native' if row==0 else 'core')+' score',ylabel='Mean21-gene maturation')
    axs[0,0].legend(fontsize=7,frameon=False)
    fig.suptitle('Proximity–maturation relationships depend on EC definition and RNA view',fontsize=16)
    write(fig,'05_graph_ring_maturation',pdf,'Pruned Delaunay graph, reached CMs only; censored/unreachable profiles are tabulated separately. Descriptive curves, not exact author graph reproduction.')

def effects(pdf):
    c=read(OUT/'vascular_maturation_coefficients.tsv.gz',index_col=None)
    fig,axs=plt.subplots(1,3,figsize=(16,6),layout='constrained')
    for ax,ex in zip(axs,['Plexus_growth_proxy','Cycling','EC_sender_Ptprm']):
        for j,stage in enumerate(STAGES):
            z=c[c.stage.eq(stage)&c.geometry.eq('proximity50')&c.adjustment.eq('receiver_conditioned')&c.outcome.eq('maturity_core_21')&c.exposure.eq(ex)].iloc[0]
            ax.errorbar(z.standardized_coefficient,j,xerr=[[z.standardized_coefficient-z.block_CI_low],[z.block_CI_high-z.standardized_coefficient]],fmt='o',c=AGES[stage],capsize=3)
        ax.axvline(0,c='gray',lw=.8);ax.set(yticks=np.arange(4),yticklabels=STAGES,xlabel='Standardized conditional coefficient (95% block CI)',title=ex.replace('_',' '));ax.invert_yaxis()
    fig.suptitle('Initial joint50µm models: small conditional associations',fontsize=16)
    write(fig,'06_conditional_growth_and_contact_associations',pdf,'Joint EC programs/families, CM sender/receiver, focal capture/area, other-cell fractions, density, spatial splines and distance proxies. No causal attribution.')
    d=read(OUT/'paired_geometry_coefficients.tsv.gz',index_col=None)
    fig,axs=plt.subplots(1,3,figsize=(16,6),layout='constrained')
    for ax,ex in zip(axs,['Plexus_growth_proxy','Cycling','EC_sender_Ptprm']):
        for j,geometry in enumerate(['proximity50','refined_contact','native_boundary']):
            z=d[d.stage.eq('P7')&d.geometry.eq(geometry)&d.exposure.eq(ex)&d.outcome.eq('maturity_core_21')].iloc[0]
            ax.errorbar(z.standardized_coefficient,j,xerr=[[max(0,z.standardized_coefficient-z.block_CI_low)],[max(0,z.block_CI_high-z.standardized_coefficient)]],fmt='o',c=plt.get_cmap('tab10')(j),capsize=3)
        ax.axvline(0,c='gray',lw=.8);ax.set(yticks=np.arange(3),yticklabels=['50µm radius','Refined boundary','Native polygon≤1µm'],xlabel='Standardized coefficient (95% block CI)',title=ex.replace('_',' '));ax.invert_yaxis()
    n=int(d[d.stage.eq('P7')].CMs.iloc[0]);fig.suptitle(f'P7: same {n:,} CMs across geometries, with source capture/coverage adjustment',fontsize=16)
    write(fig,'07_paired_geometry_depth_sensitivity',pdf,'Same focal CMs; neighboring EC nuclear/core depth and source-selection coverage included. Different graphs define different exposure variables.')

def source_audit(pdf):
    d=read(OUT/'author_EC_state_source_selection_audit.tsv',index_col=None)
    fig,axs=plt.subplots(1,2,figsize=(15,6),layout='constrained')
    z=d[d.author_state.eq('cap.EC.High.Myo')].set_index('stage').loc[STAGES];x=np.arange(4)
    for j,(field,label,color) in enumerate([('published_profiles','Public High.Myo state','#a8a8a8'),('our_broad_EC','Our broad EC','#6d9ec7'),('our_primary_EC_source','Our primary EC source','#176198')]):
        axs[0].bar(x+(j-1)*.25,z[field],width=.24,label=label,color=color)
    axs[0].set(xticks=x,xticklabels=STAGES,ylabel='Spatial profiles',title='Strict EC source selection changes the population');axs[0].legend(fontsize=8,frameon=False)
    c=read(OUT/'P7_Ptprm_source_sensitivity.tsv.gz',index_col=None);c=c[c.exposure.eq('EC_sender_Ptprm')&c.outcome.eq('maturity_core_21')]
    for i,r in enumerate(c.itertuples()):axs[1].errorbar(r.standardized_coefficient,i,xerr=[[r.standardized_coefficient-r.block_CI_low],[r.block_CI_high-r.standardized_coefficient]],fmt='o',capsize=3,color='#377eb8')
    axs[1].axvline(0,c='gray',lw=.8);axs[1].set(yticks=np.arange(len(c)),yticklabels=[x.replace('_',' ') for x in c.source_definition],xlabel='P7 sender PTPRM coefficient (95% block CI)',title='Same CMs: alternate sender definitions');axs[1].invert_yaxis()
    write(fig,'08_source_selection_sensitivity',pdf,'Public High.Myo may contain biological CM–EC programs and mixed boundaries. Exclusion is not proof of doublets; inclusion is not independent identity validation.')

def wormhole(pdf):
    b=read(OUT/'Wormhole_reserved_distance_benchmark.tsv',index_col=None)
    fig,axs=plt.subplots(1,4,figsize=(18,5),layout='constrained')
    metric=read(OUT/'Wormhole_metric_assessment.tsv',index_col=None).set_index('stage')
    for ax,stage in zip(axs,STAGES):
        z=b[b.stage.eq(stage)];ax.scatter(z.Sinkhorn_divergence,z.Wormhole_mean_squared_embedding_distance,s=18,c=AGES[stage]);lim=max(z.Sinkhorn_divergence.max(),z.Wormhole_mean_squared_embedding_distance.max());ax.plot([0,lim],[0,lim],c='gray',ls='--',lw=.8)
        r=metric.loc[stage];ax.set(xlabel='Held-out Sinkhorn divergence',ylabel='Wormhole mean squared embedding distance',title=f'{stage}: ρ={r.Wormhole_vs_Sinkhorn_Spearman:.2f}\nMean-feature baseline ρ={r.mean_vs_Sinkhorn_Spearman:.2f}')
    fig.suptitle('Wormhole trained; held-out transport fidelity does not justify general adoption',fontsize=16)
    write(fig,'09_Wormhole_distance_benchmark',pdf,'64 reserved within-age pairs/stage, non-CM gene-state clouds≤32 cells; actual encoder/decoder1200 steps. CPU benchmark, not default paper architecture.')
    p=read(OUT/'Wormhole_maturation_prediction.tsv',index_col=None)
    fig,axs=plt.subplots(1,2,figsize=(14,6),layout='constrained')
    pivot=p.pivot(index='stage',columns='predictors',values='heldout_R2').loc[STAGES]
    for j,key in enumerate(['baseline','mean_state','Wormhole']):axs[0].bar(np.arange(4)+(j-1)*.25,pivot[key],width=.24,label=key)
    axs[0].set(xticks=np.arange(4),xticklabels=STAGES,ylabel='Held-out maturation R²',title='Fixed spatial holdout: little extra predictive value');axs[0].legend(frameon=False)
    for j,key in enumerate(['mean_state','Wormhole']):axs[1].bar(np.arange(4)+(j-.5)*.3,pivot[key]-pivot.baseline,width=.28,label=key)
    axs[1].axhline(0,c='gray',lw=.8);axs[1].set(xticks=np.arange(4),xticklabels=STAGES,ylabel='Δ held-out R² relative to baseline',title='State embeddings versus technical/anatomy/composition');axs[1].legend(frameon=False)
    write(fig,'10_Wormhole_spatial_prediction',pdf,'150µm train/test buffer; sources in test clouds excluded from PCA fitting; focal maturation genes excluded from non-CM expression and normalization.')
    c=read(OUT/'cross_age_state_correspondence.tsv',index_col=None)
    fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    for ax,(a,nxt) in zip(axs,zip(STAGES[:-1],STAGES[1:])):
        z=c[c.from_stage.eq(a)&c.to_stage.eq(nxt)];m=z.pivot(index='from_state',columns='to_state',values='exact_medoid_W2_squared')
        im=ax.imshow(m,cmap='viridis_r');ax.set(xlabel=nxt+' exploratory state',ylabel=a+' exploratory state',title=a+' → '+nxt+' exact medoid OT')
        for r in z[z.nearest_state_in_next_section&z.nearest_state_in_previous_section].itertuples():ax.text(r.to_state,r.from_state,'○',ha='center',va='center',fontsize=17,color='red')
        fig.colorbar(im,ax=ax,shrink=.7,label='Empirical W² cost')
    fig.suptitle('Separate sections: exploratory state similarity, never reconstructed cell lineage',fontsize=15)
    write(fig,'11_cross_sectional_state_correspondence',pdf,'Fixed6 mean-state archetypes/age; no stability claim or probability flow. Circles: reciprocal nearest medoids. Not tracked cells, lineage conversion or true temporal rates.')

def main():
    review_table()
    with PdfPages(OUT/'cardiac_niche_coordination_report.pdf') as pdf:
        textpage(pdf,'Vascular growth, cardiac niches and CM maturation',[
          ('What the two broad contexts mean','The composition-only K2 solution separates stromal/mural-rich neighborhoods from CM/EC-rich myocardial neighborhoods. These are recurring neighbor mixtures, not two new cell types, exhaustive anatomical compartments, or a reproduction of the atlas paper’s expression-based12 niches.'),
          ('What cell concentration means','Smooth-muscle, epicardial and immune self-contacts, and early fibroblast–immune contacts, exceed local abundance-matched expectations with both native and refined boundaries. This supports spatial organization; it does not identify a secreted signal or prove an effect on CM maturation.'),
          ('How sprouting and maturation can coordinate','The2021 paper establishes perinatal expansion from a pre-existing coronary plexus and DLL4/NOTCH feedback restraint. The2026 atlas finds positive P7 vascular-proximity/maturation associations and proposes PTPRM contact signaling. Vessel growth and subsequent integration can support CM maturation without a universal negative local sprouting–maturity relationship.'),
          ('What this extension adds','Measured coronary/endocardial/Notch proxies; whole-RNA library support and coherence; native-polygon and refined-boundary comparisons; graph rings with censor flags; source/receiver controls; source-selection audits; and an actually trained, held-out Wormhole benchmark. Unmeasured sprout transcripts are never invented.'),
          ('Interpretation boundary','One spatial section per age confounds age with specimen/section. RNA expression is not oxygenation, perfusion, verified sprouting or signaling activation. Associations and cross-sectional state correspondences do not establish causal influence or track developmental ancestry.')
        ],'00_summary')
        contexts(pdf);RNA(pdf);maps(pdf);rings(pdf);effects(pdf);source_audit(pdf);wormhole(pdf)
        textpage(pdf,'What is supported, and what needs direct validation',[
          ('Supported spatial structure','Broad neighbor mixtures and selected contact enrichments survive changes in segmentation geometry. The measured growth-compatible EC proxy is more frequent among eligible ECs at P7 than P0; it remains a panel-defined state, not an angiogenesis rate.'),
          ('EC–CM maturation interpretation','The small positive P7 EC-sender PTPRM association in50µm models has intervals including zero for paired native/refined boundaries and approaches zero using the public capillary/High.Myo sender definition. It is a candidate association, not an independently established EC-specific mechanism. Growth/cycling neighborhoods do not consistently show lower overall CM maturation at bothP0 andP7.'),
          ('Method selection','Wormhole was trained on3200 source-cell clouds, assessed on800 regional held-out clouds, and benchmarked against exact/smoothed empirical transport and simple means. It fails the all-age metric gate and adds little held-out maturation prediction. Mean-state features/exact medoid transport remain descriptive comparisons; no true niche dynamics or new stable6-state taxonomy is claimed.'),
          ('Next experiments','Use independent hearts/sections and anatomically registered inner/outer myocardium; quantify vascular branches and EC EdU alongside the CM maturation panel; validate candidate EC–CM boundaries by imaging; compare CM-directed and EC-directed perturbations to resolve sender necessity. The2021 paper’s Hypoxyprobe and lineage-tracing evidence cannot be replaced by sparse RNA scores.'),
          ('Sources','Lu etal., Nature Cell Biology23,967–977(2021), doi:10.1038/s41556-021-00747-1, full uploaded13-page PDF read. Wang etal., Spatially guided in vivo single-cell functional genomics of postnatal heart, doi:10.1038/s44161-026-00861-z, full PMC13561863 text/XML read. See paper_review_and_pipeline_changes.tsv for figure/method-specific implications.')
        ],'12_interpretation_and_validation')

if __name__=='__main__':main()
