"""Maps and a biological report centered exclusively on the Xenium21-gene score."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cm21-primary-mpl')
from pathlib import Path
import json,re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM21_primary_spatial';STAGES=['P0','P7','P14','P21']
def read(p):return pd.read_csv(OUT/p,sep='\t')
def prose(s):return re.sub(r'\b(The|the|same|All|all|Only|only|Higher|Lower|higher|lower|measured|within|local|least|of|different|pooled|one|use|uses|Xenium)(?=[0-9])',r'\1 ',s)
def main():
 (OUT/'figures').mkdir(exist_ok=True)
 cohort=read('CM21_stage_cohorts.tsv');mod=read('continuous_CM21_associations.tsv');primary=mod[mod.primary_analysis]
 raw=read('all_CM21_high_low_neighborhood_descriptions.tsv');pairs=read('high_low_CM21_matched_pairs.tsv.gz');matched=read('high_low_CM21_matched_comparisons.tsv')
 data={s:read(s+'_CM21_scores_and_neighborhoods.tsv.gz') for s in STAGES}
 summaries=[];rawcorr=[]
 for s,d in data.items():
  summaries.append(dict(stage=s,CMs=len(d),score21_median=d.core_score21.median(),score21_p10=d.core_score21.quantile(.1),score21_p90=d.core_score21.quantile(.9),
   median_detected_21genes=d.core_detected_panel_genes.median(),score_capture_Spearman=spearmanr(d.core_score21,d.CM_outcome_excluded_depth).statistic))
  for col in d:
   if col.startswith('distance__'):rawcorr.append(dict(stage=s,source=col[10:],raw_score_nearness_Spearman=spearmanr(d.core_score21,-np.log1p(d[col])).statistic))
 pd.DataFrame(summaries).to_csv(OUT/'CM21_stage_score_distributions.tsv',sep='\t',index=False)
 pd.DataFrame(rawcorr).to_csv(OUT/'raw_CM21_neighborhood_correlations.tsv',sep='\t',index=False)
 def model(stage,source,view='core',adjustment='Within_local_region'):
  a=mod[mod.stage.eq(stage)&mod.source.eq(source)&mod.view.eq(view)&mod.adjustment.eq(adjustment)];assert len(a)==1;return a.iloc[0]
 def descriptive(stage,source):
  a=raw[raw.stage.eq(stage)&raw.source.eq(source)&raw.view.eq('core')].set_index('score_group')
  return a.loc['Lower quartile'],a.loc['Higher quartile']
 caprows=[]
 for s in STAGES:
  low,high=descriptive(s,'EC:Capillary-like');m=model(s,'EC:Capillary-like')
  caprows.append(f'| {s} | {low.mean_nearest_um:.1f} | {high.mean_nearest_um:.1f} | {m.beta_SD:+.3f} | {m.primary_BH_q:.3g} |')
 cohortrows='\n'.join(f'| {r.stage} | {r.QC_eligible_ventricular_CMs:,} | {r.lower_quartile:,} | {r.higher_quartile:,} |' for r in cohort.itertuples())
 g7=model('P7','FB:Gfpt2');ls14=model('P14','FB:Lsamp');ls21=model('P21','FB:Lsamp');cap7=model('P7','EC:Capillary-like');cap14=model('P14','EC:Capillary-like')
 nmatched=int(pairs[pairs.view.eq('core')].shape[0]);nall=int(cohort.QC_eligible_ventricular_CMs.sum())
 text=f'''# The21-gene CM maturation landscape and neighboring cardiac cells

## The question answered here

Where are ventricular cardiomyocytes with higher or lower21-gene maturation scores at each postnatal stage, and which non-CM subtypes surround them? Every maturation map, comparison and model in this report uses the same21 genes measured directly by Xenium. The existing RNA-supported annotations identify the neighboring cell populations.

The maps retain all {nall:,} reviewed ventricular CM profiles passing the existing QC thresholds. Anatomical classification and agreement with public labels are not requirements for inclusion. CMs in uncertain anatomical regions and valid CM sections without a visible nucleus remain in the maps. Higher and lower groups are defined separately within each age, so a higher-scoring P0 CM is not assumed to have the same absolute maturity as a higher-scoring P21 CM.

| Stage | All QC-eligible ventricular CMs | Lower21-score quartile | Higher21-score quartile |
|---|---:|---:|---:|
{cohortrows}

## Where are higher- and lower-scoring CMs?

Figure1 shows the directly measured21-gene score on each section, using the same color scale across ages. Figure2 highlights the higher and lower quartiles within each stage. Together, these maps show the age-related shift in the index and the remaining variation among CMs in the same section.

The complete cell-level tables contain each CM's coordinates,21-gene score, within-stage rank, high/low group and distances to supported non-CM sources. The interactive viewer allows selection of a stage and neighboring cell subtype, followed by inspection of individual CM scores and local source positions. All21 targets are present in the Xenium panel; individual cells need not detect every target. The median number detected per CM is10 at P0,11 at P7 and13 at P14/P21.

## Coronary capillary proximity: the raw and adjusted answers

At P14, higher-scoring CMs are closer to capillary ECs in the direct high-versus-low comparison, and the continuous-score analysis gives the same direction. At P7, the raw high-versus-low comparison does not show that pattern: the highest-scoring quartile is slightly farther from the QC-supported capillary sources. After accounting for RNA capture, cell size, CM identity and local position, a small positive capillary association appears. The statement that higher-scoring CMs are closer to capillaries at P7 is therefore conditional on those adjustments; it is not a simple unadjusted rule.

| Stage | Lower-score CMs: mean nearest capillary, µm | Higher-score CMs: mean nearest capillary, µm | Adjusted association, SD | All-source BH q |
|---|---:|---:|---:|---:|
{chr(10).join(caprows)}

Positive adjusted coefficients mean that a higher21-gene score accompanies greater capillary nearness. The adjusted coefficients are {cap7.beta_SD:+.3f} at P7 and {cap14.beta_SD:+.3f} at P14. These directions persist with deposited RNA and when anatomical labels are added while retaining unclassified regions. The P0 and P21 capillary comparisons do not pass the same all-source correction. This does not establish a biological onset or disappearance of a vascular relationship.

RNA capture is strongly related to the measured index and also varies spatially. That is why the raw and adjusted views are displayed together. Adjustment can reduce technical confounding, but it does not make one observational model a definitive account of maturation or establish that neighboring vessels cause the score differences.

## Which fibroblast subtypes accompany different21-gene scores?

The full-cohort continuous analysis identifies a small positive relationship with Gfpt2/Pi16-associated fibroblast nearness at P7: coefficient {g7.beta_SD:+.3f}, q={g7.primary_BH_q:.3f}. In contrast, Lsamp-associated fibroblast nearness relates to lower scores at P14 ({ls14.beta_SD:+.3f}, q={ls14.primary_BH_q:.3f}) and P21 ({ls21.beta_SD:+.3f}, q={ls21.primary_BH_q:.3f}). These directions persist in deposited CM RNA and the anatomy-covariate sensitivity.

These are conditional spatial associations. The raw high/low summaries do not always show the same contrast, so the results should not be translated into claims that these fibroblasts promote or suppress maturation. Nor does an Lsamp label establish an outflow-tract location.

Postn-associated fibroblast proximity does not show a corrected whole-cohort relationship at P0 or P7. Col8a1-associated fibroblast proximity does not show a corrected whole-cohort relationship at P14 or P21. Later Gfpt2/Pi16 assignments remain uncertain, and no corrected P21 association is established. A relationship in a selected compact-wall patch should not be generalized to the entire section. Unavailable reference families and inadequately supported sources remain not assessable, rather than being counted as biologically absent.

## Other neighboring cardiac cells

The five supported EC classes are evaluated where sufficient cells are available. Endocardial-like EC nearness relates to lower21-gene scores at P7 in the adjusted model. Arterial-like and venous-like sources do not show corrected associations in this full-cohort screen. Lymphatic-like sources lack sufficient support for a reliable comparison. Broad EC results overlap with their capillary subtype results and are not independent evidence for a separate mechanism.

Among other broad populations, pericyte nearness relates to slightly higher scores at P14 and P21; immune-cell nearness does so at P14. Smooth-muscle nearness relates to lower scores at P21, and epicardial nearness relates to lower scores at P7–P21 in the adjusted models. These patterns can reflect tissue organization and regional CM states. They do not establish trophic, inhibitory or signaling effects, and validated fine immune, mural and epicardial subtypes are not supplied by the present reference.

Figure3 presents every supported source rather than only selected findings. The left panel shows the raw score–nearness relationship, and the right panel shows the adjusted continuous-score relationship. The full high/low neighborhood table also gives distances and the proportion of CMs with each source within50µm. These are alternative descriptions of the same21-gene outcome, not additional maturation programs.

## What does strict high/low matching add?

The primary continuous analysis uses all eligible CMs. As a separate sensitivity, higher- and lower-score CMs were matched within CM identity and local regions for RNA capture, cell area and missing nuclear area, without matching on neighboring cells. Only {nmatched:,} disjoint pairs across the four ages satisfy the strict calipers. None of the matched distance comparisons survives correction across the source screen.

This is important uncertainty. The matched analysis does not independently confirm the source-specific findings, and its small selected subset cannot replace the all-CM question. The strongest defensible interpretation is that the21-gene maps show within-section heterogeneity and identify modest, model-dependent relationships with particular vascular and stromal surroundings. Independent-heart replication is still needed.

## Methods and scope

Only the published21 Xenium maturation targets form the CM outcome. Counts are normalized as log1p(10,000 times target count divided by all403 measured-target counts). The score is the mean of21 gene Z values using the unchanged pooled reference calibration, separately for deposited and core RNA. This preserves the score used in our earlier reports; it is not a recovery of the authors' unavailable exact Seurat AddModuleScore control set. No additional CM outcome genes, inferred gene modules or snRNA-seq maturation scores enter this analysis.

The ventricular cohort is defined by the reviewed ventricular-CM flag, broad CM identity, expression-QC pass, at least20 core transcripts and at least20 pre-existing outcome-excluded counts. No published-label concordance, conservative ventricular mask or manually annotated region is required. Missing nuclear area is retained. Existing source identities and their confidence/assessability rules are unchanged. Core RNA is a confidence-supported transcript-assignment layer, not nuclear-only truth.

The full-cohort continuous model adjusts for flexible log-capture and log-area effects, missing nuclear area, available CM identity, and200-µm local regions. Published non-CM labels on otherwise reviewed ventricular CMs are not used as CM subtype identities; those cells retain unresolved CM subtype. The anatomy sensitivity adds available image-region labels and retains unclassified cases. Within-section confidence intervals use400-µm clustered uncertainty. Core full-cohort source comparisons are corrected together; deposited RNA and anatomy models are measurement/context sensitivities.

Raw high/low summaries use within-stage quartiles of the same21-gene score. Strict matching requires the same200-µm grid and CM subtype, a maximum100-µm pair separation, capture and area differences no greater than0.25 full-stage SD, identical missing-nucleus status, and unique cells. Neighbor identities and counts do not determine matches. Matched intervals use1,999 resamples of400-µm pair-anchor clusters. The all-CM continuous model is the main inferential view because the strict high/low match has limited measurement overlap. This is an exploratory reanalysis, not a prospectively registered confirmation.

All sources, raw comparisons, adjusted models and unsuccessful matched tests are provided. The number of cells does not substitute for independent hearts: one observed Xenium section is available per age. The score measures a transcriptional index, and these spatial relationships remain associations.

Source: Wang et al., Spatially guided in vivo single-cell functional genomics of postnatal heart, Nature Cardiovascular Research, DOI10.1038/s44161-026-00861-z. https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/ Gene list and exact score reconstruction checks accompany this report.
'''
 # Keep prose spacing conventional despite compact code literals.
 for a,b in [('The21','The 21'),('lower21','lower 21'),('same21','same 21'),('measured21','measured 21'),('published21','published 21'),('all21','all 21'),('All21','All 21'),('per CM is10','per CM is 10'),('P0,11','P0, 11'),('and13','and 13'),('Figure1','Figure 1'),('Figure2','Figure 2'),('Figure3','Figure 3'),('within50','within 50'),('all403','all 403'),('of21','of 21'),('mean of21','mean of 21'),('least20','least 20'),('and200','and 200'),('use400','use 400'),('same200','same 200'),('maximum100','maximum 100'),('than0.25','than 0.25'),('use1,999','use 1,999'),('of400','of 400'),('DOI10.','DOI 10.')]:text=text.replace(a,b)
 (OUT/'Report_text.md').write_text(prose(text))
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 scores=np.concatenate([d.core_score21.to_numpy() for d in data.values()]);lo,hi=np.quantile(scores,[.02,.98])
 for mode in ['score','quartile']:
  fig,axs=plt.subplots(2,2,figsize=(12,11),layout='constrained');last=None
  for ax,s in zip(axs.ravel(),STAGES):
   d=data[s];f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{s}/atlas_annotations.tsv.gz',sep='\t')
   ax.scatter(f.x_centroid,f.y_centroid,c='#e4e4e4',s=.18,rasterized=True)
   rng=np.random.default_rng(281029);order=rng.permutation(len(d));q=d.iloc[order]
   if mode=='score':last=ax.scatter(q.x_centroid,q.y_centroid,c=q.core_score21,s=1,cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
   else:
    middle=d[d.core_score_group.eq('Middle half')];ax.scatter(middle.x_centroid,middle.y_centroid,c='#cccccc',s=.5,rasterized=True)
    q=q[q.core_score_group.ne('Middle half')];colors=q.core_score_group.map({'Lower quartile':'#2878b5','Higher quartile':'#de7928'})
    ax.scatter(q.x_centroid,q.y_centroid,c=colors,s=1,rasterized=True)
   ax.set(title=f'{s} · {len(d):,} ventricular CMs',xlabel='x (µm)',ylabel='y (µm)');ax.invert_yaxis();ax.set_aspect('equal')
  if mode=='score':fig.colorbar(last,ax=axs.ravel().tolist(),shrink=.7,label='Directly measured 21-gene maturation score');title='The same21-gene index across all four sections';stem='Figure_1_CM21_full_section_maps'
  else:
   fig.legend(handles=[Line2D([],[],color='#2878b5',marker='o',ls='',label='Lower score quartile'),Line2D([],[],color='#de7928',marker='o',ls='',label='Higher score quartile'),Line2D([],[],color='#cccccc',marker='o',ls='',label='Middle half')],loc='outside lower center',ncol=3,frameon=False)
   title='Higher and lower21-gene maturation scores within each stage';stem='Figure_2_CM21_within_stage_heterogeneity'
  fig.suptitle(title.replace('same21','same 21').replace('lower21','lower 21'),fontsize=15)
  fig.savefig(OUT/'figures'/(stem+'.png'),dpi=300);fig.savefig(OUT/'figures'/(stem+'.pdf'),dpi=300);plt.close(fig)
 order=['Endothelial','EC:Capillary-like','EC:Endocardial-like','EC:Arterial-like','EC:Venous-like','EC:Lymphatic-like','Fibroblast','FB:Postn','FB:Gfpt2','FB:Col8a1','FB:Lsamp','FB:Kirrel3','Pericyte','Smooth muscle','Immune','Epicardial']
 rawmat=pd.DataFrame(rawcorr).pivot(index='source',columns='stage',values='raw_score_nearness_Spearman').reindex(index=order,columns=STAGES)
 adjmat=primary.pivot(index='source',columns='stage',values='beta_SD').reindex(index=order,columns=STAGES)
 fig,axs=plt.subplots(1,2,figsize=(12,10),layout='constrained');cmap=plt.colormaps['coolwarm'].copy();cmap.set_bad('#eeeeee')
 labels=[v.replace('EC:','EC · ').replace('FB:','FB · ').replace('Gfpt2','Gfpt2/Pi16') for v in order]
 for ax,mat,title,limit,cbar in [(axs[0],rawmat,'a  Raw score–nearness correlation',.30,'Spearman correlation'),(axs[1],adjmat,'b  After capture and local-context adjustment',.08,'Adjusted coefficient (SD)')]:
  im=ax.imshow(np.ma.masked_invalid(mat.to_numpy()),cmap=cmap,vmin=-limit,vmax=limit,aspect='auto');ax.set(yticks=range(len(order)),yticklabels=labels,xticks=range(4),xticklabels=STAGES,title=title)
  for i in range(len(order)):
   for j in range(4):
    if not np.isfinite(mat.iloc[i,j]):ax.text(j,i,'×',ha='center',va='center',color='#777777')
  fig.colorbar(im,ax=ax,shrink=.65,label=cbar)
 for row in primary.itertuples():
  if row.primary_BH_q<.05:axs[1].text(STAGES.index(row.stage),order.index(row.source),'*',ha='center',va='center',fontsize=16)
 fig.suptitle('Which neighboring cells accompany different 21-gene CM scores?',fontsize=15)
 stem='Figure_3_CM21_nonCM_relationships';fig.savefig(OUT/'figures'/(stem+'.png'),dpi=300);fig.savefig(OUT/'figures'/(stem+'.pdf'));plt.close(fig)
 legends=[dict(stem='Figure_1_CM21_full_section_maps',title='The measured21-gene score on all four sections',legend='All131,605 QC-eligible ventricular CM profiles are mapped. Color is the unchanged core21-gene score, with one shared scale; display clipping at pooled2nd/98th percentiles does not alter analysis values. Gray points provide whole-section context. No anatomical classification or public-label concordance is required. These are measured cell scores, not imputed expression.'),
  dict(stem='Figure_2_CM21_within_stage_heterogeneity',title='Higher and lower21-gene scores within each age',legend='Blue and orange mark the lower and higher score quartiles within each stage; gray marks the middle half and section context. Quartiles are stage-relative, so an orange P0 CM is not assumed equivalent to an orange P21 CM. Rank ties are averaged. Every group is defined using the same21 measured Xenium targets.'),
  dict(stem='Figure_3_CM21_nonCM_relationships',title='Raw and adjusted relationships use the same21-gene outcome',legend='The left panel shows raw score–nearness Spearman correlations across all eligible CMs; the right shows capture-, size-, CM-identity- and local200µm-region-adjusted coefficients. Positive values associate higher21-gene scores with source nearness. Stars mark all-source BH q<0.05 only in the adjusted model. Gray crosses indicate inadequate source/reference/local support. Broad lineages and their subtypes overlap, so they are not independent mechanisms. The two panels use different association statistics and should be read by direction and context, not compared numerically. All outcomes remain the measured21-gene score; no component program is substituted.')]
 for entry in legends:
  for k in ['title','legend']:entry[k]=prose(entry[k])
 (OUT/'figure_legends.json').write_text(json.dumps(legends,indent=2)+'\n')
 (OUT/'README.md').write_text('# Direct Xenium21-gene CM maturation and neighboring cardiac cells\n\n[Word report](CM21_maturation_and_neighbors.docx) · [PDF report](CM21_maturation_and_neighbors.pdf) · [Interactive21-score viewer](CM21_spatial_viewer.html) · [Narrative](Report_text.md)\n\nAll131,605 QC-eligible ventricular CMs are retained. The sole CM maturation outcome is the directly measured21-gene Xenium score. No anatomical restriction, inferred maturation gene program or missing-gene projection defines the outcome. Existing snRNA-supported cell identities label the neighboring populations.\n\nThe full-cohort continuous analysis, raw high/low neighborhoods and strict matched sensitivities are all provided. P7 capillary proximity has a small conditional positive relationship despite the raw high/low comparison not showing closer capillaries. P14 has concordant raw and adjusted directions. All48 adjusted source comparisons are corrected together. Strict matched comparisons remain inconclusive. One observed section per stage limits replication.\n\nRun cm21_primary_analysis.py, report_cm21_primary.py and validate_cm21_primary.py with the scientific Python environment; build_cm21_viewer.py and write_cm21_primary_report.py use base Python. The21-gene panel, reconstruction checks, all-cell score/coordinate tables, source distances and confidence/coverage information are included.\n')
 p=OUT/'README.md';p.write_text(prose(p.read_text()))
 print('CM21 figures and report text written',flush=True)
if __name__=='__main__':main()
