"""Biology-first report and integrated plates for the measured 21-gene analysis."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cm21-paperstyle-mpl')
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM21_paperstyle_extension';STAGES=['P0','P7','P14','P21']
def read(p):return pd.read_csv(OUT/p,sep='\t')
def main():
 (OUT/'figures').mkdir(exist_ok=True)
 m=read('spatial_models.tsv');s=read('score_stability.tsv');c=read('fibroblast_context_interactions.tsv');raw=read('raw_vascular_comparisons.tsv')
 primary=m[m.vascular_primary];src=m[m.nonCM_BH_q.notna()];inter=c[c.context_BH_q.notna()]
 for name,frame in [('primary_vascular_results.tsv',primary),('primary_nonCM_results.tsv',src),('primary_fibroblast_context_results.tsv',inter)]:
  frame.to_csv(OUT/name,sep='\t',index=False,float_format='%.8g',na_rep='NA')
 d={st:read(st+'_CM21_extension_values.tsv.gz') for st in STAGES}
 n=sum(len(x) for x in d.values())
 scorechecks=[]
 for st,v in d.items():scorechecks.append(dict(stage=st,CMs=len(v),core_native_Spearman=spearmanr(v.native_module21,v.core_module21).statistic))
 pd.DataFrame(scorechecks).to_csv(OUT/'core_native_score_agreement.tsv',sep='\t',index=False)
 # Preserve source-definition and score disagreement, not just selected significant findings.
 robust=[]
 for (st,source),v in m[m.exposure.eq('source_nearness')].groupby(['stage','source']):
  a=v[v.view.eq('native')&v.score.eq('module_ensemble')].iloc[0]
  robust.append(dict(stage=st,source=source,primary_beta=a.beta_SD,primary_q=a.nonCM_BH_q,
   score_view_comparisons=len(v),minimum_beta=v.beta_SD.min(),maximum_beta=v.beta_SD.max(),
   same_direction_all_scores_and_views=bool((np.sign(v.beta_SD)==np.sign(a.beta_SD)).all()),
   all_CIs_exclude_zero_in_primary_direction=bool((v.CI_low.gt(0).all() if a.beta_SD>0 else v.CI_high.lt(0).all()))))
 rb=pd.DataFrame(robust);rb.to_csv(OUT/'nonCM_score_and_segmentation_robustness.tsv',sep='\t',index=False)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 def finish(fig,stem):
  fig.savefig(OUT/'figures'/(stem+'.png'),dpi=300);fig.savefig(OUT/'figures'/(stem+'.pdf'),dpi=200);plt.close(fig)
 # Figure 1: show the actual paper-style outcome rather than reusing previous score maps.
 fig,axs=plt.subplots(2,3,figsize=(13.8,9.5),layout='constrained')
 pooled=np.concatenate([v.native_module21 for v in d.values()]);lo,hi=np.quantile(pooled,[.02,.98]);last=None
 for ax,st in zip([axs[0,0],axs[0,1],axs[1,0],axs[1,1]],STAGES):
  v=d[st].sample(frac=1,random_state=18)
  last=ax.scatter(v.x_centroid,v.y_centroid,c=v.native_module21,s=.65,cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
  ax.set(title=f'{st} · {len(v):,} ventricular CMs',xlabel='x (µm)',ylabel='y (µm)');ax.set_aspect('equal');ax.invert_yaxis()
 fig.colorbar(last,ax=[axs[0,0],axs[0,1],axs[1,0],axs[1,1]],shrink=.65,label='Measured 21-gene module score')
 ax=axs[0,2]
 groups=[('Control draws',s.score.str.startswith('draw_')),('One gene omitted',s.score.str.startswith('omit_')),('Other control pools',s.score.isin(['literal100_3bins','fine24_5controls','exclude21_8bins_20controls'])),('Earlier gene-Z score',s.score.eq('earlier_geneZ'))]
 for i,(label,mask) in enumerate(groups):
  values=s.loc[mask,'Spearman_to_primary'];ax.plot([values.min(),values.max()],[i,i],lw=4,color='#237d9c');ax.scatter(values.median(),i,color='#162f48',zorder=3)
 ax.set(yticks=range(4),yticklabels=[v[0] for v in groups],xlabel='Rank agreement with primary score',xlim=(.70,1.01),title='a  Is the spatial score stable?');ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
 ax=axs[1,2]
 for i,st in enumerate(STAGES):
  a=m[m.stage.eq(st)&m.view.eq('native')&m.adjustment.eq('Within_local')&m.exposure.eq('capillary_nearness')&m.score.str.startswith('omit_')]
  p=primary[primary.stage.eq(st)&primary.exposure.eq('capillary_nearness')].iloc[0]
  ax.plot([a.beta_SD.min(),a.beta_SD.max()],[i,i],lw=4,color='#e29038');ax.scatter(p.beta_SD,i,color='#162f48',zorder=3)
 ax.axvline(0,color='#666666',lw=1);ax.set(yticks=range(4),yticklabels=STAGES,xlabel='Adjusted score–nearness coefficient (SD)',title='b  Does one gene drive the result?');ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
 fig.suptitle('The measured 21-gene maturation landscape is stable to individual-gene omission',fontsize=15)
 finish(fig,'Figure_1_CM21_score_landscape_and_stability')
 # Figure 2: separate geometry questions, and distinguish raw from conditional association.
 fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
 metrics=[('capillary_nearness','a  Nearness to capillary ECs'),('capillary_count50','b  Capillary profiles within 50 µm'),('capillary_apposition2','c  CM boundary apposition within 2 µm')]
 for ax,(metric,title) in zip(axs.ravel()[:3],metrics):
  for view,off,color in [('native',-.10,'#236c9d'),('core',.10,'#d47b28')]:
   v=m[m.view.eq(view)&m.score.eq('module_ensemble')&m.adjustment.eq('Within_local')&m.exposure.eq(metric)].set_index('stage').loc[STAGES]
   yy=np.arange(4)+off;ax.errorbar(v.beta_SD,yy,xerr=np.vstack([v.beta_SD-v.CI_low,v.CI_high-v.beta_SD]),fmt='o',color=color,ms=5,capsize=2,label='Whole-cell RNA' if view=='native' else 'Boundary-supported RNA')
   if view=='native':
    for i,st in enumerate(STAGES):
     if v.loc[st,'vascular_BH_q']<.05:ax.text(v.loc[st,'CI_high']+.003,i+off,'*',va='center',fontsize=14)
  ax.axvline(0,color='#777777',lw=1);ax.set(title=title,yticks=range(4),yticklabels=STAGES,xlabel='Adjusted association with 21-gene score (SD)');ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
 axs[0,0].legend(frameon=False,fontsize=9,loc='lower right')
 mat=raw[raw.view.eq('native')].pivot(index='exposure',columns='stage',values='raw_Spearman').reindex(index=[v[0] for v in metrics],columns=STAGES)
 ax=axs[1,1];im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-.12,vmax=.12)
 ax.set(xticks=range(4),xticklabels=STAGES,yticks=range(3),yticklabels=['Nearness','Local count','Boundary apposition'],title='d  Raw whole-cell score relationships')
 for i in range(3):
  for j in range(4):ax.text(j,i,f'{mat.iloc[i,j]:.3f}',ha='center',va='center')
 fig.colorbar(im,ax=ax,shrink=.65,label='Individual-CM Spearman correlation')
 fig.suptitle('P7 and P14 show modest vascular relationships across different spatial measures',fontsize=15)
 finish(fig,'Figure_2_CM21_capillary_geometry')
 # Figure 3: full subtype screen plus explicit tests of fibroblast context.
 sources=['Endothelial','EC:Capillary-like','EC:Endocardial-like','EC:Arterial-like','EC:Venous-like','EC:Lymphatic-like','Fibroblast','FB:Postn','FB:Gfpt2','FB:Col8a1','FB:Lsamp','FB:Kirrel3','Pericyte','Smooth muscle','Immune','Epicardial']
 families=['FB:Postn','FB:Gfpt2','FB:Col8a1','FB:Lsamp','FB:Kirrel3']
 fig,axs=plt.subplots(1,2,figsize=(12,9),layout='constrained');cmap=plt.colormaps['RdBu_r'].copy();cmap.set_bad('#eeeeee')
 for ax,data,index,key,qkey,title,lim in [(axs[0],src,sources,'source','nonCM_BH_q','a  Which neighboring subtypes\naccompany the 21-gene score?',.10),(axs[1],inter,families,'FB_family','context_BH_q','b  Does fibroblast context alter\nthe capillary relationship?',.03)]:
  mat=data.pivot(index=key,columns='stage',values='beta_SD').reindex(index=index,columns=STAGES)
  im=ax.imshow(np.ma.masked_invalid(mat),cmap=cmap,vmin=-lim,vmax=lim,aspect='auto')
  ax.set(xticks=range(4),xticklabels=STAGES,yticks=range(len(index)),yticklabels=[x.replace('EC:','EC · ').replace('FB:','FB · ').replace('Gfpt2','Gfpt2/Pi16') for x in index],title=title,box_aspect=len(index)/8)
  for i in range(len(index)):
   for j in range(4):
    if not np.isfinite(mat.iloc[i,j]):ax.text(j,i,'×',color='#888888',ha='center',va='center')
  for row in data.itertuples():
   if getattr(row,qkey)<.05:ax.text(STAGES.index(row.stage),index.index(getattr(row,key)),'*',ha='center',va='center',fontsize=15)
  fig.colorbar(im,ax=ax,shrink=.6 if key=='source' else .3,label='Adjusted coefficient (SD)')
 fig.suptitle('Subtype-specific associations do not yet establish a fibroblast-dependent vascular effect',fontsize=14)
 finish(fig,'Figure_3_CM21_nonCM_subtypes_and_FB_context')
 loo=s[s.score.str.startswith('omit_')].Spearman_to_primary.min();draw=s[s.score.str.startswith('draw_')].Spearman_to_primary.min()
 rz=s[s.score.eq('earlier_geneZ')].Spearman_to_primary
 def p(st,metric):return primary[primary.stage.eq(st)&primary.exposure.eq(metric)].iloc[0]
 def rawrow(st):return raw[raw.stage.eq(st)&raw.view.eq('native')&raw.exposure.eq('capillary_nearness')].iloc[0]
 rows=[]
 for st in STAGES:
  a=p(st,'capillary_nearness');b=p(st,'capillary_count50');e=p(st,'capillary_apposition2')
  rows.append(f'| {st} | {a.beta_SD:+.3f}; q={a.vascular_BH_q:.3g} | {b.beta_SD:+.3f}; q={b.vascular_BH_q:.3g} | {e.beta_SD:+.3f}; q={e.vascular_BH_q:.3g} |')
 rawrows=[]
 for st in STAGES:
  v=rawrow(st);rawrows.append(f'| {st} | {v.low_quartile_nearest_capillary_um:.1f} | {v.high_quartile_nearest_capillary_um:.1f} |')
 settings=json.loads((OUT/'analysis_settings.json').read_text())
 report=f'''# Cardiomyocyte maturation in its vascular and stromal surroundings

## What this analysis adds

The central question is whether ventricular cardiomyocytes with higher expression of the published 21-gene Xenium maturation panel occupy different vascular and fibroblast neighborhoods. All {n:,} previously eligible ventricular CM profiles are retained across P0, P7, P14 and P21. The positive maturation module contains only those 21 directly measured genes; other measured targets supply normalization and background controls. The existing snRNA-seq-supported annotations identify neighboring cell populations; they do not supply or impute the maturation genes.

The strongest result is at P7. Higher-scoring CMs are closer to capillary endothelial cells, have more capillary profiles nearby and have a greater fraction of their segmented boundary close to capillary profiles. At P14, nearness and boundary apposition show the same direction. These are small spatial associations, but agreement across several measurements makes the vascular relationship more credible than a single correlation between averaged neighborhood rings.

## The 21-gene pattern survives changes in individual genes and controls

The paper calculates a module score by subtracting expression-matched control genes from the mean normalized expression of the maturation genes. We reconstructed this approach with explicitly recorded settings suited to the 403-target panel and averaged 50 control draws to reduce dependence on random control selection. The authors' original controls, exact score object, anatomical masks and graph are unavailable, so this is a paper-style reconstruction rather than an exact replication.

Removing any one maturation gene leaves the cell ranking strongly similar to the full panel: the lowest rank correlation across stages and RNA views is {loo:.3f}. Each omission removes that gene from both the maturation module and the control pool. Changing the random controls also has little effect: the lowest correlation between an individual draw and the ensemble is {draw:.3f}. Alternative bin settings, including a setting permitting 100 controls per gene and a pool excluding all 21 maturation genes, preserve high rank agreement. Thus the result is not attributable to one naming gene or one favorable control draw (Fig. 1).

Changing transcript assignment has a larger effect on individual-cell rankings. Whole-cell and boundary-supported score rankings correlate at {min(v['core_native_Spearman'] for v in scorechecks):.2f}–{max(v['core_native_Spearman'] for v in scorechecks):.2f} across stages. A CM's higher- or lower-score classification therefore remains sensitive to which transcripts are assigned to it. Agreement in the direction of the vascular association across these RNA views does not mean that every cell retains the same score or rank.

There is an important difference from our earlier report. The previous score standardized each maturation gene separately before averaging. Its rank correlation with the paper-style score is only {rz.min():.2f}–{rz.max():.2f}. The two approaches therefore identify overlapping but different higher-scoring CM groups. With the paper-style score, the raw P7 comparison and the adjusted comparison now agree in direction. This explains why the earlier raw P7 result differed; it does not recover the authors' exact score or establish that one numerical index is functional maturity.

## The capillary relationship is clearest at P7 and supported at P14

The direct high-versus-low comparison uses score quartiles within each stage. At P7, the higher-scoring group is closer to capillary profiles on average. P14 shows the same direction. These are descriptive means across the section, before matching or local adjustment.

| Stage | Lower-score quartile: mean capillary distance, µm | Higher-score quartile: mean capillary distance, µm |
|---|---:|---:|
{chr(10).join(rawrows)}

The locally adjusted analysis compares CMs within 200-µm regions and the same available CM subtype, while accounting for RNA capture, cell area and missing nuclear area. A positive coefficient means a higher score with greater nearness, more capillary profiles, or more boundary apposition. The three measurements are related descriptions of tissue geometry; they do not represent independent mechanisms.

| Stage | Capillary nearness | Capillary count within 50 µm | Boundary apposition within 2 µm |
|---|---:|---:|---:|
{chr(10).join(rows)}

Coefficients are standardized associations, and q values correct the 12 whole-cell vascular comparisons together. At P7, all three pass that correction. At P14, nearness and apposition pass; the count comparison does not. The P0 and P21 comparisons do not pass the same correction. Different statistical support across ages does not demonstrate when a biological relationship begins or ends (Fig. 2).

The capillary-nearness coefficient remains positive after each individual maturation gene is omitted, with every tested control setting, and in boundary-supported RNA. At P7, alternative local grid definitions also retain a positive relationship with confidence intervals above zero. At P14, the direction remains positive, although uncertainty reaches zero in one grid sensitivity. Distance along the tissue-screened cell graph and distance between segmented boundaries give the same direction at P7 and P14. The graph avoids some connections across empty spaces but is not a reconstruction of perfused vessels.

Capillary annotation is a remaining source of uncertainty. At P7, the positive association persists when using the broader public capillary labels after EC QC. Requiring strict cell-and-nucleus subtype agreement reduces the source set from 1,464 to 268 profiles and weakens the estimate to +0.008 SD, with an interval spanning zero. At P14, the broader public-label estimate is weaker, and the strict 202-profile estimate is inconclusive with a slightly negative point estimate. The evidence is therefore strongest for the reviewed and broader P7 capillary maps; it is not invariant to every source definition. Sparse, highly selected source sets are not independent ground truth for the whole network.

## Nearness is more informative than simply counting nearby cells

More capillary profiles can reflect a generally denser tissue neighborhood. After including total nearby resolved cell count, the P7 capillary-count interval includes zero, whereas capillary nearness retains a positive relationship. The P14 nearness relationship also remains positive in this sensitivity. The present evidence therefore favors a relationship with local vascular arrangement over a claim that any increase in nearby EC number predicts a higher maturation score.

Boundary apposition provides a complementary observation: a higher score accompanies a greater fraction of the CM outline lying within 2 µm of capillary profiles at P7 and P14. These are native segmented boundaries sampled at 64 evenly spaced positions. They are not confirmed membrane contacts, and the boundaries themselves can be inaccurate. The boundary-supported RNA sensitivity challenges transcript ownership while keeping the same native geometry; it does not independently validate or replace the segmentation.

## The neighboring subtype matters, but fibroblast effects remain less settled

Capillary ECs and endocardial-like ECs do not show the same relationship. In the paper-style whole-cell analysis, endocardial-like EC nearness accompanies lower scores at P7, P14 and P21, while capillary nearness accompanies higher scores at P7 and P14. This distinction argues against treating all endothelial cells as equivalent. Endocardial-like ECs also mark different anatomical environments, so regional position remains a plausible explanation. These observations do not mean that endocardial cells inhibit maturation.

Fibroblast relationships are more dependent on measurement and score definition. Lsamp-associated fibroblast nearness accompanies lower scores at P14 and P21 in whole-cell and boundary-supported RNA, consistent with the direction in the preceding 21-gene analysis. This negative direction also persists with the cell-supported, boundary-supported and crosschecked fibroblast assignments. Col8a1-associated proximity is negatively associated in the new whole-cell score screen, but the boundary-supported RNA intervals include zero. It should therefore remain a boundary- or measurement-sensitive finding, rather than a conclusion that Col8a1 fibroblasts delay maturation. Postn and Gfpt2/Pi16 nearness do not pass the new whole-cell source-screen correction where assessable. In particular, the earlier modest P7 Gfpt2/Pi16 result is not consistently recovered across scoring choices.

The full source screen also includes pericytes, smooth muscle, immune and epicardial cells, and supported EC subtypes. Results and score/segmentation disagreements are provided together (Fig. 3). Source categories overlap with their parent lineages and are not independent biological explanations. Unknown or insufficiently supported source families remain not assessable; they are not counted as absent.

## Does fibroblast context change the capillary relationship?

We tested whether the association between capillary nearness and the 21-gene score differs as the number of nearby fibroblasts of a particular family increases. These models include capillary nearness, fibroblast context, their interaction and total local resolved cell count, alongside the technical and local-region adjustments.

None of the 12 assessable whole-cell interactions passes correction for multiple comparisons. P7 Postn context gives a positive nominal interaction, but its corrected q is {inter[(inter.stage=='P7')&(inter.FB_family=='FB:Postn')].context_BH_q.iloc[0]:.3f}. We therefore cannot conclude that a Postn-rich, Gfpt2/Pi16-rich or Col8a1-rich environment strengthens or weakens the capillary relationship. This is a direct test of the proposed combined vascular–matrix niche, and its present answer is inconclusive rather than positive.

## Biological interpretation and next experiment

The most defensible interpretation is that CM transcriptional maturation and local coronary organization are associated, with the clearest evidence in the P7 section and more source-definition-sensitive support at P14. Nearness and operational boundary apposition are more consistent than capillary abundance alone. The score is heterogeneous within a stage, and this heterogeneity relates to the identity and arrangement of neighboring cells. Capillary-label uncertainty prevents calling the vascular relationship fully resolved.

These findings make local vascular integration a reasonable hypothesis to test. They do not distinguish blood flow, oxygen delivery, contact-dependent signals, matrix organization or shared regional development. A useful next experiment would compare matched ventricular regions in independent P7 and P14 hearts, measuring the same 21-gene score alongside perfused capillaries, validated CM–EC boundaries and a functional CM measurement. Earlier P3/P5 and intermediate P10 samples would establish how these relationships develop; they cannot be supplied by interpolation from the existing four stages.

## Methods and scope

Counts are normalized as log1p(10,000 times target count divided by all 403 measured-target counts). Only the 21 published maturation genes enter the positive module. The primary paper-style score uses eight rank-balanced expression bins, 20 controls per module gene sampled without replacement, the unique control union within each draw, and the mean of 50 independently seeded draw scores. Controls may include module genes, as in the default Seurat pool. Bins use frozen pooled all-profile expression means and are fixed across stages within each RNA view. Rank binning is explicit and is not identical to Seurat's implementation. Alternative settings use three bins/100 controls, 24 bins/five controls, and eight bins/20 controls excluding the 21 panel genes. These settings challenge the reconstruction; none is asserted to recover the authors' original pool or bins. The 403-target panel alone cannot support 100 distinct controls in each of 24 balanced bins, which is a feasibility issue for that parameter combination, not evidence that the authors made an error.

The primary view is deposited whole-cell Xenium RNA, matching the paper's measurement context. Confidence-supported core RNA is a transcript-assignment sensitivity and is not nuclear-only truth. The prior all-QC ventricular cohort is unchanged: no anatomical exclusion, public-label agreement or visible-nucleus requirement is added. Non-CM source labels and their eligibility/assessability rules are inherited from the reviewed reference. Nuclear-supported EC labels and nuclear-optional multi-gene fibroblast proposals remain separate evidence layers. No new segmentation or new subtype clustering was performed here.

Local regressions use view-specific non-panel molecule counts and cell area as cubic splines, missing nuclear area, available CM subtype and 200-µm grid-by-CM-state fixed effects. Unknown CM subtype and unknown candidate anatomy are retained. Confidence intervals use 400-µm clustered uncertainty with an absorbed-group degrees-of-freedom correction. Primary vascular q values use 9,999 null-imposed clustered-score sign draws, followed by BH across 12 tests; this is not vessel relocation. Non-CM nearness and fibroblast-context tests have separate BH families of 48 and 12. All uncertainty concerns these observed sections, not independent animals. Supporting tests include alternative score settings, all 21 gene omissions, core RNA, existing candidate anatomy, local grid changes, 800-µm uncertainty blocks, 30-µm counts and alternative graph/boundary distances. Supporting comparisons are not additional discoveries.

Source-label sensitivities hold the CM score and cohort fixed while substituting EC-QC-supported public capillary labels, strict cell-and-nucleus EC calls, or cell-supported/boundary-supported/crosschecked fibroblast family calls. Their estimates and source coverage are provided without selecting a favorable source definition. The native boundaries and tissue-screened graphs are inherited, not newly expert-validated anatomical masks. Independent validation re-reads count matrices, reconstructs all 50 control draws, verifies complete cohort retention, recomputes local capillary regressions and rereads original boundary vertices for sampled CMs. One observed section per stage prevents animal-level replication, and no result establishes a causal effect or functional CM maturation. Source: Wang et al., Nature Cardiovascular Research, DOI 10.1038/s44161-026-00861-z; https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/.
'''
 (OUT/'Report_text.md').write_text(report)
 legends=[dict(stem='Figure_1_CM21_score_landscape_and_stability',title='The 21-gene score is heterogeneous within stages and stable to gene omission',legend='All 131,605 previously QC-eligible ventricular CM profiles are retained. Four maps use the same pooled color scale for the paper-style whole-cell module score; colors are clipped at the pooled second and 98th percentiles for display only. Right panels show rank agreement across scoring challenges and the range of capillary-nearness coefficients after omitting each of the 21 genes; points mark primary coefficients. Agreement is summarized across all stages and both RNA views. The earlier gene-Z score differs in weighting and yields less similar cell rankings. No new CM outcome genes are introduced.'),
 dict(stem='Figure_2_CM21_capillary_geometry',title='Capillary nearness and boundary apposition support a local vascular relationship',legend='Panels a–c show locally adjusted standardized associations with the same 21-gene outcome, with 95% within-section clustered intervals. Blue is whole-cell RNA; orange is boundary-supported RNA using the same native geometry. Stars apply only to primary whole-cell tests after correction across 12 vascular comparisons using clustered-score null tests. Panel d displays raw individual-CM Spearman correlations, not correlations between ring means; its units differ from the adjusted coefficients. Counts measure segmented profiles within 50 µm, not perfused vessel length. Boundary apposition is the fraction of 64 native CM boundary samples within 2 µm of eligible capillary polygons. Associations are small; no causal or functional inference is made.'),
 dict(stem='Figure_3_CM21_nonCM_subtypes_and_FB_context',title='Neighbor identities and fibroblast-context tests narrow the biological conclusions',legend='Panel a shows the whole-cell module-score nearness screen, with BH correction across 48 supported source comparisons. Panel b shows capillary-nearness-by-fibroblast-count interactions, with correction across 12 supported comparisons; no interaction passes correction. Positive interaction coefficients indicate a stronger positive capillary relationship with increasing fibroblast-family context, conditional on the fitted main effects and covariates. Each panel has its own color scale. Gray crosses denote insufficient or unavailable source evidence rather than biological absence. Significant primary source associations can remain sensitive to score definition and transcript assignment; all such comparisons are provided in the robustness table. Parent lineages overlap with their subtypes.')]
 legends[2]['legend']+=' Colors saturate at ±0.10 SD in panel a and ±0.03 SD in panel b; full coefficients are retained in the tables.'
 (OUT/'figure_legends.json').write_text(json.dumps(legends,indent=2)+'\n')
 (OUT/'README.md').write_text('''# Measured 21-gene CM maturation: paper-style scoring and spatial challenges

[Word report](CM21_paperstyle_vascular_stromal_report.docx) · [PDF report](CM21_paperstyle_vascular_stromal_report.pdf) · [Full narrative](Report_text.md)

The unchanged 131,605-CM cohort is scored using only the published 21 measured Xenium targets. The paper-style module ensemble is the main score; this is not an exact author replication. P7 has concordant raw and conditional capillary associations; P14 has support for nearness and operational boundary apposition. All 21 omissions and control choices are retained. Fibroblast-context interaction tests remain inconclusive. Whole-cell/core and score-dependent disagreements are preserved.

Main tables: `spatial_models.tsv`, `fibroblast_context_interactions.tsv`, `score_stability.tsv`, `nonCM_score_and_segmentation_robustness.tsv`, `source_label_sensitivity.tsv`. All-cell coordinates, scores, ranks and geometry: `P*_CM21_extension_values.tsv.gz`. Every module/draw/omission weight is in `score_gene_weights.tsv.gz`; exact 21-gene list and settings accompany it. All source labels, unknowns and reference eligibility rules are inherited. Existing boundaries are reused; new segmentation and original author masks are not claimed. One section per stage is not animal replication. P7 source-label sensitivity attenuates the strict cell-and-nucleus comparison, and P14 support is source-definition-sensitive.

Run `analysis/cm21_paperstyle_extension.py`, `analysis/cm21_source_label_sensitivity.py`, `analysis/validate_cm21_paperstyle.py` and `analysis/report_cm21_paperstyle.py` with `/workspace/heart-atlas/venv/bin/python`; run `analysis/write_cm21_paperstyle_report.py` with base Python. Validation independently rereads count matrices, rebuilds control draws, verifies omissions, reconstructs partial-regression estimates and rereads sampled original boundaries.
''')
 print('Report and three integrated plates written',flush=True)
if __name__=='__main__':main()
