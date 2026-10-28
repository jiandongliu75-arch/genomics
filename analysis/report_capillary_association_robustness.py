"""Report all targeted strengthening checks, including the inconclusive match."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/capillary-robustness-mpl')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/capillary_association_robustness'
def main():
 r=pd.read_csv(OUT/'P7_robustness_regressions.tsv',sep='\t');pairs=pd.read_csv(OUT/'P7_local_matched_outcomes.tsv',sep='\t')
 primary=r[r.primary_analysis];matched=pairs[pairs.outcome.eq('core_maturation21')].iloc[0]
 table='\n'.join(f'| {label} | {row.beta_SD:+.3f} | {row.CI_low:+.3f} to {row.CI_high:+.3f} |' for label,row in zip(
  ['Cell-center distance','Sampled CM boundary distance','Tissue-screened graph path'],primary.itertuples()))
 source=r[(r.analysis=='geometry')&(r.exposure=='near_cell_supported_capillary')&(r.outcome=='core_maturation21')].iloc[0]
 density=r[(r.analysis=='specificity_or_density')&r.extra_controls.str.contains('log_capillary_count',na=False)&r.outcome.eq('core_maturation21')].iloc[0]
 text=f'''# Strengthening the capillary–cardiomyocyte association

## Main conclusion

The additional tests strengthen a modest association between capillary proximity and CM transcriptional maturation within the P7 section. The result is similar when distance is measured from cell centers, from the sampled CM boundary, or along a tissue-screened graph. It is also stable to changes in the size and position of the local comparison regions. These checks reduce concerns that the finding depends on one geometric definition or a particular grid boundary.

They do not establish that the association is strong, that every CM subtype behaves identically, or that capillaries cause maturation. A tightly matched near–far CM comparison remains inconclusive. Independent hearts are still the most important missing evidence.

## Does the result depend on how distance is measured?

The primary comparison uses 30,485 ventricular CMs, or 30,484 with a reachable graph path. Cells are compared within the same CM subtype and 200-µm region, with flexible adjustment for RNA capture and cell size, missing nuclear area, morphology and cavity/surface distances. All three definitions give similar small positive associations:

| Distance definition | Score SD per SD of nearness | 95% spatial-cluster interval |
|---|---:|---|
{table}

The boundary analysis helps address a specific concern: larger CMs can have distant centers even when their surface is close to an EC. Its agreement with the centroid analysis makes this geometric explanation less likely to account for the whole association. Boundary distance is still approximate: it uses 64 sampled points on the deposited CM polygon and capillary territories, rather than independently validated membrane positions. The graph-path analysis avoids some straight-line shortcuts through cavities, but is not a true tissue geodesic or a measure of vascular flow.

All three associations also pass the prespecified cluster-score sign-flip sensitivity: none of 9,999 simulated sign-flip statistics was as extreme, giving a Monte Carlo P of 0.0001 and the same three-test BH q. This is a conditional statistical check using 400-µm spatial clusters, not random relocation of vessels or replication across hearts. Residual correlations across clusters and unmeasured anatomy remain possible.

## Could one region, one subtype or local tissue density explain it?

For each distance definition, I repeated the comparison using 100-, 200- and 400-µm regions and four grid positions at each size. All 36 estimates and their clustered intervals remain positive. This reduces concern that the initial result was created by the chosen 200-µm grid. Removing any single 400-µm region from the fitted partial-regression sums also leaves a positive estimate. This last check is an influence diagnostic with nuisance residuals held fixed, not independent held-out validation.

The association persists when capillary sources are defined from broad EC support plus public capillary labels without requiring nuclear subtype agreement: {source.beta_SD:+.3f} SD, with an interval of {source.CI_low:+.3f} to {source.CI_high:+.3f}. It therefore does not rely entirely on the stringent nuclear-supported capillary subset.

Additional conditional models include local cell density, capillary-source capture and capillary counts, or distances to fibroblasts, noncapillary ECs and pericytes. The capillary-nearness coefficient remains positive. With local density, source capture and capillary count held similar, it is {density.beta_SD:+.3f} SD ({density.CI_low:+.3f} to {density.CI_high:+.3f}). These comparisons make a purely generic relationship with cellular crowding less convincing. They do not show that capillaries act independently of the rest of the vascular or stromal environment; some of these controls can represent parts of the same biological relationship.

Subtype-specific point estimates are positive for working ventricular, Ankrd1-associated, Ddc-associated, cycling and Slit2-associated CMs. Several subtype intervals cross zero, especially in smaller or spatially restricted groups. A significant result in one subtype and a nonsignificant result in another does not establish a difference between subtypes. These analyses are secondary and do not justify a claim that the relationship is universal or stronger in one particular subtype.

## What does the stricter local match show?

I compared CMs within 25 µm of a QC-supported capillary center with CMs at least 50 µm away. A pair had to share CM subtype, morphology and 200-µm region, lie within 100 µm of each other, and have closely similar capture and cell area. Matching did not use the maturation score, and no CM was reused. This yielded 916 pairs across 62 spatial clusters.

The mean 21-gene difference is {matched.standardized_difference:+.3f} SD, with a spatial-block interval of {matched.CI_low:+.3f} to {matched.CI_high:+.3f}. Deposited RNA and equal-molecule scores also have positive estimates with intervals crossing zero. This comparison is compatible with the small regression estimate, but it cannot independently confirm a near–far difference. It uses a smaller selected subset and a specific distance contrast, so its wider uncertainty should be retained rather than explained away or hidden.

## What would strengthen the biological conclusion most?

First, repeat the analysis in independent P7 hearts and additional early stages using the same frozen score and distance definitions. The heart, rather than the number of cells or image fields, should be the biological replication unit. This is the most important improvement and cannot be supplied by further resampling of this section.

Second, annotate compact and trabecular myocardium, ventricular chamber, papillary regions and large-vessel surroundings independently from morphology, preferably without viewing maturation scores. Our current local comparisons and morphology proxies reduce regional confounding, but do not recover the original authors' exact anatomical masks. Prespecified analyses within these anatomical regions would make a regional developmental explanation easier to distinguish from a more local vascular relationship.

Third, validate the association against protein and CM function. Sarcomere organization, calcium handling and contraction would test whether the transcriptional association accompanies functional development. Perfusion tracing and direct oxygen measurements would characterize vascular support. PTPRM localization and cell-specific perturbations could then test a candidate mechanism. Neither mechanism should be inferred from the distance association itself.

The strengthened statement is: At P7, capillary proximity has a small association with ventricular CM transcriptional maturation that is robust to several distance definitions, local comparison scales and source-confidence rules. A tightly matched near–far comparison remains inconclusive, and independent-heart replication is needed.

## Methods and reproducibility

This follow-up retains the prior P7 strict ventricular CM cohort, 21-gene frozen-reference Xenium score, corrected RNA views and capillary definitions. No new segmentation, imputation, cell taxonomy, samples or causal inference was introduced. Core RNA is a confidence-supported assignment layer, not nuclear-only truth; anucleate CM sections remain eligible. The three primary geometry tests were chosen before examining their results. Other source definitions, grid, subtype, density and matching analyses are reported as sensitivities, without selecting favorable results as additional primary tests.

The score and exact original author masks/graph remain unavailable as author objects; our score is a directly measured gene-Z mean rather than an exact recreation of AddModuleScore. Regression uses flexible capture/area splines, within-region and CM-type centering, and 400-µm clustered uncertainty. Sign-flip inference uses null-imposed clustered partial-regression scores and assumes approximate cluster independence and symmetry; it does not preserve or randomize vessel topology. Cell-center and boundary distances are two-dimensional measurements. Graph distance sums centroid edge lengths on the previously tissue-screened Delaunay graph. Capillary count and general tissue controls define narrower conditional questions and are kept separate from the primary proximity estimate.

The local match uses fixed thresholds of at most 25 µm versus at least 50 µm from capillary centers, a 100-µm maximum pair separation, the same 200-µm grid/CM state/morphology, and at most 0.25 full-cohort SD differences in log capture and log area. Pair intervals use 1,999 resamples of 400-µm spatial clusters. They describe within-section uncertainty and do not substitute for independently sampled animals.
'''
 (OUT/'Report_text.md').write_text(text)
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 fig,ax=plt.subplots(2,2,figsize=(11,9),layout='constrained');blue='#247ba0';orange='#cd782c'
 def forest(axis,frame,labels,x='beta_SD'):
  for i,row in enumerate(frame.itertuples()):
   v=getattr(row,x);axis.errorbar(v,i,xerr=[[v-row.CI_low],[row.CI_high-v]],fmt='o',color=blue,capsize=3)
  axis.axvline(0,color='.6',lw=.8);axis.set_yticks(range(len(labels)),labels);axis.invert_yaxis()
 forest(ax[0,0],primary,['Cell centers','Sampled CM boundary','Tissue-screened path']);ax[0,0].set(title='a  Distance definitions agree',xlabel='Score SD per SD of nearness')
 grid=r[r.analysis.eq('grid_sensitivity')]
 for i,e in enumerate(primary.exposure):
  z=grid[grid.exposure.eq(e)];ax[0,1].scatter(z.beta_SD,np.full(len(z),i)+np.linspace(-.17,.17,len(z)),color=blue,s=22)
 ax[0,1].axvline(0,color='.6',lw=.8);ax[0,1].set(yticks=range(3),yticklabels=['Cell centers','CM boundary','Tissue path'],xlabel='Coefficient across 12 grid choices',title='b  Stable across comparison grids');ax[0,1].invert_yaxis()
 z=r[r.analysis.eq('CM_subtype')&r.exposure.eq('near_QC_capillary')&r.outcome.eq('core_maturation21')]
 forest(ax[1,0],z,z.CM_subtype.map({'v.CM':'Working CMs','p.CM':'Cycling CMs','Ankrd1.CM':'Ankrd1 CMs','Slit2.CM':'Slit2 CMs','Ddc.CM':'Ddc CMs'}).tolist());ax[1,0].set(title='c  Subtype estimates vary in precision',xlabel='Within-subtype coefficient (SD)')
 z=pairs[pairs.outcome.isin(['core_maturation21','native_maturation21','core_depth50'])]
 labels=[{'core_maturation21':'Core RNA','native_maturation21':'Deposited RNA','core_depth50':'Equal-molecule score'}[v] for v in z.outcome]
 forest(ax[1,1],z,labels,x='standardized_difference');ax[1,1].set(title='d  Local matching remains inconclusive',xlabel='Near minus far score (SD)')
 fig.suptitle('P7: a modest vascular association survives geometric checks',fontsize=15)
 fig.savefig(OUT/'Capillary_maturation_integrated_figure.png',dpi=300);fig.savefig(OUT/'Capillary_maturation_integrated_figure.pdf');plt.close(fig)
 (OUT/'figure_legend.txt').write_text('All panels use the previously frozen P7 ventricular cohort. (a) Three primary distance definitions give similar positive adjusted associations. (b) Estimates remain positive across three region sizes and four grid positions. (c) Subtype estimates are positive, but several intervals include zero; this does not test differences between subtypes. (d) A fixed-threshold local match yields 916 disjoint pairs; the 903 depth-eligible pairs use a fixed-50-molecule score. All three matched-score intervals include zero. Regression intervals use 400-µm spatial clusters; matched intervals use clustered resampling. No panel establishes independent-heart replication, physiological maturity or a causal capillary effect.')
 (OUT/'README.md').write_text('# Strengthening the P7 capillary–CM association\n\n[Word report](Capillary_maturation_retest.docx) · [PDF report](Capillary_maturation_retest.pdf) · [Integrated figure](Capillary_maturation_integrated_figure.pdf) · [Full narrative](Report_text.md)\n\nThe small conditional association survives alternative geometric definitions, grids and source criteria. The local near–far matched comparison remains inconclusive. One P7 section supplies no independent-heart replication. All analyses, including less supportive results, are retained.\n\nRun analysis/strengthen_capillary_association.py and analysis/report_capillary_association_robustness.py with the scientific Python environment. The existing document writer is reused with its OUT variable pointing to this result directory.\n')
if __name__=='__main__':main()
