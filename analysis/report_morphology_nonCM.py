"""Biological synthesis, source anatomy tables and integrated association plate."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/anatomy-mpl')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from morphology_nonCM_associations import assign,STAGES,BROAD
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/morphology_stratified_nonCM'
def main():
 r=pd.read_csv(OUT/'nonCM_maturation_associations.tsv',sep='\t');co=pd.read_csv(OUT/'anatomical_CM_coverage.tsv',sep='\t')
 cov=pd.read_csv(OUT/'source_coverage_and_assessability.tsv',sep='\t');prim=r[r.primary_analysis]
 def one(**kw):
  z=r
  for k,v in kw.items():z=z[z[k].eq(v)]
  assert len(z)==1,(kw,len(z));return z.iloc[0]
 def fmt(z):return f'{z.beta_SD:+.3f} ({z.CI_low:+.3f} to {z.CI_high:+.3f})'
 records=[];distribution=[]
 for s in STAGES:
  f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{s}/atlas_annotations.tsv.gz',sep='\t',index_col=0)
  fb=pd.read_csv(ROOT/f'results/fibroblast_multimethod_identity/{s}_multimethod_FB_annotations.tsv.gz',sep='\t',index_col=0)
  q=f[f.atlas_identity.isin(BROAD)].copy();a=assign(s,q[['x_centroid','y_centroid']].to_numpy());a.index=q.index
  q=q.join(a);q['source_subtype']=q.atlas_identity
  ec=q.atlas_identity.eq('Endothelial');q.loc[ec,'source_subtype']='EC:'+q.loc[ec,'EC_nuclear_type'].fillna('Unresolved')
  ids=q.index[q.atlas_identity.eq('Fibroblast')];q.loc[ids,'source_subtype']='FB:'+fb.recommended_candidate.reindex(ids).fillna('Unresolved')
  q['stage']=s;q.index.name='cell_id'
  records.append(q[['stage','x_centroid','y_centroid','atlas_identity','source_subtype','program_source_eligible','ROI','region','territory','confidence']])
  for (region,typ),g in q[q.program_source_eligible].groupby(['region','source_subtype']):distribution.append(dict(stage=s,region=region,source_subtype=typ,source_profiles=len(g)))
 pd.concat(records).to_csv(OUT/'nonCM_anatomical_assignments.tsv.gz',sep='\t',index=True)
 pd.DataFrame(distribution).to_csv(OUT/'nonCM_region_source_counts.tsv',sep='\t',index=False)
 counts=[]
 for s in STAGES:
  z=co[co.stage.eq(s)];counts.append([s,int(z.CMs.sum()),int(z.loc[z.region.eq('Compact wall'),'CMs'].sum()),int(z.loc[z.region.eq('Cavity-facing myocardium'),'CMs'].sum())])
 table='\n'.join(f'| {s} | {n:,} | {compact:,} | {cavity:,} |' for s,n,compact,cavity in counts)
 early=one(stage='P0',source='FB:Postn',scope='Annotated_regions',metric='path',outcome='Y__core_maturation21')
 cap=one(stage='P7',source='EC:Capillary-like',scope='Annotated_regions',metric='path',outcome='Y__core_maturation21')
 vein=one(stage='P14',source='EC:Venous-like',scope='Annotated_regions',metric='path',outcome='Y__core_maturation21')
 elastic=one(stage='P7',source='program__Elastic_matrix',scope='Measured_FB_programs',metric='mean_program50',outcome='Y__core_maturation21')
 text=f'''# Cardiac anatomy, neighboring cell subtypes and CM transcriptional development

## What the new analysis adds

The new analysis separates recognizable wall and cavity-facing regions before asking how nearby non-CMs relate to CM development. It retains a small positive P7 capillary estimate and identifies a P0 fibroblast association in conservative compact-wall patches. The relationships change with age, CM program and neighborhood definition. They do not support a single rule that one neighboring cell type uniformly accelerates CM maturation.

Anatomical certainty is a material limit. I annotated the available DAPI images without displaying maturation scores, RNA expression or cell-type maps during ROI construction, and froze the resulting polygons before the new model fits. The annotations identify conservative continuous compact-wall patches and visible cavity-facing myocardium. They do not establish definitive left/right ventricular masks or verified papillary muscles. Cavities are partly collapsed or obliquely cut, and the acquired antibody channels are no longer cached. This is an outcome-independent candidate annotation, not an independently blinded expert histology assessment or recovery of the authors' masks.

## What can be recognized from morphology?

There are 31 candidate ROIs across the four sections: 16 compact-wall patches, seven cavity-facing patches, two possible papillary/trabecular projections and six vessel-like surroundings. The latter two categories have low confidence and are excluded from primary tests. Only the compact and cavity-facing patches enter the main anatomical analysis (Figure 1). They contain 36,647 frozen ventricular CM profiles, about 31% of the original strict cohort; the remaining profiles are retained as unclassified or low-confidence rather than forced into anatomical categories.

| Age | Original strict CMs | Compact-wall CMs | Cavity-facing CMs |
|---|---:|---:|---:|
{table}

The cavity-facing label means visible chamber-facing or ridged myocardium. It should not be treated as a verified histological trabecular layer. At P0, no reliable cavity-facing patch was designated. Image-coordinate territories are retained for review, but anatomical LV/RV names remain unresolved. The possible papillary projections cannot be distinguished confidently from trabecular ridges in these two-dimensional images. Closed, DAPI-poor lumens define vessel-like candidates, but do not establish artery or vein identity. These small candidate regions also contain too few independent spatial clusters for reliable region-specific association tests.

A self-contained morphology-only viewer accompanies the report. It displays the DAPI images and ROI outlines without RNA or maturation scores, supports vertex editing and saves reviewed polygons as GeoJSON. This makes the anatomical interpretation concrete and reviewable. The current statistical results use the frozen version, not any subsequently edited copy.

## Fibroblasts: the strongest early relationship is broader than a subtype-specific effect

At P0, CMs closer to fibroblasts in the annotated compact-wall patches have higher 21-gene scores. The Postn-associated estimate is {fmt(early)} SD, with q={early.primary_BH_q:.4f} across the 48 primary comparisons. The metabolic-regulatory component is also higher. Deposited RNA and cell-center distances give the same direction. Targeted score checks retain the positive association after omitting Myl2, jointly omitting Cebpb/Rxrg/Lpin1, and equalizing measurement to 50 molecules per CM. This reduces concern that a single panel gene or capture depth entirely creates the result.

However, 1,946 of the 2,020 eligible P0 fibroblasts are Postn-associated. The broad-fibroblast and Postn estimates are therefore nearly identical. This is evidence for a relationship with the early Postn-dominated fibroblast environment; it cannot establish that Postn fibroblasts have a unique effect compared with other P0 fibroblast subtypes. It also does not demonstrate that fibroblasts cause CM maturation. Local matrix organization and shared regional development remain possible explanations.

At P7, the overall-score relationships with Postn-associated and Gfpt2/Pi16-associated fibroblast proximity are uncertain. Later Col8a1-associated and Gfpt2/Pi16-associated proximity also fails to provide a corrected overall-score association. Lsamp-associated candidates are tested where spatial support permits; Kirrel3 candidates lack adequate local support in these annotated patches. No general maturation-promoting role is established. Col8a1 remains one family, and Gfpt2/Pi16 is a multi-gene interpretation rather than a mandatory double-positive gate. Missing same-age reference families and inadequate local coverage are recorded as not assessable, not biological absence.

Measured fibroblast programs provide a complementary result. At P7, CMs surrounded by stronger Eln/Fbln1/Fbln2 expression have slightly lower overall maturation scores: {fmt(elastic)}, secondary q={elastic.secondary_BH_q:.3f}. The direction persists in deposited CM RNA and after controlling for informative fibroblast count and RNA capture. This is an exploratory matrix-context relationship. The genes do not measure elastin deposition, tissue stiffness or an outflow-tract location, and the finding does not identify a particular fibroblast family as causal.

## Endothelial cells: vascular proximity remains plausible, but maturation is not one uniform response

Within the annotated P7 regions, the capillary estimate remains positive: {fmt(cap)} SD. Its size is similar to the preceding within-region retest. It is less precise in this smaller anatomical subset, and q={cap.primary_BH_q:.3f} after correcting across the full 48-comparison source screen. Deposited RNA has a positive estimate with an interval crossing zero. A larger fraction of local ECs being capillary-like does not show the same clear relationship. Nearness and relative neighborhood composition ask different questions; neither should be substituted for the other to claim a stronger result.

Thus, adding morphology has not reversed the P7 capillary estimate, but the new screen does not independently confirm it across multiple animals or establish a universal vascular effect. Gene-omission checks retain a positive estimate, while equal-molecule scoring in this smaller anatomical subset gives a positive but uncertain estimate. The earlier targeted analysis and this broader screen also use different cohorts and multiplicity families. A change in statistical significance alone is not evidence that regional anatomy explains the association.

P14 venous-like EC proximity gives another candidate: {fmt(vein)}, q={vein.primary_BH_q:.3f}. It does not pass primary correction, and the local venous fraction is not clearly associated with the overall score. Arterial-like, endocardial-like and cycling-associated EC proximity do not provide corrected overall-score findings in the primary screen. Lymphatic-like cells and some stage-specific EC classes lack sufficient local support and cannot be assessed reliably.

Endocardial-like ECs illustrate why the individual CM programs matter. In the annotated P14 cohort, nearness to these cells relates to lower sarcomere RNA, including within Ankrd1-associated CMs. This differs from the earlier positive association using an endocardial neighborhood-composition measure in a broader cohort. The present fraction-based comparison is uncertain. The earlier positive relationship should therefore not be treated as an established general effect of endocardial ECs on CM structural maturation. Distance, composition, regional selection and technical definitions need to be reconciled before assigning a biological direction to that relationship. Neither result establishes greater overall CM maturity.

## Mural, immune and epicardial populations

Pericytes, smooth-muscle cells, immune cells and epicardial cells are included as supported broad sources. None provides an overall-score association surviving the primary correction. Some component-specific or region-specific estimates differ, but they remain secondary hypotheses in one section per age. The current spatial reference does not resolve validated fine immune, pericyte, smooth-muscle or epicardial subtypes, so those labels have not been invented or inferred from location alone.

The most useful biological picture is therefore a changing set of local environments. Early fibroblast-rich matrix settings, capillary proximity and particular matrix-expression programs accompany different parts of CM transcriptional development. The data identify comparisons to test, rather than a ranked list of cells that promote or inhibit maturation. Independent hearts, expert-reviewed regional anatomy and functional CM measurements are needed to establish which associations extend beyond these sections.

## Methods and interpretation limits

The ROI construction script reads only cached DAPI overviews and physical image extents. Frozen GeoJSON and image checksums record the exact annotations. Boundary and alphaSMA/vimentin TIFF caches were absent; Zenodo retrieval failed with a network proxy403 response even after network access was granted. Consequently, no newly recovered antibody images or verified vessel walls are claimed. Previously validated RNA counts, segmentation and cell identities were retained, including valid anucleate CM sections.

The CM outcome is the previously frozen mean of 21 measured Xenium gene Z values, with metabolic-regulatory, sarcomere, excitation–contraction and cycling components separate. It is not the authors' unavailable exact AddModuleScore object, and it does not measure physiological function. Separate deposited and core RNA outcomes use fixed full-stage scales. The primary source metric is minus log1p shortest weighted centroid-path distance on the frozen tissue-screened Delaunay graph; paths are not blood flow, diffusion or verified membrane contact. Cell-center nearness and 50-µm source fractions are retained as alternative questions.

Models compare CMs within image ROI, CM subtype and 200-µm grid, with flexible capture and cell-size controls, missing nuclear area and continuous cavity/surface distances. A source requires at least 30 eligible profiles, at least 100 locally supported CMs, and support in at least eight 400-µm clusters; each fitted stratum also requires adequate local source support across eight clusters. Insufficiently supported subtype/region tests are explicitly skipped. Intervals use clustered CR1 t uncertainty within the observed section. They do not incorporate independent-heart, annotation or segmentation variability.

The 48 primary tests cover supported broad sources, current EC subtypes and current fibroblast families across the four ages, using overall core CM score in the annotated regions. Broad parents and subtypes can overlap substantially; a positive subtype result is not proof of a unique effect beyond its parent lineage. All core component, alternative metric, source-rule, CM-state and program comparisons are corrected together as secondary tests; native RNA is retained as a measurement sensitivity. The measured fibroblast-program model includes informative-source count and capture; a version without those source controls is also retained. Nuclear-supported and public-plus-broad-supported EC source rules are distinguished, and unresolved identities remain in lineage denominators.

The named P0 Postn and P7 capillary comparisons receive targeted follow-up score checks, including Myl2 omission, joint Cebpb/Rxrg/Lpin1 omission and equal-50-molecule scores in the same depth-eligible cells. These 10 sensitivity tests are corrected separately and do not add primary discoveries or independent biological replication.

No new biological samples, segmentation, spatial gene imputation, lineage tracing or causal perturbation were added. Chamber, papillary and vessel-like candidates require anatomical review before they can support definitive compartment-specific conclusions. The exact original author masks remain unrecovered.

Source: Wang et al., Spatially guided in vivo single-cell functional genomics of postnatal heart, Nature Cardiovascular Research, DOI 10.1038/s44161-026-00861-z. https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/ The deposited source datasets are GSE290989 and Zenodo 14961929/14961942; the reviewed reference counts, source definitions and earlier QC remain unchanged.
'''
 (OUT/'Report_text.md').write_text(text)
 # Integrated plate: all primary findings, anatomical sampling and selected sensitivities.
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 fig=plt.figure(figsize=(12,10),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1.2,1]);a=fig.add_subplot(gs[:,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,1])
 order=['Broad:Endothelial','EC:Capillary-like','EC:Endocardial-like','EC:Arterial-like','EC:Venous-like','EC:Lymphatic-like','EC:Cycling-associated',
  'Broad:Fibroblast','FB:Postn','FB:Gfpt2','FB:Col8a1','FB:Lsamp','FB:Kirrel3','Broad:Pericyte','Broad:Smooth muscle','Broad:Immune','Broad:Epicardial']
 mat=prim.pivot(index='source',columns='stage',values='beta_SD').reindex(index=order,columns=STAGES)
 cmap=plt.colormaps['coolwarm'].copy();cmap.set_bad('#dedede');im=a.imshow(np.ma.masked_invalid(mat.to_numpy()),cmap=cmap,vmin=-.09,vmax=.09,aspect='auto')
 labels=[s.replace('Broad:','').replace('EC:','EC · ').replace('FB:','FB · ').replace('Gfpt2','Gfpt2/Pi16') for s in order]
 a.set(yticks=range(len(order)),yticklabels=labels,xticks=range(4),xticklabels=STAGES,title='a  CM score versus source nearness')
 for i in range(len(order)):
  for j in range(4):
   if not np.isfinite(mat.iloc[i,j]):a.text(j,i,'×',ha='center',va='center',color='#777777',fontsize=11)
 for row in prim.itertuples():
  if row.primary_BH_q<.05:a.text(STAGES.index(row.stage),order.index(row.source),'*',ha='center',va='center',fontsize=17)
 fig.colorbar(im,ax=a,shrink=.6,label='Adjusted coefficient (stage SD)')
 vals=np.array([x[1:] for x in counts]);compact=vals[:,1];cavity=vals[:,2];other=vals[:,0]-compact-cavity
 b.bar(STAGES,compact/vals[:,0],label='Compact wall',color='#299b86');b.bar(STAGES,cavity/vals[:,0],bottom=compact/vals[:,0],label='Cavity-facing',color='#e6a13b');b.bar(STAGES,other/vals[:,0],bottom=(compact+cavity)/vals[:,0],label='Unclassified/candidates',color='#c5c9ce')
 b.set(ylabel='Fraction of strict ventricular CMs',title='b  Image-region coverage',ylim=(0,1));b.legend(frameon=False,fontsize=8,loc='upper right')
 choices=[('P0','FB:Postn','Annotated_regions','Postn-rich fibroblasts, P0'),('P7','EC:Capillary-like','Annotated_regions','Capillary ECs, P7'),('P14','EC:Venous-like','Annotated_regions','Venous-like ECs, P14')]
 for i,(stage,source,scope,label) in enumerate(choices):
  for j,view,col in [(0,'core','#247ba0'),(1,'native','#cf7529')]:
   z=one(stage=stage,source=source,scope=scope,metric='path',outcome='Y__'+view+'_maturation21');y=i+(j-.5)*.16
   c.errorbar(z.beta_SD,y,xerr=[[z.beta_SD-z.CI_low],[z.CI_high-z.beta_SD]],fmt='o',color=col,capsize=3,label=view if i==0 else None)
 c.axvline(0,color='.6',lw=.8);c.set(yticks=range(3),yticklabels=[v[3] for v in choices],xlabel='Adjusted overall-score coefficient (SD)',title='c  CM score in two RNA views');c.invert_yaxis();c.legend(frameon=False,fontsize=9)
 fig.suptitle('CM development relates to selected vascular and matrix settings',fontsize=15)
 fig.savefig(OUT/'figures/Figure_2_nonCM_relationships.png',dpi=300);fig.savefig(OUT/'figures/Figure_2_nonCM_relationships.pdf');plt.close(fig)
 legends=[dict(stem='Figure_1_morphology_regions',title='Conservative morphology-only regions',legend='Green polygons mark compact-wall patches; orange marks cavity-facing myocardium; purple marks low-confidence papillary/trabecular projections; pink circles mark vessel-like lumens and surroundings. Uncolored tissue remains unclassified. All polygons were drawn from DAPI without displaying RNA or maturation maps during construction. These are candidate ROIs, not verified LV/RV, trabecular, papillary or arterial masks. Only green and orange patches enter the primary models.'),
  dict(stem='Figure_2_nonCM_relationships',title='Anatomy-conditioned non-CM relationships with CM development',legend='(a) All 48 primary comparisons are shown. Gray entries lack adequate reference/source/local support; stars mark primary BH q<0.05. Broad fibroblasts and Postn fibroblasts overlap almost completely at P0, so the two stars are not independent biological evidence. (b) Conservative usable anatomical patches cover about 31% of the strict CM cohort; unclassified regions and low-confidence candidates are retained outside primary tests. (c) Selected primary estimates and deposited-RNA sensitivity intervals are shown without implying all pass correction. Core is a confidence-supported RNA-assignment layer, not nuclear-only truth. Intervals describe within-section clustered uncertainty; there is one observed section per age.')]
 (OUT/'figure_legends.json').write_text(json.dumps(legends,indent=2)+'\n')
 (OUT/'README.md').write_text('# Morphology, non-CM subtypes and CM development\n\n[Word report](Morphology_nonCM_CM_development.docx) · [PDF report](Morphology_nonCM_CM_development.pdf) · [Narrative](Report_text.md) · [Morphology-only review viewer](Morphology_only_ROI_review.html) · [Frozen editable regions](morphology_ROIs.geojson)\n\nThe annotations are outcome-independent DAPI-based candidates. Definitive LV/RV, papillary and large-vessel identities remain unresolved. The usable compact/cavity-facing patches cover about 31% of strict CMs; no original author masks were recovered. The new subtype screen supports a P0 Postn-dominated fibroblast proximity association, without establishing subtype-specific causation. P7 capillary nearness remains positive but does not pass correction across the broader 48-test family.\n\nAll supported sources, alternative definitions, uncertain identities and unassessable tests are retained. Fine immune/mural/epicardial subtypes are not available and are not invented. The primary CM score is the frozen measured 21-gene score, not physiological maturity or exact author AddModuleScore.\n\nReproduce: annotate_morphology_regions.py; morphology_nonCM_associations.py; check_morphology_score_robustness.py; report_morphology_nonCM.py with the scientific Python environment; morphology_review_viewer.py and write_morphology_nonCM_report.py with base Python; validate_morphology_nonCM.py with scientific Python. acquire_anatomy_images.py is an optional retrieval helper: no antibody-channel download completed in this turn.\n')
if __name__=='__main__':main()
