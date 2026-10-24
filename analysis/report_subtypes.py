from pathlib import Path
import json,pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape,A4
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/cell_identities_subtypes';cfg=json.loads((ROOT/'analysis/subtype_criteria.json').read_text())
s=pd.read_csv(out/'subcluster_identity_evidence.tsv',sep='\t');b=pd.read_csv(out/'reviewed_broad_identities.tsv',sep='\t');counts=pd.read_csv(out/'subtype_counts_and_confidence.tsv',sep='\t')
lines=['# Cell identities and within-lineage subclustering','',
'All 73,052 QC-retained nuclei are accounted for. Broad identities were reviewed with multi-gene evidence, followed by independent subclustering of cardiomyocytes, fibroblasts, endothelial, mural and immune populations. Epicardial and neural/glial populations retain broad labels; a lineage-weak group remains unresolved.',
'', '## Main results','', '| Parent | Nuclei | Computational subclusters |','|---|---:|---:|']
for parent,g in s.groupby('parent',sort=False):lines.append(f"| {parent} | {int(g.nuclei.sum()):,} | {len(g)} |")
lines += ['', 'There are 59 computational subclusters, not 59 validated cell types. Many subclusters share a subtype-like marker profile and differ by age or cycling state. The primary Leiden resolution is 0.6, with 0.3 and 1.0 comparisons provided. Each lineage has its own PCA, neighbor graph and UMAP; coordinates between lineage-specific UMAPs are not directly comparable. No age regression or batch integration was performed.',
'', '| Parent | Proposed phenotype | Nuclei passing all evidence rules | Candidate/review nuclei | Unresolved nuclei |','|---|---|---:|---:|---:|']
for (parent,name),g in counts.groupby(['parent','proposed_subtype'],sort=False):
 sums=g.groupby('confidence').nuclei.sum()
 lines.append(f"| {parent} | {name} | {int(sums.get('provisional-supported',0)):,} | {int(sums.get('candidate/review',0)):,} | {int(sums.get('unresolved',0)):,} |")
lines += ['', 'Review flags request identity review and are not additional QC exclusions. Names explicitly containing candidate retain a review flag even if their marker/library criteria pass. Confidence describes prespecified evidence rules, not calibrated classification probability or independent biological validation. Candidate/review groups often retain coherent parent or chamber markers but fail finer enrichment or library-support requirements. For example, many ventricular-like groups share chamber markers with the rest of the CM population; lack of relative enrichment does not negate CM identity.',
'', '## Criteria used for identity','',
'1. Broad identity: at least three canonical lineage markers with mean log1p(counts normalized to 10,000) >=0.1 and >=10% of nuclei having >=1 corrected count; at least two markers permits only candidate status. Competing lineage panels must themselves have at least two supported markers. Minimum standardized panel score 0.3 and separation 0.25; an unqualified shared ion-channel gene cannot override a coherent lineage panel.',
'2. Proposed subtype: at least two supported panel genes plus the required anchor genes listed below. Detection is defined using >=1 corrected count, not any tiny positive DecontX residual. Positive evidence and exclusion-marker expression are recorded separately. Exclusion genes contribute a score penalty/check; they are not required to be completely absent.',
'3. Evidence status: >=50 nuclei; >=2 libraries each contributing >=10 nuclei; <=80% from one library; positive panel score >=0.3; score margin >=0.25; mean panel enrichment >=0.1 log units versus the remainder of the parent; positive within-library contrasts in >=2 libraries. These are reproducible analysis thresholds, not universal diagnostic cutoffs. They do not establish animal-level independence.',
'4. Cycling is a separate state: a nucleus is cycling-marker positive when >=2 of Mki67, Top2a, Cenpf, Ccnb1, Ube2c and Cdk1 have >=1 corrected count; a subcluster is cycling-enriched when >=10% meet this definition. RNA positivity does not demonstrate completed cell division.',
'5. Ddc-, Slit2- and Ankrd1-enriched CM programs are reported as separate marker states, not mutually exclusive anatomical cell types. An anchor, another supported program gene, and anchor enrichment >0.1 log units are required.',
'', '## Positive markers, anchors and interpretation','', '| Parent | Proposed identity tested | Positive marker panel | Anchor requirement | Interpretation limit |','|---|---|---|---|---|']
for parent,defs in cfg['subtypes'].items():
 for d in defs:lines.append(f"| {parent} | {d['name']} | {', '.join(d['markers'])} | >={d['anchor_min']} of {', '.join(d['anchors'])} | {d['note']} |")
lines += ['', 'Not every tested identity was detected. No supported Col8a1-enriched or Pi16/Dpt-exclusive fibroblast subtype, nodal-like CM or AV-like CM is claimed under this rule set. Their panels remain in the criteria file so lack of support is transparent.',
'', '## Important limits','',
'All labels remain provisional RNA-based interpretations. The conduction-associated CM candidates are not asserted to be anatomically or functionally verified Purkinje cells. CAR_3 is Etv1/Gja5-enriched with limited Cntn2 detection; CAR_9 has stronger Cntn2/Etv1 support. Kirrel3-enriched fibroblasts are not assigned valve location without spatial validation. Developmental Postn-associated remodeling is not equated with pathological fibrosis. Lyve1 alone does not establish lymphatic EC identity; the lymphatic candidate requires endothelial lineage plus Prox1/Flt4 and a broader lymphatic panel.',
'', 'Two libraries per age and differential QC losses (especially P7) limit abundance and developmental inference. Candidate panels can be incomplete and neighboring programs can overlap. Low detection does not prove biological absence. Cell/nucleus proportions do not equal tissue cell proportions. Descriptive marker rankings are not replicate-aware differential-expression significance tests. Spatial and/or independent reference validation is needed before treating anatomical subtype names as established.',
'', '## Data and figures','',
'- `all_cell_identity_and_subtype_annotations.csv.gz`: every retained nucleus, reviewed broad identity, subtype proposal, evidence status, global/subcluster IDs, age/library, QC, cycling and CM marker-state flags.',
'- `subcluster_identity_evidence.tsv`: each subcluster, supporting markers/anchors, library representation, positive contrasts, age composition, confidence and an explicit review reason.',
'- `subtype_identity_criteria.tsv` and `analysis/subtype_criteria.json`: full definitions, anchors, exclusion checks and thresholds.',
'- `subtype_marker_evidence.tsv`: per-subcluster mean expression and raw/corrected detection fractions for every tested marker.',
'- `candidate_panel_comparisons.tsv`: every candidate panel score, gate result and within-library contrast, including rejected candidates.',
'- `subcluster_descriptive_markers.tsv`: unbiased descriptive marker ranking within each parent, excluding mitochondrial/ribosomal genes.',
'- Parent directories contain UMAPs, marker dot plots, per-nucleus assignments, resolution labels, sample/age cross-tabs, all-gene subcluster pseudobulk counts and metadata-only H5AD objects. Metadata-only H5AD files contain no expression matrix.',
'- The full annotated count H5AD, with original integer counts and corrected fractional counts, is prepared for a GitHub Release; see the repository root for upload status.',
'', '## Broad identity review','',
'The review assigns 617 previously mixed nuclei to an endothelial profile and 178 to a pericyte profile. The 157-nucleus neural/glial group retains candidate status because of limited marker support; 103 lineage-weak nuclei remain unresolved. Earlier preliminary counts and lineage profiles used a different marker-based labeling pass and are not silently overwritten.',
'', '## Sources','',cfg['sources']['paper'],cfg['sources']['authors_annotation_code'],
'', 'The authors explicitly used additional Xenium spatial information for final annotations. Their numeric cluster labels were not transferred to this independently processed dataset.']
for parent in ['Cardiomyocyte','Fibroblast','Endothelial','Mural','Immune']:
 lines += ['', '## '+parent+' figures','', f'![{parent} subcluster UMAP]({parent}/subcluster_umap.png)','',f'![{parent} marker evidence]({parent}/subtype_marker_dotplot.png)']
(out/'README.md').write_text('\n'.join(lines)+'\n')
styles=getSampleStyleSheet();story=[Paragraph('Cell identities and subtype evidence',styles['Title']),Spacer(1,10)]
for paragraph in [lines[2], '59 computational subclusters were evaluated. Numerical clusters are not automatically distinct cell types. RNA-based identities and anatomical interpretations remain provisional.']:
 story += [Paragraph(paragraph,styles['BodyText']),Spacer(1,9)]
rows=[['Parent','Nuclei','Subclusters']]+[[parent,f'{int(g.nuclei.sum()):,}',str(len(g))] for parent,g in s.groupby('parent',sort=False)]
t=Table(rows,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.lightgrey),('FONTSIZE',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story += [t,Spacer(1,12)]
for label in ['Criteria used for identity','Important limits']:
 story.append(Paragraph(label,styles['Heading1']));start=lines.index('## '+label)+1
 for line in lines[start:]:
  if line.startswith('## '):break
  if line:story.extend([Paragraph(line.replace('&','&amp;').replace('>=','&gt;=').replace('<=','&lt;='),styles['BodyText']),Spacer(1,8)])
story.extend([PageBreak(),Paragraph('Subtype identity criteria',styles['Heading1'])])
for parent,defs in cfg['subtypes'].items():
 story.append(Paragraph(parent,styles['Heading2']))
 rows=[[Paragraph('Proposed identity',styles['BodyText']),Paragraph('Positive marker panel',styles['BodyText']),Paragraph('Required anchors',styles['BodyText'])]]
 for d in defs:
  rows.append([Paragraph(d['name'],styles['BodyText']),Paragraph(', '.join(d['markers']),styles['BodyText']),Paragraph('At least '+str(d['anchor_min'])+' of '+', '.join(d['anchors']),styles['BodyText'])])
 table=Table(rows,colWidths=[160,340,240],repeatRows=1,hAlign='LEFT');table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.lightgrey),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),6)]));story.extend([table,Spacer(1,10)])
for parent in ['Cardiomyocyte','Fibroblast','Endothelial','Mural','Immune']:
 story.extend([PageBreak(),Paragraph(parent+' subclustering',styles['Heading1']),Image(str(out/parent/'subcluster_umap.png'),width=740,height=222),Spacer(1,12)])
 # Dotplot height varies with the number of marker columns and rows.
 from PIL import Image as PILImage
 im=PILImage.open(out/parent/'subtype_marker_dotplot.png');w,h=im.size
 story.append(Image(str(out/parent/'subtype_marker_dotplot.png'),width=740,height=min(260,740*h/w)))
SimpleDocTemplate(str(out/'cell_identity_subtype_report.pdf'),pagesize=landscape(A4),leftMargin=45,rightMargin=45,topMargin=35,bottomMargin=35).build(story)
print('README and PDF report saved.')
