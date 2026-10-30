"""Combine the uploaded reports in clear biological prose with current FB data."""
from pathlib import Path
import csv,json,os,datetime
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from PIL import Image
import write_cardiac_story as style
from write_cardiac_biological_story import add_blocks

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/heart_atlas_clear_language_revision'
SOURCE=ROOT/'analysis/heart_atlas_clear_language_revision.md'
BIO='results/cardiac_biological_synthesis/figures/';FB='results/fibroblast_multimethod_identity/figures/'
FIGURES=[
('1',BIO+'Figure_1_CM_identity_and_maturation','CMs of the same subtype show different maturation-associated RNA patterns',
 'a,b, Ventricular CM subtypes at P0 and P7. c,d, Relative maturation scores in the same sections, measured using the 21-gene panel. e,f, Scores within each subtype. Cells with the same subtype label can differ in maturation-associated gene expression. A higher score does not directly measure stronger contraction or greater energy production. The labels, scores and positions are unchanged from the earlier analysis. Scale bars, 500 µm.'),
('2','results/heart_atlas_clear_language_revision/figures/Figure_2_fibroblast_matrix_patterns','Fibroblasts develop different matrix-associated features',
 'a, Matrix-associated gene expression across ages and fibroblast groups, averaged across the two RNA libraries and shown relative to the other groups. b, Marker expression in named P21 subgroups; larger dots indicate detection in more nuclei, and color shows relative expression. c, Col8a1 fibroblasts were more often near arterial structures and smooth-muscle cells than matched Gfpt2 fibroblasts in the earlier P14 comparison. This panel retains the earlier conservative labels. d, Selected early and later RNA groups resemble one another even when interpreted markers are withheld. Dashed arrows identify comparisons to follow; they do not show one group becoming another.'),
('3',FB+'Figure_3_stage_maps_and_boundary_stability','Alternative methods recover broader fibroblast candidate maps',
 'a–d, Candidate fibroblast subtypes at P0, P7, P14 and P21. e, The fraction of eligible fibroblasts receiving a candidate identity. More labels give a broader map, but do not by themselves establish greater accuracy or tissue-wide abundance. f, Agreement between whole-cell and inner-cell assignments at P21. Randomly reducing whole-cell RNA to a similar depth helps test whether disagreement reflects fewer transcripts or their position within the cell. Supported alternative nuclear identities remain flagged for review. Gray shows profiles without a confident candidate. Scale bars, 500 µm.'),
('4',FB+'Figure_2_P21_identity_and_program_maps','A Gfpt2-associated pattern extends beyond the narrow subtype',
 'a, The paper’s fibroblast labels among eligible P21 profiles. b, Candidate identities from agreement between classifiers. c, A smaller subset also supported by inner-cell RNA and by tests that hide the naming genes. d, Expression of the measured Gfpt2/Pi16/Adgrd1 pattern relative to control genes. e,f, Cells expressing Gfpt2-associated, Col8a1-associated or overlapping patterns, across all eligible fibroblasts and within the paper’s Gfpt2 group. A shared gene pattern can cross subtype boundaries. Kirrel3 candidates are not assumed to occupy valves. Gray profiles and review flags retain uncertainty. Scale bars, 500 µm.'),
('5',BIO+'Figure_4_vascular_growth_and_CM_specialization','Vascular growth and CM maturation have different relationships over time',
 'a, Coronary growth described by Lu et al. (2021); blood flow was not newly measured here. b, At P7, cycling CMs near Postn-associated fibroblasts showed lower energy-metabolism gene expression, but the relationship weakened with alternative fibroblast and neighborhood definitions. c, At P14, Ankrd1-associated CMs near endocardial-like cells showed higher sarcomere-gene expression. Intervals describe uncertainty within a section. d, Comparable CMs alongside capillary cells showed lower cell-cycle gene expression at P7 and higher energy-metabolism gene expression at P14/P21. This analysis uses a broader CM cohort. e, Growing vessels, being near a vessel and receiving a functioning blood supply are different observations. These earlier estimates have not been refitted with the newest fibroblast labels.'),
('6',BIO+'Figure_5_CM_identity_conditioned_transport','Some cellular environments resemble one another across ages',
 'a, Neighborhoods are compared around the same CM subtype and at a similar distance from the chamber lining. b, Fibroblast and endothelial composition in 96 observed neighborhoods, using the earlier strict assignments. c,d, Agreement between Wormhole and exact calculations of neighborhood similarity. e, Similarity of environments around the same or different CM subtypes within an age. f, Maturation ranks are compared after matching; they were not used to select the matches. Some matched environments contain CMs with similar relative maturation, but this varies by age and region. No individual cells were followed. The newest fibroblast labels have not been incorporated into these comparisons.'),
('7',BIO+'Figure_6_regional_developmental_model','Experiments can separate blood supply, matrix and contact effects',
 'a, A proposed sequence linking early tissue growth with later muscle specialization. Arrows indicate developmental periods, not observed cell conversion. b, Blood delivery, matrix organization and endothelial contact may each contribute to CM development. c, Intermediate-stage hearts and controlled experiments can test whether these contributions occur in a particular order or act together. Measurements should include blood flow, CM energy use, calcium handling and contraction. Lineage tracing addresses cell ancestry. These are experimental predictions, not demonstrated causal relationships.')]
SUPPLEMENT=[
('S1',FB+'Figure_1_reference_recovery_and_sensitivity','The reference tests show why more labels need not mean better assignments',
 'a, P21 Gfpt2 reference tests: recall measures how many Gfpt2 nuclei are recovered, while precision measures how many predicted Gfpt2 nuclei carry that reviewed label. Each point represents one held-out library. b, Recovery of the reviewed families by each method. c, Changes when naming genes are hidden or late-age references are pooled. d, Changes when less RNA information is available. A method can recover more Gfpt2 nuclei while also mislabeling other families. These tests assess the RNA definitions, not independently verified Xenium identities.'),
('S2',FB+'Figure_4_heldout_genes_and_published_programs','Hidden genes help check assignments, while marker patterns can overlap',
 'a, Expression of P21 naming genes after cells were grouped without allowing classifiers to use those genes. Printed values give detection percentages. b, Gfpt2- and Col8a1-associated patterns can overlap among cells carrying the paper’s Gfpt2 label. c, The number of published markers measured in the shared panel. Too few valve markers are measured to score a multi-gene pattern. d, A proposed subtype, a shared gene pattern and a cell’s developmental fate answer different questions.')]

SEGMENTATION_NOTE=[
 'To place a cell subtype on the map, we first need to know which RNA belongs to that cell. The atlas combines membrane, cell-interior and nuclear stains to identify cell outlines. This is especially useful for large cardiomyocytes, whose full shape cannot be determined by drawing a fixed area around the nucleus.',
 'The authors compared the automatic outlines with expert drawings in five P0 and five P21 image fields. Their combined-stain method performed better than nucleus expansion. Agreement was poorer at P21: approximate scores were 45–50% at P0 and around 10% at P21 under a demanding 90% overlap rule. These values describe agreement between cell outlines, not the percentage of RNA assigned correctly. A cell can fail the strict outline test while still retaining many correctly assigned interior transcripts.',
 'This matters for fibroblasts and endothelial cells beside RNA-rich cardiomyocytes. A small boundary error can introduce neighboring muscle RNA and change a subtype assignment. Conversely, a section through a large CM may contain no visible nucleus. Removing every profile without a nucleus would discard some valid muscle-cell information.',
 'We therefore compare whole-cell RNA with RNA from the inner part of the cell and the nucleus. Agreement strengthens an assignment. Sparse nuclear RNA does not by itself disprove the whole-cell label, but nuclear evidence for another cell identity requires review. Each view provides useful information; none proves the true boundary on its own.',
 'The next step is direct imaging of selected cell interfaces, using clear membrane outlines and several subtype markers in the same cells. This would help determine whether overlapping fibroblast patterns occur within individual cells or result from RNA assigned across a boundary. The detailed earlier assessment and the atlas Methods and Supplementary Figure 1 are linked in the source record.'
]

def readrows(path):
    with open(ROOT/path) as f:return list(csv.DictReader(f,delimiter='\t'))
def figure(doc,item):
    number,path,title,legend=item;p=doc.add_heading(f'Figure {number}. {title}',2);p.paragraph_format.page_break_before=True
    picture=ROOT/(path+'.png')
    with Image.open(picture) as im:w,h=im.size
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.space_after=Pt(4)
    p.add_run().add_picture(str(picture),width=Inches(min(7.0,7.10*w/h)))
    p=doc.add_paragraph(legend);p.paragraph_format.space_after=Pt(0)
    for r in p.runs:r.font.size=Pt(9)
def table_page(doc,title):
    p=doc.add_heading(title,1);p.paragraph_format.page_break_before=True
def main():
    OUT.mkdir(parents=True,exist_ok=True);text=SOURCE.read_text();title=text.splitlines()[0][2:]
    doc=Document();doc.core_properties.title=title;doc.core_properties.author='Heart atlas analysis'
    doc.core_properties.subject='Clear biological synthesis incorporating alternative fibroblast identity evidence'
    doc.core_properties.created=datetime.datetime(2026,10,26,12,tzinfo=datetime.timezone.utc);doc.core_properties.modified=doc.core_properties.created
    s=doc.sections[0];s.page_width=Inches(8.5);s.page_height=Inches(11);s.top_margin=Inches(.65);s.bottom_margin=Inches(.65);s.left_margin=Inches(.7);s.right_margin=Inches(.7)
    s.header_distance=Inches(.25);s.footer_distance=Inches(.25);s.different_first_page_header_footer=True
    n=doc.styles['Normal'];n.font.name='Cambria';n.font.size=Pt(11);n.paragraph_format.line_spacing=1.02;n.paragraph_format.space_after=Pt(4)
    for name in ['Title','Heading 1','Heading 2']:
        doc.styles[name].font.name='Arial';doc.styles[name].font.color.rgb=RGBColor.from_string('203C54')
    doc.styles['Heading 1'].font.size=Pt(14.5);doc.styles['Heading 2'].font.size=Pt(11.5)
    s.header.paragraphs[0].text='POSTNATAL HEART DEVELOPMENT | MUSCLE, MATRIX AND BLOOD VESSELS'
    for r in s.header.paragraphs[0].runs:r.font.size=Pt(8)
    footer=s.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.RIGHT;field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');footer._p.append(field)
    doc.add_paragraph('INTEGRATED BIOLOGICAL REPORT • CLEAR-LANGUAGE REVISION','Subtitle');doc.add_paragraph(title,'Title')
    doc.add_paragraph('P0–P21 mouse-heart snRNA-seq and Xenium analysis')
    doc.add_heading('Central message',2)
    doc.add_paragraph('Muscle, matrix and blood vessels develop together, but different aspects of cardiomyocyte maturation emerge at different times. Early muscle and vessel growth may be coordinated; later, vascular support and the surrounding matrix may accompany metabolic and structural development. Fibroblast subtypes share some gene patterns while differing in others.')
    doc.add_paragraph('This report brings together the cardiomyocyte, vascular and fibroblast findings. Each section describes what was observed, what it may mean biologically and what remains uncertain. Seven main figure plates and two supporting plates are included; detailed analysis rules are collected in Methods.')
    doc.add_paragraph('The new fibroblast maps have not yet been used to recalculate the earlier CM-neighbor associations or cross-age transport results. One available section per age limits developmental and causal conclusions.')
    doc.add_paragraph('Prepared 26 October 2026 • Source results and earlier reports are retained for traceability.')
    p=doc.add_heading('Abstract',1);p.paragraph_format.page_break_before=True
    body,methods=text.split('## Methods',1)
    # The first heading is already placed with an explicit page break.
    body=body.split('## Abstract',1)[1].lstrip();add_blocks(doc,body)
    for item in FIGURES:figure(doc,item)
    add_blocks(doc,'## Methods'+methods)
    table_page(doc,'What the evidence supports')
    style.table(doc,['Biological observation','What it suggests','What remains to test'],[
        ['CMs of the same subtype differ in maturation-associated expression','A subtype can contain several developmental states','Energy use, calcium handling and contraction'],
        ['Later fibroblast groups express different matrix genes','They may help organize different tissue properties','Deposited matrix proteins and tissue mechanics'],
        ['Gfpt2-associated expression occurs across several groups','A shared pattern need not define one subtype','Several markers imaged within clear cell boundaries'],
        ['EC and CM cell-cycle activity occur together at P0/P7','Muscle and blood vessels may grow together locally','Completed division and the order of events'],
        ['Capillary-adjacent CMs differ in cell-cycle expression at P7 and metabolism at P14/P21','Early growth and later vascular support may accompany different steps','Blood flow, oxygen delivery and contact signals'],
        ['P14 Ankrd1 CMs near endocardial-like cells express more sarcomere RNA','The structural relationship may differ between CM subtypes','Anatomy, mechanical loading and local signals'],
        ['Some cellular surroundings resemble one another across ages','Specific regions and states are worth following','Independent intermediate-stage hearts and lineage tracing']])
    table_page(doc,'Fibroblast identities and programs: reading the numbers')
    coverage=readrows('results/fibroblast_multimethod_identity/final_candidate_coverage.tsv')
    rows=[]
    for family in ['Postn-associated','Gfpt2-associated','Col8a1-associated','Lsamp-associated','Kirrel3-associated']:
        row=[family.replace('-associated','')]
        for stage in ['P0','P7','P14','P21']:
            r=next(r for r in coverage if r['family']==family and r['stage']==stage);row.append(f'{int(r["candidate_profiles"]):,}')
        rows.append(row)
    style.table(doc,['Candidate family','P0','P7','P14','P21'],rows)
    doc.add_paragraph('Counts refer to candidate assignments among eligible broad fibroblasts. A zero means no candidates under that reference and rule; it is not proof that every cell or marker of that kind is absent. These are not tissue-wide abundance estimates. Kirrel3 candidates are not assigned valve anatomy from their name.')
    style.table(doc,['P21 Gfpt2 evidence','Profiles','Meaning'],[
        ['Marker-pattern enrichment','1,339','Broader measured Gfpt2/Pi16/Adgrd1 signal under the module rule'],
        ['Also enriched for Col8a1 pattern','327','Overlapping RNA patterns; not proof of transition'],
        ['Public Gfpt2 profiles with Gfpt2-pattern flag','1,058 of 2,716','The public-labelled population contains varied marker evidence'],
        ['Naming-marker-hidden majority label','280','Candidate grouping made without the seven naming/anchor genes'],
        ['Cell-supported candidate','85','Reference-screened model agreement without a supported core conflict'],
        ['Cell/core and naming-marker consistency','75','A narrower supported subset; not the total biological population']])
    doc.add_heading('Which results include the new fibroblast labels?',2)
    style.table(doc,['Result','Status in this revision'],[
        ['Fibroblast identity and marker-pattern maps','Updated alternative classifiers and direct programs included'],
        ['CM-neighbor association estimates','Earlier models retained; not refitted with the newest classifiers'],
        ['P14 arterial-proximity comparison','Earlier conservative matched subset retained'],
        ['Cross-age neighborhood matching','Earlier representation and plans retained']])
    for item in SUPPLEMENT:figure(doc,item)
    table_page(doc,'How to read the cell-boundary results')
    for para in SEGMENTATION_NOTE:doc.add_paragraph(para)
    doc.add_paragraph('Sources: Wang et al. (2026), segmentation Methods and Supplementary Figure 1; the earlier detailed segmentation assessment accompanies the source record.')
    name='Postnatal_heart_integrated_clear_language_report';doc.save(OUT/(name+'.docx'))
    inventory=[]
    md=text+'\n\n## Figure plates\n\n'
    for number,path,title,legend in FIGURES+SUPPLEMENT:
        rel=os.path.relpath(ROOT/(path+'.png'),OUT)
        md+=f'### Figure {number}. {title}\n\n![Figure {number}]({rel})\n\n{legend}\n\n'
        inventory.append({'number':number,'source_stem':path,'title':title,'legend':legend,'new_plot':path.startswith('results/heart_atlas_clear_language_revision/')})
    md+='## Reading the cell-boundary results\n\n'+'\n\n'.join(SEGMENTATION_NOTE)+'\n'
    (OUT/'Report_text.md').write_text(md.rstrip()+'\n')
    (OUT/'figure_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    (OUT/'report_structure.json').write_text(json.dumps({'title':doc.core_properties.title,'main_figures':7,'supporting_figures':2,'main_text_words':len(body.split()),
      'new_analysis_models_fitted':False,'new_classifier_labels_used_to_refit_CM_associations':False,'uploaded_reports_integrated':True,
      'style_model':'User-preferred coronary–cardiomyocyte passage: observation, biological interpretation and uncertainty',
      'style_revision_scope':['abstract','introduction','results','discussion','methods','figure legends','cover','evidence tables','segmentation explanation']},indent=2)+'\n')
    print('CLEAR_INTEGRATED_REPORT_CREATED',OUT/(name+'.docx'),flush=True)
if __name__=='__main__':main()
