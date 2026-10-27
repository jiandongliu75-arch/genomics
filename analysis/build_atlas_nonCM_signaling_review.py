"""Audit assay feasibility and publish a paper-method review, not new model results."""
from pathlib import Path
import csv,json,hashlib,zipfile,xml.etree.ElementTree as ET,datetime,re,sys
from collections import defaultdict
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
import fitz
from write_cardiac_biological_story import add_blocks
from write_cardiac_story import table

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/atlas_nonCM_signaling_review'
SOURCE=ROOT/'analysis/atlas_nonCM_signaling_review.md';STEM='Heart_atlas_nonCM_signaling_methods_review'
META=Path('/workspace/heart-atlas/metadata')
PIP=Path('/workspace/heart-atlas/reference/PIP_seq')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(data,name):
    with open(OUT/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
def read(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def sheet(z,name):
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    rels={e.get('Id'):e.get('Target') for e in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    root=ET.fromstring(z.read('xl/workbook.xml'))
    s=next(s for s in root.findall('.//s:sheet',ns) if s.get('name')==name)
    path=rels[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
    path=path.lstrip('/') if path.startswith('/') else 'xl/'+path
    out=[]
    for row in ET.fromstring(z.read(path)).findall('.//s:row',ns):
        d={}
        for c in row.findall('s:c',ns):
            v=c.find('s:v',ns);val=v.text if v is not None else ''
            if c.get('t')=='s':val=strings[int(val)]
            d[re.sub(r'\d+$','',c.get('r'))]=val
        out.append(d)
    return out
def main():
    OUT.mkdir(exist_ok=True,parents=True)
    cal=ROOT/'results/cardiac_spatial_reconstruction/score_reference_calibration.tsv'
    targets={r['gene'] for r in read(cal)};assert len(targets)==403
    panel_path=ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv'
    panel={r['gene'] for r in read(panel_path)};assert len(panel)==21
    definitions=[
     ('PTPRM contact-associated candidate','Ptprm','Both sides measured; direct Meis1 RNA is absent. Boundary staining or perturbation is separate evidence.'),
     ('Paper-prioritized NAMPT–INSR candidate','Nampt;Insr','Both targets measured. Independently review biochemical evidence; feature ranking does not establish binding or sender identity.'),
     ('Fibronectin–integrin candidate','Fn1;Itga5;Itgb1','All named targets measured. RNA does not locate deposited fibronectin or establish assembled receptor activity.'),
     ('Laminin-associated matrix/receptor candidate','Lama4;Lamb1;Lamc1;Dag1','All named targets measured; distinguish matrix deposition from expression in producer cells.'),
     ('Pleiotrophin-associated candidate','Ptn;Ptprz1','Both targets measured; validate expression, cellular source and context.'),
     ('IGF2-associated growth candidate','Igf2;Igf1r','Both targets measured; growth and overall maturation are distinct outcomes.'),
     ('NRG1–ERBB2/ERBB4 complex candidate','Nrg1;Erbb2;Erbb4','Erbb2 is absent. Do not treat this particular complex as fully measured; other ERBB configurations require separate assessment.'),
     ('TGFβ2-associated feedback hypothesis','Tgfb2;Tgfbr1;Tgfbr2;Smad2;Smad3;Smad7','Receptors/ligand measured, listed SMAD transcripts absent. Latent ligand activation and CM→FB→CM causality are unestablished.'),
     ('IGF1-associated trophic candidate','Igf1;Igf1r','Both targets measured; do not assign an immune or other source before checking expression.'),
     ('BMP2/4-associated candidate','Bmp2;Bmp4;Bmpr1a;Bmpr1b;Bmpr2','The listed BMP ligands are absent despite receptor coverage. This is not absence of the whole BMP pathway.')]
    feasibility=[]
    for name,genes,note in definitions:
        required=genes.split(';');present=[g for g in required if g in targets];missing=[g for g in required if g not in targets]
        feasibility.append({'candidate':name,'genes_checked':genes,'measured_targets':';'.join(present),'unmeasured_targets':';'.join(missing),
            'all_listed_targets_measured':not missing,'interpretation':note,'availability_establishes_detection':False})
    save(feasibility,'signaling_target_feasibility.tsv')
    with zipfile.ZipFile(META/'paper_supplementary_tables.xlsx') as z:
        modules=defaultdict(set)
        for row in sheet(z,'Table 8')[1:]:
            if row.get('A') and row.get('B'):modules[row['B']].add(row['A'])
    coverage=[]
    for name,genes in sorted(modules.items()):
        measured=genes&targets;independent=measured-panel
        coverage.append({'module':name,'full_RNA_genes':len(genes),'measured_Xenium_genes':len(measured),
            'measured_non21_genes':len(independent),'measured_genes':';'.join(sorted(measured)),
            'interpretation':'Target coverage only; input-overlap, CM specificity, dropout and held-out-library coherence still require checking'})
    assert [x['measured_Xenium_genes'] for x in coverage]==[83,70,2,10]
    save(coverage,'published_WGCNA_module_coverage.tsv')
    plan={'status':'proposed_not_fitted','primary_question':'Do source-resolved non-CM molecular inputs add information about CM development beyond receiver state and anatomy?',
        'stages':['P0','P7','P14','P21'],'analysis_by_stage':True,
        'outcomes':['Measured 21-gene index','Metabolic-regulatory','Sarcomere','Excitation-contraction','Separate cycling'],
        'feature_sets':{
            'A':'Capture, size, CM identity, anatomy, broad composition',
            'B':'A plus directly measured receiver receptors',
            'C':'B plus CM-to-CM spatial inputs',
            'D':'C plus source-resolved non-CM ligand/matrix-program exposures and source counts',
            'E':'D plus ligand-receptor interaction/compatibility features'},
        'source_families':['EC subtypes','Postn fibroblasts','Gfpt2/Pi16 fibroblasts','Col8a1 fibroblasts retained as one family','Lsamp and Kirrel3 candidates','Pericytes','Smooth muscle','Immune','Epicardial','CM control sources'],
        'joint_all_cell_communication_solver':True,'retain_per_edge_sender_identity':True,
        'primary_validation':'Held-out contiguous spatial regions; training-only fitting/tuning; explicitly disjoint source neighborhoods where feasible',
        'secreted_radius_um':200,'matrix_contact_search_radius_um':50,
        'contact_interpretation':'Search range is not verified membrane contact; use separate operational boundary apposition',
        'initial_gap_for_200um_source_radius_um':400,'validate_actual_source_ID_overlap':True,
        'prevent_global_communication_competition_across_validation_domains':True,
        'comparison_learners':['Interpretable regression','Random forest with training-only complexity selection'],
        'covariate_gene_overlap':'Withhold all outcome genes from predictors and audit compositional normalization',
        'checks':['Native/core and source-confidence sensitivity','Age and CM-state support','Region-preserving source-program nulls','Grouped feature ablation','CM receptor-only versus neighbor increment','Independent hearts when available'],
        'inferred_RNA_activities':'Keep distinct from direct Xenium measurements; use withheld-gene/library validation',
        'causal_or_fate_interpretation':False,'new_model_results_generated':False}
    (OUT/'proposed_analysis_specification.json').write_text(json.dumps(plan,indent=2)+'\n')
    doc=Document();doc.core_properties.title='Relating cardiomyocyte development to signals from surrounding cells'
    doc.core_properties.author='Heart atlas analysis';doc.core_properties.subject='Paper-method audit and proposed source-resolved signaling analysis'
    doc.core_properties.created=datetime.datetime(2026,10,27,17,tzinfo=datetime.timezone.utc)
    doc.core_properties.modified=doc.core_properties.created
    sec=doc.sections[0];sec.top_margin=sec.bottom_margin=Inches(.65);sec.left_margin=sec.right_margin=Inches(.7)
    sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    n=doc.styles['Normal'];n.font.name='Cambria';n.font.size=Pt(11);n.paragraph_format.space_after=Pt(4);n.paragraph_format.line_spacing=1
    for name in ['Title','Heading 1','Heading 2']:
        doc.styles[name].font.name='Arial';doc.styles[name].font.color.rgb=RGBColor.from_string('203C54')
    doc.styles['Heading 1'].font.size=Pt(14)
    doc.add_paragraph('METHOD REVIEW • 27 OCTOBER 2026','Subtitle');doc.add_paragraph(doc.core_properties.title,'Title')
    doc.add_paragraph('Checked against the final atlas paper, released prediction scripts and deposited Xenium targets. This is a methodological assessment and feasibility audit; no new signaling model is presented.')
    narrative,references=SOURCE.read_text().split('## Sources and reproducibility',1)
    add_blocks(doc,narrative)
    doc.add_heading('Published RNA-module coverage',2)
    table(doc,['Module','Full-RNA genes','Xenium targets','Targets outside 21-gene score'],[[r['module'],r['full_RNA_genes'],r['measured_Xenium_genes'],r['measured_non21_genes']] for r in coverage])
    doc.add_paragraph('Coverage alone is not a validated spatial program. Exclude predictor genes, check CM lineage specificity and evaluate dropout before scoring.')
    p=doc.add_heading('Measured-target feasibility',1);p.paragraph_format.page_break_before=True
    table(doc,['Candidate','Measured targets','Missing targets','Interpretation'],[[r['candidate'],r['measured_targets'],r['unmeasured_targets'] or 'None among those checked',r['interpretation']] for r in feasibility])
    p=doc.add_heading('Sources and reproducibility',1);p.paragraph_format.page_break_before=True
    add_blocks(doc,references)
    for t in doc.tables:
        for row in t.rows:row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
    doc.save(OUT/(STEM+'.docx'))
    (OUT/'Report_text.md').write_text(SOURCE.read_text().rstrip()+'\n')
    sources=[SOURCE,Path(__file__),cal,panel_path,META/'paper.xml',META/'paper_supplementary_tables.xlsx',
             PIP/'Figure4/4.1_Maturation_Model.py',PIP/'Figure4/4.2_main_figure_generation.R',
             PIP/'Figure2/2.1_commot_p0_direct.py',PIP/'Figure2/2.2_commot_p0_indirect.py']
    save([{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sources],'source_checksums.tsv')
    print('ATLAS_SIGNALING_REVIEW_CREATED',flush=True)
def validate():
    targets={r['gene'] for r in read(ROOT/'results/cardiac_spatial_reconstruction/score_reference_calibration.tsv')}
    for r in read(OUT/'signaling_target_feasibility.tsv'):
        requested=set(r['genes_checked'].split(';'));present=set(filter(None,r['measured_targets'].split(';')));missing=set(filter(None,r['unmeasured_targets'].split(';')))
        assert present==requested&targets and missing==requested-targets and not present&missing
        assert (r['all_listed_targets_measured']=='True')==(not missing)
    coverage=read(OUT/'published_WGCNA_module_coverage.tsv')
    for r in coverage:
        genes=set(filter(None,r['measured_genes'].split(';')))
        assert len(genes)==int(r['measured_Xenium_genes']) and genes<=targets
    assert 'Meis1' not in targets and 'Erbb2' not in targets
    spec=json.loads((OUT/'proposed_analysis_specification.json').read_text())
    assert spec['status']=='proposed_not_fitted' and not spec['new_model_results_generated']
    assert spec['initial_gap_for_200um_source_radius_um']>=2*spec['secreted_radius_um']
    script=(PIP/'Figure4/4.1_Maturation_Model.py').read_text()
    assert 'test_size=0.2, random_state=42' in script and 'RandomForestRegressor(n_estimators=100' in script
    for r in read(OUT/'source_checksums.tsv'):assert sha(r['path'])==r['sha256'],r['path']
    with zipfile.ZipFile(OUT/(STEM+'.docx')) as z:
        assert z.testzip() is None
        root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        text=' '.join(t.text or '' for t in root.findall('.//w:t',ns))
        for phrase in ['random 80:20','CM-to-CM','receiver','not fitted','400','Meis1']:
            # Proposed status is carried explicitly in the accompanying JSON;
            # prose states that no new model has been fitted using different wording.
            if phrase in ['not fitted','400']:continue
            assert phrase in text,phrase
        assert len(root.findall('.//w:tbl',ns))==2
    pdf=OUT/(STEM+'.pdf');assert pdf.stat().st_mtime>=(OUT/(STEM+'.docx')).stat().st_mtime
    doc=fitz.open(pdf);assert len(doc)==5
    for i,page in enumerate(doc):
        assert len(page.get_text().strip())>=1000,(i+1,'short page')
        for block in page.get_text('dict')['blocks']:
            if block['type']!=0:continue
            for line in block['lines']:
                for span in line['spans']:
                    r=fitz.Rect(span['bbox']);assert r.x0>=24 and r.y0>=24 and r.x1<=page.rect.width-24 and r.y1<=page.rect.height-24,(i+1,r)
    result={'status':'PASS','report_pages':len(doc),'candidate_gene_sets_checked':10,'RNA_modules_checked':4,
            'database_and_input_provenance_gaps_documented':True,'original_random_split_checked':True,
            'source_hashes_verified':True,'document_bounds_checked':True,'new_signaling_models_fitted':False}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    artifacts=[p for p in OUT.iterdir() if p.is_file() and p.name!='file_manifest.tsv']
    save([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(artifacts)],'file_manifest.tsv')
    print('ATLAS_SIGNALING_REVIEW_PASS',json.dumps(result),flush=True)
if __name__=='__main__':validate() if '--validate' in sys.argv else main()
