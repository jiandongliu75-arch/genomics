"""Validate the consolidated Word report against existing source results."""
from pathlib import Path
import csv,json,hashlib,re,zipfile,argparse,os,subprocess,sys
import xml.etree.ElementTree as ET
import fitz
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/heart_atlas_clear_language_revision'
STEM='Postnatal_heart_integrated_clear_language_report'
SCRIPTS=['heart_atlas_clear_language_revision.md','plot_clear_heart_atlas_matrix.py','write_clear_heart_atlas_report.py','finalize_clear_heart_atlas_report.py',
         'write_cardiac_biological_story.py','write_cardiac_story.py']
DATA=['results/fibroblast_multimethod_identity/final_candidate_coverage.tsv',
 'results/fibroblast_multimethod_identity/published_Gfpt2_Col8a1_program_overlap.tsv',
 'results/fibroblast_multimethod_identity/heldout_naming_marker_profiles.tsv',
 'results/fibroblast_multimethod_identity/RNA_reference_subcluster_marker_audit.tsv',
 'results/fibroblast_multimethod_identity/numerical_validation.json',
 'results/cardiac_subtype_niches/source_subtype_coverage.tsv',
 'results/cardiac_developmental_story/fibroblast_library_matrix_programs.tsv',
 'results/cardiac_developmental_story/fibroblast_library_marker_profiles.tsv',
 'results/cardiac_developmental_story/reviewed_fibroblast_candidate_state_continuities.tsv',
 'results/fibroblast_names_and_neighborhoods/matched_family_context_sensitivity.tsv',
 'results/cardiac_biological_synthesis/nuclear_optional_FB_maturation_coefficients.tsv',
 'results/cardiac_biological_synthesis/Heart_Atlas_segmentation_assessment.md',
 'results/cardiac_biological_synthesis/literature_review_provenance.json',
 'results/cardiac_biological_synthesis/numerical_validation.json',
 'results/cardiac_spatial_reconstruction/primary_matched_component_results.tsv',
 'results/cardiac_subtype_niches/CM_maturation_rank_correspondence_audit.tsv',
 'results/cardiac_biological_synthesis/Postnatal_heart_biological_synthesis.docx',
 'results/fibroblast_multimethod_identity/Fibroblast_alternative_identity_reassessment.docx']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(path):
    with open(ROOT/path) as f:return list(csv.DictReader(f,delimiter='\t'))
def manifest(paths,name):
    with open(OUT/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
        for p in sorted(paths):p=Path(p);w.writerow({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})

def facts():
    c=[r for r in rows(DATA[0]) if r['stage']=='P21'];eligible=int(c[0]['eligible_FB_profiles']);total=sum(int(r['candidate_profiles']) for r in c)
    gf=next(r for r in c if r['family']=='Gfpt2-associated')
    p=rows(DATA[1]);allp={r['program_combination']:int(r['profiles']) for r in p if r['stage']=='P21' and r['scope']=='all_eligible'}
    pub={r['program_combination']:int(r['profiles']) for r in p if r['stage']=='P21' and r['scope']=='public_Gfpt2'}
    m=[r for r in rows(DATA[2]) if r['stage']=='P21' and r['view']=='refined' and r['grouping']=='withheld_majority' and r['identity']=='Gfpt2-associated']
    reporter={r['gene']:float(r['detection_fraction']) for r in m};n=int(m[0]['profiles'])
    r=[r for r in rows(DATA[3]) if r['stage']=='P21' and r['view']=='corrected' and r['RNA_subcluster']=='P21_F2' and r['gene']=='Gfpt2']
    detection=sum(int(x['nuclei'])*float(x['detection_fraction']) for x in r)/sum(int(x['nuclei']) for x in r)
    old=[r for r in rows(DATA[5]) if r['stage']=='P21' and r['subtype'].startswith('FB__')];oldtotal=sum(int(r['primary_source_profiles']) for r in old)
    assert (eligible,total,int(gf['candidate_profiles']),int(gf['crosschecked_profiles']))==(10867,7992,85,75)
    assert (allp['Gfpt2_only']+allp['both'],allp['both'],pub['Gfpt2_only']+pub['both'])==(1339,327,1058)
    assert n==280 and round(100*reporter['Gfpt2'])==84 and round(100*reporter['Pi16'])==91 and abs(detection-2/3)<.02
    assert oldtotal==679 and round(100*total/eligible)==74 and round(100*oldtotal/eligible)==6
    spatial=[r for r in rows(DATA[9]) if r['stage']=='P14' and r['radius_um']=='30' and r['measure']=='arterial_within_50um'];assert len(spatial)==1
    assert round(100*float(spatial[0]['target_matched_mean']))==19 and round(100*float(spatial[0]['reference_matched_mean']))==6
    matrix=[r for r in rows(DATA[6]) if r['level']=='fibroblast_family' and r['stage'] in ['P14','P21'] and r['program']=='Basement matrix' and r['group'] in ['Col8a1-associated','Gfpt2-associated']]
    values={(r['stage'],r['library'],r['view'],r['group']):float(r['mean_log1p_CPM']) for r in matrix};checks=0
    for stage,library,view,family in values:
        if family!='Col8a1-associated':continue
        assert values[(stage,library,view,family)]>values[(stage,library,view,'Gfpt2-associated')];checks+=1
    assert checks>=4
    result={'status':'PASS','P21_eligible_FB':eligible,'P21_candidate_FB':total,'P21_legacy_strict_FB':oldtotal,
       'P21_Gfpt2_program':1339,'P21_Gfpt2_and_Col8a1_programs':327,'P21_public_Gfpt2_program':1058,
       'P21_withheld_majority_Gfpt2':n,'withheld_Gfpt2_detection':reporter['Gfpt2'],'withheld_Pi16_detection':reporter['Pi16'],
       'RNA_Col8a1_Hmcn2_Gfpt2_detection':detection,'P14_near_artery_percentages_rounded':[19,6],
       'basement_matrix_direction_checks':checks,'CM_models_refitted_with_latest_classifier_labels':False,'new_biological_inference_fitted':False}
    (OUT/'data_claim_checks.json').write_text(json.dumps(result,indent=2)+'\n')

def documents():
    inventory=json.loads((OUT/'figure_inventory.json').read_text());assert len(inventory)==9
    expected=set();combined=fitz.open();sources=set()
    for item in inventory:
        stem=ROOT/item['source_stem'];png=stem.with_suffix('.png');expected.add(sha(png))
        with Image.open(png) as im:assert min(im.size)>1900 and min(im.info['dpi'])>=299;im.verify()
        ET.parse(stem.with_suffix('.svg'))
        with fitz.open(stem.with_suffix('.pdf')) as p:assert len(p)==1;combined.insert_pdf(p)
        sources|={stem.with_suffix('.'+ext) for ext in ['png','pdf','svg']}
    combined.save(OUT/'Heart_atlas_integrated_figure_plates.pdf',garbage=4,deflate=True)
    dpi=[]
    for page in combined:
        for im in page.get_image_info():
            r=fitz.Rect(im['bbox']);dpi.extend([im['width']*72/r.width,im['height']*72/r.height])
    assert min(dpi)>=299
    word=OUT/(STEM+'.docx')
    with zipfile.ZipFile(word) as z:
        assert z.testzip() is None
        assert {hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}==expected
        tree=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        text=' '.join(t.text or '' for t in tree.findall('.//w:t',ns))
        for s in ['7,992','10,867','1,339','327','1,058','2,716','280','84%','91%','75','19','50 µm','not yet been used to recalculate']:
            assert s in text,s
        assert 'outside that envelope' not in text and not any(s in text for s in ['TODO','PLACEHOLDER'])
        assert text.count('Figure S1.')==1 and text.count('Figure S2.')==1
    pdf=OUT/(STEM+'.pdf');assert pdf.stat().st_mtime>=word.stat().st_mtime
    d=fitz.open(pdf);pages=[]
    for i,page in enumerate(d):
        assert len(page.get_text().strip())>250,('blank_or_orphan_page',i+1)
        for im in page.get_image_info():
            r=fitz.Rect(im['bbox']);assert r.x0>=25 and r.y0>=25 and r.x1<=page.rect.width-25 and r.y1<=page.rect.height-25
            pages.append(i+1)
    assert len(pages)==len(set(pages))==9
    (OUT/'document_validation.json').write_text(json.dumps({'status':'PASS','rendered_pages':len(d),'main_figures':7,'supporting_figures':2,'figure_pages':pages,
       'current_png_hashes_embedded':True,'combined_minimum_raster_DPI':min(dpi),'data_claims_checked':True,
       'figure_bounds_and_no_orphan_pages_checked':True,'analysis_status_and_requested_language_checked':True},indent=2)+'\n')
    return sources

def release():
    facts();sources=documents();sources|={ROOT/p for p in DATA}|{ROOT/'analysis'/s for s in SCRIPTS}
    manifest(sources,'source_checksums.tsv')
    links=0
    for name in ['README.md','Report_text.md']:
        p=OUT/name
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http://','https://','#')):continue
            assert (p.parent/link.split('#')[0]).exists(),link;links+=1
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json']]
    manifest(paths,'file_manifest.tsv')
    (OUT/'release_validation.json').write_text(json.dumps({'status':'PASS','data_claim_and_document_checks':'PASS','relative_links_checked':links,
       'uploaded_reports_combined_in_clear_biological_language':True,'earlier_analysis_artifacts_preserved':True,
       'CM_associations_and_transport_refitted':False,'new_biological_replicates':0},indent=2)+'\n')
    print('CLEAR_HEART_REPORT_VALIDATION_PASS',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');p.add_argument('--build',action='store_true');p.add_argument('--analysis-python',default='/workspace/heart-atlas/venv/bin/python');a=p.parse_args()
    if a.resume:
        assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(OUT/name) as f:
                for r in csv.DictReader(f,delimiter='\t'):assert sha(ROOT/r['path'])==r['sha256'],r['path']
        print('Clear-language report reused after source/code/artifact hash checks');return
    if a.build:
        env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/clear-heart-mpl','XDG_CACHE_HOME':'/tmp/clear-heart-cache'}
        for command in [[a.analysis_python,'analysis/plot_clear_heart_atlas_matrix.py'],[sys.executable,'analysis/write_clear_heart_atlas_report.py'],
            ['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/(STEM+'.docx'))]]:
            subprocess.run(command,cwd=ROOT,env=env,check=True)
    release()
if __name__=='__main__':main()
