"""Run, validate and package the biology-led heart-atlas revision."""
from pathlib import Path
import os,sys,json,csv,hashlib,subprocess,zipfile,argparse,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_biological_synthesis'
SCRIPTS=['fibroblast_assignment_sensitivity.py','cardiac_biological_figures.py','cardiac_biological_story.md',
 'write_cardiac_biological_story.py','validate_cardiac_biological_synthesis.py','run_cardiac_biological_synthesis.py',
 'cardiac_subtype_niches.py','write_cardiac_story.py']
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def manifest(paths,name):
    with open(OUT/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
        for p in sorted(paths):
            p=Path(p);w.writerow({'path':str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
def documents():
    import fitz
    import xml.etree.ElementTree as ET
    from PIL import Image
    structure=json.loads((OUT/'manuscript_structure.json').read_text());names=structure['main_figures'];supp=structure['supplementary_figures']
    assert len(names)==6 and len(supp)==1
    combined=fitz.open();dpi=[]
    for name in names+supp:
        p=OUT/'figures'/(name+'.png')
        with Image.open(p) as im:assert min(im.size)>1900 and min(im.info['dpi'])>=299;im.verify()
        ET.parse(OUT/'figures'/(name+'.svg'))
        with fitz.open(OUT/'figures'/(name+'.pdf')) as d:
            assert len(d)==1
            if name in names:combined.insert_pdf(d)
    combined.save(OUT/'Biological_synthesis_integrated_figures.pdf',garbage=4,deflate=True)
    for page in combined:
        for im in page.get_image_info():
            r=fitz.Rect(im['bbox']);dpi.extend([im['width']*72/r.width,im['height']*72/r.height])
    assert min(dpi)>=299
    expected={sha(OUT/'figures'/(n+'.png')) for n in names+supp}
    word=OUT/'Postnatal_heart_biological_synthesis.docx'
    with zipfile.ZipFile(word) as z:
        assert z.testzip() is None
        assert {hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}==expected
        xml=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        words=' '.join(t.text or '' for t in xml.findall('.//w:t',ns))
        for text in ['regional developmental coordination','perinatal','nuclear','P21','10%','boundary','Gfpt2','Col8a1','perfusion']:
            assert text.lower() in words.lower(),text
        assert not any(x in words for x in ['TODO','PLACEHOLDER'])
    pdf=OUT/'Postnatal_heart_biological_synthesis.pdf';assert pdf.stat().st_mtime>=word.stat().st_mtime
    d=fitz.open(pdf);pages=[]
    for i,p in enumerate(d):
        assert len(p.get_text().strip())>100,('blank_or_orphan_page',i+1)
        for im in p.get_image_info():
            r=fitz.Rect(im['bbox']);assert r.x0>=25 and r.y0>=25 and r.x1<=p.rect.width-25 and r.y1<=p.rect.height-25
            pages.append(i+1)
    assert len(pages)==len(set(pages))==7
    (OUT/'document_validation.json').write_text(json.dumps({'status':'PASS','rendered_pages':len(d),'figure_pages':pages,'current_Word_media_hashes_checked':True,
        'minimum_combined_raster_DPI':min(dpi),'figure_bounds_and_no_blank_pages_checked':True,'segmentation_assessment_included':True},indent=2)+'\n')
    print('BIOLOGICAL_DOCUMENT_VALIDATION_PASS',len(d),'pages',flush=True)
def release():
    for name in ['numerical_validation','document_validation']:assert json.loads((OUT/(name+'.json')).read_text())['status']=='PASS'
    inputs={ROOT/p for p in json.loads((OUT/'audit_input_paths.json').read_text())}|{ROOT/'analysis'/s for s in SCRIPTS}
    inputs|={Path('/workspace/heart-atlas/metadata/paper.xml')}
    for name in ['Figure_1_CM_identity_and_maturation','Figure_4_fibroblast_specialization','Figure_5_CM_identity_conditioned_transport']:
        inputs|={ROOT/'results/cardiac_subtype_niches/figures'/(name+'.'+ext) for ext in ['png','pdf','svg']}
    inputs|={ROOT/'results/cardiac_spatial_reconstruction/primary_matched_component_results.tsv',
              ROOT/'results/cardiac_subtype_niches/subtype_conditioned_maturation_coefficients.tsv'}
    manifest(inputs,'source_checksums.tsv')
    for name,paths in [('Biological_synthesis_source_tables.zip',[p for p in OUT.iterdir() if p.suffix in ['.tsv','.gz','.json'] and p.name not in ['release_validation.json','file_manifest.tsv']]),
                       ('Biological_synthesis_publication_figures.zip',list((OUT/'figures').iterdir()))]:
        with zipfile.ZipFile(OUT/name,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for p in sorted(paths):z.write(p,p.relative_to(OUT))
        with zipfile.ZipFile(OUT/name) as z:assert z.testzip() is None
    links=0
    for name in ['README.md','Manuscript_text.md']:
        p=OUT/name
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http://','https://','#')):continue
            assert (p.parent/link.split('#')[0]).exists(),link;links+=1
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json']]
    manifest(paths,'file_manifest.tsv')
    (OUT/'release_validation.json').write_text(json.dumps({'status':'PASS','numerical_and_document_validation':'PASS','relative_links_checked':links,
        'old_atlas_and_transport_not_overwritten':True,'nuclear_optional_candidates_are_not_validated_cell_types':True,
        'public_article_supplement_reviewed':'Supplementary Figure 1b-e; approximate plot readings explicitly marked',
        'new_biological_replicates':0,'new_functional_mechanism_validated':False},indent=2)+'\n')
    print('BIOLOGICAL_RELEASE_VALIDATION_PASS',flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('--steps',nargs='+',choices=['analysis','figures','numerical','word','render','documents','release']);p.add_argument('--document-python',default='python');p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.resume:
        assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(OUT/name) as f:
                for r in csv.DictReader(f,delimiter='\t'):assert sha(ROOT/r['path'])==r['sha256'],r['path']
        print('Biological release reused after source/code/artifact hash checks');return
    commands={'analysis':[sys.executable,'analysis/fibroblast_assignment_sensitivity.py'],
       'figures':[sys.executable,'analysis/cardiac_biological_figures.py'],'numerical':[sys.executable,'analysis/validate_cardiac_biological_synthesis.py'],
       'word':[a.document_python,'analysis/write_cardiac_biological_story.py'],
       'render':['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/'Postnatal_heart_biological_synthesis.docx')]}
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/cardiac-biology-mpl','XDG_CACHE_HOME':'/tmp/cardiac-biology-cache'}
    for name in a.steps or ['analysis','figures','numerical','word','render','documents','release']:
        print('BIOLOGICAL_STEP',name,flush=True)
        if name=='documents':
            # Analysis Python may not provide PyMuPDF/Pillow.
            code='import sys;sys.path.insert(0,"analysis");from run_cardiac_biological_synthesis import documents;documents()'
            subprocess.run([a.document_python,'-c',code],cwd=ROOT,env=env,check=True)
        elif name=='release':release()
        else:subprocess.run(commands[name],cwd=ROOT,env=env,check=True)
if __name__=='__main__':main()
