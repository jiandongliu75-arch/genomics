"""Run and release the alternative FB identity/program reassessment."""
from pathlib import Path
import os,sys,argparse,subprocess,json,csv,hashlib,zipfile,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/fibroblast_multimethod_identity'
SCRIPTS=['fibroblast_multimethod_identity.py','fibroblast_multimethod_diagnostics.py','plot_fibroblast_multimethod_identity.py',
 'validate_fibroblast_multimethod_identity.py','write_fibroblast_multimethod_report.py','run_fibroblast_multimethod_identity.py',
 'fibroblast_state_common.py','write_cardiac_story.py']
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
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
    names=json.loads((OUT/'report_structure.json').read_text())['figures'];assert len(names)==4
    combined=fitz.open();expected=set();dpi=[]
    for name in names:
        p=OUT/'figures'/(name+'.png');expected.add(sha(p))
        with Image.open(p) as im:assert min(im.size)>1900 and min(im.info['dpi'])>=299;im.verify()
        ET.parse(OUT/'figures'/(name+'.svg'))
        with fitz.open(OUT/'figures'/(name+'.pdf')) as d:assert len(d)==1;combined.insert_pdf(d)
    combined.save(OUT/'Fibroblast_multimethod_integrated_figures.pdf',garbage=4,deflate=True)
    for page in combined:
        for im in page.get_image_info():
            r=fitz.Rect(im['bbox']);dpi.extend([im['width']*72/r.width,im['height']*72/r.height])
    assert min(dpi)>=299
    word=OUT/'Fibroblast_alternative_identity_reassessment.docx'
    with zipfile.ZipFile(word) as z:
        assert z.testzip() is None
        assert {hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}==expected
        tree=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        text=' '.join(t.text or '' for t in tree.findall('.//w:t',ns))
        for s in ['1,339','1,058','280','84%','91%','fate','Nuclear support','Col8a1','Lsamp']:assert s.lower() in text.lower(),s
        assert not any(s in text for s in ['TODO','PLACEHOLDER'])
    pdf=OUT/'Fibroblast_alternative_identity_reassessment.pdf';assert pdf.stat().st_mtime>=word.stat().st_mtime
    d=fitz.open(pdf);pages=[]
    for i,p in enumerate(d):
        assert len(p.get_text().strip())>250,('blank_or_orphan_page',i+1)
        for im in p.get_image_info():
            r=fitz.Rect(im['bbox']);assert r.x0>=25 and r.y0>=25 and r.x1<=p.rect.width-25 and r.y1<=p.rect.height-25
            pages.append(i+1)
    assert len(pages)==len(set(pages))==4
    (OUT/'document_validation.json').write_text(json.dumps({'status':'PASS','rendered_pages':len(d),'integrated_figure_pages':pages,'Word_embeds_current_figure_hashes':True,
      'combined_minimum_raster_DPI':min(dpi),'figure_bounds_and_no_orphan_pages_checked':True,'biological_and_fate_scope_checked':True},indent=2)+'\n')
    print('FB_MULTIMETHOD_DOCUMENT_VALIDATION_PASS',len(d),'pages',flush=True)
def release():
    for name in ['numerical_validation','document_validation']:assert json.loads((OUT/(name+'.json')).read_text())['status']=='PASS'
    sources={ROOT/p for p in json.loads((OUT/'benchmark_input_paths.json').read_text())}|{ROOT/'analysis'/s for s in SCRIPTS}
    manifest(sources,'source_checksums.tsv')
    for name,paths in [('Fibroblast_multimethod_source_data.zip',[p for p in OUT.iterdir() if p.suffix in ['.tsv','.gz','.npz','.json'] and p.name not in ['release_validation.json','file_manifest.tsv']]),
                       ('Fibroblast_multimethod_publication_figures.zip',list((OUT/'figures').iterdir()))]:
        with zipfile.ZipFile(OUT/name,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for p in sorted(paths):z.write(p,p.relative_to(OUT))
        with zipfile.ZipFile(OUT/name) as z:assert z.testzip() is None
    links=0
    for name in ['README.md','Report_text.md']:
        p=OUT/name
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http://','https://','#')):continue
            assert (p.parent/link.split('#')[0]).exists(),link;links+=1
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json']]
    assert all(p.stat().st_size<95*2**20 for p in files);manifest(files,'file_manifest.tsv')
    (OUT/'release_validation.json').write_text(json.dumps({'status':'PASS','relative_links_checked':links,'numerical_and_document_validation':'PASS',
        'legacy_atlas_and_CM_models_and_transport_preserved':True,'new_fate_or_Xenium_accuracy_established':False,'new_biological_replicates':0,
        'prediction_arrays':'variant x view x method x cell; labelled axes; confidence/margin float32'},indent=2)+'\n')
    print('FB_MULTIMETHOD_RELEASE_VALIDATION_PASS',flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('--steps',nargs='+',choices=['analysis','diagnostics','figures','numerical','word','render','documents','release']);p.add_argument('--document-python',default='python');p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.resume:
        assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(OUT/name) as f:
                for r in csv.DictReader(f,delimiter='\t'):assert sha(ROOT/r['path'])==r['sha256'],r['path']
        print('FB multimethod release reused after input/code/artifact hash checks');return
    commands={'analysis':[sys.executable,'analysis/fibroblast_multimethod_identity.py'],
      'diagnostics':[sys.executable,'analysis/fibroblast_multimethod_diagnostics.py'],'figures':[sys.executable,'analysis/plot_fibroblast_multimethod_identity.py'],
      'numerical':[sys.executable,'analysis/validate_fibroblast_multimethod_identity.py'],'word':[a.document_python,'analysis/write_fibroblast_multimethod_report.py'],
      'render':['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/'Fibroblast_alternative_identity_reassessment.docx')]}
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/fb-multimethod-mpl','XDG_CACHE_HOME':'/tmp/fb-multimethod-cache'}
    for step in a.steps or ['analysis','diagnostics','figures','numerical','word','render','documents','release']:
        print('FB_MULTIMETHOD_STEP',step,flush=True)
        if step=='documents':subprocess.run([a.document_python,'-c','import sys;sys.path.insert(0,"analysis");from run_fibroblast_multimethod_identity import documents;documents()'],cwd=ROOT,env=env,check=True)
        elif step=='release':release()
        else:subprocess.run(commands[step],cwd=ROOT,env=env,check=True)
if __name__=='__main__':main()
