"""Validate, package and hash the focused maturation/non-CM release."""
from pathlib import Path
import csv, json, re, hashlib, zipfile, argparse, subprocess, sys
import xml.etree.ElementTree as ET
import fitz
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM_maturation21_nonCM_context'
STEM='CM_maturation_and_fibroblast_development'
CODE=['assess_CM_maturation21_neighborhoods.py','validate_CM_maturation21_context.py','plot_CM_maturation21_context.py',
      'CM_maturation21_context_report.md','write_CM_maturation21_report.py','finalize_CM_maturation21_report.py',
      'cardiac_subtype_niches.py','write_cardiac_story.py','write_cardiac_biological_story.py']

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def manifest(paths,name):
    with open(OUT/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
        for p in sorted(set(paths)):w.writerow({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
def docs():
    legend=json.loads((OUT/'figure_legends.json').read_text());assert len(legend)==6
    expected=set();plates=fitz.open();min_dpi=1e9
    for r in legend:
        png=OUT/'figures'/(r['stem']+'.png');pdf=png.with_suffix('.pdf')
        with Image.open(png) as im:
            assert min(im.size)>2000 and min(im.info['dpi'])>=299;im.verify()
        expected.add(sha(png))
        with fitz.open(pdf) as d:assert len(d)==1;plates.insert_pdf(d)
    for page in plates:
        for im in page.get_image_info():
            rect=fitz.Rect(im['bbox']);min_dpi=min(min_dpi,im['width']*72/rect.width,im['height']*72/rect.height)
    assert min_dpi>=299
    plates.save(OUT/'CM_maturation_integrated_figures.pdf',garbage=4,deflate=True)
    with zipfile.ZipFile(OUT/(STEM+'.docx')) as z:
        assert z.testzip() is None
        embedded={hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}
        assert embedded==expected
        root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        text=' '.join(n.text or '' for n in root.findall('.//w:t',ns))
        for claim in ['1,108','514','235','21-gene','not clearly more mature','0.031','0.006','58%','0.281','0.102','2,716','1,339','not assessable','Hortells','TGFβ','Wormhole']:
            assert claim in text,claim
        assert len(root.findall('.//w:tbl',ns))==3
    pdf=OUT/(STEM+'.pdf');assert pdf.stat().st_mtime >= (OUT/(STEM+'.docx')).stat().st_mtime
    doc=fitz.open(pdf);figure_pages=[];short=[]
    for i,page in enumerate(doc):
        content=page.get_text().strip()
        if len(content)<250 or (not page.get_image_info() and len(content)<800 and 'References' not in content):
            short.append({'page':i+1,'text_characters':len(content),'text':content[:150]})
        for im in page.get_image_info():
            r=fitz.Rect(im['bbox']);assert r.x0>=25 and r.y0>=25 and r.x1<=page.rect.width-25 and r.y1<=page.rect.height-25,(i+1,r)
            figure_pages.append(i+1)
    assert not short,short
    assert len(figure_pages)==len(set(figure_pages))==6
    result={'status':'PASS','report_pages':len(doc),'figure_pages':figure_pages,'current_PNG_hashes_embedded':True,
            'minimum_standalone_raster_DPI':min_dpi,'figures_checked_for_bounds':True,'tables':3}
    (OUT/'document_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
def package():
    tables=[p for p in OUT.iterdir() if p.is_file() and p.name.endswith(('.tsv','.tsv.gz')) and p.name not in ['source_checksums.tsv','file_manifest.tsv','source_table_index.tsv']]
    index=[]
    with zipfile.ZipFile(OUT/'CM_maturation_source_tables.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(tables):
            info=zipfile.ZipInfo(p.name,date_time=(2026,10,27,16,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,p.read_bytes());index.append({'archive_path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
    with zipfile.ZipFile(OUT/'CM_maturation_source_tables.zip') as z:assert z.testzip() is None
    with open(OUT/'source_table_index.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['archive_path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(index)
    return len(tables)
def main():
    p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.resume:
        assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
        for file in ['source_checksums.tsv','file_manifest.tsv']:
            for r in rows(OUT/file):assert sha(ROOT/r['path'])==r['sha256'],r['path']
        print('MATURATION21_RELEASE_REUSED_AFTER_HASH_VERIFICATION');return
    analysis=json.loads((OUT/'analysis_validation.json').read_text());ind=json.loads((OUT/'independent_validation.json').read_text())
    assert analysis['status']==ind['status']=='PASS'
    doc=docs();n=package()
    inputs=set(Path(p) for p in json.loads((OUT/'input_paths.json').read_text()) if not Path(p).is_relative_to(OUT))
    inputs|={ROOT/'analysis'/p for p in CODE}
    manifest(inputs,'source_checksums.tsv')
    links=0
    for name in ['README.md','Report_text.md']:
        for link in re.findall(r'\]\(([^)]+)\)',(OUT/name).read_text()):
            if link.startswith(('http://','https://','#')):continue
            assert (OUT/link.split('#')[0]).exists(),link;links+=1
    paths={p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json']}
    manifest(paths,'file_manifest.tsv')
    result={'status':'PASS','report_pages':doc['report_pages'],'integrated_figure_plates':6,'source_tables_in_archive':n,
        'individual_per_cell_score_and_neighborhood_files':8,'relative_links_checked':links,
        'independent_neighbor_checks':ind['independent_neighbor_fraction_checks'],'independent_regression_checks':ind['independently_reconstructed_regression_coefficients'],
        'new_QC_capillary_pairs':ind['QC_capillary_pairs_checked'],'new_primary_neighborhood_tests':235,'matched_primary_tests':20,
        'new_Fibroblast_fate_trajectory':False,'exact_author_Seurat_score_reproduction_claimed':False,'new_biological_replicates':0}
    (OUT/'release_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('MATURATION21_RELEASE_PASS',json.dumps(result),flush=True)
if __name__=='__main__':main()
