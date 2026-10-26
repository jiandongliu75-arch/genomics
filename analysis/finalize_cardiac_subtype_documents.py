"""Merge publication plates and verify the revised Word/PDF layout."""
from pathlib import Path
import json,zipfile,hashlib,xml.etree.ElementTree as ET
from PIL import Image
import fitz

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_subtype_niches'
def sha(data):return hashlib.sha256(data).hexdigest()

def main():
    structure=json.loads((OUT/'manuscript_structure.json').read_text());names=structure['figures'];assert len(names)==6
    combined=fitz.open()
    for name in names:
        with fitz.open(OUT/'figures'/(name+'.pdf')) as p:assert len(p)==1;combined.insert_pdf(p)
        with Image.open(OUT/'figures'/(name+'.png')) as im:assert min(im.size)>=2000 and min(im.info['dpi'])>=299;im.verify()
        ET.parse(OUT/'figures'/(name+'.svg'))
    combined.save(OUT/'Subtype_resolved_integrated_figures.pdf',garbage=4,deflate=True)
    dpi=[]
    for p in combined:
        for info in p.get_image_info():
            rect=fitz.Rect(info['bbox']);dpi.extend([info['width']*72/rect.width,info['height']*72/rect.height])
    assert min(dpi)>=299
    expected={sha((OUT/'figures'/(name+'.png')).read_bytes()) for name in names}
    path=OUT/'Subtype_resolved_cardiac_development_manuscript.docx'
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        media={sha(z.read(n)) for n in z.namelist() if n.startswith('word/media/')};assert media==expected
        text=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        words=' '.join(x.text or '' for x in text.findall('.//w:t',ns))
        for phrase in ['within-stage','Postn','Gfpt2','Col8a1','Ankrd1','source','no two-context','maturation']:
            assert phrase.lower() in words.lower(),phrase
        assert 'results/cardiac_developmental_story and analysis/' not in words
        assert not any(x in words for x in ['TODO','PLACEHOLDER'])
    pdfpath=OUT/'Subtype_resolved_cardiac_development_manuscript.pdf';assert pdfpath.stat().st_mtime>=path.stat().st_mtime
    pdf=fitz.open(pdfpath);assert 15<=len(pdf)<=30;image_pages=[]
    for i,p in enumerate(pdf):
        assert len(p.get_text().strip())>15
        for image in p.get_image_info():
            rect=fitz.Rect(image['bbox']);assert rect.x0>=25 and rect.y0>=25 and rect.x1<=p.rect.width-25 and rect.y1<=p.rect.height-25
            assert rect.width>=400 and rect.height>=450;image_pages.append(i+1)
    assert len(image_pages)==len(set(image_pages))==6
    (OUT/'document_validation.json').write_text(json.dumps({'status':'PASS','Word_embeds_current_figure_hashes':True,'rendered_pages':len(pdf),'main_integrated_figures':6,
        'full_page_figure_pages':image_pages,'combined_figure_minimum_raster_DPI':min(dpi),'editable_SVG_and_vector_PDF_checked':True,
        'subtype_first_text_and_scope_checked':True,'word_to_pdf_layout_bounds_checked':True},indent=2)+'\n')
    print('SUBTYPE_DOCUMENT_VALIDATION_PASS',len(pdf),'pages,6 integrated subtype-first plates',flush=True)

if __name__=='__main__':main()
