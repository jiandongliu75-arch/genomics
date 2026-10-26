"""Validate current embedded figures and the rendered manuscript, not just ZIP syntax."""
from pathlib import Path
import json,zipfile,hashlib,xml.etree.ElementTree as ET
from PIL import Image
import fitz

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_developmental_story'
def digest(data):return hashlib.sha256(data).hexdigest()

def main():
    structure=json.loads((OUT/'manuscript_structure.json').read_text());figures=structure['figures']
    expected={digest((OUT/'figures'/(name+'.png')).read_bytes()) for name in figures}
    word=OUT/'Postnatal_cardiac_development_manuscript.docx'
    with zipfile.ZipFile(word) as z:
        assert z.testzip() is None
        tree=ET.fromstring(z.read('word/document.xml'))
        ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','pic':'http://schemas.openxmlformats.org/drawingml/2006/picture'}
        assert len(tree.findall('.//pic:pic',ns))==6
        images={digest(z.read(n)) for n in z.namelist() if n.startswith('word/media/')}
        assert images==expected,'Word contains stale or missing figure images'
        text=' '.join(e.text or '' for e in tree.findall('.//w:t',ns))
        assert not any(k in text for k in ['TODO','TBC','PLACEHOLDER','Insert figure here'])
        for phrase in ['Abstract','Introduction','Results','Discussion','Methods','Figure legends','33–50%','77–85%','candidate','one spatial section per age']:
            assert phrase.lower() in text.lower(),phrase
    pdfpath=OUT/'Postnatal_cardiac_development_manuscript.pdf';assert pdfpath.stat().st_mtime>=word.stat().st_mtime
    pdf=fitz.open(pdfpath);assert 15<=len(pdf)<=28
    actual_images=[];all_text=[]
    for i,page in enumerate(pdf):
        t=page.get_text();assert len(t.strip())>15,'Blank manuscript page';all_text.append(t)
        for info in page.get_image_info():
            rect=fitz.Rect(info['bbox']);assert rect.x0>=25 and rect.y0>=25 and rect.x1<=page.rect.width-25 and rect.y1<=page.rect.height-25
            assert rect.width>=400 and rect.height>=450
            actual_images.append({'page':i+1,'width_points':rect.width,'height_points':rect.height})
    assert len(actual_images)==6 and len({p['page'] for p in actual_images})==6
    assert '\ufffd' not in ''.join(all_text)
    combined=fitz.open(OUT/'Cardiac_development_integrated_figures.pdf');assert len(combined)==6
    for name in figures:
        p=OUT/'figures'/(name+'.png')
        with Image.open(p) as im:
            assert min(im.size)>=2000 and min(im.info['dpi'])>=299
            im.verify()
        vector=fitz.open(OUT/'figures'/(name+'.pdf'));assert len(vector)==1 and len(vector[0].get_text())>100
        ET.parse(OUT/'figures'/(name+'.svg'))
    check={'status':'PASS','Word_embeds_exact_current_PNG_hashes':True,'Word_embedded_figures':6,'rendered_manuscript_pages':len(pdf),
        'figure_plates_full_page_and_within_page_bounds':actual_images,'standalone_combined_figure_PDF_pages':6,
        'individual_figures_have_300dpi_PNG_vector_PDF_and_editable_SVG':True,'all_figure_PNGs_at_least2000px':True,
        'required_manuscript_sections_and_scope_present':True,'placeholder_or_replacement_characters_found':False,
        'scope':'File content, numerical-claim presence and layout bounds; visual review is additionally required'}
    (OUT/'document_validation.json').write_text(json.dumps(check,indent=2)+'\n')
    print('STORY_DOCUMENT_VALIDATION_PASS',len(pdf),'pages,6 current embedded figure plates',flush=True)

if __name__=='__main__':main()
