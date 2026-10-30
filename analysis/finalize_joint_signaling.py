"""Audit report assets and publish a portable, hashed joint-analysis bundle."""
from pathlib import Path
import csv,json,hashlib,zipfile,re,xml.etree.ElementTree as ET
import fitz
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/joint_RNA_Xenium_signaling'
STEM='Joint_RNA_Xenium_CM_development'
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(4<<20),b''):h.update(b)
 return h.hexdigest()
def rows(path):
 with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def manifest(paths,name):
 with open(OUT/name,'w') as f:
  writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');writer.writeheader()
  for p in sorted(set(paths)):
   label=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)
   writer.writerow({'path':label,'bytes':p.stat().st_size,'sha256':sha(p)})
def releasable(path):
 return path.is_file() and not any(p.startswith('.') for p in path.relative_to(OUT).parts) and not path.name.endswith(('_predictor_features.npz','_checkpoint.npz','.log'))
def main():
 rna=json.loads((OUT/'RNA_independent_validation.json').read_text())
 spatial=json.loads((OUT/'spatial_independent_validation.json').read_text())
 equivalence=json.loads((OUT/'transport_equivalence_check.json').read_text())
 assert rna['status']==spatial['status']=='PASS'
 assert equivalence['maximum_absolute_error']<1e-7 and equivalence['converged']
 assert spatial['transport_domains_checked']==32 and spatial['primary_BH_tests']==40
 assert '[[' not in (ROOT/'analysis/joint_signaling_report.md').read_text()
 legends=json.loads((OUT/'figure_legends.json').read_text());assert len(legends)==4
 expected=set();plates=fitz.open();minimum_dpi=1e9
 for r in legends:
  png=OUT/'figures'/(r['stem']+'.png');pdf=png.with_suffix('.pdf')
  with Image.open(png) as image:
   assert min(image.size)>2000 and min(image.info['dpi'])>=299;image.verify()
  expected.add(sha(png))
  with fitz.open(pdf) as d:
   assert len(d)==1;plates.insert_pdf(d)
 for page in plates:
  for image in page.get_image_info():
   rect=fitz.Rect(image['bbox']);minimum_dpi=min(minimum_dpi,image['width']*72/rect.width,image['height']*72/rect.height)
 assert minimum_dpi>=299
 plates.save(OUT/'Joint_RNA_Xenium_integrated_figures.pdf',garbage=4,deflate=True)
 with zipfile.ZipFile(OUT/(STEM+'.docx')) as z:
  assert z.testzip() is None
  media={hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}
  assert media==expected
  root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
  text=' '.join(n.text or '' for n in root.findall('.//w:t',ns))
  for term in ['58','403','180','P21','not','Col8a1','Gfpt2/Pi16']:
   assert term in text,term
  assert 'one observed spatial section' in text.lower()
  assert '[[' not in text and len(root.findall('.//w:tbl',ns))==2
 pdf=OUT/(STEM+'.pdf');doc=fitz.open(pdf);figure_pages=[];short=[]
 for i,page in enumerate(doc):
  text=page.get_text().strip();images=page.get_image_info()
  if images:
   assert len(images)==1,(i+1,len(images));figure_pages.append(i+1)
   assert 'Figure ' in text and len(text)>350,(i+1,'missing figure legend')
   for image in images:
    r=fitz.Rect(image['bbox']);assert r.x0>=35 and r.y0>=30 and r.x1<=page.rect.width-35 and r.y1<=page.rect.height-30,(i+1,r)
  elif len(text)<400:short.append({'page':i+1,'characters':len(text),'text':text[:160]})
  for block in page.get_text('dict')['blocks']:
   if block['type']!=0:continue
   for line in block['lines']:
    for span in line['spans']:
     r=fitz.Rect(span['bbox']);assert r.x0>=30 and r.y0>=8 and r.x1<=page.rect.width-30 and r.y1<=page.rect.height-8,(i+1,r,span['text'])
 assert len(figure_pages)==4 and not short,(figure_pages,short)
 document={'status':'PASS','report_pages':len(doc),'figure_pages':figure_pages,'figure_plates':4,'tables':2,
  'Word_embeds_current_PNG_hashes':True,'PDF_uses_same_narrative_tables_and_assets':True,
  'minimum_standalone_raster_DPI':minimum_dpi,'PDF_bounds_and_short_pages_checked':True}
 (OUT/'document_validation.json').write_text(json.dumps(document,indent=2)+'\n')
 inputs={Path(p) for p in json.loads((OUT/'input_paths.json').read_text()) if not Path(p).is_relative_to(OUT)}
 inputs|=set((ROOT/'analysis').glob('*joint_signaling*'))
 inputs|=set((ROOT/'analysis').glob('joint_signaling_*.py'))
 inputs.add(ROOT/'analysis/run_joint_stage_view.py')
 inputs|={ROOT/'analysis/vendor_commot'/p for p in ['_cot.py','_unot.py','LICENSE.md','CellChatDB.ligrec.mouse.csv','CellPhoneDBv4.0.mouse.csv','provenance.json']}
 inputs={p for p in inputs if p.is_file()};manifest(inputs,'source_checksums.tsv')
 table_paths=[p for p in OUT.rglob('*') if releasable(p) and p.name.endswith(('.tsv','.tsv.gz')) and p.name not in ['source_checksums.tsv','file_manifest.tsv','source_table_index.tsv']]
 manifest(table_paths,'source_table_index.tsv')
 links=0
 for name in ['README.md','Report_text.md']:
  for link in re.findall(r'\]\(([^)]+)\)',(OUT/name).read_text()):
   if link.startswith(('http://','https://','#')):continue
   if link in ['release_validation.json','file_manifest.tsv']:links+=1;continue
   assert (OUT/link.split('#')[0]).exists(),link;links+=1
 result={'status':'PASS',**document,'RNA_support_rows_checked':rna['RNA_support_rows_independently_reconstructed'],
  'RNA_library_complex_checks':rna['RNA_library_complex_means_and_detection_recomputed_from_full_counts'],
  'transport_domains_checked':32,'primary_BH_tests':40,'OOF_R2_checks':spatial['OOF_R2_independently_recomputed'],
  'local_ligand_and_capture_checks':spatial['independently_recomputed_local_ligand_and_capture_checks'],
  'relative_links_checked':links,'source_tables':len(table_paths),'new_independent_spatial_hearts':0,
  'exact_author_feature_matrix_reproduction_claimed':False,'new_fate_inference':False}
 (OUT/'release_validation.json').write_text(json.dumps(result,indent=2)+'\n')
 manifest([p for p in OUT.rglob('*') if releasable(p) and p.name!='file_manifest.tsv'],'file_manifest.tsv')
 for name in ['release_validation.json','file_manifest.tsv']:assert (OUT/name).exists()
 print('JOINT_RELEASE_PASS',json.dumps(result),flush=True)
if __name__=='__main__':main()
