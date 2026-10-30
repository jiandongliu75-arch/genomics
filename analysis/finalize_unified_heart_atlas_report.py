"""Build and validate the unified report without refitting biological models."""
from pathlib import Path
import os,csv,json,re,hashlib,zipfile,argparse,subprocess,sys,xml.etree.ElementTree as ET
import fitz
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/heart_atlas_unified_report';STEM='Postnatal_heart_unified_biological_report'
CODE=['heart_atlas_unified_report.md','plot_unified_heart_atlas.py','plot_unified_heart_matrix.py','write_unified_heart_atlas_report.py','finalize_unified_heart_atlas_report.py',
 'write_cardiac_story.py','write_cardiac_biological_story.py','write_clear_heart_atlas_report.py','write_fibroblast_CM_refit_report.py']
TABLES=[
 'cardiac_developmental_story/fibroblast_library_matrix_programs.tsv','cardiac_developmental_story/fibroblast_library_marker_profiles.tsv','cardiac_developmental_story/reviewed_fibroblast_candidate_state_continuities.tsv',
 'fibroblast_multimethod_identity/heldout_naming_marker_profiles.tsv','fibroblast_multimethod_identity/RNA_reference_subcluster_marker_audit.tsv','fibroblast_multimethod_identity/RNA_opposite_library_model_benchmark.tsv','fibroblast_multimethod_identity/published_Gfpt2_Col8a1_program_overlap.tsv',
 'cardiac_spatial_reconstruction/primary_matched_component_results.tsv','cardiac_spatial_reconstruction/nuclear_cycling_specificity_coefficients.tsv',
 'cardiac_subtype_niches/CM_identity_restricted_Wormhole_assessment.tsv','cardiac_subtype_niches/CM_maturation_rank_correspondence_audit.tsv',
 'author_atlas_baseline/descriptive_distances_to_author_EC_states.tsv','whole_heart_Ankrd1_vessels/primary_local_Ankrd1_results.tsv']
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
  for p in sorted(set(Path(p) for p in paths)):w.writerow({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
def checks():
 ref=ROOT/'results/cardiac_fibroblast_CM_refit';assert json.loads((ref/'release_validation.json').read_text())['status']=='PASS'
 data=rows(ref/'CM_fibroblast_refitted_coefficients.tsv');p=[r for r in data if r['primary_analysis']=='True'];hits=[r for r in p if float(r['primary_grid_BH_q'])<.05]
 assert len(p)==925 and len(hits)==1 and hits[0]['target']=='EC__Endocardial-like'
 coverage=rows(ref/'fibroblast_assignment_coverage.tsv');assert all(r['candidate_profiles']=='NA' for r in coverage if r['reference_assessable']=='False')
 vascular=rows(ROOT/'results/cardiac_spatial_reconstruction/primary_matched_component_results.tsv')
 for stage,outcome in [('P7','core_Z__withdrawal'),('P14','core_Z__metabolic'),('P21','core_Z__metabolic')]:assert float(next(r for r in vascular if r['stage']==stage and r['outcome']==outcome)['block_CI_low'])>0
 s=next(r for r in rows(ref/'updated_matched_fibroblast_spatial_context.tsv') if r['stage']=='P14' and r['new_label_layer']=='candidate' and r['measure']=='arterial_within_50um');assert round(100*float(s['target_matched_mean']))==11 and round(100*float(s['reference_matched_mean']))==7
 old=rows(ROOT/'results/fibroblast_multimethod_identity/heldout_naming_marker_profiles.tsv');g=[r for r in old if r['stage']=='P21' and r['view']=='refined' and r['grouping']=='withheld_majority' and r['identity']=='Gfpt2-associated'];assert int(g[0]['profiles'])==280
 info={'status':'PASS','latest_refit_primary_tests':925,'latest_refit_corrected_hits':1,'FB_corrected_hits':0,'updated_P14_arterial_percentages_rounded':[11,7],
 'earlier_capillary_comparisons_retained_separately':True,'reference_unavailability_not_a_biological_zero':True,'new_models_fitted_for_unified_document':False}
 (OUT/'data_claim_checks.json').write_text(json.dumps(info,indent=2)+'\n')
def docs():
 inventory=json.loads((OUT/'figure_inventory.json').read_text());assert len(inventory)==13
 expected=set();sources=set();plate=fitz.open();dpi=[]
 for item in inventory:
  stem=ROOT/item['source_stem'];png=stem.with_suffix('.png');expected.add(sha(png))
  with Image.open(png) as im:assert min(im.size)>1900 and min(im.info['dpi'])>=299;im.verify()
  ET.parse(stem.with_suffix('.svg'))
  with fitz.open(stem.with_suffix('.pdf')) as p:assert len(p)==1;plate.insert_pdf(p)
  sources|={stem.with_suffix('.'+ext) for ext in ['png','pdf','svg']}
 for page in plate:
  for im in page.get_image_info():
   r=fitz.Rect(im['bbox']);dpi.extend([im['width']*72/r.width,im['height']*72/r.height])
 assert min(dpi)>=299;plate.save(OUT/'Postnatal_heart_unified_figure_plates.pdf',garbage=4,deflate=True)
 word=OUT/(STEM+'.docx')
 with zipfile.ZipFile(word) as z:
  assert z.testzip() is None;assert {hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}==expected
  root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};text=' '.join(t.text or '' for t in root.findall('.//w:t',ns))
  for claim in ['925','7,992','10,867','1,339','327','280','Not assessable','11%','7%','Myl2','DecontX','Kcnd2','Wormhole','Capillary-like','Endocardial-like']:
   assert claim in text,claim
  for obsolete in ['have not yet been used to recalculate','have not been refitted','still need to be incorporated','outside that envelope','TODO','PLACEHOLDER']:assert obsolete not in text,obsolete
  assert text.count('What changed after the refit?')==1
  assert len(root.findall('.//w:tbl',ns))>=8
 pdf=OUT/(STEM+'.pdf');assert pdf.stat().st_mtime>=word.stat().st_mtime;d=fitz.open(pdf);figpages=[]
 for i,page in enumerate(d):
  assert len(page.get_text().strip())>250,('blank_or_orphan_page',i+1)
  for im in page.get_image_info():
   r=fitz.Rect(im['bbox']);assert r.x0>=25 and r.y0>=25 and r.x1<=page.rect.width-25 and r.y1<=page.rect.height-25
   figpages.append(i+1)
 assert len(figpages)==len(set(figpages))==13
 out={'status':'PASS','report_pages':len(d),'main_figures':9,'supporting_figures':4,'figure_pages':figpages,'current_figure_hashes_embedded':True,
 'minimum_standalone_raster_DPI':min(dpi),'obsolete_analysis_status_removed':True,'consolidated_evidence_tables_included':True,'bounds_and_page_content_checked':True}
 (OUT/'document_validation.json').write_text(json.dumps(out,indent=2)+'\n');return sources,out

def package():
 ref=ROOT/'results/cardiac_fibroblast_CM_refit';files=[p for p in ref.iterdir() if p.is_file() and p.name.endswith(('.tsv','.tsv.gz')) and p.name not in ['file_manifest.tsv','source_checksums.tsv']]
 files += [ROOT/'results'/p for p in TABLES];index=[]
 with zipfile.ZipFile(OUT/'Postnatal_heart_unified_source_tables.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in sorted(files):
   relative=p.relative_to(ROOT/'results');info=zipfile.ZipInfo(str(relative),date_time=(2026,10,27,12,0,0));info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,p.read_bytes())
   index.append({'archive_path':str(relative),'original_path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
 with open(OUT/'source_table_index.tsv','w') as f:
  w=csv.DictWriter(f,fieldnames=['archive_path','original_path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(index)
 with zipfile.ZipFile(OUT/'Postnatal_heart_unified_source_tables.zip') as z:assert z.testzip() is None
 return set(files),len(index)
def finalize():
 checks();figsources,doc=docs();tables,n=package()
 sources=figsources|tables|{ROOT/'analysis'/p for p in CODE}|{ROOT/'results/cardiac_fibroblast_CM_refit/Fibroblast_CM_spatial_reassessment.docx',ROOT/'results/cardiac_fibroblast_CM_refit/release_validation.json',ROOT/'results/cardiac_biological_synthesis/Heart_Atlas_segmentation_assessment.md'}
 for stage in ['P0','P7','P14','P21']:sources.add(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz')
 manifest(sources,'source_checksums.tsv');(OUT/'release_validation.json').write_text('{"status":"VALIDATING"}\n');links=0
 for name in ['README.md','Report_text.md']:
  for link in re.findall(r'\]\(([^)]+)\)',(OUT/name).read_text()):
   if link.startswith(('http://','https://','#')):continue
   assert (OUT/link.split('#')[0]).exists(),link;links+=1
 paths=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json']];manifest(paths,'file_manifest.tsv')
 result={'status':'PASS','report_pages':doc['report_pages'],'main_figures':9,'supporting_figures':4,'consolidated_source_tables':n,'relative_links_checked':links,
 'latest_FB_refit_incorporated':True,'earlier_results_versions_explicit':True,'obsolete_unfitted_status_removed':True,'new_models_fitted_for_this_document':False,'new_biological_replicates':0}
 (OUT/'release_validation.json').write_text(json.dumps(result,indent=2)+'\n');print('UNIFIED_REPORT_VALIDATION_PASS',json.dumps(result),flush=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--build',action='store_true');p.add_argument('--resume',action='store_true');p.add_argument('--analysis-python',default='/workspace/heart-atlas/venv/bin/python');a=p.parse_args()
 if a.resume:
  assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
  for name in ['source_checksums.tsv','file_manifest.tsv']:
   for row in rows(OUT/name):assert sha(ROOT/row['path'])==row['sha256'],row['path']
  print('UNIFIED_REPORT_REUSED_AFTER_HASH_CHECKS');return
 if a.build:
  env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/unified-heart-mpl','XDG_CACHE_HOME':'/tmp/clear-heart-cache'}
  for command in [[a.analysis_python,'analysis/plot_unified_heart_atlas.py'],[a.analysis_python,'analysis/plot_unified_heart_matrix.py'],[sys.executable,'analysis/write_unified_heart_atlas_report.py'],['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/(STEM+'.docx'))]]:subprocess.run(command,cwd=ROOT,env=env,check=True)
 finalize()
if __name__=='__main__':main()
