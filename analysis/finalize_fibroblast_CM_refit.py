"""Check numerical claims, rendered documents and source/release fingerprints."""
from pathlib import Path
import argparse,csv,json,hashlib,zipfile,re,io,subprocess,sys,xml.etree.ElementTree as ET
import fitz
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_fibroblast_CM_refit'
STEM='Fibroblast_CM_spatial_reassessment'
CODE=['refit_fibroblast_CM_relationships.py','fibroblast_CM_refit_checks.py','fibroblast_CM_refit_extensions.py','fibroblast_CM_gene_exclusion_sensitivity.py',
 'plot_fibroblast_CM_refit.py','write_fibroblast_CM_refit_report.py','fibroblast_CM_refit_report.md','finalize_fibroblast_CM_refit.py','run_fibroblast_CM_refit.py',
 'cardiac_subtype_niches.py','fibroblast_multimethod_identity.py','fibroblast_multimethod_diagnostics.py','fibroblast_state_common.py','name_fibroblasts_and_describe_neighbors.py','write_cardiac_story.py','write_cardiac_biological_story.py']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4<<20),b''):h.update(b)
 return h.hexdigest()
def rows(name):
 with open(OUT/name) as f:return list(csv.DictReader(f,delimiter='\t'))
def manifest(paths,name):
 with open(OUT/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
  for p in sorted(set(Path(p) for p in paths)):
   try:label=str(p.relative_to(ROOT))
   except ValueError:label=str(p)
   w.writerow({'path':label,'bytes':p.stat().st_size,'sha256':sha(p)})
def data_checks():
 assert json.loads((OUT/'numerical_validation.json').read_text())['status']=='PASS'
 data=rows('CM_fibroblast_refitted_coefficients.tsv');primary=[r for r in data if r['primary_analysis']=='True'];hit=[r for r in primary if float(r['primary_grid_BH_q'])<.05]
 assert len(primary)==925 and len(hit)==1 and hit[0]['target']=='EC__Endocardial-like' and hit[0]['CM_state']=='Ankrd1.CM' and hit[0]['stage']=='P14' and hit[0]['outcome']=='core_Z__sarcomere'
 cov=rows('fibroblast_assignment_coverage.tsv');assert all(r['table_display']=='Not assessable' and r['candidate_profiles']=='NA' for r in cov if r['reference_assessable']=='False')
 candidate=[r for r in cov if r['stage']=='P21' and r['label_layer']=='candidate' and r['reference_assessable']=='True'];assert sum(int(float(r['candidate_profiles'])) for r in candidate)==7992
 p=rows('fibroblast_program_coverage.tsv')
 for stage,program,count in [('P21','Postn.Fib',192),('P21','Gfpt2.Fib',1339)]:assert int(next(r for r in p if r['stage']==stage and r['view']=='refined' and r['program']==program)['enriched_profiles'])==count
 spatial=rows('updated_matched_fibroblast_spatial_context.tsv');s=next(r for r in spatial if r['stage']=='P14' and r['new_label_layer']=='candidate' and r['measure']=='arterial_within_50um')
 assert round(100*float(s['target_matched_mean']))==11 and round(100*float(s['reference_matched_mean']))==7
 assert max(abs(float(r['after_standardized_mean_difference'])) for r in rows('updated_spatial_matching_balance.tsv'))<.03
 t=rows('updated_frozen_transport_neighborhood_descriptors.tsv');assert len(t)==96 and sum(int(r['eligible_FB_sources'])==0 for r in t)==15
 u=[r for r in rows('updated_transport_program_comparisons.tsv') if r['stage_a']=='P7' and r['stage_b']=='P14' and r['method']=='exact']
 assert all(float(r['mean_target_minus_source'])<0 for r in u if r['program']=='Postn.Fib') and all(float(r['mean_target_minus_source'])>0 for r in u if r['program']=='Col8a1.Fib')
 no=rows('CM_gene_excluded_fibroblast_program_sensitivity.tsv');el=[r for r in no if r['stage']=='P0' and r['CM_state']=='v.CM' and r['target'].endswith('__Elastic_matrix') and r['outcome']=='core_Z__sarcomere'];assert len(el)==2 and all(float(r['CI_high'])<0 for r in el)
 extra={'status':'PASS','primary_grid_tests':len(primary),'primary_corrected_hits':len(hit),'fibroblast_corrected_hits':0,
   'reference_unavailability_not_presented_as_zero':True,'P21_candidate_FB':7992,'P21_Postn_pattern':192,'P21_Gfpt2_pattern':1339,
   'updated_P14_matched_arterial_percentages_rounded':[11,7],'matching_covariate_balance_max_abs_SD':max(abs(float(r['after_standardized_mean_difference'])) for r in rows('updated_spatial_matching_balance.tsv')),
   'transport_neighborhoods':96,'transport_without_eligible_FB':15,'P7_P14_four_context_pattern_directions_checked':True,'P0_elastic_candidate_CM_excluded_sensitivity_checked':True}
 (OUT/'report_claim_checks.json').write_text(json.dumps(extra,indent=2)+'\n')
 return extra

def documents():
 legends=json.loads((OUT/'figure_legends.json').read_text());assert len(legends)==7
 expected=set();plates=fitz.open();dpis=[]
 for item in legends:
  stem=OUT/'figures'/item['stem'];png=stem.with_suffix('.png');expected.add(sha(png))
  with Image.open(png) as im:assert min(im.size)>1900 and min(im.info['dpi'])>=299;im.verify()
  ET.parse(stem.with_suffix('.svg'))
  with fitz.open(stem.with_suffix('.pdf')) as p:
   assert len(p)==1;plates.insert_pdf(p)
 for page in plates:
  for im in page.get_image_info():
   r=fitz.Rect(im['bbox']);dpis.extend([im['width']*72/r.width,im['height']*72/r.height])
 assert min(dpis)>=299
 plates.save(OUT/'Fibroblast_CM_integrated_figures.pdf',garbage=4,deflate=True)
 word=OUT/(STEM+'.docx')
 with zipfile.ZipFile(word) as z:
  assert z.testzip() is None
  assert {hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}==expected
  tree=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};text=' '.join(t.text or '' for t in tree.findall('.//w:t',ns))
  for claim in ['7,992','10,867','925','1,339','192','Not assessable','11%','7%','Myl2','DecontX','0.05–0.12','CM genes','one available spatial section']:
   assert claim.lower() in text.lower(),claim
  assert 'outside that envelope' not in text and not any(k in text for k in ['TODO','PLACEHOLDER'])
 pdf=OUT/(STEM+'.pdf');assert pdf.stat().st_mtime>=word.stat().st_mtime;d=fitz.open(pdf);figpages=[]
 for i,page in enumerate(d):
  assert len(page.get_text().strip())>250,('blank_or_orphan',i+1)
  for im in page.get_image_info():
   r=fitz.Rect(im['bbox']);assert r.x0>=25 and r.y0>=25 and r.x1<=page.rect.width-25 and r.y1<=page.rect.height-25
   figpages.append(i+1)
 assert len(figpages)==len(set(figpages))==7
 info={'status':'PASS','report_pages':len(d),'main_figures':6,'supporting_figures':1,'figure_pages':figpages,'current_figure_hashes_embedded':True,
  'minimum_standalone_figure_raster_DPI':min(dpis),'figure_bounds_and_page_content_checked':True,'reported_numbers_and_analysis_status_checked':True}
 (OUT/'document_validation.json').write_text(json.dumps(info,indent=2)+'\n');return info

def package():
 paths=sorted([p for p in OUT.iterdir() if p.is_file() and (p.name.endswith(('.tsv','.tsv.gz')) or p.name in ['README.md','numerical_validation.json','report_claim_checks.json','report_structure.json','document_validation.json']) and p.name not in ['file_manifest.tsv','source_checksums.tsv']])
 with zipfile.ZipFile(OUT/'Fibroblast_CM_source_tables.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in paths:
   info=zipfile.ZipInfo(p.name,date_time=(2026,10,26,16,0,0));info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,p.read_bytes())
 with zipfile.ZipFile(OUT/'Fibroblast_CM_source_tables.zip') as z:assert z.testzip() is None
 return len(paths)
def main():
 p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');a=p.parse_args()
 if a.resume:
  assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
  for file in ['source_checksums.tsv','file_manifest.tsv']:
   with open(OUT/file) as f:
    for r in csv.DictReader(f,delimiter='\t'):
     path=Path(r['path']);path=path if path.is_absolute() else ROOT/path;assert sha(path)==r['sha256'],str(path)
  print('REFIT_RELEASE_REUSED_AFTER_HASH_CHECKS');return
 claims=data_checks();doc=documents();n=package()
 inputs={Path(p) for p in json.loads((OUT/'input_paths.json').read_text())}|{ROOT/'analysis'/p for p in CODE}
 for stage in ['P0','P7','P14','P21']:
  inputs.add(ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
  inputs.add(ROOT/f'results/cardiac_niche_dynamics/{stage}/nonCM_state_points.npz')
 inputs|={ROOT/'results/fibroblast_reassessment/fibroblast_QC_by_library.tsv',ROOT/'results/fibroblast_names_and_neighborhoods/matched_family_context_sensitivity.tsv',
  ROOT/'results/fibroblast_reassessment/paper_fibroblast_markers_and_coverage.tsv',ROOT/'results/fibroblast_reassessment/program_gene_coverage.tsv'}
 manifest(inputs,'source_checksums.tsv')
 (OUT/'release_validation.json').write_text(json.dumps({'status':'VALIDATING'})+'\n')
 linkcount=0
 for file in ['README.md','Report_text.md']:
  for link in re.findall(r'\]\(([^)]+)\)',(OUT/file).read_text()):
   if link.startswith(('http://','https://','#')):continue
   assert (OUT/link.split('#')[0]).exists(),link;linkcount+=1
 paths=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json']];manifest(paths,'file_manifest.tsv')
 result={'status':'PASS','report_pages':doc['report_pages'],'numerical_validation':'PASS','report_claim_checks':'PASS','document_validation':'PASS',
   'source_table_archive_files':n,'relative_links_checked':linkcount,'new_CM_models_fitted_with_latest_FB_labels':True,'transport_descriptors_updated':True,
   'transport_plans_or_encoder_retrained':False,'new_segmentation':False,'new_biological_replicates':0,'segmentation_uncertainty_fully_propagated':False}
 (OUT/'release_validation.json').write_text(json.dumps(result,indent=2)+'\n');print('REFIT_RELEASE_PASS',json.dumps(result),flush=True)
if __name__=='__main__':main()
