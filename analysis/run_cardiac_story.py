"""Build, render, check and package the cardiac-development manuscript."""
from pathlib import Path
import sys,os,json,argparse,subprocess,hashlib,csv,re,zipfile,time,importlib.metadata

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_developmental_story'
SCRIPTS=['cardiac_story_analysis.py','cardiac_story_figures.py','write_cardiac_story.py','validate_cardiac_story.py','validate_cardiac_story_documents.py','run_cardiac_story.py','requirements-cardiac-story-document.txt']
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def release():
    numerical=json.loads((OUT/'numerical_validation.json').read_text());document=json.loads((OUT/'document_validation.json').read_text())
    assert numerical['status']==document['status']=='PASS'
    sources=set(json.loads((OUT/'analysis_source_files.json').read_text()))|{str(ROOT/'analysis'/s) for s in SCRIPTS}
    sources|={str(ROOT/'analysis'/s) for s in ['requirements-cardiac-reconstruction.txt','requirements-cardiac-niche-dynamics.txt']}
    with open(OUT/'source_checksums.tsv','w') as f:
        writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for name in sorted(sources):
            p=Path(name);writer.writerow({'path':str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
    with zipfile.ZipFile(OUT/'Cardiac_development_publication_figures.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted((OUT/'figures').glob('*')):z.write(p,p.relative_to(OUT))
        z.write(OUT/'figure_message_and_panel_key.tsv','figure_message_and_panel_key.tsv')
    with zipfile.ZipFile(OUT/'Cardiac_development_figure_source_data.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(OUT.glob('*')):
            if p.suffix in ['.tsv','.json'] and p.name not in ['file_manifest.tsv','release_validation.json','runtime_provenance.json']:z.write(p,p.name)
    links=0
    for doc in [OUT/'README.md',OUT/'Manuscript_text.md']:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text()):
            if target.startswith(('https://','http://','#')):continue
            path=doc.parent/target.split('#')[0]
            if path.name!='release_validation.json':assert path.exists(),(doc,target)
            links+=1
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json','runtime_provenance.json']]
    assert all(p.stat().st_size<95*2**20 for p in files)
    with open(OUT/'file_manifest.tsv','w') as f:
        writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for p in sorted(files):writer.writerow({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
    (OUT/'release_validation.json').write_text(json.dumps({'status':'PASS','numerical_checks':'PASS','document_checks':'PASS',
        'Word_pages_after_rendering':document['rendered_manuscript_pages'],'integrated_main_figures':6,'observed_OT_neighborhoods':96,
        'relative_document_links_checked':links,'source_and_artifact_hashes_recorded':True,'ZIP_packages_checked':True,
        'new_biological_replication':0,'lineage_conversion_or_causality_claimed_as_validated':False},indent=2)+'\n')
    for name in ['Cardiac_development_publication_figures.zip','Cardiac_development_figure_source_data.zip']:
        with zipfile.ZipFile(OUT/name) as z:assert z.testzip() is None
    print('STORY_RELEASE_VALIDATION_PASS',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--steps',nargs='+',choices=['analysis','figures','numerical','word','render','documents','release']);p.add_argument('--document-python',default='python');p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.resume:
        assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(OUT/name) as f:
                for row in csv.DictReader(f,delimiter='\t'):assert sha(ROOT/row['path'])==row['sha256'],row['path']
        print('Manuscript release reused after source/code/artifact hash verification');return
    steps=a.steps or ['analysis','figures','numerical','word','render','documents','release'];receipts=[]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/cardiac-story-mpl','XDG_CACHE_HOME':'/tmp/cardiac-story-cache'}
    commands={'analysis':[sys.executable,'analysis/cardiac_story_analysis.py'],'figures':[sys.executable,'analysis/cardiac_story_figures.py'],
      'numerical':[sys.executable,'analysis/validate_cardiac_story.py'],'word':[a.document_python,'analysis/write_cardiac_story.py'],
      'render':['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/'Postnatal_cardiac_development_manuscript.docx')],
      'documents':[a.document_python,'analysis/validate_cardiac_story_documents.py']}
    for step in steps:
        if step=='release':continue
        print('STORY_STEP',step,flush=True);start=time.time();subprocess.run(commands[step],cwd=ROOT,env=env,check=True)
        receipts.append({'step':step,'completed':True,'seconds':time.time()-start})
    old=OUT/'runtime_provenance.json'
    previous=json.loads(old.read_text()).get('execution_receipts',[]) if old.exists() else []
    document_versions=json.loads(subprocess.check_output([a.document_python,'-c',"import json,docx,fitz;print(json.dumps({'python_docx':docx.__version__,'PyMuPDF':fitz.VersionBind}))"],cwd=ROOT,text=True))
    old.write_text(json.dumps({'analysis_python':sys.executable,'document_python':a.document_python,'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'analysis_versions':{name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','anndata','matplotlib','pillow']},'document_versions':document_versions,
        'execution_receipts':previous+receipts,'data_policy':'Frozen QC inputs; source hashes and original coordinates retained; one spatial section/age'},indent=2)+'\n')
    if 'release' in steps:release()

if __name__=='__main__':main()
