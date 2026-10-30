"""Run and package subtype-first cardiac neighborhood analyses."""
from pathlib import Path
import sys,os,argparse,subprocess,json,hashlib,csv,zipfile,re,time,importlib.metadata
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_subtype_niches'
SCRIPTS=['cardiac_subtype_niches.py','cardiac_subtype_transport.py','cardiac_subtype_figures.py','write_cardiac_subtype_story.py','validate_cardiac_subtypes.py','finalize_cardiac_subtype_documents.py','run_cardiac_subtype_story.py',
 'cardiac_story_analysis.py','cardiac_story_figures.py','write_cardiac_story.py','requirements-cardiac-reconstruction.txt','requirements-cardiac-story-document.txt']
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def release():
    for name in ['numerical_validation','document_validation']:assert json.loads((OUT/(name+'.json')).read_text())['status']=='PASS'
    sources=set(json.loads((OUT/'source_files.json').read_text()))|{str(ROOT/'analysis'/s) for s in SCRIPTS}
    # Include all original graph/cell inputs and saved molecular cloud sources
    # explicitly; this avoids depending on process completion order of registries.
    for stage in ['P0','P7','P14','P21']:
        for root,files in [('results/cardiac_spatial_pipeline',['atlas_annotations.tsv.gz','proximity50_graph.npz','proximity30_graph.npz','contact_graph.npz']),
                          ('results/cardiac_spatial_reconstruction',['CM_outcomes_and_anatomy.tsv.gz']),
                          ('results/cardiac_niche_dynamics',['nonCM_state_points.npz','state_cloud_membership.npz'])]:
            sources|={str(ROOT/root/stage/name) for name in files}
    with open(OUT/'source_checksums.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
        for name in sorted(sources):
            p=Path(name);w.writerow({'path':str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
    with zipfile.ZipFile(OUT/'Subtype_resolved_figure_source_tables.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(OUT.glob('*')):
            if p.suffix in ['.tsv','.json'] and p.name not in ['file_manifest.tsv','release_validation.json','runtime_provenance.json']:z.write(p,p.name)
    with zipfile.ZipFile(OUT/'Subtype_resolved_publication_figures.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted((OUT/'figures').glob('*')):z.write(p,p.relative_to(OUT))
    count=0
    for name in ['README.md','Manuscript_text.md']:
        p=OUT/name
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http://','https://','#')):continue
            target=p.parent/link.split('#')[0]
            if target.name!='release_validation.json':assert target.exists(),(name,link)
            count+=1
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json','runtime_provenance.json']]
    assert all(p.stat().st_size<95*2**20 for p in files)
    with open(OUT/'file_manifest.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
        for p in sorted(files):w.writerow({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
    (OUT/'release_validation.json').write_text(json.dumps({'status':'PASS','known_subtype_analysis_does_not_require_coarse_niche_cluster_recovery':True,
        'numerical_and_document_validation':'PASS','source_artifact_hashes_recorded':True,'relative_document_links_checked':count,'new_biological_replicates':0,
        'new_causal_mechanism_or_fate_validated':False},indent=2)+'\n')
    for p in OUT.glob('*.zip'):
        with zipfile.ZipFile(p) as z:assert z.testzip() is None
    print('SUBTYPE_RELEASE_VALIDATION_PASS',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--steps',nargs='+',choices=['analysis','transport','figures','numerical','word','render','documents','release']);p.add_argument('--document-python',default='python');p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.resume:
        assert json.loads((OUT/'release_validation.json').read_text())['status']=='PASS'
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(OUT/name) as f:
                for r in csv.DictReader(f,delimiter='\t'):assert sha(ROOT/r['path'])==r['sha256'],r['path']
        print('Subtype-focused release reused after input/code/artifact hash checks');return
    commands={'analysis':[sys.executable,'analysis/cardiac_subtype_niches.py'],'transport':[sys.executable,'analysis/cardiac_subtype_transport.py'],
        'figures':[sys.executable,'analysis/cardiac_subtype_figures.py'],'numerical':[sys.executable,'analysis/validate_cardiac_subtypes.py'],
        'word':[a.document_python,'analysis/write_cardiac_subtype_story.py'],'render':['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/'Subtype_resolved_cardiac_development_manuscript.docx')],
        'documents':[a.document_python,'analysis/finalize_cardiac_subtype_documents.py']}
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/cardiac-subtype-mpl','XDG_CACHE_HOME':'/tmp/cardiac-subtype-cache'}
    steps=a.steps or ['analysis','transport','figures','numerical','word','render','documents','release'];receipts=[]
    for name in steps:
        if name=='release':continue
        print('SUBTYPE_STEP',name,flush=True);start=time.time();subprocess.run(commands[name],cwd=ROOT,env=env,check=True);receipts.append({'step':name,'seconds':time.time()-start,'completed':True})
    path=OUT/'runtime_provenance.json';old=json.loads(path.read_text()).get('receipts',[]) if path.exists() else []
    path.write_text(json.dumps({'analysis_python':sys.executable,'document_python':a.document_python,'analysis_versions':{name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','matplotlib','anndata']},
        'starting_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'receipts':old+receipts,
        'primary_unit':'stage × CM identity × within-state maturation × known FB/EC subtypes; no two-cluster filter'},indent=2)+'\n')
    if 'release' in steps:release()

if __name__=='__main__':main()
