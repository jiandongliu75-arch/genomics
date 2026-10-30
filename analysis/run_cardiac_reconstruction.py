"""One-command paper-definition reconstruction and anatomy-led analysis."""
from pathlib import Path
import argparse,json,os,sys,subprocess,time,importlib.metadata,csv,re
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from reconstruct_cardiac_spatial import OUT,CFG,sha

STEPS=[('prepare','reconstruct_cardiac_spatial.py',['--step','all'],'python'),
 ('role_audit','audit_cardiac_source_roles.py',[],'python'),
 ('interface_sensitivity','cardiac_interface_sensitivity.py',[],'python'),
 ('RNA','cardiac_developmental_ordering.py',[],'python'),
 ('descriptive','test_cardiac_anatomy_mechanisms.py',['--step','descriptive'],'python'),
 ('R_LOESS','reconstruct_cardiac_loess.R',[str(OUT)],'R'),
 ('models','test_cardiac_anatomy_mechanisms.py',['--step','models'],'python'),
 ('matched','test_cardiac_anatomy_mechanisms.py',['--step','matched'],'python'),
 ('contact_definitions','test_cardiac_anatomy_mechanisms.py',['--step','contacts'],'python'),
 ('nuclear_cycling','check_cardiac_cycling_specificity.py',[],'python'),
 ('report','report_cardiac_reconstruction.py',[],'python'),
 ('validate','validate_cardiac_reconstruction.py',[],'python')]
SCRIPTS=[ROOT/'analysis'/name for _,name,_,_ in STEPS]+[ROOT/'analysis/run_cardiac_reconstruction.py',ROOT/'analysis/cardiac_spatial_reconstruction.json',ROOT/'analysis/requirements-cardiac-reconstruction.txt']

def release():
    from PIL import Image
    check=json.loads((OUT/'numerical_validation.json').read_text());assert check['status']=='PASS'
    text=subprocess.check_output(['pdfinfo',str(OUT/'cardiac_spatial_reconstruction_report.pdf')],text=True);assert re.search(r'Pages:\s+8\b',text)
    figures=list(OUT.glob('*.png'));assert len(figures)==10
    for path in figures:
        with Image.open(path) as im:im.verify()
    links=0;deferred_links=[]
    for doc in [OUT/'README.md',OUT/'EXPERIMENT_PLAN.md',OUT/'AUTHOR_INPUTS_NOT_PUBLIC.md']:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text()):
            if target.startswith(('http://','https://','#')):continue
            linked=doc.parent/target.split('#')[0]
            if linked==OUT/'release_validation.json':deferred_links.append(linked)
            else:assert linked.exists(),(doc,target)
            links+=1
    source={Path(p) for p in json.loads((OUT/'source_files.json').read_text())}|set(SCRIPTS)
    source|={Path('/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'),
        ROOT/'results/endothelial_QC_subtypes/snRNA/EC_annotations.tsv.gz',ROOT/'results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz',
        ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',ROOT/'analysis/cardiac_spatial_pipeline.json'}
    for stage in CFG['stages']:
        source|={ROOT/f'results/cardiac_spatial_pipeline/{stage}/{name}_graph.npz' for name in ['proximity30','proximity50','contact']}
        source|={ROOT/f'results/cardiac_niche_dynamics/{stage}/{name}.npz' for name in ['native_boundary_1um_graph','pruned_Delaunay_graph']}
        source|={ROOT/f'results/xenium_ventricular_mapping/{stage}/DAPI_overview.png',ROOT/f'results/xenium_ventricular_mapping/{stage}/morphology_display.json'}
    source|={Path('/workspace/heart-atlas/reference/PIP_seq/Figure4/4.2_main_figure_generation.R'),Path('/workspace/heart-atlas/metadata/paper.xml'),Path('/workspace/heart-atlas/data/author_webatlas/P7/obs/.zattrs')}
    sources=[]
    for path in sorted(source):
        name=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        sources.append({'path':name,'bytes':path.stat().st_size,'sha256':sha(path)})
    import pandas as pd
    pd.DataFrame(sources).to_csv(OUT/'source_checksums.tsv',sep='\t',index=False)
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','release_validation.json','runtime_provenance.json','execution_receipts.json']]
    assert all(p.stat().st_size<95*2**20 for p in files)
    pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
    (OUT/'release_validation.json').write_text(json.dumps({'status':'PASS','configuration_sha256':sha(ROOT/'analysis/cardiac_spatial_reconstruction.json'),
        'numerical_validation_status':check['status'],'PDF_pages':8,'figures':len(figures),'documentation_links_checked':links,
        'source_and_artifact_hashes_recorded':True,'GitHub_file_size_limits_passed':True,'new_observed_samples_or_spatial_replicates_generated':False,
        'exact_author_input_recovery_claimed':False},indent=2)+'\n')
    assert all(path.exists() for path in deferred_links)
    print('RELEASE VALIDATION PASS',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');p.add_argument('--steps',nargs='+',choices=[s[0] for s in STEPS]+['release']);p.add_argument('--Rscript',default='/workspace/heart-atlas/r-env/bin/Rscript');a=p.parse_args()
    if a.resume:
        status=json.loads((OUT/'release_validation.json').read_text());assert status['status']=='PASS' and status['configuration_sha256']==sha(ROOT/'analysis/cardiac_spatial_reconstruction.json')
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(OUT/name) as f:rows=list(csv.DictReader(f,delimiter='\t'))
            for row in rows:
                path=ROOT/row['path']
                if not path.exists() or sha(path)!=row['sha256']:raise RuntimeError('Validated source/code/artifact changed: '+str(path))
        print('Reusing verified reconstruction: inputs, code, configuration and artifacts match hashes.');return
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/cardiac-reconstruction-mpl'}
    provenance={'python':sys.version,'executable':sys.executable,'versions':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','scikit-learn','anndata','matplotlib','statsmodels','shapely','scikit-image','pillow']},
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'configuration_sha256':sha(ROOT/'analysis/cardiac_spatial_reconstruction.json'),
        'steps_this_invocation':a.steps or [s[0] for s in STEPS]+['release'],'Rscript':a.Rscript,
        'interpretation':'Observed spatial data at four ages; conditional RNA interpolation cannot add new observed ages, biological replication or causal evidence'}
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'runtime_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');receipts=[]
    for name,script,options,engine in STEPS:
        if a.steps and name not in a.steps:continue
        print('RECONSTRUCTION_STEP',name,flush=True);start=time.time();cmd=[sys.executable] if engine=='python' else [a.Rscript,'--vanilla']
        subprocess.run(cmd+[str(ROOT/'analysis'/script),*options],cwd=ROOT,env=env,check=True);receipts.append({'step':name,'seconds':time.time()-start,'status':'completed'})
    if not a.steps or 'release' in a.steps:release()
    (OUT/'execution_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':main()
