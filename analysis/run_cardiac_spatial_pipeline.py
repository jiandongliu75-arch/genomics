"""One-command reproducible execution with explicit optional-method boundaries."""
from pathlib import Path
import argparse, json, os, sys, subprocess, hashlib, importlib.metadata, time

ROOT=Path(__file__).resolve().parents[1]
STEPS=[('atlas','cardiac_spatial_pipeline.py',['--step','atlas']),
       ('geometry','cardiac_spatial_pipeline.py',['--step','geometry']),
       ('neighborhoods','cardiac_spatial_pipeline.py',['--step','neighborhoods']),
       ('all_lineage_QC','cardiac_spatial_pipeline_qc.py',['--step','audit']),
       ('niches','cardiac_spatial_niches.py',['--step','niches']),
       ('contact_enrichment','cardiac_spatial_niches.py',['--step','enrichment']),
       ('transport_pilot','cardiac_spatial_niches.py',['--step','transport']),
       ('models','cardiac_spatial_models.py',['--step','models']),
       ('nonlinear_prediction','cardiac_spatial_models.py',['--step','nonlinear']),
       ('QC_sensitivity','cardiac_spatial_pipeline_qc.py',['--step','sensitivity']),
       ('signaling_hypotheses','cardiac_spatial_models.py',['--step','mechanisms']),
       ('report','report_cardiac_spatial_pipeline.py',[]),
       ('validate','validate_cardiac_spatial_pipeline.py',[])]

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',default=str(ROOT/'analysis/cardiac_spatial_pipeline.json'))
    p.add_argument('--resume',action='store_true',help='Reuse a fully validated run only when configuration, inputs, code and artifacts match hashes')
    p.add_argument('--steps',nargs='+',choices=[s[0] for s in STEPS],help='Explicit partial execution; prerequisite artifacts must already exist')
    args=p.parse_args();config=Path(args.config).resolve();cfg=json.loads(config.read_text());out=ROOT/cfg['output']
    if not out.resolve().is_relative_to((ROOT/'results').resolve()):raise ValueError('Output must be a derived directory under this repository results/.')
    if cfg['niche_backend']!='sklearn_gmm':raise ValueError('Primary backend is sklearn_gmm. Use cardiac_spatial_niches.py --step cellcharter for the explicit optional benchmark.')
    if cfg['prediction_buffer_um']<2*cfg['proximity_primary_um']:raise ValueError('Primary prediction buffer must be at least twice the exposure radius.')
    if args.resume and not args.steps and (out/'validation.json').exists():
        result=json.loads((out/'validation.json').read_text())
        if result.get('status')=='PASS' and result.get('configuration_sha256')==sha(config):
            import csv
            with open(out/'file_manifest.tsv') as f:artifacts=list(csv.DictReader(f,delimiter='\t'))
            with open(out/'source_checksums.tsv') as f:inputs=list(csv.DictReader(f,delimiter='\t'))
            if all((ROOT/r['path']).exists() and sha(ROOT/r['path'])==r['sha256'] for r in artifacts+inputs):
                print('Reusing fully validated results; configuration, inputs, code and artifacts match.',flush=True);return
        raise RuntimeError('Validated-run hashes differ. Run without --resume to recompute, or select a new output directory in a copied configuration.')
    out.mkdir(parents=True,exist_ok=True)
    versions={}
    for package in ['numpy','pandas','scipy','scikit-learn','anndata','statsmodels','matplotlib','pillow','h5py']:
        versions[package]=importlib.metadata.version(package)
    provenance={'python':sys.version,'executable':sys.executable,'configuration_sha256':sha(config),'versions':versions,
        'source_repository_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'worktree_status':subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        'optional_methods_executed_by_default':False}
    (out/'runtime_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    env={**os.environ,'CARDIAC_SPATIAL_CONFIG':str(config),'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2'}
    receipts=[]
    for name,script,options in STEPS:
        if args.steps and name not in args.steps:continue
        print('PIPELINE_STEP',name,flush=True);start=time.time()
        subprocess.run([sys.executable,str(ROOT/'analysis'/script),*options],cwd=ROOT,env=env,check=True)
        receipts.append({'step':name,'script':script,'seconds':time.time()-start,'status':'completed'})
    # Receipts are runtime logs, deliberately excluded from immutable artifact hashes.
    (out/'execution_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':main()
