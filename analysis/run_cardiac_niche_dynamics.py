"""Run or verify the paper-guided extension without overwriting the v1 atlas."""
from pathlib import Path
import argparse,json,os,sys,subprocess,time,importlib.metadata
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from validate_cardiac_niche_dynamics import sha

STEPS=[('vascular_prepare','revisit_vascular_coordination.py',['--step','prepare'],False),
 ('vascular_models','revisit_vascular_coordination.py',['--step','models'],False),
 ('RNA','revisit_vascular_coordination.py',['--step','rna'],False),
 ('architecture','revisit_vascular_coordination.py',['--step','architecture'],False),
 ('source_capture','revisit_vascular_coordination.py',['--step','source_capture'],False),
 ('sensitivity','revisit_vascular_coordination.py',['--step','sensitivity'],False),
 ('state_prepare','fit_wormhole_cardiac_niches.py',['--step','prepare'],False),
 ('Wormhole_train','fit_wormhole_cardiac_niches.py',['--step','fit'],True),
 ('Wormhole_assess','fit_wormhole_cardiac_niches.py',['--step','assess'],False),
 ('report','report_cardiac_niche_dynamics.py',[],False),
 ('validate','validate_cardiac_niche_dynamics.py',[],False)]

def main():
    p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true')
    p.add_argument('--steps',nargs='+',choices=[a[0] for a in STEPS])
    p.add_argument('--wormhole-python',default=os.environ.get('WORMHOLE_PYTHON',sys.executable))
    args=p.parse_args();cfgpath=ROOT/'analysis/cardiac_niche_dynamics.json';cfg=json.loads(cfgpath.read_text());out=ROOT/cfg['output']
    assert out.resolve().is_relative_to((ROOT/'results').resolve())
    if args.resume and not args.steps:
        validation=json.loads((out/'validation.json').read_text());assert validation['status']=='PASS' and validation['configuration_sha256']==sha(cfgpath)
        import csv
        for name in ['source_checksums.tsv','file_manifest.tsv']:
            with open(out/name) as f:rows=list(csv.DictReader(f,delimiter='\t'))
            for row in rows:
                path=ROOT/row['path']
                if not path.exists() or sha(path)!=row['sha256']:raise RuntimeError('Validated input/code/artifact changed: '+str(path))
        print('Reusing verified cardiac niche extension: configuration, source, code and artifact hashes match.');return
    prior=json.loads((ROOT/cfg['upstream']/'validation.json').read_text());assert prior['status']=='PASS'
    env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2',
         'XLA_FLAGS':'--xla_cpu_multi_thread_eigen=false --xla_force_host_platform_device_count=1','MPLCONFIGDIR':'/tmp/cardiac-dynamics-mpl'}
    versions={x:importlib.metadata.version(x) for x in ['numpy','scipy','pandas','scikit-learn','anndata','matplotlib','shapely']}
    cmd="import importlib.metadata as m,json; print(json.dumps({p:m.version(p) for p in ['wassersteinwormhole','jax','jaxlib','flax','ott-jax','optax','chex','orbax-checkpoint','equinox','lineax','clu']}))"
    wormhole_versions=json.loads(subprocess.check_output([args.wormhole_python,'-c',cmd],text=True,env=env))
    out.mkdir(parents=True,exist_ok=True)
    provenance={'python':sys.version,'primary_executable':sys.executable,'primary_versions':versions,'Wormhole_executable':args.wormhole_python,
        'Wormhole_versions':wormhole_versions,'configuration_sha256':sha(cfgpath),
        'source_repository_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'steps_this_invocation':args.steps or [a[0] for a in STEPS],
        'inherited_bootstrap':{'module':'cardiac_spatial_models.joint_fit','seed':731,'spatial_block_um':400,'draws':cfg['bootstrap_draws']},
        'previous_v1_release_unchanged':True}
    (out/'runtime_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    receipts=[]
    for name,script,options,wormhole in STEPS:
        if args.steps and name not in args.steps:continue
        start=time.time();print('NICHE_EXTENSION_STEP',name,flush=True)
        subprocess.run([args.wormhole_python if wormhole else sys.executable,str(ROOT/'analysis'/script),*options],cwd=ROOT,env=env,check=True)
        receipts.append({'step':name,'seconds':time.time()-start,'status':'completed'})
    (out/'execution_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':main()
