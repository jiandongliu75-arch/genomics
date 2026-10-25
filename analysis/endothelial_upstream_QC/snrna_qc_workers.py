"""Run the remaining per-library R fits with bounded parallelism."""
from pathlib import Path
import concurrent.futures,os,subprocess,argparse
ROOT=Path('/workspace/heart-atlas')
parser=argparse.ArgumentParser();parser.add_argument('samples',nargs='+');args=parser.parse_args()
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def fit(sample):
    p=ROOT/'results/snrna_qc'/sample
    assert (p/'stage1_summary.json').exists(),sample+' initial QC incomplete'
    if (p/'decontx_R_sessionInfo.txt').exists() and (p/'decontx_metrics.csv').exists():
        return sample+' REUSED'
    with (ROOT/'logs'/f'{sample}.decontx.log').open('w') as log:
        subprocess.run([str(ROOT/'r-env/bin/Rscript'),str(ROOT/'scripts/snrna_decontx.R'),str(p)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    return sample+' R_FIT_COMPLETE'
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    for future in concurrent.futures.as_completed([executor.submit(fit,s) for s in args.samples]):
        print(future.result(),flush=True)
