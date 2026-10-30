"""Reproduce the completed refit and report, or verify a hashed release."""
from pathlib import Path
import os,argparse,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_fibroblast_CM_refit'
def main():
 p=argparse.ArgumentParser();p.add_argument('--document-python',default='python');p.add_argument('--resume',action='store_true');a=p.parse_args()
 env={**os.environ,'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/fibroblast-CM-mpl','XDG_CACHE_HOME':'/tmp/clear-heart-cache'}
 def run(command):subprocess.run(command,cwd=ROOT,env=env,check=True)
 if a.resume:run([a.document_python,'analysis/finalize_fibroblast_CM_refit.py','--resume']);return
 run([sys.executable,'analysis/refit_fibroblast_CM_relationships.py','--step','all'])
 run([sys.executable,'analysis/fibroblast_CM_refit_checks.py','--step','all'])
 run([sys.executable,'analysis/fibroblast_CM_refit_extensions.py'])
 run([sys.executable,'analysis/fibroblast_CM_gene_exclusion_sensitivity.py'])
 run([sys.executable,'analysis/plot_fibroblast_CM_refit.py'])
 run([a.document_python,'analysis/write_fibroblast_CM_refit_report.py'])
 run(['soffice','--headless','--convert-to','pdf','--outdir',str(OUT),str(OUT/'Fibroblast_CM_spatial_reassessment.docx')])
 run([a.document_python,'analysis/finalize_fibroblast_CM_refit.py'])
if __name__=='__main__':main()
