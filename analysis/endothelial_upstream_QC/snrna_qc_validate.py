from pathlib import Path
import json,hashlib,os
os.environ.setdefault('MPLCONFIGDIR','/workspace/heart-atlas/cache/matplotlib')
import numpy as np
import pandas as pd
import scanpy as sc
import h5py
ROOT=Path('/workspace/heart-atlas');OUT=ROOT/'results/snrna_qc'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
summary=pd.read_csv(OUT/'qc_summary.tsv',sep='\t')
audit=pd.read_csv(OUT/'all_barcode_decisions.csv.gz',index_col=0)
assert len(summary)==8 and len(audit)==98371 and audit.index.is_unique
assert (summary.decontx_convergence<=0.001).all()
data=sc.read_h5ad(OUT/'snrna_qc_filtered.h5ad',backed='r')
assert data.n_obs==int(summary.final_nuclei.sum())==int(audit.retained.sum())
assert set(data.obs_names)==set(audit.index[audit.retained])
assert data.obs.scdblfinder_class.eq('singlet').all()
assert (data.obs.total_counts>=1000).all()
assert data.obs.n_genes_by_counts.between(750,7500).all()
assert (data.obs.pct_counts_mt<=5).all() and (data.obs.decontx_contamination<=.5).all()
assert np.isfinite(data.obsm['X_umap']).all()
assert (OUT/'QC_REPORT.pdf').stat().st_size>10000
dimensions=list(data.shape);data.file.close()
with h5py.File(OUT/'snrna_qc_filtered.h5ad') as f:
    assert f['layers/counts_raw/data'].dtype.kind in 'iu'
    for key in ['X/data','layers/counts_raw/data']:
        values=f[key]
        for start in range(0,len(values),1_000_000):
            chunk=values[start:start+1_000_000]
            assert np.isfinite(chunk).all() and (chunk>=0).all()
for p in OUT.glob('GSM*/final_summary.json'):
    text=(ROOT/'logs'/f'{p.parent.name}.decontx.log').read_text()
    assert 'DECONTX_DONE' in text and 'Execution halted' not in text
original=json.loads((ROOT/'metadata/downloaded_files.json').read_text())
source=[]
for entry in original:
    if not entry['file'].startswith('data/snrna/'):continue
    path=ROOT/entry['file'];assert path.stat().st_size==entry['bytes']
    assert sha(path)==entry['sha256'],'Input changed: '+entry['file']
    source.append(entry)
assert len(source)==24
artifacts=[]
for p in [OUT/'snrna_qc_filtered.h5ad',OUT/'QC_REPORT.pdf',OUT/'QC_REPORT.md',OUT/'qc_summary.tsv',OUT/'all_barcode_decisions.csv.gz']:
    artifacts.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
scripts=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in sorted((ROOT/'scripts').glob('snrna*')) if p.is_file()]
result={'status':'passed','dimensions':dimensions,'original_barcodes_accounted_for':len(audit),
 'input_files_verified_unchanged':source,'outputs':artifacts,'scripts':scripts,
 'checks':['all retained nuclei pass thresholds','all retained nuclei predicted singlets',
 'eight DecontX fits converged','unique IDs and complete retention/exclusion accounting',
 'finite nonnegative corrected and raw matrices','integer original counts retained','finite QC UMAP']}
(OUT/'validation_and_provenance.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k not in ['input_files_verified_unchanged','scripts']},indent=2))
