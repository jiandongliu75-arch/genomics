from pathlib import Path
import json,hashlib,pandas as pd,numpy as np,h5py
import anndata as ad
from anndata.io import read_elem
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/cell_identities_subtypes'
full=Path('/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
ann=pd.read_csv(out/'all_cell_identity_and_subtype_annotations.csv.gz',index_col=0)
s=pd.read_csv(out/'subcluster_identity_evidence.tsv',sep='\t')
assert len(ann)==73052 and ann.index.is_unique and ann.reviewed_identity.notna().all() and ann.proposed_subtype.notna().all()
assert len(s)==59 and s.subcluster_id.is_unique and int(s.nuclei.sum())==int(ann.subcluster_id.ne('not_subclustered').sum())
passed=s.confidence.eq('provisional-supported')
for key,value in [('nuclei',50),('libraries_with_ge10_nuclei',2),('positive_within_library_contrasts',2),('panel_score',.3),('score_margin',.25),('panel_enrichment',.1)]:assert (s.loc[passed,key]>=value).all(),key
assert (s.loc[passed,'largest_library_fraction']<=.8).all()
assert ann.loc[ann.proposed_subtype.str.contains('candidate',case=False),'subtype_review_flag'].all()
for parent,g in s.groupby('parent'):
 m=ad.read_h5ad(out/parent/'subclusters_metadata_only.h5ad')
 assert m.n_vars==0 and m.n_obs==int(g.nuclei.sum()) and np.isfinite(m.obsm['X_umap']).all()
 assert m.obs_names.equals(ann.loc[m.obs_names].index)
 assert m.obs.proposed_subtype.astype(str).equals(ann.loc[m.obs_names,'proposed_subtype'].astype(str))
 assert m.obs.subtype_review_flag.equals(ann.loc[m.obs_names,'subtype_review_flag'])
for path in [out/'all_cell_annotations_metadata_only.h5ad',full]:
 with h5py.File(path) as f:
  o=read_elem(f['obs']);assert o.index.equals(ann.index)
  for col in ['reviewed_identity','subcluster_id','proposed_subtype','subtype_confidence','subtype_review_flag']:
   assert o[col].astype(str).equals(ann[col].astype(str)),str(path)+col
  if path==full:
   assert tuple(f['X'].attrs['shape'])==(73052,32285) and tuple(f['layers/counts_raw'].attrs['shape'])==(73052,32285)
   assert f['layers/counts_raw/data'].dtype.kind in 'iu'
record=json.loads((out/'full_dataset_release.json').read_text());record['sha256']=sha(full);record['bytes']=full.stat().st_size;(out/'full_dataset_release.json').write_text(json.dumps(record,indent=2))
manifest=json.loads((out/'manifest.json').read_text());manifest['criteria_sha256']=sha(ROOT/'analysis/subtype_criteria.json');manifest['scripts']=[{'file':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in sorted((ROOT/'analysis').glob('*subtype*')) if p.is_file()]
manifest['files']=[{'file':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'];assert all(f['bytes']<90_000_000 for f in manifest['files']);manifest['status']='validated after conservative nomenclature review';(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('PASSED: 73,052 identities accounted for; 59 subclusters; marker/library gates; metadata/full H5AD alignment; candidate review flags; checksums.')
print('Full H5AD:',record['bytes'],'bytes; SHA256',record['sha256'])
