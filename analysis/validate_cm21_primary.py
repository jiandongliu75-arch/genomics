"""Audit the sole21-gene outcome, complete cohort, source geometry and models."""
from pathlib import Path
import json,hashlib,zipfile,re,subprocess,shutil
import numpy as np
import pandas as pd
import h5py
from anndata.io import sparse_dataset,read_elem
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.stats import false_discovery_control,t
from sklearn.preprocessing import SplineTransformer
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM21_primary_spatial'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 method=json.loads((OUT/'methods.json').read_text());panel=pd.read_csv(OUT/'CM_maturation_21gene_panel.tsv',sep='\t').gene.tolist()
 assert len(panel)==len(set(panel))==21 and panel==method['genes'] and method['only21_CM_outcome_genes']
 assert not method['anatomical_exclusion'] and not method['published_label_concordance_required']
 cal=pd.read_csv(ROOT/'results/cardiac_spatial_reconstruction/score_reference_calibration.tsv',sep='\t').set_index(['view','scope','gene']).sort_index()
 pairs=pd.read_csv(OUT/'high_low_CM21_matched_pairs.tsv.gz',sep='\t');model=pd.read_csv(OUT/'continuous_CM21_associations.tsv',sep='\t')
 selected={'P0':'FB:Postn','P7':'EC:Capillary-like','P14':'FB:Col8a1','P21':'FB:Lsamp'}
 checks=[];total=0;inputs=set(json.loads((OUT/'input_paths.json').read_text()))
 for stage in ['P0','P7','P14','P21']:
  d=pd.read_csv(OUT/f'{stage}_CM21_scores_and_neighborhoods.tsv.gz',sep='\t',index_col=0)
  f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',sep='\t',index_col=0)
  keep=f.ventricular_CM_strict&f.atlas_identity.eq('Cardiomyocyte')&f.pass_expression_QC&f.core_transcripts.ge(20)&f.CM_outcome_excluded_depth.ge(20)
  assert d.index.equals(f.index[keep]);total+=len(d)
  assert d.index.is_unique and d.core_score21.notna().all() and d.native_score21.notna().all()
  sample=np.sort(np.random.default_rng(281029).choice(len(d),min(256,len(d)),replace=False));score_errors={}
  for view in ['core','native']:
   p=ROOT/(f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
   with h5py.File(p) as h:
    ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index;assert genes.get_indexer(panel).min()>=0
    X=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[ids.get_indexer(d.index[sample])].tocsr()
    total_count=np.asarray(X.sum(1)).ravel();xx=X[:,genes.get_indexer(panel)].toarray();ref=cal.loc[(view,'all')].reindex(panel)
    expected=((np.log1p(xx*10000/np.maximum(total_count,1)[:,None])-ref.reference_mean.to_numpy())/ref.reference_SD.to_numpy()).mean(1)
    error=float(np.max(abs(expected-d[view+'_score21'].iloc[sample])));assert error<1e-7;score_errors[view]=error
   rank=d[view+'_stage_rank'];group=np.where(rank<=.25,'Lower quartile',np.where(rank>=.75,'Higher quartile','Middle half'))
   np.testing.assert_array_equal(group,d[view+'_score_group'])
   pp=pairs[pairs.stage.eq(stage)&pairs.view.eq(view)]
   assert pp.high_CM.is_unique and pp.low_CM.is_unique and set(pp.high_CM).isdisjoint(pp.low_CM)
   hi=d.loc[pp.high_CM];lo=d.loc[pp.low_CM]
   assert hi[view+'_stage_rank'].ge(.75).all() and lo[view+'_stage_rank'].le(.25).all()
   assert hi.CM_state.to_list()==lo.CM_state.to_list() and hi.nucleus_missing.to_list()==lo.nucleus_missing.to_list()
   np.testing.assert_array_equal((hi[['x_centroid','y_centroid']].to_numpy()//200).astype(int),(lo[['x_centroid','y_centroid']].to_numpy()//200).astype(int))
   physical=np.linalg.norm(hi[['x_centroid','y_centroid']].to_numpy()-lo[['x_centroid','y_centroid']].to_numpy(),axis=1)
   assert physical.max()<100+1e-4
   for col in ['CM_outcome_excluded_depth','cell_area']:
    scale=np.log1p(d[col]).std(ddof=0);diff=(np.log1p(hi[col]).to_numpy()-np.log1p(lo[col]).to_numpy())/scale;assert abs(diff).max()<.25+1e-6
  source=selected[stage]
  if source.startswith('FB:'):
   fb=pd.read_csv(ROOT/f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',sep='\t',index_col=0)
   mask=pd.Series(False,index=f.index);mask.loc[fb.index]=fb.program_source_eligible&fb.recommended_candidate.eq(source[3:]+'-associated');mask &= f.atlas_identity.eq('Fibroblast')&f.program_source_eligible
  else:mask=f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')
  assert not d.index.isin(f.index[mask]).any()
  distance=cKDTree(f.loc[mask,['x_centroid','y_centroid']]).query(d[['x_centroid','y_centroid']])[0]
  np.testing.assert_allclose(distance,d['distance__'+source],atol=2e-4)
  G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()[f.index.get_indexer(d.index)]
  count=np.asarray(G@mask.to_numpy(float)).ravel();np.testing.assert_array_equal(count,d['neighbor_count50__'+source]);np.testing.assert_array_equal(count>0,d['presence50__'+source])
  e=-np.log1p(distance);e=(e-e.mean())/e.std();y=d.core_score21.to_numpy();y=(y-y.mean())/y.std()
  C=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(d[['CM_outcome_excluded_depth','cell_area']]))
  C=np.column_stack([C,d.nucleus_missing.astype(float),pd.get_dummies(d.CM_state,dtype=float)]);C=(C-C.mean(0))/np.maximum(C.std(0),1e-9)
  group=(d.x_centroid//200).astype(int).astype(str)+'_'+(d.y_centroid//200).astype(int).astype(str)+'::'+d.CM_state
  frame=pd.DataFrame(np.column_stack([C,e,y]),index=d.index);centered=frame-frame.groupby(group).transform('mean');X=centered.iloc[:,:-1].to_numpy();Y=centered.iloc[:,-1].to_numpy()
  coef=np.linalg.lstsq(X,Y,rcond=1e-10)[0];beta=coef[-1]
  row=model[model.stage.eq(stage)&model.source.eq(source)&model.primary_analysis].iloc[0]
  assert abs(beta-row.beta_SD)<1e-7
  residual=Y-X@coef;er=X[:,-1]-X[:,:-1]@np.linalg.lstsq(X[:,:-1],X[:,-1],rcond=1e-10)[0]
  block=(d.x_centroid//400).astype(int).astype(str)+'_'+(d.y_centroid//400).astype(int).astype(str)
  meat=pd.Series(er*residual,index=d.index).groupby(block).sum().to_numpy();ng=len(meat);rank=np.linalg.matrix_rank(X)
  se=np.sqrt(np.square(meat).sum()/(er@er)**2*ng/(ng-1)*(len(d)-1)/(len(d)-rank));CI=beta+t.ppf(.975,ng-1)*np.array([-se,se]);np.testing.assert_allclose(CI,[row.CI_low,row.CI_high],atol=1e-7)
  checks.append(dict(stage=stage,CMs=len(d),independent_21gene_errors=score_errors,independent_source_distance_and_counts_verified=True,independent_model_error=abs(beta-row.beta_SD),cluster_interval_verified=True))
 assert total==131605
 primary=model.primary_analysis;assert primary.sum()==48
 np.testing.assert_allclose(false_discovery_control(model.loc[primary,'within_section_p']),model.loc[primary,'primary_BH_q'],rtol=1e-7)
 m=pd.read_csv(OUT/'high_low_CM21_matched_comparisons.tsv',sep='\t');sel=m.distance_comparison
 np.testing.assert_allclose(false_discovery_control(m.loc[sel,'within_section_p']),m.loc[sel,'matched_distance_BH_q'],rtol=1e-7)
 assert m.loc[sel,'matched_distance_BH_q'].ge(.05).all()
 # Independently verify every matched mean contrast from the released per-cell values.
 for stage in ['P0','P7','P14','P21']:
  d=pd.read_csv(OUT/f'{stage}_CM21_scores_and_neighborhoods.tsv.gz',sep='\t',index_col=0)
  for view in ['core','native']:
   pp=pairs[pairs.stage.eq(stage)&pairs.view.eq(view)];hi=d.index.get_indexer(pp.high_CM);lo=d.index.get_indexer(pp.low_CM)
   for row in m[m.stage.eq(stage)&m.view.eq(view)].itertuples():
    col=('distance__' if row.metric=='closer_um' else 'presence50__')+row.source;v=d[col].to_numpy(float)
    point=(v[lo]-v[hi]).mean() if row.metric=='closer_um' else 100*(v[hi]-v[lo]).mean()
    # Long distances are rounded to9 significant digits in portable CSVs.
    assert abs(point-row.estimate)<1e-6
 docs=json.loads((OUT/'document_validation.json').read_text());assert docs['status']=='PASS' and len(docs['figure_pages'])==3
 html=(OUT/'CM21_spatial_viewer.html').read_text();payload=json.loads(html.split('const dataset=',1)[1].split(',limits=',1)[0]);assert sum(len(v['cm']) for v in payload.values())==131605
 if shutil.which('node'):
  js=re.search(r'<script>(.*?)</script>',html,re.S).group(1);temp=Path('/tmp/cm21_viewer_syntax.js');temp.write_text(js);subprocess.run(['node','--check',str(temp)],check=True);temp.unlink()
 pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted({Path(x) if Path(x).is_absolute() else ROOT/x for x in inputs})]).to_csv(OUT/'input_checksums.tsv',sep='\t',index=False)
 result=dict(status='PASS',sole_CM_maturation_genes=panel,number_CM_outcome_genes=21,all_QC_eligible_ventricular_CMs=total,
  no_anatomical_or_public_label_exclusion=True,stage_checks=checks,primary_full_cohort_source_tests=48,
  matched_core_pairs=int(pairs.view.eq('core').sum()),all_matched_means_independently_verified=True,
  viewer_CMs_verified=131605,viewer_JavaScript_syntax_checked=bool(shutil.which('node')),documents=docs,new_biological_replicates=0)
 (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
 files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='file_manifest.tsv')
 pd.DataFrame([dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
