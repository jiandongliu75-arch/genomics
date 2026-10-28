"""Audit masks, source support, score preservation and independent model fits."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
from scipy.stats import false_discovery_control,t
from sklearn.preprocessing import SplineTransformer
from shapely import Polygon,contains_xy
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/morphology_stratified_nonCM'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 freeze=json.loads((OUT/'anatomy_freeze.json').read_text());design=json.loads((OUT/'association_design.json').read_text())
 assert sha(OUT/'morphology_ROIs.geojson')==freeze['mask_sha256']==design['anatomy_sha256']
 assert not freeze['expression_inputs_used'] and not freeze['outcomes_loaded_by_annotation_script']
 rois=json.loads((OUT/'morphology_ROIs.geojson').read_text())['features'];assert len(rois)==31
 assert all(x['properties']['anatomical_chamber']=='Unresolved' for x in rois)
 assert sum(x['properties']['primary_eligible'] for x in rois)==23
 for source in freeze['source_images']:assert sha(ROOT/source['DAPI_path'])==source['DAPI_sha256']
 r=pd.read_csv(OUT/'nonCM_maturation_associations.tsv',sep='\t');primary=r.primary_analysis
 assert primary.sum()==design['primary_tests']==48
 np.testing.assert_allclose(false_discovery_control(r.loc[primary,'within_section_p']),r.loc[primary,'primary_BH_q'],rtol=1e-7)
 secondary=r.secondary_BH_q.notna();np.testing.assert_allclose(false_discovery_control(r.loc[secondary,'within_section_p']),r.loc[secondary,'secondary_BH_q'],rtol=1e-7)
 assert r.spatial_blocks.ge(8).all() and np.isfinite(r.beta_SD).all()
 models=[];cells=[]
 selected={'P0':'FB:Postn','P7':'EC:Capillary-like','P14':'Broad:Immune','P21':'FB:Col8a1'}
 for stage in ['P0','P7','P14','P21']:
  d=pd.read_csv(OUT/f'{stage}_CM_region_annotations_and_scores.tsv.gz',sep='\t',index_col=0)
  old=pd.read_csv(ROOT/f'results/capillary_maturation_retest/{stage}_CM_retest_values.tsv.gz',sep='\t',index_col=0);old=old[old.strict]
  assert d.index.equals(old.index)
  np.testing.assert_allclose(d.core_maturation21,old.core_maturation21,atol=1e-7)
  np.testing.assert_allclose(d.native_maturation21,old.native_maturation21,atol=1e-7)
  check=np.full(len(d),'Unclassified',dtype=object)
  for feature in rois:
   p=feature['properties']
   if p['stage']!=stage:continue
   poly=Polygon(feature['geometry']['coordinates'][0]);take=contains_xy(poly,d.x_centroid.to_numpy(),d.y_centroid.to_numpy())
   check[take]=p['ROI_id']
  assert np.array_equal(check,d.ROI.to_numpy()),stage
  assert np.array_equal(d.annotated_primary,d.confidence.isin(['high','moderate']))
  yfull=old.core_maturation21.to_numpy();np.testing.assert_allclose(d.Y__core_maturation21,(yfull-yfull.mean())/yfull.std(),atol=1e-7)
  cells.append(dict(stage=stage,strict_CMs=len(d),primary_anatomical_CMs=int(d.annotated_primary.sum()),independent_polygon_assignments_match=True))
  f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',sep='\t',index_col=0)
  source=selected[stage]
  if source=='FB:Postn' or source=='FB:Col8a1':
   fb=pd.read_csv(ROOT/f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',sep='\t',index_col=0)
   mask=pd.Series(False,index=f.index);mask.loc[fb.index]=fb.program_source_eligible&fb.recommended_candidate.eq(source[3:]+'-associated')
  elif source=='EC:Capillary-like':mask=f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')
  else:mask=f.atlas_identity.eq('Immune')&f.program_source_eligible
  g=sparse.load_npz(ROOT/f'results/cardiac_spatial_reconstruction/{stage}/tissue_screened_Delaunay_graph.npz').tocoo();xy=f[['x_centroid','y_centroid']].to_numpy()
  weighted=sparse.csr_matrix((np.sqrt(np.square(xy[g.row]-xy[g.col]).sum(1)),(g.row,g.col)),shape=g.shape)
  path=dijkstra(weighted,directed=False,indices=np.flatnonzero(mask),min_only=True)[f.index.get_indexer(d.index)]
  e=-np.log1p(path);finite=np.isfinite(e);e[finite]=(e[finite]-e[finite].mean())/e[finite].std();e[~finite]=np.nan
  z=d[d.annotated_primary&np.isfinite(e)].copy();z['e']=e[d.annotated_primary&np.isfinite(e)]
  C=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(z[['capture','cell_area']]))
  C=np.column_stack([C,z.nucleus_missing.astype(float),np.log1p(z[['distance_to_large_void_um','distance_to_outer_um']])]);C=(C-C.mean(0))/np.maximum(C.std(0),1e-8)
  groups=z.ROI+'::'+z.CM_state+'::'+(z.x_centroid//200).astype(int).astype(str)+'::'+(z.y_centroid//200).astype(int).astype(str)
  design_df=pd.DataFrame(np.column_stack([C,z.e,z.Y__core_maturation21]),index=z.index)
  centered=design_df-design_df.groupby(groups).transform('mean');X=centered.iloc[:,:-1].to_numpy();Y=centered.iloc[:,-1].to_numpy()
  beta=np.linalg.lstsq(X,Y,rcond=1e-10)[0][-1]
  row=r[primary&r.stage.eq(stage)&r.source.eq(source)].iloc[0];assert abs(beta-row.beta_SD)<1e-7
  # Independently reconstruct clustered confidence from the full-design residual.
  residual=Y-X@np.linalg.lstsq(X,Y,rcond=1e-10)[0];xr=X[:,-1]-X[:,:-1]@np.linalg.lstsq(X[:,:-1],X[:,-1],rcond=1e-10)[0]
  block=(z.x_centroid//400).astype(int).astype(str)+'_'+(z.y_centroid//400).astype(int).astype(str)
  score=pd.Series(xr*residual,index=z.index).groupby(block).sum().to_numpy();nb=len(score);rank=np.linalg.matrix_rank(X)
  se=np.sqrt((score*score).sum()/(xr@xr)**2*nb/(nb-1)*(len(z)-1)/(len(z)-rank))
  ci=beta+t.ppf(.975,nb-1)*np.array([-se,se]);np.testing.assert_allclose(ci,[row.CI_low,row.CI_high],atol=1e-7)
  models.append(dict(stage=stage,source=source,independent_coefficient_error=abs(beta-row.beta_SD),cluster_interval_verified=True))
 docs=json.loads((OUT/'document_validation.json').read_text());assert docs['status']=='PASS'
 sens=pd.read_csv(OUT/'targeted_score_definition_sensitivities.tsv',sep='\t');assert len(sens)==10
 np.testing.assert_allclose(false_discovery_control(sens.within_section_p),sens.score_sensitivity_BH_q,rtol=1e-7)
 assert sens.loc[sens.stage.eq('P0'),'CI_low'].gt(0).all()
 equal=sens[sens.stage.eq('P7')&sens.outcome.eq('S__equal_50_molecules')].iloc[0]
 assert equal.CI_low<0<equal.CI_high
 assignments=pd.read_csv(OUT/'nonCM_anatomical_assignments.tsv.gz',sep='\t')
 assert not assignments.duplicated(['stage','cell_id']).any()
 assert assignments.atlas_identity.isin(['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']).all()
 view=(OUT/'Morphology_only_ROI_review.html').read_text();assert 'const data=' in view and 'data:image/jpeg;base64,' in view
 inputs=set(json.loads((OUT/'analysis_input_paths.json').read_text()))
 inputs|={s['DAPI_path'] for s in freeze['source_images']}
 inputs|={f'results/CM_maturation21_nonCM_context/{s}/CM_maturation21_scores.tsv.gz' for s in ['P0','P7']}
 paths=[Path(p) if Path(p).is_absolute() else ROOT/p for p in inputs]
 pd.DataFrame([dict(path=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)]).to_csv(OUT/'input_checksums.tsv',sep='\t',index=False)
 result=dict(status='PASS',anatomy_sha256=freeze['mask_sha256'],ROI_count=31,primary_ROIs=23,
  all_LV_RV_unresolved=True,anatomy_cells=cells,independent_model_refits=models,primary_BH_tests=48,
  secondary_BH_tests=int(secondary.sum()),targeted_score_sensitivity_tests=10,nonCM_profiles_assigned=len(assignments),document_validation=docs,
  no_new_samples=True,no_new_segmentation=True,expert_histology_validation_completed=False)
 (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
 files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='file_manifest.tsv' and not p.name.startswith('.')]
 pd.DataFrame([dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
