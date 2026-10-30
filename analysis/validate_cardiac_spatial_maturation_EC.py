"""Independent arithmetic/geometry/source-preservation checks for spatial tests."""
from cardiac_spatial_maturation_EC import *
from PIL import Image
import subprocess
import re

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def main():
    rng=np.random.default_rng(732);checks=[]
    coefficients=pd.read_csv(OUT/'association_coefficients.tsv',sep='\t')
    source_hashes=[]
    oldhash={}
    for directory in ['EC_gene_modules','xenium_nuclear_anchored_segmentation','endothelial_QC_subtypes','whole_heart_spatial_atlas','fibroblast_spatial_consensus','endothelial_five_types','ventricular_spatial_resolution','fibroblast_Xenium_reviewed_refit']:
        f=ROOT/'results'/directory/'file_manifest.tsv'
        if not f.exists():continue
        t=pd.read_csv(f,sep='\t')
        col='path' if 'path' in t else 'file'
        for _,r in t.iterrows():
            p=str(r[col]);p=p if p.startswith(('analysis/','results/')) else 'results/'+directory+'/'+p
            oldhash[p]=r.sha256
    paths=json.loads((OUT/'source_files.json').read_text())
    for rel in paths:
        h=digest(ROOT/rel)
        if rel in oldhash:assert h==oldhash[rel],rel
        source_hashes.append({'path':rel,'sha256':h,'prior_manifest_match':h==oldhash[rel] if rel in oldhash else 'not_in_prior_manifest'})
    for s in STAGES:
        print('VALIDATE',s,flush=True)
        ec=read(f'results/cardiac_spatial_maturation_EC/{s}/EC_activity_annotations.tsv.gz')
        a=ad.read_h5ad(ROOT/f'results/EC_gene_modules/{s}/Xenium_ECs_with_gene_modules.h5ad')
        np.testing.assert_array_equal(ec.index,a.obs_names)
        np.testing.assert_allclose(ec[['x_centroid','y_centroid']],a.obs[['x_centroid','y_centroid']],atol=1e-5)
        genes=a.var_names.to_numpy(dtype=str);gi=np.isin(genes,CYCLING)
        x=a.layers['nuclear_only'];det=np.asarray((x[:,gi]>0).sum(1)).ravel()
        np.testing.assert_array_equal(ec.measured_nuclear_cycling_genes,det)
        np.testing.assert_array_equal(ec.measured_nuclear_cycling_positive,det>=2)
        z=np.log1p(x[:,gi].toarray()*1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1)[:,None]).mean(1)
        np.testing.assert_allclose(ec.measured_nuclear_cycling,z,rtol=2e-6,atol=1e-7)
        with np.load(OUT/s/'CM_panel_withheld_EC_reference_weights.npz') as b:
            assert not set(b['fitting_genes'])&set(MAT.gene.tolist()+['Angpt2']+CYCLING)
            np.testing.assert_array_equal(b['query_ids'],ec.index)
            predicted=np.einsum('nk,nkg->ng',b['neighbor_weights'],b['reference_target_expression'][b['neighbor_indices']])
            target=b['reference_gene_targets'].astype(str)
            for col,panel in [('inferred_sprouting_no_CM_panel','Apln Esm1 Angpt2 Cxcr4 Pgf Adm'.split()),
                              ('inferred_missing5_no_CM_panel','Apln Esm1 Cxcr4 Pgf Adm'.split()),
                              ('inferred_proliferation_no_CM_panel','Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split())]:
                np.testing.assert_allclose(predicted[:,np.isin(target,panel)].mean(1),ec[col],rtol=2e-6,atol=1e-7)
        cm=read(f'results/cardiac_spatial_maturation_EC/{s}/ventricular_CM_maturation.tsv.gz')
        previous=read(f'results/fibroblast_Xenium_reviewed_refit/{s}/CM_measured_outcomes.tsv.gz')
        np.testing.assert_array_equal(cm.index,previous.index)
        np.testing.assert_allclose(cm.maturation_21gene_mean,previous.maturation_21gene_mean,rtol=2e-8)
        assert cm.stage_cell_id.is_unique and cm.stage_cell_id.eq(s+'::'+cm.index.astype(str)).all()
        chosen=rng.choice(len(cm),min(64,len(cm)),replace=False)
        controls=pd.read_csv(ROOT/'results/fibroblast_Xenium_reviewed_refit/CM_maturation_fixed_controls.tsv',sep='\t')
        b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{s}/nuclear_anchored_counts.h5ad',backed='r')
        ci=b.obs_names.get_indexer(cm.index[chosen]);bg=b.var_names.to_numpy(dtype=str)
        for view,matrix in [('refined',b.X),('core',b.layers['high_confidence']),('nuclear',b.layers['nuclear_only'])]:
            counts=matrix[ci].toarray();logs=np.log1p(counts*1e4/np.maximum(counts.sum(1),1)[:,None])
            outcome=logs[:,np.isin(bg,MAT.gene)].mean(1)
            matched=outcome-logs[:,np.isin(bg,controls.loc[controls.matched_control,'gene'])].mean(1)
            np.testing.assert_allclose(cm['maturation_'+view].iloc[chosen],outcome,rtol=2e-6,atol=1e-7)
            np.testing.assert_allclose(cm['maturation_'+view+'_matched_control'].iloc[chosen],matched,rtol=2e-6,atol=1e-6)
        b.file.close()
        exy=ec[['x_centroid','y_centroid']].to_numpy();cxy=cm[['x_centroid','y_centroid']].iloc[chosen].to_numpy()
        distance=np.linalg.norm(cxy[:,None,:]-exy[None,:,:],axis=2)
        for radius in [30,50,100]:
            d=read(f'results/cardiac_spatial_maturation_EC/{s}/CM_EC_neighborhoods_{radius}um.tsv.gz')
            np.testing.assert_array_equal(d.index,cm.index)
            m=distance<=radius;n=m.sum(1)
            np.testing.assert_array_equal(n,d.eligible_EC_neighbors.iloc[chosen])
            for key in EXPOSURES[:-1]+['inferred_sprouting_no_CM_panel','inferred_proliferation_no_CM_panel']:
                means=np.divide(m@ec[key].to_numpy(float),n,out=np.full(len(n),np.nan),where=n>0)
                np.testing.assert_allclose(means,d['local_'+key].iloc[chosen],rtol=3e-6,atol=1e-7,equal_nan=True)
            np.testing.assert_allclose(distance[:,det>=2].min(1),d.nearest_cycling_EC_um.iloc[chosen],rtol=2e-6,atol=1e-4)
            np.testing.assert_array_equal(d.model_eligible,(d.eligible_EC_neighbors>=3)&(d.resolved_neighbors>=5)&d.arterial_distance_um.notna())
        restricted=read(f'results/cardiac_spatial_maturation_EC/{s}/CM_low_CM_signal_EC_neighborhoods_50um.tsv.gz')
        sample=restricted.iloc[rng.choice(len(restricted),min(64,len(restricted)),replace=False)]
        clean=ec.loc[ec.nuclear_CM_fraction.le(.05)]
        distances=np.linalg.norm(sample[['x_centroid','y_centroid']].to_numpy()[:,None,:]-clean[['x_centroid','y_centroid']].to_numpy()[None,:,:],axis=2)
        within=distances<=50;n=within.sum(1)
        np.testing.assert_array_equal(n,sample.eligible_EC_neighbors)
        for key in ['measured_nuclear_cycling','inferred_sprouting','inferred_proliferation']:
            np.testing.assert_allclose(within@clean[key].to_numpy(float)/n,sample['local_'+key],rtol=3e-6,atol=1e-7)
        positive=clean.measured_nuclear_cycling_positive.eq(1).to_numpy()
        np.testing.assert_array_equal(distances[:,positive].min(1)<=30,sample.near_cycling_EC_30um.eq(1))
        # Independently recover primary raw-unit coefficients using partial regression.
        d=read(f'results/cardiac_spatial_maturation_EC/{s}/CM_EC_neighborhoods_50um.tsv.gz')
        d=d.loc[d.model_eligible];base=design(d).to_numpy();base=np.column_stack([np.ones(len(d)),base])
        y=d.maturation_21gene_mean.to_numpy();yr=y-base@np.linalg.lstsq(base,y,rcond=1e-9)[0]
        subset=coefficients.loc[coefficients.stage.eq(s)&coefficients.radius_um.eq(50)&coefficients.model.eq('full')&coefficients.outcome.eq('maturation_21gene_mean')&coefficients.block_um.eq(400)&coefficients.exposure.isin(EXPOSURES)]
        assert len(subset)==len(EXPOSURES)
        for _,r in subset.iterrows():
            key=r.exposure if r.exposure.startswith('near_') else 'local_'+r.exposure
            e=d[key].to_numpy();er=e-base@np.linalg.lstsq(base,e,rcond=1e-9)[0]
            partial=(er@yr)/(er@er)
            np.testing.assert_allclose(partial,r.raw_unit_coefficient,rtol=3e-5,atol=1e-7)
        checks.append({'stage':s,'strict_ventricular_CMs':len(cm),'eligible_ECs':len(ec),'primary_50um_CMs':len(d),
                       'brute_force_CM_geometry_examples':len(chosen),'primary_coefficients_recovered':len(subset)})
    info=subprocess.check_output(['pdfinfo',str(OUT/'cardiac_spatial_maturation_EC_report.pdf')],text=True)
    pages=int(re.search(r'^Pages:\s+(\d+)',info,re.MULTILINE).group(1));assert pages==10
    for f in OUT.glob('*.png'):
        with Image.open(f) as im:im.verify()
    save(pd.DataFrame(source_hashes),'source_checksums.tsv')
    dump(OUT/'validation.json',{'status':'PASS','stages':checks,'PDF_pages':pages,
         'all_21_CM_targets_and_all_sprouting_targets_withheld_from_new_EC_predictors':True,
         'source_counts_coordinates_and_identity_files_match_prior_hashes_where_available':True,
         'source_CM_maturation_scores_preserved':True,'nuclear_cycling_marker_support_recalculated':True,
         'alternative_CM_boundary_scores_and_restricted_EC_neighborhoods_independently_recalculated':True,
         'heldout_module_predictions_reconstructed_from_reference_weights':True,
         'all_primary_exposure_coefficients_recovered_by_partial_regression':True,
         'scope':'Arithmetic and provenance validation, not independent biological validation'})
    rows=[]
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and p.name!='file_manifest.tsv':rows.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)})
    for name in ['cardiac_spatial_maturation_EC.py','report_cardiac_spatial_maturation_EC.py','validate_cardiac_spatial_maturation_EC.py']:
        p=ROOT/'analysis'/name;rows.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)})
    save(pd.DataFrame(rows),'file_manifest.tsv')
    print('PASS',checks,flush=True)

if __name__=='__main__':main()
