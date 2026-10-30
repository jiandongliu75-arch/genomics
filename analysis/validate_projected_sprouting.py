"""Audit observed counts, source weights, atlas coverage and inference labels."""
from project_RNA_sprouting_to_Xenium import *
from PIL import Image
import subprocess,re

def sparse_equal(a,b):
    assert a.shape==b.shape
    diff=(a-b).tocsr();diff.eliminate_zeros();assert diff.nnz==0

def main():
    coverage=[];hashes=[];all_stage_keys=[]
    baseline=pd.read_csv(OUT/'file_manifest.tsv',sep='\t').set_index('path')
    for stage in STAGES:
        folder=PROJ/stage
        srcpath=OUT/'Xenium'/stage/'EC_candidate_count_views.h5ad'
        rnapath=OUT/'snRNA'/stage/'EC_corrected_and_raw_counts.h5ad'
        source=ad.read_h5ad(srcpath);rna=ad.read_h5ad(rnapath)
        data=ad.read_h5ad(folder/'Xenium_ECs_with_inferred_sprouting.h5ad')
        selected=(source.obs.native_EC_broad&source.obs.nuclear_broad_EC_clean&~source.obs.multiple_projected_nuclear_components_review&~source.obs.geometry_unresolved_review).to_numpy()
        np.testing.assert_array_equal(data.obs_names,source.obs_names[selected])
        np.testing.assert_array_equal(data.var_names,source.var_names);assert data.n_vars==403
        assert not set(UNMEASURED)&set(data.var_names) and 'Angpt2' in data.var_names
        sparse_equal(data.X,source.layers['native'][selected])
        sparse_equal(data.layers['nuclear_only'],source.layers['nuclear_only'][selected])
        sparse_equal(data.layers['boundary_core'],source.layers['high_confidence'][selected])
        np.testing.assert_array_equal(data.obsm['spatial'],source.obs.loc[selected,['x_centroid','y_centroid']].to_numpy())
        assert data.obs.atlas_identity.to_numpy().tolist()==source.obs.loc[selected,'atlas_identity'].to_numpy().tolist()
        pred=data.obsm['inferred_sprouting_genes_RNA_log1pCP10k']
        assert list(pred.columns)==SPROUT and np.isfinite(pred.to_numpy()).all() and (pred.to_numpy()>=0).all()
        np.testing.assert_allclose(pred.mean(axis=1),data.obs.inferred_sprouting_score_6g,rtol=1e-6,atol=1e-6)
        np.testing.assert_allclose(pred[UNMEASURED].mean(axis=1),data.obs.inferred_unmeasured_sprouting_score_5g,rtol=1e-6,atol=1e-6)
        angpt2=source.var_names.get_loc('Angpt2')
        np.testing.assert_array_equal(np.asarray(data.X[:,angpt2].toarray()).ravel(),data.obs.measured_Angpt2_native_molecules)
        assert data.obs.sprouting_annotation_status.str.startswith('Exploratory').all()
        assert data.uns['sprouting_projection']['status'].startswith('Exploratory')
        # Independently reconstruct every primary gene prediction from the stored
        # RNA-source identities and weights, guarding against index misalignment.
        with np.load(folder/'nuclear_RNA_reference_neighbors.npz') as b:
            np.testing.assert_array_equal(b['query_cell_ids'].astype(str),data.obs_names)
            refs=b['reference_cell_ids'].astype(str);positions=rna.obs_names.get_indexer(refs);assert (positions>=0).all()
            targetpos=rna.var_names.get_indexer(SPROUT)
            truth=target_expression(rna.X,targetpos)[positions]
            indices=b['neighbor_indices'];weights=b['neighbor_weights']
            assert (indices>=0).all() and (indices<len(refs)).all()
            np.testing.assert_allclose(weights.sum(axis=1),1,rtol=1e-6,atol=1e-6)
            expected=np.einsum('nk,nkg->ng',weights,truth[indices])
            np.testing.assert_allclose(expected,pred.to_numpy(),rtol=1e-5,atol=1e-6)
        features=pd.read_csv(folder/'alignment_feature_genes.tsv',sep='\t')
        assert not set(features.gene)&set(UNMEASURED+CYCLING+LINEAGES['Cardiomyocyte']+SEX)
        assert features.gene.isin(source.var_names).all() and features.is_measured_sprouting_target.sum()==1
        # Reference metacells must be a disjoint partition of the correct library,
        # formed without sprouting targets. Check both raw/corrected memberships.
        predictions=pd.read_csv(folder/'opposite_library_primary_predictions.tsv.gz',sep='\t')
        for lib in sorted(rna.obs['sample'].unique()):
            expected_ids=set(rna.obs_names[rna.obs['sample'].eq(lib)])
            for layer in ['corrected','raw']:
                for size in [5,20]:
                    m=pd.read_csv(folder/f'RNA_{lib}_{layer}_metacells{size}.tsv.gz',sep='\t',index_col=0)
                    ids=[x for group in m.source_cell_ids for x in group.split(';')]
                    assert len(ids)==len(set(ids)) and set(ids)==expected_ids
                    for row in m.itertuples():assert len(row.source_cell_ids.split(';'))==row.source_nuclei
                    if layer=='corrected':
                        p=predictions[(predictions.query_library==lib)&(predictions.RNA_resolution_nuclei==size)].set_index('unit_id')
                        for row in m.itertuples():
                            pos=rna.obs_names.get_indexer(row.source_cell_ids.split(';'))
                            pooled=sparse.csr_matrix(rna.X[pos].sum(axis=0))
                            target=target_expression(pooled,rna.var_names.get_indexer(OUTCOMES))
                            for outcome,value in metrics(target).items():assert np.isclose(value[0],p.loc[row.Index,'observed_'+outcome],rtol=1e-5,atol=1e-6)
        for size,radius in [(5,75),(20,150)]:
            pools=pd.read_csv(folder/f'neighborhood{size}_inferred_sprouting.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
            assert pools.index.isin(data.obs_names).all() and pools.source_EC_profiles.eq(size).all()
            assert pools.farthest_source_distance_um.le(radius).all()
            np.testing.assert_array_equal(pools[['x_centroid','y_centroid']],data.obs.loc[pools.index,['x_centroid','y_centroid']])
            for row in pools.iloc[np.linspace(0,len(pools)-1,min(25,len(pools)),dtype=int)].itertuples():
                ids=row.source_EC_ids.split(';');assert len(ids)==len(set(ids))==size and set(ids)<=set(data.obs_names)
                xy=data.obs.loc[ids,['x_centroid','y_centroid']].to_numpy();center=np.array([row.x_centroid,row.y_centroid])
                actual=np.sqrt(((xy-center)**2).sum(axis=1)).max()
                assert actual<=radius+1e-6 and np.isclose(actual,row.farthest_source_distance_um,atol=1e-6)
        atlas=pd.read_csv(ROOT/'results/whole_heart_spatial_atlas'/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        whole=pd.read_csv(folder/'all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        np.testing.assert_array_equal(whole.index,atlas.index)
        np.testing.assert_array_equal(whole[['x_centroid','y_centroid']],atlas[['x_centroid','y_centroid']])
        assert whole.atlas_identity.to_numpy().tolist()==atlas.atlas_identity.to_numpy().tolist()
        assert int(whole.sprouting_inference_eligible.sum())==len(data)
        assert whole.loc[~whole.sprouting_inference_eligible,'inferred_sprouting_score_6g'].isna().all()
        assert whole.loc[~whole.sprouting_inference_eligible,'inferred_unmeasured_sprouting_score_5g'].isna().all()
        assert whole.stage_cell_id.is_unique
        all_stage_keys.extend(whole.stage_cell_id)
        for path in [srcpath,rnapath]:
            h=sha(path);assert h==baseline.loc[str(path.relative_to(ROOT)),'sha256'];hashes.append({'path':str(path.relative_to(ROOT)),'sha256':h})
        bench=pd.read_csv(folder/'opposite_library_projection_benchmark.tsv',sep='\t')
        assert np.isfinite(bench.Spearman.dropna()).all()
        assert (bench.effective_neighbors<=bench.neighbors).all()
        assert bench.R2_vs_query_mean.dropna().between(-100,1).all()
        gates=pd.read_csv(folder/'projection_support_gates.tsv',sep='\t')
        for row in gates.itertuples():
            d=bench[bench.RNA_resolution_nuclei.eq(row.RNA_resolution_nuclei)&bench.method.eq('assay_scaled_soft')&bench.neighbors.eq(30)&bench.outcome.eq('unmeasured_5g')]
            assert len(d)==2 and bool(d.technical_support_gate.all())==row.unmeasured_5g_both_library_gate_pass
        coverage.append({'stage':stage,'projected_ECs':len(data),'original_heart_profiles':len(whole),
            'source_candidates':len(source),'measured_gene_features':data.n_vars})
        print('PROJECTED_SPROUTING_VALIDATED',stage,flush=True)
    assert len(all_stage_keys)==len(set(all_stage_keys))==408054
    assert sum(r['projected_ECs'] for r in coverage)==29363
    assert sum(r['source_candidates'] for r in coverage)==95316
    whole=pd.read_csv(PROJ/'all_cardiac_cells_with_inferred_sprouting_layer.tsv.gz',sep='\t')
    assert len(whole)==408054 and whole.stage_cell_id.is_unique
    pdf=PROJ/'RNA_inferred_sprouting_Xenium_report.pdf'
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True);pages=int(re.search(r'^Pages:\s+(\d+)',info,re.MULTILINE).group(1));assert pages==11
    for path in PROJ.glob('*.png'):
        with Image.open(path) as image:image.verify()
    assert all(p.stat().st_size<95*1024**2 for p in PROJ.rglob('*') if p.is_file())
    write_json(PROJ/'validation.json',{'status':'PASS','coverage':coverage,'all_heart_profiles':408054,'projected_ECs':29363,
        'observed_counts_layers_gene_features_and_coordinates_exactly_preserved':True,
        'every_primary_predicted_gene_reconstructed_from_saved_RNA_weights':True,
        'unmeasured_genes_absent_from_measured_var_and_predictor_features':True,
        'RNA_metacells_disjoint_and_observed_targets_reconstructed':True,
        'spatial_neighborhood_source_membership_distances_and_original_positions_checked':True,
        'excluded_profiles_have_missing_scores_not_zero':True,'all_spatial_inference_explicitly_exploratory':True,
        'no_independent_spatial_or_biological_activation_accuracy_claim':True,'PDF_pages':pages})
    save(pd.DataFrame(hashes),PROJ/'source_counts_sha256.tsv')
    files=[p for p in PROJ.rglob('*') if p.is_file() and p.name!='file_manifest.tsv']
    files += [ROOT/'analysis'/p for p in ['project_RNA_sprouting_to_Xenium.py','report_projected_sprouting.py','validate_projected_sprouting.py']]
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]),PROJ/'file_manifest.tsv')
    print('PROJECTED_SPROUTING_VALIDATION_PASS',flush=True)

if __name__=='__main__':main()
