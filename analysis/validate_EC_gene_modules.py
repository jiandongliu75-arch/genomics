"""Check module arithmetic, preserved observed counts and honest score labels."""
from evaluate_EC_gene_modules import *
from PIL import Image
import subprocess,re

def exact_sparse(a,b):
    assert a.shape==b.shape
    diff=(a-b).tocsr();diff.eliminate_zeros();assert diff.nnz==0

def main():
    hashes=[];n=0;keys=[]
    baseline=pd.read_csv(PROJ/'file_manifest.tsv',sep='\t').set_index('path')
    for stage in STAGES:
        path=PROJ/stage/'Xenium_ECs_with_inferred_sprouting.h5ad'
        source=ad.read_h5ad(path)
        a=ad.read_h5ad(MODULE_OUT/stage/'Xenium_ECs_with_gene_modules.h5ad')
        assert a.n_vars==source.n_vars==403
        np.testing.assert_array_equal(a.obs_names,source.obs_names);np.testing.assert_array_equal(a.var_names,source.var_names)
        exact_sparse(a.X,source.X)
        for name in source.layers:exact_sparse(a.layers[name],source.layers[name])
        np.testing.assert_array_equal(a.obsm['spatial'],source.obsm['spatial'])
        np.testing.assert_array_equal(a.obs.atlas_identity,source.obs.atlas_identity)
        np.testing.assert_allclose(a.obs.inferred_sprouting_module,source.obs.inferred_sprouting_score_6g,atol=1e-6,rtol=1e-5)
        np.testing.assert_allclose(a.obs.inferred_unmeasured_sprouting_module,source.obs.inferred_unmeasured_sprouting_score_5g,atol=1e-6,rtol=1e-5)
        module_fields=[c for c in a.obs if c.startswith('inferred_module__')]
        assert len(module_fields)==9 and np.isfinite(a.obs[module_fields].to_numpy()).all()
        assert a.obs.sprouting_annotation_status.str.startswith('Exploratory').all()
        assert a.uns['gene_module_scores']['status'].startswith('Exploratory')
        rna=pd.read_csv(MODULE_OUT/stage/'RNA_reference_gene_module_scores.tsv.gz',sep='\t',index_col=0)
        with np.load(PROJ/stage/'nuclear_RNA_reference_neighbors.npz') as b:
            np.testing.assert_array_equal(b['query_cell_ids'].astype(str),a.obs_names)
            pos=rna.index.get_indexer(b['reference_cell_ids'].astype(str));assert (pos>=0).all()
            values=rna.iloc[pos].to_numpy();projected=np.einsum('nk,nkg->ng',b['neighbor_weights'],values[b['neighbor_indices']])
            np.testing.assert_allclose(projected,a.obs[module_fields].to_numpy(),rtol=1e-6,atol=1e-7)
        h=sha(path);assert h==baseline.loc[str(path.relative_to(ROOT)),'sha256']
        hashes.append({'path':str(path.relative_to(ROOT)),'sha256':h})
        keys.extend(a.obs.stage_cell_id);n+=len(a)
    assert n==len(keys)==len(set(keys))==29363
    allscores=pd.read_csv(MODULE_OUT/'all_stage_gene_module_annotations.tsv.gz',sep='\t')
    assert len(allscores)==29363 and allscores.stage_cell_id.is_unique
    equiv=pd.read_csv(MODULE_OUT/'module_aggregation_order_equivalence.tsv',sep='\t')
    assert len(equiv)==4 and equiv.module_first_vs_gene_first_max_abs_error.lt(1e-6).all()
    bench=pd.read_csv(MODULE_OUT/'module_definition_opposite_library_benchmarks.tsv',sep='\t')
    assert len(bench)==432 and set(bench.score_definition)==set(DEFINITIONS)
    assert set(bench.module)==set(MODULES) and set(bench.RNA_resolution_nuclei)=={1,5,20}
    assert bench.Spearman.dropna().between(-1,1).all()
    summary=pd.read_csv(MODULE_OUT/'module_definition_benchmark_summary.tsv',sep='\t')
    for row in summary.itertuples():
        selected=bench[(bench.stage==row.stage)&(bench.RNA_resolution_nuclei==row.RNA_resolution_nuclei)&(bench.alignment==row.alignment)&(bench.module==row.module)&(bench.score_definition==row.score_definition)]
        assert len(selected)==2 and bool(selected.technical_gate.all())==row.both_libraries_gate
    info=subprocess.check_output(['pdfinfo',str(MODULE_OUT/'EC_gene_module_comparison.pdf')],text=True)
    pages=int(re.search(r'^Pages:\s+(\d+)',info,re.MULTILINE).group(1));assert pages==5
    for p in MODULE_OUT.glob('*.png'):
        with Image.open(p) as im:im.verify()
    write_json(MODULE_OUT/'validation.json',{'status':'PASS','ECs':29363,'measured_features':403,'source_count_matrices_and_coordinates_preserved':True,
        'all_9_module_predictions_reconstructed_from_fixed_reference_weights':True,'mean_module_identical_to_previous_six_gene_score':True,
        'aggregation_order_changes_only_floating_point_rounding':True,'RNA_library_benchmark_cases':432,
        'inferred_module_definitions_not_independently_validated_activity':True,'PDF_pages':pages})
    save(pd.DataFrame(hashes),MODULE_OUT/'source_counts_sha256.tsv')
    files=[p for p in MODULE_OUT.rglob('*') if p.is_file() and p.name!='file_manifest.tsv']+[ROOT/'analysis/evaluate_EC_gene_modules.py',ROOT/'analysis/validate_EC_gene_modules.py']
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]),MODULE_OUT/'file_manifest.tsv')
    print('EC_MODULE_VALIDATION_PASS',flush=True)

if __name__=='__main__':main()
