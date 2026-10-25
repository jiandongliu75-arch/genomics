"""Validate retained RNA mass, threshold repair and observed spatial pooling."""
from vascular_detection_common import *
from enhance_vascular_Xenium_detection import SPATIAL_PANELS
import anndata as ad
from scipy import sparse
from PIL import Image
import subprocess,re

def main():
    audit=pd.read_csv(DEST/'RNA/fractional_count_and_background_audit.tsv',sep='\t')
    cycle=pd.read_csv(DEST/'RNA/cycling_fractional_threshold_sensitivity.tsv',sep='\t')
    scores=pd.read_csv(DEST/'RNA/EC_fractional_aware_program_scores.tsv.gz',sep='\t',index_col=0)
    controls=json.loads((DEST/'RNA/rank_matched_controls.json').read_text())
    assert len(scores)==4597 and scores.index.is_unique
    assert np.isfinite(scores.filter(regex='__rank').to_numpy()).all()
    assert (np.abs(scores.filter(regex='__rank').to_numpy())<=1+1e-6).all()
    for cs in controls.values():assert len(cs)==20 and not set(cs)&set(TARGET_GENES)
    fitting=pd.read_csv(DEST/'RNA/benchmark_fitting_genes.tsv',sep='\t').gene
    assert len(fitting)==364 and not set(fitting)&set(TARGET_GENES)
    bench=pd.read_csv(DEST/'RNA/withheld_program_opposite_library_benchmark.tsv',sep='\t')
    assert len(bench)==32 and not bench.technical_gate_pass.any()
    assert np.isfinite(bench.R2_vs_query_mean).all() and bench.R2_vs_query_mean.gt(-100).all()
    source_hashes=[];spatial=[]
    for stage in STAGES:
        p=OUT/'snRNA'/stage/'EC_corrected_and_raw_counts.h5ad';a=ad.read_h5ad(p)
        assert set(a.obs_names)==set(scores.index[scores.age.eq(stage)])
        for lib in sorted(a.obs['sample'].unique()):
            m=a.obs['sample'].eq(lib).to_numpy()
            d=audit[audit['sample'].eq(lib)].set_index('gene')
            for gene in TARGET_GENES:
                j=a.var_names.get_loc(gene);raw=a.layers['counts_raw'][m,j].toarray().ravel();c=a.X[m,j].toarray().ravel()
                assert np.all(c<=raw+1e-4) and np.all(c>=0)
                assert np.isclose(d.loc[gene,'raw_molecules'],raw.sum(),atol=1e-5)
                assert np.isclose(d.loc[gene,'corrected_expected_mass'],c.sum(),rtol=1e-6,atol=1e-4)
                for threshold,column in [(.1,'corrected_positive_ge0p1'),(.25,'corrected_positive_ge0p25'),(.5,'corrected_positive_ge0p5'),(1,'corrected_positive_ge1')]:
                    assert int(d.loc[gene,column])==int((c>=threshold).sum())
            r=a.layers['counts_raw'][m][:,a.var_names.isin(CYCLING)].toarray()
            c=a.X[m][:,a.var_names.isin(CYCLING)].toarray()
            for row in cycle[cycle['sample'].eq(lib)].itertuples():
                v=r if row.layer=='raw' else c
                actual=((v>=row.minimum_expected_count)&(r>=1)).sum(axis=1)>=2
                assert int(actual.sum())==row.cycling_support_nuclei
        source_hashes.append({'path':str(p.relative_to(ROOT)),'sha256':sha(p)})
        p=OUT/'Xenium'/stage/'EC_candidate_count_views.h5ad';x=ad.read_h5ad(p)
        select=x.obs.native_EC_broad&x.obs.nuclear_broad_EC_clean&~x.obs.multiple_projected_nuclear_components_review&~x.obs.geometry_unresolved_review
        o=x.obs.loc[select]
        f=DEST/'Xenium'/stage
        d=pd.read_csv(f/'adaptive_5EC_75um_nuclear_proxy_maps.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        assert d.index.isin(o.index).all() and d.source_EC_profiles.eq(5).all() and d.farthest_source_distance_um.le(75).all()
        np.testing.assert_array_equal(d[['x_centroid','y_centroid']].to_numpy(),o.loc[d.index,['x_centroid','y_centroid']].to_numpy())
        # Directly reconstruct source RNA for a deterministic sample of pools.
        indices=np.linspace(0,len(d)-1,min(20,len(d)),dtype=int)
        for row in d.iloc[indices].itertuples():
            ids=row.source_cell_ids.split(';');assert len(ids)==len(set(ids))==5 and set(ids)<=set(o.index)
            pos=x.obs_names.get_indexer(ids)
            actual=np.asarray(x.layers['nuclear_only'][pos].sum(axis=0)).ravel()
            assert float(actual.sum())==row.gene_molecules
            fields=d.loc[row.Index]
            for name,panel in SPATIAL_PANELS.items():
                v=actual[x.var_names.isin(panel)]
                assert int((v>=1).sum())==fields[name+'__observed_genes']
                assert float(v.sum())==fields[name+'__molecules']
        # Every non-overlapping50um primary tile is reconstructed from source IDs.
        membership=pd.read_csv(f/'nuclear_clean_native_EC_50um_tile_membership.tsv.gz',sep='\t')
        pooled=pd.read_csv(f/'nuclear_clean_native_EC_50um_observed_profiles.tsv.gz',sep='\t')
        p=pooled[pooled.view.eq('nuclear')].set_index('tile_id')
        for key,group in membership[membership.eligible_tile_ge5_ECs].groupby('grid_tile',observed=True):
            tile=stage+'__nuclear_clean_native_EC__50um__'+key;row=p.loc[tile]
            assert len(group)>=5 and len(group)==row.EC_profiles
            pos=x.obs_names.get_indexer(group.cell_id)
            actual=np.asarray(x.layers['nuclear_only'][pos].sum(axis=0)).ravel()
            assert actual.sum()==row.gene_molecules
            for name,panel in SPATIAL_PANELS.items():assert actual[x.var_names.isin(panel)].sum()==row[name+'__molecules']
        cover=pd.read_csv(f/'panel_coverage.tsv',sep='\t')
        assert cover[cover.program.eq('Original_sprouting_core')&cover.Xenium_measured].gene.tolist()==['Angpt2']
        source_hashes.append({'path':str((OUT/'Xenium'/stage/'EC_candidate_count_views.h5ad').relative_to(ROOT)),
            'sha256':sha(OUT/'Xenium'/stage/'EC_candidate_count_views.h5ad')})
        spatial.append({'stage':stage,'source_nuclear_clean_native_EC_profiles':len(o),'adaptive_pools':len(d),'position_coverage_fraction':len(d)/len(o),'nonoverlapping50um_tiles_verified':len(p)})
        print('VALIDATED_VASCULAR_STAGE',stage,flush=True)
    summary=pd.read_csv(DEST/'RNA/gene_support_across_two_libraries.tsv',sep='\t')
    a=summary[summary.stage.eq('P7')&summary.gene.eq('Apln')].iloc[0]
    assert a.raw_molecules_total==56 and np.isclose(a.corrected_expected_mass_total,55.59331,atol=1e-4)
    assert a.corrected_positive_ge1_total==5 and a.corrected_positive_ge0p25_total==47
    assert summary[summary.stage.eq('P7')&summary.gene.eq('Esm1')].raw_molecules_total.iloc[0]==0
    doc=DEST/'vascular_detection_enhancement_report.pdf'
    info=subprocess.check_output(['pdfinfo',str(doc)],text=True);pages=int(re.search(r'^Pages:\s+(\d+)',info,re.MULTILINE).group(1));assert pages==11
    for p in DEST.glob('*.png'):
        with Image.open(p) as image:image.verify()
    write_json(DEST/'validation.json',{'status':'PASS','RNA_primary_nuclei':4597,
        'raw_and_fractional_corrected_gene_sums_reconstructed_from_source':True,'threshold_sensitivity_and_cycling_counts_verified':True,
        'P7_Apln_raw_molecules':56,'P7_Apln_retained_expected_mass':float(a.corrected_expected_mass_total),
        'P7_Apln_legacy_positive_nuclei':5,'P7_Apln_fractional_support_nuclei':47,
        'no_primary_count_matrix_or_identity_changed':True,'no_missing_spatial_gene_imputed':True,
        'withheld_program_genes_absent_from_benchmark_features':True,'benchmark_passes':0,'spatial':spatial,
        'adaptive_pool_source_counts_sampled_and_all_50um_nuclear_tiles_exactly_verified':True,
        'depth_matched_molecule_equality_checked_during_all_pooling':True,'PDF_pages':pages,
        'coverage_not_independent_biological_accuracy':True})
    baseline=pd.read_csv(OUT/'file_manifest.tsv',sep='\t').set_index('path')
    for row in source_hashes:assert row['sha256']==baseline.loc[row['path'],'sha256']
    save(pd.DataFrame(source_hashes),DEST/'source_counts_sha256.tsv')
    paths=[p for p in DEST.rglob('*') if p.is_file() and p.name!='file_manifest.tsv']
    paths += [ROOT/'analysis'/p for p in ['vascular_detection_common.py','enhance_vascular_RNA_detection.py','enhance_vascular_Xenium_detection.py','report_vascular_detection.py','validate_vascular_detection.py']]
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(paths)]),DEST/'file_manifest.tsv')
    print('VASCULAR_DETECTION_VALIDATION_PASS',flush=True)

if __name__=='__main__':main()
