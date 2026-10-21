"""Project QC-reviewed snRNA reference onto exact shared Xenium gene symbols."""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/xenium_ventricular_mapping'
CACHE = BASE / 'results/xenium_mapping_cache'
REFERENCE = BASE / 'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'


def main():
    OUT.mkdir(exist_ok=True)
    CACHE.mkdir(exist_ok=True)
    folders = sorted((BASE / 'data/xenium').iterdir())
    panel = None
    for folder in folders:
        with h5py.File(folder / 'cell_feature_matrix.h5') as handle:
            features = handle['matrix/features']
            names = np.array([v.decode() for v in features['name'][:]])
            types = np.array([v.decode() for v in features['feature_type'][:]])
            current = names[types == 'Gene Expression']
            if panel is None:
                panel = current
            else:
                assert np.array_equal(panel, current)
    with h5py.File(REFERENCE) as handle:
        obs, var = read_elem(handle['obs']), read_elem(handle['var'])
        symbols = var.gene_symbol.astype(str).to_numpy()
        print('READ_REFERENCE_CORRECTED', flush=True)
        matrix = sparse_dataset(handle['X'])[:].tocsr().astype(np.float32)
        abundance = np.asarray(matrix.sum(axis=0)).ravel()
        matches, shared, mapping = [], [], []
        for gene in panel:
            candidates = np.flatnonzero(symbols == gene)
            feature = int(candidates[np.argmax(abundance[candidates])]) if len(candidates) else None
            mapping.append({'xenium_target': gene, 'exact_symbol_in_reference': feature is not None, 'reference_gene_id': var.index[feature] if feature is not None else 'NA', 'mapping_note': 'Exact symbol; greatest pooled abundance if duplicate features' if feature is not None else 'Excluded from transfer: no exact gene-symbol match; custom isoform targets not silently collapsed'})
            if feature is not None:
                matches.append(feature)
                shared.append(gene)
        corrected = matrix[:, matches].toarray()
        del matrix
        print('READ_REFERENCE_ORIGINAL', flush=True)
        matrix = sparse_dataset(handle['layers/counts_raw'])[:].tocsr()
        original = matrix[:, matches].toarray()
        del matrix
    assert len(panel) == 403 and 'Kcnd2' not in panel and 'Irx5' not in panel
    catalog = pd.read_csv(ROOT / 'results/ventricular_stage_subtype_catalog/stage_subtype_catalog.tsv', sep='\t')
    vent = pd.read_csv(ROOT / 'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz', index_col=0)
    ref = obs[['age', 'sample', 'reviewed_identity', 'proposed_subtype']].copy()
    ref['broad_class'] = ref.reviewed_identity.astype(str)
    ref.loc[ref.reviewed_identity.eq('Cardiomyocyte'), 'broad_class'] = 'Atrial-like CM'
    ref.loc[vent.index, 'broad_class'] = 'Ventricular-enriched CM'
    ref['ventricular_group'] = pd.NA
    ref.loc[vent.index, 'ventricular_group'] = vent.microenvironment_cluster.to_numpy()
    ref['mapping_family'] = pd.NA
    families = {
        'Mechanical-response-associated': 'Ankrd1-associated',
        'Conduction-program-enriched': 'Conduction-program candidate',
        'Regional electrical-associated': 'Regional electrical candidate',
        'Cycling': 'Cycling-associated',
        'DNA-replication/maintenance-associated': 'DNA-maintenance-associated',
        'Unassigned': 'Uncharacterized ventricular candidate',
    }
    lookup = catalog.set_index('microenvironment_cluster')
    familymap = {cluster: families.get(row.functional_family, 'Other working CM') for cluster, row in lookup.iterrows()}
    ref.loc[vent.index, 'mapping_family'] = vent.microenvironment_cluster.map(familymap).to_numpy()
    ref.to_csv(CACHE / 'reference_metadata.tsv.gz', sep='\t', na_rep='NA')
    np.savez_compressed(CACHE / 'reference_shared_panel_counts.npz', corrected=corrected, original=original, genes=np.array(shared), cell_ids=ref.index.to_numpy(dtype=str))
    pd.DataFrame(mapping).to_csv(OUT / 'panel_feature_mapping.tsv', sep='\t', index=False)
    catalog['mapping_family'] = catalog.microenvironment_cluster.map(familymap)
    catalog.to_csv(OUT / 'reference_group_catalog.tsv', sep='\t', index=False, na_rep='NA')
    ref.groupby(['age', 'broad_class'], observed=True).size().rename('reference_nuclei').to_csv(OUT / 'reference_lineage_counts.tsv', sep='\t')
    required = set('Kcnd2 Irx5 Irx3 Kcnip2 Mpped2 Etv1 Gja1 Gja5 Scn5a Scn10a Cntn2 Pcp4 Slit2 Ankrd1 Ddc Fgf14 Pdzrn4 Lrtm1 Stk39 Zbtb16 Ncam1 Pola1 Rad51b Mki67 Top2a Ppara Rreb1 Ivns1abp Myom2 Cd36 Hmgcs2 Pdk4 Tnnt2 Tnni3 Myl2 Myl4 Myl7 Hey2 Npr3 Upk3b Wt1'.split())
    for value in catalog.interpretation_markers.dropna():
        required.update(value.split(';'))
    pd.DataFrame([{'gene': gene, 'measured_by_Xenium': gene in panel, 'used_as_exact_shared_feature': gene in shared} for gene in sorted(required)]).to_csv(OUT / 'subtype_marker_panel_coverage.tsv', sep='\t', index=False)
    programs = json.loads((ROOT / 'results/ventricular_microenvironment/program_definitions.json').read_text())['programs']
    coverage = []
    for key, value in programs.items():
        present = [gene for gene in value['genes'] if gene in panel]
        absent = [gene for gene in value['genes'] if gene not in panel]
        coverage.append({'program': key, 'signature_genes': len(value['genes']), 'measured_genes': len(present), 'measured_symbols': ';'.join(present) or 'NA', 'unmeasured_symbols': ';'.join(absent) or 'NA'})
    pd.DataFrame(coverage).to_csv(OUT / 'program_panel_coverage.tsv', sep='\t', index=False)
    provenance = {'reference': str(REFERENCE), 'reference_nuclei': len(ref), 'ventricular_reference_nuclei': len(vent), 'measured_Xenium_targets': len(panel), 'exact_shared_features': len(shared), 'unmeasured_requested_marker': 'Kcnd2', 'unmeasured_regulator': 'Irx5', 'panel_identity': 'FPUV2R / mMulti_403g', 'normalization_for_transfer': 'Normalize corrected and original reference counts independently over EXACT shared panel genes; log1p(10000 * count / shared-panel total). Xenium is normalized over the same shared targets.', 'custom_targets': 'Custom splice/exon targets without exact symbol matches are retained as measured Xenium features, but excluded from reference transfer. No gene/isoform equivalence inferred.', 'source_records': ['https://zenodo.org/records/14961929', 'https://zenodo.org/records/14961942'], 'source_license': 'Zenodo CC BY 4.0; source authors Wang, Dong, Song, Qian and Liu. snRNA reference GSE290989; DOI10.1038/s44161-026-00861-z.', 'cache_npz_sha256': hashlib.sha256((CACHE / 'reference_shared_panel_counts.npz').read_bytes()).hexdigest()}
    (OUT / 'reference_preparation.json').write_text(json.dumps(provenance, indent=2))
    print(f'PREPARED {len(ref)} nuclei / {len(shared)} exact shared features. Kcnd2 and Irx5 not assayed.', flush=True)


if __name__ == '__main__':
    main()
