"""Extract explicitly reference-derived markers, including unassayed Kcnd2.

These expression values belong to snRNA nuclei, never to measured Xenium RNA.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

from run_cytospace_xenium import SOURCE, OUT, BASE, digest

GENES = 'Kcnd2 Irx5 Irx3 Etv1 Mpped2 Fgf14 Pdzrn4 Kcnip2 Kcnd3 Gja1 Scn5a Hey2 Slit2 Ankrd1 Cntn2 Gja5 Ddc Myl2'.split()


def main():
    path = BASE / 'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'
    meta = pd.read_csv(SOURCE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    output = meta[['age', 'sample', 'broad_class', 'ventricular_group', 'mapping_family']].copy()
    mapping = []
    with h5py.File(path) as handle:
        obs, var = read_elem(handle['obs']), read_elem(handle['var'])
        assert np.array_equal(obs.index, meta.index)
        symbols = var.gene_symbol.astype(str).to_numpy()
        corrected = sparse_dataset(handle['X'])[:].tocsr()
        abundance = np.asarray(corrected.sum(axis=0)).ravel()
        selected = []
        for gene in GENES:
            candidates = np.flatnonzero(symbols == gene)
            assert len(candidates)
            index = int(candidates[np.argmax(abundance[candidates])])
            selected.append(index)
            mapping.append({'gene': gene, 'reference_gene_id': var.index[index], 'duplicate_features': len(candidates),
                'selection': 'highest pooled corrected abundance among exact-symbol matches'})
        totals = np.asarray(corrected.sum(axis=1)).ravel()
        values = corrected[:, selected].toarray()
        for j, gene in enumerate(GENES):
            output[f'{gene}_corrected_count'] = values[:, j]
            output[f'{gene}_corrected_log1pCP10k'] = np.log1p(values[:, j] / totals * 10000)
        del corrected
        raw = sparse_dataset(handle['layers/counts_raw'])[:].tocsr()
        totals = np.asarray(raw.sum(axis=1)).ravel()
        values = raw[:, selected].toarray()
        for j, gene in enumerate(GENES):
            output[f'{gene}_original_count'] = values[:, j]
            output[f'{gene}_original_log1pCP10k'] = np.log1p(values[:, j] / totals * 10000)
    output.to_csv(OUT / 'snRNA_reference_marker_values_NOT_Xenium_measurements.tsv.gz', sep='\t', na_rep='NA')
    pd.DataFrame(mapping).to_csv(OUT / 'reference_marker_feature_mapping.tsv', sep='\t', index=False)
    (OUT / 'reference_marker_provenance.json').write_text(json.dumps({'source': str(path), 'source_sha256': digest(path),
        'normalization': 'log1p(10000 * corrected/original gene counts / whole-transcriptome total), independently for each count layer',
        'interpretation': 'Measured snRNA-reference values; assigning these to Xenium coordinates produces reference-derived predictions only, never measurements',
        'Kcnd2_Irx5': 'Neither gene is assayed by these Xenium panels; their predictions have no direct held-out validation'}, indent=2))


if __name__ == '__main__':
    main()
