"""Retrieve public WebAtlas annotations and verify original Xenium correspondence.

Only small annotation/coordinate arrays are fetched. Missing S3 objects may
return403; that is not treated as evidence of a network-policy block.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import itertools
import json
import time

import anndata as ad
import numcodecs
import numpy as np
import pandas as pd
import requests
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT / 'results/xenium_ventricular_mapping'
OUT = ROOT / 'results/author_atlas_baseline'
CACHE = Path('/workspace/heart-atlas/data/author_webatlas')
HOST = 'https://heart-atlas.s3.us-east-2.amazonaws.com'
STAGES = ['P0', 'P7', 'P14', 'P21']


class RemoteArrays:
    def __init__(self, stage):
        self.base = f'{HOST}/{stage.lower()}/heart-xenium-anndata.zarr/'
        self.cache = CACHE / stage
        self.records = []

    def fetch(self, key):
        local = self.cache / key
        local.parent.mkdir(parents=True, exist_ok=True)
        if not local.exists():
            r = requests.get(self.base + key, timeout=40)
            r.raise_for_status()
            local.write_bytes(r.content)
            header = {name: r.headers.get(name) for name in ['ETag', 'Last-Modified', 'Content-Type']}
            Path(str(local) + '.http.json').write_text(json.dumps(header))
        content = local.read_bytes()
        self.records.append({'key': key, 'url': self.base + key, 'bytes': len(content),
            'sha256': hashlib.sha256(content).hexdigest(),
            'response_headers': json.loads(Path(str(local) + '.http.json').read_text())})
        return content

    def array(self, path):
        m = json.loads(self.fetch(path + '/.zarray'))
        shape, chunks = m['shape'], m['chunks']
        output = np.empty(shape, dtype=np.dtype(m['dtype']))
        counts = [(s + c - 1) // c for s, c in zip(shape, chunks)]
        for index in itertools.product(*[range(n) for n in counts]):
            key = m.get('dimension_separator', '.').join(map(str, index))
            value = self.fetch(path + '/' + key)
            if m['compressor']:
                value = numcodecs.get_codec(m['compressor']).decode(value)
            for codec in reversed(m['filters'] or []):
                value = numcodecs.get_codec(codec).decode(value)
            values = value if isinstance(value, np.ndarray) else np.frombuffer(value, dtype=m['dtype'])
            values = values.reshape(chunks, order=m['order'])
            selected = tuple(slice(i*c, min((i+1)*c, s)) for i, c, s in zip(index, chunks, shape))
            crop = tuple(slice(0, sl.stop-sl.start) for sl in selected)
            output[selected] = values[crop]
        return output

    def column(self, name, axis='obs'):
        path = axis + '/' + name
        attrs = json.loads(self.fetch(path + '/.zattrs'))
        if attrs.get('encoding-type') == 'categorical':
            codes, categories = self.array(path + '/codes'), self.array(path + '/categories')
            assert ((codes >= -1) & (codes < len(categories))).all()
            return pd.Categorical.from_codes(codes, categories, ordered=attrs.get('ordered', False))
        return self.array(path)


def stage_acquisition(stage):
    reader = RemoteArrays(stage)
    folder = OUT / stage
    folder.mkdir(parents=True, exist_ok=True)
    attrs = json.loads(reader.fetch('obs/.zattrs'))
    names = reader.array('obs/' + attrs['_index'])
    desired = ['cell_id', 'x_centroid', 'y_centroid', 'transcript_counts', 'total_counts',
        'control_probe_counts', 'control_codeword_counts', 'cell_area', 'nucleus_area', 'cytoSPACE.final.anno']
    columns = [name for name in desired if name in attrs['column-order']]
    author = pd.DataFrame({name: reader.column(name) for name in columns}, index=names)
    author.index.name = 'author_row_id'
    assert len(author) and author.index.is_unique
    physical = reader.array('obsm/X_spatial')
    assert physical.shape == (len(author), 2)
    author['atlas_X_spatial_0'], author['atlas_X_spatial_1'] = physical.T
    var_attrs = json.loads(reader.fetch('var/.zattrs'))
    features = reader.array('var/' + var_attrs['_index']).astype(str)
    pd.DataFrame({'atlas_feature': features, 'TF_prefix': np.char.startswith(features, 'TF_')}).to_csv(folder / 'atlas_features.tsv', sep='\t', index=False)
    data = ad.read_h5ad(PREVIOUS / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad', backed='r')
    original = data.obs.copy()
    # Locate cells by explicit barcodes when available; numeric row indices are
    # not assumed to equal cell identities or data order.
    ids = author.cell_id.astype(str).to_numpy()
    assert len(set(ids)) == len(ids)
    # Object-array np.isin can be quadratic for hundreds of thousands of IDs.
    exact = original.index.get_indexer(ids) >= 0
    assert exact.all(), f'{stage}: {exact.sum()}/{len(ids)} explicit cell IDs match'
    source = original.loc[ids]
    np.testing.assert_allclose(author[['x_centroid', 'y_centroid']].to_numpy(), source[['x_centroid', 'y_centroid']].to_numpy(), atol=1e-7, rtol=0)
    delta = np.abs(author[['x_centroid', 'y_centroid']].to_numpy() - source[['x_centroid', 'y_centroid']].to_numpy())
    for column in ['transcript_counts', 'total_counts', 'control_probe_counts', 'control_codeword_counts']:
        assert np.array_equal(author[column].to_numpy(), source[column].to_numpy()), column
    # Keep exported display coordinates and observed µm centroids separate.
    author['stage'] = stage
    author.to_csv(folder / 'author_atlas_annotations.tsv.gz', sep='\t', compression={'method': 'gzip', 'compresslevel': 3})
    joined = original.copy()
    mapper = pd.Series(author['cytoSPACE.final.anno'].astype(str).to_numpy(), index=ids)
    joined['author_atlas_cell_state'] = mapper.reindex(joined.index).fillna('Not present in published atlas')
    joined['present_in_author_atlas'] = joined.index.isin(ids)
    joined['author_row_id'] = pd.Series(author.index, index=ids).reindex(joined.index)
    joined.to_csv(folder / 'original_cells_with_author_labels.tsv.gz', sep='\t', compression={'method': 'gzip', 'compresslevel': 3})
    summary = {'stage': stage, 'atlas_url': reader.base, 'retrieval_date_UTC': time.strftime('%Y-%m-%d'),
        'raw_Xenium_segments': len(original), 'author_atlas_segments': len(author), 'all_author_barcodes_match_raw': bool(exact.all()),
        'max_observed_centroid_difference_um': float(delta.max()), 'measured_transcript_counts_match_raw': True,
        'author_labels': author['cytoSPACE.final.anno'].value_counts().to_dict(),
        'atlas_features': len(features), 'measured_Xenium_targets': data.n_vars,
        'atlas_features_matching_measured_panel': int(np.isin(features, data.var_names).sum()),
        'Kcnd2_in_atlas_features': bool('Kcnd2' in features), 'Irx5_in_atlas_features': bool('Irx5' in features),
        'TF_Irx5_in_atlas_features': bool('TF_Irx5' in features),
        'OriginalCID_snRNA_present_in_exported_obs': 'OriginalCID_snRNA' in attrs['column-order'],
        'exported_obs_columns': attrs['column-order'],
        'interpretation': 'Retrieved final published atlas labels, not verified original CytoSPACE input labels. No original RNA reference barcode mapping is present in exported obs.',
        'spatial_array_vs_centroid': {'max_abs_same_axes': float(np.abs(physical - author[['x_centroid', 'y_centroid']].to_numpy()).max())}}
    (folder / 'acquisition_summary.json').write_text(json.dumps(summary, indent=2))
    (folder / 'http_provenance.json').write_text(json.dumps(reader.records, indent=2))
    data.file.close()
    print(f'ACQUIRED {stage}: {len(author)} author labels; {len(original)} original segments; max centroid delta {delta.max():.3g}µm', flush=True)
    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        summaries = list(pool.map(stage_acquisition, STAGES))
    (OUT / 'author_atlas_acquisition.json').write_text(json.dumps(summaries, indent=2))
    counts = [{'stage': x['stage'], 'author_cell_state': label, 'segments': n} for x in summaries for label, n in x['author_labels'].items()]
    pd.DataFrame(counts).to_csv(OUT / 'author_cell_state_counts.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
