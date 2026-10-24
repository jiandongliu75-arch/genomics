"""Audit measured panels and sample-level fibroblast signals before a focused reanalysis.

No new cell labels, trajectories, maturation scores, or spatial effects are fitted.
"""
from pathlib import Path
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/fibroblast_reassessment'
OUT.mkdir(parents=True, exist_ok=True)
INPUT = BASE / 'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'
WORKBOOK = BASE / 'metadata/paper_supplementary_tables.xlsx'

# Exact labels transcribed from Figure 4b in the user-supplied manuscript.
# Categories are the figure's categories, not a revised biological classification.
MATURATION = {
    'TFs': 'Rxrg Ppargc1a Ppargc1b Cebpa Cebpb Esrra Esrrb Esrrg Mitf'.split(),
    'Sarcomere': 'Myh6 Tnni3 Myl2'.split(),
    'Electrophysiology': 'Kcnq1 Jph1 Jph2 Bin1 Lpin1'.split(),
    'Metabolism': 'Pln Ppara Ppard Cpt1b'.split(),
}
PROGRAMS = {
    'Fibrillar_collagen': 'Col1a1 Col1a2 Col3a1 Col5a1 Col5a2'.split(),
    'Fibronectin': ['Fn1'],
    'Basement_membrane': 'Lama2 Lama4 Lama5 Lamb1 Lamb2 Lamc1 Nid1 Nid2 Col4a1 Col4a2 Hspg2'.split(),
    'Collagen_VIII_XV': 'Col8a1 Col8a2 Col15a1'.split(),
    'Collagen_processing_crosslinking': 'Lox Loxl1 Loxl2 Loxl3 Loxl4 Tgm2 Plod1 Plod2 Plod3 Serpinh1'.split(),
    'Hyaluronan_synthesis_turnover': 'Has1 Has2 Has3 Hyal1 Hyal2 Cemip'.split(),
    'Matricellular': 'Postn Tnc Thbs1 Thbs2 Sparc Cthrc1'.split(),
    'Elastic_fiber_associated': 'Eln Fbn1 Fbln1 Fbln2 Fbln5 Efemp1'.split(),
    'TGFbeta_feedback_associated': 'Smad7 Skil Pmepa1 Bambi Serpine1 Junb'.split(),
    'TGFbeta_ligands': 'Tgfb1 Tgfb2 Tgfb3'.split(),
    'TGFbeta_receptors': 'Tgfbr1 Tgfbr2 Tgfbr3'.split(),
    'Cycling': 'Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split(),
    'CM_carryover_check': 'Myh6 Tnni3 Myl2 Ttn Actn2 Pln Ryr2'.split(),
}


def table3_markers():
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(WORKBOOK) as z:
        strings = [''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        rows = []
        for row in ET.fromstring(z.read('xl/worksheets/sheet4.xml')).findall('s:sheetData/s:row', ns):
            values = {}
            for c in row:
                v = c.find('s:v', ns)
                if v is not None:
                    values[''.join(filter(str.isalpha, c.attrib['r']))] = strings[int(v.text)] if c.attrib.get('t') == 's' else v.text
            if 'Fib' in values.get('A', ''):
                rows.append({'state': values['A'], 'gene': values['B'],
                             'source_avg_log2FC': float(values['D'])})
    return pd.DataFrame(rows)


def main():
    markers = table3_markers()
    mapping = pd.read_csv(ROOT / 'results/xenium_ventricular_mapping/panel_feature_mapping.tsv', sep='\t')
    spatial_genes = set(mapping.xenium_target)
    with h5py.File(INPUT) as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        genes = set(var.gene_symbol)
        ii = np.flatnonzero(obs.reviewed_identity.eq('Fibroblast'))
        fb = obs.iloc[ii].copy()
        wanted = sorted(set(markers.gene) | {g for v in PROGRAMS.values() for g in v})
        vi = np.flatnonzero(var.gene_symbol.isin(wanted))
        symbols = var.iloc[vi].gene_symbol.to_numpy()
        matrices = {}
        for label, path in [('corrected', 'X'), ('raw', 'layers/counts_raw')]:
            full = sparse_dataset(f[path])[ii].tocsr()
            expected = fb['total_counts' if label == 'corrected' else 'raw_total_counts'].to_numpy()
            np.testing.assert_allclose(np.asarray(full.sum(axis=1)).ravel(), expected, rtol=2e-5)
            selected = full[:, vi].toarray()
            matrices[label] = pd.DataFrame(selected.T, index=symbols).groupby(level=0).sum().T
            matrices[label].index = fb.index
            del full

    markers['measured_in_Xenium'] = markers.gene.isin(spatial_genes)
    markers['present_in_snRNA_features'] = markers.gene.isin(genes)
    markers.to_csv(OUT / 'paper_fibroblast_markers_and_coverage.tsv', sep='\t', index=False)
    maturation = pd.DataFrame([{'figure_category': c, 'gene': g,
        'measured_in_Xenium': g in spatial_genes, 'present_in_snRNA_features': g in genes}
        for c, gl in MATURATION.items() for g in gl])
    assert len(maturation) == 21 and maturation.gene.is_unique and maturation.measured_in_Xenium.all()
    maturation.to_csv(OUT / 'CM_maturation_Figure4b_panel.tsv', sep='\t', index=False)
    coverage, per_gene = [], []
    panels = {**{k: list(v.gene) for k, v in markers.groupby('state')}, **PROGRAMS}
    for panel, gl in panels.items():
        present = [g for g in gl if g in spatial_genes]
        coverage.append({'panel': panel, 'genes_requested': len(gl), 'Xenium_measured': len(present),
            'Xenium_present_genes': ';'.join(present) or 'None',
            'Xenium_absent_genes': ';'.join(g for g in gl if g not in spatial_genes) or 'None'})
        per_gene.extend({'panel': panel, 'gene': g, 'Xenium_measured': g in spatial_genes,
                         'snRNA_feature_present': g in genes} for g in gl)
    pd.DataFrame(coverage).to_csv(OUT / 'program_panel_coverage.tsv', sep='\t', index=False)
    pd.DataFrame(per_gene).to_csv(OUT / 'program_gene_coverage.tsv', sep='\t', index=False)

    profile_rows = []
    for subset_name, mask in [('all_reviewed_fibroblasts', np.ones(len(fb), dtype=bool)),
            ('historical_Postn_proposal', fb.proposed_subtype.eq('Postn-associated remodeling').to_numpy())]:
        for sample, sub in fb.loc[mask].groupby('sample', observed=True):
            ids = sub.index
            for layer, counts in matrices.items():
                denom = sub['total_counts' if layer == 'corrected' else 'raw_total_counts'].to_numpy(dtype=float)
                values = counts.loc[ids]
                logs = np.log1p(values.to_numpy() * 1e4 / denom[:, None])
                for j, gene in enumerate(values.columns):
                    profile_rows.append({'subset': subset_name, 'sample': sample, 'age': str(sub.age.iloc[0]),
                        'layer': layer, 'nuclei': len(sub), 'gene': gene,
                        'pseudobulk_CPM': float(values.iloc[:, j].sum() * 1e6 / denom.sum()),
                        'mean_log1p_CP10k': float(logs[:, j].mean()),
                        'fraction_count_ge1': float((values.iloc[:, j] >= 1).mean())})
    profiles = pd.DataFrame(profile_rows)
    profiles.to_csv(OUT / 'fibroblast_marker_profiles_by_library.tsv.gz', sep='\t', index=False)
    quality = fb.groupby(['age', 'sample'], observed=True).agg(
        nuclei=('sample', 'size'), median_corrected_UMI=('total_counts', 'median'),
        median_ambient_fraction=('decontx_contamination', 'median'),
        mean_ambient_fraction=('decontx_contamination', 'mean'),
        high_contamination_fraction=('high_contamination_flag', 'mean'))
    quality.to_csv(OUT / 'fibroblast_QC_by_library.tsv', sep='\t')
    fb.groupby(['age', 'sample', 'proposed_subtype'], observed=True).size().rename('nuclei').reset_index().to_csv(
        OUT / 'historical_proposals_by_library.tsv', sep='\t', index=False)
    marker_map = {k: set(v.gene) for k, v in markers.groupby('state')}
    overlap = [{'state_1': a, 'state_2': b, 'shared_genes': ';'.join(sorted(marker_map[a] & marker_map[b]))}
        for i, a in enumerate(marker_map) for b in list(marker_map)[i+1:] if marker_map[a] & marker_map[b]]
    pd.DataFrame(overlap).to_csv(OUT / 'paper_marker_overlaps.tsv', sep='\t', index=False)
    manifest = {'status': 'PASS', 'purpose': 'Panel and existing fibroblast signal audit; not a completed new subtype or spatial analysis',
        'snRNA_fibroblasts': len(fb), 'libraries': int(fb['sample'].nunique()),
        'maturation_Figure4b_genes': len(maturation), 'all_21_maturation_genes_measured': True,
        'raw_and_corrected_per_nucleus_denominators_verified': True,
        'source_workbook_sha256': hashlib.sha256(WORKBOOK.read_bytes()).hexdigest(),
        'source_Figure4b': 'User-supplied NCVR_ms_main_021626.docx, word/media/image6.tiff; manually transcribed visible gene names',
        'published_method': 'PMC13561863 Methods: AddModuleScore, 21 genes, 100 matched control genes per module gene',
        'spatial_gene_targets': 403, 'non_gene_features': {'Unassigned Codeword': 77, 'Negative Control Codeword': 41, 'Negative Control Probe': 20},
        'limits': ['No newly fitted maturation index', 'No new cell states or trajectories',
                   'Historical proposal labels are not independent validation',
                   'Two libraries per age; do not assume nuclei are independent animals',
                   'Panel scores and ECM RNA do not measure deposited protein, stiffness, pathway causality or lineage ancestry']}
    (OUT / 'audit_validation.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(coverage).to_string(index=False), flush=True)
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    main()
