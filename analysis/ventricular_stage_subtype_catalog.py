"""Consolidate existing ventricular groups and describe cross-stage resemblance.

Aggregate-count signatures are exploratory comparisons, not matched-cell lineage
inference or biological-replicate differential expression.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/ventricular_microenvironment'
OUT = ROOT / 'results/ventricular_stage_subtype_catalog'
STAGES = ['P0', 'P7', 'P14', 'P21']

# Functional hypotheses reflect previously evaluated multi-gene programs.
# Marker-only groups explicitly retain uncertain subtype-specific functions.
CURATION = {
    'P0_M0': ('Cycling-associated CM state', 'Cycling', 'Cell-cycle progression and DNA replication; completed cell division is not demonstrated.', 'Mki67;Top2a;Ect2;Diaph3'),
    'P0_M1': ('Ddc-associated metabolic-transition candidate', 'Early metabolic/working', 'Potential adaptation of substrate use after birth, including lipid/ketone-associated expression; Ddc-specific cardiac function unresolved.', 'Ddc;Hmgcs2;Pdk4;Kcnd2'),
    'P0_M2': ('Zbtb16/Ncam1-enriched neonatal working-CM candidate', 'Marker-defined', 'Possible cell-contact and structural/electrical maturation; subtype-specific function not established.', 'Zbtb16;Ncam1;Tecrl;Vcl'),
    'P0_M3': ('Ankrd1-associated remodeling candidate A', 'Mechanical-response-associated', 'Potential sarcomere/cytoskeletal adaptation and load-responsive signaling; actual load is not measured.', 'Ankrd1;Creb5;Nrap;Tlr4'),
    'P0_M4': ('Fatty-acid-use-associated candidate', 'Early metabolic/working', 'Potential lipid uptake and fatty-acid use during postnatal metabolic adaptation; metabolic flux is not measured.', 'Cd36;Hmgcs2;Vwa8;Ntn1'),
    'P0_M5': ('Slit2-associated, conduction-program-enriched candidate', 'Conduction-program-enriched', 'Potential electrical propagation specialization within an endomyocardial/trabecular-associated program; not verified Purkinje identity.', 'Slit2;Mpped2;Etv1;Gja5;Scn10a'),
    'P0_M6': ('Ankrd1-associated remodeling candidate B', 'Mechanical-response-associated', 'Potential sarcomere and cytoskeletal coupling/remodeling; distinction from P0_M3 remains exploratory.', 'Ankrd1;Creb5;Nrap;Vcl;Lmod2'),
    'P0_M7': ('Pola1/Rad51b-associated DNA-maintenance candidate', 'DNA-replication/maintenance-associated', 'DNA-synthesis/chromatin-maintenance gene program; active proliferation or increased DNA damage is not established.', 'Pola1;Rad51b;Prim2;Hells'),
    'P7_M0': ('Ppara/Rreb1-enriched working-CM candidate', 'Metabolic-regulatory working', 'Potential metabolic regulation and contractile maturation; no qualifying standalone metabolic subtype panel.', 'Ppara;Rreb1;Ppargc1a;Kcnd2'),
    'P7_M1': ('Ankrd1-associated remodeling state', 'Mechanical-response-associated', 'Potential force-transmission apparatus remodeling and load-responsive signaling.', 'Ankrd1;Creb5;Xirp2;Nrap'),
    'P7_M2': ('Fgf14/Pdzrn4 regional electrical candidate', 'Regional electrical-associated', 'Potential electrical/RNA-processing specialization; broader conduction identity is not established.', 'Fgf14;Pdzrn4;Rbfox1;Mpped2;Scn5a'),
    'P7_M3': ('Cycling-associated CM state', 'Cycling', 'Cell-cycle progression and DNA replication; completed cell division is not demonstrated.', 'Mki67;Top2a;Ect2;Diaph3'),
    'P7_M4': ('Ivns1abp/Myom2-enriched working-CM candidate', 'Working/contractile', 'Potential contractile apparatus maturation; subtype-specific role remains uncertain.', 'Ivns1abp;Myom2;9530026P05Rik'),
    'P7_M5': ('Rad51b/Pola1-associated DNA-maintenance candidate', 'DNA-replication/maintenance-associated', 'DNA-synthesis/chromatin-maintenance gene program; active proliferation or increased DNA damage is not established.', 'Rad51b;Pola1;Prim2;Hells'),
    'P7_M6': ('Slit2-associated, conduction-program-enriched candidate', 'Conduction-program-enriched', 'Potential rapid-propagation specialization; not anatomically or functionally verified Purkinje cells.', 'Slit2;Mpped2;Etv1;Gja5;Cntn2;Scn10a'),
    'P7_M7': ('Small uncharacterized ventricular-CM candidate A', 'Unassigned', 'Subtype-specific function cannot be assigned from current evidence.', ''),
    'P7_M8': ('Small uncharacterized ventricular-CM candidate B', 'Unassigned', 'Subtype-specific function cannot be assigned from current evidence.', ''),
    'P14_M0': ('Slit2/Mpped2 regional electrical candidate', 'Regional electrical-associated', 'Potential regional repolarization specialization: relatively Irx5-enriched/Kcnd2-reduced, without the same broad conduction-marker enrichment.', 'Fgf14;Mpped2;Slit2;Rbfox1;Irx5'),
    'P14_M1': ('Lrtm1-enriched working-CM candidate', 'Working/contractile', 'Ventricular contraction and maturation background; subtype-specific function remains unresolved.', 'Lrtm1'),
    'P14_M2': ('Ankrd1-associated remodeling state', 'Mechanical-response-associated', 'Potential cytoskeletal/sarcomere remodeling and load-responsive signaling.', 'Ankrd1;Creb5;Xirp2;Tlr4;Nrap'),
    'P14_M3': ('Kcnd2/Pdzrn4-enriched repolarization candidate', 'Repolarization-associated', 'Potential relatively greater Kv4.2-related repolarization capacity; RNA does not establish potassium current.', 'Kcnd2;Pdzrn4;Fgf14;Lrtm1'),
    'P14_M4': ('Slit2-associated, conduction-program-enriched candidate', 'Conduction-program-enriched', 'Potential specialized ventricular propagation; subtype boundaries are sensitivity-dependent.', 'Slit2;Mpped2;Etv1;Gja5;Cntn2;Scn10a;Pcp4'),
    'P21_M0': ('Lrtm1/Stk39-enriched working/electrical candidate', 'Working/repolarization candidate', 'Potential ion-regulatory and working-CM specialization; broader profile resembles P14_M3 more than P14_M1.', 'Lrtm1;Stk39;Rp1;Pdk4'),
    'P21_M1': ('Marker-poor working-CM candidate', 'Working/contractile', 'General contractile/maturation background; no reproducible positive subtype-marker set or specific new function established.', ''),
    'P21_M2': ('Ankrd1-associated remodeling candidate A', 'Mechanical-response-associated', 'Potential sarcomere/cytoskeletal adaptation and load-responsive signaling.', 'Ankrd1;Creb5;Xirp2;Ankrd23;Nrap'),
    'P21_M3': ('Slit2/Mpped2 regional electrical candidate', 'Regional electrical-associated', 'Potential regional repolarization specialization: relatively Irx5-enriched/Kcnd2-reduced, with retained Gja1 and little Cntn2/Scn10a enrichment.', 'Fgf14;Mpped2;Slit2;Pdzrn4;Irx5'),
    'P21_M4': ('Ankrd1-associated remodeling candidate B', 'Mechanical-response-associated', 'Weaker Ankrd1-associated adaptation program; separation from P21_M2 and specific function remain exploratory.', 'Ankrd1;Prune2;Grm1;Stk39;Itga9'),
    'P21_M5': ('Small Slit2-associated conduction-program candidate', 'Conduction-program-enriched', 'Potential specialized ventricular propagation; small sample and historical gate dependence limit subtype confidence.', 'Slit2;Etv1;Cntn2;Scn10a;Pcp4'),
}

COMPARISONS = [
    ('Mechanical-response-associated', 'P0_M3;P0_M6', 'P7_M1', 'P14_M2', 'P21_M2;P21_M4', 'Primary all-stage series; P0_M3/M6 and P21_M2/M4 are exploratory subdivisions. P21_M2 is the closer late-stage signature match.', 'Ankrd1;Creb5;Nrap;Xirp2;cytoskeletal genes'),
    ('Conduction-program-enriched', 'P0_M5', 'P7_M6', 'P14_M4', 'P21_M5', 'Biologically coherent all-stage candidate series; P14 boundary sensitivity and P21 n=72 limit confidence. Not a lineage trajectory or verified Purkinje series.', 'Slit2;Mpped2;Etv1;Gja5;Cntn2;Scn10a'),
    ('Regional electrical-associated', '', 'P7_M2', 'P14_M0', 'P21_M3', 'Prioritize P14_M0 versus P21_M3. P7_M2 is a partial earlier counterpart; no clean independent P0 counterpart is resolved.', 'Fgf14;Rbfox1;Mpped2;Pdzrn4;Irx5/Kcnd2 axis'),
    ('Cycling', 'P0_M0', 'P7_M3', '', '', 'Strong early-stage comparison; no analogous discrete P14/P21 group resolved. This does not prove absence of rare cycling nuclei.', 'Mki67;Top2a;Ect2;Diaph3'),
    ('DNA-replication/maintenance-associated', 'P0_M7', 'P7_M5', '', '', 'Early candidate RNA-state comparison; do not equate with productive proliferation, S-phase assignment or DNA damage.', 'Pola1;Rad51b;Prim2;Hells'),
    ('Early metabolic/working arm A', 'P0_M1', 'P7_M4', '', '', 'Exploratory early working-CM comparison; Ddc does not establish a stable four-stage subtype.', 'Ivns1abp;9530026P05Rik;metabolic and contractile context'),
    ('Early metabolic/working arm B', 'P0_M4', 'P7_M0', '', '', 'Exploratory metabolic adaptation comparison; age effects and state identity are not separated fully.', 'Lipid-use and metabolic-regulatory genes'),
    ('Later repolarization/working candidate', '', 'P7_M0', 'P14_M3', 'P21_M0', 'P14_M3 versus P21_M0 is the useful late-stage pair. The P7 extension is weak and should remain exploratory.', 'Kcnd2;Pdzrn4;Lrtm1;Stk39 context'),
    ('Later working/maturation background', '', '', 'P14_M1', 'P21_M1', 'Broader gene-profile resemblance supports comparison, but both have weak subtype-specific positive markers. Avoid naming a novel functional subtype.', 'Broad relative transcriptome; few robust positive subtype markers'),
]


def main():
    OUT.mkdir(exist_ok=True)
    summary = pd.read_csv(SOURCE / 'cluster_summary.tsv', sep='\t')
    assert set(summary.microenvironment_cluster) == set(CURATION)
    catalog = summary.copy()
    fields = ['candidate_identity', 'functional_family', 'potential_function_and_limits', 'interpretation_markers']
    for i, field in enumerate(fields):
        catalog[field] = catalog.microenvironment_cluster.map(lambda cluster: CURATION[cluster][i]).replace('', pd.NA)
    totals = catalog.groupby('stage').nuclei.transform('sum')
    catalog['fraction_within_retained_stage_CMs'] = catalog.nuclei / totals
    catalog.to_csv(OUT / 'stage_subtype_catalog.tsv', sep='\t', index=False, na_rep='NA')
    comparison = pd.DataFrame(COMPARISONS, columns=['comparison_family', *STAGES, 'priority_and_limits', 'shared_program'])
    comparison.to_csv(OUT / 'cross_stage_comparison_sets.tsv', sep='\t', index=False, na_rep='NA')
    cells = pd.read_csv(SOURCE / 'cell_annotations_and_scores.csv.gz', index_col=0)
    cells = cells[['age', 'sample', 'microenvironment_cluster', 'gate_origin']].copy()
    lookup = catalog.set_index('microenvironment_cluster')
    for field in ['candidate_identity', 'functional_family', 'support_status']:
        cells[field] = cells.microenvironment_cluster.map(lookup[field])
    membership = {cluster: [] for cluster in CURATION}
    for _, row in comparison.iterrows():
        for stage in STAGES:
            for cluster in filter(None, row[stage].split(';')):
                membership[cluster].append(row.comparison_family)
    cells['proposed_comparison_memberships'] = cells.microenvironment_cluster.map(lambda cluster: ';'.join(membership[cluster]) or pd.NA)
    cells.to_csv(OUT / 'cell_catalog_and_comparison_memberships.tsv.gz', sep='\t', na_rep='NA')
    counts = []
    symbols = None
    for stage in STAGES:
        table = pd.read_csv(SOURCE / stage / 'corrected_pseudobulk.tsv.gz', sep='\t', index_col=0)
        if symbols is None:
            symbols = table.gene_symbol.astype(str)
        else:
            assert table.index.equals(symbols.index) and table.gene_symbol.astype(str).equals(symbols)
        counts.append(table.drop(columns='gene_symbol').astype(float))
    allcounts = pd.concat(counts, axis=1)
    assert set(allcounts.columns) == set(CURATION)
    # Denominators retain ALL original features, before duplicate-symbol handling.
    totals = allcounts.sum(axis=0)
    abundance = allcounts.sum(axis=1)
    mapping = pd.DataFrame({'gene_id': allcounts.index, 'gene_symbol': symbols.to_numpy(), 'pooled_corrected_counts': abundance.to_numpy()})
    chosen = mapping.sort_values('pooled_corrected_counts', ascending=False).drop_duplicates('gene_symbol').gene_id
    unique = allcounts.loc[chosen]
    gene_symbols = symbols.loc[chosen]
    expr = np.log1p(unique * 1e6 / totals)
    signatures = expr.copy()
    for stage in STAGES:
        cols = [col for col in unique if col.startswith(stage + '_')]
        rest = unique[cols].sum(axis=1).to_numpy()[:, None] - unique[cols].to_numpy()
        rest_totals = totals[cols].sum() - totals[cols].to_numpy()
        signatures[cols] -= np.log1p(np.maximum(rest, 0) * 1e6 / rest_totals[None, :])
    nuisance = gene_symbols.str.match(r'^(mt-|Rpl|Rps|Gm)') | gene_symbols.isin(['Xist', 'Tsix', 'Ddx3y', 'Eif2s3y', 'Kdm5d', 'Uty', 'Usp9y', 'Zfy1', 'Zfy2'])
    keep = ~nuisance & expr.max(axis=1).ge(np.log1p(5)) & signatures.abs().max(axis=1).ge(0.4)
    config = json.loads((SOURCE / 'program_definitions.json').read_text())['programs']
    withheld = set(g for value in config.values() for g in value['genes'])
    withheld.update(json.loads((ROOT / 'results/ventricular_additional_axes/methods.json').read_text())['genes_requested'])
    withheld.update(json.loads((ROOT / 'results/ventricular_electrical_regional/methods.json').read_text())['genes'])
    sensitivity = keep & ~gene_symbols.isin(withheld)
    primary = signatures.loc[keep].corr(method='spearman')
    excluded = signatures.loc[sensitivity].corr(method='spearman')
    features = mapping.set_index('gene_id').loc[signatures.index].copy()
    features['primary_similarity_feature'] = keep
    features['interpreted_gene_excluded_similarity_feature'] = sensitivity
    features.to_csv(OUT / 'signature_feature_mapping.tsv', sep='\t', na_rep='NA')
    signatures.loc[keep].to_csv(OUT / 'within_stage_relative_aggregate_signatures.tsv.gz', sep='\t', float_format='%.7g')
    marker = pd.read_csv(SOURCE / 'cluster_gene_markers.tsv', sep='\t')
    marker = marker[marker.libraries_with_positive_corrected_and_raw_effect.eq(2) & marker.corrected_log_difference.gt(0) & marker.raw_log_difference.gt(0)]
    sets = {cluster: set(group.gene_symbol.astype(str)) for cluster, group in marker.groupby('cluster')}
    pair_rows = []
    for source in primary.index:
        for stage in STAGES:
            if source.startswith(stage + '_'):
                continue
            target_ids = [col for col in primary if col.startswith(stage + '_')]
            rankings = primary.loc[source, target_ids].rank(ascending=False, method='min')
            for target in target_ids:
                a, b = sets.get(source, set()), sets.get(target, set())
                shared = sorted(a & b)
                reverse = [col for col in primary if col.startswith(source.split('_')[0] + '_')]
                reciprocal = primary.loc[target, reverse].idxmax() == source and rankings[target] == 1
                pair_rows.append({'source_cluster': source, 'source_stage': source.split('_')[0], 'target_cluster': target, 'target_stage': stage, 'relative_signature_spearman': primary.loc[source, target], 'interpreted_genes_excluded_spearman': excluded.loc[source, target], 'rank_within_target_stage': int(rankings[target]), 'reciprocal_top_match_between_stages': reciprocal, 'shared_two_library_positive_markers': len(shared), 'marker_jaccard': len(shared) / len(a | b) if a | b else 0., 'shared_marker_symbols': ';'.join(shared) or pd.NA})
    pd.DataFrame(pair_rows).to_csv(OUT / 'cross_stage_pairwise_similarity.tsv', sep='\t', index=False, na_rep='NA')
    primary.to_csv(OUT / 'relative_signature_similarity_matrix.tsv', sep='\t')
    excluded.to_csv(OUT / 'interpreted_genes_excluded_similarity_matrix.tsv', sep='\t')
    methods = {'nuclei': int(catalog.nuclei.sum()), 'groups': len(catalog), 'stages': STAGES, 'primary_similarity_features': int(keep.sum()), 'interpreted_genes_excluded_features': int(sensitivity.sum()), 'expression': 'log1p(CPM) of corrected-count sums by existing group; ALL original features retained in normalization denominators. Highest pooled-abundance feature selected per duplicate gene symbol.', 'signature': 'For each gene and group: aggregate log1p(CPM) minus aggregate log1p(CPM) of all other retained ventricular-enriched groups in the same stage.', 'feature_rules': 'Exclude mt-/Rpl/Rps/Gm-prefixed and listed sex-associated symbols. Require >=5 aggregate CPM in at least one group and absolute within-stage contrast>=0.4 in at least one group. Thresholds are exploratory, not calibrated.', 'similarity': 'Spearman correlation across selected relative gene signatures. Export every cross-stage pair, rank and reciprocal top match; ranks always exist even for weak similarities and do not establish identity.', 'sensitivity': 'Repeat after withholding every gene in previously interpreted microenvironment programs, 94-gene literature checks and 17-gene electrical/regional checks.', 'marker_overlap': 'Shared top80 previously selected markers with positive corrected/raw effects in both libraries; truncated discovery lists, not an independent validation set.', 'limits': ['Aggregate signatures pool libraries and can be affected by library weighting, relative stage composition, age and technical differences. This is not sample-level differential expression.', 'Previous within-library/raw/QC support is retained in catalog; new whole-transcriptome matching is corrected-count only, with interpreted-gene exclusion sensitivity.', 'No batch-adjusted counts, forced one-to-one mapping, trajectory, ancestry or functional verification.', 'A family can split/merge with development. Missing counterparts mean not resolved, not biological absence.', 'Stage-specific UMAP coordinates are not used to match groups.']}
    (OUT / 'methods.json').write_text(json.dumps(methods, indent=2))
    print(f'Cataloged {len(catalog)} groups / {catalog.nuclei.sum()} nuclei; {keep.sum()} primary and {sensitivity.sum()} excluded-panel similarity features.', flush=True)


if __name__ == '__main__':
    main()
