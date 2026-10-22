"""Assemble an evidence-qualified catalog from the published fibroblast analyses.

This is a synthesis, not a new clustering or spatial transfer. Functional readings
are hypotheses; original and reviewed annotations remain separately identifiable.
"""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fibroblast_subtype_catalog'
SOURCE = ROOT / 'results/fibroblast_states'

# Descriptions refer to marker programs, not experimentally established functions.
READINGS = {
    'P0_F1': ('Sntg1/Slit3/Ptn-associated', 'Candidate adhesion/growth-factor-associated developmental program.'),
    'P0_F2': ('Eln/Slit3-associated', 'Candidate elastic-matrix/adhesion program; no outflow-tract assignment.'),
    'P0_F3': ('Lamb1/Lamc1/Lama2-associated', 'Candidate laminin/basement-membrane assembly program.'),
    'P0_F4': ('Kirrel3/Erbb4-associated', 'Distinct adhesion/signaling-associated program; valve identity unverified.'),
    'P0_F5': ('Cycling-associated Postn', 'DNA replication/cell-cycle activity; completed division not established.'),
    'P0_F6': ('Fn1/Frem1/Hmcn1-associated', 'Candidate matrix organization/adhesion program.'),
    'P7_F1': ('Col24a1/Tnc-associated', 'Candidate developmental matrix-remodeling program.'),
    'P7_F2': ('Col15a1/Lamc1/Abca8a-associated', 'Candidate matrix/basement-membrane-associated program.'),
    'P7_F3': ('Sntg1/Ptn-associated', 'Candidate developmental adhesion/growth-factor-associated program.'),
    'P7_F4': ('Gfpt2/Sema3c-associated', 'Candidate glycosylation-related metabolic and matrix-organizing program; no hyaluronan-flux claim.'),
    'P7_F6': ('Bnc2/Cdh11/Lsamp-associated', 'Candidate adhesion-associated overlapping program.'),
    'P7_F7': ('Cycling-associated Postn', 'DNA replication/cell-cycle activity; completed division not established.'),
    'P7_F8': ('Kirrel3/Erbb4-associated', 'Distinct adhesion/signaling-associated program; valve identity unverified.'),
    'P14_F0': ('Gfpt2/Adgrd1/Cmah-associated', 'Gfpt2-family matrix/metabolic program; no hyaluronan-flux claim.'),
    'P14_F1': ('Hmcn2/Fbln2/Pcsk6-associated', 'Candidate matrix organization/processing program.'),
    'P14_F2': ('Adgrl3/Svep1/Fap-associated', 'Candidate matrix organization/adhesion program; no pathological-fibrosis claim.'),
    'P14_F3': ('Lsamp/Eln/Slit3-associated', 'Candidate elastic-matrix and adhesion/signaling program; no outflow-tract assignment.'),
    'P14_F4': ('Col8a1/Lsamp/Itgbl1 overlap', 'Overlapping matrix/adhesion program; finer identity sensitivity-dependent.'),
    'P14_F5': ('Kirrel3/Erbb4-associated', 'Distinct adhesion/signaling-associated program; valve identity unverified.'),
    'P14_F7': ('Small cycling-associated candidate', 'Cell-cycle activity in 35 nuclei; broad family unresolved.'),
    'P21_F0': ('Adgrl3/Svep1-associated late interstitial', 'Candidate matrix organization/adhesion substate within the reviewed Col8a1-related family.'),
    'P21_F1': ('Lsamp/Slit3/Rgs6-associated', 'Candidate elastic-matrix and adhesion/signaling substate.'),
    'P21_F2': ('Hmcn2/Pcsk6-associated', 'Candidate matrix organization/processing substate within the Col8a1 family.'),
    'P21_F3': ('Gfpt2/Adgrd1/Cmah-associated', 'Gfpt2-family matrix/metabolic program; no hyaluronan-flux claim.'),
    'P21_F4': ('Kirrel3/Erbb4-associated', 'Distinct adhesion/signaling-associated program; valve identity unverified.'),
    'P21_F5': ('Lsamp/Col4a5/Col4a6/Itgbl1-associated', 'Candidate elastic-matrix plus collagen-IV/adhesion substate.'),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SOURCE / 'RNA_cluster_summary.tsv', sep='\t')
    df['original_family'] = df.family
    df['reviewed_broad_family'] = df.family
    df.loc[df.cluster_id.eq('P21_F0'), 'reviewed_broad_family'] = 'Col8a1-associated'
    df['descriptive_substate'] = df.cluster_id.map(lambda x: READINGS.get(x, ('Residual lineage-program review', ''))[0])
    df['functional_interpretation_hypothesis'] = df.cluster_id.map(lambda x: READINGS.get(x, ('', 'Residual CM/mural program; excluded from fine reference and primary functional comparisons.'))[1])
    df['spatial_resolution'] = 'Family-level only where Xenium support passes; no subcluster-specific anatomical placement established.'
    df.loc[df.lineage_mixing_review_flag, 'spatial_resolution'] = 'Withheld from fine reference; not assigned a fibroblast niche.'
    df.loc[df.cluster_id.eq('P14_F7'), 'spatial_resolution'] = 'Broad family unresolved; no supported fine spatial assignment.'
    df['identity_review'] = 'Original analysis retained.'
    df.loc[df.cluster_id.eq('P21_F0'), 'identity_review'] = 'Focused P21 review supports Col8a1-related broad family; existing Xenium transfer/models not refitted.'
    cols = ['stage', 'cluster_id', 'nuclei', 'original_family', 'reviewed_broad_family',
            'descriptive_substate', 'top_markers', 'functional_interpretation_hypothesis',
            'support', 'cycling_fraction', 'lineage_mixing_review_flag', 'spatial_resolution', 'identity_review']
    df[cols].to_csv(OUT / 'all_stage_subtypes.tsv', sep='\t', index=False)
    eligible = df.loc[~df.lineage_mixing_review_flag]
    counts = eligible.pivot_table(index='reviewed_broad_family', columns='stage', values='nuclei', aggfunc='sum', fill_value=0)
    counts = counts.reindex(columns=['P0', 'P7', 'P14', 'P21'])
    counts.to_csv(OUT / 'reviewed_family_nucleus_counts.tsv', sep='\t')
    markers = pd.read_csv(ROOT / 'results/fibroblast_reassessment/paper_fibroblast_markers_and_coverage.tsv', sep='\t')
    markers.to_csv(OUT / 'paper_family_markers_and_assay_coverage.tsv', sep='\t', index=False)
    assert len(df) == 31 and int(df.nuclei.sum()) == 14985
    assert int(eligible.nuclei.sum()) == 14328
    assert int(counts.loc['Col8a1-associated', 'P21']) == 2711
    assert int(df.loc[df.lineage_mixing_review_flag, 'nuclei'].sum()) == 657
    assert set(eligible.cluster_id) == set(READINGS)
    check = {'status': 'PASS', 'computational_groups': len(df), 'source_nuclei': int(df.nuclei.sum()),
             'reference_eligible_nuclei': int(eligible.nuclei.sum()), 'lineage_review_nuclei': 657,
             'reviewed_P21_Col8a1_family_nuclei': 2711,
             'new_clustering_or_spatial_transfer': False,
             'source_analysis': 'results/fibroblast_states',
             'additive_identity_review': 'results/fibroblast_P21_identity_review'}
    (OUT / 'validation.json').write_text(json.dumps(check, indent=2) + '\n')
    text = ['# Fibroblast subtype catalog: P0–P21', '',
            'This synthesis retains five broad RNA families, within-age substates and a cycling overlay. '
            'It incorporates the focused P21_F0 identity review. Functional descriptions are marker-based '
            'hypotheses, and fine anatomical placement remains incomplete. It does not recompute a UMAP, '
            'change source counts or refit the preceding Xenium maps.', '',
            '[Family functions and spatial interpretation](INTERPRETATION.md) · '
            '[All 31 groups: markers, support, functional interpretation and spatial limits](all_stage_subtypes.tsv) · '
            '[Reviewed family counts](reviewed_family_nucleus_counts.tsv) · '
            '[Paper marker lists and assay coverage](paper_family_markers_and_assay_coverage.tsv) · '
            '[Validation](validation.json)', '',
            'All 14,985 nuclei are accounted for. The 657 nuclei in five residual CM/mural-program groups '
            'remain flagged; the 35-nucleus P14 cycling group remains family-unresolved. Thus not every '
            'computational cluster is a resolved biological subtype.', '',
            '## Stage-specific catalog', '',
            'F numbers are local to each stage. Marker lists are relative within-stage enrichments; '
            'a gene need not be exclusive or detected in every nucleus. Some rows are explicitly '
            'sensitivity-dependent in the source analysis.', '']
    for stage in ['P0', 'P7', 'P14', 'P21']:
        text += [f'### {stage}', '', '| Group | Nuclei | Reviewed family | Substate / representative enriched genes | Functional hypothesis | RNA support |',
                 '|---|---:|---|---|---|---|']
        for _, r in df.loc[df.stage.eq(stage)].iterrows():
            genes = ', '.join(str(r.top_markers).split(';')[:6]) if pd.notna(r.top_markers) else 'Not assigned'
            text.append(f'| {r.cluster_id} | {r.nuclei} | {r.reviewed_broad_family} | {r.descriptive_substate}: {genes} | {r.functional_interpretation_hypothesis} | {r.support} |')
        text.append('')
    (OUT / 'README.md').write_text('\n'.join(text).rstrip() + '\n')
    print(json.dumps(check))


if __name__ == '__main__':
    main()
