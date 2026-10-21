# Early ventricular cardiomyocyte states: within-age comparison

These comparisons address whether the separated P0/P7 ventricular groups differ beyond cycling. They use the existing annotations and expression matrices; no clusters are renamed or refitted.

## Cycling versus low-cycling groups

CAR_0 is cycling-enriched: 522/1,129 nuclei (46.2%) have at least two of Mki67, Top2a, Cenpf, Ccnb1, Ube2c and Cdk1 detected at >=1 corrected count. CAR_0 includes 325 P0 and 798 P7 nuclei (plus six P14). This is evidence of a cell-cycle transcriptional state, not proof of completed cardiomyocyte proliferation; postnatal DNA synthesis can also accompany polyploidization/binucleation. Marker-negative nuclei are not proven permanently non-proliferating.

The two dominant low-cycling P0 groups are CAR_1 (5,773 P0 nuclei) and CAR_2 (4,214). The two dominant low-cycling P7 groups are CAR_5 (1,139 P7 nuclei) and CAR_6 (10,916). Pairwise expression comparisons exclude the small number of cycling-marker-positive nuclei: P0 CAR_1=5,769, CAR_2=4,211; P7 CAR_5=1,139, CAR_6=10,916.

## Differences within P0

- **CAR_1: metabolic substrate-handling-enriched state.** Higher Hmgcs2, Cd36, Slc27a1, Pdk4 and Txnip suggest differences in lipid uptake, substrate regulation and ketogenesis-related transcription. This does not demonstrate metabolic flux, oxidative capacity, disease or greater maturity.
- **CAR_2: adhesion/cytoskeletal/calcium-associated state.** Higher Ncam1, Vcl, Smpx and Tecrl, together with Zbtb16 and Ankrd1, support a distinct regulatory/structural program. Ankrd1 expression in development is not proof of pathological stress. Ncam1 alone does not establish conduction identity.

For example, Hmgcs2 is detected at >=1 corrected count in 40.1% of CAR_1 versus 4.3% of CAR_2; Ncam1 in 62.1% versus 93.3%, respectively. These descriptive differences have the same direction in both P0 libraries and in original integer-count sensitivity checks. The selected differences also persist within Xist-positive/Y-negative nuclei, with matching directions in both libraries. Thus the P0 split is not explained solely by age, one library, sex-marker composition or ambient correction. Other technical effects, developmental state and anatomical origin remain possible; chamber-wall or regional identities are not established.

## Differences within P7

- **CAR_5: Myom2/Corin/Hmgcs2-enriched ventricular state.** Relatively higher Myom2, Corin and Hmgcs2, with Pdk4, suggest a combination of sarcomeric and metabolic regulation.
- **CAR_6: calcium-handling/contractile-associated ventricular state.** Relatively higher Srl, Jph2, Rgs6, Actn2 and Ppargc1b suggest differences in excitation–contraction and structural programs. Neither group can simply be ranked as globally more mature from these selected genes.

Selected differences persist in original counts and in both inferred sex-marker groups within the main supporting library. However, 1,115/1,139 (97.9%) P7 CAR_5 nuclei come from one library, versus only 24 from the other. Median raw counts/genes are lower in CAR_5 (2,604.5/1,329) than CAR_6 (3,475/1,560), with higher median estimated contamination (8.6% versus 4.8%). Xist/Y-linked expression composition also differs. Therefore CAR_5 is a candidate state with substantial library/technical uncertainty, not an established ventricular subtype. Two-library direction agreement for selected unstratified genes is descriptive and does not overcome the small second-library population. Matched sex-marker comparisons have adequate numbers in only one P7 library.

## Interpretation and limits

The supported description is **one cycling-enriched population plus distinct low-cycling ventricular transcriptional states**, rather than three proven functional ventricular cell types at each early age. P0 CAR_1/CAR_2 has stronger library support than the P7 CAR_5/CAR_6 split. These marker-selected descriptions are exploratory; no significance tests, pathway enrichment tests, anatomical assignments, lineage tracing or functional measurements were performed. The two libraries per age are not assumed to be independent animals. Conduction-associated CAR_3/CAR_9 comparisons across age are not performed here and do not establish that one develops into the other.

## Tables and reproducibility

- [Cycling and QC by cluster](cluster_cycling_and_qc.tsv)
- [Library and stage counts](cluster_age_library_counts.tsv)
- [Top descriptive expression differences](top_descriptive_markers.tsv)
- [P0 full gene comparison](P0_CAR_1_vs_CAR_2_all_genes.tsv.gz)
- [P7 full gene comparison](P7_CAR_5_vs_CAR_6_all_genes.tsv.gz)
- [Selected marker profiles](selected_marker_profiles.tsv)
- [Sex-marker groups by cluster and library](sex_marker_groups_by_cluster_library.tsv)
- [Detailed methods](methods.json)

Expression effects are differences of mean per-nucleus log1p(counts/total counts × 10,000), **not log2 fold changes**. Detection means >=1 corrected count, with separate original-count fractions. Comparisons use all measured genes; ranked lists exclude mitochondrial/ribosomal and generic gene-name prefixes. Within-library contrasts require >=20 nuclei in both groups. Inferred sex-marker groups require Xist detection with no Y markers or absence of Xist with >=2 detected Y markers; remaining nuclei are unclassified/ambiguous. These groups do not replace verified animal sex. Sex-matched pooled contrasts require >=30 nuclei per group. All thresholds are descriptive heuristics, not calibrated confidence probabilities.

Run `python analysis/compare_early_ventricular_states.py --input /path/to/snrna_qc_annotated_subclustered.h5ad` in the existing analysis environment.
