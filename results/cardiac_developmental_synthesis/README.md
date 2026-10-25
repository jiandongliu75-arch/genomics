# Vascular growth and the developing cardiac tissue unit

[Biological assessment report (PDF)](cardiac_developmental_biology_report.pdf) · [Developmental model](developmental_coordination_model.png) · [Cross-stage Wormhole diagnostic](cross_stage_metric_assessment.tsv)

## Overall assessment

The atlas is a sound exploratory developmental resource with a partly validated mechanistic hypothesis. Its strongest conclusion is that **CM maturation is spatially coordinated with the developing microvascular environment**. It also provides experimental support for a role of PTPRM in maturation-associated gene expression and cell-cycle regulation. It does not yet establish that capillary EC contact initiates CM maturation, that EC PTPRM is the necessary partner, or that PTPRM acts through MEIS1 to produce complete functional maturation.

The most useful biological interpretation is reciprocal coordination: **myocardial growth creates a demand for vascular expansion, while vascular access and local cell interactions help support the metabolic and structural development of muscle**. This allows EC proliferation and CM cell-cycle withdrawal to occur in the same region. The heart is assembling a working tissue unit, not synchronizing every cell type into the same proliferative state.

This assessment supersedes the biological interpretation of the earlier `cardiac_niche_dynamics` report. The numerical results remain available, but its conditional models were not a direct replication of the atlas's proximity analysis, and its original Wormhole assessment did not test cross-stage correspondence.

## What the atlas establishes, and where interpretation extends beyond evidence

| Conclusion | Assessment |
|---|---|
| Postnatal CMs acquire maturation-associated transcriptional programs with spatial heterogeneity | Well supported as a transcriptomic observation. The 21-gene index has independent developmental RNA support and gene-withholding checks. It does not measure contractile force, calcium handling, mitochondrial function or completed maturation directly. |
| At P7, higher CM maturation scores occur near capillary ECs, proliferating ECs and pericytes | Supported within the sampled tissue. The within-age comparison is a strength. These three populations also occupy related vascular structures, so their associations are not three independent demonstrations of a signaling mechanism. Anatomy, wall remodeling, cell size and local vascular density remain alternative or contributing explanations. |
| The microvascular niche contributes to CM maturation | A reasonable model supported by spatial observations and a nominated functional regulator. The data do not resolve how much comes from perfusion, contact, paracrine factors, shared anatomical development or reciprocal CM-to-EC signaling. |
| PTPRM participates in postnatal maturation-associated regulation | Supported by boundary protein localization and gene perturbation, with increased cell-cycle transcripts and reduced Meis1/maturation-related transcripts at P8. The specific validation concerns gene expression; increased cycling transcripts alone do not establish productive CM division, and the experiment does not demonstrate all dimensions of functional maturity. |
| EC PTPRM activates CM PTPRM, then MEIS1, to cause CM maturation | A plausible, incompletely tested causal chain. COMMOT and SigNET nominate spatial and regulatory compatibility; they do not measure trans-binding or establish biochemical direction. EC-specific perturbation, comparison with CM–CM contacts, and MEIS1 rescue/epistasis would distinguish the links. |
| The microvasculature is an organizing center for postnatal maturation | Useful as a working model; not established as the dominant causal hierarchy over mechanical load, hormones, nutrition, matrix remodeling or intrinsic CM programs. |

The methods have different jobs. **CytoSPACE** places RNA-derived cell states/regulatory programs into the spatial reference. The **measured Xenium maturation panel** supplies an outcome. **Graph rings** relate that outcome to neighborhood position. **CellCharter** describes recurring tissue contexts. **COMMOT** nominates expression-compatible, spatially permitted communication, and **SCENIC/SigNET** nominate downstream regulatory links. Imaging and perturbation add evidence of a different kind. Their convergence is a strength; agreement among computational scores derived from overlapping expression information is not equivalent to independent experimental replication.

The most material statistical concern is separating regional developmental organization from a local instructive effect. A compact myocardial region may simultaneously have more capillaries, more pericytes, different loading and more mature CMs. Graph rings reflect cell packing as well as physical separation; rings in growing myocardium are not fixed micrometer distances. Spatially autocorrelated cells do not provide the same replication as independent hearts. These limitations qualify causal/generalizable claims rather than invalidating the descriptive atlas.

The broad claim of vascular–CM coordination is therefore more secure than the complete proposed molecular mechanism. Sparse panel coverage, segmentation and imputed regulatory activities affect precision, but they are secondary to asking the correct biological question.

## What the perinatal vascular paper adds

Lu et al. show that a pre-existing coronary plexus expands into the developing inner myocardial wall. This begins before birth: around E18.5, VEGFR3-enriched vessels extend into trabecular bundles. Expansion is associated with regional hypoxia and preferential EC proliferation, and the inner vascular density approaches outer myocardial density by P7. The tested late-fetal/neonatal endocardial lineages do not make an appreciable contribution to this expansion; this does not negate earlier embryonic endocardial origins.

This pattern is biologically important because vascularization follows a changing tissue geometry. As the inner ventricular wall thickens and remodels, access to chamber blood across endocardial surfaces can become insufficient for the growing muscle. New intramyocardial vessels can shorten local diffusion paths and help support the increasing metabolic demands of that tissue. This is a plausible functional interpretation of the spatial pattern; the atlas itself does not measure oxygen delivery.

Higher systemic oxygen after birth and regional inner-wall hypoxia can coexist. Systemic oxygen availability rises while local diffusion distances, workload and oxygen consumption also change. The spatial response should therefore be heterogeneous, rather than a uniform switch throughout the ventricle.

DLL4/Notch regulates the pattern of this expansion: Dll4 loss enhances growth, whereas activated Notch suppresses it. More vessel growth is not automatically better vascular integration or delivery. Branch organization, lumen formation, connectivity, mural association and flow also matter. These properties cannot be read directly from an EC cycling score.

The two studies consequently support a transition from **vascular expansion in growing muscle toward a closely associated myocardial–microvascular structure**, with continued vascular adaptation as CMs enlarge and mature. They do not establish a rigid sequence in which all sprouting ends before contact or CM maturation begins.

There may also be a **change in the endothelial environment experienced by an inner-wall CM**. Early trabecular muscle has extensive association with endocardial surfaces; coronary ingrowth adds capillary and perivascular neighbors as the wall remodels. The atlas reports declining early growth-associated signaling, including NRG signaling, alongside emerging maturation-associated vascular contacts. A useful hypothesis is therefore a changing balance of endothelial influences, together with changing CM responsiveness. This is not endocardial conversion into coronary cells, nor a claim that all endocardial signals disappear or that every capillary signal is maturational. It predicts that regional endothelial context may explain some persistent differences between ventricular CM states.

## A developmental framework for P0–P21

| Interval | Tissue-level interpretation | Coordination to examine |
|---|---|---|
| Late fetal to P0 | Coronary expansion is already underway in the inner wall, alongside ventricular wall remodeling. Birth changes circulation, oxygen availability and loading. | Inner versus outer myocardium may differ in vascular access, growth demand and CM state. A single global age score will miss this organization. |
| P0–P7 | Vessels expand through growing myocardium while CMs increasingly withdraw from cycling and undergo hypertrophic growth. Fibroblast proliferation/remodeling is also prominent around P7 in the atlas. | EC proliferation can supply an expanding network to increasingly mature, enlarging CMs. Matrix remodeling can help organize vessels and transmit changing mechanical forces. These processes can overlap spatially and temporally. |
| P7–P14 | The perinatal paper describes declining inner-plexus growth activity, while the atlas identifies a marked CM regulatory-program transition between P7 and P14. | A useful hypothesis is that vascular organization and access help enable continued metabolic, structural and electrophysiological development after the initial expansion wave. Contact-dependent signaling may reinforce cell-cycle exit; it need not explain every maturation component. |
| P14–P21 | CM growth and functional specialization continue. Lower proliferation does not imply that vascular development is complete. | Determine whether established capillary/mural/matrix neighborhoods adapt to larger CMs and greater workload, rather than assuming that more mature tissue must contain more cycling ECs. |

The interval assignments summarize published observations and physiological hypotheses, not newly measured developmental trajectories. Changes in LV versus RV loading, nutrition and endocrine signaling are additional physiological drivers. They can coordinate both vessels and CMs and should not be attributed to EC contact alone.

Fibroblasts and immune cells belong in this framework because the developing tissue needs a changing matrix and organized intercellular support. The observed early fibroblast–immune contact enrichment is a candidate component of that remodeling environment. It does not establish which immune signals are active or whether they cause CM maturation. Likewise, a developmental matrix program is not automatically pathological fibrosis.

## The main mechanistic possibilities

**Permissive support:** vascular access makes sustained oxidative metabolism, growth and contraction possible. This predicts particularly strong relationships with metabolic/energetic CM programs and local vascular access. It requires perfusion/oxygen measurements for a direct physiological test.

**Instructive signaling:** a local EC/perivascular or matrix interaction changes a CM regulatory program. PTPRM is a specific candidate because protein is found at EC–CM interfaces and its perturbation changes maturation-associated expression. Its homophilic interaction requires competence on both sides. The rise in CM PTPRM after birth suggests that acquiring CM responsiveness may be as important as changing EC expression. CM–CM interactions are a competing or complementary source of the same molecular contact.

**Reciprocal demand matching:** growing, metabolically changing CMs stimulate vascular adaptation, while the resulting vascular environment supports further growth and specialization. The atlas's CM-to-EC growth-signal predictions are consistent with this reciprocal model. In that case, endothelial proliferation near mature CMs is an expected feature of coordinated development.

These mechanisms can coexist. Cell-cycle withdrawal, sarcomere organization, calcium handling and metabolic specialization should be evaluated separately, because neither a single composite score nor a MEIS1-associated program captures all of maturation.

## What Wormhole could contribute to these biological questions

Wormhole learns a fast representation of distances between distributions of cells. It can help compare neighborhood populations across samples, especially when a neighborhood contains a mixture of cell states that an average conceals. It does not itself learn developmental time, track cells or identify causal signals.

The useful cross-stage questions are:

1. **Persistence versus remodeling:** do similar CM/EC/pericyte/fibroblast neighborhoods persist from P0 to P21 while their internal programs change, or do new combinations become prominent?
2. **Expansion versus state change:** is a P7–P14 difference explained by changing numbers/proportions of cells, or by functional changes within otherwise similar cellular mixtures?
3. **Regional convergence:** do inner myocardial neighborhoods become more similar to outer myocardial neighborhoods as the coronary plexus expands? The perinatal paper gives this comparison a specific biological rationale.
4. **Asynchrony:** within an age, do some neighborhoods resemble later-stage reference neighborhoods, and are their CMs more advanced in independently measured maturation components?
5. **Coupling:** do changes in vascular, matrix and CM programs accompany one another, or do particular aspects remain uncoupled?

For a descriptive developmental map, joint CM and non-CM states can be included. For testing whether a non-CM environment relates to CM maturation, CM outcome features should be withheld and related to the inferred neighborhood afterward. Those are different representations for different questions.

A neighborhood's resemblance to a later sample must also be distinguished from anatomical identity: a persistent inner-wall state need not be an immature version of an outer-wall state. Anatomical matching and separate CM functional programs are necessary before describing a region as developmentally advanced or delayed.

Our first benchmark only tested within-age distances and within-age maturation prediction. It could not establish whether Wormhole was useful for these between-stage questions. The present extension adds **768 reserved cross-stage pairs**, including comparisons selected for similar lineage composition. Across unrelated cross-stage mixtures, Wormhole generally preserves broad distance ranking reasonably well. When composition is approximately matched, its agreement is substantially weaker, and simple mean features usually perform better. This matters biologically: the present encoder is less convincing precisely for the subtle within-neighborhood state changes we want to study. See [cross-stage assessment](cross_stage_metric_assessment.tsv).

The current input also omits within-neighborhood spatial arrangement and normalizes each cloud into a probability distribution. It therefore cannot directly identify which EC touches a CM or measure absolute vascular expansion. A useful biological representation should retain vascular density/access and contact measurements separately, use comparable anatomical compartments, and distinguish absolute abundance from cell-state composition. Our composition-matched diagnostic does not yet anatomically register the neighborhoods.

Wormhole remains a potentially useful scalable comparison tool. The current configuration is not validated for subtle developmental state alignment, and its limited performance does not mean that developmental coordination is absent. Exact transport or simpler summaries remain appropriate when they answer the same question reliably. A convincing developmental correspondence would also need stability across gene sets/count views and support from independent biological samples; stage and section remain confounded in the current spatial series.

## What should lead the next analysis

The central unit should be **a ventricular tissue neighborhood with a CM phenotype, vascular access/organization and a matrix environment**. Start with comparable anatomical regions and the authors' CM populations, then extend to the broader ventricular states. Our earlier P7 cohort contained 39,383 ventricular-like CMs, whereas the paper's reported P7 maturation cohort contains 11,086 author-labeled v.CMs; the analyses address different populations.

First determine whether vascular association concerns CM metabolism, structure, electrical development or cell-cycle withdrawal, and whether it follows regional wall remodeling. Then test whether interface-level PTPRM compatibility explains additional differences beyond vascular access, without automatically adjusting away the CM response that may mediate the mechanism. Nuclear/core masks remain valuable sensitivity views, but are not substitutes for membrane protein localization.

This strategy can produce a biologically interpretable account of how growing vessels, matrix and myocardium develop together. The decisive causal experiments would separate EC-specific contact effects from perfusion and CM-intrinsic effects, with functional CM readouts and a test of the proposed downstream mediator. The existing data can prioritize those mechanisms; additional embedding alone cannot decide among them.

## Sources and reproducibility

- [Wang et al., Spatially guided in vivo single-cell functional genomics of postnatal heart](https://doi.org/10.1038/s44161-026-00861-z), especially Fig.2–4, Extended Data7 and Methods. Full cached [PMC13561863](https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/) text/XML reviewed.
- [Lu et al., Perinatal angiogenesis from pre-existing coronary vessels via DLL4–NOTCH1 signalling](https://doi.org/10.1038/s41556-021-00747-1), Nature Cell Biology23,967–977(2021). Full uploaded PDF reviewed; it is not redistributed.
- [Puente et al., The oxygen-rich postnatal environment induces cardiomyocyte cell-cycle arrest through DNA damage response](https://doi.org/10.1016/j.cell.2014.03.032), Cell(2014), for the distinction between oxygen-related cell-cycle withdrawal and full maturation.
- [Mahmoud et al., Meis1 regulates postnatal cardiomyocyte cell cycle arrest](https://doi.org/10.1038/nature12054), Nature(2013), for the downstream regulator's established cell-cycle role.

The diagnostic reuses the saved encoder and held-out data from the earlier report. No new model training, source reassignment or original-count modification was performed. Source hashes, execution settings and [validation](cross_stage_validation.json) accompany the results. Run `analysis/cardiac_developmental_synthesis.py` in the pinned isolated Wormhole environment from the earlier release. The cross-stage pair sampling seed is 3841. Numerical validation does not validate biological correspondences or lineages.
