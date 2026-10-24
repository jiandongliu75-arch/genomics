"""Validate and summarize whole-section segmentation and fibroblast evidence."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/nuclear-anchored-mpl')
import hashlib
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from nuclear_anchored_Xenium_segmentation import ROOT,OUT,SOURCE,ATLAS,STAGES,write_json


def main():
    rows=[];tierrows=[];checks=[];mixtures=[];gate_rows=[];resolution_changes=[]
    for stage in STAGES:
        folder=OUT/stage
        seg=json.loads((folder/'segmentation_complete.json').read_text())
        count=json.loads((folder/'count_complete.json').read_text())
        proj=json.loads((folder/'fibroblast_projection_complete.json').read_text())
        rows.append({'stage':stage,**{k:seg[k] for k in ['profiles','nuclear_seeded_profiles','fallback_seeded_profiles','unseeded_profiles']},
            **{k:count[k] for k in ['all_records','Q20_gene_molecules','native_assigned','nuclear_anchored_assigned',
               'native_nuclear_molecules','nonnuclear_changed_owner','nuclear_vs_projected_owner_disagreement']},
            **{k:proj[k] for k in ['broad_fibroblasts','native_supported','nuclear_anchored_supported',
               'nuclear_only_informative','nuclear_modality_calibrated_supported','cell_nuclear_concordant','nuclear_cell_disagreement']}})
        tiers=pd.read_csv(folder/'fibroblast_evidence_tiers.tsv',sep='\t');tiers.insert(0,'stage',stage);tierrows.append(tiers)
        atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        data=ad.read_h5ad(folder/'nuclear_anchored_counts.h5ad')
        assert data.obs_names.equals(atlas.index) and data.n_vars==403
        np.testing.assert_allclose(data.obs[['x_centroid','y_centroid']],atlas[['x_centroid','y_centroid']],rtol=0,atol=1e-7)
        assert (data.X.data>=0).all()
        assert (data.layers['high_confidence']-data.X).max()<=0
        assert (data.layers['nuclear_only']-data.X).max()<=0
        for control,target in [('depth_matched','high_confidence'),('nuclear_depth_matched','nuclear_only')]:
            np.testing.assert_array_equal(np.asarray(data.layers[control].sum(1)).ravel(),np.asarray(data.layers[target].sum(1)).ravel())
        assert count['Q20_native_matrix_reconstructed_exactly']
        assert count['vendor_nuclear_ownership_preserved']==count['native_nuclear_molecules']
        audit=pd.read_csv(folder/'segmentation_profile_audit.tsv.gz',sep='\t',index_col=0)
        unresolved=audit.segmented_area_um2.eq(0)
        assert int(unresolved.sum())==seg['unseeded_profiles']
        assert int(unresolved.sum())==count['profiles_with_unresolved_2D_geometry_retaining_native_counts']
        if unresolved.any():
            original=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad',backed='r')
            ix=data.obs_names.get_indexer(audit.index[unresolved])
            np.testing.assert_array_equal(data.X[ix].toarray(),original.X[ix].toarray())
            original.file.close()
        with np.load(folder/'whole_section_segmentation.npz') as b:
            labels=b['cell_labels'];owner=b['nuclear_owner'];ids=b['cell_ids'].astype(str)
            np.testing.assert_array_equal(ids,atlas.index)
            # Check actual serialized nuclear ownership at randomly sampled pixels.
            rng=np.random.default_rng(192)
            y=rng.integers(0,labels.shape[0],200000);x=rng.integers(0,labels.shape[1],200000)
            valid=owner[y,x]>0
            assert valid.sum()>100
            np.testing.assert_array_equal(labels[y[valid],x[valid]],owner[y[valid],x[valid]])
        variants=pd.read_csv(folder/'fibroblast_variant_summary.tsv',sep='\t')
        for v in variants.variant.unique():
            f=pd.read_csv(folder/f'fibroblast_{v}_assignments.tsv.gz',sep='\t',index_col=0)
            assert len(f)==proj['broad_fibroblasts']
            selected=f.supported_family_assignment
            assert (f.loc[selected,'measured_state_markers']>=2).all()
            assert f.loc[selected,'both_library_agreement'].all() and f.loc[selected,'raw_reference_agreement'].all()
            mixtures.append({'stage':stage,'variant':v,'FB_profiles':len(f),
                'median_gene_molecules':f.total_gene_molecules.median(),
                'mixed_CM_FB_marker_review_profiles':int(f.mixed_CM_FB_marker_review.sum()),
                'mean_CM_marker_fraction':f.CM_marker_fraction.mean(),
                'median_CM_marker_fraction':f.CM_marker_fraction.median(),
                'informative_profiles':int(f.state_information_QC.sum()),
                'supported_fine_profiles':int(f.supported_family_assignment.sum())})
        for variant in list(variants.variant.unique())+['nuclear_modality_calibrated']:
            f=pd.read_csv(folder/f'fibroblast_{variant}_assignments.tsv.gz',sep='\t')
            eligible=f.state_information_QC
            gates={'weight_ge0.6':f.family_weight.ge(.6),'margin_ge0.15':f.family_margin.ge(.15),
                'at_least_two_state_markers':f.measured_state_markers.ge(2),
                'at_least_one_anchor':f.measured_state_anchors.ge(1),
                'both_RNA_libraries_agree':f.both_library_agreement,'raw_reference_agrees':f.raw_reference_agreement,
                'marker_plus_anchor':f.measured_state_markers.ge(2)&f.measured_state_anchors.ge(1)}
            for gate,passed in gates.items():
                gate_rows.append({'stage':stage,'variant':variant,'criterion':gate,
                    'informative_profiles':int(eligible.sum()),'informative_profiles_passing':int((eligible&passed).sum()),
                    'informative_profiles_failing':int((eligible&~passed).sum()),
                    'note':'Individual gate counts overlap; failures are not independent error causes'})
        evidence=pd.read_csv(folder/'fibroblast_nuclear_cell_evidence.tsv.gz',sep='\t')
        assert len(evidence)==proj['broad_fibroblasts']
        rows[-1]['conservative_recommended_supported']=int(evidence.recommended_cell_supported_family.ne('Unresolved fibroblast program').sum())
        concord=evidence.evidence_tier.eq('Cell and nuclear family concordant')
        evidence['previously_unresolved_now_cell_nuclear_concordant']=concord&evidence.native_family.eq('Unresolved fibroblast program')
        evidence['native_cell_nuclear_all_concordant']=concord&evidence.native_family.eq(evidence.nuclear_anchored_family)
        evidence['native_fine_family_change_review']=concord&evidence.native_family.ne('Unresolved fibroblast program')&evidence.native_family.ne(evidence.nuclear_anchored_family)
        rows[-1]['previously_unresolved_now_concordant']=int(evidence.previously_unresolved_now_cell_nuclear_concordant.sum())
        rows[-1]['native_cell_nuclear_all_concordant']=int(evidence.native_cell_nuclear_all_concordant.sum())
        rows[-1]['native_fine_family_change_review']=int(evidence.native_fine_family_change_review.sum())
        evidence.insert(0,'stage',stage);resolution_changes.append(evidence)
        checks.append({'stage':stage,'profiles':len(atlas),'nuclear_pixels_checked':int(valid.sum()),
                       'native_reconstruction_exact':True,'source_coordinates_unchanged':True,
                       'source_nuclear_RNA_owners_preserved':True,'matched_depths_exact':True,
                       'all_original_profiles_retained_in_count_object':True,
                       'unresolved_geometry_profiles_with_native_counts':int(unresolved.sum())})
    summary=pd.DataFrame(rows);summary.to_csv(OUT/'whole_section_summary.tsv',sep='\t',index=False)
    pd.concat(tierrows,ignore_index=True).to_csv(OUT/'fibroblast_evidence_tier_summary.tsv',sep='\t',index=False)
    pd.DataFrame(mixtures).to_csv(OUT/'fibroblast_assignment_information_summary.tsv',sep='\t',index=False)
    pd.DataFrame(gate_rows).to_csv(OUT/'fine_label_gate_coverage.tsv',sep='\t',index=False)
    pd.concat(resolution_changes,ignore_index=True).to_csv(OUT/'fibroblast_segmentation_resolution_changes.tsv.gz',sep='\t',index=False,na_rep='NA')
    write_json(OUT/'validation.json',{'status':'PASS','stages':checks,
        'total_original_profiles':int(summary.profiles.sum()),'original_records_streamed':int(summary.all_records.sum()),
        'counts_rebuilt_from_original_transcript_coordinates':True,
        'limitations':'Numerical and source-conservation validation, not independent anatomical segmentation accuracy'})
    make_report(summary)
    manifest=[]
    paths=[path for path in OUT.rglob('*') if path.is_file() and path.name!='file_manifest.tsv']
    paths += [ROOT/'analysis'/name for name in ['nuclear_anchored_Xenium_segmentation.py','project_nuclear_anchored_fibroblasts.py','report_nuclear_anchored_segmentation.py']]
    for path in sorted(paths):
        if path.is_file() and path.name!='file_manifest.tsv':
            with path.open('rb') as handle:digest=hashlib.file_digest(handle,'sha256').hexdigest()
            manifest.append({'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':digest})
    pd.DataFrame(manifest).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
    print(summary.to_string(index=False),flush=True)


def make_report(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(OUT/'nuclear_anchored_segmentation_report.pdf') as pdf:
        fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
        x=np.arange(len(summary))
        for offset,col,label in [(-.27,'native_supported','Native fine support'),(0,'nuclear_anchored_supported','Nuclear-anchored whole-cell support'),(.27,'nuclear_modality_calibrated_supported','Nuclear-only calibrated support')]:
            axes[0].bar(x+offset,summary[col],width=.27,label=label)
        axes[0].set_xticks(x,summary.stage);axes[0].set_ylabel('Broad FB profiles meeting fine-label criteria');axes[0].legend(fontsize=8)
        axes[1].bar(x,100*summary.nuclear_only_informative/summary.broad_fibroblasts)
        axes[1].set_xticks(x,summary.stage);axes[1].set_ylabel('FB profiles with informative nuclear RNA (%)')
        axes[1].set_ylim(0,100)
        fig.suptitle('Boundary refinement and nuclear RNA address different sources of uncertainty')
        fig.savefig(OUT/'all_stage_segmentation_and_nuclear_support.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        # A nested eligibility display separates ownership from fine-marker support.
        gates=pd.read_csv(OUT/'fine_label_gate_coverage.tsv',sep='\t')
        r=summary.set_index('stage').loc['P21']
        g=gates[gates.stage.eq('P21')&gates.variant.eq('nuclear_modality_calibrated')].set_index('criterion')
        names=['Broad fibroblasts','Informative nuclear RNA','Detected proposed-family anchor',
               '≥2 state markers including an anchor','All nuclear fine-label gates','Cell and nuclear family concordant']
        values=[r.broad_fibroblasts,r.nuclear_only_informative,g.loc['at_least_one_anchor','informative_profiles_passing'],
                g.loc['marker_plus_anchor','informative_profiles_passing'],r.nuclear_modality_calibrated_supported,r.cell_nuclear_concordant]
        fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True)
        ax.barh(names[::-1],values[::-1],color=['#267ea8']*5+['#999999'])
        for y,v in enumerate(values[::-1]):ax.text(v+80,y,f'{int(v):,}',va='center')
        ax.set_xlim(0,max(values)*1.13);ax.set_xlabel('P21 spatial fibroblast profiles')
        ax.set_title('P21: nuclear RNA ownership is often informative; fine-family anchors remain limiting')
        fig.savefig(OUT/'P21_nuclear_information_and_marker_gates.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        for stage in STAGES:
            for name in [f'{stage}_fibroblast_nuclear_anchored_maps.png',f'{stage}_assignment_information_comparison.png']:
                image=plt.imread(OUT/stage/name)
                fig,ax=plt.subplots(figsize=(15,7),constrained_layout=True);ax.imshow(image);ax.axis('off');pdf.savefig(fig);plt.close(fig)
    lines=['# Nuclear-anchored, antibody-guided whole-section Xenium segmentation','',
        '**Completed for P0, P7, P14 and P21:** nuclear-anchored membrane-guided refinement, molecule-level reassignment, '
        'whole-cell versus nuclear-only fibroblast mapping, paired RNA-depth controls and per-profile evidence tiers.', '',
        '[Illustrated report](nuclear_anchored_segmentation_report.pdf) · [Whole-section summary](whole_section_summary.tsv) · '
        '[Evidence tiers](fibroblast_evidence_tier_summary.tsv) · [Validation](validation.json)', '',
        '## What was done', '',
        'The deposited DAPI-derived nuclear outlines anchor a watershed on the acquired ATP1A1/CD45/E-cadherin '
        'boundary-stain image, at 0.425µm. Multiple nuclear components belonging to one source profile retain the same '
        'profile ID. Source profiles without a nucleus retain a seed in their native cell interior. The native polygon '
        'union expanded by 2µm limits the cellular territory. No RNA marker, subtype, Ankrd1 or maturation score '
        'optimizes segmentation. This refines boundaries for existing profiles; it is not independent de novo nuclear '
        'detection or discovery of missed cells/merge-split errors.', '',
        'Original Q20 gene molecules are reassigned from their actual coordinates; vendor-unassigned RNA remains '
        'unassigned in primary counts (an expanded-assignment sensitivity is supplied separately). Deposited nuclear-overlap assignments '
        'are retained rather than overturned by a coarser 2D projection. Projected nuclear ownership and disagreements '
        'are separate audit outputs. Nonnuclear RNA is assigned by the new image segmentation; projected nuclear '
        'conflicts remain unassigned unless the deposited nuclear-overlap flag supplies an owner. A high-confidence '
        'variant keeps nuclear RNA plus cytoplasmic RNA beyond a 0.85µm boundary band. Exact random depth controls '
        'match both this variant and nuclear-only RNA.', '',
        '## Results', '',
        '| Stage | All profiles | Nuclear anchors | Nucleus-free fallback anchors | Native supported FB labels | Nuclear-anchored supported | Nuclear RNA informative | Nuclear-only calibrated supported |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in summary.itertuples():
        lines.append(f'| {r.stage} | {r.profiles:,} | {r.nuclear_seeded_profiles:,} | {r.fallback_seeded_profiles:,} | '
                     f'{r.native_supported:,} | {r.nuclear_anchored_supported:,} | {r.nuclear_only_informative:,}/{r.broad_fibroblasts:,} | '
                     f'{r.nuclear_modality_calibrated_supported:,} |')
    lines += ['',
        f'Across stages, {int(summary.cell_nuclear_concordant.sum()):,} profiles have concordant supported cell/nuclear '
        f'family calls. Of these, **{int(summary.previously_unresolved_now_concordant.sum()):,} were previously unresolved** '
        '(15 P0, 114 P7, 60 P14 and 63 P21); 6,594 also agree with a supported native call. One P7 profile has '
        'concordant new cell/nuclear evidence for a family different from its native fine label and remains a '
        'cross-segmentation review case. These are additional evidence-qualified assignments, not validated anatomical truth.', '',
        '[Per-profile resolution changes](fibroblast_segmentation_resolution_changes.tsv.gz) identify the newly '
        'corroborated profiles and stable/conflicting native calls explicitly.', '',
        'All label counts use the preceding reviewed RNA hierarchy and unchanged marker/information/library gates. '
        'The native comparison exactly reproduces the preceding refit. Main variant comparisons freeze native-query '
        'gene scaling; a separately labeled nuclear-modality calibration is a sensitivity. Nuclear/cell agreement is '
        'corroborating evidence within the same profiles, not independent ground truth.', '',
        f'The unchanged fine-gate total changes from {int(summary.native_supported.sum()):,} to '
        f'{int(summary.nuclear_anchored_supported.sum()):,}. Segmentation therefore does not generally increase fine-label '
        'coverage in these data. Retain the native vendor atlas as the primary reference and use the new nuclear '
        'concordance/resolution-change flags to qualify confidence. The new count matrices and maps remain explicitly '
        'identified alternatives; these results do not justify a blanket replacement of the native boundaries.', '',
        '![All-stage summary](all_stage_segmentation_and_nuclear_support.png)', '',
        '## How to use the labels', '',
        '- **Cell and nuclear family concordant:** both assignments meet their criteria and agree; strongest available RNA/ownership corroboration, without anatomical validation.',
        '- **Cell-supported; nuclear confirmation unavailable/unresolved:** a whole-cell fine assignment remains supported, but nuclear RNA does not establish that fine identity.',
        '- **Nuclear-supported candidate; whole-cell strict criteria unmet:** useful follow-up candidate, kept separate from supported whole-cell calls.',
        '- **Nuclear-cell family disagreement or mixed CM–FB marker review:** do not force a fine family.',
        '- **Fine family unresolved:** retain the broad profile and uncertainty; no quota, smoothing or imputation creates a label.', '',
        '[Gate-coverage table](fine_label_gate_coverage.tsv) distinguishes inadequate overall RNA information from '
        'absent state markers/anchors and reference disagreements. At P21, 12,172 profiles have informative nuclear RNA, '
        'but only 1,730 of those have a detected anchor under the nuclear-calibrated proposal, and 1,494 have '
        'at least two supplied state markers including an anchor. Thus enough RNA for testing is not enough evidence for a named '
        'fine family. Gates overlap, so these are not independent causes of failure.', '',
        '![P21 nuclear information and marker gates](P21_nuclear_information_and_marker_gates.png)', '',
        'A visible nucleus helps establish ownership; it does not guarantee that the limited Xenium panel detects enough '
        'subtype markers in that nucleus. Likewise antibody segmentation does not restore unmeasured genes. Lower '
        'support after nucleus-only/core restriction must be assessed against the exact depth controls, rather than '
        'automatically interpreted as either contamination removal or failed segmentation.', '',
        '## Downloads', '',
        '| Stage | Fibroblast maps | Image gallery | Full 403-gene count object | Per-profile evidence |',
        '|---|---|---|---|---|']
    for stage in STAGES:
        lines.append(f'| {stage} | [PNG]({stage}/{stage}_fibroblast_nuclear_anchored_maps.png) | '
                     f'[8-page PDF]({stage}/{stage}_nuclear_segmentation_gallery.pdf) | '
                     f'[H5AD]({stage}/nuclear_anchored_counts.h5ad) | '
                     f'[TSV.gz]({stage}/fibroblast_nuclear_cell_evidence.tsv.gz) |')
    lines += ['', 'H5AD `X` is the nuclear-anchored whole-cell gene-count matrix. Layers contain `nuclear_only`, '
        '`high_confidence`, `depth_matched` and `nuclear_depth_matched` counts. All retain the source profile order and '
        '403 gene targets; `atlas_identity` in metadata is the original broad annotation, not a newly validated broad '
        'identity. New counts are supplied as an alternative analysis, not silently substituted into prior biological '
        'models. Existing Ankrd1/maturation spatial models have not been refitted to these new counts.', '',
        'The primary supported-label totals above count the unchanged fine gates before the additional mixed-marker '
        'and nuclear-disagreement review. The recommended cell-supported map also withholds those flagged profiles; '
        'its totals are included in `whole_section_summary.tsv`.', '',
        'One P0 endothelial profile has its entire original 2D polygon covered by another profile’s projected nuclear '
        'territory. No unique nuclear-compatible image region is forced for it: its source counts/profile are retained '
        'and geometry is explicitly flagged unresolved. No CM is removed for lacking a nucleus.', '',
        'The independent alternative uses the boundary channel and deposited nuclear outlines. The original vendor '
        'segmentation also had 18S and alphaSMA/Vimentin interior-stain inputs; this alternative is not automatically '
        'superior. Mixed-marker review rates and RNA-depth controls must be considered before replacing native '
        'boundaries. More supported labels are not an accuracy measure.', '',
        'Each stage also contains the whole-section integer segmentation/nuclear-ownership grids and physical pixel '
        'size, segmentation area/seed audit, transcript-assignment totals, all transfer variants and per-profile evidence. '
        'Nuclear-only projected counts are a separate sparse NPZ control, with the same cell/gene order as the H5AD. '
        'There is one deposited spatial section per age; cells are not independent animal replicates.', '',
        'Source: Wang, Dong, Song, Qian and Liu, deposited Xenium runs at Zenodo '
        '[14961929](https://zenodo.org/records/14961929) and [14961942](https://zenodo.org/records/14961942), CC BY4.0. '
        'The acquired stain crops and image derivatives are attributed to those authors.', '',
        'Reproduce using `analysis/acquire_xenium_boundary_inputs.py` for the boundary image and full transcript member; '
        'run `analysis/nuclear_anchored_Xenium_segmentation.py --stage STAGE --part segment` then `--part count`, '
        '`analysis/project_nuclear_anchored_fibroblasts.py --stage STAGE`, and finally '
        '`analysis/report_nuclear_anchored_segmentation.py`. Use the existing heart-atlas Python environment and '
        '`OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2`. The source DAPI/nuclear outlines, native '
        'cell polygons and reviewed RNA references must also be present. Completed stages resume from completion '
        'files; method changes require a clean derived result directory.', '',
        'Full downloaded image/transcript caches can be evicted only after their CRC provenance, SHA-256 and completed '
        'derived products are recorded. Source atlas/counts remain intact. Validation checks source reconstruction, '
        'nuclear ownership conservation, coordinates, sampled serialized nuclear pixels, count controls, label gates '
        'and profile accounting; it does not establish anatomical ground truth.', '']
    (OUT/'README.md').write_text('\n'.join(lines).rstrip()+'\n')


if __name__=='__main__':main()
