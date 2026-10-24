"""Conservative fibroblast-family annotation from existing segmentation evidence.

No new model is fitted. Agreement is a robustness criterion, not calibrated
accuracy: the views share transcripts and RNA references. Keep every broad FB
profile, original proposals and explicit reasons for withholding a fine label.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-consensus-mpl')
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
from fibroblast_state_common import ROOT,STAGES,FAMILIES,UNRESOLVED,PANELS,ANCHORS,COLORS

SOURCE=ROOT/'results/xenium_nuclear_anchored_segmentation'
OUT=ROOT/'results/fibroblast_spatial_consensus'
VIEWS=['native','nuclear_anchored','high_confidence','nuclear_only','nuclear_modality_calibrated','cytoplasmic_only']
TIERS=['A: stable cell and nuclear','B: new cell and nuclear','C: stable cell only',
       'Review: mixed lineage markers','Review: supported families conflict','Unresolved: insufficient corroboration']
TIER_COLORS=dict(zip(TIERS,['#0072b2','#56b4e9','#009e73','#d55e00','#cc79a7','#bbbbbb']))


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')


def classify(frames):
    primary=frames['nuclear_anchored'];native=frames['native'];core=frames['high_confidence']
    nuc=frames['nuclear_only'];cal=frames['nuclear_modality_calibrated']
    supported=np.column_stack([np.where(frames[v].supported_family_assignment,
                            frames[v].mapped_family,UNRESOLVED) for v in VIEWS])
    n_families=np.column_stack([(supported==fam).any(1) for fam in FAMILIES]).sum(1)
    conflict=n_families>1
    mixed=np.logical_or.reduce([frames[v].mixed_CM_FB_marker_review.to_numpy()
                              for v in ['nuclear_anchored','nuclear_only','cytoplasmic_only']])
    nuclear_agreement=(nuc.supported_family_assignment&cal.supported_family_assignment
                       &nuc.mapped_family.eq(cal.mapped_family))
    nuclear_cell=nuclear_agreement&primary.supported_family_assignment&primary.mapped_family.eq(cal.mapped_family)
    no_review=~mixed&~conflict
    stable=nuclear_cell&native.supported_family_assignment&native.mapped_family.eq(primary.mapped_family)&no_review
    new=nuclear_cell&~native.supported_family_assignment&no_review
    cell=(native.supported_family_assignment&primary.supported_family_assignment&core.supported_family_assignment
          &native.mapped_family.eq(primary.mapped_family)&core.mapped_family.eq(primary.mapped_family)&no_review&~stable&~new)
    d=pd.DataFrame(index=primary.index)
    d.index.name='cell_id'
    d['confidence_tier']=TIERS[-1]
    d.loc[conflict,'confidence_tier']=TIERS[4]
    d.loc[mixed,'confidence_tier']=TIERS[3]
    d.loc[cell,'confidence_tier']=TIERS[2]
    d.loc[new,'confidence_tier']=TIERS[1]
    d.loc[stable,'confidence_tier']=TIERS[0]
    assigned=stable|new|cell
    d['assigned_family']=np.where(assigned,primary.mapped_family,UNRESOLVED)
    d['high_confidence_family']=np.where(stable|new,primary.mapped_family,UNRESOLVED)
    d['supported_family_conflict']=conflict
    d['mixed_lineage_compartment_review']=mixed
    d['number_distinct_supported_families']=n_families
    d['both_nuclear_calibrations_agree_supported']=nuclear_agreement
    d['nuclear_RNA_informative']=nuc.state_information_QC
    for v in VIEWS:
        f=frames[v]
        for col in ['mapped_family','family_proposal','supported_family_assignment','family_weight','family_margin',
                    'measured_state_markers','measured_state_anchors','state_information_QC']:
            d[v+'__'+col]=f[col]
    return d


def main():
    OUT.mkdir(exist_ok=True)
    write_json(OUT/'annotation_rules.json',{
        'source':'Existing all-stage nuclear-anchored segmentation and reviewed RNA-family transfer; no retraining',
        'A':'Native and resegmented cell calls supported and identical; frozen-native and nuclear-modality-calibrated nuclear calls both supported and identical to the cell; no strong conflicting view or compartment mixed-marker flag',
        'B':'Previously native-unresolved; resegmented cell and both nuclear calibration calls supported and identical; no strong conflicting view or compartment mixed-marker flag',
        'C':'Native, resegmented and boundary-core cell calls supported and identical; A/B criteria unmet; no strong conflicting view or compartment mixed-marker flag',
        'conflict':'More than one family receives a supported call among native, resegmented, boundary-core, both nuclear calibrations or cytoplasmic-only views',
        'mixed':'Existing CM/FB mixed-marker review flag in resegmented whole-cell, nuclear-only or cytoplasmic-only counts; a review criterion, not proven contamination or doublets',
        'support_gates':'Inherited unchanged: weight>=0.6,margin>=0.15,>=2 measured supplied state markers including>=1 anchor,raw/both-library agreement,>=20 training molecules and>=8 training genes',
        'primary_high_confidence_map':'A+B; C supplied as a separate moderate-confidence extension',
        'no_nucleus_area_filter':True,'no_spatial_smoothing_or_label_quota':True,
        'no_cross_stage_label_borrowing':True,'no_random_thinning_votes':True,
        'interpretation':'Operational robustness tiers, not independent accuracy estimates, posterior probabilities or anatomical subtype validation'})
    allrows=[];summaries=[];markers=[];validations=[]
    for stage in STAGES:
        folder=OUT/stage;folder.mkdir(exist_ok=True)
        frames={v:pd.read_csv(SOURCE/stage/f'fibroblast_{v}_assignments.tsv.gz',sep='\t',index_col=0) for v in VIEWS}
        for v in VIEWS:assert frames[v].index.equals(frames['native'].index)
        evidence=pd.read_csv(SOURCE/stage/'fibroblast_nuclear_cell_evidence.tsv.gz',sep='\t',index_col=0)
        d=evidence[['x_centroid','y_centroid','nucleus_area']].join(classify(frames))
        d.insert(0,'stage',stage)
        d.insert(1,'broad_identity','Fibroblast')
        count=ad.read_h5ad(SOURCE/stage/'nuclear_anchored_counts.h5ad',backed='r')
        ix=count.obs_names.get_indexer(d.index);assert (ix>=0).all()
        np.testing.assert_allclose(d[['x_centroid','y_centroid']].to_numpy(),
            count.obs.iloc[ix][['x_centroid','y_centroid']].to_numpy(),rtol=0,atol=1e-7)
        genes=count.var_names.to_numpy(dtype=str)
        cell=count.X[ix].toarray()
        nuclear=count.layers['nuclear_only'][ix].toarray()
        d['cell_state_marker_genes']='';d['nuclear_state_marker_genes']=''
        d['cell_anchor_genes']='';d['nuclear_anchor_genes']=''
        for family in FAMILIES:
            selected=d.assigned_family.eq(family).to_numpy()
            indices=np.flatnonzero(selected)
            for panel,col in [(PANELS[family],'state_marker_genes'),(ANCHORS[family],'anchor_genes')]:
                gi=np.flatnonzero(np.isin(genes,panel))
                for values,prefix in [(cell,'cell'),(nuclear,'nuclear')]:
                    lists=[';'.join(genes[gi][values[j,gi]>0]) for j in indices]
                    d.loc[selected,prefix+'_'+col]=lists
            markers.append({'stage':stage,'family':family,'measured_state_markers':';'.join(g for g in PANELS[family] if g in genes),
                'measured_anchors':';'.join(g for g in ANCHORS[family] if g in genes),
                'high_confidence_profiles':int(d.high_confidence_family.eq(family).sum()),
                'all_supported_profiles':int(selected.sum())})
        count.file.close()
        d.to_csv(folder/'fibroblast_consensus_annotations.tsv.gz',sep='\t',na_rep='NA')
        allrows.append(d.reset_index())
        for tier in TIERS:
            for family in FAMILIES+[UNRESOLVED]:
                n=int((d.confidence_tier.eq(tier)&d.assigned_family.eq(family)).sum())
                if n:summaries.append({'stage':stage,'confidence_tier':tier,'family':family,'profiles':n})
        # Independently validate the evidence required for every released label.
        high=d.high_confidence_family.ne(UNRESOLVED)
        assigned=d.assigned_family.ne(UNRESOLVED)
        assert not d.loc[assigned,'supported_family_conflict'].any()
        assert not d.loc[assigned,'mixed_lineage_compartment_review'].any()
        assert d.loc[high,'both_nuclear_calibrations_agree_supported'].all()
        for v in ['nuclear_anchored','nuclear_only','nuclear_modality_calibrated']:
            assert frames[v].loc[high,'supported_family_assignment'].all()
            assert frames[v].loc[high,'mapped_family'].equals(d.loc[high,'assigned_family'])
        for cid in d.index[assigned]:
            family=d.loc[cid,'assigned_family'];j=d.index.get_loc(cid)
            assert (cell[j,np.isin(genes,PANELS[family])]>0).sum()>=2
            assert (cell[j,np.isin(genes,ANCHORS[family])]>0).sum()>=1
            if high.loc[cid]:
                assert (nuclear[j,np.isin(genes,PANELS[family])]>0).sum()>=2
                assert (nuclear[j,np.isin(genes,ANCHORS[family])]>0).sum()>=1
        c=d.confidence_tier.eq(TIERS[2])
        for v in ['native','nuclear_anchored','high_confidence']:
            assert frames[v].loc[c,'supported_family_assignment'].all()
            assert frames[v].loc[c,'mapped_family'].equals(d.loc[c,'assigned_family'])
        b=d.confidence_tier.eq(TIERS[1]);assert not frames['native'].loc[b,'supported_family_assignment'].any()
        assert d.loc[d.assigned_family.eq('Kirrel3-associated')].empty
        validations.append({'stage':stage,'broad_FB_profiles':len(d),'high_confidence':int(high.sum()),
            'new_high_confidence':int(b.sum()),'moderate_cell_supported':int(c.sum()),
            'all_supported':int(assigned.sum()),'unresolved_or_review':int((~assigned).sum()),
            'A':int(d.confidence_tier.eq(TIERS[0]).sum()),
            'mixed_review':int(d.confidence_tier.eq(TIERS[3]).sum()),
            'conflict_review':int(d.confidence_tier.eq(TIERS[4]).sum())})
        print('CONSENSUS_DONE',stage,validations[-1],flush=True)
    merged=pd.concat(allrows,ignore_index=True)
    assert len(merged)==33790 and not merged.duplicated(['stage','cell_id']).any()
    merged.to_csv(OUT/'all_stage_fibroblast_consensus.tsv.gz',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(summaries).to_csv(OUT/'family_by_tier_counts.tsv',sep='\t',index=False)
    pd.DataFrame(markers).to_csv(OUT/'identity_markers_and_assignment_counts.tsv',sep='\t',index=False)
    pd.DataFrame(validations).to_csv(OUT/'stage_confidence_summary.tsv',sep='\t',index=False)
    write_json(OUT/'validation.json',{'status':'PASS','total_profiles':len(merged),'stages':validations,
        'all_source_FB_profiles_retained':True,'all_released_calls_checked_against_required_views':True,
        'all_released_calls_have_measured_cell_marker_evidence':True,
        'scope':'Numerical/criteria verification, not independent anatomical ground truth'})
    figures(merged,pd.DataFrame(validations))
    readme(pd.DataFrame(validations))
    paths=list(OUT.rglob('*'))+[Path(__file__)]
    rows=[]
    for p in sorted(paths):
        if p.is_file() and p.name!='file_manifest.tsv':
            with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            rows.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest})
    pd.DataFrame(rows).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)


def figures(d,summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(OUT/'fibroblast_spatial_confidence_report.pdf') as pdf:
        for stage in STAGES:
            a=d.loc[d.stage.eq(stage)]
            fig,axes=plt.subplots(1,3,figsize=(18,7),constrained_layout=True)
            for ax,col,title in [(axes[0],'high_confidence_family','A+B: nuclear and cell support'),
                                  (axes[1],'assigned_family','A+B+C: including stable cell-only calls')]:
                for fam,color in COLORS.items():
                    b=a.loc[a[col].eq(fam)]
                    ax.scatter(b.x_centroid,b.y_centroid,c=color,s=1.5,rasterized=True)
                ax.set_title(title)
            for tier,color in TIER_COLORS.items():
                b=a.loc[a.confidence_tier.eq(tier)]
                axes[2].scatter(b.x_centroid,b.y_centroid,c=color,s=1.5,label=tier,rasterized=True)
            axes[2].set_title('Evidence tier; all broad FB profiles retained')
            axes[2].legend(loc='upper left',bbox_to_anchor=(1,1),fontsize=8,frameon=False,markerscale=4)
            for ax in axes:
                ax.set_xlim(a.x_centroid.min()-100,a.x_centroid.max()+100)
                ax.set_ylim(a.y_centroid.max()+100,a.y_centroid.min()-100)
                ax.set_aspect('equal');ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
            fig.legend(handles=[Line2D([],[],ls='',marker='o',color=c,label=f) for f,c in COLORS.items()],
                       loc='outside lower center',ncol=3,frameon=False)
            fig.suptitle(stage+' · Conservative integrated fibroblast-family assignment · original positions')
            fig.savefig(OUT/f'{stage}_fibroblast_confidence_map.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
        ax.bar(summary.stage,summary.high_confidence,label='A+B: nuclear and cell',color='#0072b2')
        ax.bar(summary.stage,summary.moderate_cell_supported,bottom=summary.high_confidence,label='C: stable cell only',color='#009e73')
        ax.bar(summary.stage,summary.unresolved_or_review,bottom=summary.all_supported,label='Unresolved/review',color='#bbbbbb')
        ax.set_ylabel('Broad spatial fibroblast profiles');ax.legend(frameon=False)
        ax.set_title('Confidence-qualified assignments; no forced labels')
        fig.savefig(OUT/'all_stage_confidence_summary.png',dpi=180);pdf.savefig(fig);plt.close(fig)


def readme(summary):
    lines=['# Integrated fibroblast spatial assignments with explicit confidence','',
        'This consolidates the reviewed same-age RNA references, native Xenium labels, whole-section nuclear-anchored '
        'resegmentation, boundary-core RNA, nuclear-only RNA under two calibration schemes, and cytoplasmic-only RNA. '
        'It fits no new classifier and changes no count matrix, coordinate or broad parent identity.', '',
        '[Illustrated report](fibroblast_spatial_confidence_report.pdf) · [All-profile annotations](all_stage_fibroblast_consensus.tsv.gz) · '
        '[Counts by family/tier](family_by_tier_counts.tsv) · [Rules](annotation_rules.json) · [Validation](validation.json)', '',
        '## Assignment rules', '',
        '- **A: stable cell and nuclear.** Native and resegmented cell assignments agree; both nuclear calibration schemes each pass the inherited gates for that same family.',
        '- **B: new cell and nuclear.** Native fine identity was unresolved; resegmented cell and both nuclear calibration schemes pass for the same family.',
        '- **C: stable cell only.** Native, resegmented and boundary-core cell calls all pass for the same family; the stricter A/B nuclear criterion is unmet.',
        '- **Review/unresolved.** A strong supported family conflict, CM–FB mixed-marker compartment flag, or inadequate corroboration prevents an integrated fine call.', '',
        'A/B/C all require no conflicting supported family in any of six tested views and no existing mixed CM–FB '
        'marker flag in resegmented whole-cell, nuclear-only or cytoplasmic-only RNA. A mixed-marker flag is a '
        'conservative review criterion, not proof of contamination or a doublet. Random depth controls do not vote '
        'on identity. There is no nuclear-area requirement, cross-stage label borrowing, spatial smoothing or quota.', '',
        'Each source supported call already requires weight≥0.6, margin≥0.15, at least two measured supplied state '
        'markers including an anchor, original-count and both-RNA-library agreement, and sufficient training-panel '
        'RNA (≥20 molecules and≥8 detected genes). These thresholds are unchanged. In A/B, actual measured marker '
        'and anchor genes are additionally checked in nuclear counts and published per profile.', '',
        '**The recommended primary confidence subset is A+B.** C is a separate moderate-confidence extension. '
        'These are operational robustness tiers, not calibrated probabilities, independent biological accuracy '
        'or anatomical validation: the views share RNA molecules and references. Both nuclear schemes agreeing '
        'tests sensitivity to gene scaling; it is not two independent experiments.', '',
        '## Results','',
        '| Stage | Broad FB profiles | A+B nuclear/cell-supported | B newly supported | C stable cell-only | Unresolved/review |',
        '|---|---:|---:|---:|---:|---:|']
    for r in summary.itertuples():
        lines.append(f'| {r.stage} | {r.broad_FB_profiles:,} | {r.high_confidence:,} | {r.new_high_confidence:,} | {r.moderate_cell_supported:,} | {r.unresolved_or_review:,} |')
    lines += ['',f'Total: **{int(summary.high_confidence.sum()):,} A+B**, including **{int(summary.new_high_confidence.sum()):,} B** '
        f'previously unresolved profiles, plus **{int(summary.moderate_cell_supported.sum()):,} C** cell-supported profiles. '
        'All 33,790 broad fibroblast profiles remain accounted for. A smaller released label set reflects the '
        'additional corroboration requirements; it is not a disappearance of fibroblasts.', '',
        'The preceding segmentation report found 252 previously unresolved profiles with cell/nuclear-calibrated '
        'concordance. The current B rule additionally requires nuclear frozen-calibration support, absence of '
        'cytoplasmic/core supported-family conflict, and absence of compartment mixed-marker flags. Its count '
        'therefore need not equal 252. Previously supported single-view labels remain inspectable as proposals '
        'but are not automatically promoted into this integrated annotation.', '',
        '![Confidence totals](all_stage_confidence_summary.png)', '',
        '## Family identities and limits', '',
        'Named families remain Postn-associated, Col8a1-related, Gfpt2-associated, Lsamp-associated and Kirrel3-associated '
        'RNA programs. This assigns broad fibroblast families; it does not distinguish all within-family RNA '
        'subclusters spatially, including separate P21_F0/F2 or F1/F5 locations. No Kirrel3 fine call passes the '
        'inherited two-marker gate because only one supplied valve-family marker is measured. No valve/outflow '
        'location is inferred from Kirrel3/Eln alone.', '',
        '[Measured identity markers and assigned counts](identity_markers_and_assignment_counts.tsv) documents '
        'the assay-covered panel. Per-profile `cell_state_marker_genes`, `nuclear_state_marker_genes` and anchor '
        'columns show actual detected genes for released families. These are criteria/evidence audits, not '
        'independent differential-expression discoveries or function validation.', '',
        '## Maps and data', '',
        '| Stage | Family/confidence map | Full annotations |', '|---|---|---|']
    for stage in STAGES:
        lines.append(f'| {stage} | [PNG]({stage}_fibroblast_confidence_map.png) | [TSV.gz]({stage}/fibroblast_consensus_annotations.tsv.gz) |')
    lines += ['', 'Left panels show A+B; middle panels add C; right panels show the confidence/review tier. '
        'Gray means withheld fine identity. All use original µm positions and preserve broad fibroblast parents. '
        'Use stage+cell_id as the global key. Existing native and alternative gene-count H5ADs remain in the '
        '[preceding segmentation analysis](../xenium_nuclear_anchored_segmentation/README.md). This confidence '
        'catalog does not refit previous Ankrd1/maturation association models.', '',
        'Reproduce with `OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 '
        '/workspace/heart-atlas/venv/bin/python analysis/consolidate_fibroblast_spatial_confidence.py`. '
        'Validation checks all source profile accounting/coordinates, every required supported view, absence of '
        'strong competing families, and observed cell/nuclear marker counts for every released label. '
        '[Artifact/script hashes](file_manifest.tsv) preserve provenance.', '']
    (OUT/'README.md').write_text('\n'.join(lines).rstrip()+'\n')


if __name__=='__main__':main()
