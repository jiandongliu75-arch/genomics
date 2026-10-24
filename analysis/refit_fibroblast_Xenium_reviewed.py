"""Refit all four Xenium sections with reviewed RNA labels; preserve prior results.

Parts are separate checkpoints: map, models, sensitivity, boundary, figures, validate.
No source count matrix, broad spatial identity, QC threshold or coordinate changes.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/fibroblast-refit-mpl')
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import sys
import numpy as np
import pandas as pd
import anndata as ad
from fibroblast_state_common import *

ORIGINAL = OUT
OUT = ROOT / 'results/fibroblast_Xenium_reviewed_refit'


def route(name):
    module = importlib.import_module(name)
    module.OUT = OUT
    return module


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(ORIGINAL/'RNA_cell_annotations.tsv.gz', sep='\t', index_col=0)
    ann['original_fibroblast_family'] = ann.fibroblast_family
    changed = ann.fibroblast_cluster_id.eq('P21_F0')
    assert changed.sum() == 1330 and ann.lineage_mixing_review_flag.sum() == 657
    ann.loc[changed, 'fibroblast_family'] = 'Col8a1-associated'
    ann['analysis_disposition'] = np.where(ann.fine_reference_eligible,
        'Retained in clean fibroblast analysis', 'Excluded from clean fibroblast analysis; retained for QC audit')
    ann.to_csv(OUT/'RNA_cell_annotations.tsv.gz', sep='\t', na_rep='NA')
    clean = ann.loc[ann.fine_reference_eligible]
    clean.to_csv(OUT/'RNA_clean_fibroblast_annotations.tsv.gz', sep='\t', na_rep='NA')
    ann.loc[~ann.fine_reference_eligible].to_csv(OUT/'RNA_excluded_lineage_review_nuclei.tsv.gz', sep='\t', na_rep='NA')
    summary = pd.read_csv(ORIGINAL/'RNA_cluster_summary.tsv', sep='\t')
    review = pd.read_csv(ORIGINAL/'RNA_residual_lineage_review.tsv', sep='\t')
    qc = summary.loc[summary.lineage_mixing_review_flag,
        ['stage','cluster_id','nuclei','top_markers','median_ambient_fraction','support']].merge(review, on=['stage','cluster_id'])
    save_tsv(qc, OUT/'excluded_group_QC_review.tsv')
    for stage in STAGES:
        (OUT/stage).mkdir(exist_ok=True)
    write_json(OUT/'refit_methods.json', {
        'original_results_preserved': str(ORIGINAL.relative_to(ROOT)),
        'RNA_change': 'Only the 1330 P21_F0 nuclei change broad family, unresolved -> Col8a1-associated',
        'flagged_groups': qc.cluster_id.tolist(), 'flagged_nuclei':657,
        'clean_fibroblast_nuclei':len(clean), 'source_fibroblast_nuclei':len(ann),
        'exclusion_scope':'Clean fibroblast figures/summaries and fine reference; source/QC records preserved; no automatic CM/mural reassignment',
        'mapping':'Full same-age weighted-30NN rerun: corrected, raw, both RNA libraries, and ECM-withheld reference',
        'fine_gate':'Unchanged: weight>=0.6, margin>=0.15, >=2 measured state markers, >=1 anchor, raw/both-library agreement, >=20 training molecules, >=8 training genes',
        'frozen':'Native gene counts and positions, broad atlas identities, query calibration scheme, outcome panels and model hyperparameters',
        'reference_component':'Postn where varying; otherwise Col8a1 when family weights exhaust a simplex (reviewed P21); unresolved weight is implicit reference where present',
        'limitations':'RNA-derived family interpretation; one spatial section per age; no anatomical or causal validation; no new segmentation'})


def mapping():
    prepare()
    prior_argv=sys.argv
    try:
        sys.argv=[prior_argv[0]]
        route('map_fibroblast_states_Xenium').main()
    finally:
        sys.argv=prior_argv
    rows=[];transitions=[];bench=[]
    for stage in STAGES:
        fb=pd.read_csv(OUT/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        old=pd.read_csv(ORIGINAL/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        assert fb.index.equals(old.index)
        comparison=fb[['family_proposal','mapped_family','family_weight','family_margin','both_library_agreement',
                       'raw_reference_agreement','measured_state_markers','state_information_QC']].add_prefix('reviewed_')
        comparison=comparison.join(old[['family_proposal','mapped_family','family_weight','family_margin']].add_prefix('original_'))
        comparison.to_csv(OUT/stage/'original_vs_reviewed_assignments.tsv.gz',sep='\t')
        for scope,col in [('proposal','family_proposal'),('supported_label','mapped_family')]:
            t=pd.crosstab(old[col],fb[col]).stack()
            for (a,b),n in t.items():
                if n:transitions.append({'stage':stage,'scope':scope,'original':a,'reviewed':b,'profiles':n})
        for fam in FAMILIES+[UNRESOLVED]:
            rows.append({'stage':stage,'family':fam,'broad_FB_profiles':len(fb),
                'original_proposals':int(old.family_proposal.eq(fam).sum()),'reviewed_proposals':int(fb.family_proposal.eq(fam).sum()),
                'original_supported_labels':int(old.mapped_family.eq(fam).sum()),'reviewed_supported_labels':int(fb.mapped_family.eq(fam).sum()),
                'reviewed_soft_weight_sum':fb['weight__'+fam].sum()})
        # Spatial neighbors are unchanged geometrically; reaggregate revised weights.
        context=pd.read_csv(ORIGINAL/stage/'FB_broad_neighbor_context_30um.tsv.gz',sep='\t',index_col=0)
        valid=context.resolved_other_neighbors.ge(5)&context.arterial_boundary_distance_um.notna()
        for fam in FAMILIES:
            for kind,w in [('soft_reference_weight',fb['weight__'+fam]),('supported_labels',fb.mapped_family.eq(fam).astype(float))]:
                w=w.where(valid,0)
                if w.sum()<1:continue
                for col in context.columns:
                    if col=='resolved_other_neighbors':continue
                    bench.append({'stage':stage,'family':fam,'summary':kind,'measure':col,
                        'sum_weights_or_supported_profiles':w.sum(),'weighted_mean':np.dot(context[col].fillna(0),w)/w.sum(),
                        'all_FB_mean':context.loc[valid,col].mean(),
                        'interpretation':'Reviewed-reference native-coordinate 30um neighborhood; soft weights are not calibrated identities'})
    save_tsv(pd.DataFrame(rows),OUT/'family_assignment_comparison.tsv')
    save_tsv(pd.DataFrame(transitions),OUT/'assignment_transition_counts.tsv')
    save_tsv(pd.DataFrame(bench),OUT/'FB_state_broad_spatial_context.tsv')


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.backends.backend_pdf import PdfPages
    handles=[Line2D([],[],marker='o',ls='',color=c,label=f) for f,c in COLORS.items()]
    with PdfPages(OUT/'reviewed_fibroblast_Xenium_maps.pdf') as pdf:
        ann=pd.read_csv(OUT/'RNA_clean_fibroblast_annotations.tsv.gz',sep='\t',index_col=0)
        fig,axes=plt.subplots(1,4,figsize=(19,5),constrained_layout=True)
        for ax,stage in zip(axes,STAGES):
            d=ann.loc[ann.age.eq(stage)]
            for fam,c in COLORS.items():
                a=d.loc[d.fibroblast_family.eq(fam)]
                ax.scatter(a.fibroblast_UMAP1,a.fibroblast_UMAP2,s=3,c=c,rasterized=True)
            for cid,g in d.groupby('fibroblast_cluster_id'):
                ax.text(g.fibroblast_UMAP1.median(),g.fibroblast_UMAP2.median(),cid.split('_')[-1],fontsize=9,
                        bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1))
            ax.set_title(f'{stage}: {len(d):,} retained nuclei');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
        fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
        fig.suptitle('Reviewed RNA families · residual CM/mural groups excluded · original UMAP coordinates')
        fig.savefig(OUT/'clean_fibroblast_RNA_UMAP.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        for stage in STAGES:
            fb=pd.read_csv(OUT/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
            old=pd.read_csv(ORIGINAL/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
            fig,axes=plt.subplots(1,3,figsize=(17,7),constrained_layout=True)
            for ax,d,col,title in [(axes[0],old,'mapped_family','Original supported families'),
                    (axes[1],fb,'mapped_family','Refitted supported families'),
                    (axes[2],fb,'family_proposal','Refitted proposals: include uncertain profiles')]:
                for fam,c in COLORS.items():
                    a=d.loc[d[col].eq(fam)]
                    ax.scatter(a.x_centroid,a.y_centroid,s=1.6,c=c,rasterized=True)
                ax.set_title(title);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
                ax.set_xlim(fb.x_centroid.min()-100,fb.x_centroid.max()+100)
                ax.set_ylim(fb.y_centroid.max()+100,fb.y_centroid.min()-100)
            fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
            fig.suptitle(f'{stage} · {fb.supported_family_assignment.sum():,}/{len(fb):,} profiles meet unchanged fine-label rules')
            fig.savefig(OUT/f'{stage}_reviewed_spatial_maps.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        stage='P21'
        fb=pd.read_csv(OUT/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        fig,axes=plt.subplots(1,4,figsize=(19,7),constrained_layout=True)
        for ax,fam in zip(axes,['Col8a1-associated','Lsamp-associated','Gfpt2-associated','Kirrel3-associated']):
            im=ax.scatter(fb.x_centroid,fb.y_centroid,c=fb['weight__'+fam],s=1.5,cmap='viridis',vmin=0,vmax=1,rasterized=True)
            ax.set_title(fam);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('X (µm)')
        fig.colorbar(im,ax=axes,shrink=.4,label='Reference similarity weight (not probability)')
        fig.suptitle('P21 · Reviewed family weights at original Xenium positions')
        fig.savefig(OUT/'P21_reviewed_family_weights.png',dpi=180);pdf.savefig(fig);plt.close(fig)


def validate():
    from scipy.spatial import cKDTree
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    source=pd.read_csv(ORIGINAL/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    assert ann.index.equals(source.index) and len(ann)==14985
    delta=ann.fibroblast_family.ne(source.fibroblast_family)
    assert delta.sum()==1330 and ann.loc[delta,'fibroblast_cluster_id'].eq('P21_F0').all()
    assert ann.fine_reference_eligible.equals(source.fine_reference_eligible)
    checks=[];rng=np.random.default_rng(620)
    for stage in STAGES:
        fb=pd.read_csv(OUT/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        old=pd.read_csv(ORIGINAL/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
        pd.testing.assert_frame_equal(fb[['x_centroid','y_centroid']],old[['x_centroid','y_centroid']])
        w=fb[['weight__'+f for f in FAMILIES+[UNRESOLVED]]].to_numpy()
        assert np.allclose(w.sum(1),1,atol=2e-6)
        with np.load(OUT/stage/'fibroblast_projection_calibration.npz') as b:
            ids=b['reference_ids'].astype(str)
            assert not ann.loc[ids,'lineage_mixing_review_flag'].any()
            assert not set(b['training_genes'].astype(str))&set(WITHHOLD)
        if stage!='P21':
            assert fb.family_proposal.equals(old.family_proposal)
            assert fb.mapped_family.equals(old.mapped_family)
            np.testing.assert_allclose(w,old[['weight__'+f for f in FAMILIES+[UNRESOLVED]]],atol=2e-6)
        else:
            # Since only labels change, exact kNN coarsening must match a complete rerun.
            expected=old['weight__Col8a1-associated']+old['weight__'+UNRESOLVED]
            np.testing.assert_allclose(fb['weight__Col8a1-associated'],expected,atol=2e-6)
            assert (fb['weight__'+UNRESOLVED].abs()<1e-6).all()
        supported=fb.supported_family_assignment
        assert (fb.loc[supported,'measured_state_markers']>=2).all()
        assert fb.loc[supported,'both_library_agreement'].all() and fb.loc[supported,'raw_reference_agreement'].all()
        strict=(fb.family_weight.ge(.6)&fb.family_margin.ge(.15)&fb.both_library_agreement&fb.raw_reference_agreement
                &fb.measured_state_markers.ge(2)&fb.measured_state_anchors.ge(1)&fb.state_information_QC
                &fb.family_proposal.ne(UNRESOLVED))
        assert supported.equals(strict)
        if stage=='P21':
            data=ad.read_h5ad(ROOT/'results/xenium_ventricular_mapping'/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad',backed='r')
            counts=data.X[data.obs_names.get_indexer(fb.index)].toarray()
            genes=data.var_names.to_numpy(dtype=str)
            for fam,panel in PANELS.items():
                proposed=fb.family_proposal.eq(fam).to_numpy()
                observed=(counts[proposed][:,np.isin(genes,panel)]>0).sum(axis=1)
                np.testing.assert_array_equal(observed,fb.loc[proposed,'measured_state_markers'])
            data.file.close()
        cm=pd.read_csv(OUT/stage/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        oc=pd.read_csv(ORIGINAL/stage/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        pd.testing.assert_frame_equal(cm,oc,check_exact=False,rtol=2e-6,atol=2e-6)
        tree=cKDTree(fb[['x_centroid','y_centroid']].to_numpy())
        for radius in [30,50,100]:
            env=pd.read_csv(OUT/stage/f'CM_neighborhoods_{radius}um.tsv.gz',sep='\t',index_col=0)
            for cid in rng.choice(cm.index,16,replace=False):
                ix=tree.query_ball_point(cm.loc[cid,['x_centroid','y_centroid']].to_numpy(dtype=float),radius)
                assert len(ix)==env.loc[cid,'fibroblast_neighbors']
                for fam in FAMILIES:
                    expected=fb.iloc[ix]['weight__'+fam].mean() if ix else 0
                    assert np.isclose(env.loc[cid,'family__'+fam],expected,atol=2e-6)
        checks.append({'stage':stage,'broad_fibroblasts':len(fb),'supported_fine':int(supported.sum()),
            'old_supported_fine':int(old.supported_family_assignment.sum()),'reference_nuclei':len(ids),
            'changed_proposals':int(fb.family_proposal.ne(old.family_proposal).sum()),
            'changed_supported_labels':int(fb.mapped_family.ne(old.mapped_family).sum())})
        confidence=pd.read_csv(OUT/stage/'Xenium_boundary_confidence.tsv.gz',sep='\t',index_col=0)
        assert confidence.index.equals(fb.index)
        for name in ['core_1um','depth_matched','membrane']:
            assert confidence.loc[~confidence[name+'_tested'],name+'_label_retained'].isna().all()
    for filename in ['spatial_blocked_prediction.tsv','spatial_adjusted_coefficients.tsv',
            'spatial_model_sensitivity_prediction.tsv','boundary_fibroblast_state_retention.tsv',
            'boundary_spatial_association_sensitivity.tsv','Ankrd1_vascular_pairs_after_FB_adjustment.tsv']:
        assert (OUT/filename).exists(),filename
    write_json(OUT/'validation.json',{'status':'PASS','source_nuclei':len(ann),'clean_nuclei':int(ann.fine_reference_eligible.sum()),
        'excluded_review_nuclei':int((~ann.fine_reference_eligible).sum()),'stages':checks,
        'independent_physical_neighborhood_checks':192,'all_source_coordinates_unchanged':True,
        'CM_outcomes_unchanged':True,'all_flagged_nuclei_excluded_from_reference':True,
        'P0_P7_P14_assignments_reproduced':True,'P21_full_refit_matches_expected_knn_coarsening':True})
    support=pd.DataFrame(checks).rename(columns={'broad_fibroblasts':'broad_FBs','supported_fine':'supported_fine_profiles'})
    support['fraction_supported']=support.supported_fine_profiles/support.broad_FBs
    save_tsv(support,OUT/'Xenium_fine_state_support_summary.tsv')
    manifest=[]
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and path.name!='file_manifest.tsv':
            manifest.append({'file':str(path.relative_to(OUT)),'bytes':path.stat().st_size,
                             'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    save_tsv(pd.DataFrame(manifest),OUT/'file_manifest.tsv')
    print(json.dumps(checks),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--part',required=True,
        choices=['map','models','sensitivity','boundary','figures','validate']);args=parser.parse_args()
    if args.part=='map':mapping()
    elif args.part=='models':route('compare_fibroblast_matrix_and_spatial').spatial_comparisons()
    elif args.part=='sensitivity':
        route('check_fibroblast_spatial_model_sensitivity').main()
        route('check_Ankrd1_vessel_after_fibroblast_adjustment').main()
    elif args.part=='boundary':
        route('check_fibroblast_boundary_sensitivity').main(save_profile_assignments=True)
        methods=json.loads((OUT/'boundary_methods.json').read_text())
        methods['fine_reference']='Frozen reviewed stage reference and native-query calibration; broad parent labels retained'
        write_json(OUT/'boundary_methods.json',methods)
        audit=pd.read_csv(OUT/'boundary_profile_assignments.tsv.gz',sep='\t')
        for stage in STAGES:
            fb=pd.read_csv(OUT/stage/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
            confidence=fb[['mapped_family','supported_family_assignment']].copy()
            for scope,variant,name in [('geometry','core_1um','core_1um'),('geometry','thin_to_core1','depth_matched'),('image','membrane_watershed','membrane')]:
                tested=audit.loc[audit.stage.eq(stage)&audit.scope.eq(scope)&audit.variant.eq(variant)].set_index('cell_id')
                confidence[name+'_tested']=confidence.index.isin(tested.index)
                retained=tested.mapped_family.eq(tested.native_mapped_family)
                confidence[name+'_label_retained']=retained.reindex(confidence.index).astype('boolean')
            confidence.to_csv(OUT/stage/'Xenium_boundary_confidence.tsv.gz',sep='\t',na_rep='NA')
    elif args.part=='figures':figures()
    elif args.part=='validate':validate()
