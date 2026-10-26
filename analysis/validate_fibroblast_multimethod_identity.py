"""Independent numerical checks of model inference, evidence and source IDs."""
import json
import numpy as np
import pandas as pd
from scipy.special import softmax
from scipy.stats import nbinom
from sklearn.metrics import precision_recall_fscore_support
from fibroblast_multimethod_identity import ROOT,OUT,STAGES,METHODS,FAMILIES,NAMING,EXCLUDE,refs,queries,Classifiers,read

def main():
    ann,x,raw,genes=refs();cv=read(OUT/'RNA_opposite_library_predictions.tsv.gz');bench=read(OUT/'RNA_opposite_library_model_benchmark.tsv')
    folds=read(OUT/'RNA_fold_training_genes.tsv');withheld=folds[folds.variant.eq('naming_markers_withheld')]
    assert not set(withheld.gene)&NAMING and not set(folds.gene)&EXCLUDE
    for (stage,variant,method,library),g in cv.groupby(['stage','variant','method','heldout_library']):
        assert not g.cell_id.duplicated().any()
        expected=ann.index[ann.age.eq(stage)&ann['sample'].eq(library)]
        assert set(g.cell_id)==set(expected)
        assert ann.loc[g.cell_id,'sample'].eq(library).all()
        assert np.array_equal(ann.loc[g.cell_id,'fibroblast_family'],g.true_family)
        for family in g.true_family.unique():
            tp=((g.prediction==family)&(g.true_family==family)).sum();fp=((g.prediction==family)&(g.true_family!=family)).sum();fn=((g.prediction!=family)&(g.true_family==family)).sum()
            pr=tp/max(tp+fp,1);re=tp/max(tp+fn,1);f1=2*pr*re/max(pr+re,1e-10)
            row=bench[bench.stage.eq(stage)&bench.variant.eq(variant)&bench.method.eq(method)&bench.heldout_library.eq(library)&bench.condition.eq('observed')&bench.family.eq(family)].iloc[0]
            np.testing.assert_allclose([pr,re,f1],[row.precision,row.recall,row.F1],atol=1e-7)
    rows=[];reconstructed=0
    for stage in STAGES:
        d=read(OUT/(stage+'_multimethod_FB_annotations.tsv.gz'),index_col='cell_id')
        atlas=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',index_col=0)
        assert d.index.equals(atlas.index[atlas.atlas_identity.eq('Fibroblast')])
        np.testing.assert_allclose(d[['x_centroid','y_centroid']],atlas.loc[d.index,['x_centroid','y_centroid']],rtol=0,atol=1e-9)
        with np.load(OUT/(stage+'_model_predictions.npz')) as f:
            assert np.array_equal(f['cell_ids'],d.index)
            classes=f['families'];variants=f['variants'];views=f['views'];methods=f['methods']
            assert np.array_equal(methods,METHODS)
            for v,variant in enumerate(variants):
                for w,view in enumerate(views):
                    codes=f['label_code'][v,w];support=f['reference_screen_pass'][v,w]
                    for m,method in enumerate(methods):
                        expected=classes[codes[m]];col=variant+'__'+view+'__'+method
                        np.testing.assert_array_equal(expected,d[col]);reconstructed+=len(d)
                    for screened in [False,True]:
                        result=np.full(len(d),'Unresolved',object)
                        for ci,c in enumerate(classes):
                            mask=(codes==ci)&(support if screened else True);take=mask.sum(0)>=2;result[take]=c
                        target=variant+'__'+view+('__reference_supported' if screened else '__majority')
                        np.testing.assert_array_equal(result,d[target])
        selected=d.recommended_candidate.ne('Unresolved')
        assert d.loc[selected,'program_source_eligible'].all() and not d.loc[selected,'core_counterevidence'].any()
        strict=d.crosschecked_candidate.ne('Unresolved')
        assert d.loc[strict,'recommended_candidate'].eq(d.loc[strict,'crosschecked_candidate']).all()
        assert d.loc[strict,'evidence_tier'].eq('Cell and core supported').all()
        assert d.loc[strict,'naming_withheld_same_candidate'].all()
        module=read(OUT/(stage+'_refined_published_module_scores.tsv.gz'),index_col='cell_id')
        assert module.index.equals(d.index)
        overlap=read(OUT/'published_Gfpt2_Col8a1_program_overlap.tsv');g=overlap[overlap.stage.eq(stage)]
        for scope,mask in [('all_eligible',d.program_source_eligible),('public_Gfpt2',d.program_source_eligible&d.author_atlas_cell_state.eq('Gfpt2.Fib'))]:
            z=g[g.scope.eq(scope)];assert z.profiles.sum()==mask.sum() and z.scope_profiles.eq(mask.sum()).all()
            gf=module['Gfpt2.Fib__enrichment_flag'];co=module['Col8a1.Fib__enrichment_flag']
            for name,flag in [('Gfpt2_only',gf&~co),('Col8a1_only',co&~gf),('both',gf&co),('neither',~gf&~co)]:
                assert int((flag&mask).sum())==z.loc[z.program_combination.eq(name),'profiles'].iloc[0]
        rows.append({'stage':stage,'source_FB_profiles':len(d),'eligible_sources':int(d.program_source_eligible.sum()),'cell_supported_candidates':int(selected.sum()),
                     'crosschecked_candidates':int(strict.sum()),'original_coordinates_and_IDs_preserved':True})
    # Independent NB implementation includes the common terms omitted by the
    # production score; class-normalized probabilities must still agree.
    stage='P21';take=ann.age.eq(stage).to_numpy();a=ann.loc[take];model=Classifiers(x[take],a.fibroblast_family.to_numpy(),a['sample'].to_numpy(),genes)
    atlas,q=queries(stage,genes);subset=q['refined'][[0,3,17,23,50,101,200,600]]
    p=model.predict(subset);used=subset[:,model.use];depth=np.maximum(used.sum(1),1);ll=[]
    for template in model.templates:
        mu=depth[:,None]*template
        ll.append(nbinom.logpmf(used,model.theta[None,:],model.theta[None,:]/(model.theta[None,:]+mu)).sum(1))
    independent=softmax(np.array(ll).T/np.sqrt(depth[:,None]),axis=1)
    np.testing.assert_allclose(independent,p['negative_binomial']['probabilities'],atol=2e-10)
    # Predictions cannot change merely because other query classes are added.
    extras=q['refined'][1000:1200];joint=model.predict(np.vstack([subset,extras]))
    for m in METHODS:np.testing.assert_allclose(p[m]['probabilities'],joint[m]['probabilities'][:len(subset)],atol=2e-12)
    panel=read(OUT/'published_FB_marker_panel_coverage.tsv')
    def markers(fam):return set(panel.loc[panel.author_family.eq(fam)&panel.measured_shared_panel,'gene'])
    assert 'Col8a1' in markers('Col8a1.Fib')&markers('Lsamp.Fib')
    assert {'Gfpt2','Pi16','Adgrd1'}<=markers('Gfpt2.Fib')
    contexts=read(OUT/'candidate_fibroblast_anatomical_context.tsv');assert contexts.median_um.notna().all()
    (OUT/'numerical_validation.json').write_text(json.dumps({'status':'PASS','stages':rows,'model_proposals_independently_reconstructed':reconstructed,
        'CV_precision_recall_F1_recalculated':True,'heldout_naming_genes_absent_from_training':True,'independent_NB_logpmf_cells':len(subset),
        'query_composition_invariance_checked':True,'program_overlap_counts_reconstructed':True,'new_Xenium_ground_truth_or_fate_validation':False},indent=2)+'\n')
    print('FB_MULTIMETHOD_NUMERICAL_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
