"""Spatial CM/FB/EC synthesis and within-section EC activity–CM maturation tests.

Frozen identities, coordinates and count matrices are reused. No inferred gene
is inserted into the observed Xenium matrix; no animal-level inference is made.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/endothelial-mpl')
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.preprocessing import SplineTransformer, StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/cardiac_spatial_maturation_EC'
STAGES = ['P0', 'P7', 'P14', 'P21']
MAT = pd.read_csv(ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv', sep='\t')
CYCLING = ['Mki67', 'Top2a', 'Cenpf']
GUIDANCE = ['Angpt2', 'Unc5b', 'Plxnd1', 'Itgb3', 'Lama4']
FAMILIES = ['Capillary-like', 'Endocardial-like', 'Arterial-like', 'Venous-like', 'Lymphatic-like', 'Unresolved']
EXPOSURES = ['measured_nuclear_cycling', 'measured_nuclear_cycling_positive',
             'inferred_proliferation', 'inferred_sprouting', 'inferred_missing5',
             'measured_nuclear_guidance', 'near_cycling_EC_30um']
SOURCE_FILES = []

def source(path):
    path = ROOT/path
    SOURCE_FILES.append(path)
    return path

def read(path, **kwargs):
    return pd.read_csv(source(path), sep='\t', index_col=0, **kwargs)

def save(df, name, index=False):
    df.to_csv(OUT/name, sep='\t', index=index, na_rep='NA', float_format='%.9g')

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=str)+'\n')

def simplify(s):
    return s.str.replace(' candidate', '', regex=False).replace({'Endothelial subtype unresolved':'Unresolved'})

def logscore(x, genes, panel):
    ix = np.isin(genes, panel)
    assert ix.sum() == len(panel)
    total = np.asarray(x.sum(1)).ravel()
    return np.log1p(x[:, ix].toarray()*1e4/np.maximum(total, 1)[:, None]).mean(1)

def outcome_withheld_projection(stage, a, ec):
    """Remove all21 CM outcomes and every sprouting target from EC predictors.

    This sensitivity avoids recycling CM maturation genes, including their
    contribution to the query normalization denominator, into inferred activity.
    """
    from project_RNA_sprouting_to_Xenium import fitting_genes, norm_features, target_expression, transfer, OUTCOMES, metrics, family_residual
    rna = ad.read_h5ad(source(f'results/endothelial_QC_subtypes/snRNA/{stage}/EC_corrected_and_raw_counts.h5ad'))
    fitgenes = fitting_genes()
    fitgenes = fitgenes[~np.isin(fitgenes, MAT.gene.tolist()+['Angpt2'])]
    rp = rna.var_names.get_indexer(fitgenes); qp = a.var_names.get_indexer(fitgenes)
    assert (rp>=0).all() and (qp>=0).all()
    rf = norm_features(rna.X[:,rp].toarray())
    qf = norm_features(a.layers['nuclear_only'][:,qp].toarray())
    y = target_expression(rna.X, rna.var_names.get_indexer(OUTCOMES))
    predicted, _, ix, w, sim = transfer(rf, qf, y)
    scores = metrics(predicted)
    for key, value in scores.items():
        name = {'sprouting_6g':'inferred_sprouting_no_CM_panel', 'unmeasured_5g':'inferred_missing5_no_CM_panel',
                'cycling_heldout_control':'inferred_proliferation_no_CM_panel'}[key]
        ec[name] = value
    np.savez_compressed(OUT/stage/'CM_panel_withheld_EC_reference_weights.npz', reference_ids=rna.obs_names.to_numpy(dtype=str),
        query_ids=a.obs_names.to_numpy(dtype=str), fitting_genes=fitgenes, neighbor_indices=ix, neighbor_weights=w,
        reference_gene_targets=np.array(OUTCOMES), reference_target_expression=y, nearest_similarity=sim)
    rows = []
    for lib in sorted(rna.obs['sample'].unique()):
        query = rna.obs['sample'].eq(lib).to_numpy(); train = ~query
        pred, *_ = transfer(rf[train], rf[query], y[train])
        family = rna.obs.loc[query,'EC_family'].to_numpy(dtype=str)
        for key, actual in metrics(y[query]).items():
            predicted = metrics(pred)[key]
            rows.append({'stage':stage,'query_library':lib,'module':key,'fit_genes':len(fitgenes),
                'Spearman':spearmanr(actual,predicted).statistic,
                'within_family_Spearman':spearmanr(family_residual(actual,family),family_residual(predicted,family)).statistic,
                'R2':r2_score(actual,predicted), 'CM_panel_and_all_sprouting_targets_withheld':True})
    return rows

def preparation():
    summaries, coverage, maturity, contexts, subtype_counts, benchmarks = [], [], [], [], [], []
    control = read('results/fibroblast_Xenium_reviewed_refit/CM_maturation_fixed_controls.tsv').reset_index()
    control_genes = control.loc[control.matched_control, 'gene'].tolist()
    for stage in STAGES:
        print('PREPARE', stage, flush=True)
        folder = OUT/stage; folder.mkdir(exist_ok=True)
        atlas = read(f'results/whole_heart_spatial_atlas/{stage}/all_cell_annotations.tsv.gz')
        cm = read(f'results/fibroblast_Xenium_reviewed_refit/{stage}/CM_measured_outcomes.tsv.gz')
        assert atlas.loc[cm.index, 'ventricular_CM_strict'].all()
        np.testing.assert_allclose(cm[['x_centroid', 'y_centroid']], atlas.loc[cm.index, ['x_centroid', 'y_centroid']], atol=1e-7)
        a = ad.read_h5ad(source(f'results/EC_gene_modules/{stage}/Xenium_ECs_with_gene_modules.h5ad'))
        ec = a.obs.copy(); genes = a.var_names.to_numpy(dtype=str)
        assert atlas.loc[ec.index, 'atlas_identity'].eq('Endothelial').all()
        ec['EC_type'] = simplify(ec.nuclear_supported_family)
        ec['inferred_sprouting'] = ec.inferred_sprouting_module
        ec['inferred_missing5'] = ec.inferred_unmeasured_sprouting_module
        ec['inferred_proliferation'] = ec.inferred_proliferation_module
        ec['inferred_sprouting_balanced'] = ec.inferred_sprouting_module_balanced
        benchmarks.extend(outcome_withheld_projection(stage, a, ec))
        for name, x in [('nuclear', a.layers['nuclear_only']), ('native', a.X), ('core', a.layers['boundary_core'])]:
            ec[f'measured_{name}_cycling'] = logscore(x, genes, CYCLING)
            ec[f'measured_{name}_cycling_genes'] = np.asarray((x[:, np.isin(genes, CYCLING)]>0).sum(1)).ravel()
            ec[f'measured_{name}_cycling_positive'] = ec[f'measured_{name}_cycling_genes'].ge(2).astype(float)
            ec[f'{name}_depth'] = np.asarray(x.sum(1)).ravel()
        ec['measured_nuclear_guidance'] = logscore(a.layers['nuclear_only'], genes, GUIDANCE)
        ec['log_nuclear_depth'] = np.log1p(ec.nuclear_depth)
        ec['nuclear_CM_fraction'] = ec.nuclear_only__Cardiomyocyte_fraction
        for t in FAMILIES: ec['type__'+t] = ec.EC_type.eq(t).astype(float)
        selected = ['stage', 'stage_cell_id', 'x_centroid', 'y_centroid', 'EC_type', 'confidence_tier',
                    'nuclear_CM_fraction', 'log_nuclear_depth', 'inferred_sprouting', 'inferred_missing5',
                    'inferred_proliferation', 'inferred_sprouting_balanced', 'measured_nuclear_guidance']
        selected += [c for c in ec if c.startswith('measured_') and ('cycling' in c)]
        selected += [c for c in ec if c.endswith('_no_CM_panel')]
        save(ec[list(dict.fromkeys(selected))], f'{stage}/EC_activity_annotations.tsv.gz', index=True)
        # Alternative CM ownership views, using the same frozen cells and controls.
        b = ad.read_h5ad(source(f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad'))
        bi = b.obs_names.get_indexer(cm.index); assert (bi>=0).all()
        bg = b.var_names.to_numpy(dtype=str)
        with np.load(source(f'results/fibroblast_Xenium_reviewed_refit/{stage}/CM_21gene_expression.npz')) as q:
            np.testing.assert_array_equal(q['cell_ids'], cm.index)
            logs = q['log1pCP10k']; mg = q['genes'].astype(str)
            np.testing.assert_allclose(logs.mean(1), cm.maturation_21gene_mean, rtol=2e-6)
            cm['maturation_without_Myl2'] = logs[:, mg!='Myl2'].mean(1)
            for cat, panel in MAT.groupby('figure_category'):
                cm['maturation_category__'+cat] = logs[:, np.isin(mg, panel.gene)].mean(1)
        for view, matrix in [('refined', b.X), ('core', b.layers['high_confidence']), ('nuclear', b.layers['nuclear_only'])]:
            x = matrix[bi]
            cm['maturation_'+view] = logscore(x, bg, MAT.gene.tolist())
            cm['maturation_'+view+'_matched_control'] = cm['maturation_'+view]-logscore(x, bg, control_genes)
            cm['CM_'+view+'_depth'] = np.asarray(x.sum(1)).ravel()
        del b
        cm['stage_cell_id'] = stage+'::'+cm.index.astype(str)
        fine = read(f'results/ventricular_spatial_resolution/{stage}/fine_CM_spatial_annotations.tsv.gz')
        cm['fine_group_proposal'] = fine.fine_group_proposal.reindex(cm.index)
        cm['fine_candidate_supported'] = fine.prioritized_fine_candidate.reindex(cm.index).fillna(False).astype(bool)
        cm['fine_candidate_published_CM_compatible'] = fine.prioritized_and_published_CM_compatible.reindex(cm.index).fillna(False).astype(bool)
        cm['author_ventricular_compatible'] = atlas.loc[cm.index, 'author_ventricular']
        save(cm, f'{stage}/ventricular_CM_maturation.tsv.gz', index=True)
        xy = cm[['x_centroid', 'y_centroid']].to_numpy()
        exy = ec[['x_centroid', 'y_centroid']].to_numpy()
        tree = cKDTree(exy)
        cm['nearest_eligible_EC_um'] = tree.query(xy)[0]
        cycling = ec.measured_nuclear_cycling_positive.eq(1).to_numpy()
        cm['nearest_cycling_EC_um'] = cKDTree(exy[cycling]).query(xy)[0]
        cm['nearest_no_cycling_marker_EC_um'] = cKDTree(exy[ec.measured_nuclear_cycling_genes.eq(0)]).query(xy)[0]
        cm['near_cycling_EC_30um'] = cm.nearest_cycling_EC_um.le(30).astype(float)
        edges = cKDTree(xy).sparse_distance_matrix(tree, 100, output_type='coo_matrix')
        values = EXPOSURES[:-1]+['inferred_sprouting_balanced', 'measured_native_cycling', 'measured_core_cycling',
                  'measured_native_cycling_positive', 'log_nuclear_depth', 'nuclear_CM_fraction']+['type__'+t for t in FAMILIES]
        values += [c for c in ec if c.endswith('_no_CM_panel')]
        for radius in [30, 50, 100]:
            m = edges.data<=radius
            adjacency = sparse.csr_matrix((np.ones(m.sum()), (edges.row[m], edges.col[m])), shape=(len(cm), len(ec)))
            n = np.asarray(adjacency.sum(1)).ravel()
            env = read(f'results/fibroblast_Xenium_reviewed_refit/{stage}/CM_neighborhoods_{radius}um.tsv.gz')
            d = cm.join(env[['resolved_neighbors', 'fibroblast_fraction', 'endothelial_fraction', 'mural_fraction', 'immune_fraction']])
            d['eligible_EC_neighbors'] = n
            for key in values:
                d['local_'+key] = np.divide(np.asarray(adjacency@ec[key].to_numpy(float)).ravel(), n, out=np.full(len(n), np.nan), where=n>0)
            d['model_eligible'] = (n>=3)&d.resolved_neighbors.ge(5)&d.arterial_distance_um.notna()
            save(d, f'{stage}/CM_EC_neighborhoods_{radius}um.tsv.gz', index=True)
            coverage.append({'stage':stage, 'radius_um':radius, 'ventricular_CMs':len(d), 'eligible_CMs':int(d.model_eligible.sum()),
                             'eligible_near_cycling_EC_30um':int((d.model_eligible&d.near_cycling_EC_30um.eq(1)).sum()),
                             'missing_nuclear_area_retained':int((d.model_eligible&d.nuclear_area_missing).sum())})
        for key in ['maturation_21gene_mean', 'maturation_matched_control', 'maturation_without_Myl2']:
            maturity.append({'stage':stage, 'score':key, 'CMs':len(cm), 'mean':cm[key].mean(), 'median':cm[key].median()})
        for key in ['measured_nuclear_cycling_positive', 'measured_native_cycling_positive', 'measured_core_cycling_positive']:
            summaries.append({'stage':stage, 'score':key, 'ECs':len(ec), 'positive_ECs':int(ec[key].sum()), 'fraction':ec[key].mean()})
        # Descriptive local lineage composition for supported subtype populations.
        fb = read(f'results/fibroblast_spatial_consensus/{stage}/fibroblast_consensus_annotations.tsv.gz')
        ectypes = read('results/endothelial_five_types/Xenium_five_type_annotations.tsv.gz')
        ectypes = ectypes.loc[ectypes.stage.eq(stage)].copy()
        rows = [cm[['x_centroid','y_centroid']].assign(lineage='Ventricular CM', subtype='All strict ventricular CMs')]
        for t in sorted(fb.high_confidence_family.unique()):
            if t.startswith('Unresolved'): continue
            rows.append(fb.loc[fb.high_confidence_family.eq(t), ['x_centroid','y_centroid']].assign(lineage='Fibroblast', subtype=t))
        for scope, col in [('EC strict', 'EC_type'), ('EC nuclear-supported', 'nuclear_supported_EC_type')]:
            z = ectypes.loc[ectypes.native_EC_broad&~ectypes.multiple_projected_nuclear_components_review&~ectypes.geometry_unresolved_review]
            for t in FAMILIES[:-1]:
                rows.append(z.loc[z[col].eq(t), ['x_centroid','y_centroid']].assign(lineage=scope, subtype=t))
        for group, z in cm.loc[cm.fine_candidate_supported].groupby('fine_group_proposal'):
            rows.append(z[['x_centroid','y_centroid']].assign(lineage='Ventricular CM fine candidate', subtype=group))
        resolved = atlas.loc[atlas.atlas_resolved]
        lineages = ['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
        for z in rows:
            if not len(z): continue
            label = z.lineage.iloc[0]; subtype = z.subtype.iloc[0]
            zxy = z[['x_centroid','y_centroid']].to_numpy()
            den = cKDTree(resolved[['x_centroid','y_centroid']]).query_ball_point(zxy, 30, return_length=True)-1
            valid = den>=5
            result = {'stage':stage, 'lineage':label, 'subtype':subtype, 'profiles':len(z), 'eligible_30um':int(valid.sum()), 'sparse_lt30':int(valid.sum())<30}
            for l in lineages:
                num = cKDTree(resolved.loc[resolved.atlas_identity.eq(l), ['x_centroid','y_centroid']]).query_ball_point(zxy, 30, return_length=True)
                if (l=='Cardiomyocyte' and label.startswith('Ventricular')) or (l=='Fibroblast' and label=='Fibroblast') or (l=='Endothelial' and label.startswith('EC')): num-=1
                result['fraction_'+l] = np.mean(num[valid]/den[valid]) if valid.any() else np.nan
            contexts.append(result)
            subtype_counts.append({'stage':stage, 'lineage':label, 'subtype':subtype, 'profiles':len(z)})
        del a
    save(pd.DataFrame(coverage), 'neighborhood_coverage.tsv')
    save(pd.DataFrame(maturity), 'maturation_stage_summary.tsv')
    save(pd.DataFrame(summaries), 'EC_measured_cycling_summary.tsv')
    save(pd.DataFrame(contexts), 'subtype_lineage_neighborhoods_30um.tsv')
    save(pd.DataFrame(subtype_counts), 'spatial_subtype_counts.tsv')
    save(pd.DataFrame(benchmarks), 'CM_panel_withheld_RNA_library_benchmarks.tsv')
    dump(OUT/'source_files.json', sorted(set(str(p.relative_to(ROOT)) for p in SOURCE_FILES)))

def design(d, spatial=True):
    b = pd.DataFrame(index=d.index)
    for key in ['outcome_excluded_depth','cell_area','resolved_neighbors','eligible_EC_neighbors','nearest_eligible_EC_um']:
        b['log_'+key] = np.log1p(d[key].to_numpy(float))
    for key in ['nuclear_area_missing','fibroblast_fraction','endothelial_fraction','mural_fraction','immune_fraction',
                'local_log_nuclear_depth','local_nuclear_CM_fraction','Slit2_log','cycling_log']:
        b[key] = d[key].to_numpy(float)
    # Drop the unresolved category as the composition reference; constant columns are removed in fitting.
    for t in FAMILIES[:-1]: b['local_type__'+t] = d['local_type__'+t].to_numpy(float)
    if spatial:
        for key in ['arterial_distance_um','endocardial_proxy_distance_um']:
            b['log_'+key] = np.log1p(d[key].to_numpy(float))
        xy = d[['x_centroid','y_centroid']].to_numpy()
        sp = SplineTransformer(n_knots=5, degree=3, include_bias=False).fit_transform(xy)
        for j in range(sp.shape[1]): b['spatial_'+str(j)] = sp[:,j]
        z = (xy-xy.mean(0))/np.maximum(xy.std(0),1)
        b['spatial_xy'] = z[:,0]*z[:,1]
    return b

def fit(y, b, e, xy, block=400, draws=256):
    y = np.asarray(y, float); e = np.asarray(e, float)
    ys = max(y.std(), 1e-9); es = max(e.std(), 1e-9)
    yy = (y-y.mean())/ys
    v = b.to_numpy(float); sd = v.std(0); v = v[:,sd>1e-8]; sd = sd[sd>1e-8]
    v = (v-v.mean(0))/sd
    xx = np.column_stack([np.ones(len(y)), v, (e-e.mean())/es])
    beta = np.linalg.lstsq(xx, yy, rcond=1e-9)[0]
    residual = yy-xx@beta
    blocks = np.floor(xy/block).astype(int)
    _, inv = np.unique(blocks, axis=0, return_inverse=True)
    nb = inv.max()+1
    h = np.zeros((nb, xx.shape[1], xx.shape[1])); g = np.zeros((nb, xx.shape[1]))
    for k in range(nb):
        m = inv==k; h[k] = xx[m].T@xx[m]; g[k] = xx[m].T@yy[m]
    rng = np.random.default_rng(651)
    samples = []
    for _ in range(draws):
        w = np.bincount(rng.integers(0,nb,nb), minlength=nb)
        samples.append((np.linalg.pinv(np.einsum('b,bij->ij',w,h), rcond=1e-9)@np.einsum('b,bi->i',w,g))[-1])
    low, high = np.quantile(samples, [.025,.975])
    without = xx[:,:-1]; yr = yy-without@np.linalg.lstsq(without,yy,rcond=1e-9)[0]
    er = e-without@np.linalg.lstsq(without,e,rcond=1e-9)[0]
    return {'CMs':len(y), 'spatial_blocks':nb, 'block_um':block, 'bootstrap_draws':draws,
            'beta_ySD_per_exposureSD':beta[-1], 'block_CI_low':low, 'block_CI_high':high,
            'raw_unit_coefficient':beta[-1]*ys/es, 'raw_unit_CI_low':low*ys/es, 'raw_unit_CI_high':high*ys/es,
            'raw_Spearman':spearmanr(y,e).statistic, 'adjusted_partial_Pearson':np.corrcoef(yr,er)[0,1],
            'design_rank':np.linalg.matrix_rank(xx), 'design_columns':xx.shape[1], 'residual_SD':residual.std()}

def crossvalidate(y, b, e, xy):
    blocks = np.floor(xy/400).astype(int)
    folds = ((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
    predictions = [np.full(len(y),np.nan), np.full(len(y),np.nan)]
    v = b.to_numpy(float)
    for fold in range(5):
        test=folds==fold; train=~test
        if not test.any(): continue
        pos=np.flatnonzero(train)
        train[pos[cKDTree(xy[test]).query(xy[train])[0]<100]]=False
        if train.sum()<200: continue
        for j, x in enumerate([v,np.column_stack([v,e])]):
            scaler=StandardScaler().fit(x[train]); model=Ridge(alpha=10).fit(scaler.transform(x[train]),y[train])
            predictions[j][test]=model.predict(scaler.transform(x[test]))
    valid=np.isfinite(predictions[0])&np.isfinite(predictions[1])
    if not valid.any(): return {'OOF_CMs':0}
    values=[r2_score(y[valid],p[valid]) for p in predictions]
    return {'OOF_CMs':int(valid.sum()),'baseline_OOF_R2':values[0],'plus_activity_OOF_R2':values[1],'added_OOF_R2':values[1]-values[0]}

def models():
    coefficients, predictive, contrasts = [], [], []
    for stage in STAGES:
        for radius in [30,50,100]:
            d=read(f'results/cardiac_spatial_maturation_EC/{stage}/CM_EC_neighborhoods_{radius}um.tsv.gz')
            d=d.loc[d.model_eligible].copy()
            print('MODEL',stage,radius,len(d),flush=True)
            jobs=[]
            for e in EXPOSURES:
                jobs.append((e, 'maturation_21gene_mean','full', np.ones(len(d),bool),400,256 if radius==50 else 96))
            if radius==50:
                main=['measured_nuclear_cycling','inferred_sprouting','inferred_proliferation','near_cycling_EC_30um']
                for e in main:
                    jobs.append((e,'maturation_21gene_mean','capture_density_family',np.ones(len(d),bool),400,128))
                    for outcome in ['maturation_matched_control','maturation_without_Myl2','maturation_core','maturation_core_matched_control']:
                        jobs.append((e,outcome,'full',d.CM_core_depth.ge(20).to_numpy() if 'core' in outcome else np.ones(len(d),bool),400,128))
                    jobs.append((e,'maturation_21gene_mean','full_low_EC_CM_RNA',d.local_nuclear_CM_fraction.le(.02).to_numpy(),400,128))
                    jobs.append((e,'maturation_21gene_mean','full_author_ventricular',d.author_ventricular_compatible.to_numpy(bool),400,128))
                    for block in [200,800]: jobs.append((e,'maturation_21gene_mean','full',np.ones(len(d),bool),block,128))
                for e in ['measured_native_cycling','measured_core_cycling','inferred_missing5','inferred_sprouting_balanced',
                          'inferred_sprouting_no_CM_panel','inferred_missing5_no_CM_panel','inferred_proliferation_no_CM_panel']:
                    if e not in EXPOSURES: jobs.append((e,'maturation_21gene_mean','full',np.ones(len(d),bool),400,128))
            for e,outcome,model,selected,block,draws in jobs:
                s=d.loc[selected].copy()
                key=e if e.startswith('near_') else 'local_'+e
                if len(s)<200 or s[key].std()<1e-8: continue
                b=design(s,spatial=model!='capture_density_family'); xy=s[['x_centroid','y_centroid']].to_numpy()
                result=fit(s[outcome],b,s[key],xy,block,draws)
                coefficients.append({'stage':stage,'radius_um':radius,'exposure':e,'outcome':outcome,'model':model,**result})
                if radius==50 and model=='full' and outcome=='maturation_21gene_mean' and block==400 and e in EXPOSURES:
                    predictive.append({'stage':stage,'exposure':e,**crossvalidate(s[outcome].to_numpy(),b,s[key].to_numpy(),xy)})
                    if e=='near_cycling_EC_30um':
                        split=s[key].eq(1).to_numpy()
                        contrasts.append({'stage':stage,'near_CMs':int(split.sum()),'other_CMs':int((~split).sum()),
                                          'near_mean_maturation':s.loc[split,outcome].mean(),'other_mean_maturation':s.loc[~split,outcome].mean(),
                                          'unadjusted_near_minus_other':s.loc[split,outcome].mean()-s.loc[~split,outcome].mean(),
                                          'adjusted_near_minus_other':result['raw_unit_coefficient'],
                                          'block_CI_low':result['raw_unit_CI_low'],'block_CI_high':result['raw_unit_CI_high']})
        pd.DataFrame(coefficients).to_csv(OUT/'association_coefficients.tsv',sep='\t',index=False)
        save(pd.DataFrame(predictive),'spatial_blocked_prediction.tsv')
        save(pd.DataFrame(contrasts),'near_cycling_EC_maturation_contrasts.tsv')

def extra_checks():
    a=pd.read_csv(OUT/'association_coefficients.tsv',sep='\t')
    rows=[];subtypes=[];restricted_coverage=[]
    evidence=pd.read_csv(ROOT/'results/ventricular_spatial_resolution/fine_group_evidence_summary.tsv',sep='\t').set_index('group')
    for stage in STAGES:
        d=read(f'results/cardiac_spatial_maturation_EC/{stage}/CM_EC_neighborhoods_50um.tsv.gz')
        for group,z in d.loc[d.fine_candidate_supported].groupby('fine_group_proposal'):
            subtypes.append({'stage':stage,'group':group,'candidate_identity':evidence.loc[group,'candidate_identity'],
                'functional_family':evidence.loc[group,'functional_family'],'CMs':len(z),
                'median_maturation_21gene':z.maturation_21gene_mean.median(),
                'median_endocardial_proxy_distance_um':z.endocardial_proxy_distance_um.median(),
                'median_candidate_arterial_distance_um':z.arterial_distance_um.median(),
                'fraction_within_30um_cycling_EC':z.near_cycling_EC_30um.mean(),
                'fine_labels_not_independent_of_maturation_genes':True})
        d=d.loc[d.model_eligible].copy();xy=d[['x_centroid','y_centroid']].to_numpy();b=design(d)
        for e in ['measured_nuclear_cycling','inferred_sprouting','inferred_proliferation','near_cycling_EC_30um']:
            key=e if e.startswith('near_') else 'local_'+e
            result=fit(d.maturation_21gene_mean,b.drop(columns=['Slit2_log','cycling_log']),d[key],xy,400,128)
            rows.append({'stage':stage,'radius_um':50,'exposure':e,'outcome':'maturation_21gene_mean','model':'full_without_CM_state',**result})
            for category in sorted(MAT.figure_category.unique()):
                outcome='maturation_category__'+category
                result=fit(d[outcome],b,d[key],xy,400,128)
                rows.append({'stage':stage,'radius_um':50,'exposure':e,'outcome':outcome,'model':'full',**result})
        proximity=np.log1p(d.nearest_no_cycling_marker_EC_um)-np.log1p(d.nearest_cycling_EC_um)
        result=fit(d.maturation_21gene_mean,b,proximity,xy,400,256)
        rows.append({'stage':stage,'radius_um':50,'exposure':'relative_proximity_to_cycling_EC','outcome':'maturation_21gene_mean','model':'full',**result})
        # Rebuild neighborhoods after filtering individual EC profiles, rather
        # than declaring a low average CM-marker fraction to be pure EC RNA.
        ec=read(f'results/cardiac_spatial_maturation_EC/{stage}/EC_activity_annotations.tsv.gz')
        clean=ec.loc[ec.nuclear_CM_fraction.le(.05)].copy()
        exy=clean[['x_centroid','y_centroid']].to_numpy()
        edges=cKDTree(xy).sparse_distance_matrix(cKDTree(exy),50,output_type='coo_matrix')
        mat=sparse.csr_matrix((np.ones(len(edges.data)),(edges.row,edges.col)),shape=(len(d),len(clean)))
        n=np.asarray(mat.sum(1)).ravel();s=d.copy();s['eligible_EC_neighbors']=n
        for key in ['measured_nuclear_cycling','inferred_sprouting','inferred_proliferation','log_nuclear_depth','nuclear_CM_fraction']:
            s['local_'+key]=np.divide(np.asarray(mat@clean[key].to_numpy(float)).ravel(),n,out=np.full(len(n),np.nan),where=n>0)
        for t in FAMILIES:
            s['local_type__'+t]=np.divide(np.asarray(mat@clean.EC_type.eq(t).to_numpy(float)).ravel(),n,out=np.full(len(n),np.nan),where=n>0)
        s['nearest_eligible_EC_um']=cKDTree(exy).query(xy)[0]
        positive=clean.measured_nuclear_cycling_positive.eq(1).to_numpy()
        s['nearest_cycling_EC_um']=cKDTree(exy[positive]).query(xy)[0]
        s['near_cycling_EC_30um']=s.nearest_cycling_EC_um.le(30).astype(float)
        s=s.loc[s.eligible_EC_neighbors.ge(3)].copy()
        save(s,f'{stage}/CM_low_CM_signal_EC_neighborhoods_50um.tsv.gz',index=True)
        restricted_coverage.append({'stage':stage,'ECs_before':len(ec),'ECs_CM_marker_fraction_le05':len(clean),
            'cycling_positive_ECs_remaining':int(positive.sum()),'CMs_with_3_restricted_ECs':len(s),
            'CMs_local_mean_CM_fraction_le02':int(d.local_nuclear_CM_fraction.le(.02).sum()),
            'strict_2percent_mean_model_skipped_lt200':bool((d.local_nuclear_CM_fraction<=.02).sum()<200)})
        if len(s)>=200:
            for e in ['measured_nuclear_cycling','inferred_sprouting','inferred_proliferation','near_cycling_EC_30um']:
                key=e if e.startswith('near_') else 'local_'+e
                if s[key].std()<1e-8:continue
                result=fit(s.maturation_21gene_mean,design(s),s[key],s[['x_centroid','y_centroid']].to_numpy(),400,128)
                rows.append({'stage':stage,'radius_um':50,'exposure':e,'outcome':'maturation_21gene_mean','model':'full_individual_EC_CM_markers_le05',**result})
    # Replace only this named check set, so repeated runs cannot duplicate rows.
    oldextra=a.model.isin(['full_without_CM_state','full_individual_EC_CM_markers_le05'])|a.outcome.str.startswith('maturation_category__')|a.exposure.eq('relative_proximity_to_cycling_EC')
    save(pd.concat([a.loc[~oldextra],pd.DataFrame(rows)],ignore_index=True),'association_coefficients.tsv')
    save(pd.DataFrame(subtypes),'fine_CM_candidate_spatial_summary.tsv')
    save(pd.DataFrame(restricted_coverage),'EC_CM_marker_restriction_coverage.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--part',choices=['prepare','models','extra','all'],default='all');args=p.parse_args()
    OUT.mkdir(exist_ok=True)
    if args.part in ['prepare','all']: preparation()
    if args.part in ['models','all']: models()
    if args.part in ['extra','all']: extra_checks()

if __name__=='__main__': main()
