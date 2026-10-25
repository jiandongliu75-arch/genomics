"""Independent arithmetic, graph, holdout and method-adoption checks."""
from revisit_vascular_coordination import ROOT,OUT,UP,STAGES,CFG,MAT,read,save,dump
from fit_wormhole_cardiac_niches import cloud_subset
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.optimize import linprog
from scipy.stats import spearmanr
from sklearn.preprocessing import SplineTransformer
from shapely import Polygon
from pathlib import Path
from collections import deque
import hashlib,json,subprocess,gc,importlib.metadata,re
from PIL import Image

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def bfs(graph,start,seeds,limit=10):
    if start in seeds:return 0
    seen={start};q=deque([(start,0)])
    while q:
        u,n=q.popleft()
        if n>=limit:continue
        for v in graph.indices[graph.indptr[u]:graph.indptr[u+1]]:
            if v in seeds:return n+1
            if v not in seen:seen.add(v);q.append((v,n+1))
    return np.nan

def validate():
    rng=np.random.default_rng(5129);checks=[];sources=set(json.loads((OUT/'source_files.json').read_text()))
    controls=read(OUT/'feasible_allcell_matched_controls.tsv',index_col=None)
    assert not set(controls.control_gene)&set(MAT.gene)
    assert controls.groupby('target_gene').size().eq(CFG['matched_controls_per_gene']).all()
    co=read(OUT/'paired_geometry_coefficients.tsv.gz',index_col=None)
    cover=read(OUT/'paired_geometry_cohort.tsv',index_col=None).set_index('stage')
    for stage in STAGES:
        print('VALIDATE_DYNAMICS',stage,flush=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz');q=read(OUT/stage/'measured_vascular_states.tsv.gz')
        c=read(OUT/stage/'CM_outcomes_and_context.tsv.gz');e=read(OUT/stage/'CM_vascular_exposures.tsv.gz');r=read(OUT/stage/'CM_vascular_rings.tsv.gz')
        np.testing.assert_array_equal(q.index,f.index);assert q.index.is_unique
        assert q.EC_primary.equals(f.atlas_identity.eq('Endothelial')&f.program_source_eligible)
        assert q.author_myofiber_EC.equals(q.author_atlas_cell_state.fillna('').str.rstrip('.').eq('cap.EC.High.Myo'))
        assert q.author_myofiber_EC.sum()>0
        assert c.index.equals(f.index[f.ventricular_CM_strict])
        ids=rng.choice(q.index[q.EC_primary],32,replace=False)
        a=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        assert a.obs_names.equals(q.index);genes=a.var_names.to_numpy(dtype=str);assert len(genes)==403 and 'Emcn' not in genes
        xi=a.obs_names.get_indexer(ids)
        for view,layer in [('core','high_confidence'),('nuclear','nuclear_only')]:
            x=a.layers[layer][xi].toarray();depth=x.sum(1);z=np.log1p(x*1e4/np.maximum(depth,1)[:,None])
            np.testing.assert_allclose(depth,q.loc[ids,view+'__panel_depth'])
            for gene in ['Ptprm','Flt4','Dll4','Hey1']:
                j=list(genes).index(gene);np.testing.assert_allclose(z[:,j],q.loc[ids,view+'__'+gene],rtol=3e-6,atol=1e-7)
                np.testing.assert_allclose(x[:,j],q.loc[ids,view+'_molecules__'+gene])
            for name,panel in CFG['vascular_programs'].items():
                np.testing.assert_allclose(z[:,np.isin(genes,panel)].mean(1),q.loc[ids,view+'__'+name],rtol=3e-6,atol=1e-7)
        cmids=rng.choice(c.index,32,replace=False);ci=a.obs_names.get_indexer(cmids);x=a.layers['high_confidence'][ci].toarray()
        z=np.log1p(x*1e4/np.maximum(x.sum(1),1)[:,None]);m=z[:,np.isin(genes,MAT.gene)].mean(1)
        np.testing.assert_allclose(m,c.loc[cmids,'maturity_core_21'],rtol=3e-6)
        cc=m-z[:,np.isin(genes,controls.control_gene)].mean(1)
        np.testing.assert_allclose(cc,c.loc[cmids,'maturity_allcell_matched_control'],rtol=3e-6,atol=1e-7)
        a.file.close();del a,x,z;gc.collect()
        selected=q.EC_primary&~q.EC_nuclear_type.isin(['Endocardial-like','Lymphatic-like'])&q.core_molecules__Flt4.gt(0)&q.core_molecules__Dll4.gt(0)
        assert selected.equals(q.EC_growth_compatible)
        edges=read(OUT/stage/'native_boundary_edges.tsv.gz',index_col=None);picked=edges.iloc[rng.choice(len(edges),16,replace=False)]
        ids_needed=q.index[np.unique(picked[['profile_i','profile_j']].to_numpy(int))]
        raw=next(Path('/workspace/heart-atlas/data/xenium').glob('WT_'+('P07' if stage=='P7' else stage)+'__*'))/'cell_boundaries.csv.gz'
        vv=pd.read_csv(raw,usecols=['cell_id','vertex_x','vertex_y']);vv=vv[vv.cell_id.isin(ids_needed)]
        poly={i:Polygon(v[['vertex_x','vertex_y']].to_numpy()) for i,v in vv.groupby('cell_id',sort=False)}
        graphs={threshold:sparse.load_npz(OUT/stage/f'native_boundary_{threshold:g}um_graph.npz').tocsr() for threshold in CFG['native_boundary_distances_um']}
        for row in picked.itertuples():
            i,j=int(row.profile_i),int(row.profile_j);dist=poly[q.index[i]].distance(poly[q.index[j]])
            np.testing.assert_allclose(dist,row.native_polygon_distance_um,rtol=2e-6,atol=1e-7)
            for threshold,g in graphs.items():assert bool(g[i,j])==(dist<=threshold)
        for g in graphs.values():assert g.diagonal().sum()==0 and (g-g.T).nnz==0 and g.shape==(len(q),len(q))
        for geometry,path in [('proximity50',UP/stage/'proximity50_graph.npz'),('refined_contact',UP/stage/'contact_graph.npz'),('native_boundary',OUT/stage/'native_boundary_1um_graph.npz')]:
            g=sparse.load_npz(path).tocsr()
            for cid in cmids[:12]:
                i=q.index.get_loc(cid);near=q.iloc[g.indices[g.indptr[i]:g.indptr[i+1]]];ec=near[near.EC_primary]
                expected=ec.nuclear__Ptprm.mean();np.testing.assert_allclose(expected,e.loc[cid,geometry+'__EC_sender_Ptprm'],rtol=3e-6,atol=1e-7,equal_nan=True)
                np.testing.assert_allclose(np.log1p(ec.nuclear__panel_depth).mean(),e.loc[cid,geometry+'__EC_source_log_nuclear_capture'],rtol=3e-6,atol=1e-7,equal_nan=True)
        for geometry,path in [('Delaunay',OUT/stage/'pruned_Delaunay_graph.npz'),('refined_contact',UP/stage/'contact_graph.npz')]:
            g=sparse.load_npz(path).tocsr();seeds=set(np.flatnonzero(q.EC_growth_compatible.to_numpy()))
            for cid in cmids[:6]:
                expected=bfs(g,q.index.get_loc(cid),seeds);got=r.loc[cid,geometry+'_rings__Growth_compatible']
                np.testing.assert_allclose(expected,got,equal_nan=True)
            key=geometry+'_rings__Growth_compatible';assert r[key].isna().equals(r[geometry+'_outside10__Growth_compatible'])
        with np.load(OUT/stage/'nonCM_state_points.npz') as b:
            pts=b['points'];source_ids=b['source_ids'].astype(str);focal_ids=b['focal_ids'].astype(str);train=b['training_focal_indices'];test=b['heldout_focal_indices'];pca_sources=b['PCA_training_source_indices']
        member=sparse.load_npz(OUT/stage/'state_cloud_membership.npz').tocsr()
        assert not np.isin(pca_sources,np.unique(member[test].indices)).any()
        xy=f.loc[focal_ids,['x_centroid','y_centroid']].to_numpy();assert cKDTree(xy[test]).query(xy[train])[0].min()>=CFG['holdout_buffer_um']
        assert not np.intersect1d(np.unique(member[train].indices),np.unique(member[test].indices)).size
        assert f.loc[source_ids,'program_source_eligible'].all() and not f.loc[source_ids,'atlas_identity'].eq('Cardiomyocyte').any()
        assert np.isfinite(pts).all()
        with np.load(OUT/stage/'state_niche_embeddings.npz') as b:
            assert np.isfinite(b['embedding']).all() and len(b['cell_ids'])==len(train)+len(test)
        assert co[co.stage.eq(stage)].CMs.eq(cover.loc[stage,'identical_CMs_all_three_geometries']).all()
        checks.append({'stage':stage,'profiles':len(f),'native_polygon_pairs_recomputed':16,'program_and_score_rows_recomputed':32,
            'graph_ring_BFS_checks':12,'source_overlap_train_test':0,'CMs_paired_geometry':int(cover.loc[stage,'identical_CMs_all_three_geometries']),'status':'PASS'})
    # Independently reproduce one fully adjusted coefficient via Frisch-Waugh-
    # Lovell residualization, rather than using the production block solver.
    c=read(OUT/'P7/CM_outcomes_and_context.tsv.gz');e=read(OUT/'P7/CM_vascular_exposures.tsv.gz')
    names=['Plexus_growth_proxy','Notch_response_proxy','Cycling','Capillary_fraction','Endocardial_fraction','EC_sender_Ptprm','CM_sender_Ptprm']
    take=c.CM_outcome_excluded_depth.ge(20)&c.arterial_distance_um.notna()&c.endocardial_proxy_distance_um.notna()
    for geom in ['proximity50','refined_contact','native_boundary']:take &= e[geom+'__eligible_EC_neighbors'].ge(1)&e[geom+'__resolved_neighbors'].ge(3)&e[[geom+'__'+p for p in names]].notna().all(axis=1)
    c=c.loc[take];e=e.loc[take];geom='proximity50'
    controls_num=[np.log1p(c.cell_area),np.log1p(c.CM_outcome_excluded_depth),c.nucleus_area.isna().astype(float)]
    for k in ['EC_source_log_nuclear_capture','EC_source_log_core_capture','EC_source_coverage']:controls_num.append(e[geom+'__'+k])
    controls_num.append(np.log1p(e[geom+'__resolved_neighbors']))
    for l in ['FB','SMC','Immune']:controls_num.append(e[geom+'__'+l+'_fraction_all_resolved'])
    for k in ['arterial_distance_um','endocardial_proxy_distance_um']:controls_num.append(np.log1p(c[k]))
    spl=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(c[['x_centroid','y_centroid']])
    controls_num.extend([spl[:,j] for j in range(spl.shape[1])]);controls_num.extend([c.CM_cycling,c.CM_Slit2,c.CM_Ptprm_core])
    controls_num.extend([e[geom+'__'+n] for n in names if n!='EC_sender_Ptprm'])
    b=np.column_stack(controls_num);b=b[:,b.std(0)>1e-8];b=(b-b.mean(0))/b.std(0);b=np.column_stack([np.ones(len(c)),b])
    x=e[geom+'__EC_sender_Ptprm'].to_numpy();x=(x-x.mean())/x.std();y=c.maturity_core_21.to_numpy();y=(y-y.mean())/y.std()
    xr=x-b@np.linalg.lstsq(b,x,rcond=1e-9)[0];yr=y-b@np.linalg.lstsq(b,y,rcond=1e-9)[0];estimate=xr@yr/(xr@xr)
    got=co[co.stage.eq('P7')&co.geometry.eq(geom)&co.exposure.eq('EC_sender_Ptprm')&co.outcome.eq('maturity_core_21')].iloc[0]
    np.testing.assert_allclose(estimate,got.standardized_coefficient,rtol=3e-6,atol=1e-8)
    benchmark=read(OUT/'Wormhole_reserved_distance_benchmark.tsv',index_col=None)
    metric=read(OUT/'Wormhole_metric_assessment.tsv',index_col=None).set_index('stage')
    for stage,z in benchmark.groupby('stage'):
        assert len(z)==64 and np.isfinite(z.select_dtypes('number')).all().all()
        w=spearmanr(z.Sinkhorn_divergence,z.Wormhole_mean_squared_embedding_distance).statistic
        m=spearmanr(z.Sinkhorn_divergence,z.mean_squared_feature_distance).statistic
        assert bool(metric.loc[stage,'metric_gate'])==bool(w>=.8 and w>=m-.05)
        with np.load(OUT/stage/'nonCM_state_points.npz') as b:pts=b['points'];focal=b['focal_ids'].astype(str)
        member=sparse.load_npz(OUT/stage/'state_cloud_membership.npz').tocsr()
        with np.load(OUT/'Wormhole_sample_embeddings.npz') as b:emb=b['heldout_embeddings']
        meta=read(OUT/'Wormhole_heldout_clouds.tsv',index_col=None)
        for row in z.iloc[:2].itertuples():
            clouds=[];ii=[]
            for cid in [row.cell_id_a,row.cell_id_b]:
                i=int(np.flatnonzero(focal==cid)[0]);pc,_=cloud_subset(pts,member,i,CFG['seed']);clouds.append(pc)
                ii.append(int(np.flatnonzero(meta.stage.eq(stage)&meta.cell_id.eq(cid))[0]))
            np.testing.assert_allclose(np.mean((emb[ii[0]]-emb[ii[1]])**2),row.Wormhole_mean_squared_embedding_distance,rtol=3e-6)
            a,b=clouds;na,nb=len(a),len(b);cost=cdist(a,b,'sqeuclidean')
            # Independent dense marginal construction for the small LP.
            A=np.zeros((na+nb,na*nb))
            for i in range(na):A[i,i*nb:(i+1)*nb]=1
            for j in range(nb):A[na+j,j::nb]=1
            solution=linprog(cost.ravel(),A_eq=A,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs')
            assert solution.success;np.testing.assert_allclose(solution.fun,row.exact_W2_squared,rtol=2e-6)
    gate=json.loads((OUT/'dynamic_interpretation_gate.json').read_text());assert not gate['Wormhole_metric_accepted_all_ages']
    assert not gate['tracked_cells'] and not gate['endocardial_to_coronary_lineage_inference'] and not gate['causal_influence']
    loss=read(OUT/'Wormhole_training_loss.tsv',index_col=None);assert len(loss)==1200 and np.isfinite(loss).all().all()
    assert loss.encoder_loss.tail(100).mean()<loss.encoder_loss.head(100).mean()
    feature=json.loads((OUT/'state_feature_definition.json').read_text());assert not set(feature['genes'])&set(MAT.gene)
    assert not feature['focal_CM_RNA_used'] and not feature['maturation_genes_used']
    feasibility=read(OUT/'reported_maturation_control_feasibility.tsv',index_col=None);assert not feasibility.reported100_without_replacement_feasible.any()
    pngs=sorted(OUT.glob('*.png'));assert len(pngs)==13
    for path in pngs:
        with Image.open(path) as im:im.verify()
    links=0
    for doc in [OUT/'README.md',OUT/'METHODS_REVIEW.md']:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text()):
            if target.startswith(('https://','http://','#')):continue
            assert (doc.parent/target.split('#')[0]).exists(),(doc,target)
            links+=1
    info=subprocess.check_output(['pdfinfo',str(OUT/'cardiac_niche_coordination_report.pdf')],text=True)
    assert 'Pages:           13' in info
    save(pd.DataFrame(checks),'validation_checks.tsv')
    sources|={'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',
        '/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad',
        '/workspace/heart-atlas/metadata/paper.xml','/workspace/heart-atlas/metadata/paper.txt'}
    sources|={str(p.relative_to(ROOT)) for p in (ROOT/'analysis').glob('*cardiac_niche*') if p.is_file()}
    sources|={'analysis/revisit_vascular_coordination.py','analysis/cardiac_spatial_models.py','analysis/cardiac_spatial_pipeline.py','analysis/cardiac_spatial_pipeline.json','analysis/requirements-cardiac-spatial.txt','analysis/requirements-cardiac-niche-dynamics.txt','analysis/requirements-cardiac-wormhole.txt'}
    sources|={str(UP/stage/'proximity50_graph.npz') for stage in STAGES}|{str(UP/stage/'contact_graph.npz') for stage in STAGES}|{str(UP/stage/'proximity100_graph.npz') for stage in STAGES}
    save(pd.DataFrame([{'path':str(p),'bytes':(ROOT/p).stat().st_size,'sha256':sha(ROOT/p)} for p in sorted(sources)]),'source_checksums.tsv')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','validation.json','execution_receipts.json','runtime_provenance.json']]
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]),'file_manifest.tsv')
    dump('validation.json',{'status':'PASS','stages':checks,'configuration_sha256':sha(ROOT/'analysis/cardiac_niche_dynamics.json'),
        'matched_controls_exclude_all21_targets':True,'native_polygon_distances_independently_recomputed':True,'ring_BFS_independent_check':True,
        'paired_P7_coefficient_independently_residualized':True,'Wormhole_parameters_saved_and_actual_training_loss_verified':True,
        'Wormhole_heldout_metric_gate_recomputed':True,'exact_empirical_OT_independent_dense_LP_check':True,'source_capture_and_selection_covariates_checked':True,
        'PCA_fit_and_spatial_sources_do_not_overlap_holdout':True,'focal_maturation_genes_excluded_from_state_and_normalization':True,
        'figures':len(pngs),'PDF_pages':13,'documentation_links_checked':links,'scope':'Numerical/geometry/provenance checks; no claim of validated cell boundaries, independent replication, or biological causality'})
    print('VALIDATION PASS',flush=True)

if __name__=='__main__':validate()
