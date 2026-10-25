"""Independent checks of immutable inputs, actual geometry and model arithmetic."""
from cardiac_spatial_pipeline import *
from cardiac_spatial_models import design,exposures,OUTCOMES
from scipy import ndimage
from scipy.special import logsumexp
from scipy.stats import multivariate_normal
from scipy.optimize import linprog
from scipy.spatial.distance import cdist
from PIL import Image
import re, importlib.metadata

SCRIPTS=['cardiac_spatial_pipeline.py','cardiac_spatial_niches.py','cardiac_spatial_models.py',
         'cardiac_spatial_pipeline_qc.py','report_cardiac_spatial_pipeline.py',
         'validate_cardiac_spatial_pipeline.py','run_cardiac_spatial_pipeline.py','cardiac_spatial_pipeline.json','requirements-cardiac-spatial.txt']

def validate():
    rng=np.random.default_rng(1541);checks=[];nichemeta=json.loads((OUT/'niche_methods.json').read_text())
    coef=pd.read_csv(OUT/'joint_association_coefficients.tsv.gz',sep='\t')
    sourcefiles=set(json.loads((OUT/'source_files.json').read_text()))
    sourcefiles|={'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',
                 'results/xenium_ventricular_mapping/reference_metadata.tsv.gz','results/xenium_ventricular_mapping/reference_shared_panel_counts.npz',
                 'results/ventricular_spatial_resolution/fine_group_evidence_summary.tsv'}
    prior=pd.read_csv(ROOT/'results/xenium_nuclear_anchored_segmentation/file_manifest.tsv',sep='\t').set_index('path')
    for stage in STAGES:
        print('VALIDATE_PIPELINE',stage,flush=True)
        f=frame(stage);old=read(f'results/whole_heart_spatial_atlas/{stage}/all_cell_annotations.tsv.gz')
        assert f.index.is_unique and f.stage_cell_id.is_unique
        np.testing.assert_array_equal(f.index,old.index)
        np.testing.assert_array_equal(f.atlas_identity,old.atlas_identity)
        np.testing.assert_array_equal(f.ventricular_CM_strict,old.ventricular_CM_strict)
        np.testing.assert_allclose(f[['x_centroid','y_centroid']],old[['x_centroid','y_centroid']],atol=1e-5)
        assert f.stage_cell_id.eq(stage+'::'+f.index.astype(str)).all()
        ids=f.index[f.ventricular_CM_strict];pick=rng.choice(len(ids),64,replace=False);ids=ids[pick]
        b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        gi=b.var_names.to_numpy(dtype=str);ci=b.obs_names.get_indexer(ids);x=b.layers['high_confidence'][ci].toarray()
        logs=np.log1p(x*1e4/np.maximum(x.sum(1),1)[:,None])
        np.testing.assert_allclose(logs[:,np.isin(gi,MAT.gene)].mean(1),f.loc[ids,'CM_maturation_core'],rtol=3e-6)
        for cat,z in MAT.groupby('figure_category'):
            np.testing.assert_allclose(logs[:,np.isin(gi,z.gene)].mean(1),f.loc[ids,'CM_'+cat],rtol=3e-6)
        assert not any(set(v['genes'])&set(MAT.gene) for v in CFG['programs'].values())
        for name,cfg in CFG['programs'].items():
            available=f.index[f[name].notna()]
            chosen=available[rng.choice(len(available),min(32,len(available)),replace=False)]
            index=b.obs_names.get_indexer(chosen);matrix=b.layers['nuclear_only' if cfg['view']=='nuclear' else 'high_confidence'][index].toarray()
            expected=np.log1p(matrix[:,np.isin(gi,cfg['genes'])]*1e4/np.maximum(matrix.sum(1),1)[:,None]).mean(1)
            np.testing.assert_allclose(expected,f.loc[chosen,name],rtol=3e-6,atol=1e-7)
            assert f.loc[chosen,'atlas_identity'].eq(cfg['lineage']).all()
        b.file.close();del b;gc.collect()
        with np.load(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/whole_section_segmentation.npz') as b:
            labels=b['cell_labels'];res=float(b['grid_um']);np.testing.assert_array_equal(b['cell_ids'],f.index)
        contacts=read(f"{CFG['output']}/{stage}/contact_edges.tsv.gz",index_col=None)
        g=graph(stage,'contact');assert g.shape==(len(f),len(f)) and g.diagonal().sum()==0 and (g-g.T).nnz==0
        assert g.nnz==2*len(contacts)
        # Independent bounding-box enumeration verifies complete boundaries of
        # sampled profile pairs, not merely centroid proximity.
        objects=ndimage.find_objects(labels,max_label=len(f))
        chosen=rng.choice(len(contacts),32,replace=False)
        for _,r in contacts.iloc[chosen].iterrows():
            i,j=int(r.profile_i),int(r.profile_j);si,sj=objects[i],objects[j]
            ys=slice(min(si[0].start,sj[0].start),max(si[0].stop,sj[0].stop));xs=slice(min(si[1].start,sj[1].start),max(si[1].stop,sj[1].stop))
            z=labels[ys,xs];total=0
            for a,c in [(z[:,:-1],z[:,1:]),(z[:-1],z[1:])]:total+=np.count_nonzero(((a==i+1)&(c==j+1))|((a==j+1)&(c==i+1)))
            np.testing.assert_allclose(total*res,r.shared_boundary_um,atol=1e-7)
            assert g[i,j]==1 and total*res>=CFG['contact_minimum_um']-1e-7
        edgeaudit=read(f"{CFG['output']}/{stage}/proximity_edge_audit_sample.tsv.gz",index_col=None)
        xy=f[['x_centroid','y_centroid']].to_numpy()
        edge_chosen=rng.choice(len(edgeaudit),128,replace=False)
        for _,r in edgeaudit.iloc[edge_chosen].iterrows():
            i,j=int(r.profile_i),int(r.profile_j);points=np.linspace(xy[i],xy[j],CFG['tissue_segment_samples'])
            support=[]
            for x,y in points:
                col=min(labels.shape[1]-1,max(0,int(round(x/res))));row=min(labels.shape[0]-1,max(0,int(round(y/res))))
                support.append(labels[row,col]>0)
            np.testing.assert_allclose(np.mean(support),r.tissue_support,atol=1e-7)
        for radius in [30,50,100]:
            gp=graph(stage,f'proximity{radius}');assert gp.diagonal().sum()==0 and (gp-gp.T).nnz==0
            for _,r in edgeaudit.iloc[edge_chosen[:16]].iterrows():
                expected=r.centroid_distance_um<=radius and r.tissue_support>=CFG['tissue_minimum_fraction']
                assert bool(gp[int(r.profile_i),int(r.profile_j)])==expected
            del gp
        # Actual saved neighborhood means are reconstructed from source profiles.
        d=read(f"{CFG['output']}/{stage}/proximity50_CM_neighborhoods.tsv.gz")
        gp=graph(stage,'proximity50')
        for cid in ids[:32]:
            i=f.index.get_loc(cid);neighbors=gp.indices[gp.indptr[i]:gp.indptr[i+1]];z=f.iloc[neighbors]
            assert len(neighbors)==d.loc[cid,'neighbors']
            for p in ['EC_cycling','FB_fibrillar_matrix','FB_basement_matrix']:
                expected=z[p].mean();np.testing.assert_allclose(expected,d.loc[cid,'program__'+p],rtol=4e-6,atol=1e-7,equal_nan=True)
            resolved=z[z.atlas_resolved]
            for l in LINEAGES:np.testing.assert_allclose(resolved.atlas_identity.eq(l).mean(),d.loc[cid,'fraction__'+l],atol=1e-7)
        contactd=read(f"{CFG['output']}/{stage}/contact_CM_neighborhoods.tsv.gz")
        np.testing.assert_array_equal(d.common_contact_50um_cohort,d.geometry_eligible&contactd.geometry_eligible)
        # Primary coefficient arithmetic via independently solved standardized OLS.
        data=d.loc[d.geometry_eligible];c=f.loc[data.index];e=exposures(data)
        e=e[[key for key in e if e[key].std()>1e-8 and (e[key].abs()>1e-8).sum()>=50]]
        for level in ['technical','anatomy','CM_state']:
            base=design(c,data,level);v=np.column_stack([base,e]);sd=v.std(0);use=sd>1e-8
            xx=np.column_stack([np.ones(len(v)),(v[:,use]-v[:,use].mean(0))/sd[use]])
            yy=c.CM_maturation_core.to_numpy();yy=(yy-yy.mean())/yy.std()
            estimate=np.linalg.lstsq(xx,yy,rcond=1e-9)[0];used=np.flatnonzero(use)
            for key in ['fraction__Fibroblast','program__EC_cycling','program__FB_basement_matrix']:
                pos=list(base)+list(e);j=pos.index(key);k=int(np.flatnonzero(used==j)[0])+1
                saved=coef.loc[coef.stage.eq(stage)&coef.geometry.eq('proximity50')&coef.adjustment.eq(level)&coef.outcome.eq('CM_maturation_core')&coef.exposure.eq(key)]
                assert len(saved)==1;np.testing.assert_allclose(estimate[k],saved.iloc[0].standardized_coefficient,atol=1e-7,rtol=3e-5)
        # Niche features are pure graph-derived composition, never outcome genes.
        with np.load(OUT/stage/'outcome_blind_niche_features.npz') as b:feat=b['features'];columns=b['columns'];eligible=b['eligible']
        assert feat.shape==(len(f),2*len(LINEAGES))
        assert not any(any(word in str(col) for word in ['maturation','Myl2','cycling','Slit2']) for col in columns)
        degree=np.asarray(g.sum(1)).ravel()
        for i in rng.choice(np.flatnonzero(eligible),32,replace=False):
            first=np.zeros(len(LINEAGES));second=np.zeros(len(LINEAGES));return_mass=0
            neighbors=g.indices[g.indptr[i]:g.indptr[i+1]]
            for j in neighbors:
                if f.atlas_identity.iloc[j] in LINEAGES:first[LINEAGES.index(f.atlas_identity.iloc[j])]+=1/degree[i]
                for k in g.indices[g.indptr[j]:g.indptr[j+1]]:
                    w=1/degree[i]/degree[j]
                    if k==i:return_mass+=w;continue
                    if f.atlas_identity.iloc[k] in LINEAGES:second[LINEAGES.index(f.atlas_identity.iloc[k])]+=w
            first/=first.sum();second/=second.sum()
            np.testing.assert_allclose(feat[i],np.r_[first,second],rtol=4e-6,atol=1e-7)
        with np.load(OUT/'niche_model_parameters.npz') as b:mu=b['means'];cov=b['covariances'];weights=b['weights'];sm=b['scaler_mean'];ss=b['scaler_scale']
        selected=rng.choice(len(f),128,replace=False);z=(np.sqrt(feat[selected])-sm)/ss
        logprob=np.column_stack([multivariate_normal.logpdf(z,mean=mu[k],cov=cov[k])+np.log(weights[k]) for k in range(len(mu))])
        posterior=np.exp(logprob-logsumexp(logprob,axis=1)[:,None])
        niche=read(f"{CFG['output']}/{stage}/niche_annotations.tsv.gz")
        np.testing.assert_allclose(posterior.max(1),niche.niche_maximum_posterior.iloc[selected],rtol=3e-6,atol=1e-7)
        ready=ad.read_h5ad(OUT/stage/'niche_analysis_ready.h5ad');assert ready.X is None
        assert not ready.uns['CellCharter_run'] and not ready.uns['Wasserstein_Wormhole_run']
        del ready
        # Stored point subsamples allow independent exact-distance reproduction.
        with np.load(OUT/stage/'nonCM_pointcloud_features.npz') as b:pointfeatures=b['features'];pointnames=b['feature_names']
        assert not any('CM_' in name for name in pointnames)
        with np.load(OUT/stage/'transport_pilot_point_indices.npz') as b:ia=b['point_indices_a'];ib=b['point_indices_b'];projection=b['projections']
        table=pd.read_csv(OUT/'transport_distance_pilot.tsv.gz',sep='\t');table=table.loc[table.stage.eq(stage)].reset_index(drop=True)
        for i in [0,3,7]:
            a=pointfeatures[ia[i][ia[i]>=0]];b=pointfeatures[ib[i][ib[i]>=0]]
            qa=np.quantile(a@projection,np.linspace(0,1,32),axis=0);qb=np.quantile(b@projection,np.linspace(0,1,32),axis=0)
            np.testing.assert_allclose(np.sqrt(np.mean((qa-qb)**2)),table.loc[i,'sliced_W2'],rtol=3e-6)
            na,nb=len(a),len(b);constraints=sparse.vstack([sparse.kron(sparse.eye(na),np.ones((1,nb))),sparse.kron(np.ones((1,na)),sparse.eye(nb))]).tocsr()
            result=linprog(cdist(a,b,'sqeuclidean').ravel(),A_eq=constraints,b_eq=np.r_[np.full(na,1/na),np.full(nb,1/nb)],bounds=(0,None),method='highs')
            assert result.success;np.testing.assert_allclose(np.sqrt(max(result.fun,0)),table.loc[i,'exact_empirical_W2'],rtol=3e-6,atol=1e-7)
        sourcefiles.add(f'results/xenium_nuclear_anchored_segmentation/{stage}/segmentation_profile_audit.tsv.gz')
        checks.append({'stage':stage,'profiles':len(f),'primary_CMs':len(data),'contact_pairs_independently_counted':32,
                       'barrier_samples_independently_checked':128,'niche_features_independently_checked':32,'niche_posteriors_reconstructed':128})
        del labels,objects,g,gp,f;gc.collect()
    assert sum(r['profiles'] for r in checks)==408054
    folds=pd.read_csv(OUT/'spatial_prediction_folds.tsv',sep='\t');assert (folds.minimum_train_test_separation_um>=CFG['prediction_buffer_um']-1e-6).all()
    folds=pd.read_csv(OUT/'nonlinear_prediction_folds.tsv',sep='\t');assert (folds.minimum_train_test_separation_um>=CFG['prediction_buffer_um']-1e-6).all()
    assert CFG['prediction_buffer_um']>=2*CFG['proximity_primary_um']
    sources=[]
    for rel in sorted(sourcefiles):
        h=digest(ROOT/rel)
        if rel in prior.index:assert h==prior.loc[rel,'sha256'],rel
        sources.append({'path':rel,'sha256':h,'bytes':(ROOT/rel).stat().st_size})
    save(pd.DataFrame(sources),'source_checksums.tsv')
    info=subprocess.check_output(['pdfinfo',str(OUT/'cardiac_spatial_pipeline_report.pdf')],text=True)
    pages=int(re.search(r'^Pages:\s+(\d+)',info,re.MULTILINE).group(1));assert pages==15
    for p in OUT.glob('*.png'):
        with Image.open(p) as im:im.verify()
    dump('configuration.json',CFG)
    dump('validation.json',{'status':'PASS','configuration_sha256':digest(CONFIG_PATH),'stages':checks,'PDF_pages':pages,
        'all_original_profiles_ids_coordinates_and_broad_labels_preserved':True,'source_segmentation_count_hashes_match_prior_manifest':True,
        'primary_outcomes_and_measured_program_scores_reconstructed':True,'contact_lengths_and_tissue_barriers_independently_verified':True,
        'joint_model_arithmetic_verified':True,'niche_features_and_posteriors_independently_verified':True,
        'prediction_buffer_verified':True,'transport_pilot_exact_distances_reproduced':True,
        'CellCharter_and_Wormhole_not_run_in_primary_pipeline':True,
        'biological_limits':'One section per age; arithmetic/stability checks do not establish anatomy, biological replication or causal influence'})
    rows=[]
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and p.name not in ['file_manifest.tsv','execution_receipts.json']:rows.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)})
    for script in SCRIPTS:
        p=ROOT/'analysis'/script;rows.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)})
    save(pd.DataFrame(rows),'file_manifest.tsv');print('PIPELINE_VALIDATION_PASS',checks,flush=True)

if __name__=='__main__':validate()
