"""Observed-gene EC-only spatial pooling with nuclear and boundary controls."""
from vascular_detection_common import *
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
import argparse

SPATIAL_PANELS={
    'Angpt2_single_gene':['Angpt2'],
    'Sprouting_guidance_observed_proxy':['Angpt2','Unc5b','Plxnd1','Itgb3','Lama4'],
    'EC_matrix_observed_proxy':['Mmp2','Itga5','Itgb3','Fn1'],
    'Cycling_observed':['Mki67','Top2a','Cenpf'],
}

def adaptive_pools(stage,obs,views,genes,folder):
    """Higher-coverage display pools; overlap is explicit and never replication."""
    xy=obs[['x_centroid','y_centroid']].to_numpy()
    distance,index=cKDTree(xy).query(xy,k=5,distance_upper_bound=75)
    keep=np.isfinite(distance).all(axis=1)
    indices=index[keep]
    design=sparse.csr_matrix((np.ones(indices.size),(np.repeat(np.arange(len(indices)),5),indices.ravel())),shape=(len(indices),len(obs)))
    assert np.asarray(design.sum(axis=1)).ravel().tolist()==[5]*len(indices)
    data=obs.loc[keep,['x_centroid','y_centroid']].copy()
    data['stage']=stage;data['source_EC_profiles']=5;data['farthest_source_distance_um']=distance[keep,-1]
    data['source_cell_ids']=[';'.join(obs.index[i]) for i in indices]
    assert data.farthest_source_distance_um.le(75).all()
    controls=json.loads((folder/'nuclear_clean_native_EC_50um_matched_controls.json').read_text())
    profiles=[];summaries=[];view_tables={}
    for view,counts in views.items():
        pooled=np.asarray(design@counts,dtype=np.float32);z=lognorm(pooled)
        d=data.copy();d['view']=view;d['gene_molecules']=pooled.sum(axis=1)
        cm=pooled[:,np.isin(genes,LINEAGES['Cardiomyocyte'])].sum(axis=1)
        d['CM_diagnostic_fraction']=cm/np.maximum(d.gene_molecules,1)
        for name,panel in SPATIAL_PANELS.items():
            js=pd.Index(genes).get_indexer(panel)
            score=[]
            for g,j in zip(panel,js):
                cs=pd.Index(genes).get_indexer(controls[g]);score.append(z[:,j]-z[:,cs].mean(axis=1))
            d[name+'__matched_control_excess']=np.mean(score,axis=0)
            d[name+'__observed_genes']=(pooled[:,js]>=1).sum(axis=1)
            d[name+'__molecules']=pooled[:,js].sum(axis=1)
            min_genes=1 if len(panel)==1 else 2
            summaries.append({'stage':stage,'view':view,'program':name,'source_EC_profiles':len(obs),'eligible_overlapping_pools':len(d),
                'EC_position_coverage_fraction':len(d)/len(obs),'pools_with_marker_support':int(d[name+'__observed_genes'].ge(min_genes).sum()),
                'pool_marker_support_fraction':float(d[name+'__observed_genes'].ge(min_genes).mean()),
                'interpretation':'Five-EC overlapping spatial pool within75um; display sensitivity gain at neighborhood resolution, not independent samples or cell-specific expression.'})
        view_tables[view]=d;profiles.append(d)
    np.testing.assert_array_equal(view_tables['nuclear'].gene_molecules,view_tables['nuclear_depth_control'].gene_molecules)
    np.testing.assert_array_equal(view_tables['boundary_core'].gene_molecules,view_tables['core_depth_control'].gene_molecules)
    output=view_tables['nuclear'].copy()
    for name,panel in SPATIAL_PANELS.items():
        score=output[name+'__matched_control_excess'];support=output[name+'__observed_genes'].ge(1 if len(panel)==1 else 2)
        output[name+'__higher_relative_enrichment_with_core_native_support']=support&score.gt(0)&score.ge(score.quantile(.75))&view_tables['native'][name+'__matched_control_excess'].gt(0)&view_tables['boundary_core'][name+'__matched_control_excess'].gt(0)
    pd.concat(profiles).to_csv(folder/'adaptive_5EC_75um_all_views.tsv.gz',sep='\t',na_rep='NA')
    output.to_csv(folder/'adaptive_5EC_75um_nuclear_proxy_maps.tsv.gz',sep='\t',na_rep='NA')
    save(pd.DataFrame(summaries),folder/'adaptive_pool_detection_coverage.tsv')

def tile_coherence(counts,genes,panels,cov,stage,size,view,cohort):
    z=lognorm(counts)
    rows=[]
    for name,panel in panels.items():
        js=pd.Index(genes).get_indexer(panel)
        eligible=(counts[:,js]>=1).sum(axis=0)>=10
        values=residual(z[:,js],cov)
        rs=[]
        for j in np.flatnonzero(eligible):
            other=np.flatnonzero(eligible&(np.arange(len(js))!=j))
            if len(other):rs.append(corr(values[:,j],values[:,other].mean(axis=1)))
        rows.append({'stage':stage,'tile_um':size,'cohort':cohort,'view':view,'program':name,'tiles':len(counts),
            'observed_genes':len(panel),'genes_in_ge10_tiles':int(eligible.sum()),
            'median_leave_one_gene_out_partial_Spearman':float(np.nanmedian(rs)) if len(rs) and np.isfinite(rs).any() else np.nan})
    return rows

def aggregate(stage,cohort,obs,views,genes,size,folder):
    xy=obs[['x_centroid','y_centroid']].to_numpy()
    grid=np.floor(xy/size).astype(int)
    tilekey=np.array([f'{a}_{b}' for a,b in grid])
    categories=pd.Categorical(tilekey)
    codes=categories.codes
    design=sparse.csr_matrix((np.ones(len(obs)),(codes,np.arange(len(obs)))),shape=(len(categories.categories),len(obs)))
    n=np.asarray(design.sum(axis=1)).ravel().astype(int)
    keep=n>=5
    valid=design[keep];keys=np.asarray(categories.categories)[keep]
    tilexy=np.array([[int(x) for x in k.split('_')] for k in keys])*size+size/2
    metadata=pd.DataFrame({'stage':stage,'cohort':cohort,'tile_um':size,'tile_id':stage+'__'+cohort+'__'+str(size)+'um__'+keys,
        'x_um':tilexy[:,0],'y_um':tilexy[:,1],'EC_profiles':n[keep]})
    membership=pd.DataFrame({'cell_id':obs.index,'stage':stage,'cohort':cohort,'tile_um':size,'grid_tile':tilekey,
        'eligible_tile_ge5_ECs':keep[codes]})
    pooled={v:np.asarray(valid@counts,dtype=np.float32) for v,counts in views.items()}
    assert np.array_equal(pooled['nuclear'].sum(axis=1),pooled['nuclear_depth_control'].sum(axis=1))
    assert np.array_equal(pooled['boundary_core'].sum(axis=1),pooled['core_depth_control'].sum(axis=1))
    baseline=lognorm(pooled['nuclear'])
    controls=select_controls(genes,baseline.mean(axis=0),(pooled['nuclear']>=1).mean(axis=0),20)
    write_json(folder/f'{cohort}_{size}um_matched_controls.json',controls)
    profile_rows=[];coverage=[];coherence=[];byview={}
    selected=np.asarray(valid.sum(axis=0)).ravel()>0
    for view,counts in pooled.items():
        d=metadata.copy();d['view']=view
        d['gene_molecules']=counts.sum(axis=1)
        cm=counts[:,np.isin(genes,LINEAGES['Cardiomyocyte'])].sum(axis=1)
        d['CM_diagnostic_fraction']=cm/np.maximum(d.gene_molecules,1)
        z=lognorm(counts)
        cov=np.column_stack([np.log1p(d.gene_molecules),np.log1p(d.EC_profiles),d.CM_diagnostic_fraction])
        for name,panel in SPATIAL_PANELS.items():
            js=pd.Index(genes).get_indexer(panel);assert (js>=0).all()
            detected=(counts[:,js]>=1).sum(axis=1)
            excess=[]
            for g,j in zip(panel,js):
                cj=pd.Index(genes).get_indexer(controls[g]);assert (cj>=0).all()
                excess.append(z[:,j]-z[:,cj].mean(axis=1))
            d[name+'__molecules']=counts[:,js].sum(axis=1)
            d[name+'__observed_genes']=detected
            d[name+'__mean_log1pCP10k']=z[:,js].mean(axis=1)
            d[name+'__matched_control_excess']=np.mean(excess,axis=0)
            cells=views[view][selected][:,js]
            min_genes=1 if len(panel)==1 else 2
            coverage.append({'stage':stage,'cohort':cohort,'tile_um':size,'view':view,'program':name,
                'eligible_tiles':len(d),'EC_profiles_in_eligible_tiles':int(selected.sum()),'minimum_marker_genes':min_genes,
                'individual_ECs_with_marker_support':int(((cells>=1).sum(axis=1)>=min_genes).sum()),
                'individual_EC_marker_support_fraction':float(((cells>=1).sum(axis=1)>=min_genes).mean()),
                'tiles_with_marker_support':int((detected>=min_genes).sum()),'tile_marker_support_fraction':float((detected>=min_genes).mean()),
                'limit':'Different units: pooled neighborhoods versus individual cells; increased coverage is not validated cell-level sensitivity.'})
        profile_rows.append(d);byview[view]=d
        coherence+=tile_coherence(counts,genes,SPATIAL_PANELS,cov,stage,size,view,cohort)
    consistency=[]
    if size==50:
        for name,panel in SPATIAL_PANELS.items():
            s=byview['nuclear'][name+'__matched_control_excess'].to_numpy()
            min_genes=1 if len(panel)==1 else 2
            for v in ['native','refined_cell','boundary_core','nuclear_depth_control','core_depth_control']:
                b=byview[v][name+'__matched_control_excess'].to_numpy()
                support=(byview['nuclear'][name+'__observed_genes']>=min_genes).to_numpy()&(byview[v][name+'__observed_genes']>=min_genes).to_numpy()
                cov=np.column_stack([np.log1p(byview['nuclear'].gene_molecules),metadata.EC_profiles,
                    byview['nuclear'].CM_diagnostic_fraction,byview[v].CM_diagnostic_fraction])
                a_top=(s>=np.quantile(s,.75))&support&(s>0)
                b_top=(b>=np.quantile(b,.75))&support&(b>0)
                consistency.append({'stage':stage,'cohort':cohort,'program':name,'reference_view':'nuclear','comparison_view':v,'tiles':len(s),
                    'Spearman':corr(s,b),'depth_density_CM_adjusted_Spearman':corr(residual(s,cov),residual(b,cov)),
                    'positive_topquartile_overlap_Jaccard':float((a_top&b_top).sum()/max((a_top|b_top).sum(),1)),
                    'both_views_have_multiple_observed_genes':int(support.sum())})
        # Label relative enrichment only; do not assert active angiogenesis or tip identity.
        output=byview['nuclear'].copy()
        for name,panel in SPATIAL_PANELS.items():
            core=byview['boundary_core'];native=byview['native'];score=output[name+'__matched_control_excess']
            support=(output[name+'__observed_genes']>=(1 if len(panel)==1 else 2))
            output[name+'__higher_relative_enrichment_with_core_native_support']=support&score.gt(0)&score.ge(score.quantile(.75))&core[name+'__matched_control_excess'].gt(0)&native[name+'__matched_control_excess'].gt(0)
        output.to_csv(folder/f'{cohort}_50um_nuclear_proxy_maps.tsv.gz',sep='\t',index=False,na_rep='NA')
    pd.concat(profile_rows).to_csv(folder/f'{cohort}_{size}um_observed_profiles.tsv.gz',sep='\t',index=False,na_rep='NA')
    membership.to_csv(folder/f'{cohort}_{size}um_tile_membership.tsv.gz',sep='\t',index=False)
    return coverage,coherence,consistency

def main(stage):
    folder=DEST/'Xenium'/stage;folder.mkdir(parents=True,exist_ok=True)
    source=ad.read_h5ad(OUT/'Xenium'/stage/'EC_candidate_count_views.h5ad')
    genes=source.var_names.to_numpy(dtype=str)
    # Frozen identities precede all current outcome scoring. No target score sets this cohort.
    primary=(source.obs.native_EC_broad&source.obs.nuclear_broad_EC_clean&~source.obs.multiple_projected_nuclear_components_review&~source.obs.geometry_unresolved_review).to_numpy()
    obs=source.obs.loc[primary].copy()
    values={'native':source.layers['native'][primary].toarray().astype(np.float32),
        'refined_cell':source.X[primary].toarray().astype(np.float32),
        'boundary_core':source.layers['high_confidence'][primary].toarray().astype(np.float32),
        'nuclear':source.layers['nuclear_only'][primary].toarray().astype(np.float32),
        'core_depth_control':source.layers['depth_matched'][primary].toarray().astype(np.float32),
        'nuclear_depth_control':source.layers['nuclear_depth_matched'][primary].toarray().astype(np.float32)}
    del source
    coverage=[];coherence=[];consistency=[]
    for size in [25,50,100]:
        print('SPATIAL_POOL',stage,size,len(obs),flush=True)
        a,b,c=aggregate(stage,'nuclear_clean_native_EC',obs,values,genes,size,folder)
        coverage+=a;coherence+=b;consistency+=c
    concordant=obs.EC_broad_clean.to_numpy()
    if concordant.sum()>=50:
        a,b,c=aggregate(stage,'cell_and_nuclear_clean_EC',obs.loc[concordant],{k:v[concordant] for k,v in values.items()},genes,50,folder)
        coverage+=a;coherence+=b;consistency+=c
    adaptive_pools(stage,obs,values,genes,folder)
    save(pd.DataFrame(coverage),folder/'observed_detection_coverage.tsv')
    save(pd.DataFrame(coherence),folder/'observed_gene_coherence.tsv')
    save(pd.DataFrame(consistency),folder/'cross_view_proxy_consistency.tsv')
    coverage_rows=[]
    for name,panel in PROGRAMS.items():
        for g in panel:coverage_rows.append({'program':name,'gene':g,'Xenium_measured':g in genes})
    save(pd.DataFrame(coverage_rows),folder/'panel_coverage.tsv')
    write_json(folder/'complete.json',{'status':'PASS','stage':stage,'nuclear_clean_native_EC_profiles':len(obs),
        'cell_and_nuclear_clean_EC_profiles':int(concordant.sum()),'primary_tile_um':50,'sensitivity_tile_um':[25,100],
        'minimum_EC_profiles_per_tile':5,'source_counts_and_coordinates_unchanged':True,
        'gene_measurements_added_or_imputed':False,'original_sprouting_core_observed':['Angpt2'],
        'observed_proxy_panels':SPATIAL_PANELS,
        'cohort':'Frozen native broad EC and nuclear-broad-clean evidence; no multiple projected nuclear components or unresolved geometry; independent of current outcome scores.',
        'controls':'Same profiles and non-overlapping physical tiles across native/refined/core/nuclear and exact molecule-depth controls; count equality checked after pooling.',
        'meaning':'Observed multi-gene guidance/remodeling or cycling proxies at neighborhood resolution, not full six-gene sprouting recovery, individual tip-cell identity or vessel-growth rate.'})
    print('SPATIAL_ENHANCEMENT_COMPLETE',stage,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=STAGES);a=p.parse_args();main(a.stage)
