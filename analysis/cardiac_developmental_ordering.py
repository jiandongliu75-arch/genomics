"""Observed ventricular RNA ordering and explicitly conditional gap estimates.

Interpolated times are never counted as samples or spatial replicates. Outcomes
are withheld from the RNA ordering features; separate library analyses check
whether ordering survives independently sampled RNA hearts.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-reconstruction-mpl')
from reconstruct_cardiac_spatial import ROOT,OUT,CFG,STAGES,MAT,PROGS,read,save,dump,logmatrix
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.sparse.csgraph import connected_components,dijkstra
from scipy.sparse.linalg import eigsh
from scipy.interpolate import PchipInterpolator
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn import config_context
import gc,argparse

RNA='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'
MODULES={'CM_maturation21':MAT.gene.tolist(),**{'CM_'+k:v for k,v in CFG['CM_components'].items()},
    'CM_respiratory_chain_RNA_extended':'Acadm Acadvl Hadha Hadhb Sdha Sdhb Uqcrc1 Ndufs1 Ndufv1 Cox4i1 Cox5b Atp5f1b'.split(),
    'CM_sarcomere_extended':'Ttn Actn2 Mybpc3 Actc1 Tpm1 Myl3 Fhod3 Lmod2 Myom1 Myom2 Csrp3 Nebl'.split(),
    'CM_Ca_handling_extended':'Ryr2 Atp2a2 Casq2 Slc8a1 Cacna1c Cacnb2 Trdn'.split(),
    'EC_cycling':'Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split(),
    'EC_growth_compatible':'Flt4 Dll4 Kdr Nrp1 Angpt2'.split(),
    'EC_Notch_response':'Notch1 Hey1 Hey2 Hes1 Nrarp'.split(),
    'EC_sparse_sprouting':'Apln Esm1 Angpt2 Cxcr4 Pgf Adm'.split(),
    'Hypoxia_compatible':'Slc2a1 Pdk1 Ldha Bnip3 Vegfa Pgk1 Hk2 Adm Egln3'.split(),
    **{k:v['genes'] for k,v in PROGS.items() if k.startswith('FB_')}}

def modules_and_interpolation():
    print('RNA_LIBRARY_MODULES',flush=True);a=ad.read_h5ad(RNA,backed='r');obs=a.obs;symbols=a.var.gene_symbol.astype(str).to_numpy()
    genes=sorted(set(sum(MODULES.values(),[])));positions=np.flatnonzero(np.isin(symbols,genes));xx=a.X[:,positions].toarray()
    counts=np.zeros((len(obs),len(genes)),np.float32)
    for j,pos in enumerate(positions):counts[:,genes.index(symbols[pos])]+=xx[:,j]
    total=obs.total_counts.to_numpy(float);norm=np.log1p(counts*1e4/np.maximum(total,1)[:,None])
    primary_ec=read(ROOT/'results/endothelial_QC_subtypes/snRNA/EC_annotations.tsv.gz').index
    clean_fb=read(ROOT/'results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz').index
    scopes={'Ventricular_CM':obs.reviewed_identity.eq('Cardiomyocyte')&obs.proposed_subtype.eq('Ventricular-like'),
        'QC_Endothelial':obs.index.isin(primary_ec),'QC_reassessed_Fibroblast':obs.index.isin(clean_fb)}
    rows=[];gene_rows=[]
    for scope,take in scopes.items():
        for lib,z in obs.loc[take].groupby('sample',observed=True):
            ix=obs.index.get_indexer(z.index);stage=z.age.iloc[0];cpm=counts[ix].sum(0)/max(total[ix].sum(),1)*1e6
            for name,panel in MODULES.items():
                gi=[genes.index(g) for g in panel]
                rows.append({'scope':scope,'library':lib,'stage':stage,'day':int(stage[1:]),'program':name,'nuclei':len(ix),
                    'pseudobulk_logCPM_mean':np.log1p(cpm[gi]).mean(),'per_nucleus_logCP10k_mean':norm[ix][:,gi].mean(),
                    'fraction_at_least2_genes_corrected_ge05':np.mean((counts[ix][:,gi]>=.5).sum(1)>=2),'observed':True})
                for g,j in zip(panel,gi):gene_rows.append({'scope':scope,'library':lib,'stage':stage,'program':name,'gene':g,
                    'CPM':cpm[j],'detection_fraction':np.mean(counts[ix,j]>=.001)})
    rows=pd.DataFrame(rows);save(rows,'RNA_observed_library_programs.tsv');save(pd.DataFrame(gene_rows),'RNA_program_gene_support.tsv')
    # Eight library profiles support endpoint variation; bootstrap uncertainty
    # does not observe curvature inside missing time intervals.
    rng=np.random.default_rng(CFG['seed']);days=np.array([0,7,14,21]);query=np.array(CFG['inferred_days']);estimates=[]
    for (scope,program),z in rows.groupby(['scope','program']):
        if scope=='Ventricular_CM' and not (program.startswith('CM_') or program=='Hypoxia_compatible'):continue
        if scope=='QC_Endothelial' and not (program.startswith('EC_') or program=='Hypoxia_compatible'):continue
        if scope=='QC_reassessed_Fibroblast' and not program.startswith('FB_'):continue
        vv=[z.loc[z.day.eq(day),'pseudobulk_logCPM_mean'].to_numpy() for day in days];assert all(len(x)==2 for x in vv)
        center=np.array([v.mean() for v in vv]);linear=np.interp(query,days,center);pchip=PchipInterpolator(days,center)(query)
        boot=[]
        for _ in range(CFG['RNA_pseudobulk_bootstrap_draws']):
            y=np.array([rng.choice(v,2,replace=True).mean() for v in vv]);boot.append(np.stack([np.interp(query,days,y),PchipInterpolator(days,y)(query)]))
        boot=np.array(boot);low=np.quantile(boot,.025,axis=0).min(0);high=np.quantile(boot,.975,axis=0).max(0)
        for j,day in enumerate(query):estimates.append({'scope':scope,'program':program,'inferred_day':day,'linear_estimate':linear[j],
            'PCHIP_estimate':pchip[j],'model_shape_difference':abs(linear[j]-pchip[j]),'conditional_library_and_model_low':low[j],
            'conditional_library_and_model_high':high[j],'observed_at_this_day':False,'new_biological_replicates':0,
            'interpretation':'Transcript-program interpolation; uncertainty excludes unobserved curvature and unknown animal variability, no inferred spatial section or perfusion'})
    save(pd.DataFrame(estimates),'inferred_P3_P5_P10_programs.tsv')
    # Save only outcomes needed for independent feature-withholding checks.
    c=obs.loc[scopes['Ventricular_CM'],['sample','age','replicate','total_counts','proposed_subtype']].copy();ix=obs.index.get_indexer(c.index)
    for name,panel in MODULES.items():
        if name.startswith('CM_') or name=='Hypoxia_compatible':c[name]=norm[ix][:,[genes.index(g) for g in panel]].mean(1)
    save(c,'RNA_ventricular_observed_programs.tsv.gz',index=True);a.file.close()
    dump('inferred_day_interpretation.json',{'observed_days':[0,7,14,21],'inferred_days':CFG['inferred_days'],
        'independent_RNA_libraries_per_observed_age':2,'independent_spatial_sections_per_age_available':1,
        'new_observed_samples_generated':0,'new_biological_replication_generated':0,'latent_time_is_calendar_time':False,
        'not_identifiable':'True P3/P5/P10 distributions, spatial anatomy, temporal ordering of perfusion/contact effects, functional CM maturity or causal trajectories'})

def ordering():
    a=ad.read_h5ad(RNA,backed='r');obs=a.obs;take=obs.reviewed_identity.eq('Cardiomyocyte')&obs.proposed_subtype.eq('Ventricular-like')
    o=obs.loc[take].copy();ix=obs.index.get_indexer(o.index);x=logmatrix(a.X[ix].tocsr());genes=a.var.gene_symbol.astype(str).to_numpy()
    exclude=set(sum(MODULES.values(),[]))|set('Ccnb1 Ccnb2 Cdk1 Ube2c Cdc20 Cdc25c Bub1 Bub1b Aurka Aurkb Xist Ddx3y Eif2s3y Kdm5d Uty'.split())
    allowable=~np.isin(genes,list(exclude))&~pd.Series(genes).str.match(r'^(mt-|Rpl|Rps|Mcm[0-9])').to_numpy()
    summaries=[];edges=[];outputs=[];rng=np.random.default_rng(CFG['seed'])
    for name,mask in [('pooled',np.ones(len(o),bool)),('RNA_replicate1',o.replicate.astype(str).eq('1').to_numpy()),('RNA_replicate2',o.replicate.astype(str).eq('2').to_numpy())]:
        print('RNA_ORDERING',name,int(mask.sum()),flush=True)
        if mask.sum()<500:continue
        oo=o.loc[mask];indices=np.flatnonzero(mask);fit=[]
        for lib in oo['sample'].unique():
            avail=indices[o.iloc[indices]['sample'].eq(lib).to_numpy()];fit.extend(rng.choice(avail,min(800,len(avail)),replace=False))
        fit=np.array(fit);train=x[fit];mean=np.asarray(train.mean(0)).ravel();var=np.asarray(train.power(2).mean(0)).ravel()-mean**2
        dispersion=var/np.maximum(mean,.01);bins=pd.qcut(pd.Series(mean).rank(method='first'),20,labels=False).to_numpy()
        dscore=np.zeros(len(mean))
        for b in range(20):
            v=dispersion[bins==b];dscore[bins==b]=(v-np.median(v))/max(np.median(abs(v-np.median(v))),.01)
        dscore[~allowable]=-np.inf;selected=np.argsort(dscore)[-CFG['RNA_HVGs']:]
        assert not set(genes[selected])&exclude
        pca=PCA(n_components=CFG['RNA_PCs'],random_state=CFG['seed'],svd_solver='randomized').fit(x[fit][:,selected].toarray())
        pc=np.vstack([pca.transform(x[j:j+4096][:,selected].toarray()) for j in range(0,len(o),4096)])[indices]
        # A within-age library offset is a named sensitivity; it preserves each
        # age centroid, not a forced integration of all developmental stages.
        if name=='pooled':
            raw=pc.copy()
            for age in STAGES:
                stage=oo.age.eq(age).to_numpy();center=pc[stage].mean(0)
                for lib in oo.loc[stage,'sample'].unique():
                    m=oo['sample'].eq(lib).to_numpy();pc[m]-=pc[m].mean(0)-center
        with config_context(working_memory=128):
            nn=NearestNeighbors(n_neighbors=CFG['RNA_neighbors']+1,algorithm='brute',n_jobs=2).fit(pc);dist,neigh=nn.kneighbors(pc)
        ni=neigh[:,1:];dd=dist[:,1:];rr=np.repeat(np.arange(len(oo)),CFG['RNA_neighbors']);cc=ni.ravel();vv=dd.ravel()
        dg=sparse.csr_matrix((np.maximum(vv,1e-6),(rr,cc)),shape=(len(oo),len(oo)));dg=dg.maximum(dg.T).tocsr()
        nc,labels=connected_components(dg,directed=False);largest=np.argmax(np.bincount(labels));keep=labels==largest
        ages=oo.age.to_numpy();dcomp=[]
        for comp in np.unique(labels):dcomp.append({'component':int(comp),'nuclei':int((labels==comp).sum()),'ages':';'.join(sorted(set(ages[labels==comp])))})
        for sa in STAGES:
            for sb in STAGES:
                m=ages[rr]==sa;edges.append({'ordering':name,'from_stage':sa,'to_stage':sb,'directed_edges':int(np.sum(m&(ages[cc]==sb))),
                    'fraction_of_from_stage_edges':float(np.mean(ages[cc[m]]==sb)) if m.any() else np.nan})
        degree=np.asarray(dg.astype(bool).sum(1)).ravel();root_candidates=np.flatnonzero(keep&oo.age.eq('P0').to_numpy())
        root=root_candidates[np.argmin(((pc[root_candidates]-pc[root_candidates].mean(0))**2).sum(1))]
        geodesic=dijkstra(dg,directed=False,indices=root);geodesic[~np.isfinite(geodesic)]=np.nan
        scale=np.nanquantile(geodesic,.99);geo=geodesic/max(scale,1e-8)
        # Diffusion-distance alternative uses only the largest connected graph.
        sub=dg[keep][:,keep].tocsr();radius=np.maximum(dd[:,-1],1e-4)
        weights=np.exp(-(vv**2)/(radius[rr]*radius[cc]));wg=sparse.csr_matrix((weights,(rr,cc)),shape=dg.shape);wg=wg.maximum(wg.T)[keep][:,keep]
        deg=np.asarray(wg.sum(1)).ravel();T=sparse.diags(1/np.sqrt(np.maximum(deg,1e-8)))@wg@sparse.diags(1/np.sqrt(np.maximum(deg,1e-8)))
        eigval,eigvec=eigsh(T,k=8,which='LA',tol=1e-5);order=np.argsort(eigval)[::-1];eigval=eigval[order];eigvec=eigvec[:,order]
        use=(eigval<.99999)&(eigval>0);scaled=eigvec[:,use]*(eigval[use]/np.maximum(1-eigval[use],1e-5))
        root_sub=int(np.flatnonzero(np.flatnonzero(keep)==root)[0]);dif=np.sqrt(((scaled-scaled[root_sub])**2).sum(1));dpt=np.full(len(oo),np.nan);dpt[keep]=dif/max(np.quantile(dif,.99),1e-8)
        result=oo[['sample','age','replicate','total_counts']].copy();result['ordering']=name;result['graph_geodesic_order']=geo;result['diffusion_distance_order']=dpt;result['connected_to_P0_root']=keep
        result['observed_day']=oo.age.str[1:].astype(int);outputs.append(result)
        summaries.append({'ordering':name,'nuclei':len(oo),'graph_components':nc,'largest_component_nuclei':int(keep.sum()),
            'geodesic_age_Spearman':spearmanr(result.observed_day,geo,nan_policy='omit').statistic,
            'diffusion_age_Spearman':spearmanr(result.observed_day,dpt,nan_policy='omit').statistic,
            'depth_geodesic_Spearman':spearmanr(np.log1p(result.total_counts),geo,nan_policy='omit').statistic,
            'root_observed_age':'P0','root_id':oo.index[root],'HVGs':len(selected),'outcome_genes_withheld':True,
            'interpretation':'Cross-sectional developmental ordering; graph continuity and library replication exposed, not a calibrated day or observed lineage'})
        np.savez_compressed(OUT/(name+'_RNA_ordering_model.npz'),cell_ids=oo.index.to_numpy(dtype=str),PCs=pc.astype(np.float32),HVG_symbols=np.asarray(genes[selected],dtype=str),HVG_indices=selected,
            components=pca.components_,mean=pca.mean_,root_id=oo.index[root],eigenvalues=eigval)
        dump(name+'_graph_components.json',dcomp);del pc,dg,wg,T,sub;gc.collect()
    save(pd.concat(outputs),'RNA_cell_developmental_ordering.tsv.gz',index=True);save(pd.DataFrame(summaries),'RNA_ordering_assessment.tsv');save(pd.DataFrame(edges),'RNA_cross_age_graph_edges.tsv');a.file.close()

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['modules','ordering','all'],default='all');a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if a.step in ['modules','all']:modules_and_interpolation()
    if a.step in ['ordering','all']:ordering()

if __name__=='__main__':main()
