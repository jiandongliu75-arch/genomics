"""Check shared cycling RNA using nuclear evidence and separated neighborhoods."""
from reconstruct_cardiac_spatial import ROOT,OUT,UP,CFG,STAGES,read,save,logmatrix
from test_cardiac_anatomy_mechanisms import fit,demean
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests

def run():
    rows=[];coverage=[];rng=np.random.default_rng(CFG['seed'])
    for stage in STAGES:
        print('NUCLEAR_CYCLING_CHECK',stage,flush=True)
        q=read(OUT/stage/'all_profile_observations.tsv.gz');c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');en=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz')
        a=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        x=a.layers['nuclear_only'][:].tocsr();genes=a.var_names.to_numpy(dtype=str);gi=np.isin(genes,CFG['CM_components']['cycling']);depth=np.asarray(x.sum(1)).ravel()
        cc=np.log1p(x[:,gi].toarray()*1e4/np.maximum(depth,1)[:,None]).mean(1);coherent=np.asarray((x[:,gi]>0).sum(1)).ravel()>=2
        # Primary source is nuclear-clean broad EC, not public High.Myo cells.
        source=q.atlas_identity.eq('Endothelial')&q.program_source_eligible&(depth>=20);si=np.flatnonzero(source.to_numpy())
        qc_source=cc[si];source_xy=q.iloc[si][['x_centroid','y_centroid']].to_numpy();sourceblocks=np.floor(source_xy/CFG['anatomy_block_um']).astype(int)
        _,sblock=np.unique(sourceblocks,axis=0,return_inverse=True);groups=[np.flatnonzero(sblock==k) for k in np.unique(sblock)]
        ci=q.index.get_indexer(c.index);cmdepth=depth[ci]
        c['CM_nuclear_cycling']=cc[ci];c['CM_nuclear_coherent_cycling']=coherent[ci].astype(float)
        g50=sparse.load_npz(UP/stage/'proximity50_graph.npz').tocsr()[ci][:,si];g30=sparse.load_npz(UP/stage/'proximity30_graph.npz').tocsr()[ci][:,si]
        for geometry,g in [('within50um',g50),('annulus30to50um',(g50-g30).tocsr())]:
            assert g.data.min(initial=0)>=0
            n=np.asarray(g.sum(1)).ravel();total=np.asarray(g@qc_source).ravel();ex=np.divide(total,n,out=np.full(len(n),np.nan),where=n>0)
            ecdepth=np.divide(np.asarray(g@np.log1p(depth[si])).ravel(),n,out=np.full(len(n),np.nan),where=n>0)
            valid=c.broader_ventricular_CM&(cmdepth>=20)&np.isfinite(ex)&c.distance_to_large_void_um.notna()
            d=c.loc[valid];e=pd.DataFrame({'QC_EC_nuclear_cycling':ex[valid]},index=d.index)
            b=pd.DataFrame(index=d.index);b['log_CM_nuclear_depth']=np.log1p(cmdepth[valid]);b['log_CM_area']=np.log1p(d.cell_area)
            b['log_EC_nuclear_depth']=ecdepth[valid];b['log_EC_neighbors']=np.log1p(n[valid]);b['void_fraction']=d.local_void_fraction
            b['log_void_distance']=np.log1p(d.distance_to_large_void_um);b['log_outer_distance']=np.log1p(d.distance_to_outer_um)
            b['capillary_area']=en.loc[valid,'author_capillary__nearby_cell_area_fraction_50um'];b['capillary_apposition']=en.loc[valid,'author_capillary__apposition_fraction_2um']
            fitrows=fit(d[['CM_nuclear_cycling','CM_nuclear_coherent_cycling']],b,e,d,within=True)
            # Conditional randomization of EC programs within local source blocks
            # tests finer co-localization, not signaling direction or causality.
            yy=d.CM_nuclear_cycling.to_numpy();yy=(yy-yy.mean())/max(yy.std(),1e-8);yr=demean(yy[:,None],d.anatomy_block_200um.to_numpy()).ravel()
            B=b.to_numpy();use=B.std(0)>1e-8;B=(B[:,use]-B[:,use].mean(0))/B[:,use].std(0);B=demean(B,d.anatomy_block_200um.to_numpy());B=np.column_stack([np.ones(len(d)),B]);pinv=np.linalg.pinv(B,rcond=1e-9);yr-=B@(pinv@yr)
            subg=g[valid.to_numpy()];nn=n[valid];null=[]
            for _ in range(199):
                v=qc_source.copy()
                for group in groups:v[group]=rng.permutation(v[group])
                xx=np.asarray(subg@v).ravel()/nn;xx=(xx-xx.mean())/max(xx.std(),1e-8);xx=demean(xx[:,None],d.anatomy_block_200um.to_numpy()).ravel();xx-=B@(pinv@xx)
                null.append(xx@yr/max(xx@xx,1e-12))
            observed=next(r['standardized_coefficient'] for r in fitrows if r['outcome']=='CM_nuclear_cycling')
            p=(1+np.sum(np.abs(null)>=abs(observed)))/(len(null)+1)
            for row in fitrows:rows.append({'stage':stage,'geometry':geometry,'local_source_program_permutation_p_for_continuous_score':p,**row})
            coverage.append({'stage':stage,'geometry':geometry,'nuclear_clean_EC_sources':len(si),'CMs_with_nuclear_evidence':len(d),
                'CM_fraction_two_cycling_genes_detected':d.CM_nuclear_coherent_cycling.mean(),'EC_fraction_two_cycling_genes_detected':coherent[si].mean(),
                'scope':'Observed nuclear RNA; zero/low nuclear information excluded only from this sensitivity, not removed from atlas'})
        save(c[['CM_nuclear_cycling','CM_nuclear_coherent_cycling']],f'{stage}/CM_nuclear_cycling_specificity.tsv.gz',index=True);a.file.close()
    rows=pd.DataFrame(rows);rows['within_section_FDR_BH']=multipletests(rows.approximate_block_p,method='fdr_bh')[1]
    save(rows,'nuclear_cycling_specificity_coefficients.tsv');save(pd.DataFrame(coverage),'nuclear_cycling_specificity_coverage.tsv')

if __name__=='__main__':run()
