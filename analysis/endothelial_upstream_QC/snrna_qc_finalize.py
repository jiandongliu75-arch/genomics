"""Join scDblFinder/DecontX results, recheck QC, export counts and audit tables."""
from pathlib import Path
import os
ROOT=Path('/workspace/heart-atlas')
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'cache/matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'cache/numba'))
import argparse,gzip,json,re
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy import sparse
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'results/snrna_qc'
MARKERS=['Tnnt2','Actc1','Myh6','Myh7','Ttn','Col1a1','Dcn','Lum','Pdgfra','Pecam1','Cdh5','Kdr','Ptprc','Lyz2','Rgs5','Pdgfrb','Myh11','Tagln','Wt1','Upk3b','Sox10','Plp1']

def finalize(p):
    print('FINALIZE',p.name,flush=True)
    original=sc.read_h5ad(p/'qc_pass_uncorrected.h5ad')
    dm=pd.read_csv(p/'scdblfinder_metrics.csv',index_col='cell_id')
    am=pd.read_csv(p/'decontx_metrics.csv',index_col='cell_id')
    assert set(dm.index)==set(original.obs_names)
    assert set(am.index)==set(dm.index[dm.scdblfinder_class=='singlet'])
    if (p/'decontx_shape.txt').exists():
        shape=tuple(np.loadtxt(p/'decontx_shape.txt',dtype=int))
        values=np.fromfile(p/'decontx_values.f64',dtype='<f8')
        indices=np.fromfile(p/'decontx_indices.i32',dtype='<i4')
        indptr=np.fromfile(p/'decontx_indptr.i32',dtype='<i4')
        corrected=sparse.csc_matrix((values,indices,indptr),shape=shape).T.tocsr().astype(np.float32)
    else:
        with gzip.open(p/'decontx_corrected.mtx.gz','rb') as f:
            corrected=sio.mmread(f).T.tocsr().astype(np.float32)
    assert corrected.shape==(len(am), original.n_vars)
    assert np.isfinite(corrected.data).all() and (corrected.data>=0).all()
    data=original[am.index].copy()
    raw=data.X.copy()
    difference=corrected-raw.astype(np.float32)
    assert (difference.data<=1e-4).all(), 'Correction added counts or changed feature alignment'
    del difference
    before=np.asarray(raw.sum(axis=1)).ravel()
    after=np.asarray(corrected.sum(axis=1)).ravel()
    assert np.allclose(before,am.raw_umis)
    assert np.allclose(after,am.corrected_umis,rtol=1e-5)
    assert (after<=before+0.01).all()
    data.X=corrected
    data.layers['counts_raw']=raw
    data.obs['scrublet_predicted_doublet']=data.obs.pop('predicted_doublet')
    for col in dm.columns:data.obs[col]=dm.loc[data.obs_names,col]
    data.obs['decontx_contamination']=am.loc[data.obs_names,'decontx_contamination']
    data.obs['raw_total_counts']=data.obs.total_counts
    data.obs['raw_n_genes_by_counts']=data.obs.n_genes_by_counts
    data.obs['raw_pct_counts_mt']=data.obs.pct_counts_mt
    sc.pp.calculate_qc_metrics(data,qc_vars=['mt','ribo','hemoglobin'],percent_top=None,inplace=True)
    data.obs['high_contamination_flag']=data.obs.decontx_contamination>0.3
    data.obs['severe_contamination']=data.obs.decontx_contamination>0.5
    data.obs['post_ambient_low_genes']=data.obs.n_genes_by_counts<750
    data.obs['post_ambient_low_umis']=data.obs.total_counts<1000
    data.obs['post_ambient_high_mito']=data.obs.pct_counts_mt>5
    flags=['severe_contamination','post_ambient_low_genes','post_ambient_low_umis','post_ambient_high_mito']
    data.obs['pass_final_qc']=~data.obs[flags].any(axis=1)
    data.uns['qc_method']['primary_doublet_method']='scDblFinder, per-library, clusters supplied, default library-size prior, seed 42'
    data.uns['qc_method']['doublet_method']='scDblFinder'
    data.uns['qc_method']['diagnostic_doublet_method']='Scanpy Scrublet'
    for key in ['threshold','expected_doublet_rate']:
        if key in data.uns['qc_method']:
            data.uns['qc_method']['scrublet_'+key]=data.uns['qc_method'].pop(key)
    data.uns['qc_method']['ambient_seed']=12345
    data.uns['qc_method']['ambient_method']='DecontX with library-specific clusters; no empty-droplet background'
    data.uns['qc_method']['matrix_contents']='X: unnormalized DecontX-corrected counts (fractional); layers/counts_raw: original integer counts'
    data.uns['qc_method']['severe_contamination_exclusion_threshold']=0.5
    data.obs.to_csv(p/'all_singlet_ambient_qc.csv.gz')
    # The pre-final object allows users to revisit the post-correction thresholds.
    data.write_h5ad(p/'ambient_corrected_all_singlets.h5ad',compression='gzip')
    clean=data[data.obs.pass_final_qc].copy()
    clean.write_h5ad(p/'filtered.h5ad',compression='gzip')
    audit=pd.read_csv(p/'all_barcode_qc.csv.gz',index_col=0)
    audit=audit.rename(columns={'predicted_doublet':'scrublet_predicted_doublet'})
    for col in dm.columns:audit[col]=dm[col].reindex(audit.index)
    for col in ['decontx_contamination']+flags+['pass_final_qc']:
        audit[col]=data.obs[col].reindex(audit.index)
    audit['retained']=audit.index.isin(clean.obs_names)
    reasons=[]
    for _,r in audit.iterrows():
        rr=[f for f in ['low_genes','high_genes','high_mito','low_umis'] if r[f]]
        if r.get('scdblfinder_class')=='doublet':rr.append('scDblFinder_doublet')
        rr += [f for f in flags if pd.notna(r[f]) and bool(r[f])]
        reasons.append(';'.join(rr) if rr else 'retained')
    audit['exclusion_reason']=reasons
    assert ((audit.exclusion_reason=='retained')==audit.retained).all()
    audit.to_csv(p/'all_barcode_decisions.csv.gz')
    # Marker profiles expose what correction changes; cluster IDs are not cell-type labels.
    raw_n=raw.multiply(1e4/np.maximum(before,1)[:,None]).tocsr()
    cor_n=corrected.multiply(1e4/np.maximum(after,1)[:,None]).tocsr()
    profiles=[]
    for cluster in data.obs.decontx_cluster.unique():
        mask=(data.obs.decontx_cluster==cluster).to_numpy()
        for marker in MARKERS:
            loc=np.flatnonzero(data.var.gene_symbol.to_numpy()==marker)
            if not len(loc):continue
            profiles.append({'cluster':cluster,'gene':marker,'nuclei':int(mask.sum()),
                'raw_mean_per_10k':float(raw_n[mask][:,loc].sum(axis=1).mean()),
                'corrected_mean_per_10k':float(cor_n[mask][:,loc].sum(axis=1).mean())})
    profile_frame=pd.DataFrame(profiles)
    profile_frame.to_csv(p/'marker_profiles_before_after.csv',index=False)
    cluster_order=sorted(profile_frame.cluster.unique(),key=int)
    marker_order=[g for g in MARKERS if g in profile_frame.gene.unique()]
    raw_profile=profile_frame.pivot(index='cluster',columns='gene',values='raw_mean_per_10k').reindex(index=cluster_order,columns=marker_order)
    corrected_profile=profile_frame.pivot(index='cluster',columns='gene',values='corrected_mean_per_10k').reindex(index=cluster_order,columns=marker_order)
    scale=np.maximum(raw_profile.max(axis=0),corrected_profile.max(axis=0)).clip(lower=1e-9)
    fig,axes=plt.subplots(2,1,figsize=(13,8))
    for ax,mat,title in [(axes[0],raw_profile,'Original'),(axes[1],corrected_profile,'DecontX corrected')]:
        im=ax.imshow(mat/scale,aspect='auto',vmin=0,vmax=1,cmap='viridis')
        ax.set_xticks(range(len(marker_order)),marker_order,rotation=60,ha='right')
        ax.set_yticks(range(len(cluster_order)),cluster_order)
        ax.set(ylabel='Preliminary cluster',title=title)
    fig.suptitle(p.name+' | Within-gene scaled mean expression (shared scale)')
    fig.tight_layout();fig.savefig(p/'marker_profiles_before_after.pdf');fig.savefig(p/'marker_profiles_before_after.png',dpi=150);plt.close(fig)
    s=json.loads((p/'stage1_summary.json').read_text())
    fit_log=(ROOT/'logs'/f'{p.name}.decontx.log').read_text()
    convergence=re.findall(r'Completed iteration: (\d+) \| converge: ([0-9.eE+-]+)',fit_log)
    assert convergence, 'Missing DecontX convergence diagnostics'
    iterations,delta=convergence[-1]
    assert float(delta)<=0.001, 'DecontX did not converge: '+p.name
    s.update({'scdblfinder_doublets_removed':int((dm.scdblfinder_class=='doublet').sum()),
        'scdblfinder_doublet_fraction':float((dm.scdblfinder_class=='doublet').mean()),
        'singlets_ambient_corrected':data.n_obs,'final_nuclei':clean.n_obs,
        'post_ambient_removed':int((~data.obs.pass_final_qc).sum()),
        'post_ambient_flags_overlapping':{f:int(data.obs[f].sum()) for f in flags},
        'median_decontx_contamination':float(data.obs.decontx_contamination.median()),
        'p95_decontx_contamination':float(data.obs.decontx_contamination.quantile(.95)),
        'fraction_total_umis_removed_by_decontx':float(1-after.sum()/before.sum()),
        'decontx_iterations':int(iterations),'decontx_convergence':float(delta),
        'final_retention_fraction':clean.n_obs/s['input_barcodes']})
    (p/'final_summary.json').write_text(json.dumps(s,indent=2))
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    axes[0].hist(data.obs.decontx_contamination,bins=60);axes[0].axvline(.5,c='red');axes[0].set(xlabel='DecontX estimated contamination',ylabel='Singlets')
    axes[1].scatter(before,after,s=1,c=data.obs.decontx_contamination,rasterized=True,vmin=0,vmax=1)
    lim=max(before.max(),after.max());axes[1].plot([1,lim],[1,lim],c='gray');axes[1].set(xscale='log',yscale='log',xlabel='Original UMIs',ylabel='Corrected UMIs')
    axes[2].hist(dm.loc[dm.scdblfinder_class=='singlet','scdblfinder_score'],bins=70,alpha=.7,label='Singlet')
    axes[2].hist(dm.loc[dm.scdblfinder_class=='doublet','scdblfinder_score'],bins=70,alpha=.7,label='Doublet')
    axes[2].set(xlabel='scDblFinder score',ylabel='Nuclei');axes[2].legend()
    fig.suptitle(p.name);fig.tight_layout();fig.savefig(p/'ambient_and_primary_doublets.pdf');fig.savefig(p/'ambient_and_primary_doublets.png',dpi=150);plt.close(fig)
    print('FINALIZED',json.dumps(s),flush=True)

def merge():
    paths=sorted(OUT.glob('GSM*/filtered.h5ad'))
    assert len(paths)==8,'All eight libraries must complete first'
    summaries=[json.loads((p.parent/'final_summary.json').read_text()) for p in paths]
    frame=pd.DataFrame([{k:v for k,v in s.items() if not isinstance(v,dict)} for s in summaries])
    frame.to_csv(OUT/'qc_summary.tsv',sep='\t',index=False)
    from anndata.experimental import concat_on_disk
    concat_on_disk([str(p) for p in paths],str(OUT/'snrna_qc_filtered.h5ad'),join='inner',merge='same',uns_merge='same',max_loaded_elems=20_000_000)
    import h5py
    from anndata.io import write_elem
    with h5py.File(OUT/'snrna_qc_filtered.h5ad','r+') as handle:
        write_elem(handle.require_group('uns'),'qc_method',{
            'primary_doublet_method':'scDblFinder 1.24.10, separate libraries, seed 42',
            'diagnostic_doublet_method':'Scanpy Scrublet; calls not used for exclusion',
            'ambient_method':'DecontX 1.8.0, per-library clusters, no empty-droplet background',
            'ambient_seed':12345,'min_genes':750,'max_genes':7500,'min_umis':1000,
            'max_pct_mito':5.0,'max_estimated_contamination':0.5,
            'matrix_contents':'X: unnormalized corrected fractional counts; counts_raw: original integers',
            'per_library_parameters':'See qc_summary.tsv and per-library fitted metadata'})
    obj=ad.read_h5ad(OUT/'snrna_qc_filtered.h5ad',backed='r')
    assert obj.n_obs==int(frame.final_nuclei.sum()) and obj.obs_names.is_unique
    assert 'counts_raw' in obj.layers
    obj.file.close()
    # Per-barcode decisions include excluded nuclei and all filter reasons.
    pd.concat([pd.read_csv(p.parent/'all_barcode_decisions.csv.gz',index_col=0) for p in paths]).to_csv(OUT/'all_barcode_decisions.csv.gz')
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    labels=[s['sample'].split('_',1)[1] for s in summaries]
    inp=frame.input_barcodes.to_numpy();qc=frame.initial_qc_pass.to_numpy();sing=frame.singlets_ambient_corrected.to_numpy();fin=frame.final_nuclei.to_numpy()
    bottoms=np.zeros(8)
    for values,label,color in [(fin,'Retained','#287c8e'),(sing-fin,'Post-ambient QC','#e9c46a'),(qc-sing,'Doublets','#e76f51'),(inp-qc,'Initial QC','#bdbdbd')]:
        axes[0,0].bar(labels,values,bottom=bottoms,label=label,color=color);bottoms+=values
    axes[0,0].legend(fontsize=8);axes[0,0].set(ylabel='Nuclei',title='Library retention')
    axes[0,1].bar(labels,100*frame.scdblfinder_doublet_fraction);axes[0,1].set(ylabel='% of initial QC-pass nuclei',title='scDblFinder doublets')
    axes[1,0].bar(labels,100*frame.median_decontx_contamination);axes[1,0].set(ylabel='Median estimated contamination (%)',title='DecontX estimates')
    axes[1,1].bar(labels,100*frame.fraction_total_umis_removed_by_decontx);axes[1,1].set(ylabel='% UMIs removed from singlets',title='Ambient correction')
    for ax in axes.flat:ax.tick_params(axis='x',rotation=35)
    fig.tight_layout();fig.savefig(OUT/'qc_overview.pdf');fig.savefig(OUT/'qc_overview.png',dpi=160);plt.close(fig)
    print('MERGED',int(fin.sum()),'nuclei',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--sample');parser.add_argument('--merge',action='store_true');args=parser.parse_args()
    if args.merge:merge()
    else:
        paths=[OUT/args.sample] if args.sample else sorted(OUT.glob('GSM*'))
        for p in paths:
            if not (p/'final_summary.json').exists():finalize(p)
