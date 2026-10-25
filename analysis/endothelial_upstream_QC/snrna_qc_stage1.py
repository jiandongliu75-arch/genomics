"""Library-specific QC and Scrublet; prepare singlet counts for DecontX."""
import os
from pathlib import Path
ROOT = Path('/workspace/heart-atlas')
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'cache/matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(ROOT/'cache/numba'))
import gzip, json, argparse
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy import sparse
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sc.settings.n_jobs = 4
OUT = ROOT/'results/snrna_qc'
OUT.mkdir(parents=True, exist_ok=True)

def process(prefix):
    out = OUT/prefix
    out.mkdir(exist_ok=True)
    print('START', prefix, flush=True)
    with gzip.open(ROOT/'data/snrna'/f'{prefix}_matrix.mtx.gz', 'rb') as f:
        counts = sio.mmread(f).T.tocsr().astype(np.int32)
    features = pd.read_csv(ROOT/'data/snrna'/f'{prefix}_features.tsv.gz', sep='\t', header=None)
    barcodes = pd.read_csv(ROOT/'data/snrna'/f'{prefix}_barcodes.tsv.gz', header=None)[0].astype(str)
    var = pd.DataFrame({'gene_id':features[0].values, 'gene_symbol':features[1].values}, index=features[0].values)
    var['mt'] = var.gene_symbol.str.lower().str.startswith('mt-')
    var['ribo'] = var.gene_symbol.str.match(r'^Rp[sl]\d')
    var['hemoglobin'] = var.gene_symbol.str.match(r'^Hb[ab]-')
    ad = sc.AnnData(counts, obs=pd.DataFrame(index=barcodes.values), var=var)
    ad.obs_names = prefix + ':' + ad.obs_names
    ad.obs['barcode'] = barcodes.values
    ad.obs['sample'] = prefix
    ad.obs['age'] = prefix.split('_')[1]
    ad.obs['replicate'] = prefix.split('_')[2]
    sc.pp.calculate_qc_metrics(ad, qc_vars=['mt','ribo','hemoglobin'], percent_top=None, inplace=True)
    ad.obs['low_genes'] = ad.obs.n_genes_by_counts < 750
    ad.obs['high_genes'] = ad.obs.n_genes_by_counts > 7500
    ad.obs['high_mito'] = ad.obs.pct_counts_mt > 5
    ad.obs['low_umis'] = ad.obs.total_counts < 1000
    for col in ['total_counts','n_genes_by_counts']:
        values = np.log1p(ad.obs[col].to_numpy())
        median = np.median(values)
        mad = np.median(np.abs(values - median))
        ad.obs[col+'_3mad_flag'] = np.abs(values-median)>3*mad
    flags = ['low_genes','high_genes','high_mito','low_umis']
    ad.obs['pass_initial_qc'] = ~ad.obs[flags].any(axis=1)
    work = ad[ad.obs.pass_initial_qc].copy()
    expected_rate = min(0.15, max(0.02, 0.008*ad.n_obs/1000))
    print(prefix, 'initial QC', work.n_obs, '/', ad.n_obs, 'expected doublets', expected_rate, flush=True)
    sc.pp.scrublet(work, expected_doublet_rate=expected_rate, sim_doublet_ratio=2.0,
                   n_prin_comps=30, random_state=42, verbose=True)
    threshold = float(work.uns['scrublet']['threshold'])
    ad.obs['doublet_score'] = np.nan
    ad.obs.loc[work.obs_names, 'doublet_score'] = work.obs.doublet_score
    ad.obs['predicted_doublet'] = False
    ad.obs.loc[work.obs_names, 'predicted_doublet'] = work.obs.predicted_doublet.astype(bool)
    ad.obs['pass_qc_and_doublet'] = ad.obs.pass_initial_qc & ~ad.obs.predicted_doublet
    # Retain every QC-pass nucleus: scDblFinder provides primary doublet calls in R.
    singlets = work.copy()
    ad.obs.to_csv(out/'all_barcode_qc.csv.gz')
    # Reference clusters support DecontX; these are not final biological labels.
    norm = singlets.copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    sc.pp.log1p(norm)
    sc.pp.highly_variable_genes(norm, n_top_genes=2500, flavor='seurat')
    sc.pp.pca(norm, n_comps=30, random_state=42)
    sc.pp.neighbors(norm, n_neighbors=15, n_pcs=30, random_state=42)
    sc.tl.leiden(norm, resolution=0.6, key_added='decontx_cluster', random_state=42, flavor='igraph', n_iterations=2)
    singlets.obs['decontx_cluster'] = norm.obs['decontx_cluster']
    singlets.obsm['X_pca_uncorrected'] = norm.obsm['X_pca']
    singlets.uns['qc_method'] = {'initial_thresholds':{'min_genes':750,'max_genes':7500,'max_pct_mito':5,'min_umis':1000},
        'doublet_method':'Scanpy Scrublet','expected_doublet_rate':expected_rate,'threshold':threshold,'seed':42,
        'ambient_input':'Called nuclei only; no empty droplets. DecontX uses library-specific clusters.'}
    singlets.write_h5ad(out/'qc_pass_uncorrected.h5ad', compression='gzip')
    with gzip.open(out/'decontx_input.mtx.gz', 'wb', compresslevel=1) as f:
        sio.mmwrite(f, singlets.X.T, field='integer')
    pd.DataFrame({'cell_id':singlets.obs_names, 'cluster':singlets.obs.decontx_cluster.to_numpy()}).to_csv(out/'decontx_cells.tsv', sep='\t',index=False)
    summary = {'sample':prefix,'input_barcodes':ad.n_obs,'initial_qc_pass':work.n_obs,
        'scrublet_doublets_diagnostic_only':int(work.obs.predicted_doublet.sum()),'prepared_for_scdblfinder':singlets.n_obs,
        'scrublet_threshold':threshold,'scrublet_expected_doublet_rate':expected_rate,
        'scrublet_observed_doublet_rate':float(work.obs.predicted_doublet.mean()),
        'decontx_clusters':int(norm.obs.decontx_cluster.nunique()),
        'filter_flag_counts_overlapping':{x:int(ad.obs[x].sum()) for x in flags}}
    (out/'stage1_summary.json').write_text(json.dumps(summary, indent=2))
    fig, axes = plt.subplots(1,3,figsize=(13,3.8))
    axes[0].scatter(ad.obs.total_counts, ad.obs.n_genes_by_counts,s=1,c=ad.obs.pass_initial_qc.map({True:'steelblue',False:'lightgray'}),rasterized=True)
    axes[0].set(xscale='log',yscale='log',xlabel='UMIs',ylabel='Detected genes')
    axes[1].hist(ad.obs.pct_counts_mt,bins=100);axes[1].axvline(5,c='red');axes[1].set(xlabel='Mitochondrial UMIs (%)',ylabel='Nuclei')
    axes[2].hist(work.obs.doublet_score,bins=80,alpha=.7,label='Observed')
    axes[2].hist(work.uns['scrublet']['doublet_scores_sim'],bins=80,alpha=.4,label='Simulated')
    axes[2].axvline(threshold,c='red');axes[2].set(xlabel='Scrublet score',ylabel='Count');axes[2].legend()
    fig.suptitle(prefix);fig.tight_layout();fig.savefig(out/'qc_and_doublets.pdf');fig.savefig(out/'qc_and_doublets.png',dpi=140);plt.close(fig)
    print('DONE',json.dumps(summary),flush=True)

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--sample');args=parser.parse_args()
    prefixes=[args.sample] if args.sample else [p.name.replace('_matrix.mtx.gz','') for p in sorted((ROOT/'data/snrna').glob('*matrix.mtx.gz'))]
    for prefix in prefixes:
        if (OUT/prefix/'stage1_summary.json').exists():
            print('REUSE',prefix,flush=True);continue
        process(prefix)
