"""Predeclared panels and conservative sensitivity tools; no count restoration."""
from endothelial_common import *
from scipy.stats import spearmanr, rankdata
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

DEST=ROOT/'results/vascular_program_detection'
PROGRAMS={
    'Original_sprouting_core':'Apln Esm1 Angpt2 Cxcr4 Pgf Adm'.split(),
    'Sprouting_guidance_support':'Apln Esm1 Angpt2 Cxcr4 Pgf Adm Unc5b Plxnd1 Kcne3 Fscn1 Nid2 Itgb3 Lama4'.split(),
    'EC_matrix_remodeling':'Mmp14 Mmp2 Adamts1 Adam15 Adam17 Plau Plaur Itga5 Itgb3 Fn1 Sparc Timp1 Timp2 Lgals3 Serpine1'.split(),
    'Cycling':CYCLING,
}
TARGET_GENES=list(dict.fromkeys(sum(PROGRAMS.values(),[])))
CAP='Dach1 Cd36 Rgcc Car4 Emcn Gpihbp1'.split()
LIMITS={
    'Original_sprouting_core':'Prespecified six-gene panel; Xenium only measures Angpt2, so its full score is not directly measurable spatially.',
    'Sprouting_guidance_support':'Supportive sprouting/guidance-associated panel, including multifunctional and inhibitory guidance components; not equivalent to a verified tip-cell identity or direct growth rate.',
    'EC_matrix_remodeling':'Protease/adhesion/matrix-associated RNA; shared with fibroblasts and other lineages, requiring EC identity and ownership checks. Not measured matrix degradation.',
    'Cycling':'Proliferation-associated RNA, separate from tip/stalk identity and completed cell division.',
}

def corr(a,b):
    a,b=np.asarray(a),np.asarray(b)
    keep=np.isfinite(a)&np.isfinite(b)
    if keep.sum()<10 or np.std(a[keep])<1e-10 or np.std(b[keep])<1e-10:return np.nan
    return float(spearmanr(a[keep],b[keep]).statistic)

def residual(y,cov):
    cov=np.asarray(cov,float)
    cov=(cov-cov.mean(axis=0))/np.maximum(cov.std(axis=0),1e-8)
    design=np.column_stack([np.ones(len(cov)),cov])
    return np.asarray(y)-design@np.linalg.lstsq(design,y,rcond=None)[0]

def select_controls(genes,mu,det,per_gene=20):
    excluded=set(TARGET_GENES+sum(LINEAGES.values(),[])+SEX)
    technical=pd.Series(genes).str.lower().str.startswith('mt-').to_numpy()|pd.Series(genes).str.match(r'^Rp[sl]\d').to_numpy()
    available=(det>=2/4597)&~np.isin(genes,list(excluded))&~technical
    pool=np.flatnonzero(available)
    features=np.column_stack([np.log1p(mu),np.log(np.maximum(det,1e-5))])
    features=(features-features.mean(axis=0))/np.maximum(features.std(axis=0),1e-8)
    controls={}
    for g in TARGET_GENES:
        if g not in genes:continue
        i=pd.Index(genes).get_loc(g)
        distance=((features[pool]-features[i])**2).sum(axis=1)
        controls[g]=genes[pool[np.argsort(distance,kind='stable')[:per_gene]]].tolist()
    return controls

def rank_excess(counts,genes,controls,maxrank=1500):
    """Top-rank enrichment minus expression/detection-matched controls.

    Fractional values are retained; <0.001 is only a numerical-dust floor.
    This score is an RNA enrichment diagnostic, not a pathway probability.
    """
    needed=list(dict.fromkeys(TARGET_GENES+sum(controls.values(),[])))
    needed=[g for g in needed if g in genes]
    pos=pd.Index(genes).get_indexer(needed)
    values=counts[:,pos].toarray()
    ranks=np.full(values.shape,maxrank+1,dtype=np.float32)
    for i in range(counts.shape[0]):
        a=counts.data[counts.indptr[i]:counts.indptr[i+1]]
        a=np.sort(a[a>=.001])
        valid=values[i]>=.001
        v=values[i,valid]
        ranks[i,valid]=len(a)-.5*(np.searchsorted(a,v,'left')+np.searchsorted(a,v,'right'))+.5
    enrichment=np.maximum(0,1-ranks/maxrank)
    idx={g:i for i,g in enumerate(needed)}
    rows={}
    for name,panel in PROGRAMS.items():
        gs=[g for g in panel if g in idx]
        gene_scores=[]
        for g in gs:
            cs=[idx[c] for c in controls[g]]
            gene_scores.append(enrichment[:,idx[g]]-enrichment[:,cs].mean(axis=1))
        rows[name]=np.mean(gene_scores,axis=0)
    return pd.DataFrame(rows)
