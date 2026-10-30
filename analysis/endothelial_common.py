"""Definitions and auditable utilities for endothelial RNA/spatial QC."""
from pathlib import Path
import hashlib
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/endothelial-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/endothelial-numba')
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/endothelial_QC_subtypes'
RNA_INPUT = BASE / 'results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'
STAGES = ['P0', 'P7', 'P14', 'P21']
FAMILIES = ['Capillary-like', 'Endocardial-like', 'Arterial-like', 'Venous-like candidate', 'Lymphatic-like candidate']
UNRESOLVED = 'Endothelial subtype unresolved'
PANELS = {
    FAMILIES[0]: 'Rgcc Emcn Car4 Cd36 Gpihbp1 Fabp4 Dach1 Adgrl4 Cd300lg Btnl9 Slc28a2 Dipk2b'.split(),
    FAMILIES[1]: 'Npr3 Nfatc1 Nrxn3 Nrg1 Smoc1 Col23a1 Hs3st4 Enpp6 Kcnc2'.split(),
    FAMILIES[2]: 'Gja5 Gja4 Efnb2 Sox17 Dll4 Hey1 Sema3g Mannr Ptchd1'.split(),
    FAMILIES[3]: 'Nr2f2 Ephb4 Slc38a5 Ackr1 Vwf Selp'.split(),
    FAMILIES[4]: 'Prox1 Flt4 Lyve1 Pdpn Ccl21a Nrp2 Mmrn1 Tbx1 Reln'.split(),
}
ANCHORS = {
    FAMILIES[0]: 'Rgcc Emcn Car4 Gpihbp1 Cd300lg Btnl9'.split(),
    FAMILIES[1]: 'Npr3 Nfatc1 Nrxn3'.split(),
    FAMILIES[2]: 'Gja5 Efnb2 Sox17 Dll4 Sema3g'.split(),
    FAMILIES[3]: 'Nr2f2 Ephb4 Ackr1'.split(),
    FAMILIES[4]: 'Prox1 Flt4'.split(),
}
LINEAGES = {
    'Endothelial': 'Pecam1 Cdh5 Kdr Esam Tek Flt1 Tie1 Erg Fli1'.split(),
    'Cardiomyocyte': 'Myh6 Myh7 Tnni3 Myl2 Actc1 Ttn Ryr2 Mybpc3 Tnnt2'.split(),
    'Fibroblast': 'Dcn Lum Pdgfra Tcf21 Col1a1 Col1a2 Col3a1'.split(),
    'Mural': 'Rgs5 Cspg4 Kcnj8 Abcc9 Myh11 Cnn1 Tagln Pdgfrb Notch3'.split(),
    'Immune': 'Ptprc Csf1r Adgre1 Mertk Cd3d Cd3e Cd79a Ms4a1'.split(),
    'Glial': 'Sox10 Plp1 Mpz S100b'.split(),
}
CYCLING = 'Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split()
SEX = 'Xist Tsix Ddx3y Kdm5d Eif2s3y Uty Usp9y Zfy1 Zfy2 Sry Rbmy'.split()
WITHHOLD = list(set(sum(LINEAGES.values(), []) + SEX + ['Ankrd1', 'Creb5', 'Tlr4', 'Atf3']))
COLORS = dict(zip(FAMILIES + [UNRESOLVED], ['#009e73', '#e69f00', '#0072b2', '#cc79a7', '#8b5c38', '#bdbdbd']))

def save(d, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.7g')

def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + '\n')

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def lognorm(x):
    totals = np.asarray(x.sum(axis=1)).ravel()
    if sparse.issparse(x):
        z = (sparse.diags(1e4 / np.maximum(totals, 1)) @ x).tocsr()
        z.data = np.log1p(z.data)
        return z.astype(np.float32)
    return np.log1p(1e4 * x / np.maximum(totals, 1)[:, None]).astype(np.float32)

def lineage_evidence(x, genes, modality='Xenium'):
    """Operational review flags, not experimentally proven doublets."""
    total = np.asarray(x.sum(axis=1)).ravel()
    out = pd.DataFrame(index=np.arange(x.shape[0]))
    for name, panel in LINEAGES.items():
        a = x[:, np.isin(genes, panel)]
        n = np.asarray((a >= 1).sum(axis=1)).ravel()
        mol = np.asarray(a.sum(axis=1)).ravel()
        out[name + '_genes'] = n
        out[name + '_fraction'] = mol / np.maximum(total, 1)
        out[name + '_molecules'] = mol
    # Whole-transcriptome fractions have a much larger denominator than a 403-gene panel.
    ec_min, ec_fraction, foreign_fraction = (2, .001, .02) if modality == 'snRNA' else (3, .03, .10)
    out['EC_lineage_supported'] = (out.Endothelial_genes >= ec_min) & (out.Endothelial_fraction >= ec_fraction)
    mixed = np.zeros(len(out), bool)
    for name in LINEAGES:
        if name == 'Endothelial': continue
        strong = (out[name + '_genes'] >= 3) & (out[name + '_fraction'] >= foreign_fraction) & (out[name + '_molecules'] >= 10)
        out['strong_' + name + '_program'] = strong
        mixed |= strong.to_numpy()
    out['mixed_lineage_review'] = out.EC_lineage_supported & mixed
    return out

def marker_support(counts, genes, family, minimum_anchors=1):
    a = np.asarray((counts[:, np.isin(genes, PANELS[family])] >= 1).sum(axis=1)).ravel()
    b = np.asarray((counts[:, np.isin(genes, ANCHORS[family])] >= 1).sum(axis=1)).ravel()
    ok = (a >= 2) & (b >= minimum_anchors)
    if family == FAMILIES[4]:
        # Neither Lyve1 alone nor a non-specific lymphatic panel can establish identity.
        anchors = [g for g in ['Prox1', 'Flt4'] if g in genes]
        ok &= len(anchors) == 2
        if len(anchors) == 2: ok &= np.asarray((counts[:, np.isin(genes, anchors)] >= 1).sum(axis=1)).ravel() == 2
    return a, b, ok
