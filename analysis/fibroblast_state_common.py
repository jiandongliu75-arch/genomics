"""Shared definitions for the focused fibroblast RNA and Xenium reanalysis."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/fibroblast_states'
AUDIT = ROOT / 'results/fibroblast_reassessment'
BASE = Path('/workspace/heart-atlas')
STAGES = ['P0', 'P7', 'P14', 'P21']
FAMILIES = ['Postn-associated', 'Col8a1-associated', 'Gfpt2-associated', 'Lsamp-associated', 'Kirrel3-associated']
UNRESOLVED = 'Unresolved fibroblast program'
PAPER_NAMES = {'Postn.Fib': FAMILIES[0], 'Col8a1.Fib': FAMILIES[1],
               'Gfpt2.Fib': FAMILIES[2], 'Lsamp.Fib': FAMILIES[3], 'valve.Fib': FAMILIES[4]}
ANCHORS = {FAMILIES[0]: ['Postn', 'Tnc'], FAMILIES[1]: ['Col8a1'],
           FAMILIES[2]: ['Gfpt2', 'Pi16'], FAMILIES[3]: ['Lsamp'], FAMILIES[4]: ['Kirrel3']}
SEX = 'Xist Tsix Ddx3y Kdm5d Eif2s3y Uty Usp9y Zfy1 Zfy2 Sry Rbmy'.split()
CM_CHECK = 'Myh6 Myh7 Tnni3 Myl2 Pln Actc1 Ttn Ryr2 Mybpc3 Tnnt2'.split()
MURAL_CHECK = 'Rgs5 Cspg4 Kcnj8 Abcc9 Myh11 Cnn1 Tagln'.split()
FB_CHECK = 'Dcn Lum Pdgfra Tcf21 Col1a1 Col1a2 Col3a1'.split()
CYCLING = 'Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split()
MATURATION = pd.read_csv(AUDIT / 'CM_maturation_Figure4b_panel.tsv', sep='\t').gene.tolist()
WITHHOLD = MATURATION + ['Ankrd1', 'Creb5', 'Tlr4', 'Atf3']
PROBE_MARKERS = pd.read_csv(AUDIT / 'paper_fibroblast_markers_and_coverage.tsv', sep='\t')
PANELS = {PAPER_NAMES[k]: list(v.gene) for k, v in PROBE_MARKERS.groupby('state') if k in PAPER_NAMES}
ANNOTATION_GENES = set(PROBE_MARKERS.gene)
_coverage = pd.read_csv(AUDIT / 'program_gene_coverage.tsv', sep='\t')
MATRIX_NAMES = ['Fibronectin', 'Fibrillar_collagen', 'Basement_membrane', 'Matricellular',
                'Elastic_fiber_associated', 'Collagen_processing_crosslinking', 'Hyaluronan_synthesis_turnover']
MATRIX_PANELS = {p: _coverage.loc[_coverage.panel.eq(p), 'gene'].tolist() for p in MATRIX_NAMES}
COLORS = dict(zip(FAMILIES + [UNRESOLVED], ['#e69f00', '#0072b2', '#009e73', '#cc79a7', '#8b5c38', '#bbbbbb']))


def save_tsv(data, path):
    data.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.7g')


def lognorm(counts, exclude=None):
    x = np.asarray(counts, dtype=np.float64)
    denominator = x.sum(axis=1) if exclude is None else x[:, ~exclude].sum(axis=1)
    return np.log1p(x * 1e4 / np.maximum(denominator, 1)[:, None])


def choose_markers(genes, counts, meanlog, fraction):
    """Post hoc panel evidence; an operational label is not an anatomical identity."""
    lookup = {g: i for i, g in enumerate(genes)}
    rows = []
    for fam, panel in PANELS.items():
        supported = [g for g in panel if g in lookup and fraction[lookup[g]] >= .1 and meanlog[lookup[g]] >= .1]
        anchors = [g for g in ANCHORS[fam] if g in supported]
        rows.append((fam, supported, anchors))
    return rows


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + '\n')
