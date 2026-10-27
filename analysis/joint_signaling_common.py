"""Frozen definitions for joint RNA/Xenium sender-resolved signaling analyses."""
from pathlib import Path
import os,json,hashlib
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/joint_RNA_Xenium_signaling'
STAGES=['P0','P7','P14','P21'];SEED=271028
PANEL=pd.read_csv(ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',sep='\t').gene.tolist()
COMPONENTS=json.loads((ROOT/'analysis/cardiac_spatial_reconstruction.json').read_text())['CM_components']
OUTCOME_GENES=set(PANEL)|set(COMPONENTS['cycling'])
YKEYS=['core_Z__maturation21','core_Z__metabolic','core_Z__sarcomere','core_Z__excitation_contraction','core_Z__cycling']
FB=['Postn','Gfpt2','Col8a1','Lsamp','Kirrel3']
EC=['Capillary-like','Endocardial-like','Arterial-like','Venous-like','Lymphatic-like']
GROUPS=['CM']+['EC:'+x for x in EC]+['EC:Unresolved']+['FB:'+x for x in FB]+['FB:Unresolved']+['Pericyte','Smooth muscle','Immune','Epicardial']
BROAD=['CM','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
INPUTS=set()
def read(p,**kw):
 p=Path(p);p=p if p.is_absolute() else ROOT/p;INPUTS.add(str(p));return pd.read_csv(p,sep='\t',**kw)
def save(d,name,index=False):
 OUT.mkdir(exist_ok=True,parents=True);p=OUT/name;p.parent.mkdir(exist_ok=True,parents=True)
 d.to_csv(p,sep='\t',index=index,na_rep='NA',float_format='%.8g')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4<<20),b''):h.update(b)
 return h.hexdigest()
def export_inputs():
 p=OUT/'input_paths.json';old=set(json.loads(p.read_text())) if p.exists() else set()
 p.write_text(json.dumps(sorted(old|INPUTS),indent=2)+'\n')
def canonical(x):return '_'.join(sorted(x.split('_')))
def complexes(x,genes,name):
 ix=pd.Index(genes).get_indexer(name.split('_'))
 return np.min(x[:,ix],axis=1) if (ix>=0).all() else np.zeros(len(x))
def load_database():
 path=ROOT/'analysis/vendor_commot/CellChatDB.ligrec.mouse.csv';INPUTS.add(str(path))
 d=pd.read_csv(path,index_col=0);d.columns=['ligand','receptor','pathway','signaling_type']
 d['pair_id']=[canonical(a)+'--'+canonical(b) for a,b in zip(d.ligand,d.receptor)]
 d=d.drop_duplicates('pair_id').reset_index(drop=True)
 alt=ROOT/'analysis/vendor_commot/CellPhoneDBv4.0.mouse.csv';INPUTS.add(str(alt))
 a=pd.read_csv(alt,index_col=0);pairs={canonical(x)+'--'+canonical(y) for x,y in zip(a.iloc[:,0],a.iloc[:,1])}
 d['CellPhoneDB_pair_agreement']=d.pair_id.isin(pairs)
 return d
def broad(group):return 'Endothelial' if group.startswith('EC:') else 'Fibroblast' if group.startswith('FB:') else group
def config():
 return {'stages':STAGES,'seed':SEED,'primary_database':'COMMOT 0.0.3 packaged CellChat mouse; alternate-database pair agreement annotated',
  'RNA_min_cells_per_library':30,'RNA_complex_detection_threshold':.02,'RNA_detection_definition':'All complex subunits have >=1 corrected molecule in the same nucleus', 'RNA_min_ligand_mean_logCP10k':.03,
  'RNA_min_receptor_mean_logCP10k':.03,'maximum_pairs_per_pathway':2,'maximum_global_pairs':64,
  'pair_selection_uses_spatial_outcomes':False,'primary_views':['core','native'],
  'secreted_radius_um':200,'ECM_contact_search_radius_um':50,'contact_range_is_verified_contact':False,
  'spatial_domains':4,'receiver_boundary_margin_um':200,'minimum_train_test_receiver_gap_um':400,
  'joint_solver':'Sparse collective OT component equivalent to COMMOT cot_sparse; not its four-solver ensemble',
  'joint_all_cell_types':True,'autocrine_self_edges':False,'cot_eps':.1,'cot_rho':10.,'cot_max_iterations':2000,'cot_stop':1e-3,
  'forest_trees':100,'forest_max_depth':8,'forest_min_leaf':30,'forest_max_features':.5,'maximum_training_CMs':12000,
  'ridge_alpha':100.,'ridge_minimum_molecular_scale':.05,'nested_hyperparameter_tuning':False,'feature_sets':['A','B','C','C_fine_context','D','E'],
  'new_biological_replicates':0,'new_segmentation':False,'unmeasured_genes_imputed_into_Xenium':False}
