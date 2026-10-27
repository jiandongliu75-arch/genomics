"""Independent-library RNA eligibility, cell-source evidence and CM programs."""
import gc,json,time
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from anndata.io import sparse_dataset,read_elem
from joint_signaling_common import *

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=True,parents=True);db=load_database()
 measured=set(read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').gene)
 db['all_targets_measured']=db.apply(lambda r:set((r.ligand+'_'+r.receptor).split('_'))<=measured,axis=1)
 db['outcome_gene_overlap']=db.apply(lambda r:bool(set((r.ligand+'_'+r.receptor).split('_'))&OUTCOME_GENES),axis=1)
 candidates=db[db.all_targets_measured&~db.outcome_gene_overlap].copy()
 path=Path('/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');INPUTS.add(str(path))
 programs=json.loads((ROOT/'results/ventricular_microenvironment/program_definitions.json').read_text())['programs']
 response={name:programs[name]['genes'] for name in ['Oxidative_metabolism','Fatty_acid_oxidation','Sarcomere_mechanics','HIF_response','IEG_preparation']}
 response['Cell_cycle_full']=['Mki67','Top2a','Cenpf','Ccnb1','Cdk1','Ube2c']
 required=measured|set(sum(response.values(),[]))|{'Meis1','Smad2','Smad3','Smad7'}
 with h5py.File(path) as h:
  obs=read_elem(h['obs']);var=read_elem(h['var']);symbols=var.gene_symbol.astype(str).to_numpy()
  genes=np.array(sorted(required&set(symbols)));colmap=pd.Index(genes).get_indexer(symbols)
  keep=np.flatnonzero(colmap>=0);collapse=sparse.csr_matrix((np.ones(len(keep)),(keep,colmap[keep])),shape=(len(symbols),len(genes)))
  matrices={};totals={}
  for name,key in [('corrected','X'),('raw','layers/counts_raw')]:
   print('RNA_READ',name,flush=True);x=sparse_dataset(h[key])[:].tocsr()
   totals[name]=np.asarray(x.sum(1)).ravel().astype(float);matrices[name]=(x@collapse).toarray().astype(np.float32)
   del x;gc.collect()
 np.testing.assert_allclose(totals['corrected'],obs.total_counts,rtol=2e-5)
 meta=read('results/xenium_ventricular_mapping/reference_metadata.tsv.gz',index_col=0)
 assert obs.index.equals(meta.index)
 obs['source_group']='Exclude';obs['source_broad']='Exclude'
 cm=obs.reviewed_identity.eq('Cardiomyocyte');obs.loc[cm,'source_group']='CM'
 for name in ['Pericyte','Smooth muscle','Immune','Epicardial']:
  obs.loc[obs.reviewed_identity.eq(name),'source_group']=name
 fb=read('results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz',index_col=0)
 obs.loc[fb.index,'source_group']='FB:Unresolved'
 for family in FB:obs.loc[fb.index[fb.fibroblast_family.eq(family+'-associated')],'source_group']='FB:'+family
 ec=read('results/endothelial_QC_subtypes/snRNA/EC_annotations.tsv.gz',index_col=0)
 good=ec.index[ec.primary_EC_eligible];obs.loc[good,'source_group']='EC:Unresolved'
 ec_family=ec.EC_family.str.replace(' candidate','',regex=False)
 for family in EC:obs.loc[ec.index[ec.primary_EC_eligible&ec_family.eq(family)],'source_group']='EC:'+family
 obs['source_broad']=obs.source_group.map(broad)
 obs['receiver_ventricular']=cm&~obs.proposed_subtype.fillna('').str.contains('Atrial',case=False)
 obs['CM_state']=meta.mapping_family
 save(obs[['sample','age','source_group','source_broad','receiver_ventricular','CM_state','decontx_contamination']], 'RNA_cell_source_annotations.tsv.gz',index=True)
 coverage=obs.groupby(['age','sample','source_group'],observed=True).size().reset_index(name='nuclei');save(coverage,'RNA_source_coverage.tsv')
 # Full-genome denominators; no RNA values are inserted into Xenium observations.
 log={name:np.log1p(x*1e4/np.maximum(totals[name],1)[:,None]) for name,x in matrices.items()}
 species=sorted(set(candidates.ligand)|set(candidates.receptor))
 compound={name:np.column_stack([complexes(v,genes,s) for s in species]) for name,v in log.items()}
 detection={name:np.column_stack([complexes(v,genes,s)>=1 for s in species]) for name,v in matrices.items()}
 species_index={s:i for i,s in enumerate(species)}
 expr=[];lookup={}
 groups=['CM']+['Broad:'+x for x in BROAD if x!='CM']+[g for g in GROUPS if g!='CM']
 for stage in STAGES:
  for lib in sorted(obs.loc[obs.age.eq(stage),'sample'].unique()):
   for group in groups:
    groupmask=obs.source_broad.eq(group[6:]) if group.startswith('Broad:') else obs.receiver_ventricular if group=='CM' else obs.source_group.eq(group)
    ix=(obs.age.eq(stage)&obs['sample'].eq(lib)&groupmask).to_numpy()
    if not ix.any():continue
    for view in ['corrected','raw']:
     vals=compound[view][ix];det=detection[view][ix].mean(0);means=vals.mean(0)
     for j,s in enumerate(species):
      row={'stage':stage,'sample':lib,'source_group':group,'view':view,'species':s,'nuclei':int(ix.sum()),'mean_log1p_CP10k':means[j],'fraction_complex_detected':det[j]}
      expr.append(row);lookup[(stage,lib,group,view,s)]=(int(ix.sum()),float(means[j]),float(det[j]))
 save(pd.DataFrame(expr),'RNA_complex_expression_by_library.tsv.gz')
 support=[];summary=[]
 for stage in STAGES:
  libs=sorted(obs.loc[obs.age.eq(stage),'sample'].unique());assert len(libs)==2
  for _,pair in candidates.iterrows():
   recv=[lookup.get((stage,lib,'CM','corrected',pair.receptor),(0,0,0)) for lib in libs]
   recraw=[lookup.get((stage,lib,'CM','raw',pair.receptor),(0,0,0)) for lib in libs]
   rec_ok=all(n>=30 and mu>=.03 and frac>=.02 for n,mu,frac in recv)
   groups_ok=[];strength=[]
   for group in groups:
    ligand=[lookup.get((stage,lib,group,'corrected',pair.ligand),(0,0,0)) for lib in libs]
    rawlig=[lookup.get((stage,lib,group,'raw',pair.ligand),(0,0,0)) for lib in libs]
    lig_ok=all(n>=30 and mu>=.03 and frac>=.02 for n,mu,frac in ligand)
    raw_ok=all(n>=30 and mu>=.03 and frac>=.02 for n,mu,frac in rawlig+recraw)
    qualified=rec_ok and lig_ok
    effect=min([v[1] for v in ligand]+[v[1] for v in recv]) if qualified else 0
    support.append({'stage':stage,'pair_id':pair.pair_id,'ligand':pair.ligand,'receptor':pair.receptor,'pathway':pair.pathway,
      'source_group':group,'receiver_two_library_support':rec_ok,'sender_two_library_support':lig_ok,'joint_two_library_support':qualified,
      'raw_also_supported':raw_ok,'minimum_sender_library_nuclei':min(x[0] for x in ligand),
      'minimum_sender_detection':min(x[2] for x in ligand),'minimum_receiver_detection':min(x[2] for x in recv),
      'minimum_mean_expression':effect,'interpretation':'Replicated molecular availability, not demonstrated communication'})
    if qualified:groups_ok.append(group);strength.append(effect)
   summary.append({'stage':stage,'pair_id':pair.pair_id,'supported_source_groups':';'.join(groups_ok),
       'RNA_supported':bool(groups_ok),'RNA_support_strength':max(strength,default=0)})
 support=pd.DataFrame(support);summary=pd.DataFrame(summary);save(support,'RNA_pair_source_support.tsv.gz');save(summary,'RNA_pair_stage_support.tsv')
 reverse=[]
 for stage in STAGES:
  libs=sorted(obs.loc[obs.age.eq(stage),'sample'].unique())
  for ligand in ['Tgfb1','Tgfb2']:
   receptor='Tgfbr1_Tgfbr2'
   for family in ['Broad:Fibroblast']+['FB:'+x for x in FB]:
    src=[lookup.get((stage,l,'CM','corrected',ligand),(0,0,0)) for l in libs]
    dst=[lookup.get((stage,l,family,'corrected',receptor),(0,0,0)) for l in libs]
    reverse.append({'stage':stage,'ligand':ligand,'receptor':receptor,'sender':'CM','receiver':family,
       'minimum_sender_detection':min(x[2] for x in src),'minimum_receiver_detection':min(x[2] for x in dst),
       'minimum_receiver_library_nuclei':min(x[0] for x in dst),
       'two_library_availability_support':all(n>=30 and mu>=.03 and frac>=.02 for n,mu,frac in src+dst),
       'interpretation':'Reverse-direction molecular availability only; no evidence of a fate switch or mediated CM maturation'})
 save(pd.DataFrame(reverse),'RNA_CM_to_FB_TGF_availability.tsv')
 ranking=summary.groupby('pair_id').agg(supported_stages=('RNA_supported','sum'),support_strength=('RNA_support_strength','max'))
 pool=candidates.join(ranking,on='pair_id');pool=pool[pool.supported_stages>0].sort_values(['supported_stages','support_strength','pair_id'],ascending=[False,False,True])
 selected=pool.groupby('pathway',sort=False).head(2).copy()
 # Prespecified atlas/ECM candidates may be retained when supported; no spatial outcome is inspected.
 named={('Ptprm','Ptprm'),('Nampt','Insr'),('Ptn','Ptprz1'),('Igf2','Igf1r'),('Fn1','Itga5_Itgb1'),('Tgfb2','Tgfbr1_Tgfbr2')}
 selected=pd.concat([pool[pool.apply(lambda r:(r.ligand,r.receptor) in named,axis=1)],selected]).drop_duplicates('pair_id').head(64)
 selected=selected.sort_values(['pathway','pair_id']).reset_index(drop=True)
 selected['pair_index']=np.arange(len(selected));selected['radius_um']=np.where(selected.signaling_type.eq('Secreted Signaling'),200,50)
 selected['biochemical_relationship_not_validated_here']=True
 save(selected,'selected_RNA_supported_pairs.tsv')
 db['selected_for_spatial_analysis']=db.pair_id.isin(selected.pair_id);save(db,'database_assay_coverage_and_selection.tsv')
 # Independent CM response programs contextualize signaling; spatial outcomes remain directly measured.
 allLR=set('_'.join(candidates.ligand.tolist()+candidates.receptor.tolist()).split('_'))
 definitions=[];score_rows=[];response_scores={}
 for program,panel in response.items():
  use=sorted(set(panel)&set(genes));independent=sorted(set(use)-allLR-OUTCOME_GENES)
  definitions.append({'program':program,'full_genes':';'.join(use),'independent_genes':';'.join(independent),
      'measured_independent_genes':';'.join(sorted(set(independent)&measured)),'independent_gene_count':len(independent),
      'interpretation':'Canonical RNA response/context panel; not a measured Xenium flux or functional assay'})
  if not independent:continue
  ix=pd.Index(genes).get_indexer(independent);response_scores[program]=log['corrected'][:,ix].mean(1)
 for stage in STAGES:
  for lib in sorted(obs.loc[obs.age.eq(stage),'sample'].unique()):
   for state in ['All',*sorted(obs.loc[obs.receiver_ventricular,'CM_state'].dropna().unique())]:
    ix=obs.age.eq(stage)&obs['sample'].eq(lib)&obs.receiver_ventricular
    if state!='All':ix &= obs.CM_state.eq(state)
    if ix.sum()<20:continue
    for program,values in response_scores.items():
     score_rows.append({'stage':stage,'sample':lib,'CM_state':state,'program':program,'nuclei':int(ix.sum()),'mean_log1p_CP10k':values[ix].mean()})
 save(pd.DataFrame(definitions),'RNA_independent_response_programs.tsv');save(pd.DataFrame(score_rows),'RNA_CM_programs_by_library.tsv')
 methods=config();methods.update({'RNA_reference_nuclei':len(obs),'RNA_selected_gene_count':len(genes),
    'RNA_ventricular_receiver_nuclei':int(obs.receiver_ventricular.sum()),'complete_measured_pairs':len(candidates),
    'RNA_supported_pairs_before_pathway_reduction':len(pool),'selected_pairs':len(selected),
    'RNA_full_genome_denominators':True,'RNA_duplicate_symbols_summed':True})
 (OUT/'configuration.json').write_text(json.dumps(methods,indent=2)+'\n');export_inputs()
 print('RNA_SIGNALING_SUPPORT_COMPLETE',json.dumps({k:methods[k] for k in ['RNA_reference_nuclei','RNA_ventricular_receiver_nuclei','complete_measured_pairs','RNA_supported_pairs_before_pathway_reduction','selected_pairs']}),round(time.monotonic()-start,1),flush=True)

if __name__=='__main__':main()
