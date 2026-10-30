"""Remove CM maturation/lineage genes from FB rank and control backgrounds."""
from fibroblast_CM_refit_checks import *
from fibroblast_multimethod_identity import EXCLUDE

def main():
 ann,x,raw,genes=refs();exclude=set(EXCLUDE)|{'Slit2','Ddc'};use=~np.isin(genes,list(exclude));g=genes[use]
 assert set(sum(PANELS.values(),[])).issubset(g)
 ref=x[:,use];mean=(ref/np.maximum(ref.sum(1),1)[:,None]).mean(0);rows=[];sources=[]
 for stage in STAGES:
  fb,q=queries(stage,genes);f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  d=loaded(stage,'proximity50');d=d[eligible(d)].copy()
  path=ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz';graph=sparse.load_npz(path).tocsr()[f.index.get_indexer(d.index)][:,f.index.get_indexer(fb.index)]
  for view in ['refined','core']:
   xx=q[view][:,use];depth=xx.sum(1);det=(xx>0).sum(1);good=fb.program_source_eligible.to_numpy(bool)&(depth>=20)&(det>=8)
   mod=score_modules(xx,g,PANELS,mean);layer='noCM_'+view
   denom=np.asarray(graph@good.astype(float)).ravel()
   d[layer+'__information_fraction']=denom/np.maximum(d.eligible_FB_neighbors.to_numpy(),1)
   d[layer+'__log_depth']=np.asarray(graph@(good*np.log1p(depth))).ravel()/np.maximum(denom,1)
   d[layer+'__detected_genes']=np.asarray(graph@(good*det)).ravel()/np.maximum(denom,1)
   vals=np.column_stack([mod[m+'__rank_excess'] for m in MODULES]);out=graph@(vals*good[:,None])/np.maximum(denom[:,None],1)
   for j,m in enumerate(MODULES):
    d[layer+'__'+m]=out[:,j]
    sources.append({'stage':stage,'view':view,'program':m,'retained_scoring_genes':len(g),'excluded_genes':';'.join(sorted(set(genes)&exclude)),
      'informative_profiles':int(good.sum()),'informative_with_two_markers':int((good&(mod[m+'__measured_markers'].to_numpy()>=2)).sum())})
    if sources[-1]['informative_with_two_markers']<30:continue
    a,_,_=fit_target(d.loc[denom>=1],layer+'__'+m,'program',layer,'within_region','proximity50');rows.extend(a)
  print('REFIT_CM_EXCLUDED',stage,flush=True)
 save(pd.DataFrame(rows),'CM_gene_excluded_fibroblast_program_sensitivity.tsv');save(pd.DataFrame(sources),'CM_gene_exclusion_definition.tsv');export_inputs()
if __name__=='__main__':main()
