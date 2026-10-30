"""Joint-program specificity and same-CM vascular geometry checks."""
from fibroblast_CM_refit_checks import *

def extensions():
 coverage=read(OUT/'fibroblast_program_coverage.tsv');rows=[];vrows=[];diags=[]
 for stage in STAGES:
  a=loaded(stage,'proximity50');d=a[eligible(a)].copy()
  for view in ['refined','core']:
   mods=coverage.loc[coverage.stage.eq(stage)&coverage.view.eq(view)&coverage.informative_with_two_markers.ge(30),'program'].tolist()
   for mod in mods:
    target=view+'__'+mod;others=[view+'__'+m for m in mods if m!=mod]
    r,_,di=fit_target(d,target,'program',view,'within_region','proximity50',extra_cov=others)
    rows.extend([{'other_FB_programs_adjusted':True,**x} for x in r]);diags.extend(di)
  for other in ['proximity30','contact']:
   b=loaded(stage,other);ids=a.index[eligible(a)&eligible(b)]
   for geom,z in [('proximity50',a.loc[ids]),(other,b.loc[ids])]:
    cov=[k for k in ['candidate__'+f for f in FAMILIES] if z[k].std()>1e-6]
    for target in ['EC__Capillary-like','EC__Endocardial-like']:
     r,_,di=fit_target(z,target,'vascular','candidate','within_region',geom,extra_cov=cov)
     vrows.extend([{'paired_geometry':other,'same_CM_profiles':True,**x} for x in r]);diags.extend(di)
  print('REFIT_SPECIFICITY',stage,flush=True)
 save(pd.DataFrame(rows),'joint_fibroblast_program_sensitivity.tsv');save(pd.DataFrame(vrows),'same_CM_vascular_geometry_sensitivity.tsv');save(pd.DataFrame(diags),'extension_model_designs.tsv');export_inputs()
if __name__=='__main__':extensions()
