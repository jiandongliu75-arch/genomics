"""Insert the complete held-out results into the biological report template."""
import json
from pathlib import Path
import pandas as pd
from joint_signaling_common import ROOT,OUT,STAGES

def main():
 g=pd.read_csv(OUT/'incremental_prediction_gains.tsv',sep='\t')
 primary=g[g.primary_analysis];hits=primary[primary.primary_grid_BH_q.lt(.05)&primary.delta_R2.gt(0)]
 selected=g[g.view.eq('core')&g.learner.eq('RandomForest')&g.CM_state.eq('All')&g.outcome.eq('core_Z__maturation21')&g.comparison.eq('E minus C_fine_context')].set_index('stage').reindex(STAGES)
 values='; '.join(f'{s}, {selected.loc[s,"delta_R2"]:+.4f}' for s in STAGES)
 if (selected.delta_R2<=0).all():
  opening='At all four ages, adding source-specific non-CM expression and collective ligand–receptor inputs did not improve prediction of the overall 21-gene score beyond the receptor, CM-to-CM and fine-context controls.'
 elif selected.delta_R2.abs().max()<.01:
  opening='Across ages, source-specific non-CM molecular inputs added little to prediction of the overall 21-gene score beyond the receptor, CM-to-CM and fine-context controls.'
 else:
  opening='The contribution of source-specific non-CM molecular inputs to prediction of the overall 21-gene score differed across ages.'
 results=opening+f' The changes in held-out R² were {values} (Fig. 3). These numbers describe prediction of observed RNA in unseen regions of the same section, rather than a biological effect of adding or removing a cell population.'
 if len(hits)==0:
  results+='\n\nNone of the 40 supported primary comparisons showed a positive prediction gain after within-section multiple-testing correction. The individual developmental components and source sensitivities are supplied in the figures and tables; small or negative increments should not be turned into claims that non-CMs are inactive or suppress CM development.'
 else:
  labels={'metabolic':'metabolic-regulatory','sarcomere':'sarcomere','excitation_contraction':'excitation–contraction','cycling':'cell-cycle','maturation21':'overall 21-gene'}
  entries=[]
  for r in hits.itertuples():
   label=labels[r.outcome.removeprefix('core_Z__')]
   entries.append(f'{r.stage} {label} ({r.comparison.replace("C_fine_context","fine context")}, ΔR² {r.delta_R2:+.4f}, q={r.primary_grid_BH_q:.3g})')
  results+=f'\n\nOf the 40 supported primary comparisons, {len(hits)} showed a positive gain after within-section multiple-testing correction: '+ '; '.join(entries)+'. These remain conditional prediction results from one section per age, requiring independent-heart validation.'
  for r in hits.itertuples():
   key=r.outcome.removeprefix('core_Z__')
   sensitivity=g[g.stage.eq(r.stage)&g.CM_state.eq('All')&g.comparison.eq(r.comparison)]
   native=sensitivity[sensitivity.view.eq('native')&sensitivity.learner.eq('RandomForest')&sensitivity.outcome.eq('native_Z__'+key)].iloc[0]
   linear=sensitivity[sensitivity.view.eq('core')&sensitivity.learner.eq('Ridge')&sensitivity.outcome.eq(r.outcome)].iloc[0]
   results+=f'\n\nFor the {r.stage} {labels[key]} comparison, the gain amounted to about {100*r.delta_R2:.2f} percentage points of explained variation. The deposited-RNA forest increment was {native.delta_R2:+.4f}.'
   if native.delta_R2<=0 or linear.delta_R2<=0:
    results+=' The linear sensitivity also did not support the gain.' if linear.delta_R2<=0 else ''
    results+=' The positive result therefore depends on the RNA view or learner and does not establish a reproducible source-specific mechanism.'
 results+='\n\nThe baseline already includes which cells surround a CM and where it is located. The added features ask whether measured molecular differences explain more of the CM outcome. Vascular support, matrix organization and regional development can therefore remain biologically important even when the molecular increment is small. CM state and receptor expression may also lie along a biological response pathway; conditioning on them does not estimate the total influence of neighboring cells.'
 priorities=('The next priority is to connect vascular access to specific CM functions. Compare CMs of the same identity and wall region at P7 and P14, measuring perfused capillaries together with CM metabolism, sarcomere organization and calcium handling. Add independent hearts and intermediate stages, especially P3/P5 and approximately P10. This would test whether functional vascular integration explains developmental variation better than vessel proximity alone.\n\n'
 'For fibroblasts, measure matrix protein organization around the same CMs, including early remodeling-associated matrix and later laminin-associated structures. Compare regions with similar vascular access. If a matrix difference accompanies CM structural organization while metabolic measures remain similar, prioritize a mechanical or matrix-organization experiment. If CM metabolism instead tracks perfusion, prioritize the vascular-function hypothesis. The current RNA results do not establish a single fibroblast family as the driver of the overall maturation index.\n\n'
 'Test a particular endothelial contact or secreted signal after verifying the ligand and receptor proteins in the relevant cells. For the neuregulin candidate, include the receptor partner missing from the Xenium panel. Compare direct EC–CM contact with physically separated cultures at matched oxygen and nutrient supply, and assess contraction, calcium handling and metabolism separately from cell-cycle RNA. This would distinguish a contact-dependent response from diffusible support or altered oxygen delivery.')
 source=(ROOT/'analysis/joint_signaling_report_template.md').read_text()
 assert source.count('[[RESULTS]]')==source.count('[[PRIORITIES]]')==1
 text=source.replace('[[RESULTS]]',results).replace('[[PRIORITIES]]',priorities)
 (ROOT/'analysis/joint_signaling_report.md').write_text(text)
 (OUT/'report_quantitative_claims.json').write_text(json.dumps({'overall_21gene_E_minus_fine_context_core':selected.delta_R2.to_dict(),
  'positive_primary_corrected_prediction_gains':len(hits),'primary_tests':int(primary.primary_grid_BH_q.notna().sum()),
  'interpretation':'Fixed-model conditional prediction; not causal source effects'},indent=2)+'\n')
 print('JOINT_REPORT_RESULTS_COMPOSED',len(hits),'positive primary corrected gains',flush=True)
if __name__=='__main__':main()
