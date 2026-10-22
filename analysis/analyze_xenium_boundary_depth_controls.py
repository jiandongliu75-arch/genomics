"""Count-matched random thinning controls for boundary erosion/QV sensitivity."""
import argparse
import json

import numpy as np
import pandas as pd

from design_xenium_boundary_qc import OUT,ATLAS,STAGES
from boundary_reference_projection import FrozenProjection
from analyze_xenium_boundary_sensitivity import pair_effects

CONTROLS={'thin_to_core1':('core_1um',501),'thin_to_core2':('core_2um',502),'thin_to_Q30':('native_Q30',530)}


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=STAGES,required=True);args=p.parse_args();stage=args.stage;folder=OUT/stage
    if (folder/'depth_controls_completion.json').exists():print('SKIP_DEPTH_CONTROLS',stage,flush=True);return
    design={'reason':'P0 erosion removes substantial RNA; control identity changes caused by count depth alone before interpreting boundary-specific effects',
        'scope':'Added technical controls, not outcome-optimized segmentation settings or biological replication',
        'sampling':'Within each original cell, multivariate hypergeometric sampling without replacement to exactly the perturbed gene-transcript total',
        'controls':{k:{'matches':v[0],'seed':v[1]} for k,v in CONTROLS.items()},
        'nucleus_required':False,'one_fixed_technical_draw_per_condition':True}
    path=OUT/'boundary_depth_control_design.json'
    if path.exists():assert json.loads(path.read_text())==design
    else:path.write_text(json.dumps(design,indent=2)+'\n')
    with np.load(folder/'geometry_sensitivity_counts.npz') as b:
        ids,genes=b['cell_ids'].astype(str),b['genes'].astype(str)
        native=b['native_Q20'];targets={k:b[spec[0]].sum(1) for k,spec in CONTROLS.items()}
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    model=FrozenProjection(stage,obs,genes);matrices={};effects=[];retention=[]
    for variant,(_,seed) in CONTROLS.items():
        rng=np.random.default_rng(seed);counts=np.zeros_like(native)
        for i in range(len(native)):
            n=int(targets[variant][i]);assert n<=native[i].sum()
            if n:counts[i]=rng.multivariate_hypergeometric(native[i],n)
        assert (counts<=native).all();np.testing.assert_array_equal(counts.sum(1),targets[variant])
        matrices[variant]=counts;classified=model.classify(ids,counts)
        classified.to_csv(folder/(variant+'_projection_assignments.tsv.gz'),sep='\t',na_rep='NA')
        effects.extend(pair_effects(stage,variant,'geometry_targets_and_all_frozen_pairs',ids,counts,genes,obs,classified))
        base=obs.loc[ids,'atlas_identity'];missing=obs.loc[ids,'nucleus_area'].isna().to_numpy()
        for label in sorted(base.unique()):
            for group,keep in [('all',np.ones(len(ids),bool)),('measured',~missing),('missing',missing)]:
                mask=base.eq(label).to_numpy()&keep;n=int(mask.sum())
                if n:
                    ratio=counts[mask].sum(1)/np.maximum(native[mask].sum(1),1)
                    retention.append({'stage':stage,'variant':variant,'scope':'geometry_targets_and_all_frozen_pairs',
                        'baseline_identity':label,'nuclear_area':group,'profiles':n,
                        'information_QC_pass':int(classified.information_QC_pass.to_numpy()[mask].sum()),
                        'identity_unchanged':int(np.equal(classified.new_identity.to_numpy()[mask],base.to_numpy()[mask]).sum()),
                        'mean_transcript_retained_ratio':float(ratio.mean()),'median_transcript_retained_ratio':float(np.median(ratio))})
        print('DEPTH_CONTROL',stage,variant,len(ids),flush=True)
    np.savez_compressed(folder/'depth_matched_control_counts.npz',cell_ids=ids,genes=genes,**matrices)
    pd.DataFrame(effects).to_csv(folder/'Ankrd1_depth_control_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(retention).to_csv(folder/'depth_control_identity_retention.tsv',sep='\t',index=False,na_rep='NA')
    (folder/'depth_controls_completion.json').write_text(json.dumps({'stage':stage,'profiles':len(ids),
        'exact_count_depth_matching':True,'nucleus_required':False,'conditions':list(CONTROLS)},indent=2)+'\n')
    print('DEPTH_CONTROLS_COMPLETE',stage,flush=True)


if __name__=='__main__':main()
