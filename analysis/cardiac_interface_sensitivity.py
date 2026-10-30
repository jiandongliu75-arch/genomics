"""Sample-resolution and native/refined geometry checks without new labels."""
from reconstruct_cardiac_spatial import ROOT,OUT,UP,CFG,STAGES,read,save,polys
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr
from shapely import STRtree,line_interpolate_point,get_x
import pandas as pd

def run():
    rows=[];samples=[];rng=np.random.default_rng(CFG['seed'])
    for stage in STAGES:
        print('INTERFACE_SENSITIVITY',stage,flush=True)
        q=read(OUT/stage/'all_profile_observations.tsv.gz');c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');e=read(OUT/stage/'vascular_interface_and_matrix_context.tsv.gz')
        ci=q.index.get_indexer(c.index);source=q.author_atlas_cell_state.eq('cap.EC').to_numpy(float)
        for label,path in [('refined_boundary',UP/stage/'contact_graph.npz'),('native_polygon_1um',ROOT/f'results/cardiac_niche_dynamics/{stage}/native_boundary_1um_graph.npz')]:
            g=sparse.load_npz(path).tocsr()[ci];e[label+'__capillary_contacts']=np.asarray(g@source).ravel()
            take=c.reconstructed_Figure4e_CM
            rows.append({'stage':stage,'comparison':label+'_contacts_vs_native_sampled_apposition','CMs':int(take.sum()),
                'Spearman':spearmanr(e.loc[take,label+'__capillary_contacts'],e.loc[take,'author_capillary__apposition_fraction_2um']).statistic,
                'interpretation':'Different operational2D contact definitions; neither is protein-contact ground truth'})
        pp,_=polys(stage,q.index);source_indices=np.flatnonzero(source.astype(bool));tree=STRtree(np.array([pp[i] for i in source_indices],dtype=object))
        ids=c.index[c.reconstructed_Figure4e_CM];chosen=rng.choice(ids,min(256,len(ids)),replace=False)
        boundaries=np.array([pp[q.index.get_loc(cid)].boundary for cid in chosen],dtype=object)
        for ns in [128,256]:
            pts=line_interpolate_point(boundaries[:,None],((np.arange(ns)+.5)/ns)[None,:],normalized=True)
            pairs,dd=tree.query_nearest(pts.ravel(),return_distance=True,all_matches=False)
            ds=np.full(pts.size,np.nan);ds[pairs[0]]=dd;ds=ds.reshape((len(chosen),ns));fraction=np.mean(ds<=2,axis=1)
            old=e.loc[chosen,'author_capillary__apposition_fraction_2um'].to_numpy()
            rows.append({'stage':stage,'comparison':f'64_vs_{ns}_boundary_samples','CMs':len(chosen),'Spearman':spearmanr(old,fraction).statistic,
                'mean_absolute_apposition_fraction_difference':abs(old-fraction).mean(),'interpretation':'Outcome-blinded resolution sensitivity on original native polygons'})
            for cid,value in zip(chosen,fraction):samples.append({'stage':stage,'cell_id':cid,'boundary_samples':ns,'fraction_within2um_of_author_capillary':value,
                'original64_fraction':e.loc[cid,'author_capillary__apposition_fraction_2um']})
        save(e[[k for k in e if k.startswith(('refined_boundary__','native_polygon_1um__'))]],f'{stage}/alternative_contact_counts.tsv.gz',index=True)
    save(pd.DataFrame(rows),'interface_geometry_sensitivity.tsv');save(pd.DataFrame(samples),'interface_resolution_sample.tsv')

if __name__=='__main__':run()
