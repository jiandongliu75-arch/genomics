"""Paper-definition reconstruction and independent geometry/outcome preparation."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-reconstruction-mpl')
from pathlib import Path
import json,argparse,gc,hashlib,itertools
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse,ndimage
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
from skimage.morphology import convex_hull_image,disk
from PIL import Image,ImageDraw
from shapely import Polygon,STRtree,line_interpolate_point,get_x,get_y

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'analysis/cardiac_spatial_reconstruction.json').read_text())
OUT=ROOT/CFG['output'];UP=ROOT/'results/cardiac_spatial_pipeline';STAGES=CFG['stages']
MAT=pd.read_csv(ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',sep='\t')
PROGS=json.loads((ROOT/'analysis/cardiac_spatial_pipeline.json').read_text())['programs']
def read(p,**kw):kw.setdefault('index_col',0);return pd.read_csv(p,sep='\t',**kw)
def save(d,p,index=False):
    p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,sep='\t',index=index,na_rep='NA',float_format='%.9g')
def dump(p,d):
    p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,default=str)+'\n')
def logmatrix(x):return x.multiply((1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1))[:,None]).log1p().tocsr()
def native_path(stage):return ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad'
def polys(stage,ids):
    p=next(Path('/workspace/heart-atlas/data/xenium').glob('WT_'+('P07' if stage=='P7' else stage)+'__*'))/'cell_boundaries.csv.gz'
    v=pd.read_csv(p,usecols=['cell_id','vertex_x','vertex_y']);lookup={x:i for i,x in enumerate(ids)};pp=[None]*len(ids)
    for cid,z in v.groupby('cell_id',sort=False):
        if cid in lookup and len(z)>=3:pp[lookup[cid]]=Polygon(z[['vertex_x','vertex_y']].to_numpy())
    return pp,p
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def source_audit():
    source=Path('/workspace/heart-atlas/reference/PIP_seq')
    fields=json.loads(Path('/workspace/heart-atlas/data/author_webatlas/P7/obs/.zattrs').read_text())['column-order']
    rows=[('Figure4d CM population','Exact author v.CM label; all four count totals recoverable','Recovered'),
        ('Figure4e CM population','Released six-state CM list at P7 plus four named anatomical source categories; different from Figure4d','Selection rules recovered; exact membership unavailable'),
        ('Figure4e anatomical source','ventricle_left_compact;trabecular_right;ventricle_right_compact;trabecular_left','Names recovered; source column and original masks absent from public export'),
        ('Maturation_Features1','Precomputed score read from unpublished RDS; main text standardized mean versus Methods AddModuleScore','Exact score construction/assay/pool not recoverable'),
        ('Neighbor graph / levels','Precomputed xenium_p7_metadata_cellneighbour_v2.csv; graph-builder script not released','Exact graph/levels unavailable; explicit generic6NN and tissue-screened alternatives reconstructed'),
        ('Figure4e trend plot','Mean score within each neighbor level, followed by LOESS on the means','Plotting rule recovered; input score/mask/levels not supplied'),
        ('Intermediate days / replication','Observed RNA libraries and spatial sections only at P0/P7/P14/P21','Missing ages can be interpolated only; cannot add observed replication')]
    save(pd.DataFrame(rows,columns=['item','released_evidence','reconstruction_status']),'author_reconstruction_audit.tsv')
    dump('author_definition.json',{'source_commit':'102d80e6ea5de44f594461c727a4b3c6e1bc3508','public_obs_fields':fields,
        'Figure4d_population':'v.CM','Figure4e_populations':CFG['CM_author_states'],
        'Figure4e_source_names':['ventricle_left_compact','trabecular_right','ventricle_right_compact','trabecular_left'],
        'original_mask_available':False,'original_score_available':False,'original_graph_available':False,
        'exact_Figure4e_replication_claimed':False,'independent_morphology_masks_use_CM_maturation':False,
        'script_sha256':sha(source/'Figure4/4.2_main_figure_generation.R')})

def anatomy(stage,f,pp):
    res=CFG['anatomy_grid_um'];info=json.loads((ROOT/f'results/xenium_ventricular_mapping/{stage}/morphology_display.json').read_text())
    w,h=info['extent_um'][1],info['extent_um'][2];shape=(int(np.ceil(h/res))+1,int(np.ceil(w/res))+1)
    canvas=Image.new('1',(shape[1],shape[0]));draw=ImageDraw.Draw(canvas)
    for poly in pp:
        if poly is not None and not poly.is_empty:draw.polygon([(x/res,y/res) for x,y in poly.exterior.coords],fill=1)
    tissue=np.array(canvas,bool);tissue=ndimage.binary_closing(tissue,structure=disk(2))
    select=Image.new('1',(shape[1],shape[0]));draw=ImageDraw.Draw(select)
    poly=np.array(CFG['ventricular_selector_normalized_xy'][stage])*np.array([w/res,h/res]);draw.polygon([tuple(x) for x in poly],fill=1)
    selected=np.array(select,bool);vent_tissue=tissue&selected
    lab,n=ndimage.label(vent_tissue);sizes=np.bincount(lab.ravel());sizes[0]=0
    main=lab==sizes.argmax();envelope=convex_hull_image(main)&selected
    void=envelope&~tissue;lab,n=ndimage.label(void);sizes=np.bincount(lab.ravel())*res**2
    big=np.flatnonzero(sizes>=CFG['void_minimum_area_um2']);big=big[big>0];large_void=np.isin(lab,big)
    distance_to_void=ndimage.distance_transform_edt(~large_void)*res if large_void.any() else np.full(shape,np.nan)
    distance_to_outer=ndimage.distance_transform_edt(envelope)*res
    porosity=ndimage.uniform_filter((~tissue).astype(float),size=2*int(CFG['porosity_radius_um']/res)+1)
    xy=f[['x_centroid','y_centroid']].to_numpy();ij=np.rint(xy[:,::-1]/res).astype(int);ij=np.minimum(np.maximum(ij,0),np.array(shape)-1)
    q=f[['stage_cell_id','x_centroid','y_centroid','author_atlas_cell_state']].copy()
    q['ventricular_morphology_selector']=envelope[tuple(ij.T)]&tissue[tuple(ij.T)]
    q['distance_to_large_void_um']=distance_to_void[tuple(ij.T)];q['distance_to_outer_um']=distance_to_outer[tuple(ij.T)]
    q['local_void_fraction']=porosity[tuple(ij.T)]
    q['relative_wall_geometry']=q.distance_to_large_void_um/(q.distance_to_large_void_um+q.distance_to_outer_um).clip(lower=1)
    # Thresholds use physical morphology only. No CM score or EC density enters.
    q['morphology_region']='Transition/unclassified'
    dense=q.local_void_fraction.le(.04)&q.distance_to_large_void_um.ge(80)&q.distance_to_outer_um.ge(30)
    porous=q.local_void_fraction.ge(.08)&q.distance_to_large_void_um.le(100)&q.distance_to_outer_um.ge(30)
    q.loc[dense,'morphology_region']='Compact-compatible dense myocardium'
    q.loc[porous,'morphology_region']='Cavity-associated porous myocardium'
    q.loc[~q.ventricular_morphology_selector,'morphology_region']='Outside conservative ventricular selector'
    q['anatomy_block_200um']=(np.floor(q.x_centroid/CFG['anatomy_block_um']).astype(int).astype(str)+'_'+np.floor(q.y_centroid/CFG['anatomy_block_um']).astype(int).astype(str))
    q['author_Figure4d_CM']=q.author_atlas_cell_state.eq('v.CM')
    q['author_Figure4e_label_eligible']=q.author_atlas_cell_state.isin(CFG['CM_author_states'][stage])
    q['reconstructed_Figure4e_CM']=q.author_Figure4e_label_eligible&q.ventricular_morphology_selector
    q['broader_ventricular_CM']=f.ventricular_CM_strict&q.ventricular_morphology_selector
    save(q,f'{stage}/anatomy_and_cohorts.tsv.gz',index=True)
    np.savez_compressed(OUT/stage/'morphology_masks.npz',tissue=tissue,envelope=envelope,large_void=large_void,resolution_um=res,
        distance_to_void_um=distance_to_void.astype(np.float32),distance_to_outer_um=distance_to_outer.astype(np.float32),porosity=porosity.astype(np.float32))
    return q

def prepare():
    sums={};moments={};nall=0;ncm=0;sources=[];cohorts=[];extras='Ptprm Nampt Insr Flt4 Dll4 Hey1 Hey2 Notch1 Npr3 Nrxn3'.split()
    # Count-space source files are read-only and never rewritten.
    for stage in STAGES:
        print('RECONSTRUCT_PREPARE',stage,flush=True);(OUT/stage).mkdir(parents=True,exist_ok=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz');pp,p=polys(stage,f.index);q=anatomy(stage,f,pp);del pp;gc.collect()
        a=ad.read_h5ad(native_path(stage),backed='r');assert a.obs_names.equals(f.index);genes=a.var_names.to_numpy(dtype=str)
        b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        q['native_panel_depth']=np.asarray(a.X[:].sum(1)).ravel();q['native_detected_genes']=np.asarray((a.X[:]>0).sum(1)).ravel()
        for view,x in [('native',a.X[:]),('core',b.layers['high_confidence'][:]),('nuclear',b.layers['nuclear_only'][:])]:
            x=x.tocsr();l=logmatrix(x);cm=q.author_Figure4e_label_eligible.to_numpy()
            if view in ['native','core']:
                s=np.asarray(l.sum(0)).ravel();s2=np.asarray(l.power(2).sum(0)).ravel()
                sc=np.asarray(l[cm].sum(0)).ravel();sc2=np.asarray(l[cm].power(2).sum(0)).ravel()
                moments[(stage,view)]={'all_n':len(f),'all_sum':s,'all_sum2':s2,'CM_n':cm.sum(),'CM_sum':sc,'CM_sum2':sc2}
                for key,panel in {'maturation21':MAT.gene.tolist(),**CFG['CM_components']}.items():
                    gi=np.isin(genes,panel);q[view+'_logmean__'+key]=np.asarray(l[:,gi].mean(1)).ravel()
            q[view+'_panel_depth']=np.asarray(x.sum(1)).ravel()
            for gene in extras:
                if gene in genes:
                    gi=list(genes).index(gene);q[view+'__'+gene]=l[:,gi].toarray().ravel();q[view+'_molecules__'+gene]=x[:,gi].toarray().ravel()
            if view in ['native','core']:
                for key,program in PROGS.items():
                    q[view+'__'+key]=np.asarray(l[:,np.isin(genes,program['genes'])].mean(1)).ravel()
        q['cell_area']=f.cell_area;q['nucleus_area_missing']=f.nucleus_area.isna();q['atlas_identity']=f.atlas_identity;q['program_source_eligible']=f.program_source_eligible
        q['our_ventricular_CM_strict']=f.ventricular_CM_strict;q['present_in_public_atlas']=~f.author_atlas_cell_state.eq('Not present in published atlas')
        save(q,f'{stage}/all_profile_observations.tsv.gz',index=True)
        sources.extend([native_path(stage),ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',UP/stage/'atlas_annotations.tsv.gz',p])
        cohorts.append({'stage':stage,'published_vCM_Figure4d':int(q.author_Figure4d_CM.sum()),'paper_Figure4e_label_set_all_tissue':int(q.author_Figure4e_label_eligible.sum()),
            'independent_ventricular_mask_paper_labels':int(q.reconstructed_Figure4e_CM.sum()),'our_ventricular_within_mask':int(q.broader_ventricular_CM.sum()),
            'author_source_mask_recovered':False})
        a.file.close();b.file.close();del a,b,l,x,q;gc.collect()
    cal=[]
    for view in ['native','core']:
        for scope in ['all','CM']:
            nn=sum(moments[(s,view)][scope+'_n'] for s in STAGES)
            mean=sum(moments[(s,view)][scope+'_sum'] for s in STAGES)/nn
            sd=np.sqrt(np.maximum(sum(moments[(s,view)][scope+'_sum2'] for s in STAGES)/nn-mean**2,1e-10))
            for j,gene in enumerate(genes):cal.append({'view':view,'scope':scope,'gene':gene,'reference_mean':mean[j],'reference_SD':sd[j],'reference_profiles':int(nn)})
    cal=pd.DataFrame(cal);save(cal,'score_reference_calibration.tsv');save(pd.DataFrame(cohorts),'cohort_reconstruction_counts.tsv')
    # Explicitly feasible matched-control sensitivity. Unlike the unknown author
    # assay/binning, this uses403 measured genes,3 bins and100 controls/target.
    nativecal=cal[cal.view.eq('native')&cal.scope.eq('all')].set_index('gene').loc[genes]
    bins=np.empty(len(genes),int)
    for i,part in enumerate(np.array_split(np.argsort(nativecal.reference_mean.to_numpy(),kind='stable'),3)):bins[part]=i
    rng=np.random.default_rng(CFG['seed']);mapping=[]
    for gene in MAT.gene:
        gi=list(genes).index(gene);pool=np.flatnonzero(bins==bins[gi]);chosen=rng.choice(pool,100,replace=False)
        mapping.extend({'target_gene':gene,'control_gene':genes[j],'bin':int(bins[j])} for j in chosen)
    controls=pd.DataFrame(mapping);save(controls,'matched_control_3bin_100genes.tsv')
    for stage in STAGES:
        print('RECONSTRUCT_SCORES',stage,flush=True)
        q=read(OUT/stage/'all_profile_observations.tsv.gz');take=q.author_Figure4e_label_eligible|q.our_ventricular_CM_strict
        c=q.loc[take].copy();f=read(UP/stage/'atlas_annotations.tsv.gz')
        a=ad.read_h5ad(native_path(stage),backed='r');b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r');ix=a.obs_names.get_indexer(c.index)
        for view,x in [('native',a.X[ix]),('core',b.layers['high_confidence'][ix])]:
            l=logmatrix(x);dense=l.toarray();v=cal[cal.view.eq(view)&cal.scope.eq('all')].set_index('gene').loc[genes]
            z=(dense-v.reference_mean.to_numpy())/v.reference_SD.to_numpy()
            for key,panel in {'maturation21':MAT.gene.tolist(),**CFG['CM_components']}.items():c[view+'_Z__'+key]=z[:,np.isin(genes,panel)].mean(1)
            c[view+'_Z__withdrawal']=-c[view+'_Z__cycling']
            vv=cal[cal.view.eq(view)&cal.scope.eq('CM')].set_index('gene').loc[genes]
            c[view+'_CMreference_Z__maturation21']=((dense-vv.reference_mean.to_numpy())/vv.reference_SD.to_numpy())[:,np.isin(genes,MAT.gene)].mean(1)
            m=moments[(stage,view)];mean=m['all_sum']/m['all_n'];sd=np.sqrt(np.maximum(m['all_sum2']/m['all_n']-mean**2,1e-10))
            c[view+'_within_stage_Z__maturation21']=((dense-mean)/sd)[:,np.isin(genes,MAT.gene)].mean(1)
            c[view+'_matched3_100__maturation21']=dense[:,np.isin(genes,MAT.gene)].mean(1)-dense[:,np.isin(genes,controls.control_gene)].mean(1)
        save(c,f'{stage}/CM_outcomes_and_anatomy.tsv.gz',index=True);a.file.close();b.file.close();del q,c,a,b,l,dense,z;gc.collect()
    dump('source_files.json',[str(p) for p in sorted(set(sources))]);dump('configuration.json',CFG)
    definition=[]
    for key,panel in {'maturation21':MAT.gene.tolist(),**CFG['CM_components']}.items():definition.append({'program':key,'genes':';'.join(panel),'measurement':'Direct measured403-gene Xenium RNA; native and reassigned core views','interpretation':'Transcript-program proxy, not functional measurement'})
    save(pd.DataFrame(definition),'outcome_program_definitions.tsv')

def geometry():
    audits=[];ringrows=[]
    for stage in STAGES:
        print('RECONSTRUCT_GEOMETRY',stage,flush=True)
        q=read(OUT/stage/'all_profile_observations.tsv.gz');c=read(OUT/stage/'CM_outcomes_and_anatomy.tsv.gz');xy=q[['x_centroid','y_centroid']].to_numpy();ci=q.index.get_indexer(c.index)
        valid=np.flatnonzero(q.present_in_public_atlas.to_numpy());dist,neigh=cKDTree(xy[valid]).query(xy[valid],k=CFG['paper_like_knn']+1)
        ii=np.repeat(valid,CFG['paper_like_knn']);jj=valid[neigh[:,1:].ravel()]
        graph=sparse.csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(q),len(q)));graph=graph.maximum(graph.T).tocsr()
        sparse.save_npz(OUT/stage/'paper_like_symmetric6NN_graph.npz',graph)
        mask=sparse.load_npz(ROOT/f'results/cardiac_niche_dynamics/{stage}/pruned_Delaunay_graph.npz').tocsr()
        published=sparse.diags(q.present_in_public_atlas.to_numpy(float));mask=(published@mask@published).tocsr()
        sparse.save_npz(OUT/stage/'tissue_screened_Delaunay_graph.npz',mask)
        ring=c[['author_atlas_cell_state','reconstructed_Figure4e_CM','author_Figure4d_CM','morphology_region']].copy()
        for name,source in [('cap.EC',q.author_atlas_cell_state.eq('cap.EC')),('cap.EC_plus_HighMyo',q.author_atlas_cell_state.isin(['cap.EC','cap.EC.High.Myo.'])),('p.EC',q.author_atlas_cell_state.eq('p.EC')),('Peri',q.author_atlas_cell_state.eq('Peri'))]:
            si=np.flatnonzero(source.to_numpy())
            ring['nearest_centroid__'+name]=cKDTree(xy[si]).query(xy[ci])[0]
            for label,g in [('paper_like6NN',graph),('tissue_Delaunay',mask)]:
                d=dijkstra(g,indices=si,directed=False,unweighted=True,min_only=True)[ci]
                ring[label+'__'+name]=np.where(np.isfinite(d),np.minimum(d,10),np.nan)
                ring[label+'__'+name+'__censored']=~np.isfinite(d)|(d>=10)
                audits.append({'stage':stage,'graph':label,'source':name,'source_profiles':len(si),'CMs':len(c),'CMs_unreachable':int((~np.isfinite(d)).sum()),'CMs_ge10':int((d>=10).sum())})
        save(ring,f'{stage}/reconstructed_neighbor_levels.tsv.gz',index=True)
        # Original polygon surface sampling makes geometry different from RNA
        # abundance. It remains operational2D apposition, not protein contact.
        pp,path=polys(stage,q.index);allpolys=np.array(pp,dtype=object);validpoly=np.array([p is not None and not p.is_empty for p in pp])
        frame=c[['x_centroid','y_centroid']].copy();ns=CFG['boundary_samples_per_CM'];fractions=(np.arange(ns)+.5)/ns
        cmvalid=validpoly[ci];boundaries=np.array([pp[i].boundary if validpoly[i] else None for i in ci],dtype=object)
        pts=line_interpolate_point(boundaries[:,None],fractions[None,:],normalized=True)
        flat=pts.ravel();pointvalid=~np.isnan(get_x(flat));flatvalid=flat[pointvalid]
        perimeters=np.array([pp[i].length if validpoly[i] else np.nan for i in ci]);frame['CM_native_perimeter_um']=perimeters
        groups={'author_capillary':q.author_atlas_cell_state.eq('cap.EC'),
            'author_capillary_plus_HighMyo':q.author_atlas_cell_state.isin(['cap.EC','cap.EC.High.Myo.']),
            'QC_capillary':q.atlas_identity.eq('Endothelial')&q.program_source_eligible&read(UP/stage/'atlas_annotations.tsv.gz').EC_nuclear_type.eq('Capillary-like'),
            'author_pericyte':q.author_atlas_cell_state.eq('Peri')}
        radius=sparse.load_npz(UP/stage/'proximity50_graph.npz').tocsr()[ci]
        for label,source in groups.items():
            si=np.flatnonzero(source.to_numpy()&validpoly)
            if not len(si):continue
            tree=STRtree(allpolys[si]);_,dd=tree.query_nearest(flatvalid,return_distance=True,all_matches=False)
            ds=np.full(len(flat),np.nan);ds[pointvalid]=dd;ds=ds.reshape((len(c),ns))
            frame[label+'__minimum_boundary_distance_um']=np.nanmin(ds,axis=1)
            for threshold in CFG['boundary_apposition_um']:frame[label+f'__apposition_fraction_{threshold}um']=np.where(cmvalid,np.mean(ds<=threshold,axis=1),np.nan)
            conflict=source.reindex(c.index).to_numpy(bool)
            frame['source_role_conflict__'+label]=conflict
            # An independently called CM may be an EC under a different label
            # system. Never interpret its boundary against its own source polygon.
            for key in [label+'__minimum_boundary_distance_um']+[label+f'__apposition_fraction_{v}um' for v in CFG['boundary_apposition_um']]:frame.loc[conflict,key]=np.nan
            frame[label+'__nearby_profiles_50um']=np.asarray(radius@source.to_numpy(float)).ravel()
            area=np.array([p.area if p is not None else 0 for p in pp]);frame[label+'__nearby_cell_area_fraction_50um']=np.asarray(radius@np.where(source,area,0)).ravel()/(np.pi*50**2)
        sources={'EC':q.author_atlas_cell_state.isin(['cap.EC','cap.EC.High.Myo.','p.EC']),
            'FB':q.atlas_identity.eq('Fibroblast')&q.program_source_eligible,'Pericyte':q.author_atlas_cell_state.eq('Peri')}
        for label,source in sources.items():
            n=np.asarray(radius@source.to_numpy(float)).ravel();frame[label+'__sources_50um']=n
            programs=['EC_cycling','EC_guidance_proxy'] if label=='EC' else ['FB_fibrillar_matrix','FB_basement_matrix','FB_remodeling'] if label=='FB' else ['Pericyte_vascular_support']
            for program in programs:
                for view in ['native','core']:
                    total=np.asarray(radius@np.where(source,q[view+'__'+program],0)).ravel()
                    frame[label+'__'+view+'__'+program]=np.divide(total,n,out=np.full(len(n),np.nan),where=n>0)
            capture=np.asarray(radius@np.where(source,np.log1p(q.core_panel_depth),0)).ravel()
            frame[label+'__log_core_capture']=np.divide(capture,n,out=np.full(len(n),np.nan),where=n>0)
        ecsource=sources['EC'];n=np.asarray(radius@ecsource.to_numpy(float)).ravel()
        for view in ['native','core','nuclear']:
            vals=np.asarray(radius@np.where(ecsource,q[view+'__Ptprm'],0)).ravel()
            frame['EC__'+view+'__Ptprm']=np.divide(vals,n,out=np.full(len(n),np.nan),where=n>0)
        frame['CM__core__Ptprm']=c.core__Ptprm;frame['CM_EC_Ptprm_compatibility']=np.sqrt(np.maximum(frame.EC__nuclear__Ptprm,0)*c.core__Ptprm)
        frame['resolved_neighbors_50um']=np.asarray(radius@q.atlas_identity.isin(['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']).to_numpy(float)).ravel()
        save(frame,f'{stage}/vascular_interface_and_matrix_context.tsv.gz',index=True)
        del pp,allpolys,pts,flat,frame,q,c;gc.collect()
    save(pd.DataFrame(audits),'neighbor_graph_coverage.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['audit','prepare','geometry','all'],default='all');a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if a.step in ['audit','all']:source_audit()
    if a.step in ['prepare','all']:prepare()
    if a.step in ['geometry','all']:geometry()

if __name__=='__main__':main()
