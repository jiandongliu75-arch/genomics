"""Observed Xenium sender features in source-disjoint spatial domains."""
import os,gc,json,time,argparse,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/joint-signaling-mpl')
import h5py
import numpy as np
import pandas as pd
from scipy import sparse,ndimage
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from anndata.io import sparse_dataset,read_elem
from joint_signaling_common import *
from joint_signaling_transport import collective,equivalence_check

def source_groups(f,fb):
    group=pd.Series('Exclude',index=f.index)
    for name in ['Pericyte','Smooth muscle','Immune','Epicardial']:
        group.loc[f.atlas_identity.eq(name)&f.program_source_eligible]=name
    ec=f.atlas_identity.eq('Endothelial')&f.program_source_eligible;group.loc[ec]='EC:Unresolved'
    for name in EC:group.loc[ec&f.EC_nuclear_type.eq(name)]='EC:'+name
    eligible=fb.program_source_eligible;group.loc[fb.index[eligible]]='FB:Unresolved'
    for name in FB:group.loc[fb.index[eligible&fb.recommended_candidate.eq(name+'-associated')]]='FB:'+name
    cm=f.atlas_identity.eq('Cardiomyocyte')&f.pass_expression_QC&f.core_transcripts.ge(20)
    group.loc[cm]='CM'
    return group

def graph(xy,tissue,resolution):
    t=cKDTree(xy);g=t.sparse_distance_matrix(t,200,output_type='coo_matrix')
    take=g.row!=g.col;row=g.row[take];col=g.col[take];dist=g.data[take].astype(np.float32)
    before=len(dist);keep=np.ones(before,bool);mask=ndimage.binary_dilation(tissue,iterations=1)
    # Screen straight paths against the observed tissue; no label or outcome is used.
    for lo in range(0,len(dist),250000):
        hi=min(lo+250000,len(dist));a=xy[row[lo:hi]];b=xy[col[lo:hi]];ok=np.ones(hi-lo,bool)
        for w in [.15,.3,.5,.7,.85]:
            q=np.rint(((1-w)*a+w*b)[:,::-1]/resolution).astype(int)
            q=np.clip(q,0,np.array(mask.shape)-1);ok &= mask[q[:,0],q[:,1]]
        keep[lo:hi]=ok
    row=row[keep];col=col[keep];dist=dist[keep]
    M=sparse.csr_matrix((np.maximum(dist,.001),(row,col)),shape=(len(xy),len(xy)))
    return M,{'candidate_spatial_edges':before,'retained_tissue_screened_edges':M.nnz}

def domains(stage,f,d):
    model=KMeans(n_clusters=4,n_init=20,random_state=SEED+int(stage[1:])).fit(d[['x_centroid','y_centroid']].to_numpy())
    centers=model.cluster_centers_;xy=f[['x_centroid','y_centroid']].to_numpy()
    dist2=((xy[:,None,:]-centers[None,:,:])**2).sum(2);label=dist2.argmin(1)
    margin=np.full(len(f),np.inf)
    for j in range(4):
        take=label==j
        for k in range(4):
            if j==k:continue
            value=(dist2[take,k]-dist2[take,j])/(2*np.linalg.norm(centers[j]-centers[k]))
            margin[take]=np.minimum(margin[take],value)
    return label,margin,centers

def base_features(d):
    x=pd.DataFrame(index=d.index)
    for col,label in [('CM_outcome_excluded_depth','CM_capture'),('cell_area','CM_area'),('neighbors_resolved','neighbors'),
                      ('arterial_distance_um','arterial_distance'),('endocardial_proxy_distance_um','endocardial_distance'),
                      ('distance_to_large_void_um','void_distance'),('distance_to_outer_um','outer_distance')]:x['A__log_'+label]=np.log1p(d[col])
    x['A__nuclear_area_missing']=d.nucleus_area.isna().astype(float)
    for col in ['x_centroid','y_centroid']:x['A__'+col]=d[col]/1000
    for col in d:
        if col.startswith('broad__'):x['A__'+col]=d[col]
    x=pd.concat([x,pd.get_dummies(d.CM_state,prefix='A__CM_state',dtype=float),pd.get_dummies(d.morphology_region,prefix='A__morphology',dtype=float)],axis=1)
    return x

def prepare(stage,views):
    begin=time.monotonic();cfg=config()
    f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
    d=read(f'results/cardiac_subtype_niches/{stage}_proximity50_CM_subtype_neighborhoods.tsv.gz',index_col='cell_id')
    c=read(f'results/CM_maturation21_nonCM_context/{stage}/CM_maturation21_scores.tsv.gz',index_col='cell_id')
    fb=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
    selected=read('results/joint_RNA_Xenium_signaling/selected_RNA_supported_pairs.tsv')
    support=read('results/joint_RNA_Xenium_signaling/RNA_pair_stage_support.tsv')
    pair_support=read('results/joint_RNA_Xenium_signaling/RNA_pair_source_support.tsv.gz')
    eligible_pairs=set(support.loc[support.stage.eq(stage)&support.RNA_supported,'pair_id'])
    selected=selected[selected.pair_id.isin(eligible_pairs)].reset_index(drop=True);selected['stage_pair_index']=np.arange(len(selected))
    save(selected,stage+'/stage_pairs.tsv')
    group=source_groups(f,fb);take=group.ne('Exclude')
    # Identical source cells in deposited/core comparisons, independent of outcomes.
    take &= f.core_transcripts.ge(20)
    source_ids=f.index[take];sources=f.loc[source_ids].copy();sources['source_group']=group.loc[source_ids]
    source_codes=pd.Categorical(sources.source_group,categories=GROUPS).codes
    assert (source_codes>=0).all()
    label,margin,centers=domains(stage,f,d);sources['domain']=label[f.index.get_indexer(source_ids)]
    d['domain']=label[f.index.get_indexer(d.index)];d['domain_boundary_margin_um']=margin[f.index.get_indexer(d.index)]
    d['analysis_eligible']=d.neighbors_resolved.ge(5)&d.domain_boundary_margin_um.ge(200)&d.index.isin(source_ids)
    save(d[['stage','CM_state','x_centroid','y_centroid','domain','domain_boundary_margin_um','analysis_eligible']],stage+'/receiver_domains.tsv.gz',index=True)
    save(sources[['stage','x_centroid','y_centroid','source_group','domain']],stage+'/source_domains.tsv.gz',index=True)
    (OUT/stage/'domain_centers.json').write_text(json.dumps(centers.tolist())+'\n')
    mask_path=ROOT/f'results/cardiac_spatial_reconstruction/{stage}/morphology_masks.npz';INPUTS.add(str(mask_path))
    with np.load(mask_path) as b:tissue=b['tissue'];resolution=float(b['resolution_um'])
    ligand_species=sorted(set(selected.ligand));receptor_species=sorted(set(selected.receptor))
    lig_i=np.array([ligand_species.index(x) for x in selected.ligand]);rec_i=np.array([receptor_species.index(x) for x in selected.receptor])
    radii=selected.radius_um.to_numpy(int);allgenes=sorted(set('_'.join(ligand_species+receptor_species).split('_')))
    source_names=GROUPS+['Broad:Endothelial','Broad:Fibroblast']
    supported=pair_support[pair_support.stage.eq(stage)&pair_support.joint_two_library_support]
    allowed_pairs=set(zip(supported.source_group,supported.pair_id))
    ligand_def=[];cot_def=[]
    for source in source_names:
        for ligand,radius in sorted(set(zip(selected.ligand,selected.radius_um))):
            pairs=selected[selected.ligand.eq(ligand)&selected.radius_um.eq(radius)]
            if not any((source,p) in allowed_pairs for p in pairs.pair_id):continue
            for stat in ['sum','mean']:
                ligand_def.append({'name':f'ligand__{source}__{ligand}__{radius}um__{stat}','source':source,'ligand':ligand,'radius_um':int(radius),'stat':stat})
        for p,row in selected.iterrows():
            if (source,row.pair_id) in allowed_pairs:
                cot_def.append({'name':f'cot__{source}__{row.pair_id}','source':source,'pair_id':row.pair_id,'pair_index':p,'pathway':row.pathway,'radius_um':int(row.radius_um)})
    feature_info=[]
    for row in ligand_def:feature_info.append({**row,'feature_group':'C' if row['source']=='CM' else 'D'})
    for row in cot_def:feature_info.append({**row,'feature_group':'C' if row['source']=='CM' else 'E'})
    for rec in receptor_species:feature_info.append({'name':'receptor__'+rec,'feature_group':'B','receptor':rec,'source':'Receiver CM'})
    for source in source_names:
        for radius in [50,200]:feature_info.append({'name':f'count__{source}__{radius}um','feature_group':'C_fine_context','source':source,'radius_um':radius})
    save(pd.DataFrame(feature_info),stage+'/feature_definitions.tsv')
    for view in views:
        path=(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
        INPUTS.add(str(path));print('JOINT_SPATIAL_READ',stage,view,len(sources),len(selected),flush=True)
        with h5py.File(path) as h:
            obsids=read_elem(h['obs']).index;genes=read_elem(h['var']).index
            ii=obsids.get_indexer(source_ids);assert (ii>=0).all();jj=genes.get_indexer(allgenes);assert (jj>=0).all()
            x=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[ii].tocsr()
            denom=np.asarray(x[:,~genes.isin(OUTCOME_GENES)].sum(1)).ravel()
            counts=x[:,jj].toarray().astype(np.float64);del x
        values=counts*1e4/np.maximum(denom,1)[:,None]
        S=np.column_stack([complexes(values,allgenes,s) for s in ligand_species])
        D=np.column_stack([complexes(values,allgenes,s) for s in receptor_species])
        X=base_features(d);X=X.reindex(columns=X.columns.tolist()+[r['name'] for r in feature_info],fill_value=0.)
        source_profiles=[];transport=[];cm_totals=[]
        # Each solve contains all eligible cardiac lineages within one validation domain.
        for domain in range(4):
            domain_start=time.monotonic();si=np.flatnonzero(sources.domain.eq(domain));local=sources.iloc[si]
            xy=local[['x_centroid','y_centroid']].to_numpy();M,geometry=graph(xy,tissue,resolution)
            ci=d.index.get_indexer(local.index);iscm=ci>=0;receiver_loc=np.flatnonzero(iscm);receiver_global=ci[iscm]
            print('JOINT_COT_START',stage,view,domain,len(local),M.nnz,len(receiver_loc),flush=True)
            localS=S[si];localD=D[si];sourcecode=source_codes[si]
            checkpoint=OUT/stage/f'{view}_domain{domain}_checkpoint.npz'
            fingerprint=hashlib.sha256()
            for a in [localS,localD,M.data,M.indices,M.indptr]:fingerprint.update(a.tobytes())
            fingerprint.update((''.join(X.columns)+''.join(local.index)+json.dumps(cfg,sort_keys=True)+
                sha(Path(__file__))+sha(Path(__file__).with_name('joint_signaling_transport.py'))).encode())
            fingerprint=fingerprint.hexdigest()
            if checkpoint.exists():
                with np.load(checkpoint) as z:
                    if str(z['fingerprint'])==fingerprint:
                        X.iloc[receiver_global]=z['features']
                        transport.append(json.loads(str(z['diagnostics'])))
                        source_profiles.extend(json.loads(str(z['source_profiles'])))
                        print('JOINT_COT_CHECKPOINT_RESTORED',stage,view,domain,flush=True)
                        del M,localS,localD;gc.collect();continue
            profile_start=len(source_profiles)
            membership=np.zeros((len(si),len(source_names)),np.float64)
            membership[np.arange(len(si)),sourcecode]=1
            membership[:,source_names.index('Broad:Endothelial')]=local.source_group.str.startswith('EC:').to_numpy()
            membership[:,source_names.index('Broad:Fibroblast')]=local.source_group.str.startswith('FB:').to_numpy()
            receptor_values=np.log1p(localD[receiver_loc])
            X.loc[d.index[receiver_global],['receptor__'+r for r in receptor_species]]=receptor_values
            kernel_cache={}
            for radius in [50,200]:
                K=M.copy();K.data=np.where(K.data<=radius,np.exp(-K.data/(radius/2)),0);K.eliminate_zeros()
                count=K.copy();count.data=np.ones_like(count.data)
                n=np.asarray(count.T@membership)
                for j,source in enumerate(source_names):X.loc[d.index[receiver_global],f'count__{source}__{radius}um']=np.log1p(n[receiver_loc,j])
                kernel_cache[radius]=(K,np.asarray(K.T@membership))
            for ligand,radius in sorted(set((r['ligand'],r['radius_um']) for r in ligand_def)):
                K,den=kernel_cache[radius];li=ligand_species.index(ligand)
                sums=np.asarray(K.T@(membership*localS[:,li,None]))
                for info in ligand_def:
                    if info['ligand']!=ligand or info['radius_um']!=radius:continue
                    j=source_names.index(info['source']);v=sums[:,j] if info['stat']=='sum' else sums[:,j]/np.maximum(den[:,j],1e-8)
                    X.loc[d.index[receiver_global],info['name']]=np.log1p(v[receiver_loc])
            plan,diag=collective(localS,localD,M,lig_i,rec_i,radii,maxiter=cfg['cot_max_iterations'],tolerance=cfg['cot_stop'])
            if not diag['converged']:
                print('JOINT_COT_EXTEND',stage,view,domain,diag,flush=True)
                plan,diag=collective(localS,localD,M,lig_i,rec_i,radii,maxiter=15000,tolerance=cfg['cot_stop'])
            if not diag['converged'] or diag.get('maximum_relative_supply_excess',0)>.01:
                raise RuntimeError(('Collective solver did not meet convergence/capacity criteria',stage,view,domain,diag))
            transport.append({'stage':stage,'view':view,'domain':domain,'source_profiles':len(local),'CM_profiles':len(receiver_loc),**geometry,**diag})
            if plan is not None:
                # Aggregate only after solving jointly; overlapping broad/fine descriptors do not duplicate the transport budget.
                cm_map=np.full(len(si),-1,np.int32);cm_map[receiver_loc]=np.arange(len(receiver_loc))
                edge_cm=cm_map[plan['receiver_cell']];valid=edge_cm>=0
                edge_src=sourcecode[plan['source_cell'][valid]];edge_pair=plan['pair_index'][valid];edge_mass=plan['mass'][valid]
                key=(edge_cm[valid]*len(GROUPS)+edge_src)*len(selected)+edge_pair
                incoming=np.bincount(key,weights=edge_mass,minlength=len(receiver_loc)*len(GROUPS)*len(selected)).reshape(len(receiver_loc),len(GROUPS),len(selected))
                for info in cot_def:
                    source=info['source'];p=info['pair_index']
                    if source.startswith('Broad:'):
                        take=[j for j,g in enumerate(GROUPS) if broad(g)==source[6:]];v=incoming[:,take,p].sum(1)
                    else:v=incoming[:,GROUPS.index(source),p]
                    X.loc[d.index[receiver_global],info['name']]=np.log1p(v)
                    source_profiles.append({'stage':stage,'view':view,'domain':domain,'source':source,'pair_id':info['pair_id'],
                        'pathway':info['pathway'],'received_mass_sum':v.sum(),'CM_profiles':len(v),'CMs_with_received_mass':int((v>1e-12).sum())})
            print('JOINT_COT_DONE',stage,view,domain,round(time.monotonic()-domain_start,1),diag['iterations'],flush=True)
            np.savez_compressed(checkpoint,fingerprint=np.array(fingerprint),features=X.iloc[receiver_global].to_numpy(np.float32),
                diagnostics=np.array(json.dumps(transport[-1])),source_profiles=np.array(json.dumps(source_profiles[profile_start:])))
            del M,plan,localS,localD;gc.collect()
        # Target scores use the frozen direct Xenium definition; no inferred RNA outcome is substituted.
        ykeys=[view+'_Z__'+k for k in ['maturation21','metabolic','sarcomere','excitation_contraction','cycling']]
        y=c.loc[d.index,ykeys].to_numpy(float)
        finite=np.isfinite(X.to_numpy(float)).all(1)&np.isfinite(y).all(1)
        eligible=d.analysis_eligible.to_numpy()&finite
        np.savez_compressed(OUT/stage/(view+'_predictor_features.npz'),cell_ids=d.index.to_numpy(str),columns=X.columns.to_numpy(str),
            features=X.to_numpy(np.float32),outcomes=y.astype(np.float32),outcome_names=np.array(ykeys),eligible=eligible,
            domain=d.domain.to_numpy(np.int8),CM_states=d.CM_state.to_numpy(str),xy=d[['x_centroid','y_centroid']].to_numpy(np.float32))
        save(pd.DataFrame(transport),stage+'/'+view+'_transport_diagnostics.tsv')
        save(pd.DataFrame(source_profiles),stage+'/'+view+'_source_pair_communication.tsv')
        for checkpoint in (OUT/stage).glob(view+'_domain*_checkpoint.npz'):checkpoint.unlink()
        print('JOINT_SPATIAL_VIEW_COMPLETE',stage,view,'eligible receivers',int(eligible.sum()),'features',X.shape[1],flush=True)
    export_inputs();print('JOINT_SPATIAL_COMPLETE',stage,round(time.monotonic()-begin,1),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=STAGES);p.add_argument('--views',nargs='+',choices=['core','native'],default=['core','native']);p.add_argument('--solver-check',action='store_true');p.add_argument('--resume',action='store_true');args=p.parse_args()
    if args.solver_check:
        result=equivalence_check();(OUT/'transport_equivalence_check.json').write_text(json.dumps(result,indent=2)+'\n');print('SOLVER_EQUIVALENCE_PASS',result)
    else:
        for stage in ([args.stage] if args.stage else STAGES):
            views=[v for v in args.views if not args.resume or not all((OUT/stage/(v+'_'+suffix)).exists() for suffix in ['predictor_features.npz','transport_diagnostics.tsv','source_pair_communication.tsv'])]
            if views:prepare(stage,views)

if __name__=='__main__':main()
