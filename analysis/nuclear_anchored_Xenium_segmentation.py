"""Whole-section nuclear-anchored antibody-guided boundary refinement.

Preserves deposited profile IDs, including nucleus-free section profiles. Uses
no RNA identity, subtype, Ankrd1 or maturation value to build boundaries. This
refines existing cellular territories; it is not de novo cell-number discovery.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/nuclear-anchored-mpl')
import argparse
import gc
import gzip
import hashlib
import json
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from PIL import Image,ImageDraw
from scipy import ndimage as ndi, sparse
from skimage.draw import polygon
from skimage.segmentation import watershed
import tifffile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/xenium_nuclear_anchored_segmentation'
CACHE=Path('/workspace/heart-atlas/results/xenium_boundary_inputs')
RAW=Path('/workspace/heart-atlas/data/xenium')
ATLAS=ROOT/'results/whole_heart_spatial_atlas'
SOURCE=ROOT/'results/xenium_ventricular_mapping'
STAGES=['P0','P7','P14','P21']


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')


def raw_folder(stage):
    return next(RAW.glob('WT_'+('P07' if stage=='P7' else stage)+'__*'))


def segment(stage):
    folder=OUT/stage;folder.mkdir(parents=True,exist_ok=True)
    design={'scope':'Whole-section nuclear-anchored boundary refinement retaining original profile IDs',
        'nuclear_input':'Deposited DAPI-derived nuclear boundaries; individual label_id components retain their parent cell_id',
        'boundary_input':'Acquired composite ATP1A1/CD45/E-Cadherin channel',
        'analysis_grid':'First acquired image pyramid level, 0.425um, twice native pixel size',
        'territory_mask':'Native cell-polygon union dilated by 2um, including all seed pixels',
        'cost':'Boundary-channel intensity scaled to tissue 1st/99th percentiles, Gaussian sigma0.4um',
        'watershed_compactness':0,'fallback':'Nearest free original-interior pixel, search radii24/64/256 analysis pixels from centroid',
        'transcript_QV':20,'gene_targets':403,'boundary_band':'Two grid pixels (0.85um)',
        'nuclear_assignment':'Preserve deposited overlaps_nucleus assignments; audit disagreements with 2D projected nuclear ownership; do not force originally unassigned RNA into primary counts',
        'alternatives':['Native','Deposited nuclear-only','Projected nuclear-only','Nuclear-anchored whole cell','Nuclear plus cytoplasmic interior outside0.85um boundary band','Count-matched thinning'],
        'prohibited_tuning':'No RNA subtype, lineage, Ankrd1, maturation score or desired label count sets boundaries',
        'limits':'Not independent de novo nuclear detection or cell-number/merge/split discovery; no assumption that 2D projection resolves every z overlap'}
    if (OUT/'segmentation_design.json').exists():assert json.loads((OUT/'segmentation_design.json').read_text())==design
    else:write_json(OUT/'segmentation_design.json',design)
    if (folder/'segmentation_complete.json').exists():
        print('SKIP_SEGMENT',stage,flush=True);return
    start=time.time()
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
    raw=raw_folder(stage);px=float(json.loads((raw/'experiment.xenium').read_text())['pixel_size'])
    res=px*2
    with tifffile.TiffFile(CACHE/stage/'morphology_focus_0001.ome.tif',is_ome=False) as f:
        membrane=f.pages[0].pages[0].asarray(maxworkers=2)
        assert 'ATP1A1/CD45/E-Cadherin' in f.pages[0].description
    with tifffile.TiffFile(raw/'morphology_focus_0000.ome.tif',is_ome=False) as f:
        dapi=f.pages[0].pages[0].asarray(maxworkers=2)
    assert membrane.shape==dapi.shape
    shape=membrane.shape
    print('SEGMENT_IMAGES',stage,shape,res,flush=True)
    # The native polygon union limits refinement to observed cellular territory.
    image=Image.new('I',(shape[1],shape[0]),0);draw=ImageDraw.Draw(image)
    vertices=pd.read_csv(raw/'cell_boundaries.csv.gz',usecols=['cell_id','vertex_x','vertex_y'])
    lookup={c:i+1 for i,c in enumerate(obs.index)}
    for cell,part in vertices.groupby('cell_id',sort=False):
        label=lookup.get(cell)
        if label is not None and len(part)>=3:
            xy=np.rint(part[['vertex_x','vertex_y']].to_numpy()/res).astype(int)
            draw.polygon([tuple(x) for x in xy],fill=label)
    native=np.asarray(image)
    del vertices,draw
    nuclei=pd.read_parquet(CACHE/stage/'nucleus_boundaries.parquet')
    owner=np.zeros(shape,np.int32)
    for (cell,nucleus),part in nuclei.groupby(['cell_id','label_id'],sort=False):
        label=lookup.get(cell)
        if label is None or len(part)<3:continue
        xy=part[['vertex_x','vertex_y']].to_numpy()/res
        rr,cc=polygon(xy[:,1],xy[:,0],shape)
        old=owner[rr,cc]
        empty=old==0;conflict=(old>0)&(old!=label)
        owner[rr[empty],cc[empty]]=label
        owner[rr[conflict],cc[conflict]]=-1
    del nuclei
    # Seed every unique nuclear territory; multiple nuclei of one deposited
    # profile share one ID, preserving potentially binucleated cardiomyocytes.
    markers=np.maximum(owner,0)
    seeded=np.bincount(markers.ravel(),minlength=len(obs)+1)>0
    fallback=[];unseeded=[]
    for i,(cell,r) in enumerate(obs.iterrows(),1):
        if seeded[i]:continue
        y=int(round(r.y_centroid/res));x=int(round(r.x_centroid/res))
        y=np.clip(y,0,shape[0]-1);x=np.clip(x,0,shape[1]-1)
        for radius in [24,64,256]:
            y0=max(0,y-radius);y1=min(shape[0],y+radius+1);x0=max(0,x-radius);x1=min(shape[1],x+radius+1)
            ok=(native[y0:y1,x0:x1]==i)&(markers[y0:y1,x0:x1]==0)&(owner[y0:y1,x0:x1]>=0)
            yy,xx=np.where(ok)
            if len(yy):break
        if len(yy):
            j=np.argmin((yy+y0-y)**2+(xx+x0-x)**2);markers[yy[j]+y0,xx[j]+x0]=i
            fallback.append(i)
        else:unseeded.append(i)
    mask=ndi.binary_dilation(native>0,iterations=int(np.ceil(2/res)))
    mask|=markers>0
    positive=membrane[mask];lo,hi=np.quantile(positive,[.01,.99]);del positive
    cost=np.clip((membrane.astype(np.float32)-lo)/max(hi-lo,1),0,1)
    ndi.gaussian_filter(cost,sigma=.4/res,output=cost)
    print('WATERSHED_START',stage,'nuclear',int(seeded[1:].sum()),'fallback',len(fallback),'unseeded',len(unseeded),flush=True)
    labels=watershed(cost,markers=markers,mask=mask,compactness=0).astype(np.int32)
    del markers,cost,mask;gc.collect()
    # Conflicting projected nuclear territories are not forced to a 2D owner.
    # They remain separately flagged; the deposited nuclear-overlap flag is retained.
    border=np.zeros(shape,bool)
    border[1:]|=labels[1:]!=labels[:-1];border[:-1]|=labels[1:]!=labels[:-1]
    border[:,1:]|=labels[:,1:]!=labels[:,:-1];border[:,:-1]|=labels[:,1:]!=labels[:,:-1]
    border=ndi.binary_dilation(border,iterations=2)&(labels>0)
    np.savez_compressed(folder/'whole_section_segmentation.npz',cell_ids=obs.index.to_numpy(dtype=str),
        cell_labels=labels,nuclear_owner=owner,uncertain_boundary_band=border,
        grid_um=res,x0_um=0.,y0_um=0.,source_native_shape=np.array(shape))
    areas=np.bincount(labels.ravel(),minlength=len(obs)+1)[1:]*res**2
    nuclear_area=np.bincount(np.maximum(owner,0).ravel(),minlength=len(obs)+1)[1:]*res**2
    audit=obs[['x_centroid','y_centroid','cell_area','nucleus_area','atlas_identity']].copy()
    audit['segmented_area_um2']=areas;audit['projected_unique_nuclear_area_um2']=nuclear_area
    audit['seed_type']=np.where(seeded[1:],'Deposited DAPI-derived nuclear territory','Native-profile fallback')
    audit.loc[obs.index[np.array(unseeded,dtype=int)-1],'seed_type']='No seed at analysis resolution; retained as review record'
    audit['native_nuclear_area_missing']=obs.nucleus_area.isna()
    audit.to_csv(folder/'segmentation_profile_audit.tsv.gz',sep='\t',na_rep='NA')
    make_gallery(stage,obs,native,labels,owner,membrane,dapi,res)
    write_json(folder/'segmentation_complete.json',{'stage':stage,'profiles':len(obs),'grid_um':res,
        'shape_yx':list(shape),'nuclear_seeded_profiles':int(seeded[1:].sum()),
        'fallback_seeded_profiles':len(fallback),'unseeded_profiles':len(unseeded),
        'conflicting_projected_nuclear_pixels':int((owner<0).sum()),
        'boundary_band_um':2*res,'membrane_scaling_1_99':[float(lo),float(hi)],
        'nucleus_policy':'Deposited DAPI-derived nuclei anchor boundaries; multiple nuclei per profile retain same ID; nucleus-free profiles retain native-interior fallback seed',
        'RNA_used_to_optimize_segmentation':False,
        'scope':'Whole-section profile-retaining boundary refinement; does not discover missed cells or independently resolve vendor merge/split errors',
        'seconds':time.time()-start})
    print('SEGMENT_DONE',stage,time.time()-start,flush=True)


def make_gallery(stage,obs,native,labels,owner,membrane,dapi,res):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    rois=pd.read_csv(ROOT/'results/xenium_boundary_qc/selected_ROIs.tsv',sep='\t')
    with PdfPages(OUT/stage/f'{stage}_nuclear_segmentation_gallery.pdf') as pdf:
        for roi in rois[rois.stage.eq(stage)].itertuples():
            r0=max(0,int(roi.y_min/res));r1=min(labels.shape[0],int(roi.y_max/res)+1)
            c0=max(0,int(roi.x_min/res));c1=min(labels.shape[1],int(roi.x_max/res)+1)
            sl=np.s_[r0:r1,c0:c1];ext=[c0*res,c1*res,r1*res,r0*res]
            fig,axes=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
            for ax,grid,title in [(axes[0],native[sl],'Native boundaries'),(axes[1],labels[sl],'Nuclear-anchored membrane boundaries')]:
                dna=dapi[sl].astype(float);mem=membrane[sl].astype(float)
                rgb=np.zeros((*dna.shape,3));rgb[:,:,2]=np.clip(dna/max(np.quantile(dna,.995),1),0,1)
                rgb[:,:,1]=np.clip(mem/max(np.quantile(mem,.995),1),0,1)
                ax.imshow(rgb,extent=ext)
                edge=np.zeros(grid.shape,bool);edge[1:]|=grid[1:]!=grid[:-1];edge[:,1:]|=grid[:,1:]!=grid[:,:-1]
                ax.contour(edge,levels=[.5],colors=['white'],linewidths=.35,extent=ext,origin='upper')
                ax.set_title(title)
            axes[2].imshow(owner[sl]>0,extent=ext,cmap='Blues',vmin=0,vmax=1)
            diff=(native[sl]!=labels[sl])&(labels[sl]>0)
            overlay=np.zeros((*diff.shape,4));overlay[diff]=[1,.4,0,.6]
            axes[2].imshow(overlay,extent=ext);axes[2].set_title('Nuclei (blue); changed ownership (orange)')
            for ax in axes:ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
            fig.suptitle(stage+' · '+roi.roi_id+' · images define boundaries; RNA is not an optimization target')
            fig.savefig(OUT/stage/(roi.roi_id+'_nuclear_segmentation.png'),dpi=160);pdf.savefig(fig);plt.close(fig)


def count(stage):
    folder=OUT/stage
    if (folder/'count_complete.json').exists():print('SKIP_COUNT',stage,flush=True);return
    start=time.time()
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
    source=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
    genes=source.var_names.to_numpy(dtype=str);gindex=pd.Index(genes)
    with np.load(folder/'whole_section_segmentation.npz') as b:
        ids=b['cell_ids'].astype(str);labels=b['cell_labels'];owner=b['nuclear_owner'];band=b['uncertain_boundary_band'];res=float(b['grid_um'])
    np.testing.assert_array_equal(ids,obs.index)
    seg_audit=pd.read_csv(folder/'segmentation_profile_audit.tsv.gz',sep='\t',index_col=0)
    unresolved_geometry=seg_audit.segmented_area_um2.eq(0).reindex(obs.index).to_numpy()
    names=['native','nuclear_only','projected_nuclear_only','nuclear_anchored','high_confidence','expanded_unassigned_sensitivity']
    counts={name:np.zeros((len(obs),len(genes)),np.int32) for name in names}
    totals={'all_records':0,'Q20_gene_molecules':0,'native_assigned':0,'nuclear_anchored_assigned':0,
            'native_nuclear_molecules':0,'vendor_nuclear_ownership_preserved':0,'nonnuclear_changed_owner':0,
            'nuclear_vs_projected_owner_disagreement':0,'unassigned_because_projected_nuclear_conflict':0}
    path=CACHE/stage/'transcripts.parquet'
    parquet=pq.ParquetFile(path)
    for k,batch in enumerate(parquet.iter_batches(batch_size=300000,
            columns=['cell_id','feature_name','x_location','y_location','qv','overlaps_nucleus'])):
        f=batch.to_pandas();totals['all_records']+=len(f)
        f['cell_id']=f.cell_id.map(lambda x:x.decode() if isinstance(x,bytes) else x)
        f['feature_name']=f.feature_name.map(lambda x:x.decode() if isinstance(x,bytes) else x)
        ci=obs.index.get_indexer(f.cell_id);gi=gindex.get_indexer(f.feature_name)
        ok=f.qv.ge(20).to_numpy()&(gi>=0)
        ci=ci[ok];gi=gi[ok];f=f.loc[ok]
        totals['Q20_gene_molecules']+=len(f)
        nuclear=f.overlaps_nucleus.fillna(False).astype(bool).to_numpy()&(ci>=0)
        rr=np.rint(f.y_location.to_numpy()/res).astype(int);cc=np.rint(f.x_location.to_numpy()/res).astype(int)
        inside=(rr>=0)&(cc>=0)&(rr<labels.shape[0])&(cc<labels.shape[1])
        proposed=np.full(len(f),-1,np.int32);nowner=np.full(len(f),-1,np.int32);edge=np.ones(len(f),bool)
        proposed[inside]=labels[rr[inside],cc[inside]]-1
        nowner[inside]=owner[rr[inside],cc[inside]]-1
        edge[inside]=band[rr[inside],cc[inside]]
        disagreement=nuclear&(nowner>=0)&(nowner!=ci)
        totals['nuclear_vs_projected_owner_disagreement']+=int(disagreement.sum())
        # Preserve the deposited intranuclear assignment rather than overruling
        # it from a coarser 2D projection alone. This is source evidence, not truth.
        anchored=proposed.copy();anchored[nuclear]=ci[nuclear]
        projected_conflict=nowner==-2
        anchored[projected_conflict&~nuclear]=-1
        ambiguous=np.zeros(len(ci),bool);known=ci>=0
        ambiguous[known]=unresolved_geometry[ci[known]]
        # A profile wholly covered by another nucleus in 2D cannot receive an
        # honest unique image region without violating that nuclear ownership.
        # Retain its source RNA and mark geometry unresolved instead of erasing it.
        anchored[ambiguous]=ci[ambiguous];edge[ambiguous]=True
        expanded=anchored.copy()
        # The primary refinement does not promote vendor-unassigned molecules.
        anchored[ci<0]=-1
        totals['unassigned_because_projected_nuclear_conflict']+=int((projected_conflict&~nuclear).sum())
        destinations={'native':np.where(ci>=0,ci,-1),'nuclear_only':np.where(nuclear,ci,-1),
            'projected_nuclear_only':np.where(nowner>=0,nowner,-1),'nuclear_anchored':anchored,
            'expanded_unassigned_sensitivity':expanded,
            'high_confidence':np.where(nuclear|~edge,anchored,-1)}
        for name,dest in destinations.items():
            valid=dest>=0;np.add.at(counts[name],(dest[valid],gi[valid]),1)
        totals['native_assigned']+=int((ci>=0).sum())
        totals['nuclear_anchored_assigned']+=int((anchored>=0).sum())
        totals['native_nuclear_molecules']+=int(nuclear.sum())
        totals['vendor_nuclear_ownership_preserved']+=int(np.sum(anchored[nuclear]==ci[nuclear]))
        totals['nonnuclear_changed_owner']+=int((~nuclear&(ci>=0)&(anchored!=ci)).sum())
        if k%25==0:print('COUNT_BATCH',stage,k,totals['all_records'],flush=True)
    expected=source.X.toarray().astype(np.int32)
    assert np.array_equal(counts['native'],expected),'Q20 source reconstruction must be exact'
    assert (counts['nuclear_only']<=counts['native']).all()
    assert (counts['high_confidence']<=counts['nuclear_anchored']).all()
    assert totals['native_nuclear_molecules']==totals['vendor_nuclear_ownership_preserved']
    rng=np.random.default_rng(642)
    # Paired depth control isolates spatial selection from merely losing molecules.
    thin=np.zeros_like(counts['nuclear_anchored'])
    targets=counts['high_confidence'].sum(1)
    for i,total in enumerate(targets):
        if total:thin[i]=rng.multivariate_hypergeometric(counts['nuclear_anchored'][i],int(total),method='marginals')
    counts['depth_matched']=thin
    assert np.array_equal(counts['depth_matched'].sum(1),counts['high_confidence'].sum(1))
    nuclear_thin=np.zeros_like(counts['native']);nuclear_rng=np.random.default_rng(643)
    for i,total in enumerate(counts['nuclear_only'].sum(1)):
        if total:nuclear_thin[i]=nuclear_rng.multivariate_hypergeometric(counts['native'][i],int(total),method='marginals')
    counts['nuclear_depth_matched']=nuclear_thin
    info=obs[['x_centroid','y_centroid','cell_area','nucleus_area','atlas_identity']].copy()
    info['stage']=stage
    info['geometry_unresolved_native_counts_retained']=unresolved_geometry
    obj=ad.AnnData(X=sparse.csr_matrix(counts['nuclear_anchored']),obs=info,var=pd.DataFrame(index=genes))
    for name in ['nuclear_only','high_confidence','depth_matched','nuclear_depth_matched']:obj.layers[name]=sparse.csr_matrix(counts[name])
    obj.uns['method']='Nuclear-anchored antibody-guided whole-section profile-retaining refinement; broad labels are native annotations, not revalidated identities'
    obj.uns['grid_um']=res
    obj.write_h5ad(folder/'nuclear_anchored_counts.h5ad',compression='gzip')
    # Sparse projected-nuclear control is separate from deposited nuclear overlap.
    sparse.save_npz(folder/'projected_nuclear_counts.npz',sparse.csr_matrix(counts['projected_nuclear_only']))
    sparse.save_npz(folder/'expanded_unassigned_sensitivity_counts.npz',sparse.csr_matrix(counts['expanded_unassigned_sensitivity']))
    per=pd.DataFrame({name+'_gene_molecules':values.sum(1) for name,values in counts.items()},index=obs.index)
    per.index.name='cell_id';per.to_csv(folder/'transcript_assignment_counts.tsv.gz',sep='\t')
    totals.update({'stage':stage,'profiles':len(obs),'genes':len(genes),'Q20_native_matrix_reconstructed_exactly':True,
        'profiles_with_unresolved_2D_geometry_retaining_native_counts':int(unresolved_geometry.sum()),
        'all_vendor_nuclear_transcripts_keep_native_owner':True,'seconds':time.time()-start,
        'source_transcript_file_sha256':file_digest(path),
        'interpretation':'Nuclear overlap is deposited vendor evidence; 2D polygon ownership disagreements are audited, not treated as proof of an incorrect source assignment'})
    write_json(folder/'count_complete.json',totals)
    print('COUNT_DONE',stage,totals,flush=True)


def file_digest(path):
    with path.open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def evict_cache(stage):
    folder=OUT/stage
    assert (folder/'segmentation_complete.json').exists() and (folder/'count_complete.json').exists()
    completion=json.loads((folder/'count_complete.json').read_text())
    assert completion['Q20_native_matrix_reconstructed_exactly']
    records=[]
    for name in ['morphology_focus_0001.ome.tif','transcripts.parquet']:
        path=CACHE/stage/name
        if not path.exists():continue
        provenance=json.loads((CACHE/stage/(name+'.acquisition.json')).read_text())
        assert provenance['CRC_verified'] and path.stat().st_size==provenance['bytes']
        digest=file_digest(path)
        if name=='transcripts.parquet':assert digest==completion['source_transcript_file_sha256']
        records.append({'cache_file':str(path),'bytes':path.stat().st_size,'sha256':digest,'acquisition':provenance})
    if records:
        write_json(folder/'source_cache_provenance.json',{'eviction_reason':'Conserve workspace disk after exact source reconstruction and derived whole-section outputs',
            'source_DAPI_native_polygons_and_original_count_atlas_retained':True,'members':records})
        for record in records:Path(record['cache_file']).unlink()
    print('VERIFIED_CACHE_EVICTED',stage,len(records),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=STAGES)
    p.add_argument('--part',required=True,choices=['segment','count','evict-cache']);args=p.parse_args()
    if args.part=='segment':segment(args.stage)
    elif args.part=='count':count(args.stage)
    else:evict_cache(args.stage)
