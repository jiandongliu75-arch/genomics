"""Outcome-blinded spatial atlas/niches and staged non-CM association pipeline.

Observed counts and inherited identities remain immutable. All geometry is
section-specific. Optional learned niche backends never determine CM outcomes.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-pipeline-mpl')
from pathlib import Path
import argparse, hashlib, json, time, gc, platform, subprocess
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler, SplineTransformer
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, r2_score
from sklearn.linear_model import Ridge

ROOT=Path(__file__).resolve().parents[1]
CONFIG_PATH=Path(os.environ.get('CARDIAC_SPATIAL_CONFIG',str(ROOT/'analysis/cardiac_spatial_pipeline.json'))).resolve()
CFG=json.loads(CONFIG_PATH.read_text())
OUT=ROOT/CFG['output']
STAGES=CFG['stages']
LINEAGES=['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
NC=LINEAGES[1:]
STATUSES=['Mixed-marker review','Ambiguous/review','QC/information flag']
TYPES=LINEAGES+STATUSES
MAT=pd.read_csv(ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',sep='\t')
PROGS=list(CFG['programs'])
FBF=['Postn-associated','Col8a1-associated','Gfpt2-associated','Lsamp-associated','Kirrel3-associated']
ECF=['Capillary-like','Endocardial-like','Arterial-like','Venous-like','Lymphatic-like']
MARKERS={
 'Cardiomyocyte':'Myh6 Myh7 Tnni3 Myl2 Actc1 Tnnt2'.split(),
 'Endothelial':'Pecam1 Cdh5 Esam Tek Erg Fli1 Kdr Flt1 Tie1'.split(),
 'Fibroblast':'Dcn Lum Pdgfra Tcf21'.split(),
 'Pericyte':'Rgs5 Cspg4 Pdgfrb Mcam'.split(),
 'Smooth muscle':'Myh11 Cnn1 Tagln Acta2'.split(),
 'Immune':'Ptprc Csf1r Mertk Mrc1 Cd22 Cd28 Cd8b1'.split(),
 'Epicardial':'Upk3b Wt1'.split()}

def read(rel,**kwargs):
    kwargs.setdefault('index_col',0)
    return pd.read_csv(ROOT/rel,sep='\t',**kwargs)
def save(d,rel,index=False):
    p=OUT/rel;p.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(p,sep='\t',index=index,na_rep='NA',float_format='%.9g')
def dump(rel,d):
    p=OUT/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,default=str)+'\n')
def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def score(x,genes,panel):
    use=np.isin(genes,panel);assert use.sum()==len(panel)
    den=np.asarray(x.sum(1)).ravel()
    return np.log1p(x[:,use].toarray()*1e4/np.maximum(den,1)[:,None]).mean(1)
def normalize_rows(m):
    n=np.asarray(m.sum(1)).ravel()
    return sparse.diags(1/np.maximum(n,1))@m,n
def graph(stage,name):return sparse.load_npz(OUT/stage/(name+'_graph.npz')).tocsr()
def frame(stage):return read(f"{CFG['output']}/{stage}/atlas_annotations.tsv.gz")
def simplify(s):return s.str.replace(' candidate','',regex=False).replace({'Endothelial subtype unresolved':'Unresolved'})

def atlas():
    sources=set();counts=[];program_def=[]
    for stage in STAGES:
        print('ATLAS',stage,flush=True)
        a=read(f'results/whole_heart_spatial_atlas/{stage}/all_cell_annotations.tsv.gz')
        p=ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad';sources.add(p)
        b=ad.read_h5ad(p);assert list(b.obs_names)==list(a.index)
        genes=b.var_names.to_numpy(dtype=str);core=b.layers['high_confidence'].tocsr();nuc=b.layers['nuclear_only'].tocsr()
        audit=read(f'results/xenium_nuclear_anchored_segmentation/{stage}/segmentation_profile_audit.tsv.gz')
        f=a[['stage','stage_cell_id','x_centroid','y_centroid','cell_area','nucleus_area','atlas_identity','atlas_resolved',
             'ventricular_CM_strict','author_ventricular','author_atlas_cell_state','pass_expression_QC']].copy()
        f['geometry_unseeded']=audit.seed_type.str.startswith('No seed')
        f['segmentation_seed_type']=audit.seed_type
        f['refined_area_um2']=audit.segmented_area_um2
        f['core_transcripts']=np.asarray(core.sum(1)).ravel();f['nuclear_transcripts']=np.asarray(nuc.sum(1)).ravel()
        f['core_detected_genes']=np.asarray((core>0).sum(1)).ravel()
        f['fine_identity']='Unresolved';f['fine_identity_evidence']='Unresolved'
        fb=read(f'results/fibroblast_spatial_consensus/{stage}/fibroblast_consensus_annotations.tsv.gz')
        f.loc[fb.index,'fine_identity']=fb.high_confidence_family.where(~fb.high_confidence_family.str.startswith('Unresolved'),'Unresolved')
        f.loc[fb.index,'fine_identity_evidence']=fb.confidence_tier
        ec=read(f'results/endothelial_QC_subtypes/Xenium/{stage}/EC_consensus_annotations.tsv.gz')
        e=ec.loc[ec.native_EC_broad]
        f.loc[e.index,'fine_identity']=simplify(e.fine_consensus_family).where(e.high_confidence_subtype,'Unresolved')
        f.loc[e.index,'fine_identity_evidence']=e.confidence_tier
        fine=read(f'results/ventricular_spatial_resolution/{stage}/fine_CM_spatial_annotations.tsv.gz')
        cmkeep=fine.index.intersection(f.index[f.atlas_identity.eq('Cardiomyocyte')])
        f.loc[cmkeep,'fine_identity']=fine.loc[cmkeep,'fine_group_proposal'].where(fine.loc[cmkeep,'prioritized_fine_candidate'],'Unresolved')
        f.loc[cmkeep,'fine_identity_evidence']=np.where(fine.loc[cmkeep,'prioritized_fine_candidate'],'Prioritized RNA-derived CM candidate','Unresolved')
        f['EC_nuclear_type']='Unresolved'
        f.loc[e.index,'EC_nuclear_type']=simplify(e.nuclear_supported_family)
        # Source gating is explicit, without removing profiles from the atlas.
        f['program_source_eligible']=False
        for l in NC:
            m=f.atlas_identity.eq(l).to_numpy()
            gi=np.isin(genes,MARKERS[l]);anchors=np.asarray((core[:,gi]>0).sum(1)).ravel()
            ok=m&(anchors>=1)&(f.core_transcripts.to_numpy()>=20)&~f.geometry_unseeded.to_numpy()
            if l=='Fibroblast':
                blocked=fb.mixed_lineage_compartment_review|fb.supported_family_conflict
                ok[f.index.get_indexer(fb.index[blocked])]=False
            elif l=='Endothelial':
                ids=e.index[e.nuclear_broad_EC_clean&~e.multiple_projected_nuclear_components_review&~e.geometry_unresolved_review]
                ok &= f.index.isin(ids)
            f.loc[ok,'program_source_eligible']=True
        f['source_log_capture']=np.log1p(f.core_transcripts)
        for name,cfg in CFG['programs'].items():
            panel=[g for g in cfg['genes'] if g in genes]
            assert len(panel)>=3,(name,panel)
            assert not set(panel)&set(MAT.gene)
            x=nuc if cfg['view']=='nuclear' else core
            val=score(x,genes,panel)
            take=f.atlas_identity.eq(cfg['lineage'])&f.program_source_eligible
            f[name]=np.where(take,val,np.nan)
            if stage==STAGES[0]:program_def.append({'program':name,'lineage':cfg['lineage'],'count_view':cfg['view'],
                'measured_genes':';'.join(panel),'genes':len(panel),'interpretation':cfg['meaning'],'imputed':False})
        # Fixed outcome; no maturity values participate in graph/niche construction.
        old=read(f'results/cardiac_spatial_maturation_EC/{stage}/ventricular_CM_maturation.tsv.gz')
        cmids=old.index;ix=f.index.get_indexer(cmids)
        measured=np.isin(genes,MAT.gene)
        corelog=np.log1p(core[ix].toarray()*1e4/np.maximum(f.core_transcripts.iloc[ix].to_numpy(),1)[:,None])
        f['CM_maturation_native']=np.nan;f.loc[cmids,'CM_maturation_native']=old.maturation_21gene_mean
        f['CM_maturation_core']=np.nan;f.loc[cmids,'CM_maturation_core']=corelog[:,measured].mean(1)
        z=[]
        for category,part in MAT.groupby('figure_category'):
            values=corelog[:,np.isin(genes,part.gene)].mean(1)
            f['CM_'+category]=np.nan;f.loc[cmids,'CM_'+category]=values
            z.append((values-values.mean())/max(values.std(),1e-8))
        f['CM_category_balanced']=np.nan;f.loc[cmids,'CM_category_balanced']=np.mean(z,axis=0)
        f['CM_without_Myl2']=np.nan;f.loc[cmids,'CM_without_Myl2']=corelog[:,measured&(genes!='Myl2')].mean(1)
        f['CM_cycling']=np.nan;f.loc[cmids,'CM_cycling']=score(nuc[ix],genes,['Mki67','Top2a','Cenpf'])
        f['CM_Slit2']=np.nan;f.loc[cmids,'CM_Slit2']=corelog[:,list(genes).index('Slit2')]
        f['CM_outcome_excluded_depth']=np.nan;f.loc[cmids,'CM_outcome_excluded_depth']=np.asarray(core[ix][:,~measured].sum(1)).ravel()
        for col in ['arterial_distance_um','endocardial_proxy_distance_um']:
            f[col]=np.nan;f.loc[cmids,col]=old[col]
        save(f,f'{stage}/atlas_annotations.tsv.gz',index=True)
        for l,zg in f.groupby('atlas_identity',observed=True):
            counts.append({'stage':stage,'identity':l,'profiles':len(zg),'program_eligible':int(zg.program_source_eligible.sum()),
                'fine_resolved':int(zg.fine_identity.ne('Unresolved').sum()),'nucleus_missing':int(zg.nucleus_area.isna().sum())})
        for rel in [f'results/whole_heart_spatial_atlas/{stage}/all_cell_annotations.tsv.gz',
                    f'results/fibroblast_spatial_consensus/{stage}/fibroblast_consensus_annotations.tsv.gz',
                    f'results/endothelial_QC_subtypes/Xenium/{stage}/EC_consensus_annotations.tsv.gz',
                    f'results/ventricular_spatial_resolution/{stage}/fine_CM_spatial_annotations.tsv.gz',
                    f'results/cardiac_spatial_maturation_EC/{stage}/ventricular_CM_maturation.tsv.gz',
                    f'results/xenium_nuclear_anchored_segmentation/{stage}/whole_section_segmentation.npz']:
            sources.add(ROOT/rel)
        del b,core,nuc,corelog;gc.collect()
    save(pd.DataFrame(counts),'atlas_counts.tsv');save(pd.DataFrame(program_def),'program_definitions.tsv')
    dump('source_files.json',sorted(str(p.relative_to(ROOT)) for p in sources))

def contact_pairs(labels,n):
    keys=[];weights=[]
    for start in range(0,labels.shape[0],256):
        a=labels[start:start+256]
        pairs=[(a[:,:-1],a[:,1:])]
        stop=min(start+256,labels.shape[0]-1)
        if stop>start:pairs.append((labels[start:stop],labels[start+1:stop+1]))
        for left,right in pairs:
            m=(left>0)&(right>0)&(left!=right)
            l=left[m].astype(np.int64)-1;r=right[m].astype(np.int64)-1
            k=np.minimum(l,r)*n+np.maximum(l,r)
            k,c=np.unique(k,return_counts=True);keys.append(k);weights.append(c)
    allk=np.concatenate(keys);allw=np.concatenate(weights)
    u,inv=np.unique(allk,return_inverse=True);counts=np.bincount(inv,weights=allw)
    return u//n,u%n,counts

def geometry():
    totals=[]
    for stage in STAGES:
        print('GEOMETRY',stage,flush=True)
        f=frame(stage);n=len(f);xy=f[['x_centroid','y_centroid']].to_numpy()
        with np.load(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/whole_section_segmentation.npz') as b:
            assert np.array_equal(b['cell_ids'].astype(str),f.index)
            labels=b['cell_labels'];res=float(b['grid_um'])
        i,j,pixels=contact_pairs(labels,n);length=pixels*res
        use=length>=CFG['contact_minimum_um']-1e-8
        a=sparse.csr_matrix((np.ones(use.sum()*2),(np.r_[i[use],j[use]],np.r_[j[use],i[use]])),shape=(n,n))
        sparse.save_npz(OUT/stage/'contact_graph.npz',a)
        save(pd.DataFrame({'profile_i':i[use],'profile_j':j[use],'shared_boundary_um':length[use]}),f'{stage}/contact_edges.tsv.gz')
        # A sample of the refined territory union screens straight centroid edges.
        # This is a barrier screen, not a geodesic or anatomical ground truth.
        tree=cKDTree(xy);pairs=tree.query_pairs(CFG['proximity_maximum_um'],output_type='ndarray')
        distances=np.linalg.norm(xy[pairs[:,0]]-xy[pairs[:,1]],axis=1)
        support=np.empty(len(pairs),np.float32)
        steps=np.linspace(0,1,CFG['tissue_segment_samples'])
        for start in range(0,len(pairs),50000):
            p=pairs[start:start+50000]
            samples=xy[p[:,0],None,:]*(1-steps[None,:,None])+xy[p[:,1],None,:]*steps[None,:,None]
            x=np.clip(np.rint(samples[:,:,0]/res).astype(int),0,labels.shape[1]-1)
            y=np.clip(np.rint(samples[:,:,1]/res).astype(int),0,labels.shape[0]-1)
            support[start:start+len(p)]=(labels[y,x]>0).mean(1)
        for radius in [30,50,100]:
            m=(distances<=radius)&(support>=CFG['tissue_minimum_fraction'])
            p=pairs[m]
            g=sparse.csr_matrix((np.ones(2*len(p)),(np.r_[p[:,0],p[:,1]],np.r_[p[:,1],p[:,0]])),shape=(n,n))
            sparse.save_npz(OUT/stage/f'proximity{radius}_graph.npz',g)
            totals.append({'stage':stage,'graph':f'proximity{radius}','undirected_edges':len(p),'zero_degree':int((np.asarray(g.sum(1)).ravel()==0).sum()),
                'distance_only_candidate_edges':int((distances<=radius).sum()),'barrier_rejected_edges':int(((distances<=radius)&(support<CFG['tissue_minimum_fraction'])).sum())})
        totals.append({'stage':stage,'graph':'contact','undirected_edges':int(use.sum()),'zero_degree':int((np.asarray(a.sum(1)).ravel()==0).sum()),
            'contact_minimum_um':CFG['contact_minimum_um']})
        rng=np.random.default_rng(CFG['seed']);sample=rng.choice(len(pairs),min(10000,len(pairs)),replace=False)
        save(pd.DataFrame({'profile_i':pairs[sample,0],'profile_j':pairs[sample,1],
            'centroid_distance_um':distances[sample],'tissue_support':support[sample]}),f'{stage}/proximity_edge_audit_sample.tsv.gz')
        del labels,pairs,a;gc.collect()
    save(pd.DataFrame(totals),'graph_summary.tsv')

def neighborhood(stage,g):
    f=frame(stage);m=graph(stage,g);n=len(f)
    types=pd.Categorical(f.atlas_identity,categories=TYPES).codes
    one=sparse.csr_matrix((np.ones(n),(np.arange(n),types)),shape=(n,len(TYPES)))
    num=(m@one).toarray();resolved=num[:,:len(LINEAGES)].sum(1);alln=num.sum(1)
    d=pd.DataFrame(index=f.index);d['neighbors']=alln;d['resolved_neighbors']=resolved
    d['review_fraction']=num[:,len(LINEAGES):].sum(1)/np.maximum(alln,1)
    for j,t in enumerate(LINEAGES):d['fraction__'+t]=num[:,j]/np.maximum(resolved,1)
    for l in ['Endothelial','Fibroblast']:
        source=f.atlas_identity.eq(l)&f.program_source_eligible
        count=np.asarray(m@source.to_numpy(float)).ravel()
        d[l+'_eligible_neighbors']=count
        raw=num[:,LINEAGES.index(l)];d[l+'_source_coverage']=count/np.maximum(raw,1)
        d[l+'_source_capture']=np.asarray(m@np.where(source,f.source_log_capture,0)).ravel()/np.maximum(count,1)
    for name,cfg in CFG['programs'].items():
        valid=f[name].notna().to_numpy();count=np.asarray(m@valid.astype(float)).ravel()
        total=np.asarray(m@f[name].fillna(0).to_numpy()).ravel()
        if cfg['lineage'] in ['Endothelial','Fibroblast']:
            d['program__'+name]=np.divide(total,count,out=np.full(n,np.nan),where=count>0)
        else:
            # Burden of measured source RNA per resolved profile; absence of an
            # observed source differs from lack of a measured state, audited below.
            d['program__'+name]=total/np.maximum(resolved,1)
            d['source_present__'+name]=(count>0).astype(float)
    for prefix,families,col,lineage in [('FB',FBF,'fine_identity','Fibroblast'),('EC',ECF,'fine_identity','Endothelial')]:
        total=num[:,LINEAGES.index(lineage)]
        for fam in families:
            indicator=f.atlas_identity.eq(lineage)&f[col].eq(fam)
            d[prefix+'family__'+fam]=np.asarray(m@indicator.to_numpy(float)).ravel()/np.maximum(total,1)
    return d

def neighborhoods():
    coverage=[]
    for stage in STAGES:
        print('NEIGHBORHOODS',stage,flush=True);f=frame(stage)
        cms=f.loc[f.ventricular_CM_strict].copy()
        masks=[];tables={}
        for name in ['contact','proximity30','proximity50','proximity100']:
            d=neighborhood(stage,name).loc[cms.index]
            valid=(d.resolved_neighbors>=CFG['model_minimum_neighbors'])&(d.Endothelial_eligible_neighbors>=1)&(d.Fibroblast_eligible_neighbors>=1)
            valid&=(cms.CM_outcome_excluded_depth>=CFG['model_minimum_core_transcripts'])&cms.arterial_distance_um.notna()&cms.endocardial_proxy_distance_um.notna()
            d['geometry_eligible']=valid;masks.append(valid) if name in ['contact','proximity50'] else None
            tables[name]=d
            coverage.append({'stage':stage,'geometry':name,'strict_ventricular_CMs':len(cms),'eligible_CMs':int(valid.sum())})
        common=masks[0]&masks[1]
        for name,d in tables.items():
            d['common_contact_50um_cohort']=common
            save(d,f'{stage}/{name}_CM_neighborhoods.tsv.gz',index=True)
        coverage.append({'stage':stage,'geometry':'common_contact_50um','strict_ventricular_CMs':len(cms),'eligible_CMs':int(common.sum())})
    save(pd.DataFrame(coverage),'CM_cohort_coverage.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['atlas','geometry','neighborhoods','all'],default='all');args=p.parse_args()
    OUT.mkdir(exist_ok=True);dump('configuration.json',CFG)
    for name,fn in [('atlas',atlas),('geometry',geometry),('neighborhoods',neighborhoods)]:
        if args.step in [name,'all']:fn()

if __name__=='__main__':main()
