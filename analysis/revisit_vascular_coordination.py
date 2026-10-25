"""Literature-guided measured vascular states, graph rings and contact hypotheses.

Lu2021 separates coronary plexus from endocardium and shows DLL4/NOTCH
negative feedback. Wang2026 motivates matched-control maturity, graph rings,
CM–EC PTPRM and CM–CM comparison. No unmeasured sprouting RNA is fabricated.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-dynamics-mpl')
from pathlib import Path
import json,argparse,gc,hashlib
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree,Delaunay
from scipy.sparse.csgraph import dijkstra
from sklearn.preprocessing import SplineTransformer,StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
from cardiac_spatial_pipeline import ROOT,LINEAGES,MARKERS,MAT
from cardiac_spatial_models import joint_fit

CONFIG=ROOT/'analysis/cardiac_niche_dynamics.json'
CFG=json.loads(CONFIG.read_text());OUT=ROOT/CFG['output'];UP=ROOT/CFG['upstream'];STAGES=CFG['stages']
EXTRA='Ptprm Nampt Insr Dll4 Notch1 Flt4 Prox1 Npr3 Nrxn3 Aplnr Emcn Hey1 Hey2 Kdr Nrp1 Angpt2 Mki67 Top2a Cenpf'.split()
EXTRA=list(dict.fromkeys(EXTRA+sum(CFG['vascular_programs'].values(),[])))

def read(p,**kw):kw.setdefault('index_col',0);return pd.read_csv(p,sep='\t',**kw)
def save(d,p,index=False):
    p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,sep='\t',index=index,na_rep='NA',float_format='%.9g')
def dump(p,d):
    p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,default=str)+'\n')
def logs(x):return np.log1p(x.toarray()*1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1)[:,None])
def sparse_logs(x):return x.multiply((1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1))[:,None]).log1p().tocsr()
def mean_near(g,values,source):
    count=np.asarray(g@np.asarray(source,float)).ravel()
    total=np.asarray(g@np.where(source,np.asarray(values,float),0)).ravel()
    return np.divide(total,count,out=np.full(len(count),np.nan),where=count>0),count

def native_boundary_graph(stage,f):
    """Native polygon distance is a count/geometry sensitivity, not truth."""
    from shapely import Polygon,STRtree,distance
    rawbase=Path('/workspace/heart-atlas/data/xenium')
    raw=next(rawbase.glob('WT_'+('P07' if stage=='P7' else stage)+'__*'))
    p=raw/'cell_boundaries.csv.gz';v=pd.read_csv(p,usecols=['cell_id','vertex_x','vertex_y'])
    polys=[None]*len(f);lookup={x:i for i,x in enumerate(f.index)}
    for cid,z in v.groupby('cell_id',sort=False):
        if cid in lookup and len(z)>=3:polys[lookup[cid]]=Polygon(z[['vertex_x','vertex_y']].to_numpy())
    valid=np.array([p is not None and not p.is_empty for p in polys]);ix=np.flatnonzero(valid)
    shapes=np.array([polys[i] for i in ix],dtype=object)
    pairs=STRtree(shapes).query(shapes,predicate='dwithin',distance=max(CFG['native_boundary_distances_um']))
    take=pairs[0]<pairs[1];a,b=pairs[:,take];dist=distance(shapes[a],shapes[b]);a,b=ix[a],ix[b]
    for threshold in CFG['native_boundary_distances_um']:
        keep=dist<=threshold
        g=sparse.csr_matrix((np.ones(keep.sum()*2),(np.r_[a[keep],b[keep]],np.r_[b[keep],a[keep]])),shape=(len(f),len(f)))
        sparse.save_npz(OUT/stage/f'native_boundary_{threshold:g}um_graph.npz',g)
    save(pd.DataFrame({'profile_i':a,'profile_j':b,'native_polygon_distance_um':dist}),f'{stage}/native_boundary_edges.tsv.gz')
    return str(p)

def prepare():
    summaries=[];ringaudit=[];cellmeans=[];controlcal=[];sources=[];coverage=[]
    for stage in STAGES:
        print('VASCULAR_PREPARE',stage,flush=True)
        folder=OUT/stage;folder.mkdir(parents=True,exist_ok=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz')
        a=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad')
        genes=a.var_names.to_numpy(dtype=str);core=a.layers['high_confidence'].tocsr();nuc=a.layers['nuclear_only'].tocsr()
        p=ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad'
        native=ad.read_h5ad(p,backed='r');assert np.array_equal(native.obs_names,f.index)
        order=native.var_names.get_indexer(genes);nat=native.X[:,order].tocsr()
        sources.extend([str(p.relative_to(ROOT)),str((ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad').relative_to(ROOT))])
        ann=read(ROOT/f'results/endothelial_QC_subtypes/Xenium/{stage}/EC_consensus_annotations.tsv.gz')
        e=f.atlas_identity.eq('Endothelial')
        ef=f.index.get_indexer(ann.index.intersection(f.index));assert (ef>=0).all()
        q=f[['stage_cell_id','x_centroid','y_centroid','atlas_identity','program_source_eligible','EC_nuclear_type','author_atlas_cell_state']].copy()
        q['author_state_normalized']=q.author_atlas_cell_state.fillna('').str.rstrip('.')
        q['EC_primary']=e&f.program_source_eligible
        # All broad ECs remain in the audit, including mixed myocardial programs.
        q['EC_broad']=e
        measured_extra=[g for g in EXTRA if g in genes]
        positions=pd.Index(genes).get_indexer(measured_extra)
        for name,panel in CFG['vascular_programs'].items():
            coverage.append({'stage':stage,'program':name,'requested_genes':';'.join(panel),
                'measured_genes':';'.join(g for g in panel if g in genes),'unmeasured_genes':';'.join(g for g in panel if g not in genes) or 'None'})
        for view,x in [('native',nat),('core',core),('nuclear',nuc)]:
            depth=np.asarray(x.sum(1)).ravel();values=x[:,positions].toarray()
            norm=np.log1p(values*1e4/np.maximum(depth,1)[:,None])
            for j,g in enumerate(measured_extra):q[view+'__'+g]=norm[:,j];q[view+'_molecules__'+g]=values[:,j]
            for name,panel in CFG['vascular_programs'].items():q[view+'__'+name]=norm[:,[measured_extra.index(g) for g in panel if g in measured_extra]].mean(1)
        q['EC_Flt4_Dll4_core_dual']=q.EC_primary&q.core_molecules__Flt4.gt(0)&q.core_molecules__Dll4.gt(0)&~q.EC_nuclear_type.eq('Lymphatic-like')
        q['EC_growth_compatible']=q.EC_primary&~q.EC_nuclear_type.isin(['Endocardial-like','Lymphatic-like'])&q.core_molecules__Flt4.gt(0)&q.core_molecules__Dll4.gt(0)
        q['EC_Notch_response_supported']=q.EC_primary&((q.core_molecules__Hey1>0)|(q.core_molecules__Hey2>0))
        q['author_myofiber_EC']=q.author_state_normalized.eq('cap.EC.High.Myo')
        q['EC_primary_source_exclusion']=e&~q.EC_primary
        contacts=sparse.load_npz(UP/stage/'contact_graph.npz')
        dll4_neighbors,ec_ec_n=mean_near(contacts,q.nuclear__Dll4,q.EC_primary)
        q['observed_EC_EC_neighbors']=ec_ec_n
        q['EC_EC_DLL4_NOTCH1_compatibility']=np.sqrt(np.maximum(dll4_neighbors,0)*q.nuclear__Notch1)
        q.loc[~q.EC_primary|q.observed_EC_EC_neighbors.eq(0),'EC_EC_DLL4_NOTCH1_compatibility']=np.nan
        save(q,f'{stage}/measured_vascular_states.tsv.gz',index=True)
        # All-cell and CM-only expression calibration are distinct score choices.
        lc=sparse_logs(core);ln=sparse_logs(nat)
        cellmeans.append((genes,np.asarray(lc.mean(0)).ravel(),np.asarray(ln.mean(0)).ravel()))
        cm=f.loc[f.ventricular_CM_strict].copy();ix=f.index.get_indexer(cm.index)
        matpos=np.isin(genes,MAT.gene)
        cm['maturity_native_21']=np.asarray(ln[ix][:,matpos].mean(1)).ravel()
        cm['maturity_core_21']=np.asarray(lc[ix][:,matpos].mean(1)).ravel()
        cm['CM_Ptprm_core']=q.loc[cm.index,'core__Ptprm']
        cm['CM_Insr_core']=q.loc[cm.index,'core__Insr']
        save(cm,f'{stage}/CM_outcomes_and_context.tsv.gz',index=True)
        xy=f[['x_centroid','y_centroid']].to_numpy()
        # Low-degree graph benchmark closer to the paper's neighborhood rings.
        tri=Delaunay(xy);simp=tri.simplices
        pairs=np.sort(np.vstack([simp[:,[0,1]],simp[:,[0,2]],simp[:,[1,2]]]),axis=1);pairs=np.unique(pairs,axis=0)
        maxgraph=sparse.load_npz(UP/stage/'proximity100_graph.npz').tocsr()
        ok=np.asarray(maxgraph[pairs[:,0],pairs[:,1]]).ravel()>0;pairs=pairs[ok]
        gd=sparse.csr_matrix((np.ones(2*len(pairs)),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(f),len(f)))
        sparse.save_npz(folder/'pruned_Delaunay_graph.npz',gd)
        native_source=native_boundary_graph(stage,f);sources.append(native_source)
        subtype={
            'All_broad_EC':q.EC_broad,
            'Eligible_EC':q.EC_primary,
            'Capillary_compatible':q.EC_primary&q.EC_nuclear_type.eq('Capillary-like'),
            'Endocardial_compatible':q.EC_primary&q.EC_nuclear_type.eq('Endocardial-like'),
            'Arterial_compatible':q.EC_primary&q.EC_nuclear_type.eq('Arterial-like'),
            'Growth_compatible':q.EC_growth_compatible,
            'Author_capillary':q.author_state_normalized.isin(['cap.EC','cap.EC.High.Myo']),
            'Author_proliferating_EC':q.author_atlas_cell_state.eq('p.EC')}
        ring=cm[['x_centroid','y_centroid','maturity_native_21','maturity_core_21']].copy()
        for name,take in subtype.items():
            si=np.flatnonzero(take.to_numpy())
            ring['nearest_'+name+'_um']=cKDTree(xy[si]).query(xy[ix])[0] if len(si) else np.nan
            for graphname,g in [('Delaunay',gd),('refined_contact',sparse.load_npz(UP/stage/'contact_graph.npz'))]:
                distances=dijkstra(g,directed=False,indices=si,unweighted=True,min_only=True,limit=CFG['rings_maximum']) if len(si) else np.full(len(f),np.inf)
                vals=distances[ix];ring[graphname+'_rings__'+name]=np.where(np.isfinite(vals),vals,np.nan)
                ring[graphname+'_outside10__'+name]=~np.isfinite(vals)
                ringaudit.append({'stage':stage,'graph':graphname,'EC_class':name,'seed_ECs':len(si),'CMs':len(cm),'CMs_reached_within10':int(np.isfinite(vals).sum())})
        save(ring,f'{stage}/CM_vascular_rings.tsv.gz',index=True)
        # Sender/receiver components stay separate; products alone can recycle
        # a maturation-associated receptor into a supposed non-CM effect.
        env=pd.DataFrame(index=f.index)
        for geom,g in [('proximity50',sparse.load_npz(UP/stage/'proximity50_graph.npz')),('refined_contact',sparse.load_npz(UP/stage/'contact_graph.npz')),
                       ('native_boundary',sparse.load_npz(folder/'native_boundary_1um_graph.npz'))]:
            for name in ['Plexus_growth_proxy','Notch_response_proxy','Outer_venous_compatible','Cycling']:
                vals,n=mean_near(g,q['nuclear__'+name],q.EC_primary)
                env[geom+'__'+name]=vals;env[geom+'__eligible_EC_neighbors']=n
            for name,mask in [('Capillary',subtype['Capillary_compatible']),('Endocardial',subtype['Endocardial_compatible']),('Growth',subtype['Growth_compatible'])]:
                env[geom+'__'+name+'_fraction']=np.asarray(g@mask.to_numpy(float)).ravel()/np.maximum(np.asarray(g@q.EC_primary.to_numpy(float)).ravel(),1)
            for gene in ['Ptprm','Nampt']:
                vals,n=mean_near(g,q['nuclear__'+gene],q.EC_primary);env[geom+'__EC_sender_'+gene]=vals
            vals,n=mean_near(g,q.EC_EC_DLL4_NOTCH1_compatibility.fillna(0),q.EC_EC_DLL4_NOTCH1_compatibility.notna())
            env[geom+'__EC_EC_Notch_contact_compatibility']=vals;env[geom+'__EC_EC_network_informative_neighbors']=n
            vals,n=mean_near(g,q.core__Ptprm,f.atlas_identity.eq('Cardiomyocyte'))
            env[geom+'__CM_sender_Ptprm']=vals;env[geom+'__CM_neighbors']=n
            env[geom+'__CM_receiver_Ptprm']=q.core__Ptprm
            env[geom+'__Ptprm_EC_CM_geomean']=np.sqrt(np.maximum(env[geom+'__EC_sender_Ptprm'],0)*q.core__Ptprm)
            for name,l in [('EC','Endothelial'),('FB','Fibroblast'),('SMC','Smooth muscle'),('Immune','Immune')]:
                n=np.asarray(g@f.atlas_identity.eq(l).to_numpy(float)).ravel()
                den=np.asarray(g@f.atlas_resolved.to_numpy(float)).ravel()
                env[geom+'__'+name+'_fraction_all_resolved']=n/np.maximum(den,1)
            env[geom+'__resolved_neighbors']=np.asarray(g@f.atlas_resolved.to_numpy(float)).ravel()
        save(env.loc[cm.index],f'{stage}/CM_vascular_exposures.tsv.gz',index=True)
        for key,mask in subtype.items():
            summaries.append({'stage':stage,'EC_definition':key,'EC_profiles':int(mask.sum()),'median_CM_distance_um':ring['nearest_'+key+'_um'].median()})
        sources.extend([str((UP/stage/'atlas_annotations.tsv.gz').relative_to(ROOT)),str((ROOT/f'results/endothelial_QC_subtypes/Xenium/{stage}/EC_consensus_annotations.tsv.gz').relative_to(ROOT))])
        native.file.close();del a,core,nuc,nat,lc,ln,q,env,f;gc.collect()
    save(pd.DataFrame(summaries),'vascular_stage_summary.tsv');save(pd.DataFrame(ringaudit),'ring_graph_coverage.tsv')
    save(pd.DataFrame(coverage),'vascular_program_panel_coverage.tsv')
    # Audit literal AddModuleScore control feasibility on the measured403 targets.
    genes=cellmeans[0][0];mean=np.mean([x[1] for x in cellmeans],axis=0)
    ranks=np.argsort(mean,kind='stable');bins=np.empty(len(mean),int)
    for i,part in enumerate(np.array_split(ranks,24)):bins[part]=i
    for gene in MAT.gene:
        gi=list(genes).index(gene);size=int((bins==bins[gi]).sum())
        controlcal.append({'gene':gene,'measured_pool_genes':len(genes),'default24bin_pool_size':size,'reported100_without_replacement_feasible':size>=100,
                           'scope':'Methods leave control pool/binning ambiguous; this is not evidence the authors implementation failed'})
    save(pd.DataFrame(controlcal),'reported_maturation_control_feasibility.tsv')
    # Feasible explicitly named all-cell matched-control sensitivity, not exact
    # reproduction of the insufficiently specified author100-control procedure.
    order=np.argsort(mean,kind='stable');bins=np.empty(len(mean),int)
    for i,part in enumerate(np.array_split(order,CFG['matched_control_bins'])):bins[part]=i
    rng=np.random.default_rng(CFG['seed']);controls=[];mapping=[]
    for gene in MAT.gene:
        gi=list(genes).index(gene);pool=np.flatnonzero((bins==bins[gi])&~np.isin(genes,MAT.gene))
        ix=rng.choice(pool,min(CFG['matched_controls_per_gene'],len(pool)),replace=False);controls.extend(ix)
        mapping.extend({'target_gene':gene,'control_gene':genes[j],'expression_bin':int(bins[j])} for j in ix)
    controls=np.unique(controls);save(pd.DataFrame(mapping),'feasible_allcell_matched_controls.tsv')
    for stage in STAGES:
        c=read(OUT/stage/'CM_outcomes_and_context.tsv.gz')
        a=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        idx=a.obs_names.get_indexer(c.index);core=a.layers['high_confidence'][idx];l=sparse_logs(core)
        c['maturity_allcell_matched_control']=c.maturity_core_21-np.asarray(l[:,controls].mean(1)).ravel()
        c['maturity_CM_fixed_control']=read(ROOT/f'results/cardiac_spatial_maturation_EC/{stage}/ventricular_CM_maturation.tsv.gz').maturation_core_matched_control
        save(c,f'{stage}/CM_outcomes_and_context.tsv.gz',index=True);a.file.close()
    dump('source_files.json',sorted(set(sources)));dump('configuration.json',CFG)

def models():
    coefficients=[];summaries=[]
    for stage in STAGES:
        print('VASCULAR_MODELS',stage,flush=True)
        c=read(OUT/stage/'CM_outcomes_and_context.tsv.gz');env=read(OUT/stage/'CM_vascular_exposures.tsv.gz')
        ring=read(OUT/stage/'CM_vascular_rings.tsv.gz')
        for geom in ['proximity50','refined_contact','native_boundary']:
            cols=[geom+'__'+p for p in ['Plexus_growth_proxy','Notch_response_proxy','Cycling','Capillary_fraction','Endocardial_fraction','EC_sender_Ptprm','CM_sender_Ptprm']]
            selected=(env[geom+'__eligible_EC_neighbors']>=1)&(env[geom+'__resolved_neighbors']>=3)&c.CM_outcome_excluded_depth.ge(CFG['minimum_CM_core_transcripts'])
            selected &= c.arterial_distance_um.notna()&c.endocardial_proxy_distance_um.notna()&env[cols].notna().all(axis=1)
            d=c.loc[selected];e=env.loc[selected,cols].copy();e.columns=[s.split('__',1)[1] for s in e]
            xx=d[['x_centroid','y_centroid']].to_numpy();outcomes=['maturity_core_21','maturity_native_21','maturity_allcell_matched_control','maturity_CM_fixed_control','CM_Sarcomere','CM_Metabolism']
            # Existing atlas column names use CM_ prefix; imported source includes
            # these direct sub-scores, with their original study labels.
            for level in ['technical','anatomy','CM_state','receiver_conditioned']:
                b=pd.DataFrame(index=d.index)
                for key in ['cell_area','CM_outcome_excluded_depth']:b['log_'+key]=np.log1p(d[key])
                b['nucleus_missing']=d.nucleus_area.isna().astype(float)
                b['log_neighbors']=np.log1p(env.loc[selected,geom+'__resolved_neighbors'])
                for l in ['FB','SMC','Immune']:b[l+'_fraction']=env.loc[selected,geom+'__'+l+'_fraction_all_resolved']
                if level!='technical':
                    for key in ['arterial_distance_um','endocardial_proxy_distance_um']:b['log_'+key]=np.log1p(d[key])
                    sp=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(xx)
                    for j in range(sp.shape[1]):b['xy_'+str(j)]=sp[:,j]
                if level in ['CM_state','receiver_conditioned']:b['CM_cycling']=d.CM_cycling;b['CM_Slit2']=d.CM_Slit2
                if level=='receiver_conditioned':b['CM_receiver_Ptprm']=d.CM_Ptprm_core
                if len(d)<200:continue
                ee=e[[key for key in e if e[key].std()>1e-8]]
                for r in joint_fit(d[outcomes],b,ee,xx,CFG['bootstrap_draws']):
                    coefficients.append({'stage':stage,'geometry':geom,'adjustment':level,**r})
            summaries.append({'stage':stage,'geometry':geom,'CMs':len(d)})
        for source in ['Capillary_compatible','Endocardial_compatible','Growth_compatible','All_broad_EC','Author_capillary','Author_proliferating_EC']:
            for graphname in ['Delaunay','refined_contact']:
                vals=ring[graphname+'_rings__'+source]
                bins=np.select([vals.eq(1),vals.eq(2),vals.eq(3),vals.eq(4),vals.between(5,6),vals.between(7,10)],['1','2','3','4','5-6','7-10'],default='>10/unreachable')
                for key in np.unique(bins):
                    take=bins==key
                    for outcome in ['maturity_core_21','maturity_native_21','maturity_allcell_matched_control']:
                        summaries.append({'stage':stage,'graph':graphname,'EC_definition':source,'ring':key,'outcome':outcome,'CMs':int(take.sum()),
                            'mean_maturity':c.loc[take,outcome].mean(),'median_maturity':c.loc[take,outcome].median(),
                            'scope':'Descriptive graph-ring dose response; not an independent replicated or causal test'})
    save(pd.DataFrame(coefficients),'vascular_maturation_coefficients.tsv.gz');save(pd.DataFrame(summaries),'vascular_ring_and_model_summary.tsv')

def RNA_and_source_audit():
    from scipy.stats import spearmanr
    print('RNA_FULL_TRANSCRIPTOME_PROGRAMS',flush=True)
    rna_path=Path('/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad')
    a=ad.read_h5ad(rna_path,backed='r');symbols=a.var.gene_symbol.astype(str).to_numpy()
    programs={'Sprouting_6g':'Apln Esm1 Angpt2 Cxcr4 Pgf Adm'.split(),
        'Plexus_growth_2021':CFG['vascular_programs']['Plexus_growth_proxy'],
        'Notch_response_RNA':'Notch1 Hey1 Hey2 Hes1 Nrarp'.split(),
        'Cycling_RNA':'Mki67 Top2a Cenpf Ccnb1 Ube2c Cdk1'.split(),
        'Hypoxia_compatible_RNA':'Slc2a1 Pdk1 Ldha Bnip3 Vegfa Pgk1 Hk2 Adm Egln3'.split(),
        'CM_maturation_21':MAT.gene.tolist()}
    genes=sorted(set(sum(programs.values(),[])+['Ptprm','Npr3','Aplnr','Emcn','Dll4','Flt4','Myh6','Tnni3','Actc1']))
    positions=np.flatnonzero(np.isin(symbols,genes));subset=a.X[:,positions].toarray()
    # Duplicate Ensembl-to-symbol mappings are summed, never silently selected.
    counts=np.zeros((len(a.obs),len(genes)),np.float32)
    for j,pos in enumerate(positions):counts[:,genes.index(symbols[pos])]+=subset[:,j]
    total=a.obs.total_counts.to_numpy(float);norm=np.log1p(counts*1e4/np.maximum(total,1)[:,None])
    primary=read(ROOT/'results/endothelial_QC_subtypes/snRNA/EC_annotations.tsv.gz')
    rows=[];gene_rows=[];coherence=[]
    for stage in STAGES:
        obs=a.obs;ec=obs.reviewed_identity.eq('Endothelial')&obs.age.eq(stage)
        primary_ec=ec&obs.index.isin(primary.index)
        cm=obs.reviewed_identity.eq('Cardiomyocyte')&obs.proposed_subtype.isin(['Ventricular-like','Conduction-associated CM candidate'])&obs.age.eq(stage)
        for scope,take in [('EC_primary',primary_ec),('EC_all_global_QC',ec),('Ventricular_CM',cm)]:
            for lib in sorted(obs.loc[take,'sample'].unique()):
                mask=take&obs['sample'].eq(lib);idx=np.flatnonzero(mask.to_numpy());n=len(idx)
                pooled=counts[idx].sum(0);cpm=pooled/max(total[idx].sum(),1)*1e6
                for name,panel in programs.items():
                    gi=[genes.index(g) for g in panel]
                    rows.append({'stage':stage,'library':lib,'scope':scope,'program':name,'nuclei':n,
                        'mean_per_nucleus_log1p_fullRNA_CP10k':norm[idx][:,gi].mean(),
                        'pseudobulk_mean_log1p_fullRNA_CPM':np.log1p(cpm[gi]).mean(),
                        'genes':';'.join(panel),'scope_note':'Corrected whole-transcriptome counts; no age-correlation p-values from four cross-sectional stages'})
                for j,g in enumerate(genes):gene_rows.append({'stage':stage,'library':lib,'scope':scope,'gene':g,'pseudobulk_fullRNA_CPM':cpm[j],
                    'detection_fraction_corrected_ge001':np.mean(counts[idx,j]>=.001)})
                if scope=='EC_primary':
                    fam=primary.loc[obs.index[idx],'EC_family'];dummy=pd.get_dummies(fam,dtype=float).to_numpy()
                    cov=np.column_stack([np.ones(n),np.log1p(total[idx]),dummy])
                    for name,panel in programs.items():
                        if name=='CM_maturation_21':continue
                        gi=[genes.index(g) for g in panel]
                        retained=[j for j in gi if (counts[idx,j]>=.001).sum()>=max(10,.02*n)]
                        cors=[]
                        for j in retained:
                            others=[k for k in retained if k!=j]
                            if not others:continue
                            y=norm[idx,j];z=norm[idx][:,others].mean(1)
                            yr=y-cov@np.linalg.lstsq(cov,y,rcond=1e-9)[0];zr=z-cov@np.linalg.lstsq(cov,z,rcond=1e-9)[0]
                            if np.std(yr)>1e-8 and np.std(zr)>1e-8:cors.append(spearmanr(yr,zr).statistic)
                        coherence.append({'stage':stage,'library':lib,'program':name,'eligible_genes':len(retained),'requested_genes':len(panel),
                            'median_leave_one_gene_out_partial_Spearman':np.median(cors) if cors else np.nan,
                            'interpretation':'Depth/family-adjusted descriptive coordination; sparse markers or regulatory feedback need not form one positive coexpression module'})
    save(pd.DataFrame(rows),'RNA_program_library_profiles.tsv');save(pd.DataFrame(gene_rows),'RNA_measured_and_unmeasured_gene_support.tsv')
    save(pd.DataFrame(coherence),'RNA_vascular_module_coherence.tsv');a.file.close()
    source_selection_audit()

def source_selection_audit():
    # Source selection can exclude a biologically proposed myofiber-EC state.
    audit=[]
    for stage in STAGES:
        q=read(OUT/stage/'measured_vascular_states.tsv.gz');f=read(UP/stage/'atlas_annotations.tsv.gz')
        c=read(ROOT/f'results/endothelial_QC_subtypes/Xenium/{stage}/EC_consensus_annotations.tsv.gz')
        for state in ['cap.EC.High.Myo','cap.EC','p.EC','endocardial.EC','arteriole.EC','lym.EC']:
            ids=q.index[q.author_atlas_cell_state.fillna('').str.rstrip('.').eq(state)];z=c.reindex(ids)
            audit.append({'stage':stage,'author_state':state,'published_profiles':len(ids),'our_broad_EC':int(q.loc[ids,'EC_broad'].sum()),
                'our_primary_EC_source':int(q.loc[ids,'EC_primary'].sum()),
                'EC_candidate_evidence_available':int(z.nuclear_broad_EC_clean.notna().sum()),
                'nuclear_clean_EC_support':int(z.nuclear_broad_EC_clean.fillna(False).sum()),
                'nuclear_strong_CM_program':int(z.nuclear_only__strong_Cardiomyocyte_program.fillna(False).sum()),
                'scope':'A dual myocardial/endothelial program is a review question, not automatic proof of a doublet or a verified new EC type'})
    save(pd.DataFrame(audit),'author_EC_state_source_selection_audit.tsv')

def refresh_author_labels():
    """Canonicalize a trailing punctuation difference in the public label."""
    summary=read(OUT/'vascular_stage_summary.tsv',index_col=None)
    graphcoverage=read(OUT/'ring_graph_coverage.tsv',index_col=None)
    for stage in STAGES:
        q=read(OUT/stage/'measured_vascular_states.tsv.gz');f=read(UP/stage/'atlas_annotations.tsv.gz')
        q['author_state_normalized']=q.author_atlas_cell_state.fillna('').str.rstrip('.')
        q['author_myofiber_EC']=q.author_state_normalized.eq('cap.EC.High.Myo')
        save(q,f'{stage}/measured_vascular_states.tsv.gz',index=True)
        ring=read(OUT/stage/'CM_vascular_rings.tsv.gz');xy=f[['x_centroid','y_centroid']].to_numpy();ix=f.index.get_indexer(ring.index)
        take=q.author_state_normalized.isin(['cap.EC','cap.EC.High.Myo']);seeds=np.flatnonzero(take.to_numpy())
        ring['nearest_Author_capillary_um']=cKDTree(xy[seeds]).query(xy[ix])[0]
        for name,g in [('Delaunay',sparse.load_npz(OUT/stage/'pruned_Delaunay_graph.npz')),('refined_contact',sparse.load_npz(UP/stage/'contact_graph.npz'))]:
            vals=dijkstra(g,directed=False,indices=seeds,unweighted=True,min_only=True,limit=CFG['rings_maximum'])[ix]
            ring[name+'_rings__Author_capillary']=np.where(np.isfinite(vals),vals,np.nan)
            ring[name+'_outside10__Author_capillary']=~np.isfinite(vals)
            m=graphcoverage.stage.eq(stage)&graphcoverage.graph.eq(name)&graphcoverage.EC_class.eq('Author_capillary')
            graphcoverage.loc[m,'seed_ECs']=len(seeds);graphcoverage.loc[m,'CMs_reached_within10']=int(np.isfinite(vals).sum())
        save(ring,f'{stage}/CM_vascular_rings.tsv.gz',index=True)
        m=summary.stage.eq(stage)&summary.EC_definition.eq('Author_capillary');summary.loc[m,'EC_profiles']=len(seeds)
        summary.loc[m,'median_CM_distance_um']=ring.nearest_Author_capillary_um.median()
    save(summary,'vascular_stage_summary.tsv');save(graphcoverage,'ring_graph_coverage.tsv');source_selection_audit()
    previous=read(OUT/'vascular_ring_and_model_summary.tsv',index_col=None)
    rows=previous.loc[previous.graph.isna()].to_dict('records')
    for stage in STAGES:
        c=read(OUT/stage/'CM_outcomes_and_context.tsv.gz');ring=read(OUT/stage/'CM_vascular_rings.tsv.gz')
        for source in ['Capillary_compatible','Endocardial_compatible','Growth_compatible','All_broad_EC','Author_capillary','Author_proliferating_EC']:
            for name in ['Delaunay','refined_contact']:
                vals=ring[name+'_rings__'+source]
                bins=np.select([vals.eq(1),vals.eq(2),vals.eq(3),vals.eq(4),vals.between(5,6),vals.between(7,10)],['1','2','3','4','5-6','7-10'],default='>10/unreachable')
                for key in np.unique(bins):
                    take=bins==key
                    for outcome in ['maturity_core_21','maturity_native_21','maturity_allcell_matched_control']:
                        rows.append({'stage':stage,'graph':name,'EC_definition':source,'ring':key,'outcome':outcome,'CMs':int(take.sum()),
                            'mean_maturity':c.loc[take,outcome].mean(),'median_maturity':c.loc[take,outcome].median(),
                            'scope':'Descriptive graph-ring dose response; not independent biological replicates'})
    save(pd.DataFrame(rows),'vascular_ring_and_model_summary.tsv')

def architecture_check():
    """Interpret the two coarse contexts and compare native/refined enrichment."""
    from scipy.stats import entropy
    rows=[];context=[];rng=np.random.default_rng(CFG['seed'])
    for stage in STAGES:
        print('ARCHITECTURE_CHECK',stage,flush=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz');niche=read(UP/stage/'niche_annotations.tsv.gz')
        names=niche.niche.unique()
        for name in names:
            ids=niche.index[niche.niche.eq(name)];z=f.loc[ids];g=sparse.load_npz(UP/stage/'contact_graph.npz')[f.index.get_indexer(ids)]
            count=np.asarray(g@f.atlas_resolved.to_numpy(float)).ravel()
            for lineage in LINEAGES:
                n=np.asarray(g@f.atlas_identity.eq(lineage).to_numpy(float)).ravel()
                context.append({'stage':stage,'computational_context':name,'focal_profiles':len(ids),'neighbor_type':lineage,
                    'mean_resolved_neighbor_fraction':np.mean(n/np.maximum(count,1)),
                    'meaning':'Repeated neighborhood composition, not a new cell type, two exhaustive biological niches, or proven signaling'})
        labels=f.atlas_identity.map({l:i for i,l in enumerate(LINEAGES)}).fillna(-1).to_numpy(int)
        resolved=f.atlas_resolved.to_numpy(bool);xy=f[['x_centroid','y_centroid']].to_numpy();blocks=np.floor(xy/200).astype(int)
        _,inv=np.unique(blocks,axis=0,return_inverse=True);groups=[np.flatnonzero((inv==i)&resolved) for i in np.unique(inv)]
        for geom,g in [('native_boundary_1um',sparse.load_npz(OUT/stage/'native_boundary_1um_graph.npz')),
                       ('refined_contact',sparse.load_npz(UP/stage/'contact_graph.npz'))]:
            coo=sparse.triu(g,1).tocoo();keep=resolved[coo.row]&resolved[coo.col];ii,jj=coo.row[keep],coo.col[keep]
            def tally(t):
                a=np.minimum(t[ii],t[jj]);b=np.maximum(t[ii],t[jj]);return np.bincount(a*7+b,minlength=49)
            observed=tally(labels);null=[]
            for _ in range(99):
                shuffled=labels.copy()
                for take in groups:shuffled[take]=rng.permutation(shuffled[take])
                null.append(tally(shuffled))
            null=np.array(null)
            for a in range(7):
                for b in range(a,7):
                    k=a*7+b;expected=null[:,k].mean()
                    rows.append({'stage':stage,'geometry':geom,'type_a':LINEAGES[a],'type_b':LINEAGES[b],
                        'observed_edges':observed[k],'local_shuffle_expected_edges':expected,
                        'log2_observed_over_local_null':np.log2((observed[k]+.5)/(expected+.5)),
                        'scope':'Within-section morphology/label sensitivity, not causal interaction or independent animal inference'})
    save(pd.DataFrame(rows),'native_refined_architecture_enrichment.tsv');save(pd.DataFrame(context),'coarse_context_composition_explanation.tsv')

def sensitivity():
    """Compare geometry on identical CMs and expose source-label dependence."""
    coefficients=[];coverage=[];source_coeff=[];ring_cor=[]
    outcomes=['maturity_core_21','maturity_native_21','maturity_allcell_matched_control']
    names=['Plexus_growth_proxy','Notch_response_proxy','Cycling','Capillary_fraction','Endocardial_fraction','EC_sender_Ptprm','CM_sender_Ptprm']
    def covariates(d,e,g):
        b=pd.DataFrame(index=d.index)
        for key in ['cell_area','CM_outcome_excluded_depth']:b['log_'+key]=np.log1p(d[key])
        b['nucleus_missing']=d.nucleus_area.isna().astype(float)
        for k in ['EC_source_log_nuclear_capture','EC_source_log_core_capture','EC_source_coverage']:b[k]=e[g+'__'+k]
        b['log_neighbors']=np.log1p(e[g+'__resolved_neighbors'])
        for l in ['FB','SMC','Immune']:b[l+'_fraction']=e[g+'__'+l+'_fraction_all_resolved']
        for key in ['arterial_distance_um','endocardial_proxy_distance_um']:b['log_'+key]=np.log1p(d[key])
        sp=SplineTransformer(n_knots=5,degree=3,include_bias=False).fit_transform(d[['x_centroid','y_centroid']])
        for j in range(sp.shape[1]):b['xy_'+str(j)]=sp[:,j]
        b['CM_cycling']=d.CM_cycling;b['CM_Slit2']=d.CM_Slit2;b['CM_receiver_Ptprm']=d.CM_Ptprm_core
        return b
    for stage in STAGES:
        print('PAIRED_AND_SOURCE_SENSITIVITY',stage,flush=True)
        c=read(OUT/stage/'CM_outcomes_and_context.tsv.gz');env=read(OUT/stage/'CM_vascular_exposures.tsv.gz')
        common=c.CM_outcome_excluded_depth.ge(CFG['minimum_CM_core_transcripts'])&c.arterial_distance_um.notna()&c.endocardial_proxy_distance_um.notna()
        for g in ['proximity50','refined_contact','native_boundary']:
            cols=[g+'__'+p for p in names]
            common &= env[g+'__eligible_EC_neighbors'].ge(1)&env[g+'__resolved_neighbors'].ge(3)&env[cols].notna().all(axis=1)
        d=c.loc[common];coverage.append({'stage':stage,'identical_CMs_all_three_geometries':len(d)})
        for g in ['proximity50','refined_contact','native_boundary']:
            e=env.loc[common,[g+'__'+p for p in names]].copy();e.columns=names;e=e.loc[:,e.std()>1e-8]
            for r in joint_fit(d[outcomes],covariates(d,env.loc[common],g),e,d[['x_centroid','y_centroid']].to_numpy(),CFG['bootstrap_draws']):
                coefficients.append({'stage':stage,'geometry':g,'adjustment':'receiver_conditioned','cohort':'Identical CMs for all three geometries',**r})
        # Ring profiles are descriptive; censored/unreachable distances never
        # receive an artificial numerical distance of 10 or 11.
        ring=read(OUT/stage/'CM_vascular_rings.tsv.gz')
        for key in [s for s in ring if '_rings__' in s]:
            for outcome in outcomes:
                take=ring[key].notna()&c[outcome].notna()
                x=ring.loc[take,key].to_numpy();y=c.loc[take,outcome].to_numpy()
                informative=len(x)>2 and np.std(x)>0
                ring_cor.append({'stage':stage,'distance_definition':key,'outcome':outcome,'CMs_reached_within10':int(take.sum()),
                    'unreachable_or_beyond10_CMs':int(ring[key].isna().sum()),'Pearson':np.corrcoef(x,y)[0,1] if informative else np.nan,
                    'Spearman':spearmanr(x,y).statistic if informative else np.nan,
                    'scope':'Descriptive reached-only association; no replicated or causal test'})
        if stage!='P7':continue
        q=read(OUT/stage/'measured_vascular_states.tsv.gz');f=read(UP/stage/'atlas_annotations.tsv.gz')
        g=sparse.load_npz(UP/stage/'proximity50_graph.npz');alt=pd.DataFrame(index=c.index)
        for label,source in [('primary_EC',q.EC_primary),('all_broad_EC',q.EC_broad),('author_capillary_including_HighMyo',q.author_state_normalized.isin(['cap.EC','cap.EC.High.Myo']))]:
            values,n=mean_near(g,q.nuclear__Ptprm,source)
            alt[label]=pd.Series(values,index=q.index).reindex(c.index);alt[label+'__neighbors']=pd.Series(n,index=q.index).reindex(c.index)
        cols=['proximity50__'+p for p in names]
        valid=c.CM_outcome_excluded_depth.ge(20)&c.arterial_distance_um.notna()&c.endocardial_proxy_distance_um.notna()&env[cols].notna().all(axis=1)&env.proximity50__resolved_neighbors.ge(3)
        valid &= alt[['primary_EC','all_broad_EC','author_capillary_including_HighMyo']].notna().all(axis=1)
        d=c.loc[valid];b=covariates(d,env.loc[valid],'proximity50')
        for label in ['primary_EC','all_broad_EC','author_capillary_including_HighMyo']:
            e=env.loc[valid,cols].copy();e.columns=names;e['EC_sender_Ptprm']=alt.loc[valid,label]
            for r in joint_fit(d[outcomes],b,e,d[['x_centroid','y_centroid']].to_numpy(),CFG['bootstrap_draws']):
                source_coeff.append({'stage':stage,'geometry':'proximity50','source_definition':label,
                    'scope':'Common CM cohort; sender PTPRM source changed, other primary EC exposures held fixed; public labels are not revalidated identities',**r})
        save(alt,'P7/Ptprm_source_definition_exposures.tsv.gz',index=True)
    save(pd.DataFrame(coefficients),'paired_geometry_coefficients.tsv.gz');save(pd.DataFrame(coverage),'paired_geometry_cohort.tsv')
    save(pd.DataFrame(source_coeff),'P7_Ptprm_source_sensitivity.tsv.gz');save(pd.DataFrame(ring_cor),'graph_ring_correlations.tsv')

def source_capture():
    """Expose source depth and selection coverage for geometry sensitivity."""
    for stage in STAGES:
        print('SOURCE_CAPTURE',stage,flush=True)
        q=read(OUT/stage/'measured_vascular_states.tsv.gz');e=read(OUT/stage/'CM_vascular_exposures.tsv.gz')
        a=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad',backed='r')
        for view,layer in [('core','high_confidence'),('nuclear','nuclear_only')]:
            depth=np.asarray(a.layers[layer][:].sum(1)).ravel();q[view+'__panel_depth']=depth
        for name,path in [('proximity50',UP/stage/'proximity50_graph.npz'),('refined_contact',UP/stage/'contact_graph.npz'),('native_boundary',OUT/stage/'native_boundary_1um_graph.npz')]:
            g=sparse.load_npz(path)
            for view in ['core','nuclear']:
                values,n=mean_near(g,np.log1p(q[view+'__panel_depth']),q.EC_primary)
                e[name+'__EC_source_log_'+view+'_capture']=pd.Series(values,index=q.index).reindex(e.index)
            den=np.asarray(g@q.EC_broad.to_numpy(float)).ravel()
            e[name+'__EC_source_coverage']=pd.Series(n/np.maximum(den,1),index=q.index).reindex(e.index)
        save(q,f'{stage}/measured_vascular_states.tsv.gz',index=True);save(e,f'{stage}/CM_vascular_exposures.tsv.gz',index=True);a.file.close()

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['prepare','models','rna','architecture','source_labels','source_capture','sensitivity','all'],default='all');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.step in ['prepare','all']:prepare()
    if args.step in ['models','all']:models()
    if args.step in ['rna','all']:RNA_and_source_audit()
    if args.step in ['architecture','all']:architecture_check()
    if args.step in ['source_capture','all']:source_capture()
    if args.step in ['sensitivity','all']:sensitivity()
    if args.step=='source_labels':refresh_author_labels()

if __name__=='__main__':main()
