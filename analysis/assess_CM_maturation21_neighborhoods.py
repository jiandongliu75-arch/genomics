"""Direct Xenium maturation scores, matched capillary tests and non-CM context.

Frozen cells, anatomy, source labels, graphs and pooled gene calibration are
retained. This is a transparent 21-gene standardized mean, not a claim to
recover the authors' unavailable precomputed Seurat score. No spatial RNA
imputation, new lineage tracing or animal-level replication is introduced.
"""
from pathlib import Path
import os, json, hashlib, time, argparse
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
import numpy as np
import pandas as pd
import h5py
from anndata.io import sparse_dataset, read_elem
from scipy import sparse, linalg
from scipy.stats import norm, spearmanr, pearsonr
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler, SplineTransformer
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
import cardiac_subtype_niches as prior

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/CM_maturation21_nonCM_context'
STAGES = ['P0', 'P7', 'P14', 'P21']
SEED = 271027
DRAWS = 512
DEPTH = 50
INPUTS = set()
PANEL = pd.read_csv(ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv', sep='\t').gene.tolist()
COMPONENTS = json.loads((ROOT/'analysis/cardiac_spatial_reconstruction.json').read_text())['CM_components']
REMOVALS = {'without_Cebpb_Rxrg_Lpin1': ['Cebpb', 'Rxrg', 'Lpin1'],
            'without_Myl2': ['Myl2'],
            'without_metabolic8': COMPONENTS['metabolic'],
            'without_sarcomere3': COMPONENTS['sarcomere'],
            'without_excitation5': COMPONENTS['excitation_contraction']}
PRIMARY = ['core_Z__maturation21', 'core_Z__metabolic', 'core_Z__sarcomere',
           'core_Z__excitation_contraction', 'core_Z__cycling']
SENS = ['native_Z__'+k for k in ['maturation21', 'metabolic', 'sarcomere', 'excitation_contraction', 'cycling']]
SENS += ['core_Z__'+k for k in REMOVALS]
FAMILIES = ['Postn', 'Gfpt2', 'Col8a1', 'Lsamp', 'Kirrel3']
BROAD = ['Endothelial', 'Fibroblast', 'Pericyte', 'Smooth muscle', 'Immune', 'Epicardial']

def read(p, **kwargs):
    p = Path(p); p = p if p.is_absolute() else ROOT/p
    INPUTS.add(str(p)); return pd.read_csv(p, sep='\t', **kwargs)

def save(d, name, index=False):
    OUT.mkdir(exist_ok=True, parents=True)
    d.to_csv(OUT/name, sep='\t', index=index, na_rep='NA', float_format='%.9g')

def digest(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(4<<20), b''): h.update(b)
    return h.hexdigest()

def bh(p):
    p = np.asarray(p); order = np.argsort(p); result = np.empty(len(p))
    result[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return result

def score(counts, totals, genes, cal, view):
    log = np.log1p(counts*1e4/np.maximum(totals, 1)[:,None])
    ref = cal.loc[(view, 'all')].reindex(genes)
    z = (log-ref.reference_mean.to_numpy())/ref.reference_SD.to_numpy()
    return log, z

def rarefy_panel(counts, total, rng, depth=DEPTH):
    """Exact molecule subsampling without replacement; other genes are one bin."""
    assert np.allclose(counts, np.rint(counts)) and np.allclose(total, np.rint(total))
    counts = np.rint(counts).astype(np.int64); total = np.rint(total).astype(np.int64)
    available = total >= depth
    remain = np.where(available, depth, 0).astype(np.int64)
    population = total.copy(); selected = np.zeros_like(counts)
    for j in range(counts.shape[1]):
        good = counts[:,j]; bad = population-good
        take = (remain>0)&(population>0)
        selected[take,j] = rng.hypergeometric(good[take], bad[take], remain[take])
        remain -= selected[:,j]; population -= good
    assert (selected.sum(1)<=depth).all()
    assert (selected<=counts).all()
    return selected, available

def preparation():
    cal = read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').set_index(['view','scope','gene']).sort_index()
    gene_metrics=[]; omissions=[]; summary=[]; checks=[]; depth_checks=[]
    selected_genes = PANEL + COMPONENTS['cycling']
    assert len(set(PANEL))==21 and len(selected_genes)==24
    for stage in STAGES:
        print('MATURATION21_PREPARE', stage, flush=True)
        (OUT/stage).mkdir(exist_ok=True, parents=True)
        c = read(f'results/cardiac_spatial_reconstruction/{stage}/CM_outcomes_and_anatomy.tsv.gz', index_col=0)
        f,d = prior.stage_data(stage)
        keep = ['stage_cell_id','x_centroid','y_centroid','author_atlas_cell_state', 'reconstructed_Figure4e_CM',
                'core_panel_depth','native_panel_depth','morphology_region','distance_to_large_void_um','distance_to_outer_um','relative_wall_geometry']
        export = c[keep].copy(); export['stage']=stage; export['strict_ventricular_map_cohort']=export.index.isin(d.index)
        export['CM_state'] = export.author_atlas_cell_state
        export['core_depth50_eligible']=c.core_panel_depth.ge(DEPTH)
        for view in ['native','core']:
            p = (ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad' if view=='native'
                 else ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad')
            INPUTS.add(str(p))
            with h5py.File(p) as h:
                ids=read_elem(h['obs']).index; genes=read_elem(h['var']).index
                ii=ids.get_indexer(c.index); jj=genes.get_indexer(selected_genes)
                assert (ii>=0).all() and (jj>=0).all() and len(genes)==403
                group=h['X'] if view=='native' else h['layers/high_confidence']
                xx=sparse_dataset(group)[ii].tocsr()
                total=np.asarray(xx.sum(1)).ravel().astype(float)
                counts=xx[:,jj].toarray().astype(float)
            np.testing.assert_allclose(total, c[view+'_panel_depth'], rtol=1e-7)
            log,z=score(counts,total,selected_genes,cal,view)
            full=z[:,:21].mean(1); export[view+'_Z__maturation21']=full
            export[view+'_logmean21']=log[:,:21].mean(1)
            for key,panel in COMPONENTS.items():
                export[view+'_Z__'+key]=z[:,np.isin(selected_genes,panel)].mean(1)
            for key in ['maturation21', *COMPONENTS]:
                old=c[view+'_Z__'+key].to_numpy()
                error=np.max(abs(export[view+'_Z__'+key].to_numpy()-old))
                assert error<2e-6, (stage,view,key,error)
                checks.append({'stage':stage,'check':'recomputed_'+view+'_'+key,'profiles':len(c),'maximum_absolute_error':error,'status':'PASS'})
            for gene in PANEL:
                gi=PANEL.index(gene); reduced=(21*full-z[:,gi])/20
                for scope,mask in [('all_scored_CMs',np.ones(len(c),bool)),('strict_ventricular',c.index.isin(d.index))]:
                    omissions.append({'stage':stage,'view':view,'scope':scope,'omitted_gene':gene,'profiles':int(mask.sum()),
                        'Pearson_with_full':pearsonr(full[mask],reduced[mask]).statistic,
                        'Spearman_with_full':spearmanr(full[mask],reduced[mask]).statistic,
                        'normalization_and_reference_frozen':True})
                gene_metrics.append({'stage':stage,'view':view,'gene':gene,'scope':'strict_ventricular','profiles':len(d),
                    'mean_log1p_CP10k':log[c.index.isin(d.index),gi].mean(),
                    'mean_frozen_gene_Z':z[c.index.isin(d.index),gi].mean(),
                    'fraction_detected':(counts[c.index.isin(d.index),gi]>0).mean()})
            for name,removed in REMOVALS.items():
                export[view+'_Z__'+name]=z[:,:21][:,~np.isin(PANEL,removed)].mean(1)
            estimates=[]; rng=np.random.default_rng(SEED+int(stage[1:])+(0 if view=='core' else 1000))
            for _ in range(3):
                thinned,available=rarefy_panel(counts,total,rng)
                _,zt=score(thinned,np.where(available,DEPTH,1),selected_genes,cal,view)
                m=zt[:,:21].mean(1);m[~available]=np.nan;estimates.append(m)
            mean_thinned=np.mean(estimates,axis=0)
            export[view+'_Z__depth50_mean3']=mean_thinned
            mask=c.index.isin(d.index)
            depth_checks.append({'stage':stage,'view':view,'strict_CMs':int(mask.sum()),'depth50_CMs':int((mask&available).sum()),
                'score_depth_Spearman':spearmanr(full[mask],total[mask]).statistic,
                'full_depth50_Spearman_same_CMs':spearmanr(full[mask&available],mean_thinned[mask&available]).statistic,
                'fixed_depth_molecules':DEPTH,'subsampling_repetitions':3,'subsampling_without_replacement':True})
            for scope,take in [('all_scored_CMs',np.ones(len(c),bool)),('strict_ventricular',mask)]:
                for state in ['All', *sorted(export.loc[take,'CM_state'].unique())]:
                    ix=take if state=='All' else take&export.CM_state.eq(state).to_numpy()
                    for name in ['maturation21',*COMPONENTS,*REMOVALS,'depth50_mean3']:
                        values=export.loc[ix,view+'_Z__'+name].dropna()
                        summary.append({'stage':stage,'view':view,'scope':scope,'CM_state':state,'score':name,'profiles':len(values),
                            'mean':values.mean(),'median':values.median(),'p10':values.quantile(.1),'p90':values.quantile(.9)})
        for scope,group in [('stage',None),('state','CM_state')]:
            for view in ['core','native']:
                take=export.strict_ventricular_map_cohort
                vals=export.loc[take,view+'_Z__maturation21']
                export.loc[take,view+'_'+scope+'_percentile']=(vals.rank(pct=True) if group is None else vals.groupby(export.loc[take,group]).rank(pct=True))
        export.index.name='cell_id';save(export,stage+'/CM_maturation21_scores.tsv.gz',index=True)
    save(pd.DataFrame(gene_metrics),'panel_gene_expression.tsv')
    save(pd.DataFrame(omissions),'leave_one_gene_out_correlations.tsv')
    save(pd.DataFrame(summary),'stage_state_score_summary.tsv')
    save(pd.DataFrame(checks),'independent_score_reconstruction_checks.tsv')
    save(pd.DataFrame(depth_checks),'depth_matching_checks.tsv')
    defs=[{'score':'maturation21','genes':';'.join(PANEL),'interpretation':'Direct Xenium transcriptional maturation proxy; not physiological maturity'}]
    defs += [{'score':key,'genes':';'.join(panel),'interpretation':'Measured component; cycling is separate from the 21-gene panel'} for key,panel in COMPONENTS.items()]
    defs += [{'score':key,'genes':';'.join(g for g in PANEL if g not in rem),'interpretation':'Fixed-reference omission sensitivity; no reselected controls or refitted weights'} for key,rem in REMOVALS.items()]
    save(pd.DataFrame(defs),'score_definitions.tsv')
    (OUT/'scoring_methods.json').write_text(json.dumps({'panel_genes':PANEL,'target':21,'normalization':'log1p(10000 * target counts / all 403 measured-gene counts)',
        'score':'Mean of target-gene Z values using unchanged pooled all-profile, all-age reference; separate frozen calibration for native and core views',
        'reference':'results/cardiac_spatial_reconstruction/score_reference_calibration.tsv','authors_exact_AddModuleScore_recovered':False,
        'new_control_selection':False,'one_gene_omission_reference_frozen':True,'joint_three_genes_prespecified':['Cebpb','Rxrg','Lpin1'],
        'fixed_depth':DEPTH,'depth_repetitions':3,'depth_sampling':'Exact molecule subsampling without replacement; mean score across three draws; only profiles with >=50 counts',
        'core_is_ground_truth':False,'new_biological_replicates':0},indent=2)+'\n')

def qc_capillary_pairs():
    """Refit apposition matching using QC/nuclear-supported capillary sources."""
    rows=[]
    for stage in STAGES:
        scores=read(OUT/stage/'CM_maturation21_scores.tsv.gz',index_col='cell_id')
        c=read(f'results/cardiac_spatial_reconstruction/{stage}/CM_outcomes_and_anatomy.tsv.gz',index_col=0)
        e=read(f'results/cardiac_spatial_reconstruction/{stage}/vascular_interface_and_matrix_context.tsv.gz',index_col=0)
        c=c.loc[scores.index[scores.strict_ventricular_map_cohort]].copy();e=e.loc[c.index]
        c['apposition']=e.QC_capillary__apposition_fraction_2um
        c['capillary_area']=e.QC_capillary__nearby_cell_area_fraction_50um
        x=np.column_stack([np.log1p(c.native_panel_depth),np.log1p(c.cell_area),
                           np.log1p(e.CM_native_perimeter_um),np.log1p(c.capillary_area)])
        x=(x-np.nanmean(x,axis=0))/np.maximum(np.nanstd(x,axis=0),1e-8)
        positive=c.apposition.gt(0);cut=max(.02,c.loc[positive,'apposition'].quantile(.5))
        rng=np.random.default_rng(SEED+int(stage[1:]));xy=c[['x_centroid','y_centroid']].to_numpy()
        for _,group in c.groupby(['anatomy_block_200um','author_atlas_cell_state']):
            hi=c.index.get_indexer(group.index[group.apposition.ge(cut)])
            lo=c.index.get_indexer(group.index[group.apposition.eq(0)])
            if not len(hi) or not len(lo):continue
            rng.shuffle(hi);distance=cdist(x[hi],x[lo],metric='sqeuclidean')
            physical=cdist(xy[hi],xy[lo]);capdiff=abs(x[hi,-1,None]-x[lo,-1][None,:])
            shapediff=np.max(abs(x[hi,None,:3]-x[lo,:3][None,:,:]),axis=2)
            distance[(physical>100)|(capdiff>.1)|(shapediff>.5)|(distance>4)]=np.inf
            available=np.ones(len(lo),bool)
            for k,i in enumerate(hi):
                z=np.where(available,distance[k],np.inf);j=np.argmin(z)
                if not np.isfinite(z[j]):continue
                available[j]=False;a=c.index[i];b=c.index[lo[j]]
                rows.append({'stage':stage,'high_apposition_CM':a,'zero_apposition_CM':b,
                    'anatomy_block':c.loc[a,'anatomy_block_200um'],'author_CM_state':c.loc[a,'author_atlas_cell_state'],
                    'matching_distance':z[j],'centroid_separation_um':physical[k,j],
                    'high_apposition_fraction':c.loc[a,'apposition'],'control_apposition_fraction':0,
                    'high_capillary_area':c.loc[a,'capillary_area'],'control_capillary_area':c.loc[b,'capillary_area'],
                    'source_definition':'QC broad EC, program-source eligible, nuclear-supported Capillary-like',
                    'matching_uses_maturation_score':False})
    save(pd.DataFrame(rows),'QC_capillary_refitted_pairs.tsv.gz')

def matched_capillaries():
    pairs=read('results/cardiac_spatial_reconstruction/matched_apposition_CM_pairs.tsv.gz')
    qc=read(OUT/'QC_capillary_refitted_pairs.tsv.gz')
    rows=[];pairs_out=[]
    for stage in STAGES:
        c=read(OUT/stage/'CM_maturation21_scores.tsv.gz',index_col='cell_id')
        p=pairs[pairs.stage.eq(stage)].copy()
        assert not p.high_apposition_CM.duplicated().any() and not p.zero_apposition_CM.duplicated().any()
        assert set(p.high_apposition_CM).isdisjoint(p.zero_apposition_CM)
        strict=p[p.high_apposition_CM.isin(c.index[c.strict_ventricular_map_cohort])&p.zero_apposition_CM.isin(c.index[c.strict_ventricular_map_cohort])]
        depthids=c.index[c.core_panel_depth.ge(DEPTH)&c.native_panel_depth.ge(DEPTH)]
        depthpairs=strict[strict.high_apposition_CM.isin(depthids)&strict.zero_apposition_CM.isin(depthids)]
        qcp=qc[qc.stage.eq(stage)]
        for setname,pp in [('original_matched_cohort',p),('strict_ventricular_pairs',strict),('strict_ventricular_depth50_pairs',depthpairs),('QC_capillary_strict_pairs',qcp)]:
            if len(pp)<30:continue
            ids=pp.high_apposition_CM;controls=pp.zero_apposition_CM
            assert np.array_equal(c.loc[ids,'CM_state'].to_numpy(),c.loc[controls,'CM_state'].to_numpy())
            xy=c.loc[ids,['x_centroid','y_centroid']].to_numpy();blocks=(xy//400).astype(int)
            _,inv=np.unique(blocks,axis=0,return_inverse=True);nb=inv.max()+1
            rng=np.random.default_rng(SEED+int(stage[1:]));weights=np.stack([np.bincount(rng.integers(nb,size=nb),minlength=nb) for _ in range(DRAWS)])
            for outcome in PRIMARY+SENS+['core_Z__depth50_mean3','native_Z__depth50_mean3']:
                high=c.loc[ids,outcome].to_numpy();low=c.loc[controls,outcome].to_numpy();valid=np.isfinite(high)&np.isfinite(low)
                delta=high[valid]-low[valid];codes=inv[valid]
                if len(delta)<30:continue
                sums=np.bincount(codes,weights=delta,minlength=nb);counts=np.bincount(codes,minlength=nb)
                view=outcome.split('_')[0];ref=c.loc[c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20),view+'_Z__maturation21']
                # Components retain their own full-cohort SD to match the original design.
                if outcome in PRIMARY+SENS[:5]:ref=c.loc[c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20),outcome]
                sd=ref.std();boot=(weights@sums)/np.maximum(weights@counts,1)/sd
                se=boot.std(ddof=1);point=delta.mean()/sd;lo,hi=np.quantile(boot,[.025,.975])
                rows.append({'stage':stage,'cohort':setname,'outcome':outcome,'matched_pairs':len(delta),'blocks400':len(np.unique(codes)),
                    'raw_score_difference':delta.mean(),'reference_SD':sd,'standardized_difference':point,'CI_low':lo,'CI_high':hi,
                    'within_section_p':2*norm.sf(abs(point)/max(se,1e-12)),
                    'primary_analysis':setname=='strict_ventricular_pairs' and outcome in PRIMARY,
                    'scope':'Frozen matched pairs; operational native capillary apposition; within-section uncertainty, not independent hearts'})
            pairs_out.append(pp.assign(cohort=setname))
    result=pd.DataFrame(rows);primary=result.primary_analysis
    result['primary_grid_BH_q']=np.nan;result.loc[primary,'primary_grid_BH_q']=bh(result.loc[primary,'within_section_p'])
    save(result,'matched_capillary_maturation_results.tsv');save(pd.concat(pairs_out),'matched_capillary_pairs_used.tsv.gz')

def stage_depth_audit():
    rows=[]
    for stage in STAGES:
        d=read(OUT/stage/'CM_maturation21_scores.tsv.gz',index_col='cell_id')
        d=d[d.strict_ventricular_map_cohort&d.core_panel_depth.ge(DEPTH)&d.native_panel_depth.ge(DEPTH)]
        for state in ['All',*sorted(d.CM_state.unique())]:
            z=d if state=='All' else d[d.CM_state.eq(state)]
            for view in ['core','native']:
                for score in ['maturation21','depth50_mean3','without_Cebpb_Rxrg_Lpin1']:
                    rows.append({'stage':stage,'CM_state':state,'view':view,'score':score,'same_CM_profiles':len(z),
                        'mean':z[view+'_Z__'+score].mean(),'median':z[view+'_Z__'+score].median(),
                        'scope':'Identical depth-eligible cells within age; observed ages remain separate sections'})
    save(pd.DataFrame(rows),'same_CM_depth_and_gene_omission_stage_means.tsv')
    genes=read(OUT/'panel_gene_expression.tsv');contributions=[]
    for view in ['native','core']:
        early=genes[genes.stage.eq('P0')&genes.view.eq(view)].set_index('gene')
        late=genes[genes.stage.eq('P21')&genes.view.eq(view)].set_index('gene')
        for gene in PANEL:
            contributions.append({'view':view,'gene':gene,'P21_minus_P0_gene_Z':late.loc[gene,'mean_frozen_gene_Z']-early.loc[gene,'mean_frozen_gene_Z'],
                'contribution_to_21gene_score_change':(late.loc[gene,'mean_frozen_gene_Z']-early.loc[gene,'mean_frozen_gene_Z'])/21,
                'prespecified_three_gene_omission':gene in REMOVALS['without_Cebpb_Rxrg_Lpin1'],
                'scope':'Descriptive age-section contrast in frozen strict ventricular cohorts; not a per-animal longitudinal change'})
    save(pd.DataFrame(contributions),'gene_contributions_to_stage_score_change.tsv')

def sources_and_base(stage,geometry):
    f,d=prior.stage_data(stage)
    old=read(f'results/cardiac_subtype_niches/{stage}_{geometry}_CM_subtype_neighborhoods.tsv.gz',index_col='cell_id')
    scores=read(OUT/stage/'CM_maturation21_scores.tsv.gz',index_col='cell_id')
    for col in PRIMARY+SENS+['core_Z__depth50_mean3','native_Z__depth50_mean3']:
        old[col]=scores.loc[old.index,col]
    old=old.copy()
    eligible=f.program_source_eligible.fillna(False);resolved=f.atlas_resolved.fillna(False)
    fb=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
    coverage=read('results/cardiac_fibroblast_CM_refit/fibroblast_assignment_coverage.tsv')
    assessable=set(coverage.loc[coverage.stage.eq(stage)&coverage.label_layer.eq('candidate')&coverage.reference_assessable,'family'])
    masks={};meta=[]
    for name in BROAD:masks['Broad__'+name]=f.atlas_identity.eq(name)&eligible
    for name in prior.EC:masks['EC__'+name]=f.atlas_identity.eq('Endothelial')&eligible&f.EC_nuclear_type.eq(name)
    for name in FAMILIES:
        m=pd.Series(False,index=f.index)
        m.loc[fb.index]=fb.program_source_eligible&fb.recommended_candidate.eq(name+'-associated')
        masks['FB__'+name]=m
    graph=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/{geometry}_graph.npz').tocsr()[f.index.get_indexer(old.index)]
    INPUTS.add(str(ROOT/f'results/cardiac_spatial_pipeline/{stage}/{geometry}_graph.npz'))
    all_n=np.asarray(graph@resolved.to_numpy(float)).ravel()
    fb_n=np.asarray(graph@(f.atlas_identity.eq('Fibroblast')&eligible).to_numpy(float)).ravel()
    ec_n=np.asarray(graph@(f.atlas_identity.eq('Endothelial')&eligible).to_numpy(float)).ravel()
    coords=f[['x_centroid','y_centroid']].to_numpy();cmcoords=old[['x_centroid','y_centroid']].to_numpy()
    for name,m in masks.items():
        denom=fb_n if name.startswith('FB__') else ec_n if name.startswith('EC__') else all_n
        count=np.asarray(graph@m.to_numpy(float)).ravel()
        old[name]=count/np.maximum(denom,1)
        if geometry=='proximity50':
            old['distance__'+name]=cKDTree(coords[m]).query(cmcoords)[0] if m.any() else np.nan
        available=name.split('__')[1] in assessable if name.startswith('FB__') else True
        meta.append({'stage':stage,'geometry':geometry,'source':name,'source_profiles':int(m.sum()),
                     'reference_assessable':available,'CMs_with_source':int((count>0).sum()),
                     'denominator':'All eligible broad fibroblasts' if name.startswith('FB__') else 'All eligible broad endothelial cells' if name.startswith('EC__') else 'All resolved neighbors',
                     'supported_for_model':available and int(m.sum())>=30 and int((count>0).sum())>=50})
    return f,old.copy(),pd.DataFrame(meta)

def nuisance(d,target):
    cov={'log_CM_capture':np.log1p(d.CM_outcome_excluded_depth),'log_CM_area':np.log1p(d.cell_area),
         'nucleus_missing':d.nucleus_area.isna().astype(float),'log_neighbors':np.log1p(d.neighbors_resolved),
         'log_FB_count':np.log1p(d.eligible_FB_neighbors),'log_EC_count':np.log1p(d.eligible_EC_neighbors),
         'FB_capture':d.FB_mean_log_capture,'EC_capture':d.EC_mean_log_capture}
    for col in ['arterial_distance_um','endocardial_proxy_distance_um','distance_to_large_void_um','distance_to_outer_um']:
        cov['log_'+col]=np.log1p(d[col])
    for key in ['Broad__'+b for b in BROAD]+['EC__Capillary-like','EC__Endocardial-like']:
        if key!=target:cov[key]=d[key]
    for name in sorted(d.morphology_region.unique())[1:]:cov['morphology__'+name]=d.morphology_region.eq(name).astype(float)
    return pd.DataFrame(cov,index=d.index)

def fit(d,target,geometry,mode='primary'):
    keys=PRIMARY+SENS if mode=='primary' else ['core_Z__maturation21','core_Z__depth50_mean3','core_Z__without_Cebpb_Rxrg_Lpin1']
    X=nuisance(d,target).to_numpy(float);e=d[target].to_numpy(float);Y=d[keys].to_numpy(float)
    finite=np.isfinite(X).all(1)&np.isfinite(Y).all(1)&np.isfinite(e)
    if finite.sum()<200:return []
    # Common stage scales across CM states and subset sensitivity fits.
    X=(X-np.nanmean(X,axis=0))/np.maximum(np.nanstd(X,axis=0),1e-6)
    e=(e-np.nanmean(e))/max(np.nanstd(e),1e-6)
    Y=(Y-np.nanmean(Y,axis=0))/np.maximum(np.nanstd(Y,axis=0),1e-6)
    blocks=sorted(d.loc[finite,'block400'].unique());nb=len(blocks)
    if nb<8:return []
    rng=np.random.default_rng(SEED+int(d.stage.iloc[0][1:]));w=np.stack([np.bincount(rng.integers(nb,size=nb),minlength=nb) for _ in range(DRAWS)])
    lookup={b:i for i,b in enumerate(blocks)};rows=[]
    for state in ['All_state_adjusted', 'v.CM', 'Ankrd1.CM', 'p.CM', 'Slit2.CM', 'Ddc.CM']:
        take=finite if state=='All_state_adjusted' else finite&d.CM_state.eq(state).to_numpy()
        z=d.loc[take]
        if len(z)<200 or z.block400.nunique()<8 or (z[target]>0).sum()<50:continue
        groups=(z.block200+'::'+z.CM_state).to_numpy()
        xx=prior.center(X[take],groups);ee=prior.center(e[take],groups).ravel();yy=prior.center(Y[take],groups)
        keep=np.flatnonzero(xx.std(0)>.015)
        if len(keep):
            _,r,piv=linalg.qr(xx[:,keep],mode='economic',pivoting=True)
            rank=(abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-8).sum();keep=keep[piv[:rank]]
        xx=xx[:,keep]
        residual=ee-xx@np.linalg.lstsq(xx,ee,rcond=1e-8)[0]
        if residual.std()<.02:continue
        design=np.column_stack([xx,ee]);codes=np.array([lookup[b] for b in z.block400]);p=design.shape[1]
        xtx=np.zeros((nb,p,p));xty=np.zeros((nb,p,len(keys)))
        for code in np.unique(codes):
            m=codes==code;xtx[code]=design[m].T@design[m];xty[code]=design[m].T@yy[m]
        beta=np.linalg.pinv(xtx.sum(0),rcond=1e-9)@xty.sum(0)
        boot=np.linalg.pinv(np.einsum('bn,nij->bij',w,xtx),rcond=1e-8)@np.einsum('bn,nij->bij',w,xty)
        for j,key in enumerate(keys):
            bs=boot[:,-1,j];lo,hi=np.quantile(bs,[.025,.975]);se=bs.std(ddof=1);b=beta[-1,j]
            rows.append({'stage':d.stage.iloc[0],'geometry':geometry,'analysis_set':mode,'source':target,'CM_state':state,'outcome':key,
                'coefficient':b,'CI_low':lo,'CI_high':hi,'within_section_p':2*norm.sf(abs(b)/max(se,1e-12)),
                'CM_profiles':len(z),'blocks400':z.block400.nunique(),'source_residual_SD':residual.std(),
                'condition_number':np.linalg.cond(design),
                'primary_analysis':geometry=='proximity50' and mode=='primary' and state=='All_state_adjusted' and key in PRIMARY,
                'scope':'Within 200-um region and CM identity; technical, broad-composition and anatomy adjustment; one section per age'})
    return rows

def models():
    rows=[];coverage=[];means=[]
    for stage in STAGES:
        for geometry in ['proximity50','proximity30']:
            print('MATURATION21_NEIGHBORS',stage,geometry,flush=True)
            f,d,cov=sources_and_base(stage,geometry);coverage.append(cov)
            if geometry=='proximity50':
                export=['stage','CM_state','x_centroid','y_centroid','analysis_eligible','core_Z__maturation21','native_Z__maturation21']
                export += cov.source.tolist()+['distance__'+s for s in cov.source]
                save(d[export],stage+'/CM_nonCM_neighborhoods.tsv.gz',index=True)
            z=d[d.analysis_eligible].copy()
            for source in cov.loc[cov.supported_for_model,'source']:
                rows+=fit(z,source,geometry)
                if geometry=='proximity50':rows+=fit(z[z.core_panel_depth.ge(DEPTH)].copy(),source,geometry,'depth50_same_CMs')
                if geometry=='proximity50':
                    for state,g in z.groupby('CM_state'):
                        rank=g.core_Z__maturation21.rank(pct=True)
                        for group,take in [('Lower_quartile',rank.le(.25)),('Upper_quartile',rank.ge(.75))]:
                            h=g[take]
                            means.append({'stage':stage,'CM_state':state,'score_group':group,'source':source,'CM_profiles':len(h),
                                'mean_neighbor_fraction':h[source].mean(),'fraction_with_source_neighbor':h[source].gt(0).mean(),
                                'median_nearest_source_um':h['distance__'+source].median(),'scope':'Unmatched descriptive values, not independent effects'})
    result=pd.DataFrame(rows);primary=result.primary_analysis
    result['primary_grid_BH_q']=np.nan;result.loc[primary,'primary_grid_BH_q']=bh(result.loc[primary,'within_section_p'])
    save(result,'nonCM_maturation_associations.tsv');save(pd.concat(coverage),'nonCM_source_coverage.tsv')
    save(pd.DataFrame(means),'maturation_quartile_neighborhood_descriptions.tsv')

def predictive_comparison():
    """Incremental prediction in held-out spatial blocks, without model tuning."""
    rows=[];checks=[]
    for stage in STAGES:
        print('MATURATION21_PREDICTION',stage,flush=True)
        _,d,cov=sources_and_base(stage,'proximity50')
        with np.load(ROOT/f'results/cardiac_fibroblast_CM_refit/{stage}_proximity50_neighborhood_exposures.npz') as q:
            INPUTS.add(str(ROOT/f'results/cardiac_fibroblast_CM_refit/{stage}_proximity50_neighborhood_exposures.npz'))
            assert np.array_equal(q['cell_ids'],d.index.to_numpy(str))
            for col in q.files:
                if col.startswith('refined__'):d[col]=q[col]
        d=d[d.analysis_eligible&d.refined__informative_neighbors.ge(1)].copy()
        cols=['log_CM_capture','log_CM_area','nucleus_missing','log_neighbors','log_FB_count','log_EC_count','FB_capture','EC_capture']
        base=nuisance(d,'__no_target__')
        # Avoid leaking fine EC identity into the broad baseline.
        base=base.drop(columns=['EC__Capillary-like','EC__Endocardial-like'])
        states=pd.get_dummies(d.CM_state,prefix='CM_state',dtype=float)
        base=pd.concat([base,states],axis=1)
        subtype=[s for s in cov.loc[cov.supported_for_model,'source'] if not s.startswith('Broad__')]
        programs=[c for c in d if c.startswith('refined__') and c not in ['refined__information_fraction','refined__informative_neighbors','refined__log_depth','refined__detected_genes']]
        capture=pd.concat([base[['log_CM_capture','log_CM_area','nucleus_missing']],states],axis=1)
        sets={'CM_capture_size_identity':capture,
              'anatomy_technical_broad':base,
              'plus_EC_FB_identities':pd.concat([base,d[subtype]],axis=1),
              'plus_FB_measured_programs':pd.concat([base,d[subtype+programs]],axis=1)}
        xy=d[['x_centroid','y_centroid']].to_numpy();Y=d[PRIMARY].to_numpy()
        valid=np.isfinite(Y).all(1)
        for x in sets.values():valid &= np.isfinite(x.to_numpy()).all(1)
        d=d.loc[valid];xy=xy[valid];Y=Y[valid];sets={k:v.loc[d.index].to_numpy() for k,v in sets.items()}
        blocks=(xy//400).astype(int);fold=((blocks[:,0]*73856093)^(blocks[:,1]*19349663)^271)%5
        preds={k:np.full_like(Y,np.nan) for k in sets}
        for k in range(5):
            test=fold==k;train=~test
            if not test.any():continue
            near=cKDTree(xy[test]).query(xy[train])[0]<150
            train[np.flatnonzero(train)[near]]=False
            if train.sum()<300 or test.sum()<100:continue
            # Smooth coordinate basis is fitted to training coordinates only.
            spline=SplineTransformer(n_knots=4,degree=2,include_bias=False).fit(xy[train])
            space_train=spline.transform(xy[train]);space_test=spline.transform(xy[test])
            for name,X in sets.items():
                a=X[train] if name=='CM_capture_size_identity' else np.column_stack([X[train],space_train])
                b=X[test] if name=='CM_capture_size_identity' else np.column_stack([X[test],space_test])
                scale=StandardScaler().fit(a);model=Ridge(alpha=10).fit(scale.transform(a),Y[train])
                preds[name][test]=model.predict(scale.transform(b))
            checks.append({'stage':stage,'fold':k,'training_CMs':int(train.sum()),'test_CMs':int(test.sum()),
                           'minimum_train_test_distance_um':cKDTree(xy[test]).query(xy[train])[0].min()})
        for name,pred in preds.items():
            take=np.isfinite(pred).all(1)
            if take.sum()<100:continue
            for j,outcome in enumerate(PRIMARY):
                rows.append({'stage':stage,'model':name,'outcome':outcome,'predicted_CMs':int(take.sum()),
                             'spatially_heldout_R2':r2_score(Y[take,j],pred[take,j]),'fixed_ridge_alpha':10,
                             'scope':'Buffered spatial-block prediction within a single section; not held-out animals or causal influence'})
    save(pd.DataFrame(rows),'incremental_spatial_prediction.tsv');save(pd.DataFrame(checks),'prediction_fold_checks.tsv')

def fibroblast_evidence():
    """Summarize measured state changes without introducing a forced lineage."""
    ps=read('results/fibroblast_states/RNA_corrected_matrix_pseudobulk_scores.tsv')
    ps=ps[ps.version.eq('full')&~ps.family.str.startswith('Unresolved')]
    save(ps,'fibroblast_library_matrix_scores.tsv')
    ge=read('results/fibroblast_states/RNA_within_family_stage_gene_changes.tsv')
    genes=['Fn1','Postn','Tnc','Col1a1','Col1a2','Col3a1','Col8a1','Gfpt2','Pi16','Meox1','Htra3','Lsamp','Thbs1','Vcan']
    ge=ge[ge.gene.isin(genes)].copy();ge['library_direction_agreement']=np.sign(ge.minimum_library_cross_difference)==np.sign(ge.maximum_library_cross_difference)
    save(ge,'fibroblast_matrix_gene_stage_changes.tsv')
    signature=read('results/cardiac_developmental_story/reviewed_fibroblast_candidate_state_continuities.tsv')
    save(signature,'fibroblast_candidate_expression_continuities.tsv')
    mark=read('results/fibroblast_multimethod_identity/RNA_reference_subcluster_marker_audit.tsv')
    save(mark[mark.gene.isin(genes)],'fibroblast_marker_overlap_audit.tsv')
    cov=read('results/cardiac_fibroblast_CM_refit/fibroblast_assignment_coverage.tsv')
    save(cov[cov.label_layer.eq('candidate')],'fibroblast_assignment_coverage.tsv')
    spatial=[]
    for stage in STAGES:
        fb=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
        anatomy=read(f'results/cardiac_spatial_reconstruction/{stage}/all_profile_observations.tsv.gz',index_col=0)
        physical=read(f'results/fibroblast_names_and_neighborhoods/{stage}/FB_physical_neighbors.tsv.gz',index_col='cell_id')
        for family in FAMILIES:
            available=bool(cov.loc[cov.stage.eq(stage)&cov.label_layer.eq('candidate')&cov.family.eq(family),'reference_assessable'].iloc[0])
            ids=fb.index[fb.program_source_eligible&fb.recommended_candidate.eq(family+'-associated')]
            ventricular=ids[anatomy.loc[ids,'ventricular_morphology_selector'].to_numpy(bool)]
            q=anatomy.loc[ventricular];p=physical.loc[ids]
            spatial.append({'stage':stage,'family':family,'reference_assessable':available,
                'eligible_candidate_profiles':len(ids) if available else np.nan,
                'within_ventricular_selector_profiles':len(ventricular) if available else np.nan,
                'median_void_distance_um':q.distance_to_large_void_um.median(),
                'median_outer_distance_um':q.distance_to_outer_um.median(),
                'median_relative_wall_geometry':q.relative_wall_geometry.median(),
                'fraction_within50um_arterial_structure':p.arterial_distance_um.le(50).mean() if len(p) else np.nan,
                'median_nearest_CM_um':p.nearest_Cardiomyocyte_um.median(),
                'scope':'Descriptive snapshots; void/outer geometry is not a validated endocardium-to-epicardium coordinate or tracked migration'})
    save(pd.DataFrame(spatial),'fibroblast_spatial_geometry_snapshots.tsv')
    # Col8a1 is deliberately retained as one family in every new spatial model.
    (OUT/'fibroblast_interpretation.json').write_text(json.dumps({'Col8a1_family_split':False,
        'Gfpt2_Pi16_name':'Gfpt2/Pi16-associated; multi-gene family, not obligatory double positivity',
        'new_fate_trajectory_fitted':False,'existing_continuities':'Marker-withheld cross-stage expression resemblance, not cell ancestry',
        'Thbs1_Vcan_valve_location_verified':False,'Lsamp_Htra3_outflow_location_verified':False,
        'CM_TGFbeta_causes_FB_transition_established':False,'wall_growth_excludes_transition':False},indent=2)+'\n')

def validate():
    checks=read(OUT/'independent_score_reconstruction_checks.tsv')
    assert len(checks)==40 and checks.maximum_absolute_error.max()<2e-6
    loo=read(OUT/'leave_one_gene_out_correlations.tsv');assert len(loo)==336 and loo.normalization_and_reference_frozen.all()
    assoc=read(OUT/'nonCM_maturation_associations.tsv');primary=assoc[assoc.primary_analysis]
    assert not primary.duplicated(['stage','source','outcome']).any()
    np.testing.assert_allclose(primary.primary_grid_BH_q,bh(primary.within_section_p),atol=1e-8)
    assert primary.CM_profiles.min()>=200 and primary.blocks400.min()>=8
    assert np.isfinite(assoc.coefficient).all() and (assoc.CI_low<=assoc.CI_high).all()
    assert assoc.condition_number.max()<1e8
    pairs=read(OUT/'matched_capillary_maturation_results.tsv')
    old=read('results/cardiac_spatial_reconstruction/matched_apposition_component_differences.tsv')
    comparisons=0
    for _,r in pairs[pairs.cohort.eq('original_matched_cohort')&pairs.outcome.isin(PRIMARY)].iterrows():
        rr=old[old.stage.eq(r.stage)&old.outcome.eq(r.outcome)]
        if len(rr):
            assert abs(r.standardized_difference-rr.iloc[0].mean_standardized_difference)<1e-6
            comparisons+=1
    cv=read(OUT/'prediction_fold_checks.tsv');assert len(cv)==20 and cv.minimum_train_test_distance_um.min()>=150-1e-7
    available=read(OUT/'fibroblast_assignment_coverage.tsv');assert available.loc[~available.reference_assessable,'candidate_profiles'].isna().all()
    # Independent scalar matrix/gene reconstruction checks for each stage/view.
    validation={'status':'PASS','recomputed_frozen_score_checks':len(checks),'maximum_score_reconstruction_error':checks.maximum_absolute_error.max(),
        'single_gene_omission_checks':len(loo),'old_matched_point_estimates_reproduced':comparisons,
        'new_primary_nonCM_tests':len(primary),'new_primary_nonCM_q_lt_0_05':int(primary.primary_grid_BH_q.lt(.05).sum()),
        'primary_matched_capillary_tests':int(pairs.primary_analysis.sum()),'buffered_spatial_prediction_folds':len(cv),
        'no_Col8a1_family_subdivision':True,'fate_transition_claimed':False,'new_biological_replicates':0}
    (OUT/'analysis_validation.json').write_text(json.dumps(validation,indent=2,default=float)+'\n')
    print('MATURATION21_ANALYSIS_PASS',json.dumps(validation,default=float),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['all','prepare','qc_pairs','matched','stage_audit','models','prediction','fibroblast','validate'],default='all');a=p.parse_args()
    OUT.mkdir(exist_ok=True,parents=True);start=time.monotonic()
    for name,func in [('prepare',preparation),('qc_pairs',qc_capillary_pairs),('matched',matched_capillaries),('stage_audit',stage_depth_audit),('models',models),('prediction',predictive_comparison),('fibroblast',fibroblast_evidence),('validate',validate)]:
        if a.step in ['all',name]:func()
    manifest=OUT/'input_paths.json';existing=set(json.loads(manifest.read_text())) if manifest.exists() else set()
    existing|=INPUTS|prior.SOURCES
    manifest.write_text(json.dumps(sorted(existing),indent=2)+'\n')
    print('MATURATION21_COMPLETE',a.step,round(time.monotonic()-start,1),flush=True)

if __name__=='__main__':main()
