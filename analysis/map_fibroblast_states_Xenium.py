"""Transfer reassessed fibroblast programs and compute measured CM/neighborhood data."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-state-mpl')
import argparse
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.metrics import precision_recall_fscore_support
from improve_ventricular_alignment import neighbors
from fibroblast_state_common import *

SOURCE = ROOT/'results/xenium_ventricular_mapping'
ATLAS = ROOT/'results/whole_heart_spatial_atlas'
VESSELS = ROOT/'results/whole_heart_Ankrd1_vessels'


def norm_features(counts, use):
    return np.log2(1+1e6*counts[:,use]/np.maximum(counts[:,use].sum(1),1)[:,None]).astype(np.float32)


class FibroblastProjection:
    def __init__(self,stage,qcounts,qgenes,ann,withhold_matrix=False):
        meta = pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
        with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
            rids, rgenes = b['cell_ids'].astype(str),b['genes'].astype(str)
            valid_ids=ann.index[ann.fine_reference_eligible]
            keep = meta.age.eq(stage).to_numpy() & meta.index.isin(valid_ids)
            self.rids = rids[keep]
            self.rmeta = ann.loc[self.rids]
            self.labels = self.rmeta.fibroblast_family.to_numpy(dtype=str)
            original,corrected = b['original'][keep],b['corrected'][keep]
        assert np.array_equal(rids,meta.index)
        self.classes = FAMILIES+[UNRESOLVED]
        self.qgenes = np.asarray(qgenes)
        positions = pd.Index(qgenes).get_indexer(rgenes)
        assert (positions>=0).all()
        self.positions = positions
        exclude = set(WITHHOLD)
        if withhold_matrix:exclude |= set(sum(MATRIX_PANELS.values(),[]))
        self.use = ~np.isin(rgenes,list(exclude))
        # All same-age fibroblast references compete, including unresolved families.
        self.rgenes = rgenes
        self.calibrate(qcounts)
        self.models = {}
        self.library_names = sorted(self.rmeta['sample'].unique())
        self.train_counts = corrected
        for name, counts, selection in [('primary',corrected,np.ones(len(corrected),bool)),
                ('raw_reference',original,np.ones(len(corrected),bool))]+[
                (lib,corrected,self.rmeta['sample'].eq(lib).to_numpy()) for lib in self.library_names]:
            rn = norm_features(counts[selection],self.use)
            rf = np.clip((rn-rn.mean(0))/np.maximum(rn.std(0),.25),-5,5)
            rf /= np.maximum(np.linalg.norm(rf,axis=1,keepdims=True),1e-8)
            self.models[name] = (rf.astype(np.float32),self.labels[selection])

    def calibrate(self,qcounts):
        qn = norm_features(qcounts[:,self.positions],self.use)
        self.qmu,self.qsd = qn.mean(0),np.maximum(qn.std(0),.25)

    def predict(self,qcounts,models=None):
        qn = norm_features(qcounts[:,self.positions],self.use)
        qf = np.clip((qn-self.qmu)/self.qsd,-5,5)
        qf /= np.maximum(np.linalg.norm(qf,axis=1,keepdims=True),1e-8)
        results = {}
        for name in models or self.models:
            rf,labels = self.models[name]
            ix, weights, similarity = neighbors(rf,qf.astype(np.float32))
            p = np.column_stack([(weights*(labels[ix]==fam)).sum(1) for fam in self.classes])
            assert np.allclose(p.sum(1),1,atol=1e-5)
            results[name]=(p,similarity)
        return results

    def benchmark(self):
        rows=[]
        for lib in self.library_names:
            query=self.rmeta['sample'].eq(lib).to_numpy()
            train=~query
            rn=norm_features(self.train_counts[train],self.use)
            qn=norm_features(self.train_counts[query],self.use)
            rf=np.clip((rn-rn.mean(0))/np.maximum(rn.std(0),.25),-5,5)
            qf=np.clip((qn-qn.mean(0))/np.maximum(qn.std(0),.25),-5,5)
            rf/=np.maximum(np.linalg.norm(rf,axis=1,keepdims=True),1e-8)
            qf/=np.maximum(np.linalg.norm(qf,axis=1,keepdims=True),1e-8)
            ix,w,_=neighbors(rf.astype(np.float32),qf.astype(np.float32))
            p=np.column_stack([(w*(self.labels[train][ix]==fam)).sum(1) for fam in self.classes])
            pred=np.asarray(self.classes)[p.argmax(1)]
            pr,re,f1,n=precision_recall_fscore_support(self.labels[query],pred,labels=self.classes,zero_division=0)
            for fam,a,b,c,d in zip(self.classes,pr,re,f1,n):
                rows.append({'stage':self.rmeta.age.iloc[0],'query_library':lib,'family':fam,'precision':a,'recall':b,'f1':c,'nuclei':int(d)})
        return rows


def annotate(qcounts,genes,projection,pred):
    p=pred['primary'][0]
    order=np.argsort(p,axis=1)
    top=np.asarray(projection.classes)[order[:,-1]]
    margin=p[np.arange(len(p)),order[:,-1]]-p[np.arange(len(p)),order[:,-2]]
    raw=np.asarray(projection.classes)[pred['raw_reference'][0].argmax(1)]
    both=np.ones(len(p),bool)
    for lib in projection.library_names:both&=np.asarray(projection.classes)[pred[lib][0].argmax(1)]==top
    logs=lognorm(qcounts)
    support=np.zeros(len(p),int);anchor=np.zeros(len(p),int)
    for fam,panel in PANELS.items():
        m=top==fam
        ix=np.flatnonzero(np.isin(genes,panel));ai=np.flatnonzero(np.isin(genes,ANCHORS[fam]))
        support[m]=(qcounts[m][:,ix]>0).sum(1)
        anchor[m]=(qcounts[m][:,ai]>0).sum(1)
    shared=qcounts[:,projection.positions][:,projection.use]
    informative=(shared.sum(1)>=20)&((shared>0).sum(1)>=8)
    strict=(p.max(1)>=.6)&(margin>=.15)&both&(raw==top)&(support>=2)&(anchor>=1)&informative&(top!=UNRESOLVED)
    d=pd.DataFrame({'family_proposal':top,'family_weight':p.max(1),'family_margin':margin,
        'both_library_agreement':both,'raw_reference_agreement':raw==top,'measured_state_markers':support,
        'measured_state_anchors':anchor,'state_information_QC':informative,'supported_family_assignment':strict,
        'mapped_family':np.where(strict,top,UNRESOLVED)})
    for k,fam in enumerate(projection.classes):d['weight__'+fam]=p[:,k]
    return d


def direct_programs(counts,genes):
    logs=lognorm(counts)
    d=pd.DataFrame()
    for name,panel in MATRIX_PANELS.items():
        gl=[g for g in panel if g in genes and g not in ANNOTATION_GENES]
        if not gl:continue
        d['matrix__'+name]=logs[:,np.isin(genes,gl)].mean(1)
    return d


def cm_scores(counts,genes):
    logs=lognorm(counts)
    excluded=np.isin(genes,WITHHOLD)
    without=lognorm(counts,excluded)
    mat=np.isin(genes,MATURATION)
    assert mat.sum()==21
    ank=list(genes).index('Ankrd1')
    four=np.isin(genes,['Ankrd1','Creb5','Tlr4','Atf3'])
    c=pd.DataFrame({'maturation_21gene_mean':logs[:,mat].mean(1),
        'maturation_excluded_denominator':without[:,mat].mean(1),
        'Ankrd1_log1pCP10k':np.log1p(1e4*counts[:,ank]/np.maximum(counts[:,~four].sum(1),1)),
        'Ankrd1_full_denominator':logs[:,ank],
        'outcome_excluded_depth':counts[:,~excluded].sum(1),
        'Slit2_log':logs[:,list(genes).index('Slit2')] if 'Slit2' in genes else 0,
        'cycling_log':logs[:,np.isin(genes,CYCLING)].mean(1)})
    return c,logs


def neighbors_for_cm(cm,fb,allobs,weights,programs):
    xy=cm[['x_centroid','y_centroid']].to_numpy()
    fxy=fb[['x_centroid','y_centroid']].to_numpy()
    edges=cKDTree(xy).sparse_distance_matrix(cKDTree(fxy),100,output_type='coo_matrix')
    result={}
    for radius in [30,50,100]:
        sel=edges.data<=radius
        mat=sparse.csr_matrix((np.ones(sel.sum()),(edges.row[sel],edges.col[sel])),shape=(len(cm),len(fb)))
        nf=np.asarray(mat.sum(1)).ravel()
        cols={'fibroblast_neighbors':nf}
        for j,fam in enumerate(FAMILIES):
            cols['family__'+fam]=np.asarray(mat@weights[:,j]).ravel()/np.maximum(nf,1)
            hard=fb.mapped_family.eq(fam).to_numpy(dtype=float)
            cols['hardfamily__'+fam]=np.asarray(mat@hard).ravel()/np.maximum(nf,1)
        for name in programs:
            cols[name]=np.asarray(mat@programs[name].to_numpy()).ravel()/np.maximum(nf,1)
        resolved=allobs.atlas_resolved.to_numpy(dtype=bool)
        total=cKDTree(allobs.loc[resolved,['x_centroid','y_centroid']].to_numpy()).query_ball_point(xy,radius,return_length=True)-1
        cols['resolved_neighbors']=total
        cols['fibroblast_fraction']=nf/np.maximum(total,1)
        for name,lineages in [('endothelial_fraction',['Endothelial']),('mural_fraction',['Pericyte','Smooth muscle']),('immune_fraction',['Immune'])]:
            mask=allobs.atlas_identity.isin(lineages)
            num=cKDTree(allobs.loc[mask,['x_centroid','y_centroid']].to_numpy()).query_ball_point(xy,radius,return_length=True)
            cols[name]=num/np.maximum(total,1)
        result[radius]=pd.DataFrame(cols,index=cm.index)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stages',nargs='+',default=STAGES);args=parser.parse_args()
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    benchmarks=[];calibrations={}
    for stage in args.stages:
        folder=OUT/stage
        if (folder/'spatial_preparation_complete.json').exists():
            print('SKIP_SPATIAL_PREP',stage,flush=True);continue
        print('MAP_STAGE',stage,flush=True)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        genes=data.var_names.to_numpy(dtype=str)
        atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
        np.testing.assert_array_equal(atlas.index,data.obs_names)
        fm=atlas.atlas_identity.eq('Fibroblast').to_numpy()
        fb=atlas.loc[fm].copy();fc=data.X[fm].toarray().astype(np.float32)
        model=FibroblastProjection(stage,fc,genes,ann)
        pred=model.predict(fc)
        labels=annotate(fc,genes,model,pred);labels.index=fb.index
        fb=fb.join(labels)
        # A separate fit with all matrix program genes withheld tests outcome reuse.
        alt=FibroblastProjection(stage,fc,genes,ann,withhold_matrix=True)
        altp=alt.predict(fc,models=['primary'])['primary'][0]
        for k,fam in enumerate(alt.classes):fb['matrix_withheld_weight__'+fam]=altp[:,k]
        fb['matrix_withheld_top_agrees']=np.asarray(alt.classes)[altp.argmax(1)]==fb.family_proposal.to_numpy()
        programs=direct_programs(fc,genes);programs.index=fb.index
        fb=fb.join(programs)
        fb.to_csv(folder/'Xenium_fibroblast_states.tsv.gz',sep='\t',na_rep='NA')
        np.savez_compressed(folder/'fibroblast_projection_calibration.npz',genes=genes,reference_ids=model.rids,
            training_genes=model.rgenes[model.use],query_mean=model.qmu,query_sd=model.qsd)
        benchmarks.extend(model.benchmark())
        cm=atlas.loc[atlas.ventricular_CM_strict].copy()
        cc=data.X[atlas.ventricular_CM_strict.to_numpy()].toarray().astype(np.float32)
        scores,logs=cm_scores(cc,genes);scores.index=cm.index
        cm=cm[['x_centroid','y_centroid','cell_area','nucleus_area','stage','Ankrd1_log1pCP10k']].drop(columns='Ankrd1_log1pCP10k').join(scores)
        cm['nuclear_area_missing']=cm.nucleus_area.isna()
        cm['total_gene_transcripts']=cc.sum(1)
        cm['detected_genes']=(cc>0).sum(1)
        # Fixed 403-gene calibration is finalized across ages after all stage preparations.
        np.savez_compressed(folder/'CM_maturation_gene_means.npz',genes=genes,mean_log=logs.mean(0))
        matrix=logs[:,np.isin(genes,MATURATION)].astype(np.float32)
        np.savez_compressed(folder/'CM_21gene_expression.npz',cell_ids=cm.index.to_numpy(dtype=str),genes=genes[np.isin(genes,MATURATION)],log1pCP10k=matrix)
        # Distances are original independent vascular envelopes and author endocardial proxies.
        dist=pd.read_csv(VESSELS/stage/'arterial_boundary_distances.tsv.gz',sep='\t',index_col=0)
        cm['arterial_distance_um']=dist.loc[cm.index,'independent_strict__boundary_distance_um']
        cm['arterial_segment_id']=dist.loc[cm.index,'independent_strict__boundary_segment_id']
        endo=atlas.author_atlas_cell_state.eq('endocardial.EC')
        cm['endocardial_proxy_distance_um']=cKDTree(atlas.loc[endo,['x_centroid','y_centroid']].to_numpy()).query(cm[['x_centroid','y_centroid']])[0]
        weights=np.column_stack([fb['weight__'+fam].to_numpy() for fam in FAMILIES])
        neighborhoods=neighbors_for_cm(cm,fb,atlas,weights,programs)
        for radius,table in neighborhoods.items():
            table.to_csv(folder/f'CM_neighborhoods_{radius}um.tsv.gz',sep='\t',na_rep='NA')
        cm.to_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',na_rep='NA')
        write_json(folder/'spatial_preparation_complete.json',{'stage':stage,'fibroblasts':len(fb),'strict_ventricular_CM_profiles':len(cm),
            'supported_fine_fibroblast_assignments':int(fb.supported_family_assignment.sum()),
            'gene_targets':len(genes),'maturation_targets':21,'training_genes':model.rgenes[model.use].tolist(),
            'outcome_genes_excluded_from_transfer_and_its_denominator':WITHHOLD,
            'strict_fine_assignment':'weight>=0.6,margin>=0.15,>=2 measured paper markers,>=1 anchor,raw/both-library agreement,>=20 training molecules/>=8 training genes',
            'all_broad_atlas_labels_frozen':True,'nuclear_area_gate':False})
        print('MAP_DONE',stage,len(fb),int(fb.supported_family_assignment.sum()),len(cm),flush=True)
    if benchmarks:save_tsv(pd.DataFrame(benchmarks),OUT/'RNA_opposite_library_shared_panel_recovery.tsv')
    # Equal age contribution sets expression bins; all 21 target genes are fixed in advance.
    allmeans=[]
    for stage in STAGES:
        with np.load(OUT/stage/'CM_maturation_gene_means.npz') as b:
            genes=b['genes'].astype(str);allmeans.append(b['mean_log'])
    mean=np.mean(allmeans,axis=0)
    order=np.argsort(mean,kind='stable');bins=np.empty(len(genes),int)
    for j,ix in enumerate(np.array_split(order,8)):bins[ix]=j
    rng=np.random.default_rng(42);control=[]
    for g in MATURATION:
        gi=list(genes).index(g);pool=np.flatnonzero((bins==bins[gi])&~np.isin(genes,WITHHOLD))
        control.extend(rng.choice(pool,min(10,len(pool)),replace=False).tolist())
    control=np.unique(control);controlgenes=genes[control]
    assert not set(controlgenes)&set(WITHHOLD)
    save_tsv(pd.DataFrame({'gene':genes,'expression_bin':bins,'balanced_stage_mean_log':mean,'maturation_gene':np.isin(genes,MATURATION),'matched_control':np.isin(genes,controlgenes)}),OUT/'CM_maturation_fixed_controls.tsv')
    for stage in STAGES:
        folder=OUT/stage
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad',backed='r')
        cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
        counts=data.X[data.obs_names.get_indexer(cm.index)].toarray()
        logs=lognorm(counts)
        cm['maturation_matched_control']=cm.maturation_21gene_mean.to_numpy()-logs[:,control].mean(1)
        cm.to_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',na_rep='NA')
        data.file.close()
    write_json(OUT/'Xenium_methods.json',{'native_gene_targets':403,'RNA_reference_features':381,
        'maturation_primary':'Mean log1p CP10k of the exact 21 Figure4b genes; not claimed identical to precomputed author index',
        'matched_control_sensitivity':'Fixed 8 expression bins,10 controls per target, seed42; unique control union; balanced-stage calibration; outcome genes excluded',
        'maturation_denominator_sensitivity':'21-gene mean using total excluding all maturation+Ankrd1/Creb5/Tlr4/Atf3 genes',
        'Ankrd1':'Measured log1pCP10k; denominator excludes four prior outcome genes; full-total sensitivity saved',
        'neighborhoods_um':[30,50,100],'primary_radius_um':50,
        'matrix_exposures':'Measured FB-only RNA; exact annotation marker genes withheld from tested ECM panels',
        'source_QC':'Original broad resolved atlas, frozen before fine-state/maturation fitting; inherited lineage gates use some maturation markers and are not fully outcome-independent',
        'fine_transfer':'Same-age all-fibroblast weighted 30-neighbor projection; no spatial smoothing, subtype quotas or author labels; weights not calibrated probabilities',
        'technical_recovery':'Opposite-library recovery of joint-fit RNA labels is not independent anatomical validation'})


if __name__=='__main__':main()
