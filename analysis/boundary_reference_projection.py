"""Original full-reference projection with frozen full-stage query calibration."""
import anndata as ad
import numpy as np
import pandas as pd

from reassess_whole_heart_lineages import SOURCE,OUT as FIT,PANELS,WITHHOLD,normalized,neighbors,posteriors,top
from build_whole_heart_spatial_atlas import ANCHORS


class FrozenProjection:
    def __init__(self,stage,obs,genes):
        self.stage=stage;self.obs=obs;self.genes=np.array(genes,dtype=str)
        metadata=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
        with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
            refraw,refcorrected,refgenes,allids=b['original'],b['corrected'],b['genes'].astype(str),b['cell_ids'].astype(str)
        np.testing.assert_array_equal(allids,metadata.index)
        with np.load(FIT/stage/'whole_reference_weights.npz') as b:queryids,refids=b['query_cell_ids'].astype(str),b['reference_cell_ids'].astype(str)
        self.refids=refids;self.rmeta=metadata.loc[refids]
        ri=metadata.index.get_indexer(refids);self.training=~np.isin(refgenes,WITHHOLD)
        self.shared_positions=pd.Index(self.genes).get_indexer(refgenes);assert (self.shared_positions>=0).all()
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        q=data[queryids,refgenes].X.toarray().astype(np.float32)
        qn=np.ascontiguousarray(normalized(q,self.training)[:,self.training])
        self.query_mean=qn.mean(0);self.query_sd=np.maximum(qn.std(0),.25)
        del data,q,qn
        ref=normalized(refcorrected[ri],self.training)[:,self.training]
        raw=normalized(refraw[ri],self.training)[:,self.training]
        broad=self.rmeta.reviewed_identity.astype(str).to_numpy()
        self.models={};libraries=sorted(self.rmeta['sample'].unique())
        for name,values,labels in [('primary',ref,broad),('raw_reference',raw,broad)]+[
            (lib,ref[self.rmeta['sample'].eq(lib)],broad[self.rmeta['sample'].eq(lib)]) for lib in libraries]:
            values=np.ascontiguousarray(values)
            features=np.clip((values-values.mean(0))/np.maximum(values.std(0),.25),-5,5)
            features/=np.maximum(np.linalg.norm(features,axis=1,keepdims=True),1e-8)
            self.models[name]=(features.astype(np.float32),labels)
        self.libraries=libraries

    def classify(self,ids,counts):
        ids=np.asarray(ids,dtype=str);obs=self.obs.loc[ids]
        assert counts.shape==(len(ids),len(self.genes))
        shared=counts[:,self.shared_positions].astype(np.float32)
        total=counts.sum(1);detected=(counts>0).sum(1)
        eligible=(total>=30)&(detected>=10)&(shared[:,self.training].sum(1)>=20)&((shared[:,self.training]>0).sum(1)>=8)
        eligible&=obs.cell_area.gt(0).to_numpy()&obs.negative_probe_fraction.le(.1).to_numpy()
        result=pd.DataFrame(index=pd.Index(ids,name='cell_id'))
        result['information_QC_pass']=eligible;result['total_gene_transcripts']=total;result['detected_genes']=detected
        result['baseline_identity']=obs.atlas_identity.to_numpy();result['nuclear_area_missing']=obs.nucleus_area.isna().to_numpy()
        result['new_identity']=np.where(eligible,'Ambiguous/review','QC/information flag')
        result['primary_proposal']='Not eligible';result['lineage_weight']=np.nan;result['lineage_margin']=np.nan
        result['raw_reference_agrees']=False;result['both_RNA_libraries_agree']=False
        predictions={}
        if eligible.any():
            qn=np.ascontiguousarray(normalized(shared[eligible],self.training)[:,self.training])
            features=np.clip((qn-self.query_mean)/self.query_sd,-5,5)
            features/=np.maximum(np.linalg.norm(features,axis=1,keepdims=True),1e-8)
            for model,(ref,labels) in self.models.items():
                ix,w,_=neighbors(ref,features.astype(np.float32));classes,p=posteriors(labels,ix,w)
                if model=='primary':
                    proposal,maxp,margin=top(classes,p)
                    result.loc[ids[eligible],'primary_proposal']=proposal
                    result.loc[ids[eligible],'lineage_weight']=maxp;result.loc[ids[eligible],'lineage_margin']=margin
                else:predictions[model]=np.array(classes)[p.argmax(1)]
            result.loc[ids[eligible],'raw_reference_agrees']=predictions['raw_reference']==proposal
            result.loc[ids[eligible],'both_RNA_libraries_agree']=np.logical_and.reduce([predictions[lib]==proposal for lib in self.libraries])
        panels={};canonical=np.zeros(len(ids),bool)
        for label,genes in PANELS.items():
            col=np.flatnonzero(np.isin(self.genes,genes));values=counts[:,col]
            molecules=values.sum(1);nd=(values>0).sum(1);fraction=molecules/np.maximum(total,1)
            panels[label]=(nd,molecules,fraction)
            support=(nd>=2)&(fraction>=(.02 if label=='Epicardial' else .03))
            if label=='Cardiomyocyte':support=(nd>=3)&(molecules>=10)&(fraction>=.1)
            canonical|=result.primary_proposal.eq(label).to_numpy()&support
        stable=eligible&result.lineage_weight.ge(.7)&result.lineage_margin.ge(.15)&result.raw_reference_agrees&result.both_RNA_libraries_agree&canonical
        result.loc[stable,'new_identity']=result.loc[stable,'primary_proposal']
        other=np.zeros(len(ids),bool)
        for label in ['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune']:
            nd,_,fraction=panels[label];other|=(nd>=3)&(fraction>=.1)
        nd,_,fraction=panels['Cardiomyocyte'];mixed=other&(nd>=3)&(fraction>=.1)
        result['CM_nonCM_marker_coexpression_flag']=mixed;result['CM_marker_fraction']=fraction
        for label,anchor in ANCHORS.items():
            positive=(counts[:,np.flatnonzero(np.isin(self.genes,anchor))]>0).sum(1)>0
            selected=result.new_identity.eq(label).to_numpy()
            result.loc[selected&~positive,'new_identity']='Ambiguous/review'
            result.loc[selected&mixed,'new_identity']='Mixed-marker review'
        result.loc[result.new_identity.eq('Cardiomyocyte')&mixed,'new_identity']='Mixed-marker review'
        return result
