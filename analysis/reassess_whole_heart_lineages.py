"""Whole-heart reference competition, with Ankrd1 outcomes excluded from fits."""
import argparse
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support

from improve_ventricular_alignment import ROOT,SOURCE,STAGES,normalized,unit_features,neighbors,posteriors,digest

OUT=ROOT/'results/whole_heart_Ankrd1_vessels'
WITHHOLD=['Ankrd1','Creb5','Tlr4','Atf3']
PANELS={'Cardiomyocyte':'Tnnt2 Tnni3 Myh6 Myh7 Myl2 Actc1'.split(),
 'Endothelial':'Pecam1 Cdh5 Kdr Esam Tek Vwf Flt1 Tie1 Erg Fli1'.split(),
 'Fibroblast':'Col1a1 Col1a2 Col3a1 Dcn Lum Pdgfra Tcf21'.split(),
 'Pericyte':'Rgs5 Cspg4 Pdgfrb Notch3'.split(),
 'Smooth muscle':'Myh11 Acta2 Tagln Cnn1'.split(),
 'Immune':'Ptprc Csf1r Mertk Mrc1 Fyb Cd22 Cd28 Cd8b1'.split(),
 'Epicardial':'Upk3b Wt1'.split(), 'Glial':'Plp1 Sox10 S100b Mpz'.split()}
ARTERIAL='Gja4 Gja5 Efnb2 Dll4 Hey1 Jag1'.split()


def predict(r,q,broad,leaf,scale=True):
    rf,qf=unit_features(r,q,scale)
    ix,w,similarity=neighbors(rf,qf)
    labels,p=posteriors(broad,ix,w)
    states,sp=posteriors(leaf,ix,w)
    return labels,p,states,sp,similarity,ix,w


def top(labels,p):
    order=np.argsort(p,axis=1)
    return np.array(labels)[order[:,-1]],p.max(1),p[np.arange(len(p)),order[:,-1]]-p[np.arange(len(p)),order[:,-2]]


def marker_counts(data,genes):
    ids=data.var_names.get_indexer([g for g in genes if g in data.var_names])
    if not len(ids):return np.zeros(len(data),int),np.zeros(len(data),int)
    x=data.X[:,ids]
    return np.asarray(x.sum(1)).ravel(),np.asarray((x>0).sum(1)).ravel()


def benchmark(r,rm,stage):
    rows=[]
    for library in sorted(rm['sample'].unique()):
        qi=rm['sample'].eq(library).to_numpy();ri=~qi
        for scale in [False,True]:
            labels,p,states,sp,_,_,_=predict(r[ri],r[qi],rm.iloc[np.flatnonzero(ri)].reviewed_identity.astype(str).to_numpy(),
                rm.iloc[np.flatnonzero(ri)].reference_state.to_numpy(),scale)
            for resolution,truth,classes,scores in [('lineage',rm.iloc[np.flatnonzero(qi)].reviewed_identity.astype(str).to_numpy(),labels,p),
                ('state',rm.iloc[np.flatnonzero(qi)].reference_state.to_numpy(),states,sp)]:
                predicted=np.array(classes)[scores.argmax(1)]
                allclasses=sorted(set(truth)|set(classes))
                pr,re,f1,n=precision_recall_fscore_support(truth,predicted,labels=allclasses,zero_division=0)
                for label,a,b,c,d in zip(allclasses,pr,re,f1,n):
                    rows.append({'stage':stage,'query_library':library,'method':'assay_scaled' if scale else 'Pearson',
                        'resolution':resolution,'label':label,'precision':a,'recall':b,'f1':c,'query_nuclei':int(d)})
    return pd.DataFrame(rows)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stages',nargs='+',default=STAGES,choices=STAGES);args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
    meta['reference_state']=meta.reviewed_identity+'::'+meta.proposed_subtype
    with np.load(SOURCE/'reference_shared_panel_counts.npz') as bundle:
        corrected,raw,genes,ids=bundle['corrected'],bundle['original'],bundle['genes'].astype(str),bundle['cell_ids'].astype(str)
    np.testing.assert_array_equal(ids,meta.index)
    use=~np.isin(genes,WITHHOLD)
    assert use.sum()==377
    design={'scope':'All original Xenium cells and all resolved reference lineages; no supplied spatial labels or subtype quotas',
        'outcomes_excluded_from_all_identity_fits':WITHHOLD,'training_genes':genes[use].tolist(),
        'mapping':'30-neighbor weighted projection,temperature0.10; primary separate-assay gene scaling, SDfloor0.25,clip5; Pearson comparator',
        'eligibility':'Source expression QC plus>=20 shared training molecules and>=8 training genes',
        'resolved_identity':'Weight>=0.70,margin>=0.15,canonical marker support,and agreement with raw-reference and both separate-library projections',
        'canonical_CM_support':'>=3 CM genes,>=10 CM molecules,>=10% of403-target total; no marker override of whole-reference competition',
        'mixed_review':'Resolved CM vote plus>=3 non-CM panel genes and>=10% non-CM molecules; retained as mixed/review',
        'RNA_technical_benchmark':'Opposite-library reference, all377 fitting genes; existing broad labels are not independently defined biological truth',
        'reference_sha256':digest(SOURCE/'reference_shared_panel_counts.npz')}
    if (OUT/'lineage_design.json').exists():assert json.loads((OUT/'lineage_design.json').read_text())==design
    else:(OUT/'lineage_design.json').write_text(json.dumps(design,indent=2)+'\n')
    for stage in args.stages:
        folder=OUT/stage;folder.mkdir(exist_ok=True)
        if (folder/'lineage_completion.json').exists():print('SKIP_COMPLETED',stage,flush=True);continue
        start=time.time();print('LOAD',stage,flush=True)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        allcounts=data[:,genes].X.toarray().astype(np.float32)
        qkeep=data.obs.pass_expression_QC.to_numpy()&(allcounts[:,use].sum(1)>=20)&((allcounts[:,use]>0).sum(1)>=8)
        rkeep=meta.age.eq(stage)&meta.reviewed_identity.ne('Unresolved/mixed')&(corrected[:,use].sum(1)>=20)&((corrected[:,use]>0).sum(1)>=8)&(raw[:,use].sum(1)>=20)&((raw[:,use]>0).sum(1)>=8)
        qi,ri=np.flatnonzero(qkeep),np.flatnonzero(rkeep)
        rm=meta.iloc[ri];rm.to_csv(folder/'reference_pool.tsv.gz',sep='\t',na_rep='NA')
        qn=normalized(allcounts[qi],use)[:,use]
        rn=normalized(corrected[ri],use)[:,use]
        rr=normalized(raw[ri],use)[:,use]
        broad=rm.reviewed_identity.astype(str).to_numpy();leaf=rm.reference_state.to_numpy(dtype=str)
        print('FIT_FULL',stage,len(qi),len(ri),flush=True)
        labels,p,states,sp,similarity,ix,w=predict(rn,qn,broad,leaf,True)
        primary,maxp,margin=top(labels,p)
        np.savez_compressed(folder/'whole_reference_weights.npz',query_cell_ids=data.obs_names[qi].to_numpy(dtype=str),
            reference_cell_ids=rm.index.to_numpy(dtype=str),neighbor_indices=ix,neighbor_weights=w,
            lineage_classes=np.array(labels),lineage_weights=p,state_classes=np.array(states),state_weights=sp)
        result=pd.DataFrame(index=data.obs_names[qi])
        result['whole_reference_lineage_proposal']=primary;result['lineage_weight']=maxp;result['lineage_margin']=margin
        result['nearest_reference_similarity']=similarity
        for key,reference,scaled in [('Pearson',rn,False),('raw_reference',rr,True)]:
            cl,pp,_,_,_,_,_=predict(reference,qn,broad,leaf,scaled)
            result['lineage_from_'+key]=np.array(cl)[pp.argmax(1)]
        for library in sorted(rm['sample'].unique()):
            take=rm['sample'].eq(library).to_numpy()
            cl,pp,_,_,_,_,_=predict(rn[take],qn,broad[take],leaf[take],True)
            result['lineage_from_'+library]=np.array(cl)[pp.argmax(1)]
        libraries=[c for c in result if c.startswith('lineage_from_GSM')]
        assert len(libraries)==2
        result['both_RNA_libraries_agree']=result[libraries].eq(primary,axis=0).all(axis=1)
        result['raw_reference_agrees']=result.lineage_from_raw_reference.eq(primary)
        result['Pearson_method_agrees']=result.lineage_from_Pearson.eq(primary)
        result['EC_arterial_conditional_weight']=np.divide(sp[:,states.index('Endothelial::Arterial-like')],p[:,labels.index('Endothelial')],out=np.zeros(len(qi)),where=p[:,labels.index('Endothelial')]>0)
        for label in ['Endocardial-like','Capillary-like','Venous/venular-like']:
            key='Endothelial::'+label
            result['EC_'+label+'_conditional_weight']=np.divide(sp[:,states.index(key)],p[:,labels.index('Endothelial')],out=np.zeros(len(qi)),where=p[:,labels.index('Endothelial')]>0) if key in states else np.zeros(len(qi))
        for j,label in enumerate(labels):result['lineage_weight__'+label]=p[:,j]
        obs=data.obs[['x_centroid','y_centroid','cell_area','nucleus_area','transcript_counts','negative_probe_fraction','pass_expression_QC','spatial_lineage_candidate']].copy()
        obs['eligible_whole_reference']=qkeep
        obs=obs.join(result)
        totals=np.asarray(data.X.sum(1)).ravel()
        for label,panel in PANELS.items():
            molecules,detected=marker_counts(data,panel)
            key=label.replace(' ','_')
            obs['marker_'+key+'_molecules']=molecules;obs['marker_'+key+'_genes']=detected
            obs['marker_'+key+'_fraction']=molecules/np.maximum(totals,1)
        canonical=np.zeros(len(data),bool)
        for label in labels:
            key=label.replace(' ','_')
            minimum=.02 if label=='Epicardial' else .03
            support=obs['marker_'+key+'_genes'].ge(2)&obs['marker_'+key+'_fraction'].ge(minimum)
            if label=='Cardiomyocyte':support=obs.marker_Cardiomyocyte_genes.ge(3)&obs.marker_Cardiomyocyte_molecules.ge(10)&obs.marker_Cardiomyocyte_fraction.ge(.10)
            canonical|=(obs.whole_reference_lineage_proposal.eq(label)&support).to_numpy()
        stable=qkeep&obs.lineage_weight.ge(.7)&obs.lineage_margin.ge(.15)&obs.both_RNA_libraries_agree.eq(True)&obs.raw_reference_agrees.eq(True)&canonical
        other_panels=['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune']
        mixed=np.zeros(len(data),bool)
        for label in other_panels:
            key=label.replace(' ','_')
            mixed|=(obs['marker_'+key+'_genes'].ge(3)&obs['marker_'+key+'_fraction'].ge(.10)).to_numpy()
        mixed&=obs.marker_Cardiomyocyte_genes.ge(3)&obs.marker_Cardiomyocyte_fraction.ge(.10)
        obs['mixed_CM_nonCM_marker_review']=mixed
        obs['reassessed_lineage']='Ambiguous/review'
        obs.loc[~qkeep,'reassessed_lineage']='QC/information flag'
        obs.loc[stable,'reassessed_lineage']=obs.loc[stable,'whole_reference_lineage_proposal']
        obs.loc[stable&mixed&obs.whole_reference_lineage_proposal.eq('Cardiomyocyte'),'reassessed_lineage']='Mixed CM/non-CM; review'
        obs['CM_strict']=obs.reassessed_lineage.eq('Cardiomyocyte')
        obs['CM_inclusive']=qkeep&obs.whole_reference_lineage_proposal.eq('Cardiomyocyte')&obs.lineage_weight.ge(.5)&obs.marker_Cardiomyocyte_genes.ge(2)&obs.marker_Cardiomyocyte_molecules.ge(10)
        for gene in ['Myl2','Myl4','Npr3','Ankrd1','Creb5','Tlr4','Atf3','Myh11','Acta2','Tagln','Cnn1']:
            obs['measured_'+gene]=data[:,[gene]].X.toarray().ravel() if gene in data.var_names else 0
        obs['ventricular_chamber_marker_support']=obs.measured_Myl2.ge(2)&obs.measured_Myl2.ge(obs.measured_Myl4)
        obs['ventricular_CM_strict']=obs.CM_strict&obs.ventricular_chamber_marker_support
        obs['ventricular_CM_inclusive']=obs.CM_inclusive&obs.ventricular_chamber_marker_support
        molecules,detected=marker_counts(data,ARTERIAL)
        obs['arterial_panel_molecules']=molecules;obs['arterial_panel_genes']=detected
        molecules,detected=marker_counts(data,['Gja4','Efnb2','Dll4'])
        obs['arterial_specific_genes']=detected
        obs['arterial_EC_seed']=obs.reassessed_lineage.eq('Endothelial')&obs.arterial_panel_genes.ge(2)&obs.arterial_specific_genes.ge(1)&obs.EC_arterial_conditional_weight.ge(.3)
        obs['arterial_EC_marker_seed']=qkeep&obs.whole_reference_lineage_proposal.eq('Endothelial')&obs.lineage_weight.ge(.5)&obs.marker_Endothelial_genes.ge(3)&obs.arterial_panel_genes.ge(2)&obs.arterial_specific_genes.ge(1)
        obs['SMC_seed']=obs.reassessed_lineage.eq('Smooth muscle')&obs.marker_Smooth_muscle_genes.ge(2)&(obs.measured_Myh11.gt(0)|obs.measured_Cnn1.gt(0))
        obs['SMC_marker_seed']=qkeep&obs.marker_Smooth_muscle_genes.ge(2)&(obs.measured_Myh11.gt(0)|obs.measured_Cnn1.gt(0))&obs.marker_Smooth_muscle_fraction.ge(.1)&~obs.CM_strict
        targetcols=[data.var_names.get_loc(g) for g in WITHHOLD if g in data.var_names]
        denominator=totals-np.asarray(data.X[:,targetcols].sum(1)).ravel()
        obs['outcome_excluded_transcripts']=denominator
        obs['training_detected_genes']=0;obs.loc[data.obs_names[qi],'training_detected_genes']=(allcounts[qi][:,use]>0).sum(1)
        obs['Ankrd1_log1pCP10k']=np.log1p(1e4*obs.measured_Ankrd1/np.maximum(denominator,1))
        module=np.column_stack([np.log1p(1e4*obs['measured_'+g]/np.maximum(denominator,1)) for g in ['Creb5','Tlr4','Atf3']])
        obs['response_module_log1pCP10k_mean']=module.mean(1)
        authorpath=ROOT/'results/author_atlas_baseline'/stage/'original_cells_with_author_labels.tsv.gz'
        idcol=pd.read_csv(authorpath,sep='\t',nrows=0).columns[0]
        author=pd.read_csv(authorpath,sep='\t',usecols=[idcol,'author_atlas_cell_state'],index_col=0)
        np.testing.assert_array_equal(obs.index,author.index)
        obs=obs.join(author)
        obs['author_ventricular']=obs.author_atlas_cell_state.isin(['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM'])|(obs.author_atlas_cell_state.eq('p.CM')&obs.ventricular_chamber_marker_support)
        obs['CM_consensus_ventricular']=obs.ventricular_CM_strict&obs.author_ventricular
        obs.to_csv(folder/'whole_heart_cell_annotations.tsv.gz',sep='\t',na_rep='NA')
        obs.groupby('reassessed_lineage').size().rename('cells').to_csv(folder/'lineage_counts.tsv',sep='\t')
        pd.crosstab(obs.reassessed_lineage,obs.author_atlas_cell_state).to_csv(folder/'reassessed_vs_published_lineages.tsv',sep='\t')
        print('BENCHMARK',stage,flush=True)
        benchmark(rn,rm,stage).to_csv(folder/'cross_library_technical_metrics.tsv',sep='\t',index=False)
        (folder/'lineage_completion.json').write_text(json.dumps({'stage':stage,'original_cells':len(data),'query_cells':len(qi),'reference_nuclei':len(ri),
            'resolved_lineage_cells':int(stable.sum()),'strict_CM_cells':int(obs.CM_strict.sum()),'strict_ventricular_cells':int(obs.ventricular_CM_strict.sum()),
            'arterial_EC_seeds':int(obs.arterial_EC_seed.sum()),'SMC_seeds':int(obs.SMC_seed.sum()),
            'arterial_marker_seeds':int(obs.arterial_EC_marker_seed.sum()),'SMC_marker_seeds':int(obs.SMC_marker_seed.sum()),
            'elapsed_seconds':time.time()-start,'source_query_sha256':digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad'),
            'outcome_genes_in_fit':False,'coordinates_or_counts_changed':False},indent=2)+'\n')
        print('DONE',stage,json.loads((folder/'lineage_completion.json').read_text()),flush=True)
        del data,allcounts,qn,rn,rr


if __name__=='__main__':main()
