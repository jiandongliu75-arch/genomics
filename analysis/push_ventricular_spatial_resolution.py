"""Map all existing fine CM groups and measure a held-out resolution frontier."""
import argparse
import hashlib
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support

from improve_ventricular_alignment import ROOT,SOURCE,STAGES,normalized,unit_features,neighbors,posteriors,digest

OUT=ROOT/'results/ventricular_spatial_resolution'
PREVIOUS=ROOT/'results/ventricular_alignment_improvement'
QUALITY={'min_precision':.60,'min_recall':.35,'min_f1':.50,'min_query_per_library':10,'min_reference_per_library':15}


def project_probabilities(r,q,labels,classes):
    rf,qf=unit_features(r,q,True)
    ix,w,sim=neighbors(rf,qf)
    available,p=posteriors(labels,ix,w)
    aligned=np.zeros((len(q),len(classes)),np.float32)
    for i,label in enumerate(available):aligned[:,classes.index(label)]=p[:,i]
    return aligned,sim


def scaled(values,mean,sd):
    x=np.clip((values-mean)/np.maximum(sd,.25),-5,5)
    return (x/np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-8)).astype(np.float32)


def partitions(meta):
    result=pd.Series('',index=meta.index)
    for _,part in meta.groupby(['sample','ventricular_group'],observed=True):
        ids=sorted(part.index,key=lambda s:hashlib.sha256(('heart-fine-resolution-v1:'+s).encode()).hexdigest())
        for i,cell in enumerate(ids):result.loc[cell]='calibration' if i%4<2 else ('validation' if i%4==2 else 'assessment')
    return result


def library_profiles(r,meta,classes):
    phase=partitions(meta);truth=meta.ventricular_group.astype(str).to_numpy()
    probability=np.zeros((len(meta),len(classes)),np.float32)
    source=[]
    for library in sorted(meta['sample'].unique()):
        ri=np.flatnonzero(~meta['sample'].eq(library)&phase.eq('calibration'))
        calibration=np.flatnonzero(meta['sample'].eq(library)&phase.eq('calibration'))
        qi=np.flatnonzero(meta['sample'].eq(library))
        rf=scaled(r[ri],r[ri].mean(0),r[ri].std(0))
        qf=scaled(r[qi],r[calibration].mean(0),r[calibration].std(0))
        ix,w,_=neighbors(rf,qf)
        available,p=posteriors(truth[ri],ix,w)
        for i,label in enumerate(available):probability[qi,classes.index(label)]=p[:,i]
        source.append({'query_library':library,'reference_ids':meta.index[ri].tolist(),
            'query_scaler_calibration_ids':meta.index[calibration].tolist(),
            'reference_library':sorted(meta.iloc[ri]['sample'].unique()),'assessment_or_validation_nuclei_in_reference':False})
    return probability,phase,source


def aggregate(p,classes,cut):
    return np.column_stack([p[:,[classes.index(c) for c in members]].sum(1) for members in cut])


def quality_table(p,truth,phase,meta,classes,cut,phase_name):
    mapped={g:i for i,members in enumerate(cut) for g in members}
    y=np.array([mapped[g] for g in truth]);predicted=aggregate(p,classes,cut).argmax(1)
    records=[]
    for library in sorted(meta['sample'].unique()):
        choose=(phase.to_numpy()==phase_name)&meta['sample'].eq(library).to_numpy()
        precision,recall,f1,n=precision_recall_fscore_support(y[choose],predicted[choose],labels=list(range(len(cut))),zero_division=0)
        for j,(a,b,c,d) in enumerate(zip(precision,recall,f1,n)):
            train=~meta['sample'].eq(library)&phase.eq('calibration')&meta.ventricular_group.isin(cut[j])
            records.append({'phase':phase_name,'query_library':library,'classes_at_cut':len(cut),
                'union_label':'+'.join(cut[j]),'precision':a,'recall':b,'f1':c,'query_nuclei':int(d),'reference_nuclei':int(train.sum())})
    frame=pd.DataFrame(records)
    group=frame.groupby('union_label').agg(min_precision=('precision','min'),min_recall=('recall','min'),min_f1=('f1','min'),
        min_query_nuclei=('query_nuclei','min'),min_reference_nuclei=('reference_nuclei','min'))
    group['passes']=(group.min_precision>=QUALITY['min_precision'])&(group.min_recall>=QUALITY['min_recall'])&(group.min_f1>=QUALITY['min_f1'])&(group.min_query_nuclei>=QUALITY['min_query_per_library'])&(group.min_reference_nuclei>=QUALITY['min_reference_per_library'])
    return frame,group


def hierarchy(p,truth,phase,meta,classes):
    cut=[(c,) for c in classes];cuts={};tables=[];frontier=[]
    while True:
        k=len(cut);cuts[k]=[list(m) for m in cut]
        for phase_name in ['calibration','validation','assessment']:
            table,quality=quality_table(p,truth,phase,meta,classes,cut,phase_name)
            tables.append(table)
            frontier.append({'phase':phase_name,'classes_at_cut':k,'passing_unions':int(quality.passes.sum()),
                'all_unions_pass':bool(quality.passes.all()),'minimum_precision':float(quality.min_precision.min()),
                'minimum_recall':float(quality.min_recall.min()),'minimum_f1':float(quality.min_f1.min()),
                'macro_f1':float(table.f1.mean())})
        if k==1:break
        # Merge only from calibration confusion, balanced by true class and
        # averaged across libraries. Validation/assessment labels never choose
        # the pair or the gene scaling/neighbor parameters.
        lookup={c:i for i,m in enumerate(cut) for c in m}
        y=np.array([lookup[c] for c in truth]);yh=aggregate(p,classes,cut).argmax(1)
        conf=np.zeros((k,k))
        for library in sorted(meta['sample'].unique()):
            select=phase.eq('calibration').to_numpy()&meta['sample'].eq(library).to_numpy()
            for i in range(k):
                chosen=select&(y==i)
                if chosen.any():conf[i]+=np.bincount(yh[chosen],minlength=k)/chosen.sum()
        pairs=[(float(conf[i,j]+conf[j,i]),i,j) for i in range(k) for j in range(i+1,k)]
        _,i,j=max(pairs,key=lambda row:(row[0],-row[1],-row[2]))
        merged=tuple(sorted(cut[i]+cut[j]))
        cut=sorted([m for x,m in enumerate(cut) if x not in (i,j)]+[merged])
    frontier=pd.DataFrame(frontier)
    selected=int(frontier.loc[frontier.phase.eq('validation')&frontier.all_unions_pass,'classes_at_cut'].max())
    return cuts,selected,pd.concat(tables,ignore_index=True),frontier


def panel_markers(raw,corrected,meta,genes,classes):
    records=[];selected={}
    exclude=set('Tnnt2 Actc1 Myl2 Myl4 Myl7'.split())
    rawlog=np.log1p(1e4*raw/np.maximum(raw.sum(1,keepdims=True),1))
    corrlog=np.log1p(1e4*corrected/np.maximum(corrected.sum(1,keepdims=True),1))
    for label in classes:
        effects=[];detections=[]
        for library in sorted(meta['sample'].unique()):
            inside=meta['sample'].eq(library)&meta.ventricular_group.eq(label)
            outside=meta['sample'].eq(library)&~meta.ventricular_group.eq(label)
            if not inside.any() or not outside.any():
                effects.append(np.full(len(genes),-np.inf));detections.append(np.zeros(len(genes)))
                continue
            for values in [rawlog,corrlog]:effects.append(values[inside].mean(0)-values[outside].mean(0))
            detections.append((raw[inside]>0).mean(0))
        effect=np.min(effects,axis=0);detection=np.min(detections,axis=0)
        usable=(effect>=.15)&(detection>=.15)&~np.isin(genes,list(exclude))
        order=sorted(np.flatnonzero(usable),key=lambda j:(-effect[j],genes[j]))[:6]
        selected[label]=genes[order].tolist()
        for rank,j in enumerate(order):records.append({'group':label,'gene':genes[j],'rank':rank+1,
            'minimum_raw_corrected_effect_both_libraries':float(effect[j]),'minimum_raw_detection_both_libraries':float(detection[j])})
    return selected,pd.DataFrame(records)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stages',nargs='+',default=STAGES,choices=STAGES);args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
    catalog=pd.read_csv(SOURCE/'reference_group_catalog.tsv',sep='\t').set_index('microenvironment_cluster')
    with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
        raw,corrected,genes,ids=b['original'],b['corrected'],b['genes'].astype(str),b['cell_ids'].astype(str)
    np.testing.assert_array_equal(ids,meta.index)
    design={'scope':'Existing28 stage-specific ventricular RNA groups; no new RNA clustering or changed membership',
        'reference_normalization':'log2(1+CPM),381 shared genes; same assay scaling/30-neighbor/temperature0.10 as preceding benchmark',
        'RNA_split':'Within original library and group, cell IDs hash-sorted with heart-fine-resolution-v1, ranks modulo4:0/1 calibration,2 validation,3 assessment',
        'hierarchy':'Merge maximum bidirectional class-normalized calibration confusion; no assessment labels used for merges or cut selection',
        'selected_cut':'Largest K with every class passing both-library quality criteria on validation quarter; then assess this frozen cut on final quarter',
        'quality_criteria':QUALITY,'half_gene_sets':'Sorted sha256 of heart-fine-panel-v1:gene, alternating halves',
        'panel_marker_selection':'Up to6 genes; minimum group-minus-rest mean log1pCP10k effect>=0.15 in both original libraries and raw/corrected counts; raw detection>=0.15 in both; excludes generic chamber anchors',
        'spatial_tiers':{'multi_check_candidate':'weight>=0.70,margin>=0.20,both libraries agree,original-count reference agrees,>=2 measured supporting genes,original fine class passes final RNA assessment criteria',
            'library_raw_consistent':'weight>=0.50,margin>=0.10,both libraries agree,original reference agrees,>=2 measured supporting genes',
            'exploratory':'weight>=0.35,margin>=0.05','diffuse':'remaining included cells'},
        'markers_independent_validation':False,'author_labels_used_in_fit':False,'spatial_smoothing':False,
        'reference_sha256':digest(SOURCE/'reference_shared_panel_counts.npz')}
    if (OUT/'design.json').exists():assert json.loads((OUT/'design.json').read_text())==design
    else:(OUT/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    for stage in args.stages:
        folder=OUT/stage;folder.mkdir(exist_ok=True)
        if (folder/'completion.json').exists():print('SKIP_COMPLETED',stage,flush=True);continue
        start=time.time();print('START',stage,flush=True)
        with np.load(PREVIOUS/stage/'weighted_reference_neighbors.npz') as b:
            qids,rids,classes,p=b['query_cell_ids'].astype(str),b['reference_cell_ids'].astype(str),b['fine_classes'].astype(str).tolist(),b['fine_weight_fractions']
        rindex=meta.index.get_indexer(rids);assert (rindex>=0).all()
        rm=meta.iloc[rindex];r=corrected[rindex].astype(np.float32);rr=raw[rindex].astype(np.float32)
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        assert pd.Index(qids).isin(data.obs_names).all()
        q=data[qids,genes].X.toarray().astype(np.float32)
        allgenes=np.ones(len(genes),bool)
        rn,qn=normalized(r,allgenes),normalized(q,allgenes)
        libp,phase,library_source=library_profiles(rn,rm,classes)
        cuts,chosen,quality,frontier=hierarchy(libp,rm.ventricular_group.astype(str).to_numpy(),phase,rm,classes)
        quality.to_csv(folder/'RNA_resolution_metrics.tsv',sep='\t',index=False)
        frontier.to_csv(folder/'resolution_frontier.tsv',sep='\t',index=False)
        (folder/'hierarchy.json').write_text(json.dumps({'cuts':cuts,'selected_K_from_validation':chosen,
            'selected_cut':cuts[chosen],'library_sources':library_source},indent=2)+'\n')
        np.savez_compressed(folder/'RNA_split_probabilities.npz',cell_ids=rm.index.to_numpy(dtype=str),classes=np.array(classes),
            probabilities=libp,phase=phase.to_numpy(dtype=str),known_groups=rm.ventricular_group.to_numpy(dtype=str),libraries=rm['sample'].to_numpy(dtype=str))
        print('HIERARCHY',stage,'chosen K',chosen,flush=True)
        rawp,_=project_probabilities(normalized(rr,allgenes),qn,rm.ventricular_group.astype(str).to_numpy(),classes)
        halves=sorted(genes,key=lambda g:hashlib.sha256(('heart-fine-panel-v1:'+g).encode()).hexdigest())
        half_results=[];half_genes=[]
        for groupgenes in [halves[::2],halves[1::2]]:
            use=np.isin(genes,groupgenes)
            keep_r=r[:,use].sum(1)>0;keep_q=q[:,use].sum(1)>0
            hp=np.zeros_like(p)
            hp[keep_q],_=project_probabilities(normalized(r[keep_r],use)[:,use],normalized(q[keep_q],use)[:,use],rm.iloc[np.flatnonzero(keep_r)].ventricular_group.astype(str).to_numpy(),classes)
            half_results.append(hp);half_genes.append(groupgenes)
        np.savez_compressed(folder/'fine_probabilities_and_sensitivity.npz',query_cell_ids=qids,classes=np.array(classes),
            full_panel=p,original_reference=rawp,gene_half_A=half_results[0],gene_half_B=half_results[1],
            genes_half_A=np.array(half_genes[0]),genes_half_B=np.array(half_genes[1]))
        markers,marker_table=panel_markers(rr,r,rm,genes,classes)
        marker_table.to_csv(folder/'measured_panel_supporting_markers.tsv',sep='\t',index=False)
        original_quality=quality[quality.phase.eq('assessment')&quality.classes_at_cut.eq(len(classes))].groupby('union_label').agg(
            min_precision=('precision','min'),min_recall=('recall','min'),min_f1=('f1','min'),min_query_nuclei=('query_nuclei','min'),min_reference_nuclei=('reference_nuclei','min'))
        original_quality['technical_assessment_pass']=(original_quality.min_precision>=.6)&(original_quality.min_recall>=.35)&(original_quality.min_f1>=.5)&(original_quality.min_query_nuclei>=10)&(original_quality.min_reference_nuclei>=15)
        previous=pd.read_csv(PREVIOUS/stage/'weighted_CM_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        obs=previous[['x_centroid','y_centroid','spatial_lineage_candidate','included_in_weighted_mapping','pass_expression_QC','pass_transfer_information']].copy()
        index=pd.Index(qids);best=np.array(classes)[p.argmax(1)];ordered=np.sort(p,axis=1)
        result=pd.DataFrame(index=index)
        result['fine_group_proposal']=best;result['fine_weight']=p.max(1);result['fine_margin']=ordered[:,-1]-ordered[:,-2]
        result['normalized_label_entropy']=-(p*np.log(np.maximum(p,1e-12))).sum(1)/np.log(len(classes))
        result['fine_from_original_reference']=np.array(classes)[rawp.argmax(1)]
        result['fine_from_gene_half_A']=np.where(half_results[0].sum(1)>0,np.array(classes)[half_results[0].argmax(1)],'Insufficient half-panel RNA')
        result['fine_from_gene_half_B']=np.where(half_results[1].sum(1)>0,np.array(classes)[half_results[1].argmax(1)],'Insufficient half-panel RNA')
        result['original_reference_agrees']=result.fine_from_original_reference.eq(best)
        result['both_gene_halves_agree']=result.fine_from_gene_half_A.eq(best)&result.fine_from_gene_half_B.eq(best)
        gcols=[c for c in previous if c.startswith('group_from_')]
        assert len(gcols)==2
        for c in gcols:result[c]=previous.loc[index,c]
        result['both_libraries_agree']=previous.loc[index,gcols].eq(best,axis=0).all(axis=1)
        detected=np.zeros(len(q),int)
        records=[]
        for label in classes:
            selected=best==label;gj=np.flatnonzero(np.isin(genes,markers[label]))
            if len(gj):detected[selected]=(q[selected][:,gj]>0).sum(1)
            record=catalog.loc[label].to_dict();record.update(original_quality.loc[label].to_dict())
            record.update({'stage':stage,'group':label,'proposed_cells':int(selected.sum()),'panel_supporting_genes':';'.join(markers[label]),
                'panel_supporting_gene_count':len(gj),'median_weight':float(np.median(p[selected,classes.index(label)])) if selected.any() else None,
                'raw_reference_agreement_fraction':float(result.loc[selected,'original_reference_agrees'].mean()) if selected.any() else None,
                'both_library_agreement_fraction':float(result.loc[selected,'both_libraries_agree'].mean()) if selected.any() else None,
                'both_gene_halves_agreement_fraction':float(result.loc[selected,'both_gene_halves_agree'].mean()) if selected.any() else None})
            records.append(record)
        result['measured_supporting_genes_detected']=detected
        result['fine_class_passes_RNA_assessment']=result.fine_group_proposal.map(original_quality.technical_assessment_pass)
        exploratory=(result.fine_weight>=.35)&(result.fine_margin>=.05)
        consistent=(result.fine_weight>=.5)&(result.fine_margin>=.1)&result.original_reference_agrees&result.both_libraries_agree&(detected>=2)
        high=(result.fine_weight>=.7)&(result.fine_margin>=.2)&consistent&result.fine_class_passes_RNA_assessment
        result['placement_evidence']=np.select([high,consistent,exploratory],['Multi-check candidate','Library/raw-consistent','Exploratory'],default='Diffuse')
        union_p=aggregate(p,classes,cuts[chosen]);union_labels=np.array(['+'.join(m) for m in cuts[chosen]])
        result['selected_resolution_union']=union_labels[union_p.argmax(1)]
        result['selected_resolution_union_weight']=union_p.max(1)
        for j,label in enumerate(classes):result['weight__'+label]=p[:,j]
        # The published final labels are used solely as a descriptive comparator.
        author_path=ROOT/'results/author_atlas_baseline'/stage/'original_cells_with_author_labels.tsv.gz'
        id_column=pd.read_csv(author_path,sep='\t',nrows=0).columns[0]
        author=pd.read_csv(author_path,sep='\t',index_col=0,usecols=[id_column,'author_atlas_cell_state'])
        obs=obs.join(result).join(author)
        obs.to_csv(folder/'fine_CM_spatial_annotations.tsv.gz',sep='\t',na_rep='NA')
        summary=pd.DataFrame(records)
        for tier in ['Multi-check candidate','Library/raw-consistent','Exploratory','Diffuse']:
            counts=result[result.placement_evidence.eq(tier)].fine_group_proposal.value_counts()
            summary[tier+'_cells']=summary.group.map(counts).fillna(0).astype(int)
        summary.to_csv(folder/'fine_group_evidence_summary.tsv',sep='\t',index=False)
        pd.crosstab(result.fine_group_proposal,obs.loc[index,'author_atlas_cell_state']).to_csv(folder/'fine_groups_vs_published_states.tsv',sep='\t')
        (folder/'completion.json').write_text(json.dumps({'stage':stage,'original_cells':len(data),'query_cells':len(q),'reference_nuclei':len(r),
            'fine_groups':len(classes),'selected_K_from_validation':chosen,'elapsed_seconds':time.time()-start,
            'source_query_sha256':digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad'),
            'source_neighbors_sha256':digest(PREVIOUS/stage/'weighted_reference_neighbors.npz'),
            'author_labels_used_in_fit':False,'original_counts_and_coordinates_changed':False},indent=2)+'\n')
        print('DONE',stage,'tiers',result.placement_evidence.value_counts().to_dict(),flush=True)
        del data,q,r,rr,rn,qn


if __name__=='__main__':main()
