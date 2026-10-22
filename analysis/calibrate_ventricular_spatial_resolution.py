"""Check local fine-group thresholds and audit resolution outputs."""
import json

import anndata as ad
import numpy as np
import pandas as pd

from push_ventricular_spatial_resolution import (OUT,PREVIOUS,SOURCE,ROOT,STAGES,partitions,
    normalized,unit_features,aggregate,quality_table,QUALITY,digest)

THRESHOLDS=[.35,.50,.60,.70,.80,.90]


def selective(p,truth,phase,libraries,classes,label,threshold,split):
    order=np.sort(p,axis=1);margin=order[:,-1]-order[:,-2]
    predicted=np.array(classes)[p.argmax(1)]
    rows=[]
    for library in sorted(set(libraries)):
        called=(phase==split)&(libraries==library)&(predicted==label)&(p.max(1)>=threshold)&(margin>=.10)
        n=int(called.sum());correct=int((truth[called]==label).sum())
        rows.append({'library':library,'called':n,'correct':correct,'precision':correct/n if n else np.nan})
    return rows


def passes(rows):
    return all(x['called']>=10 and x['precision']>=.75 for x in rows)


def main():
    design={'threshold_grid':THRESHOLDS,'minimum_margin':.10,'minimum_selected_per_library':10,
        'minimum_validation_precision_each_library':.75,'minimum_assessment_precision_each_library':.75,
        'selection':'Lowest threshold passing validation in both libraries; assessment evaluates that frozen threshold.',
        'spatial_local_flag':'Selected threshold and margin, both full-library projections agree, original-reference projection agrees,>=2 measured supporting genes, class passes selective RNA assessment.',
        'interpretation':'Internal selective RNA label recovery; neither a calibrated Xenium probability nor independent anatomical truth.'}
    if (OUT/'local_specificity_design.json').exists():assert json.loads((OUT/'local_specificity_design.json').read_text())==design
    else:(OUT/'local_specificity_design.json').write_text(json.dumps(design,indent=2)+'\n')
    meta=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
    with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
        raw,corr,genes,ids=b['original'],b['corrected'],b['genes'].astype(str),b['cell_ids'].astype(str)
    np.testing.assert_array_equal(ids,meta.index)
    summaries,frontiers,local_tables,audits=[],[],[],[]
    for stage in STAGES:
        folder=OUT/stage
        completion=json.loads((folder/'completion.json').read_text())
        with np.load(PREVIOUS/stage/'weighted_reference_neighbors.npz') as b:
            qids,rids=b['query_cell_ids'].astype(str),b['reference_cell_ids'].astype(str)
            original_classes=b['fine_classes'].astype(str).tolist();oldp=b['fine_weight_fractions']
        rm=meta.loc[rids];ri=meta.index.get_indexer(rids)
        with np.load(folder/'fine_probabilities_and_sensitivity.npz') as b:
            np.testing.assert_array_equal(b['query_cell_ids'],qids);classes=b['classes'].astype(str).tolist()
            assert classes==original_classes
            p,rawp,ha,hb=b['full_panel'],b['original_reference'],b['gene_half_A'],b['gene_half_B']
            ga,gb=b['genes_half_A'].astype(str),b['genes_half_B'].astype(str)
        np.testing.assert_array_equal(p,oldp)
        assert set(ga).isdisjoint(gb) and set(ga)|set(gb)==set(genes)
        for values in [p,rawp]:np.testing.assert_allclose(values.sum(1),1,atol=3e-7)
        for values in [ha,hb]:
            mass=values.sum(1);assert (np.isclose(mass,1,atol=3e-7)|np.isclose(mass,0)).all()
        with np.load(folder/'RNA_split_probabilities.npz') as b:
            np.testing.assert_array_equal(b['cell_ids'],rids)
            rp,phase,truth,libraries=b['probabilities'],b['phase'].astype(str),b['known_groups'].astype(str),b['libraries'].astype(str)
        np.testing.assert_array_equal(partitions(rm).to_numpy(),phase)
        np.testing.assert_array_equal(truth,rm.ventricular_group)
        np.testing.assert_array_equal(libraries,rm['sample'])
        np.testing.assert_allclose(rp.sum(1),1,atol=3e-7)
        hierarchy=json.loads((folder/'hierarchy.json').read_text())
        for provenance in hierarchy['library_sources']:
            train=meta.loc[provenance['reference_ids']]
            assert not train['sample'].eq(provenance['query_library']).any()
            assert set(train.index).issubset(set(rm.index[phase=='calibration']))
            assert set(provenance['query_scaler_calibration_ids']).issubset(set(rm.index[phase=='calibration']))
        assert all(set(sum(cut,[]))==set(classes) and len(sum(cut,[]))==len(classes) for cut in hierarchy['cuts'].values())
        frontier=pd.read_csv(folder/'resolution_frontier.tsv',sep='\t')
        expected_selected=int(frontier.loc[frontier.phase.eq('validation')&frontier.all_unions_pass,'classes_at_cut'].max())
        assert expected_selected==hierarchy['selected_K_from_validation']
        for k,cut in hierarchy['cuts'].items():
            for split in ['validation','assessment']:
                _,quality=quality_table(rp,truth,pd.Series(phase,index=rm.index),rm,classes,cut,split)
                row=frontier[frontier.phase.eq(split)&frontier.classes_at_cut.eq(int(k))].iloc[0]
                assert quality.passes.sum()==row.passing_unions and quality.passes.all()==row.all_unions_pass
        records=[]
        for label in classes:
            threshold=None;validation=[]
            for candidate in THRESHOLDS:
                values=selective(rp,truth,phase,libraries,classes,label,candidate,'validation')
                if passes(values):threshold=candidate;validation=values;break
            assessment=selective(rp,truth,phase,libraries,classes,label,threshold if threshold is not None else np.inf,'assessment')
            record={'stage':stage,'group':label,'selected_weight_threshold':threshold,'minimum_margin':.1,
                'validation_rule_passes':threshold is not None,'selective_RNA_assessment_pass':threshold is not None and passes(assessment),
                'validation_min_selected':min(x['called'] for x in validation) if validation else 0,
                'validation_min_precision':min(x['precision'] for x in validation) if validation else np.nan,
                'assessment_min_selected':min(x['called'] for x in assessment),
                'assessment_min_precision':min(x['precision'] for x in assessment) if all(x['called'] for x in assessment) else np.nan}
            for i,values in enumerate(assessment):
                for key in ['library','called','correct','precision']:record[f'assessment_fold{i+1}_{key}']=values[key]
            records.append(record)
        local=pd.DataFrame(records).set_index('group')
        local.to_csv(folder/'class_specific_RNA_thresholds.tsv',sep='\t',na_rep='NA')
        obs=pd.read_csv(folder/'fine_CM_spatial_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        data=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        np.testing.assert_array_equal(obs.index,data.obs_names)
        np.testing.assert_array_equal(obs[['x_centroid','y_centroid']],data.obs[['x_centroid','y_centroid']])
        assert set(obs.index[obs.included_in_weighted_mapping])==set(qids)
        result=obs.loc[qids]
        np.testing.assert_array_equal(result.fine_group_proposal,np.array(classes)[p.argmax(1)])
        np.testing.assert_allclose(result.fine_weight,p.max(1),rtol=1e-6)
        np.testing.assert_array_equal(result.fine_from_original_reference,np.array(classes)[rawp.argmax(1)])
        for values,field in [(ha,'fine_from_gene_half_A'),(hb,'fine_from_gene_half_B')]:
            expected=np.where(values.sum(1)>0,np.array(classes)[values.argmax(1)],'Insufficient half-panel RNA')
            np.testing.assert_array_equal(result[field],expected)
        q=data[qids,genes].X.toarray().astype(np.float32)
        r=corr[ri].astype(np.float32);rr=raw[ri].astype(np.float32)
        markers=pd.read_csv(folder/'measured_panel_supporting_markers.tsv',sep='\t')
        detections=np.zeros(len(q),int)
        for label in classes:
            subset=result.fine_group_proposal.eq(label).to_numpy();g=markers.loc[markers.group.eq(label),'gene']
            gi=np.flatnonzero(np.isin(genes,g))
            if len(gi):detections[subset]=(q[subset][:,gi]>0).sum(1)
        np.testing.assert_array_equal(detections,result.measured_supporting_genes_detected)
        thresholds=result.fine_group_proposal.map(local.selected_weight_threshold).fillna(np.inf)
        supported=result.fine_group_proposal.map(local.selective_RNA_assessment_pass)
        flag=supported&result.fine_weight.ge(thresholds)&result.fine_margin.ge(.1)&result.both_libraries_agree&result.original_reference_agrees&result.measured_supporting_genes_detected.ge(2)
        obs['local_RNA_threshold_candidate']=False
        obs.loc[qids,'local_RNA_threshold_candidate']=flag.to_numpy()
        obs['prioritized_fine_candidate']=False
        obs.loc[qids,'prioritized_fine_candidate']=(flag|result.placement_evidence.eq('Multi-check candidate')).to_numpy()
        obs.loc[qids,'local_RNA_threshold']=thresholds.replace(np.inf,np.nan).to_numpy()
        obs.loc[qids,'selective_RNA_assessment_pass']=supported.to_numpy()
        obs.to_csv(folder/'fine_CM_spatial_annotations.tsv.gz',sep='\t',na_rep='NA')
        group_summary=pd.read_csv(folder/'fine_group_evidence_summary.tsv',sep='\t').set_index('group')
        group_summary=group_summary.drop(columns=[c for c in local if c in group_summary and c!='stage'],errors='ignore').join(local.drop(columns='stage'))
        group_summary['stage']=stage
        group_summary['local_RNA_threshold_candidate_cells']=result.loc[flag,'fine_group_proposal'].value_counts().reindex(group_summary.index,fill_value=0)
        group_summary['prioritized_fine_candidate_cells']=result.loc[flag|result.placement_evidence.eq('Multi-check candidate'),'fine_group_proposal'].value_counts().reindex(group_summary.index,fill_value=0)
        group_summary.to_csv(folder/'fine_group_evidence_summary.tsv',sep='\t')
        # Independently fully sort sampled projection scores for all sensitivity
        # views, checking posterior sums rather than just top labels.
        sample=np.random.default_rng(982).choice(len(q),30,replace=False)
        checks=[]
        for reference,use,actual in [(rr,np.ones(len(genes),bool),rawp),(r,np.isin(genes,ga),ha),(r,np.isin(genes,gb),hb)]:
            kr=reference[:,use].sum(1)>0;kq=q[:,use].sum(1)>0
            rf,qf=unit_features(normalized(reference[kr],use)[:,use],normalized(q[kq],use)[:,use],True)
            sampled=sample[kq[sample]];qposition=np.flatnonzero(kq).searchsorted(sampled)
            score=qf[qposition]@rf.T;ix=np.argsort(-score,axis=1,kind='stable')[:,:30]
            values=np.take_along_axis(score,ix,axis=1);w=np.exp((values-values[:,:1])/.1);w/=w.sum(1,keepdims=True)
            y=rm.iloc[np.flatnonzero(kr)].ventricular_group.to_numpy(dtype=str)
            expected=np.column_stack([(w*(y[ix]==label)).sum(1) for label in classes])
            np.testing.assert_allclose(actual[sampled],expected,atol=3e-5,rtol=3e-5)
            checks.append(len(sampled))
        assert digest(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad')==completion['source_query_sha256']
        assert digest(PREVIOUS/stage/'weighted_reference_neighbors.npz')==completion['source_neighbors_sha256']
        chosen=frontier[frontier.phase.eq('assessment')&frontier.classes_at_cut.eq(expected_selected)].iloc[0]
        summaries.append({'stage':stage,'original_cells':len(obs),'mapped_ventricular_cells':len(qids),'fine_RNA_groups_tested':len(classes),
            'validation_selected_stagewide_groups':expected_selected,'selected_stagewide_cut_passes_assessment':bool(chosen.all_unions_pass),
            'fine_classes_with_selective_RNA_support':int(local.selective_RNA_assessment_pass.sum()),
            'fine_groups_with_local_candidate_cells':int(result.loc[flag,'fine_group_proposal'].nunique()),'local_candidate_cells':int(flag.sum()),
            'prioritized_fine_groups':int(result.loc[flag|result.placement_evidence.eq('Multi-check candidate'),'fine_group_proposal'].nunique()),
            'prioritized_fine_candidate_cells':int((flag|result.placement_evidence.eq('Multi-check candidate')).sum()),
            'multi_check_candidate_cells':int(result.placement_evidence.eq('Multi-check candidate').sum()),
            'library_raw_consistent_cells':int(result.placement_evidence.eq('Library/raw-consistent').sum()),
            'both_gene_halves_agree_cells':int(result.both_gene_halves_agree.sum())})
        frontiers.append(frontier.assign(stage=stage));local_tables.append(group_summary.reset_index())
        audits.append({'stage':stage,'original_rows_conserved':len(obs),'fine_proposals_checked':len(qids),
            'calibration_validation_assessment_rows_partitioned':len(rm),'no_assessment_labels_used_to_select_merges_or_cut':True,
            'sampled_sensitivity_rows_independently_recomputed':sum(checks),'all_marker_detection_counts_recomputed':len(qids),
            'source_hashes_unchanged':True})
        print('VALIDATED',stage,summaries[-1],flush=True)
        del data,q,r,rr
    pd.DataFrame(summaries).to_csv(OUT/'stage_resolution_summary.tsv',sep='\t',index=False)
    pd.concat(frontiers,ignore_index=True).to_csv(OUT/'resolution_frontier.tsv',sep='\t',index=False)
    pd.concat(local_tables,ignore_index=True).to_csv(OUT/'fine_group_evidence_summary.tsv',sep='\t',index=False)
    (OUT/'validation.json').write_text(json.dumps({'status':'PASS','stages':audits,
        'scope':'Numerical conservation, technical label recovery and sensitivity; no independent anatomical or animal-level validation.'},indent=2)+'\n')


if __name__=='__main__':main()
