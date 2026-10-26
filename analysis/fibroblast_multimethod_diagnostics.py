"""Published signatures, held-out genes and source-view diagnostics for FBs."""
from pathlib import Path
import json,zipfile,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from fibroblast_multimethod_identity import ROOT,OUT,STAGES,METHODS,FAMILIES,refs,queries,save,read,INPUTS,NAMING

PAPER=Path('/workspace/heart-atlas/metadata/paper_supplementary_tables.xlsx')
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def table(z,index,strings):
    rows=[]
    for row in ET.fromstring(z.read(f'xl/worksheets/sheet{index}.xml')).findall('m:sheetData/m:row',NS):
        d={}
        for cell in row:
            v=cell.find('m:v',NS);value='' if v is None else v.text
            if cell.attrib.get('t')=='s':value=strings[int(value)]
            d[cell.attrib['r'].rstrip('0123456789')]=value
        rows.append(d)
    return rows
def paper_tables(genes,ann):
    with zipfile.ZipFile(PAPER) as z:
        strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        markers=table(z,4,strings)[1:];author=[]
        for index,assay in [(2,'snRNA'),(3,'Xenium')]:
            for d in table(z,index,strings)[1:]:
                if '.Fib' in d.get('B',''):
                    author.append({'stage':d['A'].upper(),'assay':assay,'author_family':d['B'],'author_profiles':int(d['C']),
                                   'assay_total_profiles':int(d['D']),'scope':'Publication Supplementary Tables 1/2; different QC/cohorts from our reviewed counts'})
    author=pd.DataFrame(author);save(author,'published_fibroblast_counts.tsv')
    rows=[]
    for d in markers:
        if '.Fib' not in d.get('A',''):continue
        rows.append({'author_family':d['A'],'gene':d['B'],'log2FC':float(d['D']),'adjusted_p':float(d['E']),
                     'measured_shared_panel':d['B'] in genes})
    markers=pd.DataFrame(rows);save(markers,'published_FB_marker_panel_coverage.tsv')
    coverage=[];panels={}
    for family,g in markers.groupby('author_family'):
        selected=g[g.measured_shared_panel&g.log2FC.gt(0)&g.adjusted_p.lt(.05)].sort_values('log2FC',ascending=False).head(25)
        panels[family]=selected.gene.tolist()
        coverage.append({'author_family':family,'published_genes':len(g),'measured_upregulated_markers':len(selected),
                         'measured_marker_list':';'.join(selected.gene),'module_scored':len(selected)>=3,
                         'spatial_anatomy_established_by_module':False})
    save(pd.DataFrame(coverage),'published_FB_module_coverage.tsv')
    own=ann.groupby(['age','fibroblast_family']).size().rename('reviewed_reference_nuclei').reset_index().rename(columns={'age':'stage'})
    own['author_family']=own.fibroblast_family.str.split('-').str[0]+'.Fib'
    own=own.merge(author[author.assay.eq('snRNA')][['stage','author_family','author_profiles']],how='left',on=['stage','author_family'])
    own['same_barcode_annotation_crosswalk_available']=False
    save(own,'reviewed_vs_published_RNA_family_counts.tsv')
    return {k:v for k,v in panels.items() if len(v)>=3}

def score_modules(counts,genes,panels,reference_mean):
    r=rankdata(counts,axis=1,method='average')/len(genes);result={}
    # Controls are matched using RNA-reference abundance, never query labels.
    bins=pd.qcut(pd.Series(reference_mean).rank(method='first'),10,labels=False).to_numpy()
    rng=np.random.default_rng(261026)
    for family,panel in panels.items():
        ix=pd.Index(genes).get_indexer(panel);weights=np.zeros((len(genes),200))
        for gi in ix:
            choices=np.flatnonzero((bins==bins[gi])&~np.isin(np.arange(len(genes)),ix))
            if len(choices)==0:choices=np.flatnonzero(~np.isin(np.arange(len(genes)),ix))
            selected=rng.choice(choices,200,replace=True)
            np.add.at(weights,(selected,np.arange(200)),1/len(ix))
        observed=r[:,ix].mean(1);null=r@weights
        result[family+'__rank_excess']=observed-null.mean(1)
        result[family+'__control_percentile']=(observed[:,None]>null).mean(1)
        result[family+'__measured_markers']=(counts[:,ix]>0).sum(1)
    return pd.DataFrame(result)

def rna_reference_audit(ann,x,raw,genes):
    rows=[];names=['Gfpt2','Pi16','Adgrd1','Col8a1','Lsamp','Itgbl1','Cmklr1','Htra3','Col6a3','Postn','Tnc','Kirrel3']
    for view,values in [('corrected',x),('raw',raw)]:
        den=np.maximum(values.sum(1),1)
        for (stage,cluster,library),a in ann.groupby(['age','fibroblast_cluster_id','sample']):
            take=ann.index.isin(a.index)
            for gene in names:
                if gene not in genes:continue
                gi=list(genes).index(gene);v=values[take,gi]
                rows.append({'stage':stage,'RNA_subcluster':cluster,'family':a.fibroblast_family.iloc[0],'library':library,
                    'view':view,'nuclei':len(a),'gene':gene,'detection_fraction':float((v>0).mean()),
                    'mean_log1p_per10000_shared_panel':float(np.log1p(v*1e4/den[take]).mean()),
                    'scope':'Within reviewed RNA subcluster; overlapping marker expression is not lineage conversion'})
    save(pd.DataFrame(rows),'RNA_reference_subcluster_marker_audit.tsv')

def main():
    ann,x,raw,genes=refs();panels=paper_tables(genes,ann)
    rna_reference_audit(ann,x,raw,genes)
    module_profiles=[];anchor_profiles=[];cross=[];agreement=[];contexts=[]
    reference_mean=(x/np.maximum(x.sum(1),1)[:,None]).mean(0)
    marker_ix=pd.Index(genes).get_indexer(sorted(NAMING));marker_names=np.array(sorted(NAMING))[marker_ix>=0];marker_ix=marker_ix[marker_ix>=0]
    for stage in STAGES:
        print('FB_MODEL_DIAGNOSTICS',stage,flush=True)
        atlas,q=queries(stage,genes);d=read(OUT/(stage+'_multimethod_FB_annotations.tsv.gz'),index_col='cell_id')
        for view,counts in q.items():
            m=score_modules(counts,genes,panels,reference_mean);m.index=d.index
            m['stage']=stage;m['view']=view;m['cell_id']=d.index
            for family,panel in panels.items():
                # A descriptive enrichment flag, not a family call or p value.
                m[family+'__enrichment_flag']=m[family+'__control_percentile'].ge(.9)&m[family+'__measured_markers'].ge(2)
            m.to_csv(OUT/(stage+'_'+view+'_published_module_scores.tsv.gz'),sep='\t',index=False,na_rep='NA',float_format='%.8g')
            if view=='refined':
                m.to_csv(OUT/(stage+'_published_module_scores.tsv.gz'),sep='\t',index=False,na_rep='NA',float_format='%.8g')
                gf=m['Gfpt2.Fib__enrichment_flag'];co=m['Col8a1.Fib__enrichment_flag']
                for scope,mask in [('all_eligible',d.program_source_eligible),('public_Gfpt2',d.program_source_eligible&d.author_atlas_cell_state.eq('Gfpt2.Fib'))]:
                    for label,flag in [('Gfpt2_only',gf&~co),('Col8a1_only',co&~gf),('both',co&gf),('neither',~co&~gf)]:
                        cross.append({'stage':stage,'scope':scope,'program_combination':label,'profiles':int((mask&flag).sum()),'scope_profiles':int(mask.sum()),
                            'interpretation':'Three-/four-gene published-panel enrichment; not lineage, independent family identity or a transition'})
            labels={'public_label':d.author_atlas_cell_state,'withheld_majority':d['naming_markers_withheld__'+view+'__majority'],
                    'new_cell_candidate':d.recommended_candidate}
            logs=np.log1p(1e4*counts[:,marker_ix]/np.maximum(counts.sum(1),1)[:,None])
            for scope,label in labels.items():
                for identity in sorted(label[d.program_source_eligible].unique()):
                    take=d.program_source_eligible&label.eq(identity)
                    for family in panels:
                        module_profiles.append({'stage':stage,'view':view,'grouping':scope,'identity':identity,'profiles':int(take.sum()),'published_program':family,
                            'mean_rank_excess':m.loc[take,family+'__rank_excess'].mean(),'enriched_fraction':m.loc[take,family+'__enrichment_flag'].mean()})
                    for j,gene in enumerate(marker_names):
                        anchor_profiles.append({'stage':stage,'view':view,'grouping':scope,'identity':identity,'profiles':int(take.sum()),'gene':gene,
                            'detection_fraction':float((counts[take.to_numpy(),marker_ix[j]]>0).mean()),'mean_log1p_per10000_shared_panel':float(logs[take.to_numpy(),j].mean()),
                            'gene_was_withheld_from_named_group_classifier':scope=='withheld_majority'})
        eligible=d.program_source_eligible
        physical=read(ROOT/'results/fibroblast_names_and_neighborhoods'/stage/'FB_physical_neighbors.tsv.gz',index_col=0,
                      usecols=['cell_id','arterial_distance_um','nearest_Endothelial_um','nearest_Cardiomyocyte_um'])
        public=d[eligible&d.author_atlas_cell_state.eq('Gfpt2.Fib')]
        for col in ['shared_panel__refined__majority','shared_panel__refined__reference_supported','late_pooled_reference__refined__majority',
                    'late_pooled_reference__refined__reference_supported','recommended_candidate']:
            if col not in public:continue
            for identity,n in public[col].value_counts().items():agreement.append({'stage':stage,'comparison':col,'identity':identity,'profiles':int(n),'eligible_public_Gfpt2':len(public)})
        for family in FAMILIES:
            selected=eligible&d.recommended_candidate.eq(family)
            if not selected.any():continue
            for measure in ['arterial_distance_um','nearest_Endothelial_um','nearest_Cardiomyocyte_um']:
                z=physical.loc[d.index[selected],measure].dropna()
                assert len(z)>0
                contexts.append({'stage':stage,'family':family,'profiles':int(selected.sum()),'measure':measure,'median_um':z.median(),
                    'within50um_fraction':z.le(50).mean(),'scope':'Descriptive candidate geometry; not vessel identity validation or a matched causal contrast'})
        d['naming_withheld_same_candidate']=d.recommended_candidate.eq(d.naming_markers_withheld__refined__majority)&d.recommended_candidate.ne('Unresolved')
        d['late_pooled_same_candidate']=d.recommended_candidate.eq(d.get('late_pooled_reference__refined__majority',pd.Series('Unavailable',index=d.index)))&d.recommended_candidate.ne('Unresolved')
        d['boundary_supported_candidate']=d.recommended_candidate.where(d.evidence_tier.eq('Cell and core supported'),'Unresolved')
        d['crosschecked_candidate']=d.boundary_supported_candidate.where(d.naming_withheld_same_candidate,'Unresolved')
        d[['x_centroid','y_centroid']]=atlas.loc[d.index,['x_centroid','y_centroid']]
        d.to_csv(OUT/(stage+'_multimethod_FB_annotations.tsv.gz'),sep='\t',na_rep='NA')
    save(pd.DataFrame(module_profiles),'published_module_profiles_by_assignment.tsv')
    save(pd.DataFrame(anchor_profiles),'heldout_naming_marker_profiles.tsv')
    save(pd.DataFrame(cross),'published_Gfpt2_Col8a1_program_overlap.tsv')
    save(pd.DataFrame(agreement),'public_Gfpt2_assignment_comparison.tsv')
    save(pd.DataFrame(contexts),'candidate_fibroblast_anatomical_context.tsv')
    path=OUT/'benchmark_input_paths.json';inputs=set(json.loads(path.read_text()))|INPUTS|{str(PAPER)}
    path.write_text(json.dumps(sorted(inputs),indent=2)+'\n')
    print('FB_MULTIMETHOD_DIAGNOSTICS_COMPLETE',flush=True)
if __name__=='__main__':main()
