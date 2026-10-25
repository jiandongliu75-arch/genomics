"""Primary EC presentation at the five requested vascular identity categories.

Presentation-only consolidation: source QC, identities, confidence, graphs,
expression counts and physical/UMAP coordinates remain unchanged.
"""
from endothelial_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import subprocess
import re

SOURCE=OUT
DEST=ROOT/'results/endothelial_five_types'
TYPES=['Capillary-like','Endocardial-like','Arterial-like','Venous-like','Lymphatic-like']
LABELS=dict(zip(FAMILIES+[UNRESOLVED],TYPES+['Unresolved']))
PALETTE={LABELS[f]:c for f,c in COLORS.items()}

def relabel(series):
    assert series.isin(LABELS).all(),set(series)-set(LABELS)
    return series.map(LABELS)

def legend(ax,labels,where='outside'):
    handles=[Line2D([],[],ls='',marker='o',color=c,markersize=5,label=f'{t} (n={int(labels.eq(t).sum()):,})') for t,c in PALETTE.items()]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1,1),fontsize=9,frameon=False)

def scatter(ax,xy,labels,size=5):
    # Review/unresolved profiles form the background rather than hide accepted labels.
    for t in ['Unresolved']+TYPES:
        m=labels.eq(t).to_numpy()
        ax.scatter(xy[m,0],xy[m,1],s=size,color=PALETTE[t],alpha=.75,edgecolors='none',rasterized=True)

def finish(fig,name,pdf):
    fig.savefig(DEST/(name+'.png'),dpi=180,bbox_inches='tight')
    pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def main():
    DEST.mkdir(exist_ok=True)
    ann=pd.read_csv(SOURCE/'snRNA/EC_annotations.tsv.gz',sep='\t',index_col=0)
    qc=pd.read_csv(SOURCE/'snRNA/all_endothelial_QC_decisions.tsv.gz',sep='\t',index_col=0)
    joint=pd.read_csv(SOURCE/'snRNA/all_stage_unintegrated_UMAP.tsv.gz',sep='\t',index_col=0)
    spatial=pd.read_csv(SOURCE/'Xenium/all_stage_EC_consensus.tsv.gz',sep='\t',index_col=0)
    rna=qc[['age','sample','primary_EC_eligible','QC_disposition','cycling_marker_count','cycling_state']].copy()
    rna['EC_type']='Unresolved'
    rna.loc[ann.index,'EC_type']=relabel(ann.EC_family)
    rna['identity_status']='QC review / unassigned'
    rna.loc[ann.index,'identity_status']=np.where(ann.EC_family.eq(UNRESOLVED),'Vascular identity unresolved',
        np.where(ann.EC_family.isin(FAMILIES[3:]),'Provisional candidate identity','RNA-based identity'))
    rna['cycling_flag']=rna.cycling_marker_count.ge(2)
    rna['source_group_support']=ann.EC_grouping_support.reindex(rna.index)
    rna['source_group_audit_only']=ann.EC_cluster_id.reindex(rna.index)
    for c in ['EC_UMAP1','EC_UMAP2']:rna[c]=ann[c].reindex(rna.index)
    rna['joint_UMAP1']=joint.UMAP1.reindex(rna.index);rna['joint_UMAP2']=joint.UMAP2.reindex(rna.index)
    rna.to_csv(DEST/'snRNA_five_type_annotations.tsv.gz',sep='\t',na_rep='NA')
    x=spatial[['stage','x_centroid','y_centroid','atlas_identity','native_EC_broad','confidence_tier',
        'high_confidence_subtype','assigned_including_moderate','nuclear_only_subtype_candidate',
        'nuclear_rescue_with_mixed_cell_RNA','multiple_projected_nuclear_components_review','geometry_unresolved_review']].copy()
    x['EC_type']=relabel(spatial.fine_consensus_family).where(spatial.high_confidence_subtype,'Unresolved')
    x['EC_type_including_cell_supported']=relabel(spatial.fine_consensus_family)
    x['nuclear_supported_EC_type']=relabel(spatial.nuclear_supported_family)
    x['additional_nuclear_only_candidate_type']=x.nuclear_supported_EC_type.where(x.nuclear_only_subtype_candidate,'Unresolved')
    x.insert(0,'stage_cell_id',[s+'::'+c for s,c in zip(x.stage,x.index)])
    x['EC_type_is_provisional']=x.EC_type.isin(TYPES[3:])
    x['nuclear_type_is_provisional']=x.nuclear_supported_EC_type.isin(TYPES[3:])
    x['identity_status']=np.where(x.EC_type.eq('Unresolved'),'Unresolved under strict cell+nuclear criteria',
        np.where(x.EC_type.isin(TYPES[3:]),'Provisional candidate identity; strict technical support','Strict cell+nuclear technical support'))
    x.to_csv(DEST/'Xenium_five_type_annotations.tsv.gz',sep='\t',na_rep='NA')
    key=pd.read_csv(SOURCE/'endothelial_subtype_name_key.tsv',sep='\t')
    key['subtype_name']=relabel(key.subtype_name)
    key['interpretation_status']=['RNA-based identity','RNA-based identity','RNA-based identity; stage-dependent support','Provisional candidate','Provisional candidate']
    save(key,DEST/'five_type_marker_key.tsv')
    rows=[]
    for stage in STAGES:
        for t in TYPES+['Unresolved']:
            r=rna[rna.age.eq(stage)&rna.primary_EC_eligible]
            rows.append({'assay':'snRNA','scope':'primary_QC_endothelial','stage':stage,'EC_type':t,'cells':int(r.EC_type.eq(t).sum())})
            s=x[x.stage.eq(stage)]
            for scope,m in [('previously_labeled_EC',s.native_EC_broad),('all_EC_candidates',pd.Series(True,index=s.index))]:
                for suffix,col in [('strict_cell_nuclear','EC_type'),('including_cell_supported','EC_type_including_cell_supported')]:
                    rows.append({'assay':'Xenium','scope':scope+'__'+suffix,'stage':stage,'EC_type':t,'cells':int(s.loc[m,col].eq(t).sum())})
                # This is a separate, disjoint candidate layer, not part of the strict classification.
                if t!='Unresolved':
                    rows.append({'assay':'Xenium','scope':scope+'__additional_nuclear_only_candidates','stage':stage,'EC_type':t,
                        'cells':int(s.loc[m,'additional_nuclear_only_candidate_type'].eq(t).sum())})
    counts=pd.DataFrame(rows);save(counts,DEST/'five_type_counts_by_stage.tsv')
    rna_counts=counts[counts.assay.eq('snRNA')].pivot(index='EC_type',columns='stage',values='cells').reindex(TYPES+['Unresolved'])[STAGES]
    save(rna_counts.reset_index(),DEST/'snRNA_five_type_stage_counts.tsv')
    # Verify the requested level of presentation preserves every source decision.
    assert len(rna)==5573 and int(rna.primary_EC_eligible.sum())==4597
    np.testing.assert_array_equal(rna.loc[ann.index,'EC_type'],relabel(ann.EC_family))
    assert x.high_confidence_subtype.equals(spatial.high_confidence_subtype)
    assert x.assigned_including_moderate.equals(spatial.assigned_including_moderate)
    np.testing.assert_array_equal(x[['x_centroid','y_centroid']],spatial[['x_centroid','y_centroid']])
    assert len(x)==95316 and x.stage_cell_id.is_unique and x.EC_type.isin(TYPES+['Unresolved']).all()
    assert not x.loc[~x.high_confidence_subtype,'EC_type'].ne('Unresolved').any()
    assert rna_counts.sum().sum()==4597
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
    with PdfPages(DEST/'five_endothelial_types_report.pdf') as pdf:
        fig=plt.figure(figsize=(12,8));fig.text(.06,.92,'Five endothelial types · primary classification',fontsize=20,weight='bold')
        lines=['Capillary-like · Endocardial-like · Arterial-like · Venous-like · Lymphatic-like',
            'The five vascular identity categories are the primary labels in this report and its tables.',
            'Venous-like and lymphatic-like interpretations remain provisional; confidence is a separate field.',
            'Cycling is a separate RNA annotation, not a sixth vascular type.',
            'Unresolved means insufficient evidence for a vascular subtype; it is not a sixth biological type.',
            'Prior QC, reference mapping, expression counts, UMAP coordinates and physical positions are preserved.',
            'RNA: 4,597 primary endothelial nuclei across P0, P7, P14 and P21; 976 QC-review nuclei retained.',
            'Xenium: 69,617 previously labeled EC profiles plus uncertain EC proposals; 95,316 candidates total.',
            'Spatial panels separate strict cell+nuclear agreement from nuclear-only evidence.',
            'Zero supported calls at a stage do not establish biological absence.']
        for i,line in enumerate(lines):fig.text(.06,.82-i*.064,line,fontsize=11)
        pdf.savefig(fig);plt.close(fig)
        a=rna.loc[joint.index];xy=joint[['UMAP1','UMAP2']].to_numpy()
        fig,axes=plt.subplots(1,2,figsize=(17,7),layout='constrained')
        scatter(axes[0],xy,a.EC_type);legend(axes[0],a.EC_type);axes[0].set_title('Five EC types')
        for stage,c in zip(STAGES,['#0072b2','#e69f00','#009e73','#cc3311']):
            m=a.age.eq(stage).to_numpy();axes[1].scatter(xy[m,0],xy[m,1],s=5,color=c,label=stage,rasterized=True)
        axes[1].legend(frameon=False);axes[1].set_title('Developmental stage · same UMAP')
        for ax in axes:ax.set(xlabel='UMAP 1',ylabel='UMAP 2',xticks=[],yticks=[])
        fig.suptitle('Endothelial nuclei · five vascular identity categories\nExisting all-stage UMAP; no age/batch integration; candidate status retained separately',fontsize=15)
        finish(fig,'all_stage_five_type_UMAP',pdf)
        for stage in STAGES:
            a=rna[rna.age.eq(stage)&rna.primary_EC_eligible];xy=a[['EC_UMAP1','EC_UMAP2']].to_numpy()
            fig,axes=plt.subplots(1,2,figsize=(16,7),layout='constrained')
            scatter(axes[0],xy,a.EC_type,size=8);legend(axes[0],a.EC_type);axes[0].set_title('Five EC types')
            for lib,c in zip(sorted(a['sample'].unique()),['#0072b2','#d55e00']):
                m=a['sample'].eq(lib).to_numpy();axes[1].scatter(xy[m,0],xy[m,1],s=7,color=c,label=lib,rasterized=True)
            axes[1].legend(fontsize=9,frameon=False);axes[1].set_title('RNA libraries · same coordinates')
            for ax in axes:ax.set(xlabel='UMAP 1',ylabel='UMAP 2',xticks=[],yticks=[])
            fig.suptitle(f'{stage} · {len(a):,} primary endothelial nuclei\nWithin-stage library adjustment; original UMAP preserved',fontsize=15)
            finish(fig,stage+'_five_type_RNA_UMAP',pdf)
        for stage in STAGES:
            a=x[x.stage.eq(stage)]
            atlas=pd.read_csv(ROOT/'results/whole_heart_spatial_atlas'/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0)
            fig,axes=plt.subplots(1,2,figsize=(17,8),layout='constrained')
            for ax,col,title in [(axes[0],'EC_type','Strict cell+nuclear agreement'),(axes[1],'nuclear_supported_EC_type','Nuclear-supported identities · separate evidence')]:
                ax.scatter(atlas.x_centroid,atlas.y_centroid,s=.15,color='#e6e6e6',rasterized=True)
                scatter(ax,a[['x_centroid','y_centroid']].to_numpy(),a[col],size=2.5)
                legend(ax,a[col]);ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='X (µm)',ylabel='Y (µm)');ax.set_title(title,fontsize=12)
            fig.suptitle(f'{stage} Xenium · five EC types\nGray: subtype unresolved under the displayed evidence criteria; source coordinates preserved',fontsize=15)
            finish(fig,stage+'_five_type_Xenium_map',pdf)
        fig,ax=plt.subplots(figsize=(12,6));ax.axis('off')
        table=ax.table(cellText=rna_counts.astype(int).to_numpy(),rowLabels=rna_counts.index,colLabels=STAGES,cellLoc='center',loc='center')
        table.auto_set_font_size(False);table.set_fontsize(13);table.scale(1,2)
        ax.set_title('Primary snRNA-seq counts · five EC types',fontsize=18,pad=20)
        fig.text(.08,.08,'Unresolved profiles remain in the dataset. Capture/QC selection and limited replication prevent biological abundance claims.',fontsize=10)
        pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    md=['| EC type | P0 | P7 | P14 | P21 |','|---|---:|---:|---:|---:|']
    for t,d in rna_counts.iterrows():md.append('| '+t+' | '+' | '.join(f'{int(v):,}' for v in d)+' |')
    (DEST/'README.md').write_text('''# Five endothelial types: primary presentation

The primary EC labels are **Capillary-like, Endocardial-like, Arterial-like, Venous-like and Lymphatic-like**. Finer gene-enrichment groups remain available only in the preceding audit. Cycling is a separate RNA annotation. Venous-like/lymphatic-like candidate status and all spatial confidence tiers remain explicit metadata. Unresolved profiles retain that label and are not forced into a type.

[Simplified report, 11 pages](five_endothelial_types_report.pdf) · [Marker key](five_type_marker_key.tsv) · [Stage counts](five_type_counts_by_stage.tsv) · [RNA annotations including QC-review nuclei](snRNA_five_type_annotations.tsv.gz) · [Xenium annotations](Xenium_five_type_annotations.tsv.gz)

## Primary snRNA classification

'''+ '\n'.join(md)+'''

Counts are QC-selected nuclei, not unbiased tissue abundances. Zero calls do not establish absence. In particular, arterial identity is unresolved at P7 under the current conservative model; no type was added to fill a stage quota. Small lymphatic-like groups and the P14 venous-like group remain provisional.

![Same RNA UMAP colored by five types and stages](all_stage_five_type_UMAP.png)

## Spatial evidence

`EC_type` is the main five-type classification under strict cell+nuclear agreement. `EC_type_including_cell_supported` adds the preceding moderate boundary-stable calls. `nuclear_supported_EC_type` and `additional_nuclear_only_candidate_type` expose nuclear evidence separately; mixed cell-body profiles are not silently promoted to clean whole-cell calls. Confidence, merge review, original atlas identity and source coordinates remain in the table.

The two spatial panels display strict cell+nuclear labels and nuclear-supported labels separately. Original gray review profiles remain in the background. These evidence levels do not measure independently validated accuracy.

## Scope and reproducibility

This is a presentation-only consolidation of the [completed endothelial QC and subclustering](../endothelial_QC_subtypes/README.md), source commit `d9f02610105a90ee821c2d93e05d9f77643e0e9f`. The existing 4,597 primary RNA nuclei, 976 RNA QC-review nuclei and 95,316 Xenium candidates (including all 69,617 original EC profiles) are accounted for. No model is retrained, RNA molecule removed, physical coordinate moved or UMAP recomputed. Source marker criteria and confidence remain unchanged. `Unresolved` is an assignment status, not a sixth biological EC type.

Full expression matrices and the detailed QC/segmentation evidence remain in the preceding results directory. Candidate labels are shortened for the main display, with their provisional status preserved in metadata and the marker key. The 24 finer RNA group IDs are retained only as audit references, not used as the main EC taxonomy.

Run `python analysis/present_five_endothelial_types.py` in the existing analysis environment. The script checks source identity correspondence, coordinates, coverage, category counts and preserved evidence flags; `validation.json` and `file_manifest.tsv` record the checks and SHA256 hashes.
''')
    pdfinfo=subprocess.check_output(['pdfinfo',str(DEST/'five_endothelial_types_report.pdf')],text=True)
    pages=int(re.search(r'^Pages:\s+(\d+)',pdfinfo,re.MULTILINE).group(1));assert pages==11
    write_json(DEST/'validation.json',{'status':'PASS','presentation_only':True,'primary_categories':TYPES,
        'RNA_total':len(rna),'RNA_primary':int(rna.primary_EC_eligible.sum()),'Xenium_candidates':len(x),'Xenium_previously_labeled_EC':int(x.native_EC_broad.sum()),
        'source_classifications_and_confidence_preserved':True,'source_coordinates_preserved':True,'RNA_counts_not_modified':True,
        'unresolved_preserved':True,'cycling_separate':True,'venous_lymphatic_candidate_status_preserved':True,'PDF_pages':pages})
    files=[p for p in DEST.iterdir() if p.is_file() and p.name!='file_manifest.tsv']+[Path(__file__)]
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]),DEST/'file_manifest.tsv')
    print(rna_counts.to_string(),flush=True)
    print('FIVE_EC_TYPES_COMPLETE',DEST,flush=True)

if __name__=='__main__':main()
