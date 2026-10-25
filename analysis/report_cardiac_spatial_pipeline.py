"""Reviewable atlas, spatial-distribution discoveries and association report."""
from cardiac_spatial_pipeline import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from scipy import ndimage
import textwrap

COLORS=dict(zip(TYPES,['#ce6a79','#377eb8','#37a36c','#9467bd','#d98b27','#c73932','#e5b828','#999999','#d0d0d0','#ededed']))
FBPAL=dict(zip(FBF,['#e69f00','#0072b2','#009e73','#cc79a7','#8b5c38']))
ECPAL=dict(zip(ECF,['#0072b2','#e69f00','#cc3311','#009e73','#cc79a7']))

def legend(ax,palette,fs=8):
    ax.legend(handles=[Line2D([],[],ls='',marker='o',color=c,markersize=4,label=k) for k,c in palette.items()],
              loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=fs)
def points(ax,f,c,s=.7,**kw):return ax.scatter(f.x_centroid,f.y_centroid,c=c,s=s,linewidths=0,rasterized=True,**kw)
def bg(ax,f):points(ax,f,'#ededed',.15)
def finish(ax):ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='X (µm)',ylabel='Y (µm)')
def writefig(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def cover(pdf):
    fig=plt.figure(figsize=(12,9));fig.text(.06,.94,'Cardiac spatial pipeline · integrated atlas and maturation',fontsize=19,weight='bold')
    lines=[
      'Scope:408,054 original spatial profiles, seven resolved cardiac lineages, four ages; unresolved/mixed profiles retained.',
      'Geometry: antibody-guided segmentation contacts plus30/50/100µm centroid graphs screened against the tissue-territory mask. Contact is an operational shared boundary, not independently verified membrane contact.',
      'QC: inherited broad identities, strict FB/EC subtype evidence, all-lineage nuclear/core marker audit, and a separately named conservative-mask sensitivity. Missing nuclear area alone never removes a CM.',
      'Niches: stage-balanced graph-neighbor composition, no focal RNA or maturation inputs, review burden excluded from clustering features, seed stability and assignment uncertainty exposed.',
      'Maturation: measured21-gene core/native indices plus four panel categories, equal-category relative score andMyl2-withheld score. No physiological maturity calibration is asserted.',
      'Non-CMs: measured FB matrix, EC cycling/guidance, immune, mural, pericyte and epicardial proxy programs; inferred sprouting is excluded from primary mechanism models.',
      'Inference: joint composition/program models, one common CM cohort per geometry, technical → anatomy → CM-state adjustment, paired contact/proximity comparisons,512 primary block-bootstrap draws.',
      'Prediction: spatial block holdout with150µm buffer for the primary50µm neighborhoods; compare baseline, composition, programs and joint predictors.',
      'Method assessment: CellCharter is optional niche discovery; Wasserstein Wormhole is optional distribution embedding. This release uses neither package. It exports compatible inputs and runs a small exact/sliced OT distance pilot.',
      'Biological scope: spatial association and candidate signaling hypotheses. One section per age cannot establish causal influence, longitudinal transitions or independently replicated age effects.'
    ]
    y=.85
    for line in lines:
        w=textwrap.fill(line,122);fig.text(.06,y,w,fontsize=10,va='top');y-=.024*(w.count('\n')+1)+.012
    pdf.savefig(fig);plt.close(fig)

def stage_maps(stage,pdf):
    f=frame(stage);n=read(f"{CFG['output']}/{stage}/niche_annotations.tsv.gz")
    q=read(f"{CFG['output']}/{stage}/all_lineage_count_view_QC.tsv.gz")
    fig,axes=plt.subplots(2,2,figsize=(19,15),layout='constrained')
    ax=axes[0,0]
    for t in STATUSES+LINEAGES:points(ax,f.loc[f.atlas_identity.eq(t)],COLORS[t],.7)
    legend(ax,COLORS);ax.set_title('All seven lineages · inherited broad labels; review profiles retained')
    ax=axes[0,1];bg(ax,f);cm=f.loc[f.ventricular_CM_strict]
    im=points(ax,cm,cm.CM_maturation_core,1,cmap='viridis',vmin=1.5,vmax=3.7)
    fig.colorbar(im,ax=ax,shrink=.55,label='Measured21-gene core mean log1p(panel CP10k)')
    ax.set_title('Ventricular CM maturation · boundary-core RNA')
    ax=axes[1,0];bg(ax,f)
    aliases=json.loads((OUT/'niche_methods.json').read_text())['aliases'];palette={name:matplotlib.colors.to_hex(plt.get_cmap('tab20')(i)) for i,name in enumerate(aliases.values())}
    for name,c in palette.items():points(ax,f.loc[n.niche.eq(name)],c,.8)
    legend(ax,palette,7);ax.set_title('Outcome-blinded computational contexts\nLow-information/uncertain assignments remain gray')
    ax=axes[1,1];points(ax,f,'#c8d5df',.25);points(ax,f.loc[q.conservative_QC_review],'#c73932',.65)
    legend(ax,{'No additional review flag':'#c8d5df','Count-view lineage review':'#c73932'})
    ax.set_title('All-lineage QC · conservative review flags\nFlag does not prove doublet or contamination')
    for ax in axes.ravel():finish(ax)
    fig.suptitle(stage+' · integrated spatial atlas, maturation, contexts and confidence',fontsize=18)
    writefig(fig,stage+'_integrated_atlas',pdf)
    fig,axes=plt.subplots(1,3,figsize=(23,8),layout='constrained')
    for ax in axes:bg(ax,f)
    evidence=pd.read_csv(ROOT/'results/ventricular_spatial_resolution/fine_group_evidence_summary.tsv',sep='\t').set_index('group')
    names=sorted(set(f.loc[f.atlas_identity.eq('Cardiomyocyte'),'fine_identity'])-{'Unresolved'})
    palette={}
    for i,name in enumerate(names):
        c=matplotlib.colors.to_hex(plt.get_cmap('tab10')(i));z=f.loc[f.atlas_identity.eq('Cardiomyocyte')&f.fine_identity.eq(name)]
        points(axes[0],z,c,1.7);palette[name+': '+evidence.loc[name,'mapping_family']]=c
    legend(axes[0],palette,7);axes[0].set_title('Prioritized CM candidates · finer labels not forced')
    for name,c in FBPAL.items():points(axes[1],f.loc[f.atlas_identity.eq('Fibroblast')&f.fine_identity.eq(name)],c,2)
    legend(axes[1],FBPAL);axes[1].set_title('Strict FB cell+nuclear family assignments')
    for name,c in ECPAL.items():points(axes[2],f.loc[f.atlas_identity.eq('Endothelial')&f.fine_identity.eq(name)],c,2)
    legend(axes[2],ECPAL);axes[2].set_title('Strict EC cell+nuclear subtype assignments')
    for ax in axes:finish(ax)
    fig.suptitle(stage+' · supported finer identities on original coordinates\nGray includes unresolved fine identity and other cardiac profiles; sparse categories are not absent-cell claims',fontsize=16)
    writefig(fig,stage+'_supported_subtypes',pdf)

def contact_examples(pdf):
    fig,axes=plt.subplots(2,2,figsize=(15,14),layout='constrained')
    for stage,ax in zip(STAGES,axes.ravel()):
        f=frame(stage);xy=f[['x_centroid','y_centroid']].to_numpy();tile=np.floor(xy/200).astype(int)
        codes=pd.Categorical(f.atlas_identity,categories=TYPES).codes
        _,inv=np.unique(tile,axis=0,return_inverse=True)
        count=np.zeros((inv.max()+1,len(TYPES)))
        np.add.at(count,(inv,codes),1)
        prob=count[:,:len(LINEAGES)]/np.maximum(count[:,:len(LINEAGES)].sum(1),1)[:,None]
        entropy=-np.sum(prob*np.log(np.maximum(prob,1e-8)),axis=1)
        entropy[count[:,:len(LINEAGES)].sum(1)<40]=-1
        chosen=int(np.argmax(entropy));center=xy[inv==chosen].mean(0);width=240
        with np.load(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/whole_section_segmentation.npz') as b:
            labels=b['cell_labels'];res=float(b['grid_um'])
        x0=max(0,int((center[0]-width/2)/res));x1=min(labels.shape[1],int((center[0]+width/2)/res))
        y0=max(0,int((center[1]-width/2)/res));y1=min(labels.shape[0],int((center[1]+width/2)/res))
        crop=labels[y0:y1,x0:x1]
        rgba=np.vstack([np.array([1,1,1,1]),np.array([matplotlib.colors.to_rgba(COLORS[t]) for t in f.atlas_identity])])
        im=rgba[crop];border=np.zeros(crop.shape,bool);border[:,1:]|=crop[:,1:]!=crop[:,:-1];border[1:]|=crop[1:]!=crop[:-1]
        im[border]=[.2,.2,.2,1]
        ax.imshow(im,extent=[x0*res,x1*res,y1*res,y0*res],interpolation='nearest');ax.set_title(stage+' · diversity-selected ROI');ax.set(xlabel='X (µm)',ylabel='Y (µm)')
        save(f.loc[(f.x_centroid.between(x0*res,x1*res))&(f.y_centroid.between(y0*res,y1*res)),['stage_cell_id','atlas_identity','x_centroid','y_centroid']],f'{stage}/contact_ROI_profiles.tsv.gz',index=True)
        dump(f'{stage}/contact_ROI_selection.json',{'x0_um':x0*res,'x1_um':x1*res,'y0_um':y0*res,'y1_um':y1*res,
            'selection':'Maximum resolved-lineage entropy in200um tiles with>=40 resolved profiles; no maturation outcome used'})
        del labels;gc.collect()
    fig.suptitle('Segmentation-aware spatial context · boundaries and labels\nAntibody-guided refinement retains deposited IDs; neither these boundaries nor labels are independent anatomical truth',fontsize=15)
    writefig(fig,'segmentation_contact_examples',pdf)

def enrichments(pdf):
    a=pd.read_csv(OUT/'contact_enrichment_local_null.tsv',sep='\t')
    fig,axes=plt.subplots(2,2,figsize=(15,12),layout='constrained')
    for stage,ax in zip(STAGES,axes.ravel()):
        mat=np.zeros((len(LINEAGES),len(LINEAGES)))
        for _,r in a.loc[a.stage.eq(stage)].iterrows():
            i=LINEAGES.index(r.cell_type_1);j=LINEAGES.index(r.cell_type_2);mat[i,j]=mat[j,i]=r.log2_observed_over_local_null
        im=ax.imshow(mat,cmap='RdBu_r',vmin=-1.5,vmax=1.5)
        ax.set_xticks(range(len(LINEAGES)),LINEAGES,rotation=40,ha='right',fontsize=8);ax.set_yticks(range(len(LINEAGES)),LINEAGES,fontsize=8);ax.set_title(stage)
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.55,label='log2 observed contacts / within200µm shuffled expectation')
    fig.suptitle('Cell-distribution hypotheses after preserving local composition\n199 within-section label permutations; conditional co-occurrence, not communication or animal replication',fontsize=15)
    writefig(fig,'contact_enrichment_local_context',pdf)

def associations(pdf):
    a=pd.read_csv(OUT/'joint_association_coefficients.tsv.gz',sep='\t')
    keys=['fraction__Fibroblast','fraction__Endothelial','program__FB_fibrillar_matrix','program__FB_basement_matrix','program__FB_elastic_matrix','program__FB_remodeling','program__EC_cycling','program__EC_guidance_proxy']
    fig,axes=plt.subplots(1,3,figsize=(20,9),layout='constrained',sharey=True)
    colors=['#0072b2','#e69f00','#009e73','#cc3311']
    for level,ax in zip(['technical','anatomy','CM_state'],axes):
        for j,stage in enumerate(STAGES):
            z=a.loc[a.stage.eq(stage)&a.geometry.eq('proximity50')&a.adjustment.eq(level)&a.outcome.eq('CM_maturation_core')].set_index('exposure').reindex(keys)
            y=np.arange(len(keys))+(j-1.5)*.16
            # Bootstrap intervals need not contain the point estimate; draw segments directly.
            ax.hlines(y,z.block_CI_low,z.block_CI_high,color=colors[j],lw=1.5)
            ax.plot(z.standardized_coefficient,y,'o',ms=4,color=colors[j],label=stage)
        ax.axvline(0,color='#999999',lw=1);ax.grid(axis='x',alpha=.2);ax.set_title(level+' adjustment');ax.set_xlabel('Maturation SD per exposure SD');ax.legend(frameon=False,fontsize=8)
    axes[0].set_yticks(range(len(keys)),[k.replace('fraction__','Neighborhood ').replace('program__','').replace('_',' ') for k in keys],fontsize=9);axes[0].invert_yaxis()
    fig.suptitle('Joint FB/EC and other non-CM models · fixed cohort across adjustment levels\nPrimary50µm tissue-screened neighborhoods;512 spatial-block draws; conditioning target changes across columns',fontsize=15)
    writefig(fig,'staged_nonCM_maturation_associations',pdf)
    pred=pd.read_csv(OUT/'spatial_prediction.tsv',sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for level,c in zip(['technical','anatomy','CM_state'],['#0072b2','#d55e00','#009e73']):
        z=pred.loc[pred.adjustment.eq(level)].pivot(index='stage',columns='predictors',values='OOF_R2').loc[STAGES]
        axes[0].plot(STAGES,z.baseline,'o--',c=c,label=level+' baseline');axes[0].plot(STAGES,z.joint,'o-',c=c,label=level+' + non-CMs')
        axes[1].plot(STAGES,z.joint-z.baseline,'o-',c=c,label=level)
    axes[0].set(ylabel='Held-out R²',title='Spatially buffered maturation prediction');axes[1].set(ylabel='Added held-out R²',title='Increment beyond each baseline');axes[1].axhline(0,c='#999999',lw=1)
    if (OUT/'nonlinear_prediction.tsv').exists():
        z=pd.read_csv(OUT/'nonlinear_prediction.tsv',sep='\t').pivot(index='stage',columns='predictors',values='OOF_R2').loc[STAGES]
        axes[1].plot(STAGES,z.joint-z.baseline,'s--',c='black',label='Nonlinear anatomy-adjusted')
    axes[0].legend(fontsize=8,frameon=False);axes[1].legend(frameon=False,fontsize=9)
    fig.suptitle('Prediction is distinct from association and causation\nFive400µm-block folds;150µm exclusion buffer; core21-gene outcome',fontsize=14)
    writefig(fig,'nonCM_spatial_prediction',pdf)

def methods_comparison(pdf):
    sel=pd.read_csv(OUT/'niche_model_selection.tsv',sep='\t')
    transport=pd.read_csv(OUT/'transport_pilot_summary.tsv',sep='\t').set_index('stage').loc[STAGES]
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    axes[0].plot(sel.k,sel.median_seed_ARI,'o-',label='Median seed ARI');axes[0].plot(sel.k,sel.minimum_seed_ARI,'s--',label='Minimum seed ARI');axes[0].axhline(.8,c='#999999',ls=':');axes[0].legend(frameon=False)
    axes[0].set(xlabel='Number of computational contexts',ylabel='Seed agreement (ARI)',title='Fine contexts require demonstrated stability')
    axes[1].plot(STAGES,transport.SW_vs_exact_Spearman,'o-',label='Sliced OT versus exact OT');axes[1].plot(STAGES,transport.mean_vs_exact_Spearman,'s-',label='Mean features versus exact OT')
    axes[1].set(ylabel='Spearman correlation,64 pairs/age',title='Pilot before training a transformer');axes[1].legend(fontsize=9,frameon=False)
    fig.suptitle('Method choice is evidence-based\nCellCharter/Wormhole optional; pilot uses observed non-CM type/program distributions and tests distances only',fontsize=14)
    writefig(fig,'niche_and_transport_method_assessment',pdf)

def findings(pdf):
    a=pd.read_csv(OUT/'candidate_relationships.tsv',sep='\t')
    selected=a.loc[a.priority_within_section_candidate]
    cov=pd.read_csv(OUT/'model_cohort_coverage.tsv',sep='\t')
    n=json.loads((OUT/'niche_methods.json').read_text())
    fig=plt.figure(figsize=(12,9));fig.text(.06,.94,'Findings and interpretation gates',fontsize=21,weight='bold')
    lines=[f'All408,054 original profiles remain accessible. Seven broad resolved lineages and conservative finer CM/FB/EC candidates are carried forward with explicit uncertainty.',
        f'Niche selection:K={n["chosen_k"]}. {n["status"]}. These are graph-context classes, not new cell identities or manually validated cardiac compartments.',
        'Contact-enrichment patterns include concentrated smooth-muscle, epicardial and immune neighborhoods. Local-shuffle comparisons distinguish these from a simple global abundance contrast.',
        'Joint primary models use measured non-CM RNA on one CM cohort at each age. EC cycling/guidance do not show a consistent negative relationship with the overall maturation panel acrossP0/P7.',
        'Higher local fibroblast fraction has a small negative maturation association at later stages. It is a relative-composition contrast, not proof that fibroblasts suppress maturation; individual FB matrix programs add weaker evidence.',
        f'{len(selected)} stage/exposure relationships meet the predeclared operational robustness screen (effect≥0.05SD, interval direction, native/balanced/Myl2-withheld score agreement, paired geometry agreement, technical-model agreement, low collinearity and conservative-QC support).',
        'The screen prioritizes follow-up within these sections. It is not a multiplicity-corrected discovery procedure or independently replicated validation. All tested effects and failed gates are exported.',
        'Candidate signaling tables use a small curated ligand/receptor-component list and opposite-library RNA expression support. Receptor complexes, protein abundance, activation and causal effects are not established.',
        'CellCharter could provide a niche-backend sensitivity. Wormhole could scale comparisons of richer neighborhood distributions. Neither fixes missing genes, uncertain segmentation, age/batch confounding or lack of biological replication.',
        'Next experimental priorities: independent age-matched sections/animals, histological chamber/compact-trabecular annotations, selected boundary review, and targeted perturbation of a nominated non-CM pathway with CM maturation readouts.'
    ]
    y=.85
    for line in lines:
        w=textwrap.fill(line,120);fig.text(.06,y,w,fontsize=10,va='top');y-=.024*(w.count('\n')+1)+.012
    pdf.savefig(fig);plt.close(fig)
    return selected

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    with PdfPages(OUT/'cardiac_spatial_pipeline_report.pdf') as pdf:
        cover(pdf)
        for stage in STAGES:stage_maps(stage,pdf)
        contact_examples(pdf);enrichments(pdf);associations(pdf);methods_comparison(pdf);selected=findings(pdf)
    print('REPORT_COMPLETE',len(selected),'priority candidates',flush=True)

if __name__=='__main__':main()
