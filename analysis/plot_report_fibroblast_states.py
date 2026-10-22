"""Create readable RNA/state, matrix, spatial, prediction and boundary figures."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibroblast-state-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib import patheffects
from fibroblast_state_common import *

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
STAGECOLORS=['#0072b2','#e69f00','#009e73','#cc6677']


def save(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    pdf.savefig(fig,bbox_inches='tight');plt.close(fig)


def clean_spatial(ax):
    ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')


def main():
    ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
    clusters=pd.read_csv(OUT/'RNA_cluster_summary.tsv',sep='\t')
    models=pd.read_csv(OUT/'spatial_blocked_prediction.tsv',sep='\t')
    coefs=pd.read_csv(OUT/'spatial_adjusted_coefficients.tsv',sep='\t')
    means=pd.read_csv(OUT/'RNA_corrected_matrix_pseudobulk_scores.tsv',sep='\t')
    with PdfPages(OUT/'fibroblast_states_matrix_spatial_report.pdf') as pdf:
        fig=plt.figure(figsize=(11.7,8.3));fig.patch.set_facecolor('white')
        fig.text(.07,.89,'Fibroblast states and the postnatal CM environment',fontsize=22,weight='bold')
        text=(
            'Completed reanalysis: stage-specific RNA discovery, matrix comparisons, and Xenium association tests.\n\n'
            f'{len(ann):,} original QC-retained fibroblast nuclei · eight libraries · four ages\n'
            f'{len(clusters)} computational groups; {int(ann.lineage_mixing_review_flag.sum()):,} nuclei in residual lineage-review groups\n'
            'No nuclei are deleted from the original data. Lineage-review groups are withheld from the fine reference.\n\n'
            'CM maturation uses the exact 21 genes displayed in manuscript Figure 4b.\n'
            'Primary score: mean log1p(CP10k); fixed matched-control and denominator sensitivities.\n'
            'This is an independent score using the intended panel, not an exact reproduction of the author index.\n\n'
            'Main findings\n'
            '• Early Postn-associated groups contain several distinct transcriptional programs.\n'
            '• A Gfpt2-associated group is already present at P7 and resembles the later group.\n'
            '• Basement-membrane and other matrix programs differ by state as well as stage.\n'
            '• Fine individual Xenium labels are much less well supported at P14/P21.\n'
            '• Fibroblast features add little to spatially held-out CM maturation prediction.\n'
            '• Modest Ankrd1 associations depend on the exposure/boundary definition.\n\n'
            'One spatial section per age. No lineage ancestry, matrix mechanics, or causal signaling is established.\n'
            'Spatial-block intervals describe variation within the observed section, not independent animals.'
        )
        fig.text(.07,.79,text,fontsize=12,va='top',linespacing=1.55)
        pdf.savefig(fig);plt.close(fig)
        # RNA family-colored UMAPs; numerical groups label the exact computed partitions.
        fig,axes=plt.subplots(1,4,figsize=(18,5),constrained_layout=True)
        for ax,stage in zip(axes,STAGES):
            a=ann[ann.age.eq(stage)]
            for fam,color in COLORS.items():
                m=a.fibroblast_family.eq(fam)
                ax.scatter(a.loc[m,'fibroblast_UMAP1'],a.loc[m,'fibroblast_UMAP2'],s=2,color=color,rasterized=True)
            for cid,group in a.groupby('fibroblast_cluster_id',observed=True):
                pos=group[['fibroblast_UMAP1','fibroblast_UMAP2']].median()
                label=cid.replace(stage+'_','')+('*' if group.lineage_mixing_review_flag.iloc[0] else '')
                ax.text(*pos,label,fontsize=11,weight='bold',ha='center',va='center',path_effects=[patheffects.withStroke(linewidth=3,foreground='white')])
            ax.set_title(stage+f' · n={len(a):,}',fontsize=14);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
            ax.set_xticks([]);ax.set_yticks([])
        handles=[Line2D([],[],marker='o',linestyle='',color=c,label=f,markersize=7) for f,c in COLORS.items()]
        fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
        fig.suptitle('Fibroblast RNA groups within each stage\n* residual lineage-mixing review; independent UMAP coordinates between stages',fontsize=16)
        save(fig,'RNA_fibroblast_states_UMAP',pdf)
        # Composition is explicitly a sampled-nucleus description.
        counts=pd.crosstab(ann.loc[ann.fine_reference_eligible,'age'],ann.loc[ann.fine_reference_eligible,'fibroblast_family']).reindex(STAGES).reindex(columns=FAMILIES+[UNRESOLVED],fill_value=0)
        save_tsv(counts.reset_index(),OUT/'RNA_family_counts_by_stage.tsv')
        fig,ax=plt.subplots(figsize=(11,5),constrained_layout=True)
        bottom=np.zeros(4)
        for fam in FAMILIES+[UNRESOLVED]:
            v=counts[fam].to_numpy()/counts.sum(1).to_numpy()
            ax.bar(STAGES,v,bottom=bottom,color=COLORS[fam],label=fam);bottom+=v
        ax.set_ylabel('Fraction of reviewed RNA reference nuclei');ax.set_ylim(0,1)
        ax.set_title('Stage-associated changes in the sampled fibroblast compartment\nCounts reflect nucleus recovery and operational labels; they do not establish population expansion or ancestry')
        ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',frameon=False)
        save(fig,'RNA_fibroblast_family_composition',pdf)
        # Matrix profiles, using only non-annotation genes.
        selected=means[means.version.eq('annotation_markers_withheld') & means.program.isin(MATRIX_NAMES)]
        table=selected.groupby(['age','family','program'],observed=True).score_mean_log1p_CPM.mean().unstack('program')
        keys=sorted(table.index,key=lambda p:(STAGES.index(p[0]),(FAMILIES+[UNRESOLVED]).index(p[1])))
        table=table.reindex(keys).reindex(columns=MATRIX_NAMES)
        z=(table-table.mean(0))/table.std(0)
        fig,ax=plt.subplots(figsize=(12,8),constrained_layout=True)
        im=ax.imshow(z.to_numpy(),cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
        ax.set_yticks(range(len(table)),[a+' · '+b for a,b in table.index],fontsize=10)
        ax.set_xticks(range(len(table.columns)),[s.replace('_','\n') for s in table.columns],fontsize=10)
        ax.set_title('Matrix RNA varies with fibroblast program and age\nAnnotation genes withheld; two library pseudobulks averaged; each program standardized separately')
        fig.colorbar(im,ax=ax,label='Relative program score (z)',shrink=.7)
        save(fig,'RNA_matrix_programs_by_state_stage',pdf)
        # Whole-transcriptome stage matches.
        sim=pd.read_csv(OUT/'RNA_cross_stage_signature_similarity.tsv',sep='\t',index_col=0)
        ids=sorted(sim.index,key=lambda k:(STAGES.index(k.split('_')[0]),int(k.split('F')[-1])))
        sim=sim.loc[ids,ids]
        fig,ax=plt.subplots(figsize=(11,9),constrained_layout=True)
        im=ax.imshow(sim,cmap='RdBu_r',vmin=-.5,vmax=1)
        ax.set_xticks(range(len(ids)),ids,rotation=90,fontsize=8);ax.set_yticks(range(len(ids)),ids,fontsize=8)
        ax.set_title('RNA resemblance with interpreted gene panels withheld\nSame-stage background removed; correlation does not demonstrate ancestry')
        fig.colorbar(im,ax=ax,label='Pearson correlation',shrink=.7)
        save(fig,'RNA_cross_stage_program_similarity',pdf)
        # Maturation distributions and spatial maps with native coordinates.
        maturation=[];strictcounts=[]
        for stage in STAGES:
            folder=OUT/stage
            fb=pd.read_csv(folder/'Xenium_fibroblast_states.tsv.gz',sep='\t',index_col=0)
            cm=pd.read_csv(folder/'CM_measured_outcomes.tsv.gz',sep='\t',index_col=0)
            for outcome in ['maturation_21gene_mean','maturation_matched_control']:
                maturation.append({'stage':stage,'outcome':outcome,'profiles':len(cm),'mean':cm[outcome].mean(),'median':cm[outcome].median(),
                    'p25':cm[outcome].quantile(.25),'p75':cm[outcome].quantile(.75)})
            strictcounts.append({'stage':stage,'broad_FBs':len(fb),'supported_fine_profiles':int(fb.supported_family_assignment.sum()),
                'fraction_supported':fb.supported_family_assignment.mean(),'raw_reference_agreement':fb.raw_reference_agreement.mean(),
                'both_library_agreement':fb.both_library_agreement.mean()})
            fig,axes=plt.subplots(1,4,figsize=(18,6),constrained_layout=True)
            for fam,color in COLORS.items():
                m=fb.mapped_family.eq(fam)
                axes[0].scatter(fb.loc[m,'x_centroid'],fb.loc[m,'y_centroid'],s=.8,c=color,rasterized=True)
            axes[0].set_title(f'Fibroblast labels meeting evidence rules\n{fb.supported_family_assignment.sum():,}/{len(fb):,}; gray = unresolved')
            for ax,frame,col,title,cmap,lim in [
                    (axes[1],cm,'maturation_21gene_mean','CM maturation: intended 21-gene panel','viridis',(1.5,4)),
                    (axes[2],cm,'Ankrd1_log1pCP10k','Measured CM Ankrd1','magma',(0,8)),
                    (axes[3],fb,'matrix__Basement_membrane','Measured FB basement-membrane RNA','viridis',(.7,4.2))]:
                im=ax.scatter(frame.x_centroid,frame.y_centroid,c=frame[col],s=.7,cmap=cmap,vmin=lim[0],vmax=lim[1],rasterized=True)
                ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.45,label='log-normalized RNA score')
            allxy=pd.concat([fb[['x_centroid','y_centroid']],cm[['x_centroid','y_centroid']]])
            for ax in axes:
                ax.set_xlim(allxy.x_centroid.min()-100,allxy.x_centroid.max()+100)
                ax.set_ylim(allxy.y_centroid.min()-100,allxy.y_centroid.max()+100)
                clean_spatial(ax)
            fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=9)
            fig.suptitle(stage+' · Original Xenium positions; no nucleus requirement or spatial label smoothing',fontsize=16)
            save(fig,stage+'_fibroblast_CM_spatial_maps',pdf)
            # Soft proposals are separate from supported labels, not silently upgraded.
            fig,ax=plt.subplots(figsize=(9,8),constrained_layout=True)
            for fam,color in COLORS.items():
                m=fb.family_proposal.eq(fam)
                ax.scatter(fb.loc[m,'x_centroid'],fb.loc[m,'y_centroid'],s=1,c=color,label=fam,alpha=.65,rasterized=True)
            clean_spatial(ax);ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False,markerscale=4)
            ax.set_title(stage+' · Fine-reference proposals for all broad fibroblasts\nIncludes uncertain profiles; colors are hypotheses, not verified fine identities')
            save(fig,stage+'_fibroblast_soft_proposals',pdf)
        save_tsv(pd.DataFrame(maturation),OUT/'CM_maturation_stage_summary.tsv')
        save_tsv(pd.DataFrame(strictcounts),OUT/'Xenium_fine_state_support_summary.tsv')
        fig,axes=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
        for ax,outcome,title in zip(axes,['maturation_21gene_mean','maturation_matched_control'],['21-gene mean log1p(CP10k)','Same panel minus fixed matched controls']):
            vals=[pd.read_csv(OUT/s/'CM_measured_outcomes.tsv.gz',sep='\t')[outcome].to_numpy() for s in STAGES]
            box=ax.boxplot(vals,tick_labels=STAGES,showfliers=False,patch_artist=True)
            for patch,color in zip(box['boxes'],STAGECOLORS):patch.set_facecolor(color);patch.set_alpha(.6)
            ax.set_title(title);ax.set_ylabel('Measured panel score')
        fig.suptitle('CM maturation distributions across the four observed sections\nOne section per age; no animal-level test or exact author-index reproduction',fontsize=14)
        save(fig,'CM_21gene_maturation_distributions',pdf)
        # Main spatial result is prediction on separated spatial blocks, not an in-sample p-value.
        fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
        for ax,outcome,title in zip(axes,['maturation_21gene_mean','Ankrd1_log1pCP10k'],['CM maturation','CM Ankrd1']):
            for j,model in enumerate(['families','matrix','both']):
                sub=models[models.outcome.eq(outcome)&models.model.eq(model)].set_index('stage').reindex(STAGES)
                ax.bar(np.arange(4)+(j-1)*.23,sub.incremental_OOF_R2,width=.21,label=model)
            ax.axhline(0,color='black',lw=.8);ax.set_xticks(range(4),STAGES)
            ax.set_title(title+('\nP0 complete model R² remains below zero' if outcome=='Ankrd1_log1pCP10k' else ''))
            ax.set_ylabel('Added spatially held-out R²');ax.legend(frameon=False)
        fig.suptitle('Additional prediction from fibroblast states / matrix programs\n50µm neighborhoods; 400µm blocks and 100µm training buffer; baseline includes vessels, coverage, composition and XY',fontsize=13)
        save(fig,'spatial_prediction_increment',pdf)
        fig,axes=plt.subplots(1,2,figsize=(14,7),constrained_layout=True)
        predictors=['matrix__'+p for p in MATRIX_NAMES[:5]]
        for ax,outcome,title in zip(axes,['maturation_21gene_mean','Ankrd1_log1pCP10k'],['CM maturation','CM Ankrd1']):
            for j,stage in enumerate(STAGES):
                sub=coefs[coefs.stage.eq(stage)&coefs.radius_um.eq(50)&coefs.model.eq('matrix')&coefs.outcome.eq(outcome)].set_index('predictor').reindex(predictors)
                y=np.arange(5)+(j-1.5)*.13
                # Intervals are bootstrap quantiles and need not be centered on the point estimate.
                ax.hlines(y,sub.spatial_block_bootstrap_low,sub.spatial_block_bootstrap_high,color=STAGECOLORS[j],lw=1.2)
                ax.scatter(sub.standardized_coefficient,y,c=STAGECOLORS[j],s=20,label=stage,zorder=3)
            ax.axvline(0,color='black',lw=.8);ax.set_yticks(range(5),[p.replace('_',' ') for p in MATRIX_NAMES[:5]])
            ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('Adjusted standardized coefficient');ax.legend(frameon=False)
        fig.suptitle('Conditional associations with measured neighboring fibroblast matrix RNA\n95% spatial-block bootstrap intervals describe the same section; no causal or animal-level interpretation',fontsize=14)
        save(fig,'spatial_matrix_association_coefficients',pdf)
        boundary=pd.read_csv(OUT/'boundary_spatial_association_sensitivity.tsv',sep='\t')
        fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
        for ax,outcome,title in zip(axes,['maturation_21gene_mean','Ankrd1_log1pCP10k'],['CM maturation','CM Ankrd1']):
            sub=boundary[boundary.outcome.eq(outcome)&boundary.predictor.eq('matrix__Basement_membrane')&boundary.scope.eq('image')]
            for variant,color in [('native_Q20','#0072b2'),('membrane_watershed','#d55e00')]:
                v=sub[sub.variant.eq(variant)].set_index('stage').reindex(STAGES)
                ax.plot(STAGES,v.partial_correlation,'o-',label=variant,color=color)
            ax.axhline(0,c='black',lw=.8);ax.set_title(title);ax.set_ylabel('Partial correlation with FB basement-membrane RNA');ax.legend(frameon=False)
        fig.suptitle('Boundary sensitivity on the identical image-covered CM subsets\n≥90% of 50µm FB neighbors audited; alternative segmentation is not anatomical truth',fontsize=14)
        save(fig,'boundary_matrix_association_sensitivity',pdf)
        pairs=pd.read_csv(OUT/'Ankrd1_vascular_pairs_after_FB_adjustment.tsv',sep='\t')
        fig,ax=plt.subplots(figsize=(11,6),constrained_layout=True)
        for j,(model,label,color) in enumerate([
                ('frozen_pairs_unadjusted','Same-pair native contrast','#333333'),
                ('matched_covariates_plus_CM_maturation','Covariates + CM maturation','#0072b2'),
                ('plus_FB_families_and_matrix','Also FB families + matrix','#d55e00')]):
            sub=pairs[pairs.model.eq(model)].set_index('stage').reindex(STAGES)
            x=np.arange(4)+(j-1)*.15
            ax.vlines(x,sub.segment_bootstrap_low,sub.segment_bootstrap_high,color=color,lw=1.5)
            ax.scatter(x,sub.conditional_near_minus_far_log_expression,c=color,s=30,label=label)
        ax.axhline(0,c='black',lw=.8);ax.set_xticks(range(4),[s+'\nn='+str(int(pairs[pairs.stage.eq(s)].same_matched_pairs.iloc[0])) for s in STAGES])
        ax.set_ylabel('Near minus far Ankrd1 log-expression contrast')
        ax.set_title('Frozen vascular pairs after conditioning on fibroblast environments\nSame neighborhood-supported subsets; segment-bootstrap intervals; P21 estimates unstable')
        ax.legend(frameon=False)
        save(fig,'Ankrd1_vascular_contrast_FB_adjustment',pdf)
    print('REPORT_WRITTEN',OUT/'fibroblast_states_matrix_spatial_report.pdf',flush=True)


if __name__=='__main__':main()
