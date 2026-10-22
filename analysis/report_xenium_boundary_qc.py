"""Summary figures/report for completed molecule and segmentation sensitivity."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from design_xenium_boundary_qc import OUT,STAGES

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
STAGE_COLORS=['#3579b5','#e68a26','#3b9a58','#bb4d79']
NAMES={'native_Q20':'Native Q20','core_1um':'Core >1 µm','core_2um':'Core >2 µm','dilation_2um':'Add unassigned ≤2 µm',
    'native_Q30':'Native Q30','thin_to_core1':'Random: core1 depth','thin_to_core2':'Random: core2 depth','thin_to_Q30':'Random: Q30 depth',
    'native_raster_Q20':'Native raster control','centroid_Voronoi':'Centroid Voronoi','membrane_watershed':'Membrane watershed'}


def table(name):return pd.read_csv(OUT/(name+'.tsv'),sep='\t')


def footer(fig,text):fig.text(.5,.018,text,ha='center',va='bottom',fontsize=9,linespacing=1.4)


def save(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=180,facecolor='white');fig.savefig(OUT/(name+'.pdf'),facecolor='white')
    pdf.savefig(fig,facecolor='white');plt.close(fig)


def overview(pdf):
    qc=table('whole_run_molecule_QC').set_index('stage').loc[STAGES]
    facing=table('CM_marker_boundary_localization')
    full=table('primary_Ankrd1_assignment_sensitivity');image=table('image_ROI_Ankrd1_assignment_sensitivity')
    fig,axes=plt.subplots(2,2,figsize=(14,11));x=np.arange(4)
    ax=axes[0,0]
    ax.bar(x-.18,qc.gene_Q20_fraction,.36,color='#527db4',label='QV≥20 / all decoded gene RNA')
    ax.bar(x+.18,qc.assigned_Q20_gene_fraction,.36,color='#2d9e88',label='Assigned / QV≥20 gene RNA')
    ax.set_xticks(x,STAGES);ax.set_ylim(0,1);ax.set_ylabel('Fraction');ax.legend(fontsize=8,frameon=False)
    ax.set_title('Whole-run transcript quality and assignment',weight='bold')
    ax=axes[0,1]
    for lineage,style,color in [('Endothelial','-','#008fa5'),('Fibroblast','-','#ad43b4')]:
        p=facing[facing.lineage.eq(lineage)].pivot(index='stage',columns='region',values='CM_marker_fraction').loc[STAGES]
        ratio=p.opposite_facing_shell/p.deep_core
        ax.plot(x,ratio,'o'+style,c=color,label=lineage+': facing / core')
        ax.plot(x,p.opposite_facing_shell/p.away_shell,'s:',c=color,label=lineage+': facing / other border',alpha=.8)
    ax.axhline(1,c='black',ls=':',lw=1);ax.set_xticks(x,STAGES);ax.set_ylim(bottom=0)
    ax.set_ylabel('CM-facing marker fraction / comparison fraction')
    ax.set_title('Boundary enrichment has imperfect direction specificity',weight='bold');ax.legend(frameon=False,fontsize=8,loc='lower left')
    ax=axes[1,0]
    for variant in ['native_Q20','core_1um','core_2um','dilation_2um','native_Q30']:
        p=full[full.variant.eq(variant)].set_index('stage').loc[STAGES]
        ax.plot(x,p.near_far_mean_Ankrd1_CP10k_ratio,'o-',label=NAMES[variant],lw=1.4)
    ax.axhline(1,c='black',ls=':',lw=1);ax.set_xticks(x,STAGES);ax.set_ylabel('Near/control mean measured Ankrd1 CP10k ratio')
    ax.set_title('All 3,104 original matched CM pairs',weight='bold');ax.legend(fontsize=8,frameon=False)
    ax=axes[1,1]
    for offset,variant,color,label in [(-.18,'native_Q20','#a4b9d6','Original cells'),(.18,'membrane_watershed','#527db4','Membrane alternative')]:
        p=image[image.variant.eq(variant)].set_index('stage').loc[STAGES]
        ax.bar(x+offset,p.near_far_mean_Ankrd1_CP10k_ratio,.36,color=color,label=label)
        if variant=='native_Q20':
            for i,s in enumerate(STAGES):ax.text(i,3.3,f'n={int(p.loc[s,"pairs"])}',ha='center',fontsize=9)
    ax.axhline(1,c='black',ls=':',lw=1);ax.set_xticks(x,STAGES);ax.set_ylim(0,3.6)
    ax.set_ylabel('Near/control Ankrd1 ratio in image-covered pairs')
    ax.set_title('Image comparisons use the same smaller subset',weight='bold');ax.legend(fontsize=8,frameon=False)
    fig.suptitle('Xenium transcript-boundary audit and segmentation sensitivity',fontsize=19,weight='bold')
    footer(fig,'Molecule reconstruction matches every deposited count matrix exactly at QV≥20. No nucleus is required for a CM profile.\nBoundary-localized marker RNA is compatible with mixing, not proof of contamination. Alternative segmentation is not validated ground truth.\nOne section per stage; ratios are within-section contrasts. Full-pair and image-region subsets must not be treated as the same cohort.')
    fig.subplots_adjust(left=.07,right=.97,top=.91,bottom=.13,hspace=.35,wspace=.28)
    save(fig,'boundary_QC_key_results',pdf)


def identities(pdf):
    data=pd.concat([table('identity_and_information_retention'),table('depth_control_identity_retention')],ignore_index=True)
    row_keys=[(s,c) for s in STAGES for c in ['Cardiomyocyte','Endothelial','Fibroblast','Smooth muscle']]
    geometry=['core_1um','thin_to_core1','core_2um','thin_to_core2','dilation_2um','native_Q30','thin_to_Q30']
    morphology=['native_raster_Q20','centroid_Voronoi','membrane_watershed']
    fig,axes=plt.subplots(1,2,figsize=(19,12),gridspec_kw={'width_ratios':[2.2,1]})
    for ax,columns,title in zip(axes,[geometry,morphology],['Boundary/QV changes and count-matched controls','Image-covered cells: boundary alternatives']):
        values=np.full((len(row_keys),len(columns)),np.nan);labels=[]
        for i,(s,lineage) in enumerate(row_keys):
            n=None
            for j,variant in enumerate(columns):
                part=data[data.stage.eq(s)&data.variant.eq(variant)&data.baseline_identity.eq(lineage)&data.nuclear_area.eq('all')]
                if len(part):
                    r=part.iloc[0];values[i,j]=r.identity_unchanged/r.profiles;n=int(r.profiles)
            short={'Cardiomyocyte':'CM','Endothelial':'EC','Fibroblast':'Fib','Smooth muscle':'SMC'}[lineage]
            labels.append(f'{s} · {short} (n={n})')
        img=ax.imshow(values,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
        ax.set_xticks(range(len(columns)),[NAMES[c] for c in columns],rotation=55,ha='right',fontsize=9)
        ax.set_yticks(range(len(row_keys)),labels,fontsize=9)
        for i in range(len(row_keys)):
            for j in range(len(columns)):
                if np.isfinite(values[i,j]):ax.text(j,i,f'{values[i,j]:.0%}',ha='center',va='center',fontsize=9,color='white' if values[i,j]>.70 else 'black')
        for cut in [3.5,7.5,11.5]:ax.axhline(cut,c='white',lw=2)
        ax.set_title(title,weight='bold',fontsize=12)
    fig.suptitle('Broad-identity stability depends on boundaries and measurement depth',fontsize=18,weight='bold')
    cax=fig.add_axes([.94,.32,.013,.45]);fig.colorbar(img,cax=cax,label='Fraction retaining the original broad label')
    footer(fig,'Label retention is not an accuracy estimate or segmentation-error rate. Lost labels include uncertainty and insufficient RNA, not only lineage switches.\nRandom controls sample original transcripts without replacement to exactly match each cell’s eroded/Q30 transcript total (one fixed technical draw).\nThe two panels use different cohorts; image comparisons require adequate crop coverage. All original seeds, including anucleate profiles, are retained.')
    fig.subplots_adjust(left=.12,right=.90,top=.91,bottom=.23,wspace=.44)
    save(fig,'boundary_identity_and_depth_controls',pdf)


def anucleate(pdf):
    data=table('identity_and_information_retention');variants=['native_Q20','core_1um','core_2um','dilation_2um','native_Q30']
    fig,axes=plt.subplots(2,2,figsize=(15,11))
    for ax,s in zip(axes.ravel(),STAGES):
        for status,color,marker in [('measured','#527db4','o'),('missing','#c64f84','s')]:
            p=data[data.stage.eq(s)&data.baseline_identity.eq('Cardiomyocyte')&data.nuclear_area.eq(status)].set_index('variant').loc[variants]
            n=int(p.profiles.iloc[0]);ax.plot(range(5),p.identity_unchanged/p.profiles,marker=marker,c=color,label=f'Nuclear area {status} (n={n:,})')
            ax.plot(range(5),p.information_QC_pass/p.profiles,c=color,ls=':',alpha=.7)
        ax.set_xticks(range(5),['Native','Core1','Core2','Dilate2','Q30']);ax.set_ylim(0,1.05);ax.set_title(s,weight='bold')
        ax.set_ylabel('Fraction of original CM profiles');ax.legend(fontsize=8,frameon=False,loc='lower left')
    fig.suptitle('CM profiles without measured nuclear area are retained',fontsize=19,weight='bold')
    footer(fig,'Solid: same CM identity; dotted: sufficient expression information. No nucleus criterion enters classification or alternative seed placement.\nLarge RNA losses under erosion are stress tests, not recommendations to discard legitimate anucleate CM profiles.\nNuclear-area groups differ biologically and technically; these curves do not estimate a causal effect of nuclear presence.')
    fig.subplots_adjust(left=.07,right=.97,top=.91,bottom=.13,hspace=.28,wspace=.23)
    save(fig,'anucleate_CM_retention',pdf)


def vessels(pdf):
    d=table('selected_vessel_envelope_sensitivity')
    fig,axes=plt.subplots(1,2,figsize=(15,8))
    original=d[d.variant.eq('native_ROI_Q20')&d.pair_subset.eq('same_original_pairs')].set_index('roi_id')
    new=d[d.variant.eq('membrane_watershed')&d.pair_subset.eq('same_original_pairs')].set_index('roi_id')
    endpoints=[]
    for i,(roi,row) in enumerate(new.iterrows()):
        color=STAGE_COLORS[STAGES.index(row.stage)]
        axes[0].scatter(row.envelope_IoU,row.fraction_original_distance_pair_categories_preserved,c=color,s=55)
        axes[0].annotate(str(i+1),(row.envelope_IoU,row.fraction_original_distance_pair_categories_preserved),xytext=(5,5),textcoords='offset points',fontsize=9)
        axes[1].plot([0,1],[original.loc[roi,'Ankrd1_near_far_ratio'],row.Ankrd1_near_far_ratio],'o-',c=color,alpha=.8)
        endpoints.append((row.Ankrd1_near_far_ratio,roi,color,i+1))
    previous=-np.inf
    for value,roi,color,index in sorted(endpoints):
        y=max(value,previous+.20);previous=y
        axes[1].annotate(f'{index}: {roi}',(1,value),xytext=(1.12,y),fontsize=8,color=color,
            arrowprops={'arrowstyle':'-','color':color,'lw':.6})
    axes[0].set_xlim(.55,.95);axes[0].set_ylim(.75,1.0);axes[0].set_xlabel('Native / alternative vascular-envelope IoU')
    axes[0].set_ylabel('Fraction of original pairs keeping near/far distance categories')
    axes[0].set_title('Eight selected arterial regions: geometry changes',weight='bold')
    axes[1].set_xticks([0,1],['Native RNA/boundaries','Membrane alternative']);axes[1].set_xlim(-.15,1.65)
    axes[1].axhline(1,c='black',ls=':',lw=1);axes[1].set_ylabel('Near/control Ankrd1 ratio in the same local pairs')
    axes[1].set_title('RNA reassignment and vessel-envelope sensitivity',weight='bold')
    footer(fig,'Same original arterial EC/SMC memberships; envelopes rebuilt from each alternative’s cell polygons, buffered3 µm.\nComparisons use crop-covered original pairs. Additional distance-preserved-pair results are in the table.\nThese are eight selected vessel regions, not eight independent animals or a complete reannotation of the vasculature.')
    fig.subplots_adjust(left=.07,right=.96,top=.90,bottom=.16,wspace=.31)
    save(fig,'selected_vessel_boundary_sensitivity',pdf)


def main():
    with PdfPages(OUT/'Xenium_transcript_boundary_QC_report.pdf') as pdf:
        overview(pdf);identities(pdf);anucleate(pdf);vessels(pdf)
        # Every ROI remains in the per-stage8-page galleries; this report selects4/stage.
        for stage in STAGES:
            for suffix in ['Endothelial_flagged','Fibroblast_flagged','anucleate_CM_1','arterial_1']:
                path=OUT/stage/(stage+'_'+suffix+'_transcript_boundary_inspection.png')
                image=plt.imread(path);fig=plt.figure(figsize=(18,6.5));ax=fig.add_axes([0,0,1,1]);ax.imshow(image);ax.axis('off')
                pdf.savefig(fig,facecolor='white');plt.close(fig)
    print('BOUNDARY_QC_REPORT_COMPLETE',flush=True)


if __name__=='__main__':main()
