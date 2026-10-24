"""Render whole-heart lineage reassessment and measured Ankrd1 spatial contrasts."""
import json
import importlib.metadata
import os
import sys
from datetime import datetime,timezone
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

from reassess_whole_heart_lineages import OUT,ROOT,STAGES,digest
from plot_xenium_ventricular import background,frame

COLORS=dict(zip(STAGES,['#3579b5','#e68a26','#3b9a58','#bb4d79']))
LINEAGES={'Cardiomyocyte':'#466bb0','Endothelial':'#29998e','Fibroblast':'#b45c99',
    'Pericyte':'#977449','Smooth muscle':'#d5a829','Immune':'#cf4c4c','Epicardial':'#8263b4',
    'Glial':'#7a9c43','Ambiguous/review':'#b8b8b8','Mixed CM/non-CM; review':'#777777',
    'QC/information flag':'#e0e0e0'}
DEFS=['independent_strict','independent_marker','published_comparator']
DLABELS=['Independent strict','Marker sensitivity','Published comparator']
COHORTS=['whole_reference_strict','whole_reference_inclusive','annotation_consensus',
         'published_ventricular','strict_Pearson_agreement','previous_CM_gate']
CLABELS=['Strict ventricular','Inclusive ventricular','Annotation consensus',
         'Published ventricular','Strict + Pearson agreement','Previous CM gate']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})


def table(name):
    return pd.read_csv(OUT/name,sep='\t',float_precision='round_trip')


def save(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=160,facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),facecolor='white')
    pdf.savefig(fig,facecolor='white');plt.close(fig)


def outline(ax,geom,color='#222222',lw=.8,ls='-',alpha=1):
    if geom.geom_type=='Polygon':
        for ring in [geom.exterior,*geom.interiors]:
            xy=np.asarray(ring.coords);ax.plot(xy[:,0],xy[:,1],color=color,lw=lw,ls=ls,alpha=alpha)
    elif hasattr(geom,'geoms'):
        for g in geom.geoms:outline(ax,g,color,lw,ls,alpha)


def footer(fig,text):
    fig.text(.5,.018,text,ha='center',va='bottom',fontsize=9,linespacing=1.4)


def summary(pdf):
    primary=table('primary_local_Ankrd1_results.tsv').set_index('stage').loc[STAGES]
    profiles=table('distance_expression_profiles.tsv')
    effects=table('segment_effects.tsv')
    fig,axes=plt.subplots(2,2,figsize=(13,10))
    ax=axes[0,0]
    values=primary.near_far_Ankrd1_CP10k_ratio.to_numpy()
    ax.bar(STAGES,values,color=[COLORS[s] for s in STAGES],width=.6)
    ax.axhline(1,c='black',lw=1,ls='--');ax.set_ylim(0,1.95)
    for i,(s,row) in enumerate(primary.iterrows()):
        ax.text(i,values[i]+.04,f'{values[i]:.2f}×\n{row.matched_pairs:,} pairs',ha='center',fontsize=10)
    ax.set_title('Matched Ankrd1 expression: near / local control',weight='bold')
    ax.set_ylabel('Ratio of mean CP10k expression')
    ax=axes[0,1]
    for i,s in enumerate(STAGES):
        selected=effects[effects.stage.eq(s)&effects.definition.eq('independent_strict')&effects.cohort.eq('whole_reference_strict')]
        y=selected.Ankrd1_difference.to_numpy();jitter=np.random.default_rng(17).uniform(-.16,.16,len(y))
        ax.scatter(i+jitter,y,color=COLORS[s],s=18,alpha=.7)
        ax.plot([i-.25,i+.25],[y.mean()]*2,c='black',lw=2)
        ax.text(i,1.65,f'{len(y)} segments',ha='center',fontsize=9)
    ax.axhline(0,c='black',lw=1,ls='--');ax.set_xticks(range(4),STAGES)
    ax.set_title('Each dot: one candidate vascular segment',weight='bold')
    ax.set_ylabel('Mean paired Ankrd1 log1p(CP10k) difference')
    ax=axes[1,0]
    for s in STAGES:
        p=profiles[profiles.stage.eq(s)&profiles.definition.eq('independent_strict')&profiles.cohort.eq('whole_reference_strict')].iloc[:6]
        ax.plot((p.distance_low_um+p.distance_high_um)/2,p.mean_Ankrd1_log1pCP10k,'o-',color=COLORS[s],label=s)
    ax.set_title('Unmatched descriptive distance profiles',weight='bold')
    ax.set_xlabel('Distance to vascular-cell envelope (µm)');ax.set_ylabel('Mean Ankrd1 log1p(CP10k)');ax.legend(frameon=False)
    ax=axes[1,1]
    x=np.arange(4);ax.bar(x-.17,primary.near_Ankrd1_high_fraction,width=.34,color=[COLORS[s] for s in STAGES],label='Near, ≤30 µm')
    ax.bar(x+.17,primary.far_Ankrd1_high_fraction,width=.34,color=[COLORS[s] for s in STAGES],alpha=.35,label='Control, 60–150 µm')
    ax.set_xticks(x,STAGES);ax.set_ylim(0,.62);ax.set_ylabel('Fraction above stage-specific upper quartile')
    ax.set_title('Measured Ankrd1-high fraction in matched pairs',weight='bold');ax.legend(frameon=False,fontsize=9)
    fig.suptitle('Ankrd1 and candidate arterial neighborhoods: P0–P21',fontsize=18,weight='bold')
    footer(fig,'Primary: strict whole-reference ventricular CMs; independently defined arterial EC + SMC segments.\nNear: 0–30 µm; local controls: 60–150 µm; 1:1 covariate matching within the same segment.\nOne section per age. Segment dots are exploratory within-section summaries, not independent animals.')
    fig.subplots_adjust(top=.90,bottom=.13,hspace=.40,wspace=.27)
    save(fig,'Ankrd1_arterial_neighborhood_summary',pdf)


def stage_maps(stage,pdf):
    folder=OUT/stage;obs=pd.read_csv(folder/'whole_heart_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    fs=json.loads((folder/'arterial_segments.geojson').read_text())['features']
    features=[f for f in fs if f['properties']['definition']=='independent_strict' and f['properties']['primary_small_segment']]
    bg,extent,bounds=background(stage,obs)
    fig,axes=plt.subplots(1,3,figsize=(19,7.5))
    frame(axes[0],bg,extent,bounds,'Competing whole-heart lineage assignments')
    handles=[]
    for label,group in obs.groupby('reassessed_lineage',sort=True):
        color=LINEAGES.get(label,'#dddddd')
        axes[0].scatter(group.x_centroid,group.y_centroid,s=.55,c=color,edgecolors='none',rasterized=True)
        handles.append(Line2D([],[],marker='o',linestyle='',color=color,markersize=4,label=f'{label}: {len(group):,}'))
    axes[0].legend(handles=handles,fontsize=7,loc='upper left',framealpha=.9)
    frame(axes[1],bg,extent,bounds,f'{len(features)} candidate small arterial segments')
    axes[1].scatter(obs.x_centroid,obs.y_centroid,s=.15,c='#cccccc',edgecolors='none',rasterized=True)
    for key,color,label in [('arterial_EC_ids','#1678bf','Arterial EC'),('SMC_ids','#d5423c','SMC')]:
        ids=[cell for f in features for cell in f['properties'][key]];o=obs.loc[ids]
        axes[1].scatter(o.x_centroid,o.y_centroid,s=3,c=color,edgecolors='none',label=label,rasterized=True)
    for f in features:outline(axes[1],shape(f['geometry']),lw=.55)
    axes[1].legend(loc='upper left',fontsize=8)
    frame(axes[2],bg,extent,bounds,'Measured Ankrd1 in strict ventricular CMs')
    cm=obs[obs.ventricular_CM_strict].sort_values('Ankrd1_log1pCP10k')
    pts=axes[2].scatter(cm.x_centroid,cm.y_centroid,c=cm.Ankrd1_log1pCP10k,s=.7,vmin=0,vmax=7,cmap='magma_r',edgecolors='none',rasterized=True)
    fig.colorbar(pts,ax=axes[2],shrink=.55,pad=.02,label='Ankrd1 log1p(CP10k)')
    fig.suptitle(f'{stage}: whole-heart reassessment and vascular context',fontsize=19,weight='bold')
    footer(fig,'Original cell centroids and DAPI coordinates. All original cells retained; uncertain identities remain flagged.\nArterial candidates require arterial EC + surrounding smooth muscle; outlines are vascular-cell envelopes, not traced vessel walls.\nAnkrd1, Creb5, Tlr4 and Atf3 were excluded from all new identity fits and fitting size-factor denominators.')
    fig.subplots_adjust(left=.04,right=.98,top=.90,bottom=.14,wspace=.24)
    save(fig,f'{stage}_whole_heart_arterial_maps',pdf)

    pairs=pd.read_csv(folder/'default_contrast_matched_pairs.tsv.gz',sep='\t')
    chosen=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')]
    top=chosen.groupby('segment_id').size().sort_values(ascending=False).head(4)
    byid={f['properties']['segment_id']:f for f in features}
    fig,axes=plt.subplots(2,2,figsize=(13,12))
    for ax,(sid,n) in zip(axes.ravel(),top.items()):
        f=byid[sid];g=shape(f['geometry']);p=f['properties']
        x0,y0,x1,y1=g.bounds;b=[x0-165,x1+165,y0-165,y1+165]
        frame(ax,bg,extent,b,f'{sid.split("_")[-1]}: {n} matched pairs; {p["arterial_ECs"]} EC + {p["SMC_cells"]} SMC')
        local=cm[cm.x_centroid.between(b[0],b[1])&cm.y_centroid.between(b[2],b[3])]
        pts=ax.scatter(local.x_centroid,local.y_centroid,c=local.Ankrd1_log1pCP10k,s=16,vmin=0,vmax=7,cmap='magma_r',edgecolors='none',rasterized=True)
        for k,c,m in [('arterial_EC_ids','#1678bf','s'),('SMC_ids','#20a69c','^')]:
            o=obs.loc[p[k]];ax.scatter(o.x_centroid,o.y_centroid,s=30,c=c,marker=m,edgecolors='white',lw=.3)
        outline(ax,g,lw=1.3);outline(ax,g.buffer(30),'#2266aa',1.1,'--');outline(ax,g.buffer(60),'#555555',.8,':');outline(ax,g.buffer(150),'#555555',.8,':')
        fig.colorbar(pts,ax=ax,shrink=.7,pad=.02,label='Ankrd1 log1p(CP10k)')
    for ax in axes.ravel()[len(top):]:ax.set_visible(False)
    fig.suptitle(f'{stage}: candidate vessel geometry and nearby ventricular CMs',fontsize=17,weight='bold')
    footer(fig,'Examples selected by most matched pairs, without selecting on Ankrd1 effects. Square: arterial EC; triangle: SMC.\nBlack: vascular-cell envelope; blue dashed: 30 µm; gray dotted: 60 and 150 µm.\nCM color shows directly measured expression, without smoothing or imputation. Geometry requires anatomical review.')
    fig.subplots_adjust(left=.07,right=.95,top=.91,bottom=.12,hspace=.29,wspace=.30)
    save(fig,f'{stage}_arterial_neighborhood_zooms',pdf)


def sensitivities(pdf):
    data=table('local_Ankrd1_contrasts.tsv')
    fig,axes=plt.subplots(1,3,figsize=(17,9))
    norm=TwoSlopeNorm(vmin=-.7,vcenter=0,vmax=.7)
    for ax,definition,label in zip(axes,DEFS,DLABELS):
        sel=data[data.definition.eq(definition)&data.distance_metric.eq('boundary')&data.near_limit_um.eq(30)]
        a=sel.pivot(index='cohort',columns='stage',values='equal_segment_mean_Ankrd1_log_difference').reindex(index=COHORTS,columns=STAGES)
        ratio=sel.pivot(index='cohort',columns='stage',values='near_far_Ankrd1_CP10k_ratio').reindex(index=COHORTS,columns=STAGES)
        n=sel.pivot(index='cohort',columns='stage',values='matched_pairs').reindex(index=COHORTS,columns=STAGES)
        img=ax.imshow(a,norm=norm,cmap='RdBu_r',aspect='auto')
        ax.set_xticks(range(4),STAGES);ax.set_yticks(range(6),CLABELS if ax is axes[0] else [])
        for i in range(6):
            for j in range(4):ax.text(j,i,f'{ratio.iloc[i,j]:.2f}×\nn={int(n.iloc[i,j]):,}',ha='center',va='center',fontsize=9,color='white' if abs(a.iloc[i,j])>.45 else 'black')
        ax.set_title(label,weight='bold')
    fig.colorbar(img,ax=axes.tolist(),shrink=.65,pad=.025,label='Equal-segment paired log-expression difference')
    fig.suptitle('All identity and vessel-definition sensitivities',fontsize=18,weight='bold')
    footer(fig,'Text: near/control Ankrd1 mean CP10k ratio and matched-pair count. Color: equal-segment mean paired log difference.\nDefault 30 vs 60–150 µm contrast. Different definitions select different cells and structures.\nThese are sensitivity analyses of the same four sections, not independent replications.')
    fig.subplots_adjust(left=.18,right=.86,top=.89,bottom=.14,wspace=.12)
    # Explicit colorbar position avoids overlap when the plotting area is adjusted.
    fig.axes[-1].set_position([.89,.25,.018,.52])
    save(fig,'identity_vessel_sensitivity',pdf)

    fig,axes=plt.subplots(1,2,figsize=(14,7))
    rows=[('boundary',20,'20 vs 40–100 µm'),('boundary',30,'30 vs 60–150 µm'),('boundary',50,'50 vs 100–200 µm'),
          ('unbuffered',30,'Unbuffered cells: 30 vs 60–150'),('hull',30,'Convex hull: 30 vs 60–150'),('EC_centroid',30,'EC centroid: 30 vs 60–150')]
    for j,(metric,near,label) in enumerate(rows):
        sel=data[data.definition.eq('independent_strict')&data.cohort.eq('whole_reference_strict')&data.distance_metric.eq(metric)&data.near_limit_um.eq(near)].set_index('stage').loc[STAGES]
        for ax,key in zip(axes,['near_far_Ankrd1_CP10k_ratio','equal_segment_mean_Ankrd1_log_difference']):
            ax.plot(range(4),sel[key],marker='o',label=label,alpha=.85)
    for ax in axes:ax.set_xticks(range(4),STAGES);ax.grid(axis='y',alpha=.2)
    axes[0].axhline(1,c='black',lw=.8,ls='--');axes[1].axhline(0,c='black',lw=.8,ls='--')
    axes[0].set_ylabel('Ratio of near/control mean Ankrd1 CP10k');axes[1].set_ylabel('Equal-segment mean paired log difference')
    axes[0].legend(fontsize=8,frameon=False,loc='upper left')
    fig.suptitle('Distance and geometry sensitivity: strict ventricular CMs',fontsize=17,weight='bold')
    footer(fig,'All alternatives shown; vascular membership and small-segment restrictions held fixed.\nDistances to envelopes, hulls and EC centroids are different biological approximations; matching cohorts can change.')
    fig.subplots_adjust(left=.07,right=.97,top=.89,bottom=.16,wspace=.25)
    save(fig,'distance_geometry_sensitivity',pdf)


def concordance(pdf):
    d=table('CM_lineage_concordance.tsv');fig,ax=plt.subplots(figsize=(11,7))
    for offset,cohort,color,label in [(-.25,'previous_CM_gate','#b8b8b8','Previous CM gate'),(0,'new_strict_CM','#436cad','New strict CM'),(.25,'new_strict_ventricular','#24958c','New strict ventricular')]:
        p=d[d.cohort.eq(cohort)].set_index('stage').loc[STAGES]
        x=np.arange(4)+offset;ax.bar(x,p.published_CM_fraction,width=.24,color=color,label=label)
        for pos,(_,r) in zip(x,p.iterrows()):ax.text(pos,r.published_CM_fraction+.015,f'{r.published_CM_fraction:.0%}',ha='center',fontsize=9)
    ax.set_xticks(range(4),STAGES);ax.set_ylim(0,1.13);ax.set_ylabel('Fraction with a published cardiomyocyte annotation')
    ax.legend(loc='lower left',frameon=False);ax.set_title('Whole-reference CM assignments agree more often with published CM labels',fontsize=14,weight='bold')
    footer(fig,'Agreement is not independently established accuracy; cohort membership and several analysis choices changed.\nPublished spatial labels were not used to fit the new primary lineage assignments. Ambiguous cells were retained.\nThe spatial outcome analysis also changes vessel geometry and local controls, so its revision cannot be attributed to lineage mapping alone.')
    fig.subplots_adjust(left=.1,right=.97,top=.88,bottom=.19)
    save(fig,'CM_lineage_concordance',pdf)


def main():
    assert json.loads((OUT/'validation.json').read_text())['status']=='PASS'
    with PdfPages(OUT/'whole_heart_Ankrd1_arterial_report.pdf') as pdf:
        summary(pdf)
        for stage in STAGES:stage_maps(stage,pdf)
        sensitivities(pdf);concordance(pdf)
    (OUT/'runtime_versions.json').write_text(json.dumps({'python':sys.version,'recorded_utc':datetime.now(timezone.utc).isoformat(),
        'packages':{k:importlib.metadata.version(k) for k in ['anndata','numpy','pandas','scipy','scikit-learn','matplotlib','shapely']}},indent=2)+'\n')
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='file_manifest.tsv')
    files+=list((ROOT/'analysis').glob('*whole_heart*Ankrd1*vessels.py'))
    files+=[ROOT/'analysis'/n for n in ['reassess_whole_heart_lineages.py','define_heart_arterial_structures.py',
        'test_Ankrd1_arterial_neighborhoods.py','requirements-vessel-analysis.txt','improve_ventricular_alignment.py',
        'plot_xenium_ventricular.py']]
    pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(set(files))]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
    print('REPORT_COMPLETE',OUT,flush=True)


if __name__=='__main__':main()
