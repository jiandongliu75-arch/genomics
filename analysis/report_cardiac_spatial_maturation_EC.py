"""Figures and synthesis for the frozen cardiac atlas and EC-neighbor analysis."""
from cardiac_spatial_maturation_EC import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import textwrap

PALETTE = {'Cardiomyocyte':'#ce6a79','Endothelial':'#377eb8','Fibroblast':'#37a36c'}
ECPAL = dict(zip(FAMILIES, ['#0072b2','#e69f00','#cc3311','#009e73','#cc79a7','#c5c5c5']))
FBPAL = dict(zip(['Postn-associated','Col8a1-associated','Gfpt2-associated','Lsamp-associated','Kirrel3-associated'],
                ['#e69f00','#0072b2','#009e73','#cc79a7','#8b5c38']))
NAMES = {'measured_nuclear_cycling':'Measured nuclear cycling module',
         'measured_nuclear_cycling_positive':'Fraction of ECs with ≥2 nuclear cycling genes',
         'near_cycling_EC_30um':'Within 30 µm of a ≥2-marker nuclear cycling EC',
         'inferred_proliferation':'RNA-inferred proliferation module',
         'inferred_sprouting':'RNA-inferred six-gene sprouting module',
         'inferred_missing5':'RNA-inferred unmeasured-five component',
         'measured_nuclear_guidance':'Measured nuclear guidance proxy (5 genes)'}

def legend(ax, items, size=8):
    ax.legend(handles=[Line2D([],[],ls='',marker='o',markersize=4,color=c,label=k) for k,c in items],
              loc='upper left', bbox_to_anchor=(1,1), frameon=False, fontsize=size)

def finish_map(ax, atlas):
    ax.set_aspect('equal'); ax.invert_yaxis(); ax.set(xlabel='X (µm)',ylabel='Y (µm)')
    ax.set_xlim(atlas.x_centroid.min()-50, atlas.x_centroid.max()+50)
    ax.set_ylim(atlas.y_centroid.max()+50, atlas.y_centroid.min()-50)
    ax.tick_params(labelsize=8)

def points(ax, d, color, size=.8, **kwargs):
    if isinstance(color,tuple): color=matplotlib.colors.to_hex(color)
    return ax.scatter(d.x_centroid, d.y_centroid, c=color, s=size, linewidths=0, rasterized=True, **kwargs)

def background(ax, atlas):
    points(ax,atlas,'#ededed',.15)

def write_figure(fig,name,pdf):
    fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight')
    pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def spatial(stage,pdf):
    atlas=read(f'results/whole_heart_spatial_atlas/{stage}/all_cell_annotations.tsv.gz')
    cm=read(f'results/cardiac_spatial_maturation_EC/{stage}/ventricular_CM_maturation.tsv.gz')
    ec=read(f'results/cardiac_spatial_maturation_EC/{stage}/EC_activity_annotations.tsv.gz')
    fb=read(f'results/fibroblast_spatial_consensus/{stage}/fibroblast_consensus_annotations.tsv.gz')
    evidence=read('results/ventricular_spatial_resolution/fine_group_evidence_summary.tsv')
    fig,axes=plt.subplots(2,3,figsize=(23,13),layout='constrained')
    for ax in axes.ravel(): background(ax,atlas)
    ax=axes[0,0]
    for t,c in PALETTE.items(): points(ax,atlas.loc[atlas.atlas_identity.eq(t)],c,.6)
    legend(ax,PALETTE.items());ax.set_title('Three principal lineages · frozen broad identities')
    ax=axes[0,1]
    m=points(ax,cm,cm.maturation_21gene_mean,1,vmin=1.8,vmax=3.7,cmap='viridis')
    fig.colorbar(m,ax=ax,shrink=.55,label='Measured 21-gene mean log1p(CP10k)')
    ax.set_title(f'Ventricular CM maturation · {len(cm):,} strict profiles')
    ax=axes[0,2];points(ax,cm,'#c8c8c8',.5)
    groups=sorted(cm.loc[cm.fine_candidate_supported,'fine_group_proposal'].dropna().unique())
    colors=plt.get_cmap('tab10').colors
    items=[]
    for i,g in enumerate(groups):
        z=cm.loc[cm.fine_candidate_supported&cm.fine_group_proposal.eq(g)]
        points(ax,z,colors[i%10],1.5)
        name=evidence.loc[g,'mapping_family'] if g in evidence.index else 'RNA candidate'
        items.append((g+': '+name,colors[i%10]))
    legend(ax,items,7);ax.set_title('Prioritized fine ventricular candidates\nGray CMs retain unresolved fine identity')
    ax=axes[1,0];points(ax,atlas.loc[atlas.atlas_identity.eq('Fibroblast')],'#c8c8c8',.8)
    for t,c in FBPAL.items(): points(ax,fb.loc[fb.high_confidence_family.eq(t)],c,2)
    legend(ax,FBPAL.items());ax.set_title('Fibroblast families · strict cell+nuclear support')
    ax=axes[1,1]
    for t in ['Unresolved']+FAMILIES[:-1]: points(ax,ec.loc[ec.EC_type.eq(t)],ECPAL[t],1.5)
    legend(ax,ECPAL.items());ax.set_title('Eligible ECs · nuclear-supported five types\nSeparate from stricter whole-cell subtype calls')
    ax=axes[1,2]
    points(ax,ec,'#a6bddb',1)
    points(ax,ec.loc[ec.measured_nuclear_cycling_positive.eq(1)],'#d73027',3)
    legend(ax,[('Eligible EC','#a6bddb'),('≥2 nuclear cycling genes','#d73027')])
    ax.set_title('Directly measured EC proliferation-marker support\nMki67 / Top2a / Cenpf; not confirmed division')
    for ax in axes.ravel(): finish_map(ax,atlas)
    fig.suptitle(stage+' · cardiac lineages, supported subtypes and ventricular maturation\nNative section coordinates; gray fine identities are retained; confidence-selected coverage is incomplete',fontsize=17)
    write_figure(fig,stage+'_cardiac_spatial_overview',pdf)

def effects(pdf):
    a=pd.read_csv(OUT/'association_coefficients.tsv',sep='\t')
    primary=a.loc[a.radius_um.eq(50)&a.outcome.eq('maturation_21gene_mean')&a.model.eq('full')&a.block_um.eq(400)]
    fig,axes=plt.subplots(1,2,figsize=(16,8),layout='constrained',sharex=True)
    panels=[['measured_nuclear_cycling','measured_nuclear_cycling_positive','inferred_proliferation','near_cycling_EC_30um'],
            ['inferred_sprouting','inferred_missing5','measured_nuclear_guidance']]
    stages=STAGES;colors=['#0072b2','#e69f00','#009e73','#cc3311']
    for ax,keys in zip(axes,panels):
        for j,s in enumerate(stages):
            d=primary.loc[primary.stage.eq(s)].set_index('exposure').reindex(keys)
            y=np.arange(len(keys))+(j-1.5)*.16
            ax.errorbar(d.beta_ySD_per_exposureSD,y,xerr=[d.beta_ySD_per_exposureSD-d.block_CI_low,d.block_CI_high-d.beta_ySD_per_exposureSD],
                        fmt='o',ms=4,color=colors[j],capsize=2,label=s)
        ax.axvline(0,color='#999999',lw=1);ax.set_yticks(range(len(keys)),[NAMES[k] for k in keys],fontsize=9)
        ax.invert_yaxis();ax.set_xlabel('Adjusted maturation SD per EC-exposure SD\nNegative = lower panel maturation near greater activity')
        ax.grid(axis='x',alpha=.2);ax.legend(frameon=False)
    axes[0].set_title('Proliferation · direct and inferred evidence')
    axes[1].set_title('Sprouting-associated signals · exploratory inference')
    fig.suptitle('50-µm EC neighborhoods · within-section adjusted associations\n95% intervals resample 400-µm spatial blocks; these are not animal-replicate intervals',fontsize=15)
    write_figure(fig,'EC_activity_CM_maturation_associations',pdf)
    fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    for ax,stage in zip(axes.ravel(),STAGES):
        d=read(f'results/cardiac_spatial_maturation_EC/{stage}/CM_EC_neighborhoods_50um.tsv.gz')
        d=d.loc[d.model_eligible].copy()
        b=design(d);base=np.column_stack([np.ones(len(d)),b.to_numpy(float)])
        y=d.maturation_21gene_mean.to_numpy();res=y-base@np.linalg.lstsq(base,y,rcond=1e-9)[0]
        for e,c in [('measured_nuclear_cycling','#0072b2'),('inferred_sprouting','#d55e00')]:
            val=d['local_'+e].to_numpy();er=val-base@np.linalg.lstsq(base,val,rcond=1e-9)[0]
            bins=pd.qcut(er,5,labels=False,duplicates='drop')
            mean=[np.mean(res[bins==i]) for i in range(5)]
            ax.plot(range(1,6),mean,'o-',color=c,label=NAMES[e])
        ax.axhline(0,color='#aaaaaa',lw=1);ax.set(xlabel='Quintile of activity residual after covariate adjustment',ylabel='Mean maturation residual (log-expression units)',title=stage)
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('CM maturation across local EC activity gradients\nDescriptive residual bins; overlapping neighborhoods are not independent replicates',fontsize=15)
    write_figure(fig,'maturation_across_EC_activity_gradients',pdf)

def maturity_and_context(pdf):
    a=pd.read_csv(OUT/'maturation_stage_summary.tsv',sep='\t')
    e=pd.read_csv(OUT/'EC_measured_cycling_summary.tsv',sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for key,label in [('maturation_21gene_mean','21-gene mean'),('maturation_matched_control','21-gene mean minus fixed matched controls')]:
        d=a.loc[a.score.eq(key)].set_index('stage').loc[STAGES]
        axes[0].plot(STAGES,d['mean'],'o-',label=label)
    axes[0].set(title='Ventricular maturation panel across sections',ylabel='Mean score');axes[0].legend(fontsize=9,frameon=False)
    for key,label in [('measured_nuclear_cycling_positive','Nuclear ≥2 genes'),('measured_native_cycling_positive','Native whole profile ≥2 genes'),('measured_core_cycling_positive','Boundary core ≥2 genes')]:
        d=e.loc[e.score.eq(key)].set_index('stage').loc[STAGES]
        axes[1].plot(STAGES,100*d.fraction,'o-',label=label)
    axes[1].set(title='EC proliferation-marker support',ylabel='% of selected eligible EC profiles');axes[1].legend(fontsize=9,frameon=False)
    fig.suptitle('Across-age trends and within-age proximity are different questions\nOne deposited section per age; capture and selection affect these descriptive fractions',fontsize=14)
    write_figure(fig,'maturation_and_EC_cycling_by_stage',pdf)
    context=pd.read_csv(OUT/'subtype_lineage_neighborhoods_30um.tsv',sep='\t')
    context=context.loc[context.lineage.isin(['Fibroblast','EC strict','EC nuclear-supported'])&~context.sparse_lt30].copy()
    context['row']=context.stage+' '+context.lineage+' '+context.subtype
    cols=['fraction_'+k for k in ['Cardiomyocyte','Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']]
    fig,ax=plt.subplots(figsize=(12,max(6,.27*len(context)+2)),layout='constrained')
    im=ax.imshow(context[cols],vmin=0,vmax=.65,cmap='YlGnBu',aspect='auto')
    ax.set_yticks(range(len(context)),context.row,fontsize=8)
    ax.set_xticks(range(len(cols)),[k.replace('fraction_','') for k in cols],rotation=30,ha='right')
    fig.colorbar(im,ax=ax,shrink=.6,label='Mean resolved-neighbor fraction (30 µm)')
    ax.set_title('Supported subtype neighborhoods · groups with ≥30 eligible profiles\nEC strict and nuclear-supported are distinct evidence views, not additive populations',fontsize=13)
    write_figure(fig,'subtype_neighborhood_composition',pdf)

def findings(pdf):
    a=pd.read_csv(OUT/'near_cycling_EC_maturation_contrasts.tsv',sep='\t').set_index('stage').loc[STAGES]
    m=pd.read_csv(OUT/'maturation_stage_summary.tsv',sep='\t')
    m=m.loc[m.score.eq('maturation_21gene_mean')].set_index('stage').loc[STAGES]
    fig=plt.figure(figsize=(13,9))
    fig.text(.06,.93,'Findings: spatial organization and the proximity hypothesis',fontsize=19,weight='bold')
    fig.text(.06,.865,'P0/P7: no convincing reduction in the overall 21-gene maturation score near cycling ECs.',fontsize=13,weight='bold')
    ax=fig.add_axes([.06,.54,.88,.25]);ax.axis('off')
    rows=[]
    for s in STAGES:
        r=a.loc[s]
        rows.append([s,f'{m.loc[s,"mean"]:.3f}',f'{int(r.near_CMs):,}',f'{int(r.other_CMs):,}',
            f'{r.adjusted_near_minus_other:+.4f}',f'[{r.block_CI_low:+.4f}, {r.block_CI_high:+.4f}]'])
    t=ax.table(cellText=rows,colLabels=['Stage','Mean 21-gene score','Near cycling EC','Other CMs','Adjusted near−other','95% spatial-block interval'],
               loc='center',cellLoc='center',colWidths=[.07,.17,.14,.14,.19,.29])
    t.auto_set_font_size(False);t.set_fontsize(10);t.scale(1,2)
    lines=[
      'Near: CM centroid≤30µm from an eligible EC with≥2 distinct nuclear Mki67/Top2a/Cenpf genes. Comparison: other eligible CMs in the same section; both groups have≥3 eligible ECs within 50µm. Effects use log-expression score units.',
      'Direct nuclear cycling, inferred proliferation and relative distance to cycling versus marker-negative ECs give no convincing overall immaturity association at P0/P7. Boundary-core and fixed-control score checks support that conclusion.',
      'Early inferred sprouting has a tiny positive association with maturation (~0.02 outcome SD/exposure SD), not the predicted negative direction. Missing-gene recovery is weak; this is not a verified sprouting niche.',
      'An exploratory P0 sarcomere-only component is modestly lower with greater local nuclear cycling (~−0.045 SD); the full panel and P7 do not show that pattern. It is not evidence for globally less mature CMs.',
      'P21 has a small negative near-cycling contrast in native counts. It attenuates with CM boundary-core counts and is uncertain in the author-compatible subset; this is boundary-sensitive rather than a firm late niche.',
      'Spatially: capillary-like ECs are interspersed with CMs; endocardial-like ECs favor lining-associated territories; arterial-like ECs have many more smooth-muscle neighbors. Small venous/lymphatic candidates remain provisional.',
      'Early Postn fibroblasts occupy CM-containing interstitial territories. Later Col8a1 and Gfpt2 families have distinct stromal/vascular contexts; separate P21 RNA substates remain unresolved spatially. Conduction-program CM candidates lie closer to endocardial proxies.',
      'Across-age maturation increases while EC cycling-marker support decreases. That developmental co-trend does not establish a within-age local relationship or causality. One section/age limits generalization.'
    ]
    y=.53
    for line in lines:
        w=textwrap.fill(line,145);fig.text(.06,y,w,fontsize=9.5,va='top');y-=.022*(w.count('\n')+1)+.011
    pdf.savefig(fig);plt.close(fig)

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    with PdfPages(OUT/'cardiac_spatial_maturation_EC_report.pdf') as pdf:
        fig=plt.figure(figsize=(12,8));fig.text(.07,.92,'Cardiac spatial organization and CM maturation',fontsize=21,weight='bold')
        lines=[
          'Question: are ventricular CMs near proliferating or sprouting-associated ECs less mature, particularly at P0/P7?',
          'Three principal lineages: cardiomyocytes, fibroblasts and endothelial cells; all other cardiac profiles remain in the background.',
          'Maturation: measured mean log1p(CP10k) of the exact 21 Figure4b genes, with fixed-control, Myl2-withheld and boundary-core sensitivities.',
          'EC proliferation: directly measured Mki67/Top2a/Cenpf nuclear RNA, multi-marker support and the separate RNA-inferred proliferation module.',
          'Sprouting: exploratory RNA-inferred six-gene module; only Angpt2 of that signature is measured by Xenium. The missing-five component is tracked separately.',
          'A new alignment sensitivity excludes all 21 CM maturation genes and all sprouting targets from EC predictors and their normalization.',
          'Primary neighborhood: at least3 eligible ECs within 50µm and at least5 resolved cardiac neighbors;30/100µm are sensitivities.',
          'Adjustment: CM capture/area, CM cycling/Slit2, EC capture/CM-marker fraction, EC subtype composition, local lineage density/composition, XY splines and vascular/endocardial proxies.',
          'Spatial uncertainty:400µm block bootstrap;200/800µm sensitivities; buffered spatial prediction. No independent animal-level significance is claimed.',
          'Confidence-qualified fibroblast and EC families are displayed without assigning unsupported fine states or treating gray profiles as absent.',
          'One Xenium section per age; centroid proximity does not demonstrate contact, causation, active sprouting or completed cell division.'
        ]
        y=.84
        for line in lines:
            wrapped=textwrap.fill(line,116);fig.text(.07,y,wrapped,fontsize=10,va='top');y-=.028*(wrapped.count('\n')+1)+.013
        pdf.savefig(fig);plt.close(fig)
        findings(pdf)
        for s in STAGES: spatial(s,pdf)
        effects(pdf);maturity_and_context(pdf)

if __name__=='__main__': main()
