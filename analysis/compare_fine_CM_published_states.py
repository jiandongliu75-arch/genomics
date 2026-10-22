"""Describe agreement and lineage conflicts with published final labels."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from push_ventricular_spatial_resolution import OUT,STAGES

CM=['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM','a.CM','av.CM']


def main():
    overlaps=[];compatibility=[]
    for stage in STAGES:
        path=OUT/stage/'fine_CM_spatial_annotations.tsv.gz'
        obs=pd.read_csv(path,sep='\t',index_col=0,float_precision='round_trip')
        obs['published_CM_compatible']=obs.author_atlas_cell_state.isin(CM)
        obs['published_nonCM_discordant']=~obs.published_CM_compatible&obs.author_atlas_cell_state.ne('Not present in published atlas')
        obs['prioritized_and_published_CM_compatible']=obs.prioritized_fine_candidate&obs.published_CM_compatible
        obs.to_csv(path,sep='\t',na_rep='NA')
        for subset,mask in [('All proposals',obs.included_in_weighted_mapping),('Prioritized',obs.prioritized_fine_candidate)]:
            table=pd.crosstab(obs.loc[mask,'fine_group_proposal'],obs.loc[mask,'author_atlas_cell_state'])
            table.to_csv(OUT/stage/('prioritized_vs_published_states.tsv' if subset=='Prioritized' else 'fine_groups_vs_published_states.tsv'),sep='\t')
            for group,row in table.iterrows():
                total=int(row.sum());top=row.sort_values(ascending=False).head(4)
                compat=int(row.reindex(CM,fill_value=0).sum())
                record={'stage':stage,'subset':subset,'group':group,'cells':total,'published_CM_cells':compat,
                    'published_CM_fraction':compat/total,'top_author_states':'; '.join(f'{k} {v} ({v/total:.1%})' for k,v in top.items())}
                for name in ['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM']:
                    record[name+'_fraction_all']=int(row.get(name,0))/total
                    record[name+'_fraction_given_published_CM']=int(row.get(name,0))/compat if compat else np.nan
                overlaps.append(record)
            compatibility.append({'stage':stage,'subset':subset,'cells':int(mask.sum()),
                'published_CM_cells':int((mask&obs.published_CM_compatible).sum()),
                'published_nonCM_cells':int((mask&obs.published_nonCM_discordant).sum()),
                'absent_from_published_atlas':int((mask&obs.author_atlas_cell_state.eq('Not present in published atlas')).sum())})
    overlaps=pd.DataFrame(overlaps);overlaps.to_csv(OUT/'published_state_overlap_summary.tsv',sep='\t',index=False)
    pd.DataFrame(compatibility).to_csv(OUT/'published_broad_lineage_compatibility.tsv',sep='\t',index=False)
    part=overlaps[overlaps.subset.eq('Prioritized')].sort_values(['stage','group'])
    labels=['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM','Atrial/AV','Published non-CM/absent']
    colors=['#6489ad','#e2903b','#269e9c','#334c81','#49a257','#bf79aa','#b8b8b8']
    fig,ax=plt.subplots(figsize=(13,7));left=np.zeros(len(part))
    values=[]
    for row in part.itertuples():
        obs=pd.read_csv(OUT/row.stage/'prioritized_vs_published_states.tsv',sep='\t',index_col=0).loc[row.group]
        fractions=[obs.get(x,0)/row.cells for x in CM[:5]]
        fractions.append((obs.get('a.CM',0)+obs.get('av.CM',0))/row.cells)
        fractions.append(1-sum(fractions));values.append(fractions)
    values=np.array(values)
    for j,(label,color) in enumerate(zip(labels,colors)):
        ax.barh(np.arange(len(part)),values[:,j],left=left,color=color,label=label,height=.7);left+=values[:,j]
    ax.set_yticks(np.arange(len(part)),[f'{r.group} (n={r.cells:,})' for r in part.itertuples()]);ax.invert_yaxis()
    ax.set_xlim(0,1);ax.set_xlabel('Fraction of all prioritized cells, including lineage conflicts')
    ax.set_title('Our prioritized fine placements versus the published final state labels',fontsize=15,weight='bold')
    ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),frameon=False,fontsize=9)
    fig.text(.5,.035,'These are annotation-overlap proportions, not accuracy estimates. Gray explicitly retains cells the authors labeled non-CM or omitted.\n'
        'Author labels were not used for fitting or selecting our candidates. Agreement is strongest for Slit2/conduction-associated groups; lineage conflicts require review.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.16,right=.77,bottom=.17,top=.9)
    fig.savefig(OUT/'comparison_with_published_states.png',dpi=180,facecolor='white')
    fig.savefig(OUT/'comparison_with_published_states.pdf',facecolor='white');plt.close(fig)
    print(part[['stage','group','cells','published_CM_fraction','top_author_states']].to_string(index=False))


if __name__=='__main__':main()
