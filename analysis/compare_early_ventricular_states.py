"""Descriptive within-age comparisons; nuclei are not independent biological replicates."""
from pathlib import Path
import argparse,json
import h5py,numpy as np,pandas as pd
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('--input',default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');a=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/cardiomyocyte_early_states';out.mkdir(exist_ok=True)
with h5py.File(a.input) as f:
 obs=read_elem(f['obs']);var=read_elem(f['var']);idx=np.flatnonzero(obs.reviewed_identity.eq('Cardiomyocyte').to_numpy());o=obs.iloc[idx].copy()
 print('Reading CM matrices',flush=True)
 x=sparse_dataset(f['X'])[idx,:].tocsr().astype(np.float32)
 raw=sparse_dataset(f['layers/counts_raw'])[idx,:].tocsr().astype(np.float32)
genes=var.gene_symbol.astype(str).to_numpy();ids=var.index.to_numpy()
def normalized(x):
 z=sparse.diags(1e4/np.asarray(x.sum(axis=1)).ravel())@x;z=z.tocsr();z.data=np.log1p(z.data);return z
log=normalized(x);rlog=normalized(raw)
ygenes=['Ddx3y','Kdm5d','Eif2s3y','Uty']
yindex=np.flatnonzero(np.isin(genes,ygenes));xindex=np.flatnonzero(genes=='Xist')
ydet=np.asarray((raw[:,yindex]>=1).sum(axis=1)).ravel();xdet=np.asarray((raw[:,xindex]>=1).sum(axis=1)).ravel()>0
sex=np.full(len(o),'unclassified/ambiguous',dtype=object);sex[xdet&(ydet==0)]='Xist-positive/Y-negative';sex[(~xdet)&(ydet>=2)]='Y-positive/Xist-negative'
o['inferred_sex_marker_group']=sex
o.groupby(['subcluster_id','age','sample','inferred_sex_marker_group'],observed=True).size().rename('nuclei').to_csv(out/'sex_marker_groups_by_cluster_library.tsv',sep='\t')
def mean(z,m):return np.asarray(z[m].mean(axis=0)).ravel()
def frac(z,m):return np.asarray((z[m]>=1).mean(axis=0)).ravel()
quality=[]
for cid,g in o.groupby('subcluster_id',observed=True):
 row={'subcluster_id':cid,'nuclei':len(g),'cycling_positive_n':int(g.cycling_state.eq('cycling-marker positive').sum()),'cycling_positive_fraction':float(g.cycling_state.eq('cycling-marker positive').mean()),'largest_library_fraction':float(g['sample'].value_counts().max()/len(g))}
 for col in ['raw_total_counts','raw_n_genes_by_counts','decontx_contamination','raw_pct_counts_mt','scdblfinder_score']:
  row['median_'+col]=float(g[col].median())
 quality.append(row)
pd.DataFrame(quality).to_csv(out/'cluster_cycling_and_qc.tsv',sep='\t',index=False)
o.groupby(['subcluster_id','age','sample'],observed=True).size().rename('nuclei').to_csv(out/'cluster_age_library_counts.tsv',sep='\t')
comparisons=[]
for age,c1,c2 in [('P0','CAR_1','CAR_2'),('P7','CAR_5','CAR_6')]:
 base=o.age.eq(age)&o.cycling_state.ne('cycling-marker positive');m1=(base&o.subcluster_id.eq(c1)).to_numpy();m2=(base&o.subcluster_id.eq(c2)).to_numpy()
 table=pd.DataFrame({'gene_id':ids,'gene_symbol':genes,'mean_log1p_10k_A':mean(log,m1),'mean_log1p_10k_B':mean(log,m2),'corrected_fraction_ge1_A':frac(x,m1),'corrected_fraction_ge1_B':frac(x,m2),'raw_fraction_ge1_A':frac(raw,m1),'raw_fraction_ge1_B':frac(raw,m2)})
 table['difference_A_minus_B']=table.mean_log1p_10k_A-table.mean_log1p_10k_B
 table['raw_log_difference_A_minus_B']=mean(rlog,m1)-mean(rlog,m2)
 libs=[]
 for lib in o.loc[base,'sample'].unique():
  s1=m1&o['sample'].eq(lib).to_numpy();s2=m2&o['sample'].eq(lib).to_numpy()
  if min(s1.sum(),s2.sum())>=20:
   table[lib+'_difference']=mean(log,s1)-mean(log,s2);libs.append(lib)
 table['concordant_libraries']=sum((np.sign(table[lib+'_difference'])==np.sign(table.difference_A_minus_B)).astype(int) for lib in libs)
 table['eligible_libraries']=len(libs)
 for sexgroup in ['Xist-positive/Y-negative','Y-positive/Xist-negative']:
  sm=o.inferred_sex_marker_group.eq(sexgroup).to_numpy();a1=m1&sm;a2=m2&sm
  if min(a1.sum(),a2.sum())>=30:
   table[sexgroup+'_difference']=mean(log,a1)-mean(log,a2)
   sexlibs=[]
   for lib in o.loc[base,'sample'].unique():
    b1=a1&o['sample'].eq(lib).to_numpy();b2=a2&o['sample'].eq(lib).to_numpy()
    if min(b1.sum(),b2.sum())>=20:
     col=sexgroup+'_'+lib+'_difference';table[col]=mean(log,b1)-mean(log,b2);sexlibs.append(col)
   table[sexgroup+'_eligible_libraries']=len(sexlibs)
   table[sexgroup+'_concordant_libraries']=sum((np.sign(table[col])==np.sign(table.difference_A_minus_B)).astype(int) for col in sexlibs) if sexlibs else 0
   print(age,sexgroup,'n',int(a1.sum()),int(a2.sum()),'libraries',len(sexlibs),flush=True)
 table.to_csv(out/f'{age}_{c1}_vs_{c2}_all_genes.tsv.gz',sep='\t',index=False)
 valid=~table.gene_symbol.str.match(r'^(mt-|Rpl|Rps|Gm\d|\d)')
 for direction,mask in [(c1,valid&(table.difference_A_minus_B>.25)&(table.corrected_fraction_ge1_A>=.2)),(c2,valid&(table.difference_A_minus_B<-.25)&(table.corrected_fraction_ge1_B>=.2))]:
  t=table.loc[mask].copy();t['absolute_difference']=t.difference_A_minus_B.abs();t=t.sort_values('absolute_difference',ascending=False).head(60);t.insert(0,'enriched_cluster',direction);t.insert(0,'age',age);comparisons.append(t)
 print(age,c1,int(m1.sum()),c2,int(m2.sum()),'libraries',libs,flush=True)
pd.concat(comparisons).to_csv(out/'top_descriptive_markers.tsv',sep='\t',index=False)
PANEL=['Mki67','Top2a','Cenpf','Ccnb1','Ube2c','Cdk1','Pcna','Mcm2','Mcm5','Tnni1','Tnni3','Myh6','Myh7','Myl2','Myl3','Ncam1','Txnip','Hmgcs2','Pdk4','Cpt1b','Ppargc1a','Acadm','Cox6a2','Ryr2','Atp2a2','Pln','Casq2','Ankrd1','Nppa','Nppb','Myom2','Myom1','Sntb1','Ddc','Slit2','Xist','Ddx3y','Kdm5d','Eif2s3y','Uty','Fos','Jun','Hspa1a','Acta1','Col1a1','Pecam1','Ptprc']
rows=[]
for cid in sorted(o.subcluster_id.unique(),key=lambda v:int(v.split('_')[1])):
 m=o.subcluster_id.eq(cid).to_numpy();means=mean(log,m);fractions=frac(x,m);rfrac=frac(raw,m)
 for g in PANEL:
  for j in np.flatnonzero(genes==g):rows.append({'subcluster_id':cid,'gene':g,'mean_log1p_10k':means[j],'fraction_corrected_ge1':fractions[j],'fraction_raw_ge1':rfrac[j]})
pd.DataFrame(rows).to_csv(out/'selected_marker_profiles.tsv',sep='\t',index=False)
(out/'methods.json').write_text(json.dumps({'input':a.input,'comparison':'P0 CAR_1 versus CAR_2; P7 CAR_5 versus CAR_6; exclude cycling-marker-positive cells from pairwise comparisons only','expression':'Per-nucleus log1p(corrected counts / total corrected counts * 10000); original integer counts normalized independently for sensitivity check','effect':'Difference of mean log-normalized expression; not log2 fold change','detection':'Fraction with at least one corrected or original integer count','library_checks':'Within-library expression differences calculated when both groups have at least 20 retained nuclei; descriptive, not animal-level inference','sex_sensitivity':'Original integer counts: female-marker-like if Xist >=1 and no detected Y markers; male-marker-like if Xist undetected and >=2 of Ddx3y/Kdm5d/Eif2s3y/Uty detected. Others unclassified/ambiguous. Inferred marker groups, not verified animal sex. Matched-marker-group contrasts require >=30 cells per cluster; library contrasts >=20 each.', 'limitations':'Only two libraries per age; animal independence unverified. No significance tests. Cycling markers do not establish completed cell division.'},indent=2))
print('Finished descriptive comparisons',flush=True)
