"""Endothelial reference competition and nuclear/boundary sensitivity, all ages.

Reuses the already completed, image-guided whole-section segmentation. No
droplet-specific decontamination model or artificial spatial doublet prior.
"""
import argparse
import anndata as ad
import scanpy as sc
from sklearn.metrics import precision_recall_fscore_support, adjusted_rand_score
from endothelial_common import *
from improve_ventricular_alignment import neighbors

SOURCE=ROOT/'results/xenium_ventricular_mapping'
ATLAS=ROOT/'results/whole_heart_spatial_atlas'
SEG=ROOT/'results/xenium_nuclear_anchored_segmentation'

def norm(x):
    return np.log2(1+1e6*x/np.maximum(x.sum(axis=1),1)[:,None]).astype(np.float32)

def scale(x,mu,sd):
    f=np.clip((x-mu)/sd,-5,5)
    return (f/np.maximum(np.linalg.norm(f,axis=1,keepdims=True),1e-8)).astype(np.float32)

class Projection:
    def __init__(self,stage,qgenes,qnative,qglobal,ann):
        meta=pd.read_csv(SOURCE/'reference_metadata.tsv.gz',sep='\t',index_col=0)
        with np.load(SOURCE/'reference_shared_panel_counts.npz') as b:
            assert np.array_equal(b['cell_ids'].astype(str),meta.index)
            self.genes=b['genes'].astype(str)
            take=meta.age.eq(stage).to_numpy() & meta.reviewed_identity.ne('Unresolved/mixed').to_numpy()
            take &= (meta.reviewed_identity.ne('Endothelial') | meta.index.isin(ann.index)).to_numpy()
            rm=meta.loc[take].copy();raw=b['original'][take].astype(np.float32);corrected=b['corrected'][take].astype(np.float32)
        self.position=pd.Index(qgenes).get_indexer(self.genes);assert (self.position>=0).all()
        # CM program is diagnostic only and cannot define a new endothelial subtype.
        self.use=~np.isin(self.genes,SEX+LINEAGES['Cardiomyocyte']+['Ankrd1','Creb5','Tlr4','Atf3'])
        self.libraries=sorted(rm['sample'].unique())
        self.classes=FAMILIES+[UNRESOLVED]
        ec=rm.index.isin(ann.index)
        self.rmeta=ann.loc[rm.index[ec]].copy()
        self.rlabels=self.rmeta.EC_family.to_numpy(dtype=str)
        self.ec_corrected=corrected[ec]
        self.models={};self.broad={}
        for name,counts,mask in [('primary',corrected,np.ones(len(rm),bool)),('raw_reference',raw,np.ones(len(rm),bool))]+[
                (lib,corrected,rm['sample'].eq(lib).to_numpy()) for lib in self.libraries]:
            # Fine classification uses all eligible same-age nuclei (no quotas).
            sel=mask & ec
            rn=norm(counts[sel][:,self.use])
            labels=ann.loc[rm.index[sel],'EC_family'].to_numpy(dtype=str)
            self.models[name]=(scale(rn,rn.mean(axis=0),np.maximum(rn.std(axis=0),.25)),labels)
            # Broad competition retains non-EC alternatives. Library x previous
            # RNA subtype centroids summarize all nuclei; centroid counts do not
            # act as priors: score each lineage by its best matching centroid.
            bn=norm(counts[mask][:,self.use]);bm=rm.loc[mask].copy()
            mu,sd=bn.mean(axis=0),np.maximum(bn.std(axis=0),.25)
            feat=scale(bn,mu,sd)
            bm['prototype']=bm['sample'].astype(str)+'::'+bm.reviewed_identity.astype(str)+'::'+bm.proposed_subtype.astype(str)
            centroids=[];blabels=[]
            for _,ii in bm.groupby('prototype',observed=True).indices.items():
                if len(ii)<5:continue
                c=feat[ii].mean(axis=0);c/=max(np.linalg.norm(c),1e-8)
                centroids.append(c);blabels.append(str(bm.iloc[ii[0]].reviewed_identity))
            self.broad[name]=(np.asarray(centroids,dtype=np.float32),np.asarray(blabels))
        self.calibrate(qnative,qglobal)

    def calibrate(self,q,qglobal=None):
        a=norm(q[:,self.position][:,self.use])
        self.mu,self.sd=a.mean(axis=0),np.maximum(a.std(axis=0),.25)
        if qglobal is not None:
            b=norm(qglobal[:,self.position][:,self.use])
            self.bmu,self.bsd=b.mean(axis=0),np.maximum(b.std(axis=0),.25)

    def predict(self,q):
        qf=scale(norm(q[:,self.position][:,self.use]),self.mu,self.sd)
        bq=scale(norm(q[:,self.position][:,self.use]),self.bmu,self.bsd)
        predictions={};broad={}
        for name,(r,labels) in self.models.items():
            ix,w,sim=neighbors(r,qf)
            p=np.column_stack([(w*(labels[ix]==fam)).sum(axis=1) for fam in self.classes])
            predictions[name]=(p,sim)
            bc,bl=self.broad[name];score=bq@bc.T
            ecscore=score[:,bl=='Endothelial'].max(axis=1)
            other=score[:,bl!='Endothelial'].max(axis=1)
            broad[name]=(ecscore>=other+.05,ecscore-other,bl[score.argmax(axis=1)])
        return predictions,broad

    def benchmark(self,stage):
        rows=[]
        for lib in self.libraries:
            q=self.rmeta['sample'].eq(lib).to_numpy();r=~q
            a,b=norm(self.ec_corrected[r][:,self.use]),norm(self.ec_corrected[q][:,self.use])
            rf,qf=scale(a,a.mean(axis=0),np.maximum(a.std(axis=0),.25)),scale(b,b.mean(axis=0),np.maximum(b.std(axis=0),.25))
            ix,w,_=neighbors(rf,qf)
            p=np.column_stack([(w*(self.rlabels[r][ix]==c)).sum(axis=1) for c in self.classes])
            pred=np.asarray(self.classes)[p.argmax(axis=1)]
            prec,rec,f1,n=precision_recall_fscore_support(self.rlabels[q],pred,labels=self.classes,zero_division=0)
            for c,a,b,f,k in zip(self.classes,prec,rec,f1,n):
                rows.append({'stage':stage,'query_library':lib,'family':c,'precision':a,'recall':b,'f1':f,'nuclei':k})
        return rows

def annotate(q,genes,model,pred,broad,nuclear=False):
    p=pred['primary'][0];order=np.argsort(p,axis=1)
    top=np.asarray(model.classes)[order[:,-1]]
    margin=p[np.arange(len(p)),order[:,-1]]-p[np.arange(len(p)),order[:,-2]]
    raw=np.asarray(model.classes)[pred['raw_reference'][0].argmax(axis=1)]
    both=np.ones(len(q),bool)
    for lib in model.libraries:both &= np.asarray(model.classes)[pred[lib][0].argmax(axis=1)]==top
    markers=np.zeros(len(q),int);anchors=np.zeros(len(q),int);supported=np.zeros(len(q),bool)
    for f in FAMILIES:
        m=top==f;a,b,ok=marker_support(q,genes,f)
        markers[m],anchors[m],supported[m]=a[m],b[m],ok[m]
    ev=lineage_evidence(q,genes)
    if nuclear:ev['EC_lineage_supported']=ev.Endothelial_genes.ge(2)&ev.Endothelial_fraction.ge(.03)
    foreign=np.column_stack([ev['strong_'+c+'_program'] for c in LINEAGES if c!='Endothelial']).any(axis=1)
    ev['mixed_lineage_review']=ev.EC_lineage_supported & foreign
    fit=q[:,model.position][:,model.use]
    informative=(fit.sum(axis=1)>=20)&((fit>=1).sum(axis=1)>=8)
    broadpass=informative & ev.EC_lineage_supported.to_numpy()
    for v in broad.values():broadpass &= v[0]
    strict=informative & (p.max(axis=1)>=.6) & (margin>=.15) & both & (raw==top) & supported & (top!=UNRESOLVED) & broadpass & ~ev.mixed_lineage_review.to_numpy()
    d=ev.copy()
    d['family_proposal']=top;d['family_weight']=p.max(axis=1);d['family_margin']=margin
    d['both_library_agreement']=both;d['raw_reference_agreement']=raw==top
    d['marker_genes']=markers;d['anchor_genes']=anchors;d['information_QC']=informative
    d['broad_EC_competition_pass']=broadpass
    d['broad_EC_cosine_margin']=broad['primary'][1]
    d['broad_lineage_proposal']=broad['primary'][2]
    d['fine_supported']=strict;d['mapped_family']=np.where(strict,top,UNRESOLVED)
    for k,c in enumerate(model.classes):d['weight__'+c]=p[:,k]
    d['total_gene_molecules']=q.sum(axis=1);d['detected_genes']=(q>=1).sum(axis=1)
    return d

def cluster_views(stage,ids,counts,frames,genes,folder,extra_eligible=None):
    eligible=frames['nuclear_anchored'].broad_EC_competition_pass & ~frames['nuclear_anchored'].mixed_lineage_review
    if extra_eligible is not None:eligible &= np.asarray(extra_eligible)
    ii=np.flatnonzero(eligible);o=frames['nuclear_anchored'].loc[ids[ii]].copy()
    q=ad.AnnData(sparse.csr_matrix((len(o),0)),obs=o)
    use=~np.isin(genes,WITHHOLD)
    assert use.sum()>=100
    sensitivity=[]
    for variant in ['nuclear_anchored','high_confidence','nuclear_only','nuclear_depth_matched']:
        z=lognorm(counts[variant][ii][:,use])
        a=ad.AnnData(z,var=pd.DataFrame(index=genes[use]))
        sc.pp.scale(a,zero_center=True,max_value=10)
        sc.pp.pca(a,n_comps=20,svd_solver='arpack',random_state=42)
        q.obsm['X_'+variant]=a.obsm['X_pca']
        sc.pp.neighbors(q,n_neighbors=20,use_rep='X_'+variant,key_added=variant,random_state=42)
        sc.tl.leiden(q,neighbors_key=variant,resolution=.5,random_state=42,flavor='igraph',n_iterations=3,key_added=variant+'_cluster')
    sc.tl.umap(q,neighbors_key='nuclear_anchored',random_state=42,min_dist=.4,maxiter=250)
    for seed in [7,123]:
        sc.tl.leiden(q,neighbors_key='nuclear_anchored',resolution=.5,random_state=seed,flavor='igraph',n_iterations=3,key_added='seed_'+str(seed))
    for res in [.25,.8]:
        sc.tl.leiden(q,neighbors_key='nuclear_anchored',resolution=res,random_state=42,flavor='igraph',n_iterations=3,key_added='resolution_'+str(res))
    label=q.obs.nuclear_anchored_cluster.to_numpy(dtype=str)
    z=lognorm(counts['nuclear_anchored'][ii]);rows=[];mr=[]
    for cl in sorted(set(label),key=int):
        m=label==cl;mean=z[m].mean(axis=0);rest=z[~m].mean(axis=0)
        frac=(counts['nuclear_anchored'][ii[m]]>=1).mean(axis=0)
        rank=np.flatnonzero((mean-rest>.2)&(frac>=.15)&use)
        rank=rank[np.argsort((mean-rest)[rank])[::-1]][:30]
        top=genes[rank[:2]];alias='/'.join(top)+' EC program' if len(top) else 'Marker-poor EC program'
        cid=stage+'_XEC'+cl
        vc=q.obs.loc[m,'mapped_family'].value_counts();best=vc.index[0];share=vc.iloc[0]/m.sum()
        row={'stage':stage,'cluster_id':cid,'gene_name':alias,'profiles':int(m.sum()),'most_common_individual_label':best,
            'most_common_label_fraction':share,'marker_genes':';'.join(genes[rank[:12]])}
        for key in ['high_confidence_cluster','nuclear_only_cluster','nuclear_depth_matched_cluster','seed_7','seed_123','resolution_0.25','resolution_0.8']:
            other=q.obs[key].to_numpy(dtype=str)
            row[key+'_best_jaccard']=max(np.sum(m&(other==c))/np.sum(m|(other==c)) for c in set(other))
        q.obs.loc[m,'Xenium_cluster_id']=cid;q.obs.loc[m,'Xenium_gene_name']=alias
        rows.append(row)
        for j in rank:mr.append({'stage':stage,'cluster_id':cid,'gene':genes[j],'mean_log_difference':mean[j]-rest[j],'fraction_ge1':frac[j]})
    for key in ['high_confidence_cluster','nuclear_only_cluster','nuclear_depth_matched_cluster','seed_7','seed_123','resolution_0.25','resolution_0.8']:
        sensitivity.append({'stage':stage,'comparison':key,'adjusted_rand_index':adjusted_rand_score(label,q.obs[key]),'clusters':q.obs[key].nunique()})
    q.obs['Xenium_UMAP1'],q.obs['Xenium_UMAP2']=q.obsm['X_umap'][:,0],q.obsm['X_umap'][:,1]
    q.write_h5ad(folder/'EC_de_novo_embeddings_metadata_only.h5ad',compression='gzip')
    q.obs.to_csv(folder/'EC_de_novo_clusters.tsv.gz',sep='\t',na_rep='NA')
    save(pd.DataFrame(rows),folder/'EC_de_novo_cluster_summary.tsv')
    save(pd.DataFrame(mr),folder/'EC_de_novo_cluster_markers.tsv')
    save(pd.DataFrame(sensitivity),folder/'EC_de_novo_cluster_sensitivity.tsv')
    return q.obs

def main(stage):
    folder=OUT/'Xenium'/stage;folder.mkdir(parents=True,exist_ok=True)
    atlas=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    # Include previously uncertain EC proposals; no cell is automatically rescued.
    mask=(atlas.atlas_identity.eq('Endothelial')|atlas.whole_reference_lineage_proposal.eq('Endothelial')).to_numpy()
    ids=atlas.index[mask]
    data=ad.read_h5ad(SEG/stage/'nuclear_anchored_counts.h5ad',backed='r')
    native=ad.read_h5ad(SOURCE/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad',backed='r')
    np.testing.assert_array_equal(data.obs_names,atlas.index);np.testing.assert_array_equal(native.obs_names,atlas.index)
    np.testing.assert_array_equal(data.var_names,native.var_names)
    genes=data.var_names.to_numpy(dtype=str);assert len(genes)==403
    global_native=native.X[:].toarray().astype(np.float32)
    global_nuclear=data.layers['nuclear_only'].toarray().astype(np.float32)
    values={'native':global_native[mask],'nuclear_anchored':data.X[mask].toarray().astype(np.float32)}
    for key in ['high_confidence','depth_matched','nuclear_only','nuclear_depth_matched']:
        values[key]=data.layers[key][mask].toarray().astype(np.float32)
    geometry=data.obs.loc[ids,'geometry_unresolved_native_counts_retained'].to_numpy()
    data.file.close();native.file.close()
    values['cytoplasmic_only']=values['nuclear_anchored']-values['nuclear_only'];assert (values['cytoplasmic_only']>=0).all()
    ann=pd.read_csv(OUT/'snRNA/EC_annotations.tsv.gz',sep='\t',index_col=0)
    model=Projection(stage,genes,values['native'],global_native,ann)
    del global_native
    save(pd.DataFrame(model.benchmark(stage)),folder/'opposite_RNA_library_benchmark.tsv')
    frames={};summary=[]
    for variant,counts in values.items():
        print('XENIUM_EC_VIEW',stage,variant,len(ids),flush=True)
        pred,broad=model.predict(counts)
        f=annotate(counts,genes,model,pred,broad,variant=='nuclear_only');f.index=ids
        f['geometry_unresolved_review']=geometry
        f.loc[geometry,'fine_supported']=False;f.loc[geometry,'mapped_family']=UNRESOLVED
        f.to_csv(folder/f'EC_{variant}_assignments.tsv.gz',sep='\t',na_rep='NA')
        frames[variant]=f
    # Nuclear-specific assay scaling is a sensitivity view, not independent truth.
    informative=frames['nuclear_only'].information_QC.to_numpy()
    if informative.sum()>=30:model.calibrate(values['nuclear_only'][informative],global_nuclear)
    del global_nuclear
    pred,broad=model.predict(values['nuclear_only'])
    f=annotate(values['nuclear_only'],genes,model,pred,broad,True);f.index=ids
    f['geometry_unresolved_review']=geometry;f.loc[geometry,'fine_supported']=False;f.loc[geometry,'mapped_family']=UNRESOLVED
    f.to_csv(folder/'EC_nuclear_calibrated_assignments.tsv.gz',sep='\t',na_rep='NA');frames['nuclear_calibrated']=f
    for key,f in frames.items():
        for family in FAMILIES+[UNRESOLVED]:
            summary.append({'stage':stage,'variant':key,'family':family,'candidate_EC_profiles':len(f),'profiles':f.mapped_family.eq(family).sum(),
                'information_QC':f.information_QC.sum(),'broad_EC_competition_pass':f.broad_EC_competition_pass.sum(),'mixed_lineage_review':f.mixed_lineage_review.sum()})
    save(pd.DataFrame(summary),folder/'view_assignment_summary.tsv')
    out=atlas.loc[ids,['x_centroid','y_centroid','cell_area','nucleus_area','atlas_identity','whole_reference_lineage_proposal','author_atlas_cell_state']].copy()
    out['geometry_unresolved_review']=geometry
    for key,f in frames.items():
        for col in ['mapped_family','fine_supported','mixed_lineage_review','broad_EC_competition_pass','total_gene_molecules','detected_genes','Cardiomyocyte_fraction','strong_Cardiomyocyte_program','Endothelial_genes']:
            out[key+'__'+col]=f[col]
    primary=frames['nuclear_anchored'];nuc=frames['nuclear_calibrated'];old=frames['native'];core=frames['high_confidence']
    strong_labels=pd.concat([frames[k].mapped_family.rename(k).where(frames[k].fine_supported) for k in ['native','nuclear_anchored','high_confidence','nuclear_only','nuclear_calibrated']],axis=1)
    conflict=strong_labels.nunique(axis=1)>1
    mixed=primary.mixed_lineage_review|nuc.mixed_lineage_review
    out['fine_consensus_family']=UNRESOLVED
    out['confidence_tier']='Fine subtype unresolved'
    new=primary.fine_supported
    nuclear=new & nuc.fine_supported & primary.mapped_family.eq(nuc.mapped_family)
    # Both nuclear calibrations must agree for the highest-confidence release.
    nuclear &= frames['nuclear_only'].fine_supported & primary.mapped_family.eq(frames['nuclear_only'].mapped_family)
    stable=old.fine_supported & old.mapped_family.eq(primary.mapped_family)
    coreagree=core.fine_supported & core.mapped_family.eq(primary.mapped_family)
    tiers=[(new & stable & coreagree,'C: Cell-boundary stable; nuclear unconfirmed'),
           (nuclear & ~old.fine_supported,'B: New cell+nuclear supported'),
           (nuclear & stable,'A: Native, new cell and nuclear concordant')]
    for mask_t,tier in tiers:
        out.loc[mask_t,'fine_consensus_family']=primary.loc[mask_t,'mapped_family'];out.loc[mask_t,'confidence_tier']=tier
    out.loc[conflict,'confidence_tier']='Conflicting supported subtypes; review'
    out.loc[mixed,'confidence_tier']='Mixed EC/foreign program; review'
    out.loc[geometry,'confidence_tier']='Unresolved projected nuclear geometry; review'
    out.loc[conflict|mixed|geometry,'fine_consensus_family']=UNRESOLVED
    out['high_confidence_subtype']=out.confidence_tier.str.startswith(('A:','B:'))
    out['assigned_including_moderate']=out.confidence_tier.str.startswith(('A:','B:','C:'))
    # Broad evidence and fine subtype uncertainty are separate fields.
    out['EC_broad_clean']=primary.broad_EC_competition_pass & ~mixed & ~geometry
    out['native_EC_broad']=out.atlas_identity.eq('Endothelial')
    out['new_broad_EC_candidate']=~out.native_EC_broad & out.EC_broad_clean
    out['boundary_sensitive_CM_program']=old.strong_Cardiomyocyte_program & ~core.strong_Cardiomyocyte_program
    out['nuclear_supported_CM_program']=nuc.strong_Cardiomyocyte_program & nuc.EC_lineage_supported
    out['nuclear_fraction']=values['nuclear_only'].sum(axis=1)/np.maximum(values['nuclear_anchored'].sum(axis=1),1)
    # Depth-controlled molecular-program comparison. Thinning changes detection
    # and variance; persistence/decline is diagnostic evidence, not proof of source.
    program=[]
    for key in values:
        ev=frames[key]
        for cohort,m in [('all_candidates',np.ones(len(ids),bool)),('native_strong_CM_program',old.strong_Cardiomyocyte_program.to_numpy())]:
            if not m.any():continue
            program.append({'stage':stage,'cohort':cohort,'view':key,'profiles':int(m.sum()),
                'mean_CM_fraction':ev.loc[m,'Cardiomyocyte_fraction'].mean(),'median_CM_fraction':ev.loc[m,'Cardiomyocyte_fraction'].median(),
                'mean_EC_fraction':ev.loc[m,'Endothelial_fraction'].mean(),'strong_CM_fraction':ev.loc[m,'strong_Cardiomyocyte_program'].mean(),
                'median_molecules':ev.loc[m,'total_gene_molecules'].median(),'median_genes':ev.loc[m,'detected_genes'].median()})
    save(pd.DataFrame(program),folder/'mixed_RNA_depth_control_comparison.tsv')
    clusters=cluster_views(stage,ids,values,frames,genes,folder)
    for c in ['Xenium_cluster_id','Xenium_gene_name','Xenium_UMAP1','Xenium_UMAP2']:
        out[c]=clusters[c].reindex(out.index)
    out.to_csv(folder/'EC_consensus_annotations.tsv.gz',sep='\t',na_rep='NA')
    counts=ad.AnnData(sparse.csr_matrix(values['nuclear_anchored']),obs=out.copy(),var=pd.DataFrame(index=genes))
    for key in ['native','high_confidence','nuclear_only','nuclear_depth_matched','depth_matched']:
        counts.layers[key]=sparse.csr_matrix(values[key])
    counts.uns['source']='Existing whole-section DAPI-anchored, acquired boundary-channel watershed; original profiles retained; QV>=20 gene transcripts; no RNA-based boundary tuning.'
    counts.write_h5ad(folder/'EC_candidate_count_views.h5ad',compression='gzip')
    write_json(folder/'complete.json',{'stage':stage,'candidate_profiles':len(out),'native_broad_EC':int(out.native_EC_broad.sum()),
        'broad_clean_after_segmentation':int(out.EC_broad_clean.sum()),'new_broad_candidates':int(out.new_broad_EC_candidate.sum()),
        'high_confidence_subtypes':int(out.high_confidence_subtype.sum()),'assigned_including_moderate':int(out.assigned_including_moderate.sum()),
        'mixed_review':int(mixed.sum()),'supported_subtype_conflicts':int(conflict.sum()),'coordinates_unchanged':True,
        'segmentation_source_sha256':sha(SEG/stage/'nuclear_anchored_counts.h5ad'),
        'nuclear_calibration':'Separate nuclear scaling among information-eligible candidates; primary nuclear view retains frozen native calibration',
        'spatial_doublets':'Mixed/merged profiles assessed by competing lineage programs and ownership sensitivity, not scDblFinder/Scrublet',
        'ambient_RNA':'No droplet ambient model or global spatial decontamination; compare actual nuclear, cytoplasmic, boundary-core counts plus exact molecule-depth controls',
        'de_novo_clustering':'403 measured genes minus listed foreign-lineage/sex/outcome diagnostics; within each section, no spatial coordinates or RNA subtype labels in the graph; resolution .5,20PC,20NN,seed42',
        'broad_competition':'Best same-age RNA subtype x library centroid per lineage, cosine EC margin>=.05 and consistent corrected/raw/both library models; marker support and information gates also required',
        'fine_gates':'soft30 weighted RNA neighbors; weight>=.6,margin>=.15,raw/both-library agreement,>=2 subtype markers including>=1 anchor; lymphatic needs both Prox1 and Flt4;>=20 fitting molecules,>=8 genes; no mixed flag; broad competition pass',
        'consensus':'A/B: new cell + both nuclear calibrations, A also native; C native+new+boundary core; supported subtype conflicts and strong mixed programs excluded; all candidate profiles retained in audit.'})
    print('XENIUM_EC_COMPLETE',stage,out.confidence_tier.value_counts().to_dict(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=STAGES);a=p.parse_args();main(a.stage)
