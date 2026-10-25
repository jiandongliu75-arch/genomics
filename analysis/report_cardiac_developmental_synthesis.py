"""A biology-first assessment, with cross-stage metrics kept in an appendix."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cardiac-synthesis-mpl')
from cardiac_developmental_synthesis import OUT,SOURCE,ROOT,sha
import textwrap,json,subprocess,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from PIL import Image

plt.rcParams.update({'font.size':11,'pdf.fonttype':42})

def wrap(s,n=111):return textwrap.fill(s,n,break_long_words=False,break_on_hyphens=False)
def page(pdf,title,sections,name):
    fig=plt.figure(figsize=(12,10));fig.text(.06,.94,title,fontsize=20,weight='bold');y=.86
    total_lines=sum(wrap(body).count('\n')+1 for _,body in sections)
    step=min(.0205,(.79-.051*len(sections))/total_lines);fontsize=min(10.8,step*10*72/1.35)
    for heading,body in sections:
        fig.text(.06,y,heading,fontsize=12,weight='bold',va='top');y-=.033
        body=wrap(body);fig.text(.06,y,body,fontsize=fontsize,va='top',linespacing=1.35)
        y-=step*(body.count('\n')+1)+.018
    assert y>=.065,(name,y)
    fig.text(.06,.025,'Sources: Wang et al., Nature Cardiovascular Research (2026), Fig.2–4/Extended Data7; Lu et al., Nature Cell Biology (2021).',fontsize=8)
    fig.savefig(OUT/(name+'.png'),dpi=150,bbox_inches='tight');pdf.savefig(fig);plt.close(fig)

def box(ax,x,y,w,h,title,body,color):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012',facecolor=color,edgecolor='#385263',lw=1);ax.add_patch(p)
    ax.text(x+w/2,y+h-.025,title,ha='center',va='top',fontsize=12,weight='bold')
    ax.text(x+w/2,y+h-.068,textwrap.fill(body,max(22,int(w*115))),ha='center',va='top',fontsize=10,linespacing=1.3)

def developmental_model(pdf):
    fig,ax=plt.subplots(figsize=(12,9));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.02,.98,'A growing ventricle assembles a working tissue unit',fontsize=19,weight='bold',va='top')
    phases=[('Late fetal / P0','Inner-wall vascular expansion is already underway; birth changes circulation and load.'),
        ('P0–P7','Vascular expansion, CM growth and matrix remodeling overlap; inner vascular density catches up.'),
        ('P7–P14','Plexus growth activity declines while CM regulatory programs undergo a marked transition.'),
        ('P14–P21','Established neighborhoods adapt as CMs enlarge and continue functional specialization.')]
    for i,(title,body) in enumerate(phases):box(ax,.02+i*.245,.66,.215,.25,title,body,'#e7f0f5')
    ax.text(.5,.61,'Developmental windows overlap; this is not a measured sequence of individual cells.',ha='center',fontsize=10,color='#555555')
    box(ax,.02,.26,.29,.25,'Cardiomyocytes','Growth and workload increase; structural, metabolic and electrical programs develop.','#f6e1e4')
    box(ax,.68,.26,.29,.25,'Microvasculature','Network expansion and organization support access; EC/perivascular contacts provide candidate cues.','#ddecf5')
    box(ax,.35,.02,.29,.20,'Matrix / supporting cells','Fibroblast remodeling and mural association help organize the changing tissue.','#e4efdf')
    for a,b,label,y in [((.33,.45),(.66,.45),'CM demand and growth signals',.475),((.66,.32),(.33,.32),'Delivery and candidate local cues',.345)]:
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=17,color='#456c7c',lw=1.6));ax.text(.5,y,label,ha='center',fontsize=9)
    ax.add_patch(FancyArrowPatch((.49,.23),(.23,.26),arrowstyle='<->',mutation_scale=15,lw=1.3,color='#607b55'))
    ax.add_patch(FancyArrowPatch((.51,.23),(.77,.26),arrowstyle='<->',mutation_scale=15,lw=1.3,color='#607b55'))
    ax.text(.5,.565,'Shared context: wall geometry, local oxygen demand, loading, nutrition and endocrine maturation',ha='center',fontsize=10)
    fig.text(.04,-.015,'Conceptual synthesis. Timed observations come from the two papers; feedback arrows denote mechanisms to test, not measured causal pathways.',fontsize=8)
    fig.savefig(OUT/'developmental_coordination_model.png',dpi=160,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)

def report():
    with PdfPages(OUT/'cardiac_developmental_biology_report.pdf') as pdf:
        page(pdf,'How sound is the atlas’s biological conclusion?',[
          ('Assessment','The atlas is sound as an exploratory developmental resource and provides a partly validated mechanistic hypothesis. Its strongest conclusion is that CM maturation is spatially coordinated with a developing microvascular environment. It provides experimental support for PTPRM-dependent maturation-associated regulation; it does not establish the complete EC-contact-to-MEIS1 causal chain or a dominant vascular control hierarchy.'),
          ('Why the spatial relationship is persuasive','Within P7 tissue, higher CM maturation scores occur near capillary ECs, proliferating ECs and pericytes. Independent developmental RNA and gene-withholding checks support the transcript-based index. The convergence of these vascular-associated populations points to a developing microvascular unit, although they share anatomy and are not three independent mechanistic tests.'),
          ('What the molecular evidence adds','Spatial communication analysis nominates PTPRM; regulatory inference nominates MEIS1; protein staining places PTPRM at EC–CM interfaces; Ptprm perturbation increases cycling transcripts and reduces Meis1 and maturation-related transcripts. This is stronger than association alone. It still leaves EC-specific necessity, the MEIS1-mediated link and full functional maturation unresolved.'),
          ('The important alternative explanations','Wall remodeling, local workload, cell size and vascular density can change together. A regional association may include both real vascular support and shared developmental drivers. Spatial cells are not independent hearts, graph rings are not fixed physical distances, and transcript indices do not measure force, calcium handling or perfusion. These limitations qualify causal claims rather than invalidate the atlas.'),
          ('The biological model worth retaining','Myocardial growth promotes vascular demand and adaptation; vascular access and local interactions can then support further muscle growth and specialization. Reciprocal coordination is better supported than a one-direction claim that ECs initiate every component of CM maturation.')
        ],'01_atlas_assessment')
        developmental_model(pdf)
        page(pdf,'What could coordinate vascular and CM development?',[
          ('Growth follows a changing ventricular geometry','The perinatal study shows pre-existing coronary vessels entering the developing inner myocardium, beginning before birth, with regional hypoxia and preferential EC proliferation. Inner vessel density approaches outer-wall density by P7. A plausible interpretation is that growing, remodeling muscle increasingly needs intramyocardial supply as diffusion paths and demand change. Higher systemic oxygen and local inner-wall hypoxia can coexist.'),
          ('Permissive support and instructive signals can coexist','Vascular access can support oxidative metabolism, contractile work and growth. Local EC/perivascular/matrix interactions may also change CM regulatory programs. PTPRM is an interface-level candidate requiring competence on both sides. Its rise in CMs after birth suggests that CM responsiveness may be as important as EC transcript abundance. Oxygen-related cell-cycle withdrawal alone is not complete functional maturation.'),
          ('Why cycling ECs can neighbor more mature CMs','EC proliferation expands a network supporting muscle that is leaving the cell cycle and growing by hypertrophy. The two cell types need not share the same proliferative state. A cycling EC can be part of an already integrated, expanding vascular bed. DLL4/Notch restrains and patterns expansion; more proliferation does not automatically mean better vascular delivery.'),
          ('A changing endothelial environment','Coronary ingrowth adds capillary/perivascular neighbors to inner-wall muscle associated with endocardial surfaces. Declining early growth-associated signals and emerging vascular contacts motivate testing a changing mix of endothelial influences and CM responsiveness. This is not conversion of endocardium into coronary ECs. Waning plexus growth around P7–P14 overlaps the atlas’s marked CM regulatory transition; the papers do not prove a strict causal sequence.'),
          ('The other cardiac cells are part of the mechanism','Fibroblast matrix remodeling and mural association can organize a growing vascular–muscle structure and its mechanical environment. Early fibroblast–immune contacts nominate another remodeling relationship. Their functions and signals need direct testing. The atlas should be read as a developing tissue system, not an isolated EC-to-CM pathway.')
        ],'03_coordination_biology')
        page(pdf,'Wormhole: what could cross-stage comparison reveal?',[
          ('The biological questions','Do similar cellular neighborhoods persist while their internal states mature? Does a stage difference reflect expansion of a cell population or a state change within a similar mixture? Do inner and outer myocardial neighborhoods converge during vascular expansion? Are some neighborhoods developmentally advanced or delayed within an age? Do vascular, matrix and CM programs change together or become uncoupled?'),
          ('What a distribution method offers','Wormhole can accelerate comparisons of cell-state distributions, including mixtures hidden by an average. A descriptive developmental map may include CM and non-CM states. A test of non-CM relationships to CM maturation should withhold CM outcomes and examine them afterward. Those are different questions and should use deliberately different representations.'),
          ('What our additional cross-stage check says','The previous within-age benchmark could not decide between-stage usefulness. The new test evaluates 768 reserved cross-stage pairs, including pairs with similar lineage composition. Broad mixture differences are recovered more reliably than subtle differences after composition matching; simple means usually perform better. The current encoder is therefore not validated for the state-level developmental correspondences of greatest interest.'),
          ('What must be represented biologically','The present cloud representation omits internal spatial arrangement and absolute cell abundance. It cannot directly identify EC–CM contact or vascular expansion. Useful stage comparisons need comparable anatomical compartments and separate measures of vascular access, density/contact, cell mixture and cell state. The current composition match does not also match anatomy.'),
          ('What to prioritize','First resolve stage- and region-specific coordination of CM metabolism, structure, electrical development and cell-cycle withdrawal with vascular organization and matrix state. Use distribution comparison to identify conserved, changing or asynchronous neighborhoods. Wormhole is a possible computational tool for that task; it does not supply tracked cells, signaling direction or causal temporal trajectories.')
        ],'04_Wormhole_biological_questions')

def appendix():
    m=pd.read_csv(OUT/'cross_stage_metric_assessment.tsv',sep='\t');pairs=m[['stage_a','stage_b']].drop_duplicates().to_records(index=False)
    fig,axs=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for ax,mode in zip(axs,['random_cross_stage','composition_matched_cross_stage']):
        z=m[m.pair_selection.eq(mode)]
        ax.plot(np.arange(len(z)),z.Wormhole_vs_Sinkhorn_Spearman,'o-',label='Wormhole',c='#b64e50')
        ax.plot(np.arange(len(z)),z.mean_vs_Sinkhorn_Spearman,'s-',label='Mean state features',c='#377eb8')
        ax.set(xticks=np.arange(len(z)),xticklabels=[a+'→'+b for a,b in zip(z.stage_a,z.stage_b)],ylim=(.45,1),ylabel='Rank agreement with held-out Sinkhorn distance',title='Broad mixture differences' if mode=='random_cross_stage' else 'Similar lineage composition');ax.tick_params(axis='x',rotation=25);ax.legend(frameon=False)
    fig.suptitle('Diagnostic appendix: can the existing encoder compare stages?',fontsize=16)
    fig.text(.04,-.06,'64 pairs per stage contrast and selection. No anatomical matching or independent sample replication; metric agreement does not validate biological correspondences.',fontsize=8)
    fig.savefig(OUT/'cross_stage_Wormhole_diagnostic.png',dpi=150,bbox_inches='tight');plt.close(fig)

def validate_report():
    p=OUT/'cardiac_developmental_biology_report.pdf';info=subprocess.check_output(['pdfinfo',str(p)],text=True);assert re.search(r'Pages:\s+4\b',info)
    for f in OUT.glob('*.png'):
        with Image.open(f) as im:im.verify()
    links=0
    for target in re.findall(r'\]\(([^)]+)\)',(OUT/'README.md').read_text()):
        if target.startswith(('http://','https://','#')):continue
        assert (OUT/target).exists(),target;links+=1
    artifacts=[f for f in OUT.rglob('*') if f.is_file() and f.name not in ['file_manifest.tsv','report_validation.json']]
    pd.DataFrame([{'path':str(f.relative_to(ROOT)),'bytes':f.stat().st_size,'sha256':sha(f)} for f in sorted(artifacts)]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
    scripts=[ROOT/'analysis/cardiac_developmental_synthesis.py',ROOT/'analysis/report_cardiac_developmental_synthesis.py']
    (OUT/'report_validation.json').write_text(json.dumps({'status':'PASS','PDF_pages':4,'PNG_integrity_checked':True,'documentation_links_checked':links,
        'script_SHA256':{str(p.relative_to(ROOT)):sha(p) for p in scripts},'scope':'Document/artifact integrity and numerical checks, not biological validation'},indent=2)+'\n')

if __name__=='__main__':report();appendix();validate_report()
