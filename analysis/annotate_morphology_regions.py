"""Conservative, outcome-independent image ROIs, frozen before new model fits.

Only acquired DAPI and physical extents are read. These are analyst-drawn
candidate regions, not expert-validated histological masks or a recovery of
author LV/RV masks. No RNA, cell identity or maturation field enters this step.
"""
from pathlib import Path
import json,hashlib,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/anatomy-mpl')
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Circle
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/morphology_stratified_nonCM'
STAGES=['P0','P7','P14','P21']
# Normalized image coordinates were drawn while inspecting DAPI alone.
ROIS={
 'P0':[
  ('C1','Compact wall','Main ventricular mass','high',[(.48,.17),(.64,.15),(.76,.18),(.80,.26),(.73,.30),(.59,.26),(.46,.26)]),
  ('C2','Compact wall','Main ventricular mass','high',[(.19,.39),(.27,.32),(.32,.37),(.25,.50),(.25,.60),(.20,.64),(.16,.55)]),
  ('C3','Compact wall','Main ventricular mass','high',[(.77,.30),(.82,.33),(.82,.48),(.76,.59),(.69,.62),(.68,.55),(.75,.46)]),
  ('C4','Compact wall','Main ventricular mass','high',[(.29,.65),(.39,.68),(.53,.70),(.61,.66),(.66,.69),(.59,.76),(.44,.77),(.31,.73)]),
 ],
 'P7':[
  ('C1','Compact wall','Main ventricular mass','high',[(.46,.18),(.55,.15),(.64,.17),(.71,.22),(.69,.27),(.59,.25),(.50,.24)]),
  ('C2','Compact wall','Main ventricular mass','high',[(.73,.31),(.81,.40),(.87,.51),(.90,.63),(.87,.72),(.80,.72),(.79,.60),(.74,.47),(.69,.38)]),
  ('C3','Compact wall','Main ventricular mass','high',[(.43,.73),(.55,.77),(.70,.79),(.80,.76),(.83,.82),(.71,.86),(.53,.86),(.43,.81)]),
  ('C4','Compact wall','Main ventricular mass','high',[(.32,.46),(.38,.46),(.39,.58),(.42,.68),(.40,.75),(.34,.68),(.31,.56)]),
  ('T1','Cavity-facing myocardium','Main cavity margin','moderate',[(.42,.27),(.47,.25),(.52,.25),(.53,.28),(.49,.32),(.45,.35),(.41,.34)]),
  ('T2','Cavity-facing myocardium','Main cavity margin','moderate',[(.48,.34),(.54,.32),(.59,.34),(.63,.36),(.66,.35),(.68,.39),(.59,.41),(.52,.39)]),
  ('P1','Papillary/trabecular projection candidate','Main cavity margin','low',[(.565,.305),(.600,.300),(.620,.325),(.608,.348),(.579,.336)]),
 ],
 'P14':[
  ('C1','Compact wall','Image-left ventricular wall','high',[(.18,.39),(.23,.31),(.28,.28),(.29,.32),(.25,.43),(.25,.51),(.21,.52),(.17,.47)]),
  ('C2','Compact wall','Image-right ventricular wall','high',[(.50,.11),(.65,.08),(.77,.10),(.85,.14),(.84,.19),(.73,.17),(.61,.16),(.52,.19)]),
  ('C3','Compact wall','Image-right ventricular wall','high',[(.80,.24),(.88,.24),(.90,.38),(.88,.51),(.83,.58),(.76,.59),(.78,.48),(.79,.36)]),
  ('C4','Compact wall','Image-right ventricular wall','high',[(.58,.60),(.67,.63),(.77,.60),(.82,.57),(.83,.63),(.74,.68),(.62,.70),(.55,.67)]),
  ('T1','Cavity-facing myocardium','Image-left cavity margin','moderate',[(.30,.34),(.34,.33),(.37,.40),(.36,.47),(.33,.52),(.30,.49),(.31,.42)]),
  ('T2','Cavity-facing myocardium','Image-right cavity margin','moderate',[(.52,.40),(.55,.35),(.58,.34),(.60,.39),(.57,.44),(.55,.49),(.51,.53),(.48,.51)]),
  ('T3','Cavity-facing myocardium','Image-right cavity margin','moderate',[(.65,.23),(.70,.20),(.73,.22),(.68,.28),(.64,.33),(.61,.34),(.62,.28)]),
  ('P1','Papillary/trabecular projection candidate','Image-right cavity margin','low',[(.535,.445),(.565,.440),(.576,.467),(.559,.493),(.535,.482)]),
 ],
 'P21':[
  ('C1','Compact wall','Image-left ventricular mass','high',[(.13,.27),(.20,.20),(.24,.24),(.20,.38),(.20,.51),(.24,.60),(.20,.64),(.14,.57),(.09,.43)]),
  ('C2','Compact wall','Superior ventricular wall','high',[(.33,.09),(.47,.055),(.59,.055),(.70,.10),(.72,.16),(.62,.14),(.49,.12),(.36,.16)]),
  ('C3','Compact wall','Image-right ventricular mass','high',[(.82,.28),(.88,.32),(.91,.43),(.93,.57),(.90,.64),(.86,.65),(.86,.55),(.84,.43),(.79,.36)]),
  ('C4','Compact wall','Image-left ventricular mass','high',[(.25,.65),(.34,.69),(.42,.68),(.48,.66),(.48,.72),(.40,.76),(.30,.75),(.23,.70)]),
  ('T1','Cavity-facing myocardium','Image-right cavity margin','moderate',[(.74,.54),(.79,.53),(.82,.58),(.81,.63),(.77,.66),(.72,.64),(.74,.59)]),
  ('T2','Cavity-facing myocardium','Central basal cavity margin','moderate',[(.51,.65),(.55,.64),(.58,.67),(.59,.72),(.55,.75),(.52,.71)]),
 ]}
VESSELS={'P0':[],'P7':[(.55,.171,45),(.755,.339,25)],'P14':[(.50,.615,27),(.415,.225,25)],'P21':[(.435,.415,35),(.485,.355,35)]}
COLORS={'Compact wall':'#299b86','Cavity-facing myocardium':'#e6a13b','Papillary/trabecular projection candidate':'#9a66bd','Large-vessel-like surroundings':'#db667c'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True);(OUT/'figures').mkdir(exist_ok=True)
 features=[];sources=[];fig,axs=plt.subplots(2,2,figsize=(12,13),layout='constrained')
 for ax,stage in zip(axs.ravel(),STAGES):
  path=ROOT/f'results/xenium_ventricular_mapping/{stage}/DAPI_overview.png'
  info_path=path.with_name('morphology_display.json');info=json.loads(info_path.read_text());image=np.array(Image.open(path).convert('L'))
  extent=info['extent_um'];width=extent[1];height=extent[2]
  ax.imshow(image,cmap='gray',extent=extent,vmin=0,vmax=210)
  sources.append(dict(stage=stage,DAPI_path=str(path.relative_to(ROOT)),DAPI_sha256=sha(path),display_metadata_sha256=sha(info_path),extent_um=extent))
  for rid,region,territory,confidence,xy in ROIS[stage]:
   points=np.array(xy)*[width,height];points=np.vstack([points,points[0]])
   properties=dict(stage=stage,ROI_id=stage+'_'+rid,region=region,image_territory=territory,confidence=confidence,
    anatomical_chamber='Unresolved',histological_trabecular_verified=False,papillary_verified=False,
    primary_eligible=confidence in ['high','moderate'],annotation_basis='DAPI architecture only; no expression or cell labels viewed during ROI construction')
   features.append(dict(type='Feature',properties=properties,geometry=dict(type='Polygon',coordinates=[points.tolist()])))
   ax.add_patch(Polygon(points,closed=True,facecolor=COLORS[region],edgecolor=COLORS[region],alpha=.28,lw=1))
   center=points[:-1].mean(0);ax.text(*center,rid,color='white',fontsize=8,ha='center',bbox=dict(facecolor='black',alpha=.5,pad=1))
  for i,(x,y,lumen_radius) in enumerate(VESSELS[stage],1):
   center=np.array([x*width,y*height]);radius=lumen_radius+50
   theta=np.linspace(0,2*np.pi,65);points=center+radius*np.column_stack([np.cos(theta),np.sin(theta)])
   properties=dict(stage=stage,ROI_id=f'{stage}_V{i}',region='Large-vessel-like surroundings',image_territory='Image-defined vessel-like lumen',
    confidence='low',anatomical_chamber='Unresolved',primary_eligible=False,lumen_center_um=center.tolist(),
    lumen_radius_um=lumen_radius,outer_radius_um=radius,vessel_identity_verified=False,
    annotation_basis='Closed DAPI-poor lumen and surrounding tissue; cannot establish artery/vein or exclude nonvascular void')
   features.append(dict(type='Feature',properties=properties,geometry=dict(type='Polygon',coordinates=[points.tolist()])))
   ax.add_patch(Circle(center,radius,fill=False,edgecolor=COLORS['Large-vessel-like surroundings'],lw=1.3));ax.text(*center,'V'+str(i),color='white',fontsize=8)
  ax.set(title=stage,xlabel='x (µm)',ylabel='y (µm)')
 fig.suptitle('Morphology-only annotations: conservative wall and cavity-facing regions',fontsize=15)
 fig.savefig(OUT/'figures/Figure_1_morphology_regions.png',dpi=240);fig.savefig(OUT/'figures/Figure_1_morphology_regions.pdf');plt.close(fig)
 obj=dict(type='FeatureCollection',coordinate_system='Xenium micrometers; y increases downward',features=features)
 p=OUT/'morphology_ROIs.geojson';p.write_text(json.dumps(obj,indent=2)+'\n')
 design=dict(mask_sha256=sha(p),source_images=sources,expression_inputs_used=False,outcomes_loaded_by_annotation_script=False,
  outcome_independent_not_independently_blinded_expert_review=True,
  primary_regions=['Compact wall','Cavity-facing myocardium'],primary_region_meanings={'Compact wall':'Conservative continuous dense ventricular wall patches',
   'Cavity-facing myocardium':'Visible cavity margin/ridged wall; not definitively histological trabecular myocardium'},
  chamber_policy='Image-coordinate territories retained; definitive LV/RV not assigned without orientation and landmarks',
  papillary_policy='Low-confidence papillary/trabecular projection candidates; excluded from primary tests',
  vessel_policy='DAPI-only vessel-like lumens; no artery/vein verification; sensitivity only',
  unannotated_areas='Unclassified rather than forced into compact/trabecular/chamber classes',
  annotation_is_original_author_mask=False,new_biological_replicates=0,
  unavailable_images='Boundary and alphaSMA/vimentin cache missing; public Zenodo retrieval returned proxy403 even after network permission was granted')
 (OUT/'anatomy_freeze.json').write_text(json.dumps(design,indent=2)+'\n')
 print('FROZEN_MORPHOLOGY_MASK',design['mask_sha256'],len(features),'ROIs',flush=True)
if __name__=='__main__':main()
