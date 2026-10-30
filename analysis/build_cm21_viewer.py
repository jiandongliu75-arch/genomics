"""Offline interactive map: only the measured21-gene CM score and source labels."""
from pathlib import Path
import csv,gzip,json,math
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM21_primary_spatial'
def rows(p):
 opener=gzip.open if str(p).endswith('.gz') else open
 with opener(p,'rt') as f:yield from csv.DictReader(f,delimiter='\t')
def yes(x):return x=='True'
def main():
 stages={};scores=[]
 coverage={(r['stage'],r['source']):r for r in rows(OUT/'nonCM_source_coverage.tsv')}
 for s in ['P0','P7','P14','P21']:
  cm=[]
  for r in rows(OUT/f'{s}_CM21_scores_and_neighborhoods.tsv.gz'):
   score=float(r['core_score21']);scores.append(score)
   cm.append([r['cell_id'],round(float(r['x_centroid']),1),round(float(r['y_centroid']),1),round(score,6),round(float(r['core_stage_rank'])*100,3),{'Lower quartile':0,'Middle half':1,'Higher quartile':2}[r['core_score_group']]])
  fb={r['cell_id']:r for r in rows(ROOT/f'results/fibroblast_multimethod_identity/{s}_multimethod_FB_annotations.tsv.gz')}
  sources={k[1]:[] for k in coverage if k[0]==s}
  for r in rows(ROOT/f'results/cardiac_spatial_pipeline/{s}/atlas_annotations.tsv.gz'):
   broad=r['atlas_identity']
   if broad not in sources or not yes(r['program_source_eligible']):continue
   xy=[round(float(r['x_centroid']),1),round(float(r['y_centroid']),1)];sources[broad].append(xy)
   if broad=='Endothelial':
    key='EC:'+r['EC_nuclear_type']
    if key in sources:sources[key].append(xy)
   if broad=='Fibroblast':
    a=fb.get(r['cell_id']);key='FB:'+a['recommended_candidate'].replace('-associated','') if a else ''
    if a and yes(a['program_source_eligible']) and key in sources:sources[key].append(xy)
  info=json.loads((ROOT/f'results/xenium_ventricular_mapping/{s}/morphology_display.json').read_text())
  stages[s]=dict(cm=cm,sources=sources,width=info['extent_um'][1],height=info['extent_um'][2],status={k:coverage[(s,k)]['support_status'] for k in sources})
 scores.sort();lo=scores[int(.02*len(scores))];hi=scores[int(.98*len(scores))]
 page='''<!doctype html><meta charset="utf-8"><title>21-gene CM maturation map</title>
<style>body{font:15px system-ui;color:#20384b;margin:20px;background:#f3f6f8}h2{margin-bottom:6px}.controls{display:flex;flex-wrap:wrap;gap:12px;align-items:center}select,button{font:inherit;padding:6px}main{display:grid;grid-template-columns:minmax(500px,1fr) 300px;gap:15px;margin-top:14px}canvas{background:white;width:100%;height:76vh;touch-action:none;border:1px solid #cad4dc}aside{background:white;padding:15px;line-height:1.5}.key{padding:8px;background:#fff;margin:8px 0}.muted{font-size:13px;color:#526779}#bar{display:inline-block;width:180px;height:12px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)}</style>
<h2>Where are higher- and lower-scoring ventricular cardiomyocytes?</h2><p>The sole maturation measure is the directly measured Xenium21-gene score. All QC-eligible ventricular CMs are retained.</p>
<div class="controls"><label>Stage <select id="stage"></select></label><label>CM color <select id="mode"><option value="score">21-gene score, common scale</option><option value="rank">Rank within this stage</option><option value="groups">Lower/higher stage quartiles</option></select></label><label>Neighbor overlay <select id="source"></select></label><button id="reset">Reset view</button><button id="save">Save map PNG</button></div>
<div class="key" id="legend"></div><main><canvas id="map"></canvas><aside><b id="count"></b><p>Wheel to zoom; drag to pan. Click a CM to inspect its score and nearby selected source cells.</p><div id="details">Select a cardiomyocyte.</div><hr><div id="sourceinfo"></div><p class="muted">Magenta outlines show the selected non-CM source. Viewer distances/counts use rounded cell centers; formal analyses use the full-precision tables and tissue-screened neighborhoods. Scores are transcriptional indices; neighbor associations do not establish causation.</p></aside></main>
<script>const dataset=__DATA__,limits=__LIMITS__;const stage=document.getElementById('stage'),mode=document.getElementById('mode'),source=document.getElementById('source'),canvas=document.getElementById('map'),ctx=canvas.getContext('2d');let scale=1,tx=0,ty=0,drag=null,selected=null;
for(const s of Object.keys(dataset)){let o=document.createElement('option');o.textContent=s;stage.append(o)}
const labels=k=>k.replace('EC:','EC · ').replace('FB:','Fibroblast · ').replace('Gfpt2','Gfpt2/Pi16');
function sources(){const old=source.value;source.replaceChildren();let n=document.createElement('option');n.value='';n.textContent='None';source.append(n);for(const k of Object.keys(dataset[stage.value].sources)){let o=document.createElement('option');o.value=k;o.textContent=labels(k);source.append(o)}if(dataset[stage.value].sources[old])source.value=old;else source.value='EC:Capillary-like'}
function color(t){t=Math.max(0,Math.min(1,t));let a=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]],p=t*4,i=Math.min(3,Math.floor(p)),f=p-i;return 'rgb('+a[i].map((v,j)=>Math.round(v*(1-f)+a[i+1][j]*f)).join(',')+')'}
function resize(){const r=canvas.getBoundingClientRect();canvas.width=Math.round(r.width*devicePixelRatio);canvas.height=Math.round(r.height*devicePixelRatio);reset()}
function reset(){let d=dataset[stage.value],w=canvas.width,h=canvas.height;scale=.95*Math.min(w/d.width,h/d.height);tx=(w-d.width*scale)/2;ty=(h-d.height*scale)/2;draw()}
function draw(){const d=dataset[stage.value];ctx.clearRect(0,0,canvas.width,canvas.height);const radius=Math.max(.6,Math.min(3,Math.sqrt(scale)))*devicePixelRatio;for(const c of d.cm){const x=tx+c[1]*scale,y=ty+c[2]*scale;if(x<0||y<0||x>canvas.width||y>canvas.height)continue;ctx.fillStyle=mode.value==='groups'?['#2878b5','#cccccc','#de7928'][c[5]]:color(mode.value==='rank'?c[4]/100:(c[3]-limits[0])/(limits[1]-limits[0]));ctx.fillRect(x-radius/2,y-radius/2,radius,radius)}
let points=d.sources[source.value]||[];ctx.strokeStyle='#e3499b';ctx.lineWidth=devicePixelRatio;for(const p of points){const x=tx+p[0]*scale,y=ty+p[1]*scale;if(x<0||y<0||x>canvas.width||y>canvas.height)continue;ctx.beginPath();ctx.arc(x,y,2*devicePixelRatio,0,Math.PI*2);ctx.stroke()}
if(selected){ctx.strokeStyle='black';ctx.lineWidth=2*devicePixelRatio;ctx.beginPath();ctx.arc(tx+selected[1]*scale,ty+selected[2]*scale,6*devicePixelRatio,0,Math.PI*2);ctx.stroke()}
document.getElementById('count').textContent=stage.value+': '+d.cm.length.toLocaleString()+' ventricular CMs';document.getElementById('legend').innerHTML=mode.value==='groups'?'<span style="color:#2878b5">■ Lower quartile</span> &nbsp; <span style="color:#de7928">■ Higher quartile</span> &nbsp; <span style="color:#999">■ Middle half</span>':'<span id="bar"></span> &nbsp; '+(mode.value==='rank'?'0–100th percentile within this stage':limits.map(x=>x.toFixed(2)).join(' to ')+'; shared21-gene score scale');document.getElementById('sourceinfo').textContent=source.value?labels(source.value)+': '+points.length.toLocaleString()+' eligible source profiles. '+d.status[source.value]:'';details()}
function details(){if(!selected)return;let points=dataset[stage.value].sources[source.value]||[],nearest=Infinity,n=0;for(const p of points){let distance=Math.hypot(p[0]-selected[1],p[1]-selected[2]);nearest=Math.min(nearest,distance);if(distance<=50)n++}document.getElementById('details').textContent='CM '+selected[0]+' ·21-gene score '+selected[3].toFixed(4)+' · stage percentile '+selected[4].toFixed(1)+(source.value?' · nearest '+labels(source.value)+': '+(Number.isFinite(nearest)?nearest.toFixed(1)+' µm':'not assessable')+' · sources within50µm: '+n:'')}
canvas.onpointerdown=e=>{drag={x:e.clientX,y:e.clientY,tx,ty,moved:false};canvas.setPointerCapture(e.pointerId)};canvas.onpointermove=e=>{if(!drag)return;let dx=(e.clientX-drag.x)*devicePixelRatio,dy=(e.clientY-drag.y)*devicePixelRatio;if(Math.hypot(dx,dy)>4)drag.moved=true;tx=drag.tx+dx;ty=drag.ty+dy;draw()};canvas.onpointerup=e=>{if(!drag)return;if(!drag.moved){const r=canvas.getBoundingClientRect(),x=((e.clientX-r.left)*devicePixelRatio-tx)/scale,y=((e.clientY-r.top)*devicePixelRatio-ty)/scale;let best=Infinity,c=null;for(const p of dataset[stage.value].cm){let distance=Math.hypot(p[1]-x,p[2]-y);if(distance<best){best=distance;c=p}}if(best*scale<=15*devicePixelRatio)selected=c;draw()}drag=null};
canvas.addEventListener('wheel',e=>{e.preventDefault();const r=canvas.getBoundingClientRect(),x=(e.clientX-r.left)*devicePixelRatio,y=(e.clientY-r.top)*devicePixelRatio,f=Math.exp(-e.deltaY*.001);tx=x-(x-tx)*f;ty=y-(y-ty)*f;scale*=f;draw()},{passive:false});stage.onchange=()=>{selected=null;sources();reset()};mode.onchange=draw;source.onchange=draw;document.getElementById('reset').onclick=reset;document.getElementById('save').onclick=()=>{let a=document.createElement('a');a.download='CM21_'+stage.value+'_map.png';a.href=canvas.toDataURL();a.click()};window.onresize=resize;sources();resize();</script>'''
 page=page.replace('__DATA__',json.dumps(stages,separators=(',',':'))).replace('__LIMITS__',json.dumps([lo,hi]))
 page=page.replace('Xenium21','Xenium 21').replace('shared21','shared 21').replace('·21','· 21').replace('within50','within 50')
 (OUT/'CM21_spatial_viewer.html').write_text(page);print('21-gene viewer written',len(page),'bytes')
if __name__=='__main__':main()
