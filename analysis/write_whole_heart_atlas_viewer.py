"""Self-contained, offline canvas viewer of all original Xenium cell positions."""
import base64
import gzip
import json

import pandas as pd

from build_whole_heart_spatial_atlas import OUT,SOURCE,STAGES,LINEAGES,REVIEW
from plot_whole_heart_spatial_atlas import COLORS


def main():
    stages={};images={}
    for stage in STAGES:
        obs=pd.read_csv(OUT/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        info=json.loads((SOURCE/stage/'morphology_display.json').read_text())
        stages[stage]={'ids':obs.index.tolist(),'x':obs.x_centroid.tolist(),'y':obs.y_centroid.tolist(),
            'labels':pd.Index(LINEAGES+REVIEW).get_indexer(obs.atlas_identity).tolist(),
            'weights':[None if pd.isna(v) else float(v) for v in obs.lineage_weight],
            'reasons':obs.atlas_review_reason.fillna('').tolist(),
            'extent':info['extent_um']}
        images[stage]='data:image/png;base64,'+base64.b64encode((SOURCE/stage/'DAPI_overview.png').read_bytes()).decode()
    payload={'stages':stages,'names':LINEAGES+REVIEW,'colors':[COLORS[k] for k in LINEAGES+REVIEW],'images':images}
    encoded=base64.b64encode(gzip.compress(json.dumps(payload,separators=(',',':'),allow_nan=False).encode(),mtime=0)).decode()
    html=r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Postnatal heart: all-cell Xenium atlas</title>
<style>
*{box-sizing:border-box}body{font:14px system-ui,sans-serif;color:#202631;margin:0;background:#fafbfc}header{padding:14px 20px;border-bottom:1px solid #ddd;background:white}h1{font-size:22px;margin:0 0 6px}.note{max-width:1050px;color:#535c68;line-height:1.5;font-size:13px}main{display:flex;height:calc(100vh - 125px);min-height:520px}aside{width:285px;flex-shrink:0;padding:16px;overflow:auto;background:white;border-right:1px solid #ddd}label.row{display:flex;align-items:center;gap:7px;margin:9px 0}.dot{width:10px;height:10px;border-radius:50%;display:inline-block;flex-shrink:0}.controls{display:flex;gap:8px;flex-wrap:wrap;margin:8px 0 18px}button,select{font:inherit;padding:6px 9px;border:1px solid #bfc6cf;border-radius:5px;background:white;cursor:pointer}button:hover{background:#eef2f7}#wrap{position:relative;flex:1;min-width:0}canvas{width:100%;height:100%;display:block;touch-action:none;cursor:crosshair}#tip{position:absolute;pointer-events:none;display:none;background:#fffef1;border:1px solid #888;border-radius:4px;padding:8px;white-space:pre-line;font-size:12px;box-shadow:0 2px 6px #bbb;max-width:330px}#status{font-size:12px;color:#5b6570;margin:10px 0;line-height:1.5}#scale{position:absolute;bottom:10px;left:12px;background:#ffffffe6;padding:7px;border-radius:4px;font-size:12px}#loading{position:absolute;top:35%;left:20%;padding:20px;background:#fff;font-size:18px}#details{font-size:12px;line-height:1.5;color:#555;white-space:pre-wrap}a{color:#2764a0}@media(max-width:720px){aside{width:200px;padding:10px}header{padding:10px}h1{font-size:18px}main{height:calc(100vh - 160px)}}
</style>
<header><h1>Postnatal heart: all cardiac cells on Xenium</h1><div class="note">P0–P21 · Native section coordinates · Broad reference-assisted identities, with lineage-anchor and mixed-marker review. Gray profiles remain uncertain. One section per age; these maps do not establish fine subtypes or replicated developmental changes.</div></header>
<main><aside><label for="stage"><b>Stage</b></label> <select id="stage"><option>P0</option><option>P7</option><option>P14</option><option>P21</option></select>
<div class="controls"><button id="reset">Reset view</button><button id="all">All cells</button><button id="noncm">Non-CMs</button><button id="resolved">Resolved</button></div>
<div id="legend"></div><hr><label class="row"><input type="checkbox" id="dapi"> Show DAPI background</label>
<label class="row">Point size <input id="size" type="range" min="1" max="5" value="2"></label>
<button id="csv">Download current-stage annotations</button><div id="status"></div><div class="note">Scroll to zoom; drag to pan. Hover over a cell for its ID and review status. Stage changes reset the view. Neighboring sections are not registered.</div><hr><div id="details">Click a cell to pin its information here.</div></aside>
<div id="wrap"><canvas id="canvas"></canvas><div id="tip"></div><div id="scale"></div><div id="loading">Loading embedded atlas…</div></div></main>
<script>
const encoded='__PAYLOAD__';
(async()=>{try{
const binary=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));
const stream=new Blob([binary]).stream().pipeThrough(new DecompressionStream('gzip'));
const atlas=JSON.parse(await new Response(stream).text());
const canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d'),wrap=document.getElementById('wrap');
const tip=document.getElementById('tip'),status=document.getElementById('status'),stageSelect=document.getElementById('stage');
let stage='P0',data=atlas.stages[stage],enabled=atlas.names.map(()=>true),W=1,H=1,cx=0,cy=0,scale=1,drag=null,hover=-1,raf=0;
const images={};for(const s of Object.keys(atlas.images)){const img=new Image();img.src=atlas.images[s];img.onload=()=>{if(s===stage)requestDraw()};images[s]=img;}
const boxes=[];atlas.names.forEach((name,i)=>{const label=document.createElement('label');label.className='row';const box=document.createElement('input');box.type='checkbox';box.checked=true;box.onchange=()=>{enabled[i]=box.checked;requestDraw()};const dot=document.createElement('span');dot.className='dot';dot.style.background=atlas.colors[i];const text=document.createElement('span');label.append(box,dot,text);document.getElementById('legend').append(label);boxes.push({box,text});});
function counts(){const n=atlas.names.map(()=>0);data.labels.forEach(i=>n[i]++);boxes.forEach((r,i)=>{r.text.textContent=atlas.names[i]+' ('+n[i].toLocaleString()+')';r.box.checked=enabled[i];});}
function reset(){const xmin=Math.min(...data.x.slice(0,10000));let x0=xmin,x1=xmin,y0=data.y[0],y1=y0;for(let i=0;i<data.x.length;i++){x0=Math.min(x0,data.x[i]);x1=Math.max(x1,data.x[i]);y0=Math.min(y0,data.y[i]);y1=Math.max(y1,data.y[i]);}cx=(x0+x1)/2;cy=(y0+y1)/2;scale=Math.min((W-60)/(x1-x0+200),(H-60)/(y1-y0+200));requestDraw();}
function requestDraw(){if(!raf)raf=requestAnimationFrame(draw);}
function draw(){raf=0;ctx.fillStyle='#fff';ctx.fillRect(0,0,W,H);if(document.getElementById('dapi').checked&&images[stage].complete){const e=data.extent;ctx.globalAlpha=.22;ctx.drawImage(images[stage],W/2+(e[0]-cx)*scale,H/2+(Math.min(e[2],e[3])-cy)*scale,(e[1]-e[0])*scale,Math.abs(e[2]-e[3])*scale);ctx.globalAlpha=1;}
 const size=+document.getElementById('size').value;let visible=0;
 // Review profiles first, so resolved lineages remain visible above them.
 const order=[7,8,9,0,1,2,3,4,5,6];for(const group of order){if(!enabled[group])continue;ctx.fillStyle=atlas.colors[group];for(let i=0;i<data.x.length;i++){if(data.labels[i]!==group)continue;const x=W/2+(data.x[i]-cx)*scale,y=H/2+(data.y[i]-cy)*scale;if(x<0||y<0||x>W||y>H)continue;ctx.fillRect(x-size/2,y-size/2,size,size);visible++;}}
 const length=100*scale<80?500:100;ctx.strokeStyle='#111';ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(28,H-45);ctx.lineTo(28+length*scale,H-45);ctx.stroke();ctx.fillStyle='#111';ctx.font='12px system-ui';ctx.fillText(length+' µm',28,H-54);
 status.textContent=stage+': '+data.x.length.toLocaleString()+' original cells; '+visible.toLocaleString()+' displayed in the current view.';document.getElementById('scale').textContent='Native µm coordinates; y increases downward. All points represent original cells.';}
function resize(){const oldW=W,oldH=H;W=wrap.clientWidth;H=wrap.clientHeight;const ratio=window.devicePixelRatio||1;canvas.width=W*ratio;canvas.height=H*ratio;ctx.setTransform(ratio,0,0,ratio,0,0);if(oldW===1)reset();else scale*=Math.min(W/oldW,H/oldH);requestDraw();}
new ResizeObserver(resize).observe(wrap);
function setEnabled(fn){enabled=atlas.names.map((_,i)=>fn(i));counts();requestDraw();}
document.getElementById('all').onclick=()=>setEnabled(()=>true);document.getElementById('noncm').onclick=()=>setEnabled(i=>i>0&&i<7);document.getElementById('resolved').onclick=()=>setEnabled(i=>i<7);document.getElementById('reset').onclick=reset;
document.getElementById('dapi').onchange=requestDraw;document.getElementById('size').oninput=requestDraw;
stageSelect.onchange=()=>{stage=stageSelect.value;data=atlas.stages[stage];hover=-1;tip.style.display='none';counts();reset();};
canvas.addEventListener('wheel',e=>{e.preventDefault();const r=canvas.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;const nx=cx+(x-W/2)/scale,ny=cy+(y-H/2)/scale;scale=Math.max(.02,Math.min(20,scale*Math.exp(-e.deltaY*.0015)));cx=nx-(x-W/2)/scale;cy=ny-(y-H/2)/scale;requestDraw();},{passive:false});
function information(i){return stage+' · '+data.ids[i]+'\n'+atlas.names[data.labels[i]]+'\nX '+data.x[i].toFixed(2)+' µm; Y '+data.y[i].toFixed(2)+' µm\nReference lineage weight: '+(data.weights[i]===null?'not eligible':data.weights[i].toFixed(3))+' (not a calibrated probability)'+(data.reasons[i]?'\nReview: '+data.reasons[i]:'');}
canvas.onpointerdown=e=>{drag={x:e.clientX,y:e.clientY,cx,cy,moved:false};canvas.setPointerCapture(e.pointerId);tip.style.display='none';};
canvas.onpointerup=e=>{if(drag&&!drag.moved&&hover>=0)document.getElementById('details').textContent=information(hover);drag=null;};
canvas.onpointermove=e=>{if(drag){const dx=e.clientX-drag.x,dy=e.clientY-drag.y;if(Math.abs(dx)+Math.abs(dy)>3)drag.moved=true;cx=drag.cx-dx/scale;cy=drag.cy-dy/scale;requestDraw();return;}const r=canvas.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;const nx=cx+(x-W/2)/scale,ny=cy+(y-H/2)/scale;let best=(9/scale)**2,chosen=-1;for(let i=0;i<data.x.length;i++){if(!enabled[data.labels[i]])continue;const d=(data.x[i]-nx)**2+(data.y[i]-ny)**2;if(d<best){best=d;chosen=i;}}hover=chosen;if(chosen<0){tip.style.display='none';return;}tip.textContent=information(chosen);tip.style.display='block';tip.style.left=Math.max(0,Math.min(x+15,W-335))+'px';tip.style.top=Math.max(0,Math.min(y+15,H-155))+'px';};
canvas.onpointerleave=()=>{tip.style.display='none';};
document.getElementById('csv').onclick=()=>{const quote=x=>'"'+String(x).replaceAll('"','""')+'"';const rows=['stage,cell_id,x_um,y_um,atlas_identity,lineage_weight,review_reason'];for(let i=0;i<data.x.length;i++)rows.push([stage,data.ids[i],data.x[i],data.y[i],atlas.names[data.labels[i]],data.weights[i]??'',data.reasons[i]].map(quote).join(','));const url=URL.createObjectURL(new Blob([rows.join('\n')],{type:'text/csv'}));const a=document.createElement('a');a.href=url;a.download=stage+'_all_cardiac_cells.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
counts();resize();document.getElementById('loading').remove();
}catch(error){document.getElementById('loading').textContent='Unable to load atlas: '+error.message+'. Open this downloaded HTML in a current browser supporting DecompressionStream.';}})();
</script></html>'''.replace('__PAYLOAD__',encoded)
    (OUT/'interactive_whole_heart_atlas.html').write_text(html)
    (OUT/'interactive_viewer_data_summary.json').write_text(json.dumps({'stages':{s:{'original_cells':len(d['ids']),'unique_ids':len(set(d['ids']))} for s,d in stages.items()},
        'offline':True,'external_resources':False,'coordinate_rounding':False,'categories':LINEAGES+REVIEW,
        'browser_requirement':'Current browser with DecompressionStream(gzip), canvas and ResizeObserver'},indent=2)+'\n')
    print('OFFLINE_VIEWER_WRITTEN',len(html),'bytes',flush=True)


if __name__=='__main__':main()
