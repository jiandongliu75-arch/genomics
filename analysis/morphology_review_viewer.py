"""Self-contained morphology-only viewer and editable ROI export; no RNA data."""
from pathlib import Path
import json,base64,io
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/morphology_stratified_nonCM'
def main():
 obj=json.loads((OUT/'morphology_ROIs.geojson').read_text());images={}
 for stage in ['P0','P7','P14','P21']:
  p=ROOT/f'results/xenium_ventricular_mapping/{stage}/DAPI_overview.png'
  im=Image.open(p).convert('L');im.thumbnail((1800,1800));b=io.BytesIO();im.save(b,format='JPEG',quality=88)
  info=json.loads(p.with_name('morphology_display.json').read_text())
  images[stage]={'image':'data:image/jpeg;base64,'+base64.b64encode(b.getvalue()).decode(),'width':info['extent_um'][1],'height':info['extent_um'][2]}
 html='''<!doctype html><meta charset="utf-8"><title>Heart morphology ROI review</title>
<style>body{font:15px system-ui;margin:20px;background:#f2f5f7;color:#20384b}main{display:grid;grid-template-columns:minmax(600px,1fr) 320px;gap:16px}svg{background:#000;width:100%;max-height:85vh}button,select,input{font:inherit;margin:5px;padding:6px}aside{background:white;padding:16px}small{display:block;line-height:1.5}.selected{stroke:white!important;stroke-width:12!important}textarea{width:95%;height:70px}</style>
<h2>Review cardiac anatomy using morphology alone</h2><p>No maturation scores, RNA expression or cell identities are displayed. These are conservative candidate ROIs. Unresolved chamber and papillary identities remain unresolved until reviewed.</p>
<label>Section <select id="stage"></select></label><label><input id="overlay" type="checkbox" checked>Show regions</label><button id="export">Download edited GeoJSON</button>
<main><svg id="canvas" xmlns="http://www.w3.org/2000/svg"></svg><aside><h3 id="name">Select a region</h3><small>Click an outline to select it. Drag its white vertices to adjust the boundary. Coordinates are in Xenium micrometers, with y increasing downward.</small><p>Region</p><input id="region"><p>Chamber interpretation</p><select id="chamber"><option>Unresolved</option><option>LV candidate</option><option>RV candidate</option><option>LV expert-confirmed</option><option>RV expert-confirmed</option><option>Other</option></select><p>Confidence</p><select id="confidence"><option>high</option><option>moderate</option><option>low</option><option>none</option></select><p>Review note</p><textarea id="note"></textarea><button id="apply">Apply interpretation</button><small>Edits affect this downloaded review copy. They do not change the frozen analysis masks or recompute results. Expert confirmation should be documented before interpreting these candidates as validated histological compartments.</small></aside></main>
<script>const data=__GEOJSON__,images=__IMAGES__;const svg=document.getElementById('canvas'),stage=document.getElementById('stage');let selected=null,drag=null;
for(const s of Object.keys(images)){let o=document.createElement('option');o.textContent=s;stage.append(o)}
const colors={'Compact wall':'#299b86','Cavity-facing myocardium':'#e6a13b','Papillary/trabecular projection candidate':'#9a66bd','Large-vessel-like surroundings':'#db667c'};
function el(n,attrs){let x=document.createElementNS('http://www.w3.org/2000/svg',n);for(const[k,v]of Object.entries(attrs))x.setAttribute(k,v);return x}
function choose(f){selected=f;document.getElementById('name').textContent=f.properties.ROI_id;document.getElementById('region').value=f.properties.region;document.getElementById('chamber').value=f.properties.anatomical_chamber||'Unresolved';document.getElementById('confidence').value=f.properties.confidence;document.getElementById('note').value=f.properties.review_note||'';render()}
function render(){const m=images[stage.value];svg.setAttribute('viewBox',`0 0 ${m.width} ${m.height}`);svg.replaceChildren();svg.append(el('image',{href:m.image,width:m.width,height:m.height}));if(!document.getElementById('overlay').checked)return;
for(const f of data.features.filter(f=>f.properties.stage===stage.value)){let pts=f.geometry.coordinates[0],c=colors[f.properties.region]||'#cccccc';let p=el('polygon',{points:pts.map(p=>p.join(',')).join(' '),fill:c,'fill-opacity':.17,stroke:c,'stroke-width':7});p.onclick=()=>choose(f);if(f===selected)p.classList.add('selected');svg.append(p);let q=pts.slice(0,-1).reduce((a,p)=>[a[0]+p[0]/(pts.length-1),a[1]+p[1]/(pts.length-1)],[0,0]);let t=el('text',{x:q[0],y:q[1],fill:'white','font-size':m.width/85});t.textContent=f.properties.ROI_id.split('_').slice(1).join('_');t.style.pointerEvents='none';svg.append(t);
if(f===selected)pts.slice(0,-1).forEach((p,i)=>{let v=el('circle',{cx:p[0],cy:p[1],r:m.width/350,fill:'white'});v.onpointerdown=e=>{e.preventDefault();drag={f,i};svg.setPointerCapture(e.pointerId)};svg.append(v)})}}
svg.onpointermove=e=>{if(!drag)return;let p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;p=p.matrixTransform(svg.getScreenCTM().inverse());let a=drag.f.geometry.coordinates[0];a[drag.i]=[p.x,p.y];if(drag.i===0)a[a.length-1]=[p.x,p.y];render()};svg.onpointerup=()=>drag=null;
stage.onchange=()=>{selected=null;render()};document.getElementById('overlay').onchange=render;
document.getElementById('apply').onclick=()=>{if(!selected)return;selected.properties.region=document.getElementById('region').value;selected.properties.anatomical_chamber=document.getElementById('chamber').value;selected.properties.confidence=document.getElementById('confidence').value;selected.properties.review_note=document.getElementById('note').value;selected.properties.reviewed_in_morphology_only_viewer=true;render()};
document.getElementById('export').onclick=()=>{const b=new Blob([JSON.stringify(data,null,2)],{type:'application/geo+json'}),a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='reviewed_heart_morphology_ROIs.geojson';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)};render();</script>'''
 html=html.replace('__GEOJSON__',json.dumps(obj)).replace('__IMAGES__',json.dumps(images))
 (OUT/'Morphology_only_ROI_review.html').write_text(html)
 print('Morphology review HTML written; no RNA/outcome data embedded')
if __name__=='__main__':main()
