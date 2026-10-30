"""Recover acquired morphology overviews without retaining large source TIFFs.

ZIP members are fully streamed/CRC checked. Sparse temporary TIFFs retain only
the header and tail; all decoded overview tiles must lie within retained bytes.
No expression matrices or cell labels are loaded by this acquisition step.
"""
from pathlib import Path
import sys,json,struct,zlib,argparse,os
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import requests,tifffile
from remotezip import RemoteZip
from PIL import Image
from acquire_xenium_boundary_inputs import discover,ranged
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/morphology_stratified_nonCM'

def acquire(stage,channel):
 folder=OUT/'morphology_inputs';folder.mkdir(parents=True,exist_ok=True)
 output=folder/f'{stage}_channel{channel}_overview.npz';receipt=output.with_suffix('.json')
 if output.exists() and receipt.exists():return output
 archive,record=discover(stage);url=archive['links']['self'];name=f'morphology_focus_000{channel}.ome.tif'
 with RemoteZip(url) as z:
  infos=[x for x in z.infolist() if not x.filename.startswith('__MACOSX/') and Path(x.filename).name==name]
  assert len(infos)==1;info=infos[0]
 header=ranged(requests.Session(),url,info.header_offset,info.header_offset+29)
 fn,extra=struct.unpack_from('<HH',header,26);start=info.header_offset+30+fn+extra
 assert info.compress_type in [0,8]
 decoder=zlib.decompressobj(-15) if info.compress_type==8 else None
 temp=folder/f'.{stage}_{channel}.sparse.tif';keep_head=65536;keep_tail=48<<20;size=info.file_size
 positions=list(range(0,info.compress_size,16<<20));written=0;crc=0
 print('MORPHOLOGY_STREAM',stage,channel,info.compress_size,flush=True)
 try:
  with temp.open('wb') as target,ThreadPoolExecutor(max_workers=4) as pool:
   target.truncate(size)
   def fetch(pos):return ranged(requests.Session(),url,start+pos,start+min(pos+(16<<20),info.compress_size)-1)
   pending={pos:pool.submit(fetch,pos) for pos in positions[:4]}
   for k,pos in enumerate(positions):
    compressed=pending.pop(pos).result()
    if k+4<len(positions):pending[positions[k+4]]=pool.submit(fetch,positions[k+4])
    plain=decoder.decompress(compressed) if decoder else compressed
    end=written+len(plain);crc=zlib.crc32(plain,crc)
    for lo,hi in [(0,min(keep_head,size)),(max(keep_head,size-keep_tail),size)]:
     a=max(written,lo);b=min(end,hi)
     if b>a:target.seek(a);target.write(plain[a-written:b-written])
    written=end
   if decoder:assert decoder.eof and not decoder.unused_data and decoder.flush()==b''
  assert written==size and crc==info.CRC
  with tifffile.TiffFile(temp,is_ome=False) as f:
   page=f.pages[0].pages[2] # full-resolution / 8
   assert all(a>=size-keep_tail and a+b<=size for a,b in zip(page.dataoffsets,page.databytecounts))
   image=page.asarray(maxworkers=2);shape=f.pages[0].shape
  np.savez_compressed(output,image=image,pixel_um=np.array(1.7),full_shape_yx=np.array(shape))
  receipt.write_text(json.dumps(dict(stage=stage,channel=channel,channel_name='ATP1A1/CD45/E-Cadherin' if channel==1 else 'alphaSMA/Vimentin',
   record=record,archive=archive['key'],member=info.filename,url=url,source_bytes=size,compressed_bytes=info.compress_size,
   stream_CRC32=crc,expected_CRC32=info.CRC,whole_stream_CRC_verified=True,overview_level=3,
   overview_shape_yx=list(image.shape),overview_pixel_um=1.7,all_decoded_tiles_in_retained_tail=True,
   expression_or_maturation_inputs_used=False),indent=2)+'\n')
  print('MORPHOLOGY_OVERVIEW_VERIFIED',stage,channel,image.shape,flush=True)
 finally:
  if temp.exists():temp.unlink()
 return output

def display(stage):
 from matplotlib import pyplot as plt
 folder=OUT/'morphology_inputs';(OUT/'figures').mkdir(parents=True,exist_ok=True)
 raw=next(Path('/workspace/heart-atlas/data/xenium').glob('WT_'+('P07' if stage=='P7' else stage)+'__*'))
 with tifffile.TiffFile(raw/'morphology_focus_0000.ome.tif',is_ome=False) as f:dapi=f.pages[0].pages[2].asarray(maxworkers=2)
 membrane=np.load(folder/f'{stage}_channel1_overview.npz')['image'];stroma=np.load(folder/f'{stage}_channel3_overview.npz')['image']
 assert dapi.shape==membrane.shape==stroma.shape
 def scale(x):
  lo,hi=np.quantile(x[x>0],[.01,.997]);return np.sqrt(np.clip((x.astype(float)-lo)/max(hi-lo,1),0,1))
 rgb=np.stack([scale(stroma),scale(membrane),scale(dapi)],axis=2)
 np.savez_compressed(folder/f'{stage}_morphology_only.npz',DAPI=dapi,membrane=membrane,stroma=stroma,pixel_um=np.array(1.7))
 fig,axs=plt.subplots(1,2,figsize=(14,8),layout='constrained')
 extent=[0,rgb.shape[1]*1.7,rgb.shape[0]*1.7,0]
 for ax,img,title in zip(axs,[scale(dapi),rgb],['DAPI','DAPI (blue), boundary (green), αSMA/vimentin (red)']):
  ax.imshow(img,cmap='gray',extent=extent);ax.set(title=f'{stage}: {title}',xlabel='x (µm)',ylabel='y (µm)')
 fig.savefig(OUT/'figures'/f'{stage}_morphology_only_review.png',dpi=180);plt.close(fig)
 for ch in [1,3]:(folder/f'{stage}_channel{ch}_overview.npz').unlink()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=['P0','P7','P14','P21']);a=p.parse_args()
 os.environ.setdefault('MPLCONFIGDIR','/tmp/anatomy-mpl')
 for c in [1,3]:acquire(a.stage,c)
 display(a.stage)
