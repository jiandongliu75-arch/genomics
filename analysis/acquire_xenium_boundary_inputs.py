"""Acquire exact ZIP members with resumable bounded HTTP range decompression."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import struct
import time
import zlib

import requests
from remotezip import RemoteZip

BASE=Path('/workspace/heart-atlas')
CACHE=BASE/'results/xenium_boundary_inputs'
STAGES=['P0','P7','P14','P21']


def discover(stage):
    record='14961929' if stage in ['P0','P7'] else '14961942'
    data=json.loads((BASE/f'metadata/zenodo_{record}.json').read_text())
    alias='P07' if stage=='P7' else stage
    archive=next(f for f in data['files'] if f['key'].startswith('WT_'+alias+'__'))
    return archive,record


def ranged(session,url,start,end):
    for attempt in range(6):
        try:
            response=session.get(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},timeout=(30,150))
            response.raise_for_status()
            assert response.status_code==206,(response.status_code,response.headers.get('Content-Range'))
            assert response.headers['Content-Range'].startswith(f'bytes {start}-'),response.headers
            assert len(response.content)==end-start+1,(len(response.content),start,end)
            return response.content
        except Exception:
            if attempt==5:raise
            time.sleep(min(8,1+attempt))


def extract(stage,name):
    archive,record=discover(stage);folder=CACHE/stage;folder.mkdir(parents=True,exist_ok=True)
    url=archive['links']['self']
    with RemoteZip(url) as remote:
        candidates=[i for i in remote.infolist() if not i.filename.startswith('__MACOSX/') and Path(i.filename).name==name]
        assert len(candidates)==1,(stage,name,len(candidates))
        info=candidates[0]
    path=folder/name;meta=folder/(name+'.acquisition.json')
    if path.exists() and meta.exists() and path.stat().st_size==info.file_size:
        print('CACHED',stage,name,path.stat().st_size,flush=True);return path
    session=requests.Session();header=ranged(session,url,info.header_offset,info.header_offset+29)
    assert header[:4]==b'PK\x03\x04'
    fn,extra=struct.unpack_from('<HH',header,26);offset=info.header_offset+30+fn+extra
    assert info.compress_type in [0,8]
    decoder=zlib.decompressobj(-15) if info.compress_type==8 else None
    temp=path.with_suffix(path.suffix+'.partial');crc=0;written=0;block=16*1024*1024
    print('ACQUIRE',stage,name,'compressed',info.compress_size,'bytes',flush=True)
    positions=list(range(0,info.compress_size,block))
    workers=6
    with temp.open('wb') as target,ThreadPoolExecutor(max_workers=workers) as pool:
        def fetch(position):
            size=min(block,info.compress_size-position)
            return ranged(requests.Session(),url,offset+position,offset+position+size-1)
        pending={p:pool.submit(fetch,p) for p in positions[:workers]}
        for index,position in enumerate(positions):
            size=min(block,info.compress_size-position)
            compressed=pending.pop(position).result()
            if index+workers<len(positions):
                upcoming=positions[index+workers];pending[upcoming]=pool.submit(fetch,upcoming)
            plain=decoder.decompress(compressed) if decoder else compressed
            target.write(plain);written+=len(plain);crc=zlib.crc32(plain,crc)
            if position==0 or (position//block)%8==0 or position+size==info.compress_size:
                print('PROGRESS',stage,name,position+size,info.compress_size,flush=True)
        if decoder:
            tail=decoder.flush();target.write(tail);written+=len(tail);crc=zlib.crc32(tail,crc)
            assert decoder.eof and not decoder.unused_data
    assert written==info.file_size and crc==info.CRC,(written,info.file_size,crc,info.CRC)
    temp.replace(path)
    meta.write_text(json.dumps({'stage':stage,'record':record,'archive':archive['key'],'url':url,'member':info.filename,
        'compressed_bytes':info.compress_size,'bytes':written,'zip_crc32':crc,'bounded_range_bytes':block,
        'CRC_verified':True},indent=2)+'\n')
    print('ACQUIRED_CRC_PASS',stage,name,written,flush=True);return path


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=STAGES,required=True);p.add_argument('--members',nargs='+',required=True);args=p.parse_args()
    for member in args.members:extract(args.stage,member)


if __name__=='__main__':main()
