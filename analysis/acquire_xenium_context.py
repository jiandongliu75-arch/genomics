"""Download morphology and segmentation context by ZIP member HTTP ranges."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import json
import zlib

from remotezip import RemoteZip

BASE = Path('/workspace/heart-atlas')


def acquire(task):
    archive, record = task
    folder = BASE / 'data/xenium' / archive['key'].removesuffix('.zip')
    wanted = {'morphology_focus_0000.ome.tif', 'cell_boundaries.csv.gz'}
    results = []
    with RemoteZip(archive['links']['self']) as remote:
        for info in remote.infolist():
            if info.filename.startswith('__MACOSX/') or Path(info.filename).name not in wanted:
                continue
            path = folder / Path(info.filename).name
            if path.exists() and path.stat().st_size == info.file_size:
                status = 'existing'
            else:
                print(f'DOWNLOAD {folder.name} {path.name} {info.file_size}', flush=True)
                data = remote.read(info.filename)  # ZipFile checks extracted member CRC.
                assert len(data) == info.file_size and zlib.crc32(data) == info.CRC
                temp = path.with_name(path.name + '.partial')
                temp.write_bytes(data)
                temp.replace(path)
                status = 'downloaded_CRC_verified'
                del data
            results.append({'record': record, 'archive': archive['key'], 'member': info.filename, 'local': str(path), 'bytes': info.file_size, 'zip_crc32': info.CRC, 'status': status})
            print(f'CONTEXT_READY {folder.name} {path.name}', flush=True)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    tasks = []
    for record in ['14961929', '14961942']:
        data = json.loads((BASE / f'metadata/zenodo_{record}.json').read_text())
        tasks.extend((archive, record) for archive in data['files'])
    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for results in pool.map(acquire, tasks):
            records.extend(results)
    (BASE / 'metadata/xenium_context_downloads.json').write_text(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
