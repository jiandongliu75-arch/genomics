"""Read only the acquired DAPI TIFF page and its pyramid, preserving µm axes."""
from pathlib import Path
import json
import hashlib

import numpy as np
from PIL import Image
import tifffile

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/xenium_ventricular_mapping'

for stage in ['P0', 'P7', 'P14', 'P21']:
    pattern = 'WT_P07__*' if stage == 'P7' else f'WT_{stage}__*'
    source = next((BASE / 'data/xenium').glob(pattern))
    experiment = json.loads((source / 'experiment.xenium').read_text())
    path = source / 'morphology_focus_0000.ome.tif'
    folder = OUT / stage
    folder.mkdir(exist_ok=True)
    with tifffile.TiffFile(path) as handle:
        page = handle.pages[0]
        full = page.shape
        pages = [page, *list(page.pages)]
        level = next(i for i, candidate in enumerate(pages) if max(candidate.shape) <= 3500)
        array = pages[level].asarray()
        # Do not load the OME multi-file series; only DAPI has been acquired.
        assert array.ndim == 2 and 'Name="DAPI"' in handle.ome_metadata
    positive = array[array > 0]
    low, high = np.quantile(positive, [.02, .998])
    display = np.sqrt(np.clip((array.astype(float) - low) / max(high - low, 1), 0, 1))
    image = np.round(display * 255).astype(np.uint8)
    Image.fromarray(image).save(folder / 'DAPI_overview.png')
    px = float(experiment['pixel_size'])
    info = {'stage': stage, 'source_section': source.name, 'source': str(path), 'channel': 'DAPI only', 'full_resolution_shape_yx': list(full), 'pyramid_level': level, 'overview_shape_yx': list(image.shape), 'pixel_size_um_at_full_resolution': px, 'extent_um': [0, full[1] * px, full[0] * px, 0], 'display_transform': '2nd–99.8th percentile clipping of nonzero pyramid intensity, square-root display contrast, uint8. Analysis counts and coordinates unchanged.', 'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    (folder / 'morphology_display.json').write_text(json.dumps(info, indent=2))
    print('DAPI_READY', stage, full, image.shape, flush=True)
