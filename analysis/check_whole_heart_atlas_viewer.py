"""Optional offline browser verification (Python Playwright + local Chromium)."""
from pathlib import Path
import json
import os

import numpy as np
import pandas as pd
from playwright.sync_api import sync_playwright


def main():
    out=Path(__file__).resolve().parents[1]/'results/whole_heart_spatial_atlas';checks=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(executable_path=os.environ.get('ATLAS_CHROMIUM_PATH','/usr/bin/chromium'),
            headless=True,args=['--no-sandbox','--disable-dev-shm-usage'])
        context=browser.new_context(viewport={'width':1440,'height':1000},accept_downloads=True,offline=True)
        page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.set_content((out/'interactive_whole_heart_atlas.html').read_text(),wait_until='load',timeout=60000)
        page.locator('#loading').wait_for(state='detached',timeout=60000)
        page.wait_for_function("document.querySelector('#status').textContent.includes('47,611 displayed')")
        checks.append('P0: all47611 original cells displayed')
        page.locator('#noncm').click()
        page.wait_for_function("document.querySelector('#status').textContent.includes('10,306 displayed')")
        checks.append('P0: non-CM filter displays10306 cells')
        page.locator('#stage').select_option('P21')
        page.wait_for_function("document.querySelector('#status').textContent.includes('47,304 displayed')")
        checks.append('P21 switch preserves non-CM filter and displays47304 cells')
        page.locator('#all').click()
        page.wait_for_function("document.querySelector('#status').textContent.includes('143,944 displayed')")
        page.locator('#dapi').check();page.wait_for_timeout(300)
        page.locator('#resolved').click()
        page.wait_for_function("document.querySelector('#status').textContent.includes('108,511 displayed')")
        checks.append('P21: resolved filter displays108511 cells; embedded DAPI toggles offline')
        page.screenshot(path=str(out/'interactive_atlas_preview.png'),full_page=True)
        box=page.locator('#canvas').bounding_box()
        page.mouse.move(box['x']+box['width']/2,box['y']+box['height']/2)
        page.mouse.wheel(0,-800);page.wait_for_timeout(250)
        assert '108,511 displayed' not in page.locator('#status').inner_text()
        page.locator('#reset').click()
        page.wait_for_function("document.querySelector('#status').textContent.includes('108,511 displayed')")
        checks.append('Zoom changes visible count; reset restores complete section')
        with page.expect_download() as download_info:page.locator('#csv').click()
        download_info.value.save_as('/tmp/P21_atlas_browser_export.csv')
        downloaded=pd.read_csv('/tmp/P21_atlas_browser_export.csv',float_precision='round_trip')
        original=pd.read_csv(out/'P21/all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        assert len(downloaded)==143944
        assert downloaded.cell_id.tolist()==original.index.tolist()
        assert downloaded.atlas_identity.tolist()==original.atlas_identity.tolist()
        np.testing.assert_array_equal(downloaded[['x_um','y_um']],original[['x_centroid','y_centroid']])
        checks.append('CSV download retains all143944 P21 IDs, exact coordinates and identities')
        assert not errors,errors
        result={'status':'PASS','browser':browser.version,'offline':True,'checks':checks,'page_errors':errors}
        browser.close()
    (out/'viewer_browser_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
