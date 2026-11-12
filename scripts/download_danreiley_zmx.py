#!/usr/bin/env python3
"""
Downloads all Zemax .zmx files from Dan Reiley's lens design exchange site:
https://sites.google.com/site/danreiley/a-file-exchange-site-for-lens-designs
"""
import urllib.request
import re
import os
import sys
import time
import concurrent.futures

CATEGORIES = {
    'photographic-lenses-prime': '1P38j7dtEsnPvlGWBaxrLd4tA98hyiyds',
    'photographic-lenses-zoom': '1x8_9hb72-Gyvmlx-DE026UXR-zdPACB0',
    'projectors': '1rz0OD5erCPJG_8NLyLry5CSfbS4Cy32a',
    'microscope-objectives': '1Gi4ZnWE4G_m-pWJbeuRQC7O9uXqP6WAO',
    'telescopes': '1leEyQualluLTwUfJSNr2TcgR5S1i0xoR',
    'eyepieces': '1F4GBcEohFNWh9kLPcS7aDgur_DHOADFO',
    'endoscopes': '1eXA_cBzkEosQz7-Q1JuPc2gIcz0G3BNw',
    'scan-lenses': '13IjbkaoEkGwKN3MnEK-Yub9vNpucPHH4',
    'spectro': '18fU4UPQz2InGtwYzlBUIF91Q1SV2fIf7',
}

DEST_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'test_data', 'dan_reiley'))

def get_folder_entries(cat, fid):
    url = f'https://drive.google.com/embeddedfolderview?id={fid}#list'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    html = urllib.request.urlopen(req).read().decode('utf-8', errors='ignore')
    entries = re.findall(r'id=\"entry-([^\"]+)\".*?<div class=\"flip-entry-title\">([^<]+)</div>', html, re.DOTALL)
    zmx_entries = [(cat, fid_entry, name.strip()) for fid_entry, name in entries if name.strip().lower().endswith('.zmx')]
    return zmx_entries

def download_file(item):
    cat, fid, name = item
    cat_dir = os.path.join(DEST_DIR, cat)
    os.makedirs(cat_dir, exist_ok=True)
    dest_path = os.path.join(cat_dir, name)
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 100:
        return True, name, 'already exists'
    
    url = f'https://drive.google.com/uc?export=download&id={fid}'
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            data = urllib.request.urlopen(req, timeout=15).read()
            if len(data) > 0:
                with open(dest_path, 'wb') as f:
                    f.write(data)
                return True, name, f'{len(data)} bytes'
        except Exception as e:
            time.sleep(1 + attempt)
    return False, name, 'failed after 3 attempts'

def main():
    os.makedirs(DEST_DIR, exist_ok=True)
    print(f'Fetching file lists for {len(CATEGORIES)} categories...')
    all_tasks = []
    for cat, fid in CATEGORIES.items():
        try:
            tasks = get_folder_entries(cat, fid)
            print(f'  {cat}: {len(tasks)} files')
            all_tasks.extend(tasks)
        except Exception as e:
            print(f'  Error indexing {cat}: {e}', file=sys.stderr)
            
    print(f'Total .zmx files found: {len(all_tasks)}')
    print(f'Starting concurrent download to {DEST_DIR}...')
    
    start_time = time.time()
    success_count = 0
    fail_count = 0
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(download_file, task): task for task in all_tasks}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            task = futures[future]
            success, name, info = future.result()
            if success:
                success_count += 1
            else:
                fail_count += 1
                print(f'Failed: {task[0]}/{name}: {info}', file=sys.stderr)
            if i % 100 == 0 or i == len(all_tasks):
                print(f'Progress: {i}/{len(all_tasks)} ({success_count} ok, {fail_count} failed)')
                
    elapsed = time.time() - start_time
    print(f'Done in {elapsed:.1f}s. Successfully downloaded: {success_count}, Failed: {fail_count}')

if __name__ == '__main__':
    main()
