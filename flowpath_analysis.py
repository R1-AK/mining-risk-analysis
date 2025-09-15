# -*- coding: utf-8 -*-
"""
Refactored Flowpath Analysis Module (from Notebook).
"""
import os
import math
import time
import zipfile
from pathlib import Path
import sys
from datetime import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import warp
from pysheds.grid import Grid
from shapely.geometry import Point
from tqdm import tqdm
import ee
import requests

def run_flowpath_analysis(mining_shapefile_path, gee_credentials, parameters, risk_indicators):
    """
    Main function to execute the GEE-based flowpath analysis workflow.
    """
    print("--- FLOWPATH ANALYSIS BACKEND (GEE METHOD) ---")
    
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Use a fixed output directory for final results
    OUTPUT_DIR = Path('./outputs').resolve()
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Use a temporary directory for intermediate files
    TEMP_DIR = Path(f'./temp_flowpath_{run_timestamp}').resolve()

    SERVICE_ACCOUNT = gee_credentials["SERVICE_ACCOUNT"]
    KEY_FILE = gee_credentials["SERVICE_ACCOUNT_KEY_FILE"]
    PROJECT_ID = gee_credentials["PROJECT_ID"]
    SCALE = parameters["SCALE"]
    RADIUS_M = parameters["RADIUS_M"]
    BANDS = [(0, 1_000, '0_1km'), (1_000, 5_000, '1_5km'), (5_000, 10_000, '5_10km'), (10_000, 20_000, '10_20km')]
    DEM_EE = 'WWF/HydroSHEDS/03CONDEM'
    WORLDPOP_2020 = 'WorldPop/GP/100m/pop'
    LC_2019 = 'COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019'
    LC_BAND = 'discrete_classification'
    NO_DATA = -32768

    def make_square(ee_geom, radius_m):
        proj = 'EPSG:3395'
        pt_m = ee_geom.transform(proj, 1)
        coords = pt_m.coordinates()
        x, y, r = ee.Number(coords.get(0)), ee.Number(coords.get(1)), ee.Number(radius_m)
        rect_m = ee.Geometry.Rectangle([x.subtract(r), y.subtract(r), x.add(r), y.add(r)], proj, False)
        return rect_m.transform('EPSG:4326', 1)

    def robust_ee_download(img, region, out_path, scale):
        out_path = Path(out_path)
        if out_path.exists(): return
        for attempt in range(6):
            try:
                url = img.getDownloadURL({'region': region, 'scale': scale, 'format': 'GEO_TIFF'})
                r = requests.get(url, stream=True, timeout=300)
                r.raise_for_status()
                tmp_path = out_path.with_suffix('.part')
                with open(tmp_path, 'wb') as f:
                    for chunk in r.iter_content(1 << 20): f.write(chunk)
                if zipfile.is_zipfile(tmp_path):
                    with zipfile.ZipFile(tmp_path, 'r') as zf:
                        tif_name = next((n for n in zf.namelist() if n.lower().endswith('.tif')), None)
                        if not tif_name: raise ValueError("Zip file from GEE did not contain a .tif file.")
                        with zf.open(tif_name) as zt, open(out_path, 'wb') as out: out.write(zt.read())
                    tmp_path.unlink()
                else:
                    tmp_path.rename(out_path)
                return
            except Exception as e:
                if attempt == 5: raise e
                time.sleep(2 ** attempt)

    def initialize_environment_local():
        print("Initializing environment...")
        for d in ['dem_tiles', 'source_masks', 'flowpath_tiles', 'pop_tiles', 'lc_tiles']:
            (TEMP_DIR / d).mkdir(parents=True, exist_ok=True)
        try:
            credentials = ee.ServiceAccountCredentials(SERVICE_ACCOUNT, KEY_FILE)
            ee.Initialize(credentials, project=PROJECT_ID)
            print('Earth Engine authentication successful.')
            return True
        except Exception as e:
            raise Exception(f"Fatal: Earth Engine initialization failed. Check credentials. Error: {e}")

    def create_lease_centroids_local(path):
        print("\nCreating centroids for each lease...")
        leases = gpd.read_file(path).to_crs('EPSG:4326')
        leases = leases.reset_index(drop=True).reset_index().rename(columns={'index': 'lease_id'})
        def get_accurate_centroid(geom):
            lon, lat = geom.representative_point().x, geom.representative_point().y
            utm_epsg = 32600 + int((lon + 180) / 6) + 1
            if lat < 0: utm_epsg += 100
            return gpd.GeoSeries([geom], crs='EPSG:4326').to_crs(utm_epsg).centroid.to_crs('EPSG:4326').iloc[0]
        centroids = leases.copy()
        centroids['geometry'] = centroids['geometry'].apply(get_accurate_centroid)
        print(f"   Processed {len(centroids)} centroids.")
        return leases, centroids

    def download_data_local(centroids_gdf, indicators):
        print("\nDownloading required data tiles from GEE...")
        centroids_ee = ee.FeatureCollection(centroids_gdf.__geo_interface__)
        feat_list = centroids_ee.toList(centroids_ee.size()).getInfo()
        dem_img = ee.Image(DEM_EE).select(0).updateMask(ee.Image(DEM_EE).select(0).gt(0))
        pop_img = ee.ImageCollection(WORLDPOP_2020).filter(ee.Filter.eq('year', 2020)).mosaic()
        lc_img = ee.Image(LC_2019).select(LC_BAND)
        for f in tqdm(feat_list, desc="Downloading GEE data"):
            lease_id = f['properties']['lease_id']
            ee_geom = ee.Geometry.Point(f['geometry']['coordinates'])
            region = make_square(ee_geom, RADIUS_M)
            robust_ee_download(dem_img, region, TEMP_DIR / 'dem_tiles' / f'dem_{lease_id}.tif', scale=SCALE)
            source_mask = ee.Image().toByte().paint(ee.Feature(f), 1).selfMask()
            robust_ee_download(source_mask, region, TEMP_DIR / 'source_masks' / f'source_{lease_id}.tif', scale=SCALE)
            if 'population' in indicators:
                robust_ee_download(pop_img, region, TEMP_DIR / 'pop_tiles' / f'pop_{lease_id}.tif', scale=100)
            if 'forest' in indicators or 'agriculture' in indicators:
                robust_ee_download(lc_img, region, TEMP_DIR / 'lc_tiles' / f'lc_{lease_id}.tif', scale=100)

    def process_flowpaths_local(centroids_gdf):
        print("\nCalculating individual flowpaths...")
        for _, lease in tqdm(list(centroids_gdf.iterrows()), desc="Processing flowpaths"):
            lease_id = lease.lease_id
            dem_path = TEMP_DIR / 'dem_tiles' / f'dem_{lease_id}.tif'
            source_path = TEMP_DIR / 'source_masks' / f'source_{lease_id}.tif'
            out_path = TEMP_DIR / 'flowpath_tiles' / f'flowpath_{lease_id}.tif'
            if not dem_path.exists() or not source_path.exists() or out_path.exists(): continue
            try:
                with rasterio.open(source_path) as src_src:
                    src_arr = src_src.read(1)
                idx = np.argwhere(src_arr == 1)
                if len(idx) != 1: continue
                row, col = idx[0]
                grid = Grid.from_raster(str(dem_path))
                dem = grid.read_raster(str(dem_path))
                dem_f = grid.fill_pits(dem); dem_f = grid.fill_depressions(dem_f); dem_f = grid.resolve_flats(dem_f)
                fdir = grid.flowdir(dem_f); fdir_arr = np.asarray(fdir)
                stepmap = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
                path = np.full(fdir_arr.shape, NO_DATA, dtype=np.float32)
                r, c = row, col
                dist_prev = 0.0
                while fdir_arr[r, c] in stepmap:
                    path[r, c] = dist_prev
                    dr, dc = stepmap[fdir_arr[r, c]]
                    nr, nc = r + dr, c + dc
                    if not (0 <= nr < path.shape[0] and 0 <= nc < path.shape[1]): break
                    dist_prev += float(math.hypot(dr, dc) * SCALE)
                    r, c = nr, nc
                with rasterio.open(dem_path) as dem_src:
                    profile = dem_src.profile.copy()
                profile.update(dtype='float32', nodata=NO_DATA, compress='lzw')
                with rasterio.open(out_path, 'w', **profile) as dst:
                    dst.write(path, 1)
            except Exception as e:
                print(f"   Warning: Could not delineate flowpath for lease {lease_id}. Error: {e}")

    def run_impact_analysis_local(leases_gdf, indicators):
        print("\nRunning impact analysis...")
        results_data = []
        for _, lease in tqdm(list(leases_gdf.iterrows()), desc="Analyzing lease impacts"):
            lease_id = lease.lease_id
            flow_path = TEMP_DIR / 'flowpath_tiles' / f'flowpath_{lease_id}.tif'
            if not flow_path.exists(): continue
            with rasterio.open(flow_path) as flow_src:
                dist_img = flow_src.read(1)
                dist_img[dist_img == NO_DATA] = np.nan
                dst_transform, dst_crs, dst_shape = flow_src.transform, flow_src.crs, flow_src.shape
            lease_result = {'lease_id': lease_id}
            if 'population' in indicators:
                pop_path = TEMP_DIR / 'pop_tiles' / f'pop_{lease_id}.tif'
                if pop_path.exists():
                    with rasterio.open(pop_path) as pop_src:
                        pop_arr_reproj = np.zeros(dst_shape, dtype=np.float32)
                        warp.reproject(pop_src.read(1), pop_arr_reproj, src_transform=pop_src.transform, src_crs=pop_src.crs, dst_transform=dst_transform, dst_crs=dst_crs)
                        pop_rescaled = pop_arr_reproj / ((100 / SCALE) ** 2)
                        for qi, qf, tag in BANDS:
                            m = (dist_img >= qi) & (dist_img < qf)
                            lease_result[f'people_{tag}'] = float(np.nansum(pop_rescaled[m]))
            if 'forest' in indicators or 'agriculture' in indicators:
                lc_path = TEMP_DIR / 'lc_tiles' / f'lc_{lease_id}.tif'
                if lc_path.exists():
                    with rasterio.open(lc_path) as lc_src:
                        lc_arr_reproj = np.zeros(dst_shape, dtype=np.uint8)
                        warp.reproject(lc_src.read(1), lc_arr_reproj, src_transform=lc_src.transform, src_crs=lc_src.crs, dst_transform=dst_transform, dst_crs=dst_crs, resampling=warp.Resampling.nearest)
                        area_km2 = (SCALE ** 2) / 1e6
                        if 'forest' in indicators:
                            forest = ((lc_arr_reproj >= 111) & (lc_arr_reproj <= 126))
                            for qi, qf, tag in BANDS:
                                lease_result[f'forest_km2_{tag}'] = float(np.sum(forest[(dist_img >= qi) & (dist_img < qf)]) * area_km2)
                        if 'agriculture' in indicators:
                            agri = (lc_arr_reproj == 40)
                            for qi, qf, tag in BANDS:
                                lease_result[f'agri_km2_{tag}'] = float(np.sum(agri[(dist_img >= qi) & (dist_img < qf)]) * area_km2)
            results_data.append(lease_result)
        return pd.DataFrame(results_data)

    if not initialize_environment_local():
        return "GEE initialization failed."
    leases_orig, centroids_gdf = create_lease_centroids_local(mining_shapefile_path)
    download_data_local(centroids_gdf, risk_indicators)
    process_flowpaths_local(centroids_gdf)
    results_df = run_impact_analysis_local(centroids_gdf, risk_indicators)
    
    # Cleanup temporary directory
    import shutil
    shutil.rmtree(TEMP_DIR)

    if results_df.empty:
        return "No results generated."
    print("\nFinalizing and saving results...")
    
    result_filename = f"flowpath_analysis_summary_{run_timestamp}.xlsx"
    results_path = OUTPUT_DIR / result_filename
    results_df.to_excel(results_path, index=False)

    return result_filename
