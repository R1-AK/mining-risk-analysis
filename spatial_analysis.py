# -*- coding: utf-8 -*-
"""
Comprehensive Mining Impact Analysis.
"""

import os
import sys
import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union
from tqdm import tqdm
import rasterio
from rasterio.mask import mask
import numpy as np
from datetime import datetime
from pathlib import Path

def run_spatial_analysis(mining_filepath, buffer_kms, risk_data_config):
    """
    Main function to execute the spatial (buffer) analysis workflow.
    """
    print("--- SPATIAL ANALYSIS BACKEND ---")

    PROJECTED_CRS = 'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],UNIT["Meter",1.0]]'
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Use a fixed output directory name
    output_dir = Path('./outputs').resolve()
    output_dir.mkdir(exist_ok=True)

    all_results = []

    print(f"Loading mining data from: {mining_filepath}")
    mining_gdf = gpd.read_file(mining_filepath).to_crs(PROJECTED_CRS)
    dissolved_mining_geom = unary_union(mining_gdf.geometry)

    print("Creating non-cumulative ring buffers...")
    buffer_distances_m = [km * 1000 for km in buffer_kms]
    all_cumulative_buffers = {dist: dissolved_mining_geom.buffer(dist) for dist in buffer_distances_m}
    sorted_dists = sorted(buffer_distances_m)

    rings = {f"0-{sorted_dists[0] / 1000} km": all_cumulative_buffers[sorted_dists[0]]}
    for i in range(1, len(sorted_dists)):
        inner_dist_m, outer_dist_m = sorted_dists[i - 1], sorted_dists[i]
        rings[f"{inner_dist_m / 1000}-{outer_dist_m / 1000} km"] = all_cumulative_buffers[outer_dist_m].difference(
            all_cumulative_buffers[inner_dist_m])

    for config in risk_data_config:
        print(f"\nProcessing dataset: {config['name']} ({config['type']})")
        is_biodiversity = config.get('subcategory') == 'Terrestrial Biodiversity'
        if is_biodiversity:
            results = _analyze_species_ring(rings, config, PROJECTED_CRS)
        elif config['type'] == 'VECTOR':
            results = _analyze_vector_ring(rings, config, PROJECTED_CRS)
        elif config['type'] == 'RASTER':
            results = _analyze_raster_ring(rings, config, PROJECTED_CRS)
        all_results.extend(results)

    if not all_results:
        return "Analysis complete, but no results were generated."

    results_df = pd.DataFrame(all_results)
    
    # Create a unique filename instead of a unique directory
    result_filename = f"spatial_analysis_summary_{run_timestamp}.xlsx"
    results_path = output_dir / result_filename
    results_df.to_excel(results_path, index=False)

    print(f"\nAnalysis complete. Results saved to: {results_path}")
    
    # Return ONLY the unique filename
    return result_filename


# --- (All helper functions like _get_subcategory, _filter_species_data, etc. remain the same) ---

def _get_subcategory(config):
    if config.get('subcategory') == 'Others':
        return config.get('other_subcategory', 'Others (unspecified)')
    return config.get('subcategory', 'General')

def _filter_species_data(config, target_crs):
    try:
        species_data = gpd.read_file(config['path'])
        taxon_name = config['name'].lower()
        if 'reptile' in taxon_name:
            species_id_col = 'Binomial'
            if 'Group_' in species_data.columns:
                species_data = species_data[species_data['Group_'].isin(['lizard', 'turtle', 'snake'])]
        else:
            species_id_col = 'sci_name'
            if 'category' in species_data.columns:
                species_data = species_data[species_data['category'].isin(['EN', 'CR', 'VU'])]
            if 'presence' in species_data.columns:
                species_data = species_data[species_data['presence'].isin([1, 2, 3])]
            if 'origin' in species_data.columns:
                species_data = species_data[species_data['origin'].isin([1, 2, 3])]
        if species_id_col not in species_data.columns: return gpd.GeoDataFrame()
        species_data.rename(columns={species_id_col: 'species_name'}, inplace=True)
        return species_data[['species_name', 'geometry']].to_crs(target_crs)
    except Exception as e:
        print(f"   Warning: Could not filter {config['name']}. Error: {e}")
        return gpd.GeoDataFrame()

def _analyze_species_ring(rings, config, target_crs):
    all_results = []
    species_gdf = _filter_species_data(config, target_crs)
    total_unique_species = species_gdf['species_name'].nunique()
    if total_unique_species == 0:
        print(f"   No threatened species found for {config['name']} after filtering. Skipping.")
        return all_results
    print(f"   Found {total_unique_species} unique threatened species for analysis.")
    subcategory_name = _get_subcategory(config)
    for ring_name, ring_geom in tqdm(rings.items(), desc=f"   Analyzing rings for {config['name']}"):
        if ring_geom.is_empty: continue
        buffer_gdf = gpd.GeoDataFrame(geometry=[ring_geom], crs=target_crs)
        possible_matches = species_gdf.iloc[list(species_gdf.sindex.intersection(buffer_gdf.total_bounds))]
        num_impacted = 0
        if not possible_matches.empty:
            impacted_gdf = gpd.overlay(possible_matches, buffer_gdf, how='intersection')
            num_impacted = impacted_gdf['species_name'].nunique()
        percent_impacted = (num_impacted / total_unique_species) * 100 if total_unique_species > 0 else 0
        all_results.append({'Main Category': config.get('main_category', 'N/A'), 'Subcategory': subcategory_name, 'Dataset': config['name'], 'Impact Zone': ring_name, 'Class': 'Count of Species', 'Value': num_impacted, 'Unit': 'Number of Species'})
        all_results.append({'Main Category': config.get('main_category', 'N/A'), 'Subcategory': subcategory_name, 'Dataset': config['name'], 'Impact Zone': ring_name, 'Class': 'Percentage of Species', 'Value': percent_impacted, 'Unit': '%'})
    return all_results

def _analyze_vector_ring(rings, config, target_crs):
    results = []
    vector_gdf = gpd.read_file(config['path']).to_crs(target_crs)
    subcategory_name = _get_subcategory(config)
    for ring_name, ring_geom in tqdm(rings.items(), desc=f"   Analyzing rings for {config['name']}"):
        if ring_geom.is_empty: continue
        total_area_m2 = ring_geom.area
        clipped = gpd.clip(vector_gdf, ring_geom)
        if not clipped.empty:
            if config['processing'] == 'CATEGORICAL' and 'column' in config:
                col = config['column']
                clipped['area'] = clipped.geometry.area
                class_areas = clipped.groupby(col)['area'].sum()
                for class_value, area_m2 in class_areas.items():
                    results.append({'Main Category': config.get('main_category', 'N/A'), 'Subcategory': subcategory_name, 'Dataset': config['name'], 'Impact Zone': ring_name, 'Class': class_value, 'Value': (area_m2 / total_area_m2) * 100 if total_area_m2 > 0 else 0, 'Unit': 'Area Overlap (%)'})
            else:
                overlap_m2 = clipped.area.sum()
                results.append({'Main Category': config.get('main_category', 'N/A'), 'Subcategory': subcategory_name, 'Dataset': config['name'], 'Impact Zone': ring_name, 'Class': 'Total', 'Value': (overlap_m2 / total_area_m2) * 100 if total_area_m2 > 0 else 0, 'Unit': 'Area Overlap (%)'})
    return results

def _analyze_raster_ring(rings, config, target_crs):
    results = []
    subcategory_name = _get_subcategory(config)
    with rasterio.open(config['path']) as src:
        for ring_name, ring_geom in tqdm(rings.items(), desc=f"   Analyzing rings for {config['name']}"):
            if ring_geom.is_empty: continue
            total_area_m2 = ring_geom.area
            try:
                buffer_reproj = gpd.GeoDataFrame(geometry=[ring_geom], crs=target_crs).to_crs(src.crs)
                masked, _ = mask(src, buffer_reproj.geometry, crop=True, nodata=src.nodata)
                masked = masked[0]
                if config['processing'] == 'CATEGORICAL' and 'classes' in config and config['classes']:
                    class_defs = [item.split(',') for item in config['classes'].split(';')]
                    for min_val, max_val, name in class_defs:
                        class_mask = (masked >= int(min_val)) & (masked <= int(max_val))
                        pixel_count = np.sum(class_mask)
                        pixel_area = abs(src.transform[0] * src.transform[4])
                        area_m2 = pixel_count * pixel_area
                        results.append({'Main Category': config.get('main_category', 'N/A'), 'Subcategory': subcategory_name, 'Dataset': config['name'], 'Impact Zone': ring_name, 'Class': name, 'Value': (area_m2 / total_area_m2) * 100 if total_area_m2 > 0 else 0, 'Unit': 'Area Coverage (%)'})
                else:
                    valid_pixels = masked[(masked != src.nodata) & (~np.isnan(masked))]
                    if valid_pixels.size > 0:
                        results.append({'Main Category': config.get('main_category', 'N/A'), 'Subcategory': subcategory_name, 'Dataset': config['name'], 'Impact Zone': ring_name, 'Class': 'Average Value', 'Value': valid_pixels.mean(), 'Unit': 'Mean Pixel Value'})
            except (ValueError, IndexError):
                pass
    return results
