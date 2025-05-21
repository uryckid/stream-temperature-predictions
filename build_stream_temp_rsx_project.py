#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 14:30:58 2025

@author: dawn urycki
"""

'''
This script does stuff with stuff
'''

import os
import datetime
import time
import pandas as pd
import shutil
import sqlite3
import zipfile
import glob
import geopandas as gpd
from shapely import Point, LineString, Polygon
from shapely import force_2d
from rsxml.project_xml import (
   Project,
   MetaData,
   Meta,
   ProjectBounds,
   Coords,
   BoundingBox,
   Realization,
   Geopackage,
   GeopackageLayer,
   GeoPackageDatasetTypes,
   Dataset
)


import sys


def create_proj_dir(huc6, new_dir, ref_path):
    path = os.path.join(ref_path, f'{huc6}/{new_dir}')
    
    if not os.path.exists(path):
        os.makedirs(path)
        print(f'Project directory created: {os.path.relpath(path, ref_path)}')
    else: print(f'Project directory exists: {os.path.relpath(path, ref_path)}')
    return(path)

def add_file(filepath, dest, ref_path = None):
    filename =  filepath.split("/")[-1]
    shutil.copy(filepath, dest)
    return(print(f'\n{filename} successfully added to {os.path.relpath(dest, ref_path)}.'))
    
def add_database(df, db_name, db_path, zip_file = True):
    connection= sqlite3.connect(db_path)
    df.to_sql(db_name, connection, if_exists = 'replace')
    connection.close()
    print('Database created successfully')
    if zip_file:
        zip_path = db_path + '.zip'
        try:
            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(db_path, os.path.basename(db_path))
            os.remove(db_path)
            print(f"Successfully zipped '{db_name}' to '{os.path.relpath(zip_path)}'")
        except:
            print(f"\nERROR zipping database {db_name}.")
    return

def create_database_file(file, date_col, proj_path, db_type, compression = None, add_cols = None, zip_file = True, overwrite = False):
    df = pd.read_csv(file, parse_dates = [date_col], compression = compression)
    df.columns = df.columns.str.split(r'\.').str[-1].str.strip()
    df.rename({'COMID': 'comid', 'Huc10': 'HUC10'}, axis = 1, inplace = True)
    out_path = os.path.join(proj_path, f'daily_{db_type}.db')
    if not (os.path.exists(out_path)) or (overwrite!=True):
        if db_type == 'stream_temperature':
            try:
                dropcols = ['antec_air_temp', 'std_mean_flow']
                df_trim = df.drop(dropcols, axis = 1, errors = 'ignore')
                df_out = df_trim[['comid', 'date', 'stream_temp']]
                print(f'\n\nStream temperature file processed: {len(df_out.comid.unique())} comids')
            except: print('\n\nStream temperature file FAILED.')
       
        else: 
            try:
                df = df.loc[:,~df.columns.duplicated()].copy()
                df_out = pd.merge(df, add_cols, how = 'left', on=['comid','date'])
                print(f'\n\nTemporal covariates file processed: {len(df_out.comid.unique())} comids.')
            except: print('\n\nTemporal covariates file FAILED.')
       
        out_table = df_out.fillna(-999).round(4)
        print('Writing to database...')
        add_database(out_table, db_type, db_path = out_path, zip_file = zip_file)
    return(df)

def check_comids(comids, data, layer_name):
    try: no_data = set(comids).difference(set(data.reset_index().comid.unique()))
    except: no_data = set(comids).difference(set(data.reset_index().featureid.unique()))
    if no_data==set(): 
        return(print(f'\n\n{layer_name} successfully loaded for {len(comids)} comids'))
    else: 
        return(print(f'\n\n{layer_name} missing data for {round((len(no_data) / len(comids)), 2)} of comids:', 
                     [str(n) for n in no_data]))

def layer_to_gpkg(comids, data, path, layer_name, epsg = 'EPSG:4326', check=True):
    if check: check_comids(comids, data, layer_name)
    try: 
        data.to_crs(epsg = epsg, inplace = True)
    except: pass

    if not os.path.exists(path):
        data.round(4).to_file(path, layer = layer_name, driver = 'GPKG')
        print('\n\nGeopackage successfully initiated.')
    else: gpd.GeoDataFrame(data.round(4)).to_file(path, layer = layer_name)
    
    return(print(f'\n{layer_name} added to {os.path.relpath(path)}.'))
    
def load_spatial_covariates(file):
    df = pd.read_csv(file, index_col = 'COMID')
    df.columns = df.columns.str.split(r'\.').str[-1].str.strip()
    df.index.rename('comid', inplace = True)
    return(df)

def find_containing_huc12(h12s, wings):
    h12s = h12s.dissolve(by='HUC12')
    
    wings_list = []
    h12s_list = []
    fail_list = []
    for i, wing in wings.iterrows():
        cont_h12 = h12s[h12s.geometry.contains(wing.geometry.representative_point())]
        if len(cont_h12) == 1:
            wings_list.append(i)
            h12s_list.append(cont_h12.index[0])
        else: 
            fail_list.append(i)
    if fail_list == []: print(f'\n\n{len(wings_list)} contributing areas successfully matched to HUC 12s.')
    else: print(f'\n\nERROR finding HUC 12s containing {len(fail_list)} contributing areas: {fail_list}')
    return(pd.DataFrame(zip(wings_list, h12s_list), columns = ['featureid', 'HUC12']))
        
def assign_season(doy):
    """Assigns a season based on the day of year (doy)."""
    if 1 <= doy <= 79:
        return "winter"
    elif 80 <= doy <= 171:
        return "spring"
    elif 172 <= doy <= 263:
        return "summer"
    elif 264 <= doy <= 356:
        return "fall"
    else:
        return "winter" 
    
def assign_water_year(date):
    """Assigns water year (Oct-Sept) based on date."""
    if date.month >= 10:
        return (date.year +1)
    else: return date.year
    
def getSeasonalAnomalies(indf, baseline_period: range, anomaly_periods: list, cols: list):
    print(f'\n\nCalculating seasonal anomalies: {cols}')
    try:
        cur_huc_anoms = []
        cur_comids = list(indf.comid.unique())
        for j in cur_comids:
            #Assign water year and season
            cur_df = indf[indf.comid == j][['date']+ cols]
            cur_ens_med = cur_df.groupby('date').median() # sets index to date
            cur_ens_med['doy'] = cur_ens_med.index.dayofyear
            cur_ens_med['season'] = cur_ens_med['doy'].apply(assign_season)
            for col in cols: cur_ens_med[f'n_{col}'] = cur_ens_med[col].notnull().astype(int)
            new_cols = cols
            
            if 'SWE' in cols:
                #Use water year for swe
                cur_swe = cur_ens_med.copy()
                cur_swe['WY'] = cur_swe.index.to_series().apply(assign_water_year)
                
                # Reclassify end of Sept as summer (pull it into summer of WY[y-1])
                cur_swe.loc[cur_swe.index.month == 9, 'season'] = 'summer'  
                swe_cum = cur_swe[['SWE', 'WY', 'season', 'n_SWE']].groupby(['WY', 'season']).sum()
                
                new_cols = [x.replace('SWE', 'SwS') for x in cols]
                
            cur_anomalies = {}
            
            #Determine baseline (median)
            cur_met = cur_ens_med.copy()
            baseline = cur_met[cur_met.index.year.isin(baseline_period)].groupby('season').median()
            for col in cols: 
                baseline[f'n_{col}'] = cur_met[cur_met.index.year.isin(baseline_period)].groupby('season').sum()[f'n_{col}'] # n

            if 'SWE' in cols:
                swe_df = swe_cum.loc[pd.IndexSlice[baseline_period], :, :]                
                baseline['SwS'] = swe_df.groupby(level=1).median()['SWE']/1000 # mm to m
                baseline['n_SwS'] = swe_df.groupby(level=1).sum()['n_SWE'] # n
            
            
            baseline_trim = pd.concat([baseline[[col, f'n_{col}']] for col in new_cols], axis = 1)
            
            cur_anomalies[f'{baseline_period[0]}s'] = baseline_trim
    
            #Calculate anomalies
            for per in anomaly_periods:                
                cur_anoms = cur_ens_med[cur_ens_med.index.year.isin(per)].groupby('season').median().subtract(baseline)
                for col in cols: 
                    cur_anoms[f'n_{col}'] = cur_ens_med[cur_ens_med.index.year.isin(per)].groupby('season').sum()[f'n_{col}'] # n
                        
                if 'SWE' in cols:
                    swe = swe_cum.loc[pd.IndexSlice[per], :, :]                
                    cur_anoms['SwS'] = ((swe.groupby(level=1).median()['SWE'])/1000).subtract(baseline['SwS'])
                    cur_anoms['n_SwS'] = swe.groupby(level=1).sum()['n_SWE'] # n
                
                cur_anoms_trim = pd.concat([cur_anoms[[col, f'n_{col}']] for col in new_cols], axis = 1)
                
                cur_anomalies[f'{per[0]}s'] = cur_anoms_trim
            
            cur_seas_anomalies = pd.concat(cur_anomalies, axis = 1)
            cur_huc_anoms.append(cur_seas_anomalies)
        print('Done.')
        res = (pd.concat(cur_huc_anoms, keys = cur_comids))
        res.index.names = ['comid', 'season']
        return res.reset_index()
    except: 
        return(print('Error calculating anomalies.'))
    
def collapse_levels(indf, level_name):
    outdf = indf[['spring', 'summer', 'fall', 'winter']]
    outdf.columns = [f'{level_name}_{c[:2]}' for c in outdf.columns]
    return(outdf)
    
def parse_epochs(indf: pd.DataFrame, epochs: list):
    epoch_dict = {}
    idx = indf[['comid', 'season']]
    idx.columns = idx.columns.droplevel(1)
    df = pd.DataFrame(idx)
    for e in epochs:
        vardf = pd.concat([df, indf.xs(e, axis = 1).copy()], axis = 1)
        n_cols = [c for c in vardf.columns if c[:2] == 'n_']
        pv_cols = [c for c in vardf.columns[2:] if c not in n_cols]
        varpv = vardf.pivot(index = 'comid', columns = 'season', values = pv_cols)
        varnpv = vardf.pivot(index = 'comid', columns = 'season', values = n_cols)
        frames = []
        n_frames = []
        for col in pv_cols:
            curdf = varpv.xs(col, axis = 1)
            frames.append(collapse_levels(curdf, col))
            ndf = varnpv.xs(f'n_{col}', axis = 1)
            n_frames.append(collapse_levels(ndf, f'n_{col}'))
        epoch_dict[e] = pd.concat([pd.concat(frames, axis = 1), pd.concat(n_frames, axis = 1)], axis = 1)
    return(epoch_dict)

def area_weighted_mean(metric_df, weights, agg_units):
    weighted_metrics = metric_df.apply(lambda x: x*metric_df.index.map(weights))
    agg_units_expanded = pd.Series(weighted_metrics.index.map(agg_units), index = weighted_metrics.index, name = agg_units.name)
    by_agg_unit = pd.concat([weighted_metrics, agg_units_expanded], axis = 1)
    weighted_mean = by_agg_unit.groupby([agg_units.name]).sum()
    return(weighted_mean)
    
def create_index(dbpath, layer, ixname):
    conn = sqlite3.connect(dbpath)
    with conn:
        curs = conn.cursor()
        try:
            curs.execute(f'''CREATE INDEX ix_{layer}_{ixname} 
                         ON {layer}({ixname});''')
            conn.commit()
            return(print(f'\nix created: {ixname}'))
        except: return(print('\nFailed to create index on {ixname}'))

def add_view(dbpath, layer_name, feature_table, f_index, a_index, shape = 'polygons'):
    conn = sqlite3.connect(dbpath)
    with conn:
        curs = conn.cursor()
        try: 
            curs.execute(f'''
                         CREATE VIEW vw_{layer_name}_{shape} AS 
                         SELECT f.fid, f.geom, a.* 
                         FROM {feature_table} AS f INNER JOIN {layer_name} a
                         ON f.{f_index} = a.{a_index};''')
            
            curs.execute(f'''
                         INSERT INTO gpkg_contents (table_name, data_type, identifier, min_x, min_y, max_x, max_y, srs_id) 
                         SELECT 'vw_{layer_name}_{shape}', data_type, 'vw_{layer_name}_{shape}', min_x, min_y, max_x, max_y, srs_id 
                         FROM gpkg_contents WHERE table_name = '{feature_table}';''')
            
            curs.execute(f'''
                         INSERT INTO gpkg_geometry_columns (table_name, column_name, geometry_type_name, srs_id, z, m)
                         SELECT 'vw_{layer_name}_{shape}', column_name, geometry_type_name, srs_id, z, m
                         FROM gpkg_geometry_columns
                         WHERE table_name = '{feature_table}';''')         

            conn.commit()
            print(f'\nSuccesfully added spatial view: vw_{layer_name}_{shape}')        
            
        except: pass
    return

def get_datasets(output_gpkg: str):
   """
   Returns a list of the datasets from the output GeoPackage.
   """

   conn = sqlite3.connect(output_gpkg)
   conn.enable_load_extension(True)
   curs = conn.cursor()
   
   # Get the names of all the tables in the database
   curs.execute("SELECT table_name FROM gpkg_contents")
   datasets = [GeopackageLayer(
       lyr_name=row[0],
       ds_type=GeoPackageDatasetTypes.VECTOR,
       name=row[0]
   ) for row in curs.fetchall()]
   return datasets
    
#%%
# USER INPUTS
proj_type = 'Stream_Temperature_Retrospective_'
timeframe = '1990-2021_'

# Required
curhuc = '1701010108' #huc10
huc8 = curhuc[:8]
proj_crs = 'EPSG:9674' # NAD83 / USFS R6 Albers projection, Oregon and Washington

# Optional
os.chdir(os.path.dirname(__file__))
curdir = os.getcwd()
datadir = os.path.join(curdir, 'DATA/')
geojson_dir = glob.glob(os.path.join(datadir, '*HUC17_project_bounds/' ))[0]
overwrite = False

# Define anomaly periods
periods = {'2000s': list(range(2000, 2010)),
           '2010s': list(range(2010, 2020)),
           '2020s': list(range(2020, 2030)), 
           '1990s': list(range(1990, 2000)), 
           '2050s': list(range(2050, 2060)), 
           '2080s': list(range(2080, 2090))}

anom_covs = ['air_temp_ws', 'NWM_flow_log', 'SWE']


#%% Define directories 

outdir = os.path.join(datadir, 'Outputs/retro/')
tempdir = os.path.join(datadir, 'preds_retro/')
nhddir = os.path.join(datadir, 'NHDPlusPN/NHDPlus17/')

#%% Load in geometry files

hydro_files = glob.glob(os.path.join(nhddir, r'NHDSnapshot/Hydrography/*.dbf'))

flowlines = gpd.read_file(hydro_files[0])
flowlines['geometry'] = flowlines['geometry'].apply(lambda x: force_2d(x))
flowlines.rename(columns = str.lower, inplace=True)

spatial_covs_file = os.path.join(datadir, 'spatial_data.csv')
spatial_covs = load_spatial_covariates(spatial_covs_file)

contrib_area_file = os.path.join(nhddir, 'NHDPlusCatchment/Catchment.dbf')
contrib_areas = gpd.read_file(contrib_area_file)
contrib_areas.rename({'FEATUREID': 'featureid', 'AreaSqKM': 'area_sqkm'}, axis = 1, inplace = True)

huc12s_file = os.path.join(nhddir, 'WBDSnapshot/WBD/WBD_Subwatershed.dbf')
huc12s = gpd.read_file(huc12s_file, columns = ['OBJECTID', 'HUC_12', 'HUC_8','HU_12_NAME','HUC_10','HU_10_NAME','geometry'])
huc12s.rename({'HUC_12': 'HUC12', 'HU_12_NAME': 'name'}, axis = 1, inplace = True)
h10name = huc12s.loc[huc12s.HUC_10 == curhuc]['HU_10_NAME'].mode()[0]
h10name = h10name.replace(" ","_")
#%% Create project directory and add daily temperature and covariate databases

proj_name = proj_type + timeframe + h10name
proj_path = create_proj_dir(huc8[:6], proj_name, outdir)
log_path = os.path.join(proj_path, 'processing.log')
terminal_stdout = sys.stdout
with open(log_path, 'w') as f:
    sys.stdout = f
    #sys.stdout = open(log_path,'w')
    start = datetime.datetime.now()
    print(start.strftime("%Y-%m-%d %H:%M:%S"))
    print(f'\nBuilding {proj_type} for {h10name}')
    
    # Add README
    readme_path = os.path.join('DATA', 'preds_retro', 'readme.txt')
    add_file(readme_path, proj_path, outdir)
    
    # Add covariate metadata
    add_file(os.path.join('DATA', 'covariate_metadata.csv'), proj_path, outdir)
    
    # Temperature
    temp_file = os.path.join(tempdir, 'predictions_temperature', f'st_pred_{curhuc}.csv')
    temp_df = create_database_file(temp_file, date_col = 'tim.date', proj_path = proj_path, 
                                   db_type = 'stream_temperature', add_cols = None, compression = None, 
                                   overwrite = overwrite)
    
    # Covariates
    cov_file = os.path.join(tempdir, 'predictions_covariates', 'cov_csvs', f'{curhuc}_covs.zip')
    cov_df = create_database_file(cov_file, date_col = 'date', proj_path = proj_path,
                                  db_type = 'covariates', add_cols = temp_df[['comid', 'date', 'antec_air_temp', 'std_mean_flow']], 
                                  compression = 'zip', overwrite = overwrite)
    
    #%% Build geopackage
    
    gpkg_filename = 'seasonal_anomalies_spatial_covariates.gpkg'
    gpkg_path = os.path.join(proj_path, gpkg_filename)
    
    cur_comids = list(set(cov_df.comid.unique()).union(set(temp_df.comid.unique())))
    cur_comids.sort()
    
    # Layer: Contributing area
    cur_layer = 'contributing_area'
    cur_contrib = contrib_areas.set_index('featureid').loc[cur_comids][['area_sqkm', 'geometry']]
    centroid = cur_contrib.union_all().centroid
    bounding_rect = cur_contrib.total_bounds
    cur_huc12s = huc12s.loc[huc12s.HUC_8 == huc8] # Constrain search area
    matched_h12s = find_containing_huc12(cur_huc12s, cur_contrib)
    cur_contrib['HUC12'] = cur_contrib.index.map(matched_h12s.set_index('featureid')['HUC12'])
    lyr2 = cur_contrib[['HUC12', 'area_sqkm', 'geometry']]
    layer_to_gpkg(cur_comids, lyr2, gpkg_path, cur_layer)
    for ix in ['featureid', 'HUC12']:
        create_index(gpkg_path, cur_layer, ix)
    
    
    # Layer: Flowlines
    cur_layer = 'flowlines'
    cur_flowlines = flowlines.set_index('comid').loc[cur_comids][['fcode', 'geometry']]
    cur_flowlines['HUC12'] = cur_flowlines.index.map(matched_h12s.set_index('featureid')['HUC12'])
    newcols = [c for c in cur_flowlines.columns if c != 'geometry'] + ['geometry']
    cur_flowlines = cur_flowlines[newcols]
    layer_to_gpkg(cur_comids, cur_flowlines, gpkg_path, cur_layer)
    for ix in ['comid', 'fcode', 'HUC12']:
        create_index(gpkg_path, cur_layer, ix)
      
        
    # Layer 2: HUC 12 geometries and area
    cur_layer = 'HUC12_boundaries'
    cur_boundaries = cur_huc12s.set_index('HUC12').loc[matched_h12s['HUC12'].unique()][['name','geometry']]
    area_sqkm = (cur_boundaries.to_crs(proj_crs).area).multiply(1e-6)
    lyr3 = pd.concat([area_sqkm.rename('area_sqkm'), cur_boundaries], axis =1)
    layer_to_gpkg(cur_comids, lyr3, gpkg_path, cur_layer, check=False)
    create_index(gpkg_path, cur_layer, 'HUC12')
    
    # Layer 3: Fcode + spatial covariates
    cur_layer = 'spatial_covariates'
    cur_spatial = spatial_covs.loc[cur_comids]
    lyr1 = pd.concat([cur_flowlines['fcode'], cur_spatial], axis = 1)
    layer_to_gpkg(cur_comids, lyr1, gpkg_path, cur_layer)
    for ix in ['comid', 'fcode']:
        create_index(gpkg_path, cur_layer, ix)
    
    
    # Calculate fraction of HUC_12 area
    cur_area = lyr2[['HUC12', 'area_sqkm']]
    total_area = cur_area.groupby('HUC12').sum()
    cur_area['h12_area'] = cur_area.HUC12.map(
        pd.Series(total_area.area_sqkm.values, index = total_area.index))
    frac_total_area = cur_area.area_sqkm.divide(cur_area.h12_area)
    
    # Layers 4-6: Seasonal temperature anomalies (7-9 by HUC12)
    cur_layer = 'stream_temperature'
    epochs = ['2010s', '2000s', '1990s']
    temp_anomalies = getSeasonalAnomalies(temp_df, baseline_period = periods[epochs[0]], 
                                        anomaly_periods = [periods[key] for key in (epochs[1], epochs[2])], cols = ['stream_temp'])
    temp_epoch_frames = parse_epochs(temp_anomalies, epochs = epochs)
    
    
    for e in epochs:
        cur_anomalies = temp_epoch_frames[e]
        layer_name = cur_layer + f'_{e}' if e=='2010s' else cur_layer + f'_anomalies_{e}'
        layer_to_gpkg(cur_comids, cur_anomalies, gpkg_path, layer_name)
        create_index(gpkg_path, layer_name, 'comid')
        add_view(gpkg_path, layer_name, 'contributing_area', 'featureid', 'comid')
        add_view(gpkg_path, layer_name, 'flowlines', 'comid', 'comid', shape = 'lines')
    
       
        #cur_awm = area_weighted_mean(cur_anomalies[[c for c in cur_anomalies.columns if c[:2]!='n_']], frac_total_area, cur_contrib['HUC12'])
        cur_awm = area_weighted_mean(cur_anomalies, frac_total_area, cur_contrib['HUC12'])
        layer_name = cur_layer + f'_{e}' + '_byHUC12' if e=='2010s' else cur_layer + f'_anomalies_{e}' + '_byHUC12'
        layer_to_gpkg(cur_comids, cur_awm, gpkg_path, layer_name, check = False)
        create_index(gpkg_path, layer_name, 'HUC12')
        add_view(gpkg_path, layer_name, 'HUC12_boundaries', 'HUC12', 'HUC12')
    
    # Layers 10-12: Seasonal covariate anomalies (13-15 by HUC12)
    cur_layer = 'temporal_covariates'
    cov_anomalies = getSeasonalAnomalies(cov_df, baseline_period = periods[epochs[0]], 
                                        anomaly_periods = [periods[key] for key in (epochs[1], epochs[2])], cols = anom_covs)
    cov_epoch_frames = parse_epochs(cov_anomalies, epochs = epochs)
    
    
    for e in epochs:
        cur_anomalies = cov_epoch_frames[e]
        layer_name = cur_layer +f'_{e}' if e=='2010s' else cur_layer + f'_anomalies_{e}'
        layer_to_gpkg(cur_comids, cur_anomalies, gpkg_path, layer_name)
        create_index(gpkg_path, layer_name, 'comid')
        add_view(gpkg_path, layer_name, 'contributing_area', 'featureid', 'comid')
        add_view(gpkg_path, layer_name, 'flowlines', 'comid', 'comid', shape = 'lines')
    
        
        cur_awm = area_weighted_mean(cur_anomalies, frac_total_area, cur_contrib['HUC12'])
        layer_name = cur_layer +f'_{e}' + '_byHUC12' if e=='2010s' else cur_layer + f'_anomalies_{e}' + '_byHUC12'
        layer_to_gpkg(cur_comids, cur_awm, gpkg_path, layer_name, check = False)
        create_index(gpkg_path, layer_name, 'HUC12')
        add_view(gpkg_path, layer_name, 'HUC12_boundaries', 'HUC12', 'HUC12')
    
    #%% Import geojson
    
    print('\n\nImporting project bounds as project_bounds.geojson ...')
    proj_bounds = os.path.join(geojson_dir, f'{curhuc}_bounds.geojson')
    dest_filename = 'project_bounds.geojson'
    try:
        bounds = shutil.copy(proj_bounds, os.path.join(proj_path, dest_filename))
        print(f'\n{dest_filename} successfully added to {os.path.relpath(proj_path, outdir)}')
    except: print(f'Failed to add {dest_filename}')
    sys.stdout =terminal_stdout
#%% Build project rs.xml

rs_project = Project(
       project_name = proj_name,
       description= f"""Predicted stream temperatures and covariates for {len(cur_comids)} reaches for each day between 1/1/1990 and 9/30/2021, 
       produced from a statistical model described in Siegel et al. (2023).""",
       meta_data=MetaData([
           Meta('Date Created',  str(datetime.datetime.now().isoformat()), type='isodate', ext=None),
           Meta('HUC', curhuc),
           Meta('Hydrologic Unit Code', curhuc),
           Meta('Watershed Name', h10name),
           Meta('reach IDs (NHDv2)', cur_comids),
       ]),
       bounds=ProjectBounds(
           Coords(centroid.x, centroid.y),
           BoundingBox(bounding_rect[0], bounding_rect[1], bounding_rect[2], bounding_rect[3]),
           os.path.basename(proj_bounds)
       ),
       realizations=[Realization(
           name='Realization1',
           xml_id='REALIZATION1',
          date_created=datetime.datetime.now(),
          product_version='1.0.0',
          datasets=[
              Geopackage(
                  name= gpkg_filename,
                  xml_id='OUTPUT',
                  path=os.path.relpath(gpkg_path, proj_path),
                  layers=get_datasets(gpkg_path)
              ),
              Dataset(
                  xml_id='README',
                  name='README',
                  description='Description of project and data files',
                  path=os.path.relpath(readme_path, proj_path),
              ),
              Dataset(
                  xml_id='LOG',
                  name='Log File',
                  description='Processing log file',
                  path=os.path.relpath(log_path, proj_path),
              ),
          ]
      )]
  )

merged_project_xml = os.path.join(proj_path, 'project.rs.xml')
rs_project.write(merged_project_xml)
print(f'Project XML file written to {merged_project_xml}')

#%%
runtime = datetime.datetime.now() - start
print(f'\n\n\n----------  DONE. Runtime {runtime.seconds // 3600}:{runtime.seconds // 60}:{runtime.seconds % 60} hours  ---------')

