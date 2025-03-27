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
import pandas as pd
import shutil
import sqlite3
import zipfile
import glob
import geopandas as gpd
from shapely import Point, LineString, Polygon
from shapely import force_2d


def create_proj_dir(huc10):
    path = os.path.join(outdir, f'retro/{huc10[:6]}/{huc10}')
    if not os.path.exists(path):
        os.makedirs(path)
        print('Project directory created.')
    else: print('Project directory exists.')
    return(path)

def add_cov_metadata(path):
    cov_metadata = os.path.join('DATA', 'covariate_metadata.csv')
    shutil.copy(cov_metadata, os.path.join(path, 'covariate_metadata.csv'))
    return(print('Covariate metadata successfully added.'))
    
def add_database(df, db_name, db_path, zip_file = True):
    connection= sqlite3.connect(db_path)
    df.to_sql(db_name, connection, if_exists = 'replace')
    connection.close()
    print('Database created succesfully')
    if zip_file:
        zip_path = db_path + '.zip'
        try:
            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(db_path, os.path.basename(db_path))
            os.remove(db_path)
            print(f"Successfully zipped '{db_name}' to '{zip_path}'")
        except:
            print(f"Error zipping database {db_name}.")
    return

def create_database_file(file, date_col, proj_path, db_type, compression = None, add_cols = None, zip_file = True, overwrite = False):
    df = pd.read_csv(file, parse_dates = [date_col], compression = compression)
    df.columns = df.columns.str.split(r'\.').str[-1].str.strip()
    out_path = os.path.join(proj_path, f'daily_{db_type}.db')
    if not (os.path.exists(out_path)) or (overwrite==True):
        if db_type == 'stream_temperature':
            try:
                dropcols = ['antec_air_temp', 'std_mean_flow']
                df_trim = df.drop(dropcols, axis = 1, errors = 'ignore')
                df_out = df_trim[['COMID', 'date', 'stream_temp']]
                print('Stream temperature file processed.')
            except: print('Stream temperature file failed.')
       
        else: 
            try:
                df = df.loc[:,~df.columns.duplicated()].copy()
                df_out = pd.merge(df, add_cols, how = 'left', on=['COMID','date'])
                print('Temporal covariates file processed.')
            except: print('Daily covariates file failed.')
       
        out_table = df_out.fillna(-999).round(4)
        print('Writing to database...')
        add_database(out_table, db_type, db_path = out_path, zip_file = zip_file)
    return(df)

def check_comids(comids, data, layer_name):
    try: no_data = set(comids).difference(set(data.reset_index().COMID.unique()))
    except: no_data = set(comids).difference(set(data.reset_index().FEATUREID.unique()))
    if no_data==set(): 
        return(print(f'{layer_name} successfully loaded for all comids'))
    else: 
        return(print('{layer_name} missing for some comids:', no_data))

def layer_to_gpkg(comids, data, path, layer_name, epsg = 'EPSG:4326'):
    check_comids(comids, data, layer_name)
    try: 
        data.to_crs(epsg = epsg, inplace = True)
    except: pass

    if not os.path.exists(path):
        data.to_file(path, layer = layer_name, driver = 'GPKG')
        print('Geopackage successfully initiated.')
    else: gpd.GeoDataFrame(data).to_file(path, layer = layer_name)
    
    return(print(f'{layer_name} added to {path}.'))
    
def load_spatial_covariates(file):
    df = pd.read_csv(file, index_col = 'COMID')
    df.columns = df.columns.str.split(r'\.').str[-1].str.strip()
    return(df)

def find_containing_huc12(h12s, wings):
    h12s = h12s.dissolve(by='HUC_12')
    
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
    if fail_list == []: print('All contributing areas successfully matched to HUC 12s.')
    else: print(f'Error finding HUC 12s containing {fail_list}')
    return(pd.DataFrame(zip(wings_list, h12s_list), columns = ['FEATUREID', 'HUC_12']))
        
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
    print('Calculating seasonal anomalies...')
    try:
        cur_huc_anoms = []
        cur_comids = list(indf.COMID.unique())
        for j in cur_comids:
            #Assign water year and season
            cur_df = indf[indf.COMID == j][['date']+ cols]
            cur_ens_med = cur_df.groupby('date').median() # sets index to date
            cur_ens_med['doy'] = cur_ens_med.index.dayofyear
            cur_ens_med['season'] = cur_ens_med['doy'].apply(assign_season)
            for col in cols: cur_ens_med[f'n_{col}'] = cur_ens_med[col].notnull().astype(int)
    
            if 'cov.SWE' in cols:
                #Use water year for swe
                cur_swe = cur_ens_med.copy()
                cur_swe['WY'] = cur_swe.index.to_series().apply(assign_water_year)
        
                # Reclassify end of Sept as summer (pull it into summer of WY[y-1])
                cur_swe.loc[cur_swe.index.month == 9, 'season'] = 'summer'        
                swe_cum = cur_swe[['cov.SWE', 'WY', 'season', 'n_cov.SWE']].groupby(['WY', 'season']).sum()
    
            
            cur_anomalies = {}
            
            #Determine baseline (median)
            cur_met = cur_ens_med.copy()
            baseline = cur_met[cur_met.index.year.isin(baseline_period)].groupby('season').median()
            for col in cols: 
                baseline[f'n_{col}'] = cur_met[cur_met.index.year.isin(baseline_period)].groupby('season').sum()[f'n_{col}'] # n
                
            if 'cov.SWE' in cols:
                swe_df = swe_cum.loc[pd.IndexSlice[baseline_period], :, :]                
                baseline['cov.SWE'] = swe_df.groupby(level=1).median()['cov.SWE']
                baseline['n_cov.SWE'] = swe_df.groupby(level=1).sum()['n_cov.SWE'] # n
            
            
            baseline_trim = pd.concat([baseline[[col, f'n_{col}']] for col in cols], axis = 1)
            
            cur_anomalies[f'{baseline_period[0]}s'] = baseline_trim
    
            #Calculate anomalies
            for per in anomaly_periods:                
                cur_anoms = cur_ens_med[cur_ens_med.index.year.isin(per)].groupby('season').median().subtract(baseline)
                for col in cols: 
                    cur_anoms[f'n_{col}'] = cur_ens_med[cur_ens_med.index.year.isin(per)].groupby('season').sum()[f'n_{col}'] # n
    
                if 'cov.SWE' in cols:
                    swe = swe_cum.loc[pd.IndexSlice[per], :, :]                
                    cur_anoms['cov.SWE'] = swe.groupby(level=1).median()['cov.SWE'].subtract(baseline['cov.SWE'])
                    cur_anoms['n_cov.SWE'] = swe.groupby(level=1).sum()['n_cov.SWE'] # n
                
                cur_anoms_trim = pd.concat([cur_anoms[[col, f'n_{col}']] for col in cols], axis = 1)
                
                cur_anomalies[f'{per[0]}s'] = cur_anoms_trim
            
            cur_seas_anomalies = pd.concat(cur_anomalies, axis = 1)
            cur_huc_anoms.append(cur_seas_anomalies)
        print('Done.')
        res = (pd.concat(cur_huc_anoms, keys = cur_comids))
        res.index.names = ['COMID', 'season']
        return res.reset_index()
    except: 
        return(print('Error calculating anomalies.'))
    
def parse_epochs(indf: pd.DataFrame, epochs: list):
    indf.reset_index()
    epoch_dict = {}
    idx = indf[['COMID', 'season']]
    idx.columns = idx.columns.droplevel(1)
    df = pd.DataFrame(idx)
    for e in epochs:
        epoch_dict[e] = pd.concat([df, indf.xs(e, axis = 1).copy()], axis = 1)
    return(epoch_dict)
#%%
# USER INPUTS
# Required
curhuc = '1701010107' #huc10
huc8 = curhuc[:8]

# Optional
datadir = os.path.join((os.path.dirname(__file__)), 'DATA/')
overwrite = False

# Define anomaly periods
periods = {'2000s': list(range(2000, 2010)),
           '2010s': list(range(2010, 2020)),
           '2020s': list(range(2020, 2030)), 
           '1990s': list(range(1990, 2000)), 
           '2050s': list(range(2050, 2060)), 
           '2080s': list(range(2080, 2090))}

#if name=='__main__':
#%% Define directories 

tempdir = os.path.join(datadir, 'preds_retro/')
outdir = os.path.join(datadir, 'Outputs/')
nhddir = os.path.join(datadir, 'NHDPlusPN/NHDPlus17/')

#%% Load in geometry files

hydro_files = glob.glob(os.path.join(nhddir, r'NHDSnapshot/Hydrography/*.dbf'))

flowlines = gpd.read_file(hydro_files[0])
flowlines['geometry'] = flowlines['geometry'].apply(lambda x: force_2d(x))

spatial_covs_file = os.path.join(datadir, 'spatial_data.csv')
spatial_covs = load_spatial_covariates(spatial_covs_file)

contrib_area_file = os.path.join(nhddir, 'NHDPlusCatchment/Catchment.dbf')
contrib_areas = gpd.read_file(contrib_area_file)

huc12s_file = os.path.join(nhddir, 'WBDSnapshot/WBD/WBD_Subwatershed.dbf')
huc12s = gpd.read_file(huc12s_file, columns = ['OBJECTID', 'HUC_12', 'HUC_8','geometry'])
#%%

# Create project diorectory
proj_path = create_proj_dir(curhuc) 
                            
#%%
# Add daily temperature and covariate retrospective databases

#Temperature
temp_file = os.path.join(tempdir, 'predictions_temperature', f'st_pred_{curhuc}.csv')
temp_df = create_database_file(temp_file, date_col = 'tim.date', proj_path = proj_path, 
                               db_type = 'stream_temperature', add_cols = None, compression = None, 
                               overwrite = overwrite)

#Covariates
add_cov_metadata(proj_path)

cov_file = os.path.join(tempdir, 'predictions_covariates', 'cov_csvs', f'{curhuc}_covs.zip')
cov_df = create_database_file(cov_file, date_col = 'date', proj_path = proj_path,
                              db_type = 'covariates', add_cols = temp_df[['COMID', 'date', 'antec_air_temp', 'std_mean_flow']], 
                              compression = 'zip', overwrite = overwrite)

#%%
# Build geopackage
gpkg_path = os.path.join(proj_path, 'seasonal_anomalies_spatial_covariates.gpkg')

cur_comids = list(set(cov_df.COMID.unique()).union(set(temp_df.COMID.unique())))
cur_comids.sort()

# Layer 0: Flowlines
cur_layer = 'flowlines'
cur_flowlines = flowlines.set_index('COMID').loc[cur_comids][['FCODE', 'geometry']]
layer_to_gpkg(cur_comids, cur_flowlines, gpkg_path, cur_layer)

# Layer 1: Fcode + spatial covariates
cur_layer = 'spatial_covarites'
cur_spatial = spatial_covs.loc[cur_comids]
lyr1 = pd.concat([cur_flowlines['FCODE'], cur_spatial], axis = 1)
layer_to_gpkg(cur_comids, lyr1, gpkg_path, cur_layer)

# Layer 2: Contributing area
cur_layer = 'contriubting_area'
cur_contrib = contrib_areas.set_index('FEATUREID').loc[cur_comids][['AreaSqKM', 'geometry']]
cur_huc12s = huc12s.loc[huc12s.HUC_8 == huc8]
matched_h12s = find_containing_huc12(cur_huc12s, cur_contrib)
cur_contrib['HUC_12'] = cur_contrib.index.map(matched_h12s.set_index('FEATUREID')['HUC_12'])
lyr2 = cur_contrib[['AreaSqKM', 'HUC_12', 'geometry']]
layer_to_gpkg(cur_comids, lyr2, gpkg_path, cur_layer)

# Layers 3-5 Seasonal temperature anomalies
cur_layer = 'stream_temperature'
epochs = ['2010s', '2000s', '1990s']
temp_anomalies = getSeasonalAnomalies(temp_df, baseline_period = periods[epochs[0]], 
                                    anomaly_periods = [periods[key] for key in (epochs[1], epochs[2])], cols = ['stream_temp'])
temp_epoch_frames = parse_epochs(temp_anomalies, epochs = epochs)
for e in epochs:
    cur_anomalies = temp_epoch_frames[e]
    layer_name = cur_layer +f'_{e}' if e=='2010s' else cur_layer + f'_anomalies_{e}'
    layer_to_gpkg(cur_comids, cur_anomalies, gpkg_path, layer_name)

# Layers 6-8 Seasonal covariate anomalies
cur_layer = 'temporal_covariates'
cov_anomalies = getSeasonalAnomalies(cov_df, baseline_period = periods[epochs[0]], 
                                    anomaly_periods = [periods[key] for key in (epochs[1], epochs[2])], cols = ['air_temp_ws', 'NWM_flow_log', 'SWE'])
cov_epoch_frames = parse_epochs(cov_anomalies, epochs = epochs)
for e in epochs:
    cur_anomalies = cov_epoch_frames[e]
    layer_name = cur_layer +f'_{e}' if e=='2010s' else cur_layer + f'_anomalies_{e}'
    layer_to_gpkg(cur_comids, cur_anomalies, gpkg_path, layer_name)
    
# Layer 9: HUC 12 geometries and area

