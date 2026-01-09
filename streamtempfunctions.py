#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 14:30:58 2025

@author: dawn urycki
"""

'''
This script does stuff with stuff
'''

import geopandas as gpd
import os
import numpy as np
import pandas as pd
import xarray as xr
import dask.dataframe as dd
import shutil
import sqlite3
import gc
import pyarrow as pa
import pyarrow.parquet as pq
import zipfile
from os.path import basename

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

def create_proj_dir(huc6, new_dir, ref_path):
    path = os.path.join(ref_path, f'{huc6}/{new_dir}')
    
    if not os.path.exists(path):
        os.makedirs(path)
        print(f'Project directory created: {os.path.relpath(path, ref_path)}')
    else: 
        print(f'Project directory exists: {os.path.relpath(path, ref_path)}')
        #path = path + '_1'
        #os.makedirs(path)
    return(path)

def add_file(filepath, dest, ref_path = None):
    filename =  filepath.split("/")[-1]
    shutil.copy(filepath, dest)
    outpath = os.path.abspath(os.path.join(dest, filename))
    print(f'\n{filename} successfully added to {os.path.relpath(dest, ref_path)}.')
    return outpath

def add_dataset(df, ds_name, ds_path, ref_path = None):
    try: 
        df_filled = df.drop('GCM', axis =1).fillna(-999).round(4)
        df_filled['GCM'] = df['GCM']
    except:
        df_filled = df.copy()
        for col in df.columns: 
            if col not in ['GCM', 'lookup']:
                df.loc[:,col] = df.loc[:,col].fillna(-999).round(4)
    try:
        df.to_parquet(ds_path, compute = True, engine = 'pyarrow', write_index = False, compression = 'gzip')
        shutil.make_archive(f'{ds_path}', 'zip', ds_path)
        shutil.rmtree(ds_path)
        print(f"Successfully created '{os.path.relpath(ds_path, ref_path)}'")
    except:
        try:
            try: df['lookup'] = df['lookup'].astype(str)
            except: pass
            
            writer = None
            for i in range(0, len(df), 10000):
                chunk_df = df.iloc[i:i + 10000]
                table = pa.Table.from_pandas(chunk_df)

                if writer is None:
                    # Create the Parquet writer with GZIP compression and desired row group size
                    writer = pq.ParquetWriter(ds_path, table.schema, compression='gzip')
                # Write the current chunk (table) to the Parquet file
                writer.write_table(table)

            # Close the writer to finalize the Parquet file
            if writer is not None:
                writer.close()
            
            #df.to_parquet(ds_path, engine = 'pyarrow', index = False, compression = 'gzip')
            print(f"Successfully created '{os.path.relpath(ds_path, ref_path)}'")
            df = df.loc[df.GCM=='ens_mean']
        
        except: print(f"\nERROR creating {ds_path}.")
    
    df_out = df[df.GCM=='ens_mean'].replace(-999, np.nan)
    df_out = df_out.compute()
    
    return df_out

def create_database_file(files, cur_comids, date_col, proj_path, db_type, 
                         add_cols = None, overwrite = False, compression = None,
                         covs_to_count = ['antec_air_temp', 'NWM_flow_log', 'SWE']):
    frames = []
    for file in files:
        if not os.path.exists(file):
            print(f'File not found: {os.path.relpath(file)}')
            if (frames == []) & (file == files[-1]):
                print('NO DATA.')
                return pd.DataFrame(), [], []
            continue
        
        
        print(f'\n\nRetrieving data from {os.path.relpath(file)}')
        
        try:
            if file[-3:]=='.nc': 
                ds = xr.open_dataset(file)
                df = ds[['COMID','date', 'GCM', *ds.data_vars]].to_dask_dataframe()
                                   # ens_mean = ds.mean(dim='GCM')
                # ens_mean = ens_mean.assign_coords(GCM='ens_mean')
                # ds_comb = xr.concat([ds, ens_mean], dim = 'GCM')
                # #df_mean = ens_mean[['COMID','date', 'GCM', *ens_mean.data_vars]].to_dataframe().reset_index()
                # cur_df = ds_comb[['COMID','date', 'GCM', *ds_comb.data_vars]].to_dataframe().reset_index()
                
            elif file[-5:] == '.parq':
                df = dd.read_parquet(file, filters = [('COMID', 'in', (tuple(cur_comids)))], engine = 'pyarrow')
                df.columns = [c.strip('cov.') for c in df.columns]
                df['date'] = dd.to_datetime(df['tim.year'] * 1000 + df['tim.doy'], format='%Y%j')
                #ds = df.to_xarray()
            
            print(f'Calculating ensemble mean {db_type}...')
            
            try: df = df.drop('lookup', axis = 1)
            except: pass
        
            ens_mean = df.groupby(['COMID', 'date']).mean().reset_index()
            ens_mean['GCM'] = 'ens_mean'
            df['GCM'] = df['GCM'].astype('str')
            cur_df = dd.concat([df, ens_mean], axis = 0)
            frames.append(cur_df)
            print('Done.')
               
        except: print('ERROR. Ensemble mean could not be calculated.')            
        
    df = dd.concat(frames, axis = 0)
    gc.collect()           

    try:
        if not pd.api.types.is_datetime64_any_dtype(df['tim.date']): 
            df = df.rename(columns = {'tim.date': 'r_date'})    
    except:pass
    
    df.columns = [col.split('.')[-1] for col in list(df.columns)]
    df = df.rename(columns = {'COMID': 'comid', 'Huc10': 'HUC10', 'Flow_log': 'flow_log', 'Flow': 'flow'})
    out_path = os.path.join(proj_path, f'daily_{db_type}.parquet')

    if db_type == 'stream_temperature':
        try:
            datapoints = np.array(df.stream_temp.count(axis = 0).compute())
            print(f'''\n\nStream temperature file processed: {len(df.comid.unique().compute())} COMIDs, {df.shape[0].compute()} records, {datapoints} non-null stream temperature values.''')
        except: print('\n\nStream temperature file FAILED.')
   
    else: 
        try:
            ens_mean = ens_mean.rename(columns = {'COMID': 'comid', 'Huc10': 'HUC10', 'Flow_log': 'flow_log', 'Flow': 'flow'})
            datapoints = ens_mean[covs_to_count].count(axis=0).compute().values * 10
            print(f'\n\nTemporal covariates file processed: {len(ens_mean['comid'].unique().compute())} COMIDs;\nNumber covariate predictions {covs_to_count}: {datapoints}.')
        except: 
            print('\n\nTemporal covariates file FAILED.')
            return pd.DataFrame(), [], []
        
    if datapoints.sum() > 0:
        print('Writing to database...')
        try: 
            df = add_dataset(df, db_type, ds_path = out_path, ref_path = proj_path)
            print(f"Successfully created {os.path.relpath(out_path+'.zip', proj_path)}.")
        except: print(f"Failed to write {os.path.relpath(out_path+'.zip', proj_path)}.")
    return df, datapoints, out_path + '.zip'  


def check_comids(comids, data, layer_name):
    n_cols = [col for col in data.columns if col[0:2] == 'n_']
    if n_cols != []: data_comids = set(data[data[n_cols].sum(axis =1) != 0].index)
    else: data_comids = set(data.index)
    
    return(print(f'\n\n{layer_name} successfully loaded for {len(comids)} features'))

def layer_to_gpkg(comids, data, path, rel_path, layer_name, epsg = '4326', check=True):
    if check: check_comids(comids, data, layer_name)
    try: 
        data.to_crs(epsg = epsg, inplace = True)
    except: pass

    if not os.path.exists(path):
        data.fillna(-999).infer_objects(copy=False).round(4).to_file(path, layer = layer_name, driver = 'GPKG')
        print('\n\nGeopackage successfully initiated.')
    else: gpd.GeoDataFrame(data.fillna(-999).infer_objects(copy=False).round(4)).to_file(path, layer = layer_name)
    
    return(print(f'\n{layer_name} added to {os.path.relpath(path, rel_path)}.'))
    
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
    
def getSeasonalAnomalies(indf, baseline_period: range, anomaly_periods: list, cols: list): # Add comid to groupby list to avoid looping over comids?
    print('\n\nCalculating seasonal anomalies...')
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
                swe_cum = cur_swe[['SWE', 'WY', 'season', 'n_SWE']].groupby(['WY', 'season']).sum(min_count=1)
                
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
        print(f'\t{new_cols} anomalies done.')
        res = (pd.concat(cur_huc_anoms, keys = cur_comids))
        res.index.names = ['comid', 'season']
        return res.reset_index()
    except: 
        return(print('Error calculating anomalies.')) 
    
def collapse_levels(indf, level_name):
    outdf = indf[['spring', 'summer', 'fall', 'winter']]
    #value_if_true if condition else value_if_false
    outdf.columns = [f'{c[:2]}_{level_name}' if level_name[:2] != 'n_' else f'n_{c[:2]}_{level_name[2:]}' for c in outdf.columns]
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

def area_weighted_mean(metric_df, area_df, agg = 'HUC12'):
    outframes = []
    metrics = [col for col in metric_df.columns if col[:2] != 'n_']
    for m in metrics:
        cur_met = pd.concat([metric_df[[m, f'n_{m}']].dropna(), area_df], join = 'inner', axis = 1)
        # Calculate fraction of HUC_12 area
        total_area = area_df.loc[cur_met.index].groupby(agg).sum()
        cur_met['h12_area'] = cur_met.HUC12.map(total_area.area_sqkm)
        cur_met['weight'] = cur_met.area_sqkm.divide(cur_met.h12_area)
    
        metric_comps = cur_met[[m, f'n_{m}']].mul(cur_met.weight, axis  =0)
        metric_comps[agg] = metric_comps.index.map(area_df[agg])
        outframes.append(metric_comps.groupby(agg).sum(min_count=1)[[m, f'n_{m}']])
    weighted_mean = pd.concat(outframes, axis = 1)
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
        except: return(print(f'\nFailed to create index on {ixname}'))

def add_view(dbpath, layer_name, feature_table, f_index, a_index, shape = 'polygons'):
    conn = sqlite3.connect(dbpath)
    with conn:
        curs = conn.cursor()
        try: 
            curs.execute(f'''
                         CREATE VIEW vw_{layer_name}_{shape} AS 
                         SELECT a.*, f.geom 
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
       
def check_bounds(rsx_bounds, cur_geoms):    
    if len(gpd.overlay(rsx_bounds, cur_geoms.to_crs(rsx_bounds.crs).dissolve(), how = 'intersection')) < 1:
        print('Not overlapping')
        return(False)
    else: 
        print('Overlapping')
        return(True)
        
def lyr_name_to_sentence(text):
   if text[0:3] == 'vw_': text = text[3:]
   return text[0].upper() + text[1:].replace('_', ' ')

def get_datasets(output_gpkg: str):
   """
   Returns a list of the datasets from the output GeoPackage.
   """

   conn = sqlite3.connect(output_gpkg)
   #conn.enable_load_extension(True)
   curs = conn.cursor()
   
   # Get the names of all the tables in the database
   curs.execute("SELECT table_name FROM gpkg_contents")
   datasets = [GeopackageLayer(
       lyr_name=row[0],
       ds_type=GeoPackageDatasetTypes.VECTOR,
       name=lyr_name_to_sentence(row[0])
   ) for row in curs.fetchall()]
   return datasets