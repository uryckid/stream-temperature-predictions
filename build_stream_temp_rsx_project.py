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
import pandas as pd
import shutil
import glob
import geopandas as gpd
from shapely import force_2d
from collections import Counter

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
import numpy as np
import sys
import subprocess

import streamtempfunctions as stp   


#%%
# USER INPUTS
proj_type = 'Stream_Temperature_Retrospective_'
timeframe = '1990-2021_'

# Required
proj_crs = 'EPSG:9674' # NAD83 / USFS R6 Albers projection, Oregon and Washington

# Optional
os.chdir(os.path.dirname(__file__))
curdir = os.getcwd()
datadir = os.path.join(curdir, 'DATA/')
geojson_dir = glob.glob(os.path.join(datadir, '*HUC17_project_bounds/' ))[0]
overwrite = True

# Define anomaly periods
periods = {'2000s': list(range(2000, 2010)),
           '2010s': list(range(2010, 2020)),
           '2020s': list(range(2020, 2030)), 
           '1990s': list(range(1990, 2000)), 
           '2050s': list(range(2050, 2060)), 
           '2080s': list(range(2080, 2090))}

anom_covs = ['air_temp_ws', 'NWM_flow_log', 'SWE']
anom_covs = ['antec_air_temp', 'NWM_flow_log', 'SWE']



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
spatial_covs = stp.load_spatial_covariates(spatial_covs_file)

contrib_area_file = os.path.join(nhddir, 'NHDPlusCatchment/Catchment.dbf')
contrib_areas = gpd.read_file(contrib_area_file)
contrib_areas.rename({'FEATUREID': 'featureid', 'AreaSqKM': 'area_sqkm'}, axis = 1, inplace = True)
#%%
huc12s_file = os.path.join(nhddir, 'WBDSnapshot/WBD/WBD_Subwatershed.dbf')
huc12s = gpd.read_file(huc12s_file, columns = ['OBJECTID', 'HUC_12', 'HUC_8','HU_12_NAME','HUC_10','HU_10_NAME','geometry'])
huc12s.rename({'HUC_12': 'HUC12', 'HU_12_NAME': 'name'}, axis = 1, inplace = True)

xwalk = pd.read_csv(os.path.join(datadir, 'COMID_to_HUC12.csv'), index_col = 'COMID')
xwalk.fillna({'Huc12':0}, inplace=True)
xwalk.Huc12 = xwalk.Huc12.astype(int)
xwalk_usgs = pd.read_csv(os.path.join(datadir, 'CrosswalkTable_NHDplus_HU12.csv'), index_col='FEATUREID')
xwalk['Huc12_usgs'] = xwalk_usgs['HUC_12'].astype(int)
h10s_af = list({str(h)[:10] for h in xwalk.Huc12.unique()})
h10s_usgs = list({str(h)[:10] for h in xwalk.Huc12_usgs.unique()})
h10s_usgs.sort()
#conflicts = xwalk[xwalk.Huc12 != xwalk.Huc12_usgs]

for curhuc in [h10s_usgs[5]]: 
    print(curhuc)
    
    cur_comids = list(xwalk.loc[xwalk.Huc12_usgs.astype(str).str.contains(curhuc)].index)
    cur_comids.sort()
    
    alth10s = list({str(c)[:10] for c in xwalk.loc[cur_comids].Huc12.unique()})
    legacy_h10s = [h for h in alth10s if h != curhuc]
        
    h10name_hr = huc12s.loc[huc12s.HUC_10 == curhuc]['HU_10_NAME'].mode()[0]
    h10name = h10name_hr.replace(" ","_")
    proj_name_hr = f'{proj_type.replace('_', ' ')}({timeframe[:-1]}) for {h10name_hr}'
    #%% Create project directory and add daily temperature and covariate databases
    proj_name = proj_type + timeframe + h10name
    proj_path = stp.create_proj_dir(curhuc[:6], proj_name, outdir)
    log_path = os.path.join(proj_path, 'processing.log')
    terminal_stdout = sys.stdout
    with open(log_path, 'w') as f:
        #sys.stdout = f
        #sys.stdout = open(log_path,'w')
        start = datetime.datetime.now()
        print(start.strftime("%Y-%m-%d %H:%M:%S"))
        print(f'\nBuilding {proj_type[:-1]} for {h10name} ({curhuc})')
        
        # Add README
        readme_path = os.path.join(outdir, 'README.html')
        readme_out_path = stp.add_file(readme_path, proj_path, outdir)
        
        # # Add covariate metadata
        # cov_meta_path = os.path.join('DATA', 'covariate_metadata.csv')
        # cov_meta_out_path = stp.add_file(cov_meta_path, proj_path, outdir)
        
        # Temperature
        temp_files = [os.path.join(tempdir, 'predictions_temperature', f'st_pred_{h}.csv') for h in alth10s]
        temp_df, no_st_values, st_db_path = stp.create_database_file(temp_files, cur_comids, date_col = 'tim.date', proj_path = proj_path, 
                                       db_type = 'stream_temperature', add_cols = None, compression = None, 
                                       overwrite = overwrite)
        
        # Covariates
        cov_files = [os.path.join(tempdir, 'predictions_covariates', 'cov_csvs', f'{h}_covs.zip') for h in alth10s]
        cov_df, no_cov_values, cov_db_path = stp.create_database_file(cov_files, cur_comids, date_col = 'date', proj_path = proj_path,
                                      db_type = 'covariates', add_cols = None, compression = 'zip', overwrite = overwrite)

        #%% Build geopackage

        gpkg_filename = 'seasonal_anomalies_spatial_covariates.gpkg'
        gpkg_path = os.path.join(proj_path, gpkg_filename)
    
        
        # Layer: Contributing area
        cur_layer = 'contributing_area'
        cur_contrib = contrib_areas.set_index('featureid').loc[cur_comids][['area_sqkm', 'geometry']]
        centroid = cur_contrib.union_all().centroid
        bounding_rect = cur_contrib.total_bounds
        cur_contrib['HUC12'] = cur_contrib.index.map(xwalk.Huc12_usgs)
        lyr2 = cur_contrib[['HUC12', 'area_sqkm', 'geometry']]
        stp.layer_to_gpkg(cur_comids, lyr2, gpkg_path, cur_layer)
        for ix in ['featureid', 'HUC12']:
            stp.create_index(gpkg_path, cur_layer, ix)
        
        
        # Layer: Flowlines
        cur_layer = 'flowlines'
        cur_flowlines = flowlines.set_index('comid').loc[cur_comids][['fcode', 'geometry']]
        cur_flowlines['HUC12'] = cur_flowlines.index.map(xwalk.Huc12_usgs)
        newcols = [c for c in cur_flowlines.columns if c != 'geometry'] + ['geometry']
        cur_flowlines = cur_flowlines[newcols]
        stp.layer_to_gpkg(cur_comids, cur_flowlines, gpkg_path, cur_layer)
        for ix in ['comid', 'fcode', 'HUC12']:
            stp.create_index(gpkg_path, cur_layer, ix)
          
            
        # Layer 2: HUC 12 geometries and area
        cur_layer = 'HUC12_boundaries'
        cur_boundaries = huc12s.set_index('HUC12').loc[cur_contrib.HUC12.unique().astype(str)][['name','geometry']]
        area_sqkm = (cur_boundaries.to_crs(proj_crs).area).multiply(1e-6)
        lyr3 = pd.concat([area_sqkm.rename('area_sqkm'), cur_boundaries], axis =1)
        stp.layer_to_gpkg(cur_comids, lyr3, gpkg_path, cur_layer, check=False)
        stp.create_index(gpkg_path, cur_layer, 'HUC12')
        
        # Layer 3: Fcode + spatial covariates
        cur_layer = 'spatial_covariates'
        cur_spatial = spatial_covs.loc[cur_comids]
        lyr1 = pd.concat([cur_flowlines['fcode'], cur_spatial], axis = 1)
        stp.layer_to_gpkg(cur_comids, lyr1, gpkg_path, cur_layer)
        for ix in ['comid', 'fcode']:
            stp.create_index(gpkg_path, cur_layer, ix)
        
        
        # Layers 4-6: Seasonal temperature anomalies (7-9 by HUC12)
        cur_layer = 'stream_temperature'
        epochs = ['2010s', '2000s', '1990s']
        temp_anomalies = stp.getSeasonalAnomalies(temp_df, baseline_period = periods[epochs[0]], 
                                            anomaly_periods = [periods[key] for key in (epochs[1], epochs[2])], cols = ['stream_temp'])
        temp_epoch_frames = stp.parse_epochs(temp_anomalies, epochs = epochs)
        
        
        for e in epochs:
            cur_anomalies = temp_epoch_frames[e]
            layer_name = cur_layer + f'_{e}' if e=='2010s' else cur_layer + f'_anomalies_{e}'
            stp.layer_to_gpkg(cur_comids, cur_anomalies, gpkg_path, layer_name)
            stp.create_index(gpkg_path, layer_name, 'comid')
            stp.add_view(gpkg_path, layer_name, 'contributing_area', 'featureid', 'comid')
            stp.add_view(gpkg_path, layer_name, 'flowlines', 'comid', 'comid', shape = 'lines')
        
           
            cur_awm = stp.area_weighted_mean(cur_anomalies,  lyr2[['HUC12', 'area_sqkm']])
            layer_name = cur_layer + f'_{e}' + '_byHUC12' if e=='2010s' else cur_layer + f'_anomalies_{e}' + '_byHUC12'
            stp.layer_to_gpkg(cur_comids, cur_awm, gpkg_path, layer_name, check=False)
            stp.create_index(gpkg_path, layer_name, 'HUC12')
            stp.add_view(gpkg_path, layer_name, 'HUC12_boundaries', 'HUC12', 'HUC12')

#%%        
        # Layers 10-12: Seasonal covariate anomalies (13-15 by HUC12)
        cur_layer = 'temporal_covariates'
        cov_anomalies = stp.getSeasonalAnomalies(cov_df, baseline_period = periods[epochs[0]], 
                                            anomaly_periods = [periods[key] for key in (epochs[1], epochs[2])], cols = anom_covs)
        cov_epoch_frames = stp.parse_epochs(cov_anomalies, epochs = epochs)
        
        
        for e in epochs:
            cur_anomalies = cov_epoch_frames[e]
            layer_name = cur_layer +f'_{e}' if e=='2010s' else cur_layer + f'_anomalies_{e}'
            stp.layer_to_gpkg(cur_comids, cur_anomalies, gpkg_path, layer_name)
            stp.create_index(gpkg_path, layer_name, 'comid')
            stp.add_view(gpkg_path, layer_name, 'contributing_area', 'featureid', 'comid')
            stp.add_view(gpkg_path, layer_name, 'flowlines', 'comid', 'comid', shape = 'lines')
        
            
            cur_awm = stp.area_weighted_mean(cur_anomalies, lyr2[['HUC12', 'area_sqkm']])
            layer_name = cur_layer +f'_{e}' + '_byHUC12' if e=='2010s' else cur_layer + f'_anomalies_{e}' + '_byHUC12'
            stp.layer_to_gpkg(cur_comids, cur_awm, gpkg_path, layer_name, check=False)
            stp.create_index(gpkg_path, layer_name, 'HUC12')
            stp.add_view(gpkg_path, layer_name, 'HUC12_boundaries', 'HUC12', 'HUC12')
        out_anoms = list(dict.fromkeys([c[:-3] for c in cur_anomalies.columns if c[:2]!= 'n_']))
        
        
        #%% Import geojson
        
        print('\nChecking geometries in HUC 10 project boundaries...')
        alth10s.append(curhuc)
        
        for h10 in alth10s:
            try: proj_bounds = os.path.join(geojson_dir, f'{h10}_bounds.geojson') 
            except: continue
            print(f'Checking RSX project bounds for {h10}...')
            if stp.check_bounds(gpd.read_file(proj_bounds), cur_contrib): break
        
        print('\n\nImporting project bounds as project_bounds.geojson ...')
        print('Done.')
                    
        dest_filename = 'project_bounds.geojson'
        try:
            bounds = shutil.copy(proj_bounds, os.path.join(proj_path, dest_filename))
            print(f'\n{dest_filename} successfully added to {os.path.relpath(proj_path, outdir)}')
        except: print(f'Failed to add {dest_filename}')
        
    
        #%% Build project rs.xml
        
        rs_project = Project(
               name = proj_name_hr,
               project_type = 'streamtemp',
               description= f"""Predicted stream temperatures and covariates for {len(cur_comids)} reaches for each day between 1/1/1990 and 9/30/2021, 
               produced from a statistical model described in Siegel et al. (2023).""",
               meta_data=MetaData([
                   Meta('Date Created',  str(datetime.datetime.now().isoformat()), type='isodate', ext=None),
                   Meta('Hydrologic Unit Code (HUC)', curhuc),
                   Meta('ModelVersion', '1.0.0'),
                   Meta('Model Documentation','https://doi.org/10.1371/journal.pwat.0000119', type='url'),
                   Meta('README', 'https://noaa-nwfsc.github.io/stream-temperature-predictions/', type='url'),
                   Meta('Watershed Name', h10name.replace('_', ' ')),
                   Meta('Legacy HUC ID(s)', legacy_h10s),
                   Meta('Number stream temperature predictions', [int(no_st_values)]),
                   Meta(f'Number covariate predictions {out_anoms}', [int(v) for v in no_cov_values]),
                   Meta('Feature IDs (NHDv2 COMIDs)', [int(c) for c in cur_comids])
               ]),
               bounds=ProjectBounds(
                   Coords(centroid.x, centroid.y),
                   BoundingBox(bounding_rect[0], bounding_rect[1], bounding_rect[2], bounding_rect[3]),
                   os.path.relpath(bounds, proj_path)
               ),
               realizations=[Realization(
                   name='Realization1',
                   xml_id='REALIZATION1',
                  date_created=datetime.datetime.now(),
                  product_version='1.0.0',
                  datasets=[
                      Geopackage(
                          name= 'Seasonal anomalies and spatial covariates', 
                          xml_id='OUTPUT',
                          path=os.path.relpath(gpkg_path, proj_path),
                          layers=stp.get_datasets(gpkg_path)
                      ),
                      Dataset(
                          xml_id='README',
                          name='README',
                          ds_type = 'File',
                          description='Project description and metadata',
                          path=os.path.relpath(readme_out_path, proj_path),
                      ),
                      # Dataset(
                      #     xml_id='COVARIATE_METADATA',
                      #     name='Covariate metadata',
                      #     ds_type = 'File',
                      #     description='Temporal and spatial covariates considered or used in the final temperature predication model, with their symbology, definitions, and sources.',
                      #     path=os.path.relpath(cov_meta_out_path, proj_path),
                      # ),
                      Dataset(
                          xml_id='DAILY_COVARIATES',
                          name='Daily covariate predictions',
                          ds_type = 'File',
                          description='Zipped database of the full set of daily temporal covariate values, by reach, for covariates used in the temperature prediction model, as well as many others that were computed but not used in the final model',
                          path=os.path.relpath(cov_db_path, proj_path),
                      ),
                      Dataset(
                          xml_id='LOG',
                          name='Log File',
                          ds_type = 'File',
                          description='Processing log file',
                          path=os.path.relpath(log_path, proj_path),
                      ),
                  ]
              )]
          )
                          
        if no_st_values > 0:
            rs_project.realizations[0].datasets.append(
                Dataset(
                    xml_id='DAILY_STREAM_TEMPERATURE',
                    name='Daily stream temperature predictions',
                    ds_type = 'File',
                    description='Zipped database of the full set of daily stream temperature predictions (1990-2021) for each stream reach (i.e., each comid in the 1:100,000 National Hydrography Dataset).',
                    path=os.path.relpath(st_db_path, proj_path),
                ),)
        
        merged_project_xml = os.path.join(proj_path, 'project.rs.xml')
        rs_project.write(merged_project_xml)
        print(f'Project XML file written to {os.path.relpath(merged_project_xml, proj_path)}')
        #%% Finish
        runtime = datetime.datetime.now() - start
        print(f'\n\n\n----------  DONE. Runtime: {runtime.seconds // 60}:{runtime.seconds % 60:02d} minutes  ---------')
        sys.stdout = terminal_stdout
        print(f'\n\n\n----------  DONE. Runtime: {runtime.seconds // 60}:{runtime.seconds % 60:02d} minutes  ---------')
        stop
    try:
        res = subprocess.run(['rscli','upload','--org', 'fe204da2-8c52-4165-9e90-f4c9807ac57a','--visibility','PRIVATE', 
                          f'{proj_path}', '--verbose'], input = 'Y\r', capture_output=True, text = True)
        with open(os.path.join(outdir, 'upload_log.txt'), 'a') as f:
            if res.stdout[-14:-6]=='COMPLETE':
                f.write(f'\n@ {curhuc} - {h10name_hr} SUCCESS\n {res.stdout}')
            else: f.write(f'\n@ {curhuc} - {h10name_hr} ERROR\n {res.stdout}')
    
    except subprocess.CalledProcessError as e:
        print(e)
        with open(os.path.join(outdir, 'upload_log_retro.txt'), 'a') as f:
            f.write(f'\n@ {curhuc} - {h10name_hr} SUBPROCESS ERROR \n{e}')
