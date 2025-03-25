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
import sqlite3

os.path.dirname(__file__)
os.chdir(os.path.dirname(__file__))

def create_proj_dir(huc10):
    path = os.path.join(outdir, f'retro/{huc10[:6]}/{huc10}')
    if not os.path.exists(path):
        os.makedirs(path)
    return(path)

def add_cov_metadata():
    import shutil
    cov_metadata = os.path.join('DATA', 'covariate_metadata.csv')
    shutil.copy(cov_metadata, os.path.join(path, 'covariate_metadata.csv'))
    return

def add_database_file(file, date_col, db_type, add_cols = None, compression = None):
    df = pd.read_csv(file, parse_dates = [date_col], compression = compression)
    print(df.columns)
    dropcols = ['COMID', 'tim.date', 'cov.antec_air_temp', 'cov.std_mean_flow']
    outcols = df[dropcols]
    df.drop(dropcols[2:], axis = 1, inplace = True, errors = 'ignore')
    df.columns = df.columns.str.split('\.').str[-1].str.strip()
    print(df.columns)
    try: df_out = df[['COMID', 'date', 'stream_temp']]
    except: df_out = pd.merge(df, add_cols, how = 'left', on=['COMID','date'])
    out_table = df_out.fillna(-999).round(4)
    connection= sqlite3.connect(os.path.join(path, f'daily_{db_type}.db'))
    out_table.to_sql('daily_{db_type}', connection, if_exists = 'replace')
    connection.close()
    return(outcols)
    
    

tempdir = 'DATA/preds_retro'
outdir = 'DATA/Outputs/'

#%%
curhuc = '1701010107'

# Create project diorectory
path = create_proj_dir(curhuc) 

add_cov_metadata()

# Load in geometry files
                            
#%%
# Add daily temperature and covariate retrospective databases

#Temperature
temp_file = os.path.join(tempdir, 'predictions_temperature', f'st_pred_{curhuc}.csv')
outcols = add_database_file(temp_file, date_col = 'tim.date', db_type = 'stream_temperature', add_cols = None, compression = None)

#temp_out = temp_df[['COMID', 'date', 'stream_temp']]
#temp_out = temp_out.fillna(-999).round(4)

#connection= sqlite3.connect(os.path.join(path, 'daily_temperature_retrospective.db'))
#temp_out.to_sql('daily_stream_temperature', connection, if_exists = 'replace')
#connection.close()

#Covariates
cov_file = os.path.join(tempdir, 'predictions_covariates', 'cov_csvs', f'{curhuc}_covs.zip')
res = add_database_file(cov_file, date_col = 'date', db_type = 'covariates', 
                        add_cols = outcols.rename(columns={'tim.date':'date'}), compression = 'zip')
#cov_df = pd.read_csv(cov_file, compression = 'zip', parse_dates = ['date'])
#cov_df.drop(['tim.date','cov.antec_air_temp', 'cov.std_mean_flow'], axis = 1, inplace = True)

#cov_df_all = pd.merge(cov_df, temp_df.drop(['lookup', 'stream_temp'], axis =1), how = 'left', on=['COMID','date'])
#cov_df_all.columns = cov_df_all.columns.str.split('\.').str[-1].str.strip()
#cov_out = cov_df_all.fillna(-999).round(4)

#connection= sqlite3.connect(os.path.join(path, 'daily_covariate_retrospective.db'))
#cov_out.to_sql('daily_covariate_temperature', connection, if_exists = 'replace')
#connection.close()

#%%
# Build geopackage



