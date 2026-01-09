#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 10:40:55 2025

@author: dawn.urycki
"""

import pandas as pd
import json
import glob
import sys
import subprocess
import os
import shutil
import time

projects_dir = 'DATA/Outputs/gcms/'
gpkgdir = os.path.join(projects_dir, 'geopackages')

for i in range(6):
    print(i)
    
    old_proj = pd.read_table('DATA/Outputs/gcms/update_readme.txt', header= None)
    uploaded.columns = ['huc10', 'hname', 'id']
    
    #built.to_json(f'{projects_dir}upload_log.json', orient = 'records', indent = 4)
    built = pd.concat([built_clean.set_index('huc10'), uploaded.set_index('huc10')['id']], axis = 1, join = 'outer').reset_index()
    #log = pd.read_json(f'{projects_dir}upload_log.json') 
        
    for i, row in built.iterrows():
        print (i, row.huc10, row.status, row.id)
        if (row.status == 'COMPLETE') and (len(str(row.id)) != 36):
                paths = glob.glob(f'{projects_dir}{str(row.huc10)[:6]}/Stream_Temperature_Projections_1950-2099_{row.hname.replace(' ', '_')}')
                print(f'uploading {paths}')
                try:
                    path = paths[0]
                    res = subprocess.run(['rscli','upload','--org', 'fe204da2-8c52-4165-9e90-f4c9807ac57a','--visibility','PUBLIC', 
                                      f'{path}', '--verbose'], input = 'Y\r', capture_output=True, text = True)
                    ware_id = res.stdout.split('Project ID: ', 1)[1][:36]
                    with open(os.path.join(projects_dir, 'upload_log_proj.txt'), 'a') as f:
                        if res.stdout[-14:-6]=='COMPLETE':
                            f.write(f'\n@ {row.huc10} - {row.hname} SUCCESS\n {res.stdout}')
                            with open(f'{projects_dir}upload_list.txt', 'a') as j:
                                j.write(f'\n{row.huc10}\t{row.hname}\t{ware_id}')
                        
                            outdir = os.path.join(gpkgdir, f'{str(row.huc10)}_anomalies.gpkg')
                            try:
                                shutil.copy2(os.path.join(path, 'seasonal_anomalies_spatial_covariates.gpkg'), outdir)
                                try: 
                                    shutil.copy2(os.path.join(path, 'processing.log'), os.path.join(gpkgdir, 'logs', f'{str(row.huc10)}_processing.log'))
                                    try:
                                        shutil.rmtree(path)
                                        print(f"Directory '{path}' and all its contents have been removed.")
                                    except OSError as e:
                                        print(f"Error: {path} : {e.strerror}")
                                except Exception as e:
                                    print(f"A processing.log error occurred: {e}")
                                    continue
                            except FileNotFoundError:
                                print("Source file geopackage not found.")
                                continue
                            except Exception as e:
                                print(f"An geopackage error occurred: {e}")
                              
                            
                        else: f.write(f'\n@ {row.huc10} - {row.hname} ERROR\n {res.stdout}')   
                except subprocess.CalledProcessError as e:
                     print(e)
                     with open(os.path.join(projects_dir, 'upload_log_proj.txt'), 'a') as f:
                         f.write(f'\n@ {row.huc10} - {row.hname} SUBPROCESS ERROR \n{e}')
                     continue 
                except Exception as e:
                    with open(os.path.join(projects_dir, 'upload_log_proj.txt'), 'a') as f:
                    	f.write(f'\n@ {row.huc10} - {row.hname} ERROR: \n{e}')
                    print(f"An error occurred: {e}")
                    continue   
        else: pass
    print('sleeping')    
    time.sleep(60*60*4)

#%%
# xml_files = glob.glob(projects_dir+'/*/*/project.rs.xml')

# import xml.etree.ElementTree as ET
# for file_path in xml_files:
#     try:
#         tree = ET.parse(file_path)
#         root = tree.getroot()

#         # Example: Fix a specific typo in text content
#         # Iterate through all elements to find text
#         for elem in root.iter():
#             if elem.text:
#                 # Simple replace for exact typo
#                 if "/.zip" in elem.text:
#                     elem.text = elem.text.replace("/.zip", ".zip")
#         tree.write(file_path, encoding='utf-8', xml_declaration=True)
#         print(f"Fixed: {file_path}")
        
#     except ET.ParseError as e:
#         print(f"Error parsing {file_path}: {e}")
#     except Exception as e:
#         print(f"An unexpected error occurred with {file_path}: {e}")  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        