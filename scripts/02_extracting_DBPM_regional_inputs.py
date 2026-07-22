#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
import numpy as np
import pandas as pd
from glob import glob
from distributed import Client
from multiprocessing import Process, freeze_support

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    # Defining folder where masks with all FAO regions are stored
    mask_folder = '/g/data/vf71/shared_resources/fao_lme_masks'
    
    # Define location of DBPM inputs
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/'
    
    # Define variables for which data will be extracted
    vars_int = ['input-w20m', 'expc-bot', 'simask', 'tob', 'tos', 'deptho', 
                'areacello', 'phyc', 'phypico', 'thetao', 'bio-weights',
                'thkcello']
    
    # Define resolutions
    resolutions = ['1deg', '025deg']
    
    #Defining stable spin and spinup periods
    stable_spin = pd.date_range('1741-01', end = '1840-12', freq = 'MS')
    spinup_period = pd.date_range('1841-01', end = '1960-12', freq = 'MS')
    
    for res in resolutions:
        if res == '1deg':
            mask_all = xr.open_dataarray(os.path.join(
                mask_folder, 
                'gfdl-mom6-cobalt2_fao-major_lme_60arcmin_global_fixed.nc'))
        elif res == '025deg':
            mask_all = xr.open_dataarray(os.path.join(
                mask_folder, 
                'gfdl-mom6-cobalt2_fao-major_lme_15arcmin_global_fixed.nc'))
        
        # Getting region codes included in mask
        fao_lme_id = (np.unique(mask_all.values[np.isfinite(mask_all.values)]).
            astype(int))
        
        #Define GFDL folder
        file_list = glob(os.path.join(
            base_dir, 'global_gridded_zarr', res, '*'))

        #List all files to be extracted
        for aoi in fao_lme_id:
            gfdl_out = os.path.join(
                base_dir, 'fao_lme_inputs', f'fao_lme-{aoi}', 'gridded', res)
            os.makedirs(gfdl_out, exist_ok = True)
            mask = xr.where(mask_all == aoi, 1, np.nan)
            for dv in vars_int:
                #Ignoring files not needed
                file_dv = [f for f in file_list if f'_{dv}_' in f]
                
                #Extracting data for FAO area
                for f in file_dv:
                    #Create file path to save outputs
                    f_out = os.path.basename(f).replace('global', 
                                                        f'fao_lme-{aoi}')
                    f_out = os.path.join(gfdl_out, f_out)
                    #Apply function
                    if aoi in [1, 54, 65, 161, 171, 181, 188]:
                        cross_dateline = True
                    else:
                        cross_dateline = False
                    uf.extract_gfdl(f, mask, f_out, 
                                    cross_dateline = cross_dateline)

            # Vertically integrate phytoplankton inputs up to threshold depth
            [weight_file] = glob(os.path.join(gfdl_out, '*_bio-weights_*'))
            for exp in ['ctrlclim', 'obsclim']:
                base_fn = (weight_file.replace('_bio-weights_', '_var_').
                    replace('_ctrlclim_', f'_{exp}_'))
                phyc, phypico, temp_ocean = uf.integrating_inputs(
                    gfdl_out, exp, thresh_depth = 200, averaging = 'custom', 
                    weights = weight_file)
                #Save outputs
                phyc.to_zarr(base_fn.replace('_var_', '_phyc-vint-weighted_'),
                             consolidated = True, mode = 'w')
                phypico.to_zarr(base_fn.replace(
                    '_var_', '_phypico-vint-weighted_'), 
                                consolidated = True, mode = 'w')
                temp_ocean.to_zarr(base_fn.replace(
                    '_var_', '_ocean-temp-weighted_'), 
                                   consolidated = True, mode = 'w')

                #Calculate phytoplankton size distribution and export ratio
                sphy, lphy, er = uf.getExportRatio(gfdl_out, exp)
                #Save outputs
                sphy.to_zarr(base_fn.replace('_var_', '_sphy_'), 
                             consolidated = True, mode = 'w')
                lphy.to_zarr(base_fn.replace('_var_', '_lphy_'),
                             consolidated = True, mode = 'w')
                er.to_zarr(base_fn.replace('_var_', '_er_'),
                           consolidated = True, mode = 'w')

            #The spinup period goes from 1841 and 1960. It is created by repeating 
            #inputs from "ctrlclim" experiment between 1961 and 1980
            #The stable spinup period goes from 1741 to 1840. It is created by 
            #repeating the mean for the year 1841 in the spinup period
            for dv in ['sphy', 'lphy', 'er', 'ocean-temp-weighted']:
                [f_in] = glob(os.path.join(
                    gfdl_out, f'gfdl-*ctrlclim_{dv}_*monthly*'))
                f_out = (f_in.replace('ctrlclim', 'spinup').replace('1961', '1841').
                    replace('2010', '1960'))
                uf.gridded_spinup(f_in, '1961-01', '1980-12', spinup_period, 
                                  file_out = f_out)
    
                # The stable spinup period 
                fout_stable = (f_in.replace('ctrlclim', 'stable-spin').
                    replace('1961', '1741').replace('2010', '1840'))
                uf.gridded_spinup(f_out, '1841-01', '1841-12', stable_spin,
                                  mean_spinup = True, file_out = fout_stable)
            
            # Calculating slope and intercept
            lphy_files = sorted(glob(os.path.join(gfdl_out, '*_lphy_*')))
            sphy_files = sorted(glob(os.path.join(gfdl_out, '*_sphy_*')))
            for l, s in zip(lphy_files, sphy_files):
                intercept, slope = uf.GetPPIntSlope(sphy_file = s, lphy_file = l)
                #Save outputs
                intercept.to_zarr(l.replace('_lphy_', '_intercept_'), 
                                  consolidated = True, mode = 'w')
                slope.to_zarr(l.replace('_lphy_', '_slope_'), 
                              consolidated = True, mode = 'w')
