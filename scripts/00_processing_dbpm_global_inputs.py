#!/usr/bin/env python3

# Loading libraries
import os
import numpy as np
from glob import glob
import useful_functions as uf
import xarray as xr
import dask
from distributed import Client
from multiprocessing import Process, freeze_support

# This script pre-processes outputs from GFDL-MOM6-COBALT2 (referred to as GFDL 
# from hereon) prior to their use as forcings in DBPM. GFDL outputs are  
# available at two horizontal resolutions (1 deg and 0.25 deg).
# Here, GFDL outputs are stored as zarr files to process data faster and 
# phytoplankton intercept and slope, as well as export ratio are calculated.

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    #Base folder where GFDL outputs are stored 
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/global_inputs'

    #Define location of area of grid cell
    grid_dir = '/g/data/vf71/shared_resources/grid_cell_vars_ESMs/isimip3a'
    
    #Define experiments and resolution
    exp_name = ['ctrlclim', 'obsclim']
    resolutions = ['1deg', '025deg']

    #Define variables of interest
    dbpm_var = ['phyc', 'phypico', 'siconc', 'deptho', 'expc-bot', 'tob', 'tos',
                'thetao', 'mlotst-0125']
    
    #Loop through experiments and resolutions
    for res in resolutions:
        #Define output folder
        gfdl_out = f'/g/data/vf71/fishmip_inputs/ISIMIP3a/global_gridded_zarr/{res}'
        os.makedirs(gfdl_out, exist_ok = True)
        
        #Store degrees to arcmin
        if res == '1deg':
            arc_res = '60arcmin'
        elif res == '025deg':
            arc_res = '15arcmin'

        # Process area of grid cell files
        [area_file] = glob(os.path.join(grid_dir, 
                                        f'gfdl*_areacello_{arc_res}*.nc'))
        f_out = os.path.basename(area_file).replace('.nc', '.zarr')
        f_out = os.path.join(gfdl_out, f_out)
        uf.netcdf_to_zarr(area_file, f_out)
        #Load file area to use as land mask for sea ice processing step
        area = xr.open_zarr(f_out)['cellareao']

        # Process thickness of grid cell files
        [thick_file] = glob(os.path.join(grid_dir,
                                         f'gfdl*_thkcello_{arc_res}*.nc'))
        f_out = os.path.basename(thick_file).replace('.nc', '.zarr')
        f_out = os.path.join(gfdl_out, f_out)
        uf.netcdf_to_zarr(thick_file, f_out)

        #Process all other variables
        for exp in exp_name:
            #Define folder containing netCDF files
            gfdl_folder = os.path.join(base_dir, exp, res)
            base_fn = f'gfdl-mom6-cobalt2_{exp}_var_{arc_res}_global_monthly_1961_2010.zarr'
            for var in dbpm_var:
                [gfdl_file] = glob(os.path.join(gfdl_folder, f'*clim_{var}_*'))
                f_out = os.path.basename(gfdl_file).replace('.nc', '.zarr')
                f_out = os.path.join(gfdl_out, f_out)
                #Apply function
                uf.netcdf_to_zarr(gfdl_file, f_out)

                if var == 'phyc' and exp == 'ctrlclim':
                    depths = (xr.open_dataarray(thick_file).drop_vars('time').
                        squeeze().fillna(0))
                    phyc = (xr.open_dataarray(f_out).
                        sel(time = slice('1961', '1980')).mean('time').fillna(0))
                    weights = phyc*depths
                    weights.name = 'weights'
                    weights = weights.drop_attrs()
                    weights = weights.assign_attrs({
                        'long_name':
                        'Climatological mean (1961-1980) biomass weighting'})
                    fn_weights = f_out.replace('_phyc_', '_bio-weights_')
                    weights.drop_encoding().to_zarr(
                        fn_weights, consolidated = True, mode = 'w')
                
            # Transforming expc-bot (mol m-2 s-1) to input_w (gWW m-3 yr-1)
            input_w = uf.detrital_input_seafloor(gfdl_out, exp,
                                                 benthic_habitat_depth = 20)
            # Save outputs
            input_w.drop_encoding().to_zarr(
                os.path.join(gfdl_out, 
                             base_fn.replace('_var_', '_input-w20m_')), 
                consolidated = True, mode = 'w')

            # Vertically integrate phytoplankton inputs up to threshold depth
            [weight_file] = glob(os.path.join(gfdl_out, '*_bio-weights_*'))
            phyc, phypico, temp_ocean = uf.integrating_inputs(
                gfdl_out, exp, thresh_depth = 200, 
                averaging = 'custom', weights = weight_file)
            #Save outputs
            phyc.to_zarr(
                os.path.join(gfdl_out, base_fn.replace(
                    '_var_', '_phyc-vint-weighted_')), 
                consolidated = True, mode = 'w')
            phypico.to_zarr(
                os.path.join(gfdl_out, base_fn.replace(
                    '_var_', '_phypico-vint-weighted_')), 
                consolidated = True, mode = 'w')
            temp_ocean.to_zarr(
                os.path.join(gfdl_out, base_fn.replace(
                    '_var_', '_ocean-temp-weighted_')), 
                consolidated = True, mode = 'w')

            #Calculate phytoplankton size distribution and export ratio
            sphy, lphy, er = uf.getExportRatio(gfdl_out, exp)
            #Save outputs
            sphy.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_sphy_')), 
                consolidated = True, mode = 'w')
            lphy.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_lphy_')),
                consolidated = True, mode = 'w')
            er.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_er_')),
                consolidated = True, mode = 'w')
