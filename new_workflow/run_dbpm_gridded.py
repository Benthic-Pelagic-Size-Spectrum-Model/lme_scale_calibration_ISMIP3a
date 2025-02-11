#!/usr/bin/env python3

#Loading libraries
import os
from glob import glob
import xarray as xr
import pandas as pd
import useful_functions as uf
from dask.distributed import Client
from multiprocessing import Process, freeze_support

if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1)

    #Name of region and model resolution
    region = 'fao-48'
    model_res = '025deg'
    
    #Defining input and output folders
    base_folder = '/g/data/vf71/la6889/dbpm_inputs/weddell/'
    gridded_folder = os.path.join(base_folder, 'gridded_params', model_res)
    out_folder = os.path.join(base_folder, 'run_fishing', model_res)
    os.makedirs(out_folder, exist_ok = True) 
    
    #Loading fixed DBPM parameters
    depth = xr.open_zarr(glob(os.path.join(base_folder, 'gridded', 
                                           model_res, '*obsclim_deptho_*'))[0])['deptho']
    log10_size_bins_mat = xr.open_zarr('outputs/log10_size_bins_matrix.zarr/')['size_bins']
    
    ## Loading predator, detritivores and detritus initialisation data
    #predators = xr.open_zarr(glob(os.path.join(gridded_folder, 
    #                                           'predators_*'))[0])['predators'] 
    #detritivores = xr.open_zarr(glob(os.path.join(gridded_folder, 
    #                                             'detritivores_*'))[0])['detritivores']
    #detritus = xr.open_zarr(glob(os.path.join(gridded_folder, 'detritus_*'))[0])['detritus']
    
    # If picking up from a specific time step:
    # Load predator results for the month before you want to run DBPM
    init_time = '1851-03'
    #Timestep from when to restart DBPM 
    subset_time = (pd.Timestamp(init_time)+pd.DateOffset(months = 1)).strftime('%Y-%m')
    #Timestep from when to add init effort data
    effort_time = (pd.Timestamp(init_time)+pd.DateOffset(months = 2)).strftime('%Y-%m')
    predators = xr.open_dataarray(glob(os.path.join(out_folder, f'predators_*_{init_time}.nc'))[0])
    detritivores = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                       f'detritivores_*_{init_time}.nc'))[0])
    detritus = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                   f'detritus_*_{init_time}.nc'))[0])

    # Dynamic data
    #Intercept plankton size spectrum
    ui0 = xr.open_zarr(glob(os.path.join(gridded_folder, 'ui0_spinup*'))[0])['ui0']
    
    #Slope plankton size spectrum
    slope = xr.open_zarr(glob(os.path.join(base_folder, 'gridded', 
                                           model_res, '*spinup_slope_*'))[0])['slope']
    #Temperature effects
    pel_tempeffect = xr.open_zarr(glob(
        os.path.join(gridded_folder, 'pel-temp-eff_spinup*'))[0])['pel_temp_eff']
    
    ben_tempeffect = xr.open_zarr(glob(
        os.path.join(gridded_folder, 'ben-temp-eff_spinup*'))[0])['ben_temp_eff']

    #Sinking rate
    sinking_rate = xr.open_zarr(
        glob(os.path.join(base_folder, 'gridded', model_res, '*_spinup_er_*'))[0])['export_ratio']

    # Loading effort
    effort = xr.open_zarr(glob(os.path.join(gridded_folder, 'effort_spinup*'))[0])['effort']

    #If running from a specific point in time, then subset data from the month DBPM is meant to 
    #begin (i.e., one month after predators data)
    ui0 = ui0.sel(time = slice(subset_time, None))
    slope = slope.sel(time = slice(subset_time, None))
    pel_tempeffect = pel_tempeffect.sel(time = slice(subset_time, None))
    ben_tempeffect = ben_tempeffect.sel(time = slice(subset_time, None))
    sinking_rate = sinking_rate.sel(time = slice(subset_time, None))
    #Effort   
    #Load effort for time step DBPM starts
    e_start = xr.open_dataarray(glob(os.path.join(out_folder, f'effort_*_{subset_time}.nc'))[0])
    effort = effort.sel(time = slice(effort_time, None))
    #Combine both data arrays
    effort = xr.concat([e_start, effort], dim = 'time')
    
    #Running DBPM
    uf.gridded_sizemodel(gridded_folder, predators, detritivores, detritus, pel_tempeffect, 
                         ben_tempeffect, effort, ui0, sinking_rate, slope, depth, 
                         log10_size_bins_mat, region = region, model_res = model_res, 
                         out_folder = out_folder)
