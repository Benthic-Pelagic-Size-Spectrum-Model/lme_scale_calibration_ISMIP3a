#!/usr/bin/env python3

#Loading libraries
import os
from glob import glob
import xarray as xr
import pandas as pd
import useful_functions as uf
import dask
from distributed import Client
from multiprocessing import Process, freeze_support

if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client()

    #Do not print warning about large chunks
    #dask.config.set({'array.slicing.split_large_chunks': False})
    #Reduce size of chunks to 100 MB
    #dask.config.set({'array.chunk-size': '100 MB'})

    ## Name of region and model resolution ----
    region = 'fao-58'
    model_res = '1deg'
    
    ## Defining input and output folders ----
    base_folder = '/g/data/vf71/la6889/dbpm_inputs/east_antarctica/'
    gridded_folder = os.path.join(base_folder, 'gridded_params', model_res)
    out_folder = os.path.join(base_folder, 'run_fishing', model_res)
    os.makedirs(out_folder, exist_ok = True) 

    ## If starting DBPM run from a specific time step ----
    # Character: Year and month from when DBPM initialisation values should be loaded
    # If starting model for the first time, it should be set to None
    init_time = '1960-10'
    
    ## Loading fixed DBPM parameters ----
    ds_fixed = uf.loading_dbpm_inputs(gridded_folder)
    #Adding additional fixed DBPM parameters to dataset
    #Depth
    depth = xr.open_zarr(glob(os.path.join(base_folder, 'gridded', 
                                           model_res, '*obsclim_deptho_*'))[0])['deptho']
    ds_fixed['depth'] = depth
    #Size bins in log10
    log10_size_bins_mat = xr.open_zarr('outputs/log10_size_bins_matrix.zarr/')['size_bins']
    ds_fixed['log10_size_bins_mat'] = log10_size_bins_mat
    #Removing datarrays added to fixed inputs
    del depth, log10_size_bins_mat

    # #Scatter fixed dataset across workers
    # ds_fixed_scattered = client.scatter(ds_fixed)
    # #Complete scattering before use ("dask future")
    # ds_fixed_fut = ds_fixed_scattered.result()

    
    ## Loading predator, detritivores and detritus initialisation data ----
    if init_time is None:
        predators = xr.open_zarr(glob(os.path.join(gridded_folder, 
                                                  'predators_*'))[0])['predators'] 
        detritivores = xr.open_zarr(glob(os.path.join(gridded_folder, 
                                                    'detritivores_*'))[0])['detritivores']
        detritus = xr.open_zarr(glob(os.path.join(gridded_folder, 
                                                  'detritus_*'))[0])['detritus']
    else:
        predators = xr.open_dataarray(
            glob(os.path.join(out_folder, f'predators_*_{init_time}.nc'))[0])
        detritivores = xr.open_dataarray(
            glob(os.path.join(out_folder, f'detritivores_*_{init_time}.nc'))[0])
        detritus = xr.open_dataarray(
            glob(os.path.join(out_folder, f'detritus_*_{init_time}.nc'))[0])

    #Create dataset for predator, detritivores and detritus initialisation data
    ds_init = xr.Dataset(data_vars = {'predators': predators, 
                                      'detritivores': detritivores, 
                                      'detritus': detritus})
    # #Scatter initialisation dataset across workers
    # ds_init_scattered = client.scatter(ds_init)
    # #Complete scattering before use ("dask future")
    # ds_init_fut = ds_init_scattered.result()

    ## Loading dynamic data ----
    #Spinup data is loaded if init_time is None or if the init_time year is less than 1960
    if init_time is None or pd.Timestamp(init_time).year < 1960:
        #Intercept plankton size spectrum
        ui0 = xr.open_zarr(glob(os.path.join(gridded_folder,
                                             'ui0_spinup*'))[0])['ui0']
        #Slope plankton size spectrum
        slope = xr.open_zarr(glob(
            os.path.join(base_folder, 'gridded', 
                         model_res, '*spinup_slope_*'))[0])['slope']
        #Temperature effects
        pel_tempeffect = xr.open_zarr(glob(
            os.path.join(gridded_folder, 
                         'pel-temp-eff_spinup*'))[0])['pel_temp_eff']
        ben_tempeffect = xr.open_zarr(glob(
            os.path.join(gridded_folder, 
                         'ben-temp-eff_spinup*'))[0])['ben_temp_eff']
        #Sinking rate
        sinking_rate = xr.open_zarr(glob(
            os.path.join(base_folder, 'gridded', model_res,
                         '*_spinup_er_*'))[0])['export_ratio']
        # Loading effort
        effort = xr.open_zarr(glob(os.path.join(gridded_folder, 
                                                'effort_spinup*'))[0])['effort']
    #Spinup data plus obsclim are loaded if init_time is 1960
    elif pd.Timestamp(init_time).year == 1960:
        exp = ['spinup', 'obsclim']
        #Intercept plankton size spectrum
        ui0 = xr.open_mfdataset(glob(os.path.join(gridded_folder, 'ui0_*')), 
                                engine = 'zarr')['ui0']
        #Slope plankton size spectrum
        slope = xr.open_mfdataset([f for ex in exp for f in glob(
            os.path.join(base_folder, 'gridded', model_res, f'*{ex}_slope_*'))], 
                                  engine = 'zarr')['slope']
        #Temperature effects
        pel_tempeffect = xr.open_mfdataset(glob(
            os.path.join(gridded_folder, 'pel-temp-eff_*')),
                                           engine = 'zarr')['pel_temp_eff']
        ben_tempeffect = xr.open_mfdataset(glob(
            os.path.join(gridded_folder, 'ben-temp-eff_*')),
                                           engine = 'zarr')['ben_temp_eff']
        #Sinking rate
        sinking_rate = xr.open_mfdataset([f for ex in exp for f in glob(
            os.path.join(base_folder, 'gridded', model_res, f'*{ex}_er_*'))], 
                                         engine = 'zarr')['export_ratio']
        # Loading effort
        effort = xr.open_mfdataset(glob(os.path.join(gridded_folder, 'effort_*')),
                                   engine = 'zarr')['effort']
    #Obsclim data loaded from 1961 onwards:
    else:
        #Intercept plankton size spectrum
        ui0 = xr.open_zarr(glob(
            os.path.join(gridded_folder, 'ui0_[0-9]*'))[0])['ui0']
        #Slope plankton size spectrum
        slope = xr.open_zarr(glob(
            os.path.join(base_folder, 'gridded', 
                         model_res, '*obsclim_slope_*'))[0])['slope']
        #Temperature effects
        pel_tempeffect = xr.open_zarr(glob(
            os.path.join(gridded_folder, 'pel-temp-eff_[0-9]*'))[0])['pel_temp_eff']
        ben_tempeffect = xr.open_zarr(glob(
            os.path.join(gridded_folder, 'ben-temp-eff_[0-9]*'))[0])['ben_temp_eff']
        #Sinking rate
        sinking_rate = xr.open_zarr(glob(
            os.path.join(base_folder, 'gridded', model_res,
                         '*_obsclim_er_*'))[0])['export_ratio']
        # Loading effort
        effort = xr.open_zarr(glob(
            os.path.join(gridded_folder, 'effort_[0-9]*'))[0])['effort']

    #If running from a specific point in time, then data is subsetted from the month after
    #init_time
    if init_time is not None:
        #Timestep from when to restart DBPM 
        subset_time = (pd.Timestamp(init_time)+
                       pd.DateOffset(months = 1)).strftime('%Y-%m')
        #Timestep from when to add init effort data
        effort_time = (pd.Timestamp(init_time)+
                       pd.DateOffset(months = 2)).strftime('%Y-%m')
        
        #begin (i.e., one month after predators data)
        ui0 = ui0.sel(time = slice(subset_time, None))
        slope = slope.sel(time = slice(subset_time, None))
        pel_tempeffect = pel_tempeffect.sel(time = slice(subset_time, None))
        ben_tempeffect = ben_tempeffect.sel(time = slice(subset_time, None))
        sinking_rate = sinking_rate.sel(time = slice(subset_time, None))
        #Effort   
        #Load effort for time step DBPM starts
        e_start = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                      f'effort_*_{subset_time}.nc'))[0])
        effort = effort.sel(time = slice(effort_time, None))
        #Combine both data arrays
        effort = xr.concat([e_start, effort], dim = 'time')

    #Create a single dataset for dynamic inputs
    ds_dynamic = xr.Dataset(data_vars = {'ui0': ui0, 'slope': slope,
                                         'pel_tempeffect': pel_tempeffect,
                                         'ben_tempeffect': ben_tempeffect, 
                                         'sinking_rate': sinking_rate,
                                         'effort': effort})
    # #Scatter initialisation dataset across workers
    # ds_dynamic_scattered = client.scatter(ds_dynamic)
    # #Complete scattering before use ("dask future")
    # ds_dynamic_fut = ds_dynamic_scattered.result()

    
    ## Running spatial DBPM ----
    uf.gridded_sizemodel(gridded_folder, ds_fixed, ds_init, ds_dynamic, 
                         region = region, model_res = model_res, out_folder = out_folder)

    #Closing client
    # client.close()
