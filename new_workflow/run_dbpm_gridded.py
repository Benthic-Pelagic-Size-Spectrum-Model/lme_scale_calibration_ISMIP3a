#!/usr/bin/env python3

#Loading libraries
import os
from glob import glob
import xarray as xr
import useful_functions as uf
from dask.distributed import Client
from multiprocessing import Process, freeze_support

if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1)

    #Name of region and model resolution
    region = 'fao-48'
    model_res = '1deg'
    
    #Defining input and output folders
    base_folder = '/g/data/vf71/la6889/dbpm_inputs/weddell'
    gridded_folder = os.path.join(base_folder, 'gridded_params', model_res)
    out_folder = os.path.join(base_folder, 'run_fishing', model_res)
    os.makedirs(out_folder, exist_ok = True) 
    
    #Loading fixed DBPM parameters
    depth = xr.open_zarr(
         glob(os.path.join(base_folder, 'gridded', model_res, '*obsclim_deptho_*'))[0])['deptho']
    log10_size_bins_mat = xr.open_zarr('outputs/log10_size_bins_matrix.zarr/')['size_bins']
    
    ## Loading initial predator data
    predators = xr.open_zarr(glob(os.path.join(gridded_folder, 
                                               'predators_spinup*1845*'))[0])['predators'] 
    
    #Load predator results for the month before you want to run DBPM
    p = xr.open_dataarray(glob(os.path.join(out_folder, 'predators_*_1844-12.nc'))[0])
    # #Subset predator data for the period you want to run the model
    # predators = predators.sel(time = slice('1845-01', '1849-12'))
    #Combine both data arrays
    predators = xr.concat([p, predators], dim = 'time')
    
    ## Load detritus and detritivores
    #Initialisiation files
    # detritivores = xr.open_zarr(glob(os.path.join(gridded_folder, 
    #                                               'detritivores_*'))[0])['detritivores']
    # detritus = xr.open_zarr(glob(os.path.join(gridded_folder, 'detritus_*'))[0])['detritus']

    #If picking up from a specific time step
    detritivores = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                       'detritivores_*_1844-12.nc'))[0])
    detritus = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                   'detritus_*_1844-12.nc'))[0])
    
    ## Loading effort
    effort = xr.open_zarr(glob(os.path.join(gridded_folder, 'effort_spinup*'))[0])['effort']
    
    #Load effort for time step DBPM starts
    e_start = xr.open_dataarray(glob(os.path.join(out_folder, 'effort_*_1845-01.nc'))[0])
    #Subset effort data from the month after DBPM begins
    effort = effort.sel(time = slice('1845-02', '1849-12'))
    #Combine both data arrays
    effort = xr.concat([e_start, effort], dim = 'time')
    
    #Temperature effects
    pel_tempeffect = xr.open_zarr(glob(
        os.path.join(gridded_folder, 'pel-temp-eff_spinup*'))[0])['pel_temp_eff']
    
    ben_tempeffect = xr.open_zarr(glob(
        os.path.join(gridded_folder, 'ben-temp-eff_spinup*'))[0])['ben_temp_eff']
    
    pel_tempeffect = pel_tempeffect.sel(time = slice('1845-01', '1849-12'))
    ben_tempeffect = ben_tempeffect.sel(time = slice('1845-01', '1849-12'))
    
    #Sinking rate
    sinking_rate = xr.open_zarr(
        glob(os.path.join(base_folder, 'gridded', model_res, '*_spinup_er_*'))[0])['export_ratio']
    
    sinking_rate = sinking_rate.sel(time = slice('1845-01', '1849-12'))
    
    #Running DBPM
    uf.gridded_sizemodel(gridded_folder, predators, detritivores, detritus, pel_tempeffect, 
                         ben_tempeffect, effort, sinking_rate, depth, log10_size_bins_mat,
                         region = region, model_res = model_res, out_folder = out_folder)