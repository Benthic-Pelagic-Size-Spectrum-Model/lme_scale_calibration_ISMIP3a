#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
from glob import glob

# Defining base folder
base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_lme_inputs/'
# Getting list of FAO regions
fao_lme_code = [f for f in os.listdir(base_dir) if 'fao_lme' in f]

# Variables to be processed
vars_int = ['tob', 'tos', 'ocean-temp-weighted', 'er', 'simask', 'lphy', 
            'sphy', 'expc-bot', 'input-w20m']

# Experiments
exp_name = ['obsclim', 'ctrlclim', 'spinup', 'stable-spin']

# Apply data processing workflow to all regions - High resolution only
for aoi in fao_lme_code:
    # Defining input and output folders
    gridded_folder = os.path.join(base_dir, aoi, 'gridded', '025deg')
    gfdl_out = os.path.join(base_dir, aoi, 'monthly_weighted')
    os.makedirs(gfdl_out, exist_ok = True)

    # Load area of grid cells
    area = xr.open_zarr(glob(os.path.join(
        gridded_folder, '*area*'))[0])['cellareao'].fillna(0)

    # Crearting biomass-area weights
    # Calculating long-term (1961-1980) mean total phytoplankton 
    # biomass from "ctrlclim" experiment
    lphy = xr.open_zarr(glob(os.path.join(
        gridded_folder, '*ctrlclim_lphy*'))[0])['lphy']
    sphy = xr.open_zarr(glob(os.path.join(
        gridded_folder, '*ctrlclim_sphy*'))[0])['sphy']
    totphy = ((lphy+sphy).sel(time = slice('1961', '1980')).
        mean('time')).fillna(0)
    # Weighting by phytoplankton biomass per area of grid cell
    weights = totphy*area
    
    for exp in exp_name:
        try:
            depth = xr.open_zarr(glob(os.path.join(
                gridded_folder, f'*{exp}_deptho*'))[0])['deptho']
        except:
            depth = xr.open_zarr(glob(os.path.join(
                gridded_folder, '*ctrlclim_deptho*'))[0])['deptho']
        
        area_weighted_depth = (depth.weighted(area).
            mean(('lat', 'lon')).values)
        
        region_int = aoi.replace('-', ' ').upper()

        # Getting a list of all files contained in the LME/FAO folder
        all_fn = glob(os.path.join(gridded_folder, f'gfdl*{exp}*'))

        exp_fn = []
        for var in vars_int:
            [var_fn] = [fn for fn in all_fn if f'_{var}_' in fn]
            exp_fn.append(var_fn)

        weighted_inputs = uf.weighted_mean_timestep(
            exp_fn, weights, area, region_int)

        (weighted_inputs['intercept'], 
         weighted_inputs['slope']) = uf.GetPPIntSlope(
             sphy_file = weighted_inputs['sphy'].values, 
             lphy_file = weighted_inputs['lphy'].values)
        
        weighted_inputs['depth'] = area_weighted_depth
        weighted_inputs['depth_m_bio_weighted'] = (depth.
            weighted(weights).mean(('lat', 'lon')).values)
        
        #Saving data
        start_yr = weighted_inputs.year.min()
        end_yr = weighted_inputs.year.max()
        
        weighted_inputs.to_parquet(os.path.join(
            gfdl_out, 
            f'{exp}_dbpm_clim-inputs_{aoi}_{start_yr}-{end_yr}.parquet'),
                                   index = False)
