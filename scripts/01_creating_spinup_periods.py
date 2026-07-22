#!/usr/bin/env python3

# Loading libraries
import os
from glob import glob
import numpy as np
import xarray as xr
import useful_functions as uf
import pandas as pd
import dask
from distributed import Client
from multiprocessing import Process, freeze_support

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    #Base folder where GFDL outputs are stored 
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a'

    #Define resolutions 
    resolutions = ['1deg', '025deg']

    #Define variables for which a spinup period will be created
    dynamic_vars = ['input-w20m', 'expc-bot', 'tob', 'tos', 'siconc']

    #Defining stable spin and spinup periods
    stable_spin = pd.date_range('1741-01', end = '1840-12', freq = 'MS')
    spinup_period = pd.date_range('1841-01', end = '1960-12', freq = 'MS')
   
    #Loop through experiments and resolutions
    for res in resolutions:
        #Define GFDL folder where sea ice masks are stored
        gfdl_folder = os.path.join(base_dir, 'global_gridded_zarr', res)

        [area_file] = glob(os.path.join(gfdl_folder, '*areacello*'))

        #The spinup period goes from 1841 and 1960. It is created by repeating 
        #inputs from "ctrlclim" experiment between 1961 and 1980
        #The stable spinup period goes from 1741 to 1840. It is created by 
        #repeating the mean for the year 1841 in the spinup period
        for dv in dynamic_vars:
            [f_in] = glob(os.path.join(gfdl_folder, 
                                       f'gfdl-*ctrlclim_{dv}_*monthly*'))
            f_out = (f_in.replace('ctrlclim', 'spinup').replace('1961', '1841').
                replace('2010', '1960'))
            uf.gridded_spinup(f_in, '1961-01', '1980-12', spinup_period, 
                              file_out = f_out)

            # The stable spinup period 
            fout_stable = (f_in.replace('ctrlclim', 'stable-spin').
                replace('1961', '1741').replace('2010', '1840'))
            uf.gridded_spinup(f_out, '1841-01', '1841-12', stable_spin,
                              mean_spinup = True, file_out = fout_stable)

        #Creating sea ice masks
        si_files = glob(os.path.join(gfdl_folder, '*siconc*'))
        for f in si_files:
            si_mask = uf.sea_ice_masks(f, area_file)
            mask_fn = f.replace('_siconc_', '_simask_')
            si_mask.to_zarr(mask_fn, consolidated = True, mode = 'w')

