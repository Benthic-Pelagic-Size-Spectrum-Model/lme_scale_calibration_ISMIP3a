#!/usr/bin/env python3

# Library of useful functions
# Authors: Denisse Fierro Arcos
# Date last update: 2024-10-10

# Loading libraries
import xarray as xr
import re
import numpy as np
import os
from glob import glob


#Transforming netCDF files to zarr
def netcdf_to_zarr(file_path, path_out):
    '''
    Inputs:
    - file_path (character) File path where GFDL output is located
    - path_out (character) File path where outputs should be stored as zarr files

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    da = xr.open_dataarray(file_path).chunk({'lat': '50MB', 'lon': '50MB'})

    #Save results
    da.to_zarr(path_out, consolidated = True, mode = 'w')


## Extracting GFDL outputs for region of interest using boolean mask
def extract_gfdl(file_path, mask, path_out):
    '''
    Inputs:
    - file_path (character) File path where GFDL zarr file is located
    - mask (boolean data array) Grid cells within region of interest should be identified
    as 1.
    - path_out (character) File path where outputs should be stored as zarr files

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    da = xr.open_zarr(file_path)

    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    da = da[var]
    
    #Apply mask and remove rows where all grid cells are empty to reduce data array size
    da = da.where(mask == 1, drop = 'all')

    #Rechunking data
    if 'time' in da.dims:
        da = da.chunk({'time': '50MB', 'lat': '50MB', 'lon': '50MB'})
    else:
        da = da.chunk({'lat': '50MB', 'lon': '50MB'})

    #Save results
    da.to_zarr(path_out, consolidated = True, mode = 'w')


## Calculating area weighted means
def weighted_mean_timestep(file_paths, weights, region):
    '''
    Inputs:
    - file_paths (list) File paths pointing to zarr files from which weighted means 
    will be calculated and stored in a single data frame
    - weights (data array) Contains the weights to be used when calculating weighted mean.
    It should NOT include NaN, zeroes (0) should be used instead.
    - region (character) Name of the region to be recorded in data frame

    Outputs:
    - df (pandas data frame) containing area weighted mean per time step
    '''
    
    #Loading all zarr files into a single dataset
    da = xr.open_mfdataset(file_paths, engine = 'zarr')
    #Fix time format if needed
    try:
        new_time = da.indexes['time'].to_datetimeindex()
        da['time'] = new_time
    except:
        pass

    #Apply weights
    da_weighted = da.weighted(weights)
    #Calculate weighted mean
    da_weighted_mean = da_weighted.mean(('lat', 'lon'))
    
    #Transform to data frame
    df = da_weighted_mean.to_dataframe().reset_index()

    #Getting names of variables in data frame
    col_names = [i for i in df.columns if i != 'time']
    #Getting name of experiment from file path
    [exp] = re.findall('cobalt2_(.*?)_', file_paths[0])

    #If no depth file is included, add variable and leave it empty
    if 'deptho' not in col_names:
        df['depth_m'] = np.nan
    else:
        col_names.remove('deptho')
        df = df.rename(columns = {'deptho': 'depth_m'})
    
    #Add metadata to data frame
    df['tot_area_m2'] = weights.values.sum()
    df['year'] = df.apply(lambda x: x.time.year, axis = 1)
    df['month'] = df.apply(lambda x: x.time.strftime('%B'), axis = 1)
    df['region'] = region
    df['scenario'] = exp

    #Rearrange columns 
    names = ['region', 'scenario', 'time', 'year', 'month', 
             'depth_m', 'tot_area_m2'] + col_names
    df = df[names]

    return df


#Calculating export ratio
def getExportRatio(folder_gridded_data, gfdl_exp):
    '''
    Inputs:
    - folder_gridded_data (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - gfdl_exp (character) Select 'ctrl_clim' or 'obs_clim'

    Outputs:
    - sphy (data array) Contains small phytoplankton. This is data frame simply
    renames 'phypico-vint' data
    - lphy (data array) Contains large phytoplankton: difference between 'phyc_vint'
    and 'phypico_vint'
    - er (data array) Contains export ratio
    '''

    #Get list of files in experiment
    file_list = glob(os.path.join(folder_gridded_data, f'*_{gfdl_exp}_*'))
    
    #Load sea surface temperature
    tos = xr.open_zarr([f for f in file_list if '_tos_' in f][0])['tos']
    
    #load depth
    depth = xr.open_zarr([f for f in file_list if '_deptho_' in f][0])['deptho']
    
    #Load phypico-vint
    sphy = xr.open_zarr([f for f in file_list if '_phypico-vint_' in f][0])['phypico-vint']
    #Rename phypico-vint to sphy
    sphy.name = 'sphy'

    #Load phyc-vint
    ptotal = xr.open_zarr([f for f in file_list if '_phyc-vint_' in f][0])['phyc-vint']

    #Calculate large phytoplankton
    lphy = ptotal-sphy
    #Give correct name to large phytoplankton dataset
    lphy.name = 'lphy'

    #Calculate phytoplankton size ratios
    plarge = lphy/ptotal
    psmall = sphy/ptotal

    #Calculate export ration
    er = (np.exp(-0.032*tos)*((0.14*psmall)+(0.74*(plarge)))+
          (0.0228*(plarge)*(depth*0.004)))/(1+(depth*0.004))
    #If values are negative, assign a value of 0
    er = xr.where(er < 0, 0, er)
    #If values are above 1, assign a value of 1
    er = xr.where(er > 1, 1, er)
    er.name = 'export_ratio'
    
    return sphy, lphy, er


#Calculate slope and intercept
def GetPPIntSlope(file_paths, mmin = 10**(-14.25), mmid = 10**(-10.184), 
                  mmax = 10**(-5.25), output = 'both'):
    '''
    Inputs:
    - file_paths (list) File paths pointing to zarr files from which slope and
    intercept will be calculated
    - mmin (numeric)  Default is 10**(-14.25). ????
    - mmid (numeric)  Default is 10**(-10.184). ????
    - mmax (numeric)  Default is 10**(-5.25). ????
    - output (character) Default is 'both'. Select what outputs should be returned. 
    Choose from 'both', 'slope', or 'intercept'

    Outputs:
    - (Data array) - Depends on value of 'output' parameter. 
    '''
        
    #Load depth
    depth = xr.open_zarr([f for f in file_paths if '_deptho_' in f][0])['deptho']
    
    #load large phytoplankton
    lphy = xr.open_zarr([f for f in file_paths if '_lphy_' in f][0])['lphy']
    
    #Load small phytoplankton
    sphy = xr.open_zarr([f for f in file_paths if '_sphy_' in f][0])['sphy']
    
    #Convert sphy and lphy from mol C / m^3 to g C / m^3
    sphy = sphy*12.0107
    lphy = lphy*12.0107

    #From Appendix of Barnes 2010 JPR, used in Woodworth-Jefcoats et al 2013
    #GCB. The scaling of biomass with body mass can be described as B=aM^b
    #the exponent b (also equivalent to the slope b in a log B vs log M 
    #relationship) can be assumed:
    #0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
    # most studies seem to suggest N~M^-1, so can assume that and test 
    #sensitivity of our results to this assumption. 
      
    #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
    #in log10 gww
    midsmall = np.log10((mmin+mmid)/2) 
    midlarge = np.log10((mmid+mmax)/2)

    #convert to log10 (gww/size class median size) for log10 abundance
    small = np.log10((sphy*10)/(10**midsmall))
    #convert to log10 (gww/size class median size) for log10 abundance
    large = np.log10((lphy*10)/(10**midlarge))

    #Calculating lope
    b = (small-large)/(midsmall-midlarge)
    b.name = 'slope'

    #a is really log10(a), same a when small, midsmall are used
    a = large-(b*midlarge)
    a.name = 'intercept'

    # a could be used directly to replace 10^pp in sizemodel()
    if output == 'slope':
        return b
    if output == 'intercept':
        return a
    if output == 'both':
        return a, b


# Calculate spinup from gridded data
def gridded_spinup(file_path, start_spin, end_spin, spinup_period,
                  **kwargs):
    '''
    Inputs:
    - file_path (character) File path pointing to zarr file from which spinup
    will be calculated
    - start_spin (character or numeric) Year or date spinup starts
    - end_spin (character or numeric) Year or date spinup ends
    - spinup_period (pandas Datetime array) New time labels for spinup period.
    Must be a multiple of spinup range (i.e., difference between start and end
    spin)
    **Optional**: 
    - file_out (character) File path to save results

    Outputs:
    spinup_da (data array) Spinup data containing information within spinup
    range
    '''
    
    #Loading data
    da = xr.open_zarr(file_path)
    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    #Select period to be used for spinup
    da = da[var].sel(time = slice(str(start_spin), str(end_spin)))

    #Create spinup data array
    spinup_da = xr.concat([da, da], dim = 'time')
    while len(spinup_da.time) < len(spinup_period):
        spinup_da = xr.concat([spinup_da, da], dim = 'time')
    spinup_da['time'] = spinup_period
    
    #Updating chunks
    spinup_da = spinup_da.chunk({'time': '50MB', 'lat': 100, 'lon': 240})
    spinup_da = spinup_da.drop_encoding()

    #Save result if path is provided
    if kwargs.get('file_out', None):
        f_out = kwargs.get('file_out')
        #Ensure folder in file path exists
        os.makedirs(os.path.dirname(f_out), exist_ok = True)
        spinup_da.to_zarr(f_out, consolidated = True, mode = 'w')
    
    return spinup_da