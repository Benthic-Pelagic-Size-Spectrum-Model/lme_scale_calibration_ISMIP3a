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
import pandas as pd


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
                   mean_spinup = False, **kwargs):
    '''
    Inputs:
    - file_path (character) File path pointing to zarr file from which spinup
    will be calculated
    - start_spin (character or numeric) Year or date spinup starts
    - end_spin (character or numeric) Year or date spinup ends
    - spinup_period (pandas Datetime array) New time labels for spinup period.
    Must be a multiple of spinup range (i.e., difference between start and end
    spin)
    - mean_spinup (boolean) Default is False. If set to True, then the spinup
    period will be based on the mean value over the spin period.
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
    #If spinup should be created based on mean values over spinup period
    if mean_spinup:
        da = da.mean('time')

    #Create spinup data array
    spinup_da = [da] * len(spinup_period)
    spinup_da = xr.concat(spinup_da, dim = 'time')
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

### DPBM functions ----
# Build a lookup table for diet preference. Looks at all combinations of predator 
# and prey body size: diet preference (in the predator spectrum only)
def phi_f(q, log10_pred_prey_ratio, log_prey_pref):
    phi = np.where(q > 0, 
                   np.exp(-(q-log10_pred_prey_ratio)*(q-log10_pred_prey_ratio)/
                          (2*log_prey_pref**2))/(log_prey_pref*np.sqrt(2.0*np.pi)),
                   0) 
    return phi


# Function to build lookup tables for (constant) growth
# Considers components which remain constant
def gphi_f(pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref):
    gphi = 10**(-pred_prey_matrix)*phi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                                         log_prey_pref)
    return gphi


# Function to build lookup tables for (constant) mortality
def mphi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref, 
           metabolic_req_pred):
    mphi = (10**(metabolic_req_pred*rev_pred_prey_matrix)*
            phi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref))
    return mphi


# Function to build lookup table for components of 10^(alpha*x)
def expax_f(log10_size_bins, metabolic_req_pred):
    expax = list(10**(np.array(log10_size_bins)*metabolic_req_pred))
    return expax


# Create predator-prey combination matrix
# The square matrix holds the log(predatorsize/preysize) for all combinations of
# predator-prey sizes
def pred_prey_matrix(log10_size_bins):
    '''
    Inputs:
    - log10_size_bins (numpy array) Containing the log of all sizes of predator and
    prey items
    
    Outputs:
    ppm (numpy array) Contains all combinations of predator-prey sizes
    '''
    
    ppm = np.empty([len(log10_size_bins), len(log10_size_bins)])
    for i, j in enumerate(log10_size_bins):
        ppm[:, i] = log10_size_bins[i] - log10_size_bins
    return ppm


# Initialising matrices using size classes and time as dimensions
def init_da(log10_size_bins, time):
    '''
    Inputs:
    - log10_size_bins (numpy array) Containing the log of all sizes of predator and
    prey items
    - time (numpy array) Containing dates to be included in final data array
    
    Outputs:
    da (data array) Containing zeroes with size classes and time as dimensions
    '''
    data_start = np.zeros((len(log10_size_bins), len(time)))
    da = xr.DataArray(data = data_start, dims = ['size_class', 'time'], 
                      coords = {'size_class': log10_size_bins, 'time': time})
    return da


# Run model per grid cell or averaged over an area ------
def sizemodel(params, dbpm_input, ERSEM_det_input = False, temp_effect = True,
              use_init = False):
    '''
    Inputs:
    - params (dictionary) Output from 'sizemodel' function
    - dbpm_input (data frame) xxx
    - ERSEM_det_input (boolean) Default is False. xxxx
    - temp_effect (boolean) Default is True. xxxx
    - use_init (boolean) Default is False. xxxx

    Outputs:
    - return_list (dictionary) xxxx
    '''
    # time series of intercept of plankton size spectrum (estimated from GCM, 
    # biogeophysical model output or satellite data).		
    ui0 = xr.DataArray(data = 10**np.array(params['int_phy_zoo']),
                       dims = 'time', coords = {'time': dbpm_input.time})
    slope_phy_zoo = xr.DataArray(data = params['slope_phy_zoo'], dims = 'time',
                                 coords = {'time': dbpm_input.time})
    effort = xr.DataArray(data = params['effort'], dims = 'time', 
                          coords = {'time': dbpm_input.time})
    [numb_size_bins] = params['numb_size_bins']
    #Size bin index
    size_bin_index = np.arange(0, numb_size_bins)
    #Size bins
    log10_size_bins = np.array(params['log10_size_bins'])
    size_bin_vals = 10**log10_size_bins
    [log10_pred_prey_ratio] = params['log10_pred_prey_ratio']
    [log_prey_pref] = params['log_prey_pref']
    [metabolic_req_pred] = params['metabolic_req_pred']
    [metabolic_req_detritivore] = params['metabolic_req_detritivore']
    #Index for minimum detritivore size
    [ind_min_detritivore_size] = params['ind_min_detritivore_size']
    #Index for minimum predator size
    [ind_min_pred_size] = params['ind_min_pred_size']
    #Index for minimum size of detritivore fished
    [ind_min_fish_det] = params['ind_min_fish_det']
    #Index for minimum size of predator fished
    [ind_min_fish_pred] = params['ind_min_fish_pred']
    #Time
    time = pd.date_range(dbpm_input.time.min()-pd.DateOffset(months = 1),
                         dbpm_input.time.max(), freq = 'MS')
    #Sinking rate
    sinking_rate = xr.DataArray(data = params['sinking_rate'], dims = 'time', 
                                coords = {'time': dbpm_input.time})

    #data arrays for recording the two size spectra (V - det, U - pred)
    predators = init_da(log10_size_bins, time)
    detritivores = init_da(log10_size_bins, time)
    
    #vector to hold detritus biomass density (g.m-2)
    detritus = xr.DataArray(data = np.zeros(len(time)),
                            dims = 'time',
                            coords = {'time': time})
    
    #data arrays for keeping track of growth (GG_v, GG_u) and reproduction (R_v, R_u) 
    #from ingested food:
    growth_int_pred = growth_det = reprod_pred = reprod_det = init_da(log10_size_bins, time)
    # reprod_pred = init_da(log10_size_bins, time)
    # growth_det = init_da(log10_size_bins, time)
    # growth_int_pred = init_da(log10_size_bins, time)
    
    #data arrays for keeping track of predation mortality (PM_v, PM_u)
    pred_mort_pred = pred_mort_det = init_da(log10_size_bins, time)
    # pred_mort_pred = init_da(log10_size_bins, time)
    
    #data arrays for keeping track of catches (Y_v, Y_u)
    catch_pred = catch_det = init_da(log10_size_bins, time)
    # catch_pred = init_da(log10_size_bins, time)
    
    #data arrays for keeping track of total mortality (Z_v, Z_u)
    tot_mort_pred = tot_mort_det = init_da(log10_size_bins, time)
    # tot_mort_pred = init_da(log10_size_bins, time)
    
    #data arrays to hold fishing mortality rates at each size class at time 
    #(Fvec_v, Fvec_u)
    fishing_mort_pred = fishing_mort_det = init_da(log10_size_bins, time)
    # fishing_mort_pred = init_da(log10_size_bins, time)
    
    #data arrays for keeping track of senescence mortality (SM_v, SM_u) and other 
    #(intrinsic) mortality (OM_v, OM_u)
    other_mort_det = senes_mort_pred = senes_mort_det = np.zeros(len(log10_size_bins))
    # senes_mort_pred = np.zeros(len(log10_size_bins))
    # other_mort_det = np.zeros(len(log10_size_bins))
    other_mort_pred = np.zeros(len(log10_size_bins))

    #Build look up tables
    #lookup tables for terms in the integrals which remain constant over time (gphi, mphi)
    constant_growth = gphi_f(pred_prey_matrix(log10_size_bins), log10_pred_prey_ratio, 
                             log_prey_pref)
    constant_mortality = mphi_f(-pred_prey_matrix(log10_size_bins), log10_pred_prey_ratio, 
                                log_prey_pref, metabolic_req_pred)
    
    #lookup table for components of 10^(metabolic_req_pred*log10_size_bins) (expax)
    met_req_log10_size_bins = expax_f(log10_size_bins, metabolic_req_pred)

    #Numerical integration
    # set up with the initial values from param - same for all grid cells
    #(phyto+zoo)plankton size spectrum (U) 
    predators.isel(size_class = slice(None, 120)).loc[{'time': predators.time.min()}] = \
    params['plank_pred_sizes'][:120]
    
    # set initial detritivore spectrum (V)
    detritivores.isel(size_class = slice((ind_min_detritivore_size-1), 120)).\
    loc[{'time': detritivores.time.min()}] = \
    np.array(params['detritivore_sizes'])[ind_min_detritivore_size-1:120]
    
    # set initial detritus biomass density (g.m^-3) (W)
    detritus.loc[{'time': detritus.time.min()}] = params['init_detritus'][0]

    if use_init:
        # set up with the initial values from previous run
        predators.isel(size_class = slice((ind_min_pred_size-1), None)).\
        loc[{'time': predators.time.min()}] = \
        params['plank_pred_sizes'][(ind_min_pred_size-1):]
        
        # set initial detritivore spectrum from previous run
        detritivores.isel(size_class = slice((ind_min_detritivore_size-1), None)).\
        loc[{'time': detritivores.time.min()}] = \
        np.array(params['detritivore_sizes'])[(ind_min_detritivore_size-1):]

    #intrinsic natural mortality (OM.u, OM.v)
    other_mort_pred = other_mort_det = params['natural_mort']*10**(-0.25*log10_size_bins)
    # other_mort_pred = params['natural_mort']*10**(-0.25*log10_size_bins)
    
    #senescence mortality rate to limit large fish from building up in the system
    #same function as in Law et al 2008, with chosen parameters gives similar M2 values
    #as in Hall et al. 2006 (SM.u, SM.v)
    senes_mort_pred = senes_mort_det = (params['const_senescence_mort']*
                                        10**(params['exp_senescence_mort']*
                                             (log10_size_bins - params['size_senescence'])))

    #Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!) (FVec.u, FVec.v)
    # from Benoit & Rochet 2004 
    # here fish_mort_pred and fish_mort_pred= fixed catchability term for predators and 
    # detritivores to be estimated along with ind_min_det and ind_min_fish_pred
    fishing_mort_pred.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).\
    loc[{'time': fishing_mort_pred.time.min()}] = (params['fish_mort_pred']*
                                                   effort.sel(time = dbpm_input.time.min()).values)
    
    fishing_mort_det.isel(size_class = slice(int(ind_min_fish_det-1), -1)).\
    loc[{'time': fishing_mort_det.time.min()}] = (params['fish_mort_detritivore']*
                                                  effort.sel(time = dbpm_input.time.min()).values)

    #output fisheries catches per yr at size - predators (Y_u)
    catch_pred.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).\
    loc[{'time': catch_pred.time.min()}] = \
    (fishing_mort_pred.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).
     sel(time = fishing_mort_pred.time.min())*
     predators.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).
     sel(time = predators.time.min())*size_bin_vals[int(ind_min_fish_pred-1): -1])
    
    #output fisheries catches per yr at size - detritivores (Y_v)
    catch_det.isel(size_class = slice(int(ind_min_fish_det-1), -1)).\
    loc[{'time': catch_pred.time.min()}] = \
    (fishing_mort_det.isel(size_class = slice(int(ind_min_fish_det-1), -1)).
     sel(time = fishing_mort_pred.time.min())*
     detritivores.isel(size_class = slice(int(ind_min_fish_det-1), -1)).
     sel(time = predators.time.min())*size_bin_vals[int(ind_min_fish_det-1): -1])

    if temp_effect:
        #Adding time dimension to temperature effect for pelagic group
        pel_tempeffect = xr.DataArray(data = np.exp(params['c1']-params['activation_energy']/
                                                    (params['boltzmann']*
                                                     (np.array(params['sea_surf_temp'])+273))),
                                      dims = ['time'],
                                      coords = {'time': dbpm_input.time})
        
        #Adding time dimension to temperature effect for benthic group
        ben_tempeffect = xr.DataArray(data = np.exp(params['c1']-params['activation_energy']/
                                                    (params['boltzmann']*
                                                     (np.array(params['sea_floor_temp'])+273))),
                                      dims = ['time'],
                                      coords = {'time': dbpm_input.time})
    else:
        pel_tempeffect = 1
        ben_tempeffect = 1
    
    #To be applied to feeding rates for pelagics and benthic groups
    feed_multiplier = params['hr_volume_search']*10**(log10_size_bins*metabolic_req_pred)
    
    # Ingested food
    growth_prop = 1-np.array(params['defecate_prop'])
    # High quality food
    high_prop = 1-np.array(params['def_low'])

    #iteration over time, N [days]
    for i in range(0, len(dbpm_input)):
        t = dbpm_input.time[i]
        ts = time[i]
        # Calculate Growth and Mortality
        # feeding rates pelagics yr-1 (f_pel)
        pred_growth = np.dot((predators.sel(time = ts)*params['log_size_increase']), 
                             constant_growth)
        feed_rate_pel = (pel_tempeffect.sel(time = t).values*\
                         (((feed_multiplier*params['pref_pelagic'])*pred_growth)/
                          (1+params['handling']*((feed_multiplier*params['pref_pelagic'])
                                                 *pred_growth))))
        
        # feeding rates benthics yr-1 (f_ben)
        detrit_growth = np.dot((detritivores.sel(time = ts)*params['log_size_increase']), 
                               constant_growth)
        feed_rate_bent = (pel_tempeffect.sel(time = t).values*\
                         (((feed_multiplier*params['pref_benthos'])*detrit_growth)/
                          (1+params['handling']*((feed_multiplier*params['pref_benthos'])
                                                 *detrit_growth))))
        
        # feeding rates detritivores yr-1 (f_det)
        detritus_multiplier = ((1/size_bin_vals)*params['hr_vol_filter_benthos']*
                               (10**(log10_size_bins*metabolic_req_detritivore)*
                               detritus.sel(time = ts).values))
        feed_rate_det = (ben_tempeffect.sel(time = t).values*detritus_multiplier/
                         (1+params['handling']*detritus_multiplier))
        
        # Predator growth integral yr-1 (GG_u)
        growth_int_pred.loc[{'time': ts}] = (growth_prop*np.array(params['growth_pred'])*
                                             feed_rate_pel+high_prop*
                                             np.array(params['growth_detritivore'])*feed_rate_bent)
        
        # Reproduction yr-1 (R_u)
        if params['dynamic_reproduction']:
            reprod_pred.loc[{'time': ts}] = (growth_prop*(1-(np.array(params['growth_pred'])+
                                                             params['energy_pred']))*
                                             feed_rate_pel+growth_prop*
                                             (1-(np.array(params['growth_detritivore'])+
                                                 params['energy_detritivore']))*feed_rate_bent)
        
        # Predator death integrals 
        #Satiation level of predator for pelagic prey
        divisor = feed_multiplier*params['pref_pelagic']*pred_growth
        sat_pel = [fp/fpt if fp > 0 else 0 for fp, fpt in zip(feed_rate_pel, divisor)]
        
        # yr-1 (PM_u)
        pred_mort_pred.loc[{'time': ts}] = (np.array(params['pref_pelagic'])*
                                            params['hr_volume_search']*met_req_log10_size_bins*
                                           np.dot((predators.sel(time = ts)*sat_pel*
                                                   params['log_size_increase']), 
                                                  constant_mortality))
        
        # yr-1 (Z_u)
        tot_mort_pred.loc[{'time': ts}] = (pred_mort_pred.sel(time = ts)+
                                           pel_tempeffect.sel(time = t).values*
                                           other_mort_pred+senes_mort_pred+
                                           fishing_mort_pred.sel(time = ts))
        
        # Benthos growth integral yr-1 (GG_v)
        growth_det.loc[{'time': ts}] = (np.array(high_prop)*params['growth_detritus']*
                                        feed_rate_det)
        
        #reproduction yr-1 (R_v)
        if params['dynamic_reproduction']:
            reprod_det.loc[{'time': ts}] = (high_prop*(1-(np.array(params['growth_detritus'])+
                                                      params['energy_detritivore']))*
                                        feed_rate_det)
        
        # Benthos death integral
        # Satiation level of predator for benthic prey 
        divisor = (params['hr_volume_search']*10**(log10_size_bins*
                                                   metabolic_req_detritivore)*
                    params['pref_benthos']*detrit_growth)
        sat_ben = [fb/fpt if fb > 0 else 0 for fb, fpt in zip(feed_rate_bent, divisor)]
        
        # yr-1 (PM_v)
        new_mort = ((np.array(params['pref_benthos'])*params['hr_volume_search']*
                     met_req_log10_size_bins)*
                    np.dot((predators.sel(time = ts)*sat_ben*params['log_size_increase']),
                           constant_mortality))
        pred_mort_det.loc[{'time': ts}] = [nm if sb > 0 else 0 for sb, nm in zip(sat_ben, new_mort)]
        
        # yr-1 (Z_v)
        tot_mort_det.loc[{'time': ts}] = (pred_mort_det.sel(time = ts)+
                                          ben_tempeffect.sel(time = t).values*
                                          other_mort_det+senes_mort_det+
                                          fishing_mort_det.sel(time = ts))
        
        #total biomass density eaten by pred (g.m-2.yr-1)
        eatenbypred = (size_bin_vals*feed_rate_pel*predators.sel(time = ts)+size_bin_vals*
                       feed_rate_bent*predators.sel(time = ts))
        
        #detritus output (g.m-2.yr-1)
        # losses from detritivore scavenging/filtering only:
        output_w = (size_bin_vals*feed_rate_det*detritivores.sel(time = ts)*
                    params['log_size_increase']).sum().values
        
        #total biomass density defecated by pred (g.m-2.yr-1)
        defbypred = (params['defecate_prop']*(feed_rate_pel)*size_bin_vals*
                     predators.sel(time = ts)+params['def_low']*feed_rate_bent*
                     size_bin_vals*predators.sel(time = ts))
        
        # Increment values of detritus, predators & detritivores for next time step  
        #Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of 
        #detritus pool and solve for detritus biomass density in next time step 
        if not ERSEM_det_input:
            #considering pelagic faeces as input as well as dead bodies from both 
            #pelagic and benthic communities and phytodetritus (dying sinking
            #phytoplankton)
            if params['detritus_coupling']:
                # pelagic spectrum inputs (sinking dead bodies and faeces) - export 
                # ratio used for "sinking rate" + benthic spectrum inputs (dead stuff
                # already on/in seafloor)
                input_w = (sinking_rate.sel(time = t)* 
                   ((defbypred.isel(size_class = slice(ind_min_pred_size-1, 
                                                       numb_size_bins))*
                     params['log_size_increase']).sum()+
                    (pel_tempeffect.sel(time = t).values*other_mort_pred*
                     predators.sel(time = ts)*size_bin_vals*
                     params['log_size_increase']).sum()+ 
                    (pel_tempeffect.sel(time = t).values*senes_mort_pred*
                     predators.sel(time = ts)*size_bin_vals*
                     params['log_size_increase']).sum())+
                   ((ben_tempeffect.sel(time = t).values*other_mort_det*
                     detritivores.sel(time = ts)*size_bin_vals*
                     params['log_size_increase']).sum()+ 
                    (ben_tempeffect.sel(time = t).values*senes_mort_det*
                     detritivores.sel(time = ts)*size_bin_vals*
                     params['log_size_increase']).sum()))
            else:
                input_w = ((ben_tempeffect.sel(time = t).values*other_mort_det*
                            detritivores.sel(time = ts)*size_bin_vals*
                            params['log_size_increase']).sum()+
                           (ben_tempeffect.sel(time = t).values*senes_mort_det*
                            detritivores.sel(time = ts)*size_bin_vals*
                            params['log_size_increase']).sum())
        
            # get burial rate from Dunne et al. 2007 equation 3
            burial = input_w*(0.013+0.53*input_w**2/(7+input_w)**2)
        
            # losses from detritivory + burial rate (not including remineralisation
            # bc that goes to p.p. after sediment, we are using realised p.p. as
            # inputs to the model) 
            dW = input_w-(output_w+burial)
        
            #biomass density of detritus g.m-2
            detritus.loc[{'time': t}] = ((detritus.sel(time = ts).values+
                                          dW.values*params['timesteps_years']).item())
        else:
            detritus.loc[{'time': t}] = detritus.sel(time = ts)
        
        # Pelagic Predator Density (nos.m-2)- solve for time + timesteps_years using
        # implicit time Euler upwind finite difference (help from Ken Andersen 
        # and Richard Law)
        
        # Matrix setup for implicit differencing 
        Ai_u = np.zeros(len(log10_size_bins))
        Bi_u = np.zeros(len(log10_size_bins))
        Si_u = np.zeros(len(log10_size_bins))
        
        Ai_u[1:] = ((1/np.log(10))*
                    -growth_int_pred.isel(size_class = slice(None, -1)).sel(time = ts)*
                    params['timesteps_years']/params['log_size_increase'])
        Bi_u[1:] = (1+(1/np.log(10))*
                    growth_int_pred.isel(size_class = slice(1, None)).sel(time = ts)*
                    params['timesteps_years']/params['log_size_increase']+
                    tot_mort_pred.isel(size_class = slice(1, None)).sel(time = ts)*
                    params['timesteps_years'])
        Si_u[1:] = predators.isel(size_class = slice(1, None)).sel(time = ts)
        
        # Boundary condition at upstream end 
        Ai_u[(ind_min_pred_size-1)] = 0
        Bi_u[(ind_min_pred_size-1)] = 1
        Si_u[(ind_min_pred_size-1)] = predators.isel(size_class = (ind_min_pred_size-1)).\
        sel(time = ts)
        
        # Invert matrix
        #recruitment at smallest consumer mass
        #continuation of plankton hold constant  
        if use_init:
            predators.isel(size_class = slice(None, ind_min_pred_size)).\
            loc[{'time': t}] = (ui0.sel(time = t).values*
                                (10**(slope_phy_zoo.sel(time = t).values*
                                      log10_size_bins))[:ind_min_pred_size])
        else:
            predators.isel(size_class = slice(None, ind_min_pred_size)).\
            loc[{'time': t}] = (ui0.sel(time = t+pd.DateOffset(months = 1)).values*
                                (10**(slope_phy_zoo.sel(time = t+pd.DateOffset(months = 1)).values*
                                      log10_size_bins))[:ind_min_pred_size])
        
        # apply transfer efficiency of 10% *plankton density at same size
        # reproduction from energy allocation
        if params['dynamic_reproduction']:
            predators.isel(size_class = ind_min_pred_size-1).\
            loc[{'time': t}] = (predators.isel(size_class = ind_min_pred_size-1).\
                                sel(time = ts).values+
                                ((reprod_pred.isel(size_class = slice(ind_min_pred_size,
                                                                      None)).sel(time = ts)*
                                  size_bin_vals[ind_min_pred_size:]*
                                  predators.isel(size_class = slice(ind_min_pred_size, 
                                                                    None)).sel(time = ts)*
                                  params['log_size_increase']).sum().values*
                                 params['timesteps_years'])/
                                (np.array(params['log_size_increase'])*
                                 size_bin_vals[ind_min_pred_size-1])-
                                (np.array(params['timesteps_years'])/
                                 params['log_size_increase'])*(1/np.log(10))*
                                (growth_int_pred.isel(size_class = ind_min_pred_size-1).\
                                 sel(time = ts)*
                                 predators.isel(size_class = ind_min_pred_size-1).\
                                 sel(time = ts)).values-params['timesteps_years']*
                                ((tot_mort_pred.isel(size_class = ind_min_pred_size-1).\
                                  sel(time = ts)*
                                  predators.isel(size_class = ind_min_pred_size-1).\
                                  sel(time = ts)).values)).item()
        #main loop calculation
        for j in range(ind_min_pred_size, numb_size_bins):
            predators.isel(size_class = j).\
            loc[{'time': t}] = ((Si_u[j]-Ai_u[j]*predators.isel(size_class = j-1).\
                                 sel(time = t))/Bi_u[j])

        # Benthic Detritivore Density (nos.m-2) 
        Ai_v = np.zeros(len(log10_size_bins))
        Bi_v = np.zeros(len(log10_size_bins))
        Si_v = np.zeros(len(log10_size_bins))
        
        #shorthand for matrix referencing
        idx = np.arange(ind_min_detritivore_size, numb_size_bins)
        
        Ai_v[idx] = (((1/np.log(10))*
                      -growth_det.isel(size_class = slice(ind_min_detritivore_size-1, -1)).\
                      sel(time = ts)).values*params['timesteps_years']/
                     params['log_size_increase'])
        Bi_v[idx] = (1+(1/np.log(10))*
                     growth_det.isel(size_class = slice(ind_min_detritivore_size, None)).\
                     sel(time = ts).values*params['timesteps_years']/
                     params['log_size_increase']+
                     tot_mort_det.isel(size_class = slice(ind_min_detritivore_size, None)).\
                     sel(time = ts)*params['timesteps_years'])
        Si_v[idx] = (detritivores.isel(size_class = slice(ind_min_detritivore_size, None)).\
                     sel(time = ts))
        
        #boundary condition at upstream end
        Ai_v[ind_min_detritivore_size-1] = 0
        Bi_v[ind_min_detritivore_size-1] = 1
        Si_v[ind_min_detritivore_size-1] = (detritivores.\
                                            isel(size_class = ind_min_detritivore_size-1).\
                                            sel(time = ts))
        
        #invert matrix
        #recruitment at smallest detritivore mass  
        #hold constant continuation of plankton with sinking rate multiplier 
        (detritivores.isel(size_class = slice(None, 
                                             ind_min_detritivore_size)).\
         loc[{'time': t}]) = (detritivores.isel(size_class = slice(None,
                                                                  ind_min_detritivore_size)).\
                                               sel(time = ts))
        
        # apply a very low of transfer efficiency 1%* total biomass of detritus
        #divided by minimum size
        if params['dynamic_reproduction']:
            (detritivores.\
             isel(size_class = ind_min_detritivore_size-1).\
             loc[{'time': t}]) = (detritivores.isel(size_class = ind_min_detritivore_size-1).\
                                  sel(time = ts).values+
                                  ((reprod_det.isel(size_class = slice(ind_min_detritivore_size,
                                                                       None)).sel(time = ts)*
                                    size_bin_vals[ind_min_detritivore_size:]*
                                    detritivores.isel(size_class = slice(ind_min_detritivore_size,
                                                                         None)).sel(time = ts))\
                                   .values*params['log_size_increase']).sum()*
                                  np.array(params['timesteps_years'])/
                                  (np.array(params['log_size_increase'])*
                                   size_bin_vals[ind_min_detritivore_size-1])-
                                  (np.array(params['timesteps_years'])/
                                   params['log_size_increase'])*(1/np.log(10))*
                                  (growth_det.isel(size_class = ind_min_detritivore_size-1).\
                                   sel(time = ts)*
                                   detritivores.isel(size_class = ind_min_detritivore_size-1).\
                                   sel(time = ts)).values-
                                  (params['timesteps_years']*
                                   (tot_mort_det.isel(size_class = ind_min_detritivore_size-1).\
                                    sel(time = ts)*
                                    detritivores.isel(size_class = ind_min_detritivore_size-1).\
                                    sel(time = ts)).values)).item()

        #loop calculation
        for j in range(ind_min_pred_size, numb_size_bins):
            detritivores.isel(size_class = j).\
            loc[{'time': t}] = ((Si_v[j]-Ai_v[j]*detritivores.isel(size_class = j-1).\
                                 sel(time = t))/Bi_v[j])

        # increment fishing 
        (fishing_mort_pred.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).\
         loc[{'time': t}]) = (params['fish_mort_pred']*
                              effort.sel(time = t+pd.DateOffset(months = 1)).values)
        (fishing_mort_det.isel(size_class = slice(int(ind_min_fish_det-1), -1)).\
         loc[{'time': t}]) = (params['fish_mort_detritivore']*
                              effort.sel(time = t+pd.DateOffset(months = 1)).values)
        
        #output fisheries catches per yr at size
        catch_pred.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).\
         loc[{'time': t}] = ((fishing_mort_pred.isel(size_class = slice(int(ind_min_fish_pred-1), 
                                                                        -1)).\
                              sel(time = t)*
                              predators.isel(size_class = slice(int(ind_min_fish_pred-1), -1)).\
                              sel(time = t)).values*size_bin_vals[int(ind_min_fish_pred-1): -1])
        
        #output fisheries catches per yr at size
        catch_det.isel(size_class = slice(int(ind_min_fish_det-1), -1)).\
         loc[{'time': t}] = ((fishing_mort_det.isel(size_class = slice(int(ind_min_fish_det-1), 
                                                                       -1)).\
                              sel(time = t)*
                              detritivores.isel(size_class = slice(int(ind_min_fish_det-1), -1)).\
                              sel(time = t)).values*size_bin_vals[int(ind_min_fish_det-1): -1])

    return_list = {'predators': predators,
               'growth_int_pred': growth_int_pred,
               'pred_mort_pred': pred_mort_pred,
               'detritivores': detritivores,
               'growth_det': growth_det,
               'pred_mort_det': pred_mort_det,
               'catch_pred': catch_pred,
               'catch_det': catch_det,
               'detritus': detritus,
               'params': params}