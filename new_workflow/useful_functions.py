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
def extract_gfdl(file_path, mask, path_out, drop = 'all'):
    '''
    Inputs:
    - file_path (character) File path where GFDL zarr file is located
    - mask (boolean data array) Grid cells within region of interest should be identified
    as 1.
    - path_out (character) File path where outputs should be stored as zarr files
    - drop (character) Default is 'all'. It drops all grid cells where condition is not met. 
    If a dimension name is provided, it will drop rows if ALL grid cells along dimension
    do not meet the condition.

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    da = xr.open_zarr(file_path)

    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    da = da[var]
    
    #Apply mask and remove rows where all grid cells are empty to reduce data array size
    if drop == 'all':
        da = da.where(mask == 1, drop = True)
    else:
        da = da.where(mask == 1).dropna(drop, how = 'all')

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


#Format gridded_params in Python friendly way
def gridded_param_python(gridded_params):
    '''
    Inputs:
    - gridded_params (dictionary) Containing the output from the `sizeparam` function in R.
    Stored as a `json` file.

    Outputs:
    griddded_python (dictionary) Containing gridded_params in a Python friendly format. New
    entries calculated as needed, and entries not used in the DBPM gridded model were removed.
    '''

    gridded_python = {'timesteps_years': gridded_params['timesteps_years'][0],
                      'numb_time_steps': gridded_params['numb_time_steps'][0],
                      'effort': gridded_params['effort'],
                      'fish_mort_pred': gridded_params['fish_mort_pred'][0],
                      'fish_mort_detritivore': gridded_params['fish_mort_detritivore'][0],
                      'hr_volume_search': gridded_params['hr_volume_search'][0],
                      'detritus_coupling': gridded_params['detritus_coupling'][0],
                      'log10_pred_prey_ratio': gridded_params['log10_pred_prey_ratio'][0],
                      'log_prey_pref': gridded_params['log_prey_pref'][0],
                      'hr_vol_filter_benthos': gridded_params['hr_vol_filter_benthos'][0],
                      'metabolic_req_pred': gridded_params['metabolic_req_pred'][0],
                      'metabolic_req_detritivore': gridded_params['metabolic_req_detritivore'][0],
                      'defecate_prop': gridded_params['defecate_prop'][0],
                      'growth_prop': 1-(gridded_params['defecate_prop'][0]),
                      'def_low': gridded_params['def_low'][0],
                      'high_prop': 1-(gridded_params['def_low'][0]),
                      'growth_pred': gridded_params['growth_pred'][0],
                      'growth_detritivore': gridded_params['growth_detritivore'][0],
                      'growth_detritus': gridded_params['growth_detritus'][0],
                      'energy_pred': gridded_params['energy_pred'][0],
                      'energy_detritivore': gridded_params['energy_detritivore'][0],
                      'handling': gridded_params['handling'][0],
                      'dynamic_reproduction': gridded_params['dynamic_reproduction'][0],
                      'c1': gridded_params['c1'][0],
                      'activation_energy': gridded_params['activation_energy'][0],
                      'boltzmann': gridded_params['boltzmann'][0],
                      'natural_mort': gridded_params['natural_mort'][0],
                      'size_senescence': gridded_params['size_senescence'][0],
                      'exp_senescence_mort': gridded_params['exp_senescence_mort'][0],
                      'const_senescence_mort': gridded_params['const_senescence_mort'][0],
                      'log_size_increase': gridded_params['log_size_increase'][0],
                      'log10_size_bins': gridded_params['log10_size_bins'],
                      'numb_size_bins':  gridded_params['numb_size_bins'][0],
                      'ind_min_pred_size': (gridded_params['ind_min_pred_size'][0])-1,
                      'ind_min_detritivore_size': (gridded_params['ind_min_detritivore_size'][0])-1,
                      'idx_new': list(range(gridded_params['ind_min_detritivore_size'][0]+1,
                                            gridded_params['numb_size_bins'][0])),
                      'ind_min_fish_pred': int(gridded_params['ind_min_fish_pred'][0]-1),
                      'ind_min_fish_det': int(gridded_params['ind_min_fish_det'][0]-1),
                      'idx': (np.array(gridded_params['idx'])-1).tolist(),
                      'init_pred': gridded_params['init_pred'],
                      'init_detritivores': gridded_params['init_detritivores'],
                      'init_detritus': gridded_params['init_detritus'][0]
                     }
    
    return gridded_python


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
    expax = 10**(log10_size_bins*metabolic_req_pred)
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


# Gravity model
# Redistribute total effort across grid cells according to proportion of biomass in that 
# grid cell using graivity model, Walters & Bonfil, 1999, Gelchu & Pauly 2007 ideal free
# distribution - Blanchard et al 2008
def gravitymodel(effort, prop_b, depth, n_iter):
    '''
    Inputs:
    - effort (Data array) Fishing effort for a single time step
    - prop_b (Data array) Proportion of total fishable biomass for each grid cell at a 
    single time step
    - depth (Data array) Bathymetry of the area of interest
    - n_iter (integer) Number of iterations needed to redistribute fishing effort
    
    Outputs:
    eff (data array) Containing redistributed fishing effort
    '''

    eff = effort
    d = (1-(depth/depth.max()))

    #Initialise loop
    i = 1
    while(i <= n_iter):
        suit = prop_b*d
        rel_suit = suit/suit.sum()
        neweffort = eff+(rel_suit*eff)
        mult = eff.sum()/neweffort.sum()
        eff = neweffort*mult
        i += i

    return eff


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
    
    #data arrays for keeping track of growth (GG_v, GG_u) and reproduction (R_v, R_u) 
    #from ingested food:
    reprod_det = xr.zeros_like(predators)
    reprod_pred = xr.zeros_like(predators)
    growth_det = xr.zeros_like(predators)
    growth_int_pred = xr.zeros_like(predators)
    
    #data arrays for keeping track of predation mortality (PM_v, PM_u)
    pred_mort_det = xr.zeros_like(predators)
    pred_mort_pred = xr.zeros_like(predators)
        
    #data arrays for keeping track of total mortality (Z_v, Z_u)
    tot_mort_det = xr.zeros_like(predators)
    tot_mort_pred = xr.zeros_like(predators)
    
    #To be applied to feeding rates for pelagics and benthic groups
    feed_mult_pel = (params['hr_volume_search']*
                     10**(log10_size_bins_mat*metabolic_req_pred)*
                     params['pref_pelagic'])

    feed_mult_ben = (params['hr_volume_search']*
                     10**(log10_size_bins_mat*metabolic_req_pred)*
                     params['pref_benthos'])

    #iteration over time, N [days]
    for i in range(0, numb_time_steps):
        t = dbpm_input_time[i]
        ts = time[i]
        #Select relevant timestep from predators and detritivores
        pred_short = predators.sel(time = ts)
        detrit_short = detritivores.sel(time = ts)
        
        # Calculate Growth and Mortality
        # feeding rates pelagics yr-1 (f_pel)
        pred_growth = (feed_mult_pel*np.dot((pred_short*log_size_increase),
                                            constant_growth))
        feed_rate_pel = (pel_tempeffect.sel(time = t)*
                         (pred_growth/(1+params['handling']*pred_growth)))
        
        # feeding rates benthics yr-1 (f_ben)
        detrit_growth = (feed_mult_ben*np.dot((detrit_short*log_size_increase),
                                              constant_growth))
        feed_rate_bent = (pel_tempeffect.sel(time = t)*
                          (detrit_growth/(1+params['handling']*detrit_growth)))
        
        # feeding rates detritivores yr-1 (f_det)
        detritus_multiplier = ((1/size_bin_vals)*params['hr_vol_filter_benthos']*
                               10**(log10_size_bins_mat*metabolic_req_detritivore)*
                               detritus.sel(time = ts))
        feed_rate_det = (ben_tempeffect.sel(time = t)*detritus_multiplier/
                         (1+params['handling']*detritus_multiplier))
        
        # Predator growth integral yr-1 (GG_u)
        growth_int_pred.loc[{'time': ts}] = (growth_prop*params['growth_pred']*
                                             feed_rate_pel+high_prop*
                                             params['growth_detritivore']*feed_rate_bent)
        
        # Reproduction yr-1 (R_u)
        if params['dynamic_reproduction']:
            reprod_pred.loc[{'time': ts}] = (growth_prop*
                                             (1-(np.array(params['growth_pred'])+
                                                 params['energy_pred']))*
                                             feed_rate_pel+growth_prop*
                                             (1-(np.array(params['growth_detritivore'])+
                                                 params['energy_detritivore']))*
                                             feed_rate_bent)
        
        # Predator death integrals 
        #Satiation level of predator for pelagic prey
        sat_pel = xr.where(feed_rate_pel > 0, feed_rate_pel/pred_growth, 0)
        
        # yr-1 (PM_u)
        pred_mort_pred.loc[{'time': ts}] = ((params['pref_pelagic']*
                                             met_req_log10_size_bins*
                                             params['hr_volume_search'])*
                                             np.dot((pred_short*sat_pel*
                                                     log_size_increase), 
                                                    constant_mortality))
        
        # yr-1 (Z_u)
        tot_mort_pred.loc[{'time': ts}] = (pred_mort_pred.sel(time = ts)+
                                           pel_tempeffect.sel(time = t)*
                                           other_mort_pred+senes_mort_pred+
                                           fishing_mort_pred.sel(time = t))
        
        # Benthos growth integral yr-1 (GG_v)
        growth_det.loc[{'time': ts}] = high_prop*params['growth_detritus']*feed_rate_det
        
        #reproduction yr-1 (R_v)
        if params['dynamic_reproduction']:
            reprod_det.loc[{'time': ts}] = (high_prop*
                                            (1-(np.array(params['growth_detritus'])+
                                                params['energy_detritivore']))*
                                            feed_rate_det)
        
        # Benthos death integral
        # Satiation level of predator for benthic prey 
        divisor = ((params['hr_volume_search']*
                    10**(log10_size_bins_mat*metabolic_req_detritivore)*
                    params['pref_benthos'])*
                   np.dot((detrit_short*log_size_increase), 
                          constant_growth))
        sat_ben = xr.where(feed_rate_bent > 0, feed_rate_bent/divisor, 0)
        
        # yr-1 (PM_v)
        pred_mort_det.loc[{'time': ts}] = xr.where(sat_ben > 0, 
                                                   ((params['pref_benthos']*
                                                     met_req_log10_size_bins*
                                                     params['hr_volume_search'])*
                                                    np.dot((pred_short*sat_ben*
                                                            log_size_increase),
                                                           constant_mortality)), 0)
        
        # yr-1 (Z_v)
        tot_mort_det.loc[{'time': ts}] = (pred_mort_det.sel(time = ts)+
                                          ben_tempeffect.sel(time = t)*
                                          other_mort_det+senes_mort_det+
                                          fishing_mort_det.sel(time = t))
        
        #detritus output (g.m-2.yr-1)
        # losses from detritivore scavenging/filtering only:
        output_w = (size_bin_vals*feed_rate_det*detrit_short*log_size_increase).sum()
        
        #total biomass density defecated by pred (g.m-2.yr-1)
        defbypred = (params['defecate_prop']*feed_rate_pel*size_bin_vals*
                     pred_short+params['def_low']*feed_rate_bent*size_bin_vals*
                     pred_short)
        
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
                           ((defbypred.isel(size_class = slice(ind_min_pred_size,
                                                               numb_size_bins))*
                             log_size_increase).sum()+
                            (pel_tempeffect.sel(time = t)*other_mort_pred*
                             pred_short*size_bin_vals*log_size_increase).sum()+
                            (pel_tempeffect.sel(time = t)*senes_mort_pred*
                             pred_short*size_bin_vals*log_size_increase).sum())+
                           ((ben_tempeffect.sel(time = t)*other_mort_det*
                             detrit_short*size_bin_vals*log_size_increase).sum()+
                            (ben_tempeffect.sel(time = t)*senes_mort_det*
                             detrit_short*size_bin_vals*log_size_increase).sum()))
            else:
                input_w = ((ben_tempeffect.sel(time = t)*other_mort_det*
                            detrit_short*size_bin_vals*log_size_increase).sum()+
                           (ben_tempeffect.sel(time = t)*senes_mort_det*detrit_short*
                            size_bin_vals*log_size_increase).sum())
        
            # get burial rate from Dunne et al. 2007 equation 3
            burial = input_w*(0.013+0.53*input_w**2/(7+input_w)**2)
        
            # losses from detritivory + burial rate (not including remineralisation
            # bc that goes to p.p. after sediment, we are using realised p.p. as
            # inputs to the model) 
            dW = input_w-(output_w+burial)
        
            #biomass density of detritus g.m-2
            detritus.loc[{'time': t}] = detritus.sel(time = ts)+dW*timesteps_years
        else:
            detritus.loc[{'time': t}] = detritus.sel(time = ts)
        
        # Pelagic Predator Density (nos.m-2)- solve for time + timesteps_years using
        # implicit time Euler upwind finite difference (help from Ken Andersen 
        # and Richard Law)
        
        # Matrix setup for implicit differencing 
        Ai_u = np.zeros(numb_size_bins)
        Bi_u = np.zeros(numb_size_bins)
        Si_u = np.zeros(numb_size_bins)
        
        Ai_u[idx] = ((1/np.log(10))*
                     -growth_int_pred.isel(size_class = slice(None, -1)).
                     sel(time = ts)*timesteps_years/log_size_increase)
        Bi_u[idx] = (1+(1/np.log(10))*
                     growth_int_pred.isel(size_class = slice(1, None)).
                     sel(time = ts)*timesteps_years/log_size_increase+
                     tot_mort_pred.isel(size_class = slice(1, None)).
                     sel(time = ts)*timesteps_years)
        Si_u[idx] = pred_short.isel(size_class = slice(1, None))
        
        # Boundary condition at upstream end
        Ai_u[ind_min_pred_size] = 0
        Bi_u[ind_min_pred_size] = 1
        Si_u[ind_min_pred_size] = pred_short.isel(size_class = ind_min_pred_size)
        
        # apply transfer efficiency of 10% *plankton density at same size
        # reproduction from energy allocation
        if params['dynamic_reproduction']:
            predators.isel(size_class = ind_min_pred_size).\
            loc[{'time': t}] = (pred_short.isel(size_class = ind_min_pred_size)+
                                ((reprod_pred.sel(time = ts)*size_bin_vals*
                                  pred_short*log_size_increase).
                                 isel(size_class = slice(ind_min_pred_size+1, None)).
                                 sum()*timesteps_years)/
                                (log_size_increase*
                                 size_bin_vals.isel(size_class = ind_min_pred_size))-
                                (timesteps_years/log_size_increase)*(1/np.log(10))*
                                ((growth_int_pred.sel(time = ts)*pred_short).
                                 isel(size_class = ind_min_pred_size))-
                                timesteps_years*
                                ((tot_mort_pred.sel(time = ts)*pred_short).
                                 isel(size_class = ind_min_pred_size)))

        #main loop calculation
        for j in range((ind_min_pred_size+1), numb_size_bins):
            predators.isel(size_class = j).\
            loc[{'time': t}] = ((Si_u[j]-Ai_u[j]*predators.sel(time = t).
                                 isel(size_class = j-1))/Bi_u[j])

        # Benthic Detritivore Density (nos.m-2) 
        Ai_v = np.zeros(numb_size_bins)
        Bi_v = np.zeros(numb_size_bins)
        Si_v = np.zeros(numb_size_bins)
        
        #shorthand for matrix referencing
        Ai_v[idx_new] = ((1/np.log(10))*
                         -growth_det.isel(size_class = idx_new-1).sel(time = ts)*
                         timesteps_years/log_size_increase)
        Bi_v[idx_new] = (1+(1/np.log(10))*
                         growth_det.isel(size_class = idx_new).sel(time = ts)*
                         timesteps_years/log_size_increase+
                         tot_mort_det.isel(size_class = idx_new).
                         sel(time = ts)*timesteps_years)
        Si_v[idx_new] = detrit_short.isel(size_class = idx_new)
        
        #boundary condition at upstream end
        Ai_v[ind_min_detritivore_size] = 0
        Bi_v[ind_min_detritivore_size] = 1
        Si_v[ind_min_detritivore_size] = (detrit_short.
                                          isel(size_class = ind_min_detritivore_size))
        
        #invert matrix
        #recruitment at smallest detritivore mass  
        #hold constant continuation of plankton with sinking rate multiplier 
        (detritivores.isel(size_class = slice(None, ind_min_detritivore_size+1)).
         loc[{'time': t}]) = (detrit_short.
                              isel(size_class = slice(None, 
                                                      ind_min_detritivore_size+1)))
        
        # apply a very low of transfer efficiency 1%* total biomass of detritus
        #divided by minimum size
        if params['dynamic_reproduction']:
            (detritivores.isel(size_class = ind_min_detritivore_size).
             loc[{'time': t}]) = (detrit_short.
                                  isel(size_class = ind_min_detritivore_size)+
                                  ((reprod_det.sel(time = ts)*size_bin_vals*
                                    detrit_short*log_size_increase).
                                   isel(size_class = idx_new).sum()*timesteps_years)/
                                  (log_size_increase*size_bin_vals.
                                   isel(size_class = ind_min_detritivore_size))-
                                  (timesteps_years/log_size_increase)*(1/np.log(10))*
                                  ((growth_det.sel(time = ts)*detrit_short).
                                   isel(size_class = ind_min_detritivore_size))-
                                  timesteps_years*
                                  ((tot_mort_det.sel(time = ts)*detrit_short).
                                   isel(size_class = ind_min_detritivore_size)))

        #loop calculation
        for j in idx_new:
            detritivores.isel(size_class = j).\
            loc[{'time': t}] = ((Si_v[j]-Ai_v[j]*detritivores.isel(size_class = j-1).
                                 sel(time = t))/Bi_v[j])

    #output fisheries catches per yr at size
    (catch_pred.isel(size_class = slice(ind_min_fish_pred, None)).
     loc[:,dbpm_input_time]) = ((fishing_mort_pred*predators*size_bin_vals).
                               isel(size_class = slice(ind_min_fish_pred, None)))
    
    #output fisheries catches per yr at size
    (catch_det.isel(size_class = slice(ind_min_fish_det, None)).
     loc[:,dbpm_input_time]) = ((fishing_mort_det*detritivores*size_bin_vals).
                               isel(size_class = slice(ind_min_fish_det, None)))

    # Subsetting predator and detritivore results to include relevant size 
    # classes
    predators = predators.isel(size_class = slice(ind_min_pred_size, None))
    detritivores = detritivores.isel(size_class = slice(ind_min_detritivore_size,
                                                        None))

    return_list = {'predators': predators, 'growth_int_pred': growth_int_pred,
                   'pred_mort_pred': pred_mort_pred, 'detritivores': detritivores,
                   'growth_det': growth_det, 'pred_mort_det': pred_mort_det,
                   'catch_pred': catch_pred, 'catch_det': catch_det, 
                   'detritus': detritus, 'params': params}
    return return_list
