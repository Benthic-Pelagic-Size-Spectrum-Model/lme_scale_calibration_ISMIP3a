###############################################################################
# Supporting DBPM functions
# Functions have been adapted from previous DBPM work
# 
# Edited by: Denisse Fierro Arcos
# Date of update: 2024-10-16


# Loading libraries -------------------------------------------------------
library(arrow)
library(janitor)
library(dplyr)
library(lhs)


dbpm_inputs <- read_parquet("/g/data/vf71/la6889/dbpm_inputs/east_antarctica/monthly_weighted/dbpm_clim-fish-inputs_fao-58_1841-2010.parquet")

#Making function reproducible
set.seed(1234)

#Construct a hypercube with random numbers
#num_iter defines number of rows in hypercube
#columns represent five specific parameters needed
fishing_params <- data.frame(randomLHS(100, 5))
colnames(fishing_params) <- c("fmort_u", "fmort_v", "fminx_u", "fminx_v", 
                              "search_vol")

#Adjust range of mi size params, others go from 0-1
fishing_params <- fishing_params |> 
  mutate(fminx_u = fminx_u*2, 
         fminx_v = fminx_v*2,
         # adjust range of search vol, others go from 0-1
         search_vol = search_vol+0.001)


test <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3, 
          xmin_consumer_v = -3, tstepspryr = 12)


# Getting DBPM model parameters ready -------------------------------------
sizeparam <- function(dbpm_inputs, fishing_params, dx = 0.1, xmin = -12, 
                      xmin_consumer_u = -7, xmin_consumer_v = -7, xmax = 6, 
                      tstepspryr = 48, Ngrid = NA, use_init = FALSE, 
                      u_initial = NA, v_initial = NA, w_initial = NA, 
                      equilibrium = FALSE){
  
  #Inputs:
  # - dbpm_inputs (data frame) Containing climate and fishing data needed as 
  # inputs for DBPM. This is produced by script "03_processing_effort_fishing_inputs.R"
  # - fishing_params (data frame) Containing fishing parameters: "f_u", "f_v",
  # "f_minu", "f_minv", and "search_vol"
  # - dx (numeric) Default value is 0.1. Size increment after discretization for 
  # integration (log body weight)
  # - xmin (numeric) Default value is -12. Minimum log10 body size of plankton
  # - xmin_consumer_u (numeric) Default value is -7. Minimum log10 body size in 
  # dynamics predators
  # - xmin_consumer_v (numeric) Default value is -7. Minimum log10 body size in
  # dynamics benthic detritivores
  # - xmax (numeric). Default value is 6. Maximum log10 body size of predators
  # - tstepspryr (numeric) Default value is 48. Number of time steps to include 
  # within a year
  # - Ngrid (numeric) Optional. Number of grid cells.
  # - use_init (boolean). Default value is FALSE. If set to TRUE, the 
  # initialisation values for predators (U_initial), detritivores (V.initial) 
  # and detritus (W.initial) will be used
  # - u_initial (numeric). Optional. Default is NA. Initialisation value for
  # predators. If provided, 'use_init' must be set to TRUE
  # - v_initial (numeric). Optional. Default is NA. Initialisation value for
  # detritivores. If provided, 'use_init' must be set to TRUE
  # - w_initial (numeric). Optional. Default is NA. Initialisation value for
  # detritus If provided, 'use_init' must be set to TRUE
  # - equilibrium (boolean). Default value is FALSE.
  #
  # Outputs:
  # er (numeric vector) Export ratio
  
  
  # Creating empty list to store DBPM parameters:
  param <- list()
  
  # depth
  param$depth <- dbpm_inputs$depth
  
  # number of years to run model (tmax)
  param$n_years <- length(unique(dbpm_inputs$year))
  
  # discretisation of year(delta_t)
  param$timesteps_years <- (1/tstepspryr)
  
  # number of time bins (Neq)
  param$numb_time_steps <- param$n_years*tstepspryr
  
  # export ratio (er)
  param$export_ratio <- dbpm_inputs$export_ratio
  # fraction of sinking detritus reaching the seafloor (from export ratio input)
  # (sinking.rate)
  param$sinking_rate <- param$export_ratio 
  
  # plankton parameters
  # intercept of phyto-zooplankton spectrum (pp)
  param$int_phy_zoo <- dbpm_inputs$intercept
  # slope of phyto-zooplankton spectrum (r.plank)
  param$slope_phy_zoo <- dbpm_inputs$slope
  
  # temperature parameters 
  # sea-surface temperature - degrees Celsius (sst)
  param$sea_surf_temp <- dbpm_inputs$tos
  # near sea-floor temperature - degrees Celsius (sft)
  param$sea_floor_temp <- dbpm_inputs$tob
  
  # get rescaled effort
  param$effort <- dbpm_inputs$nom_active_area_m2_relative
  
  # fishing parameters 
  # fishing mortality rate for predators (Fmort.u)
  param$fish_mort_pred <- fishing_params$fmort_u
  # fishing mortality rate for detritivores (Fmort.v)
  param$fish_mort_detritivore <- fishing_params$fmort_v
  # minimum log10 body size fished for predators (min.fishing.size.u)
  param$min_fishing_size_pred <- fishing_params$fminx_u 
  # minimum log10 body size fished for detritivores (min.fishing.size.v)
  param$min_fishing_size_detritivore <- fishing_params$fminx_v 
  
  # Coordinates
  if("lat" %in% names(dbpm_inputs)){
    param$lat <- dbpm_inputs$lat
    param$lon <- dbpm_inputs$lon
  }else{
    param$lat <- NA
    param$lon <- NA
  }
  
  # Benthic-pelagic coupling parameters
  
  # originally 640, but if 64 then using Quest-fish default of 64 hourly rate
  # volume searched constant m3.yr-1 for fish. need to check this value, its
  # quite large. (A.u)
  param$hr_volume_search <- fishing_params$search_vol
  
  # set predator coupling to benthos, depth dependent - 0.75 above 500 m, 0.5
  # between 500-1800 and 0 below 	1800m (suggestions of values from Clive
  # Trueman based on stable isotope work, and proportion of biomass, 	Rockall 
  # Trough studies) (pref.ben)
  param$pref_benthos <- 0.8*exp(-1/250*param$depth)
  
  # preference for pelagic prey (pref.pel)
  param$pref_pelagic <- 1-param$pref_benthos 
  
  # detritus coupling on? 1 = yes, 0 = no (det.coupling)
  param$detritus_coupling <- TRUE
  
  # feeding and energy budget parameters
  # Mean log10 predator prey mass ratio 100:1. (q0)
  param$log10_pred_prey_ratio <- 2.0
  
  # 0.6 to 1.0 for log normal prey preference function. (sd.q)
  param$log_prey_pref <- 1.0
  
  # hourly rate volume filtered constant m3*yr-1 for benthos. this value yields 
  # believable growth curve. Approximately 10 times less than A.u (A.v)
  param$hr_vol_filter_benthos <- param$hr_volume_search*0.1    
  
  # exponent for metabolic requirements plus swimming for predators (Ware et al 
  # 1978) (alpha.u)
  param$metabolic_req_pred <- 0.82  
  
  # exponent for whole organism basal (sedentary) metabolic rate (see growth.v) 
  # from Peters (1983) and Brown et al. (2004) for detritivores (alpha.v)
  param$metabolic_req_detritivore <- 0.75
  
  # fraction of ingested food that is defecated (Peters,1983) (def.high)
  param$defecate_prop <- 0.3  
  
  # low <- low quality (K) food, high <- high quality (K) food (def.low)
  param$def_low <- 0.5
  
  # net growth conversion efficiency for organisms in the "predator" spectrum
  # from Ware (1978) (K.u)
  param$growth_pred <- 0.3
  
  # net growth conversion efficiency for organisms in the "detritivore" 
  # spectrum (K.v)
  param$growth_detritivore <- 0.2
  # net growth conversion efficiency for detritus (K.d)
  param$growth_detritus <- param$growth_detritivore
  
  #fraction of energy required for maintenance & activity etc. (AM.u and AM.v)
  param$energy_pred <- 0.5
  param$energy_detritivore <- 0.7
  
  # if handling time is > 0 Type II functional response, if = 0 linear (no 
  # predator satiation) - 5.7e-7
  param$handling <- 0
  
  # dynamic reproduction on = 1, off = 0 (repro.on)
  param$dynamic_reproduction <- TRUE
  
  # constant used in Jennings et al. 2008 Proc B to standardize 
  # metabolism-temperature effects for Boltzmann equation. Derived from Simon's
  # fit to Andy Clarke's data
  param$c1 <- 25.22 
  
  # activation energy, eV (E)
  param$activation_energy <- 0.63
  
  # Boltzmann's constant (Boltzmann)
  param$boltzmann <- 8.62*10^-5
  
  # "other" mortality parameters
  
  # residual natural mortality (mu0)
  param$natural_mort <- 0.2
  
  # size at senescence (xs)
  param$size_senescence <- 3
  
  # exponent for senescence mortality (p.s)
  param$exp_senescence_mort <- 0.3
  
  # constant for senescence mortality (k.sm)
  param$const_senescence_mort <- 0.2
  
  # Parameters for numerical integration (size & time discretisation)
  # size increment after discretization for integration (log body weight) (dx)
  param$log_size_increase <- dx
  
  # minimum log10 body size of plankton (xmin)
  param$min_log10_plankton <- xmin
  
  # minimum log10 body size in dynamics predators (x1)
  param$min_log10_pred <- xmin_consumer_u
  
  # minimum log10 body size in dynamics benthic detritivores (x1.det)
  param$min_log10_detritivore <- xmin_consumer_v
  
  # maximum log10 body size of predators (xmax)
  param$max_log10_pred <- xmax
  
  # maximum log10 body size before senescence kicks in (departure form 
  # linearity) (xmax2)
  param$max_log10_senescense <- 4.0
  
  # Log10 size bins (x)
  param$log10_size_bins <- seq(xmin, xmax, dx)
  
  # Number of log10 size bins (Nx)
  param$numb_size_bins <- length(param$log10_size_bins)
  
  # index for minimum log10 predator size (ref)
  param$ind_pred <- which(log10_size_bins == param$min_log10_pred)
    #((param$min_log10_pred-param$min_log10_plankton)/dx)+1
  
  # index for minimum log10 detritivore size (ref.det)
  param$ind_detritivore <- which(log10_size_bins == param$min_log10_detritivore)
    #((param$min_log10_detritivore-param$min_log10_plankton)/dx)+1 
  
  # index in F vector corresponding to smallest size fished in U (Fref.u)
  param$ind_min_pred <- ((param$min_fishing_size_pred-param$min_log10_plankton)/
                           dx)+1
  
  # index in F vector corresponding to smallest size fished in V (Fref.v)
  param$ind_min_det <- ((param$min_fishing_size_detritivore-
                           param$min_log10_plankton)/dx)+1
  
  #short hand for matrix indexing (idx)
  param$idx <- 2:param$numb_size_bins
  
  # number of grid cells (Ngrid)
  param$numb_grid_cells <- Ngrid 
  
  # (phyto+zoo)plankton + pelagic predator size spectrum (U.init)
  param$plank_pred_sizes <- 10^param$int_phy_zoo[1]*
    10^(param$slope_phy_zoo[1]*param$log10_size_bins)
  
  # set initial detritivore spectrum (V.init)
  param$detritivore_sizes <- param$sinking_rate[1]*10^param$int_phy_zoo[1] * 
    10^(param$slope_phy_zoo[1]*param$log10_size_bins)
  
  # abritrary initial value for detritus (W.init)
  param$init_detritus <- 0.00001
  
  if(use_init == TRUE){
    #(U.init)
    param$init_pred <- u_initial
    #(V.init)
    param$init_detritivores <- v_initial
    #(W.init)
    param$init_detritus <- w_initial
  }
  
  param$equilibrium <- equilibrium
  
  return(param)
}


