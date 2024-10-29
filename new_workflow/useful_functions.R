###############################################################################
# Supporting DBPM functions
# Functions have been adapted from previous DBPM work
# 
# Edited by: Denisse Fierro Arcos
# Date of update: 2024-10-23


# Loading libraries -------------------------------------------------------
# library(dplyr)


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
  
  # discretisation of year(delta.t)
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
  param$ind_min_pred_size <- which(param$log10_size_bins == param$min_log10_pred)
    #((param$min_log10_pred-param$min_log10_plankton)/dx)+1
  
  # index for minimum log10 detritivore size (ref.det)
  param$ind_min_detritivore_size <- which(param$log10_size_bins == 
                                            param$min_log10_detritivore)
    #((param$min_log10_detritivore-param$min_log10_plankton)/dx)+1 
  
  # index in F vector corresponding to smallest size fished in U (Fref.u)
  param$ind_min_fish_pred <- ((param$min_fishing_size_pred-param$min_log10_plankton)/
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


# Build a lookup table for diet preference ------
# Looks at all combinations of predator and prey body size: diet preference 
# (in the predator spectrum only)
phi_f <- function(q, log10_pred_prey_ratio, log_prey_pref){
  phi <- ifelse(q > 0, 
                exp(-(q-log10_pred_prey_ratio)*(q-log10_pred_prey_ratio)/
                      (2*log_prey_pref*log_prey_pref))/
                  (log_prey_pref*sqrt(2.0*pi)),
                0) 
  return(phi)
}

# Build lookup tables for (constant) growth ------
# Considers components which remain constant
gphi_f <- function(pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref){
  return(10^(-pred_prey_matrix)*phi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                                      log_prey_pref))
}	

# Build lookup tables for (constant) mortality ------
mphi_f <- function(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref,
                   metabolic_req_pred){
  return(10^(metabolic_req_pred*rev_pred_prey_matrix)*
           phi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref))
}	

# Build lookup table for components of 10^(alpha*x) ------
expax_f <- function(log10_size_bins, metabolic_req_pred){
  return(10^(metabolic_req_pred*log10_size_bins)) 
}


# Run model per grid cell or averaged over an area ------
sizemodel <- function(params, ERSEM_det_input = F, temp_effect = T, 
                      use_init = F){
  with(params,{
    #--------------------------------------------------------------------------
    # Model for a dynamical ecosystem comprised of: two functionally distinct 
    # size spectra (predators and detritivores), size structured primary 
    # producers and an unstructured detritus resource pool. 
    # time implicit upwind difference discretization algorithm (from Ken 
    # Andersen)
    # (Numerical Recipes and Richards notes, see "C://...standardised test 
    # implicit.txt" for basic example)
    # U represents abundance density (m-2) of "fish" and V abundance density 
    # (m-2) of "detritivores"
    # W is the pool of detritus (expressed in biomass density, g.m-2) - not 
    # size-based
    # Fish feed on other fish and benthic detritivores with size preference 
    # function 
    # Benthic detritivores feed on detritus produced from pelagic spectrum 
    # Senescence Mortality also included.  Options to include dynamic 
    # reproduction and predator handling time (but not currently used).
    #
    # Code modified for global fishing mortality rate application. 
    # JLB 17/02/2014
    # Code modified to include temperature scaling on senescence and detrital 
    # flux. RFH 18/06/2020
    # -------------------------------------------------------------------------
    
    # Input parameters to vary:
    # time series of intercept of plankton size spectrum (estimated from GCM, 
    # biogeophysical model output or satellite data).		
    ui0 <- 10^int_phy_zoo  
    # time series of slope of plankton size spectrum (estimated from GCM, 
    #biogeophysical model output or satellite data).  	
    slope_phy_zoo <- slope_phy_zoo
    # time series of temperature in water column
    sea_surf_temp <- sea_surf_temp
    #time series of temperature near seabed
    sea_floor_temp <- sea_floor_temp
    #time series of export ratio (read in sizeparam)
    sinking_rate <- sinking_rate
    
    
    #--------------------------------------------------------------------------
    # Initialising matrices
    #--------------------------------------------------------------------------
    
    #q1 is a square matrix holding the log(predatorsize/preysize) for all 
    #combinations of sizes (q1)
    pred_prey_matrix  <- matrix(NA, numb_size_bins, numb_size_bins)
    for(i in 1:numb_size_bins){
      pred_prey_matrix[, i] <- log10_size_bins[i] - log10_size_bins
    }
    
    #q2 is the reverse matrix holding the log(preysize/predatorsize) for all 
    #combinations of sizes (q2)
    rev_pred_prey_matrix <- matrix(-pred_prey_matrix, numb_size_bins, 
                                   numb_size_bins)
    
    #matrix for recording the two size spectra (V - det, U - pred)
    detritivores <- predators <- array(0, c(numb_size_bins, numb_time_steps+1))
    
    #vector to hold detritus biomass density (g.m-2) (W)
    detritus <- array(0, numb_time_steps+1)
    
    #matrix for keeping track of growth (GG_v, GG_u) and reproduction (R_v, R_u) 
    #from ingested food:
    reprod_det <- reprod_pred <- array(0, c(numb_size_bins, numb_time_steps+1))
    growth_det <- growth_pred <- array(0, c(numb_size_bins, numb_time_steps+1))
    
    #matrix for keeping track of predation mortality (PM_v, PM_u)
    pred_mort_det <- array(0, c(numb_size_bins, numb_time_steps+1))
    pred_mort_pred <- array(0, c(numb_size_bins, numb_time_steps+1))
    
    #matrix for keeping track of catches (Y_v, Y_u)
    catch_det <- catch_pred <- array(0, c(numb_size_bins, numb_time_steps+1))  
    
    #matrix for keeping track of total mortality (Z_v, Z_u)
    tot_mort_det <- tot_mort_pred <- array(0, c(numb_size_bins,
                                                numb_time_steps+1))
    
    #matrix for keeping track of senescence mortality (SM_v, SM_u) and other 
    #(intrinsic) mortality (OM_v, OM_u)
    senes_mort_det <- senes_mort_pred <- array(0, numb_size_bins)
    other_mort_det <- other_mort_pred <- array(0, numb_size_bins)
    
    #empty vector to hold fishing mortality rates at each size class at time 
    #(Fvec_v, Fvec_u)
    fishing_mort_det <- fishing_mort_pred <- array(0, c(numb_size_bins, 
                                                        numb_time_steps+1))
    
    #lookup tables for terms in the integrals which remain constant over time
    #(gphi, mphi)
    constant_growth <- gphi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                              log_prey_pref)
    constant_mortality <- mphi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, 
                                 log_prey_pref, metabolic_req_pred)
    
    #lookup table for components of 10^(metabolic_req_pred*log10_size_bins) 
    #(expax)
    met_req_log10_size_bins <- expax_f(log10_size_bins, metabolic_req_pred)
    
    #--------------------------------------------------------------------------
    # Numerical integration
    #--------------------------------------------------------------------------
    
    # set up with the initial values from param
    #(phyto+zoo)plankton size spectrum  
    predators[1:(ind_min_pred_size-1), 1] <- plank_pred_sizes[1:(ind_min_pred_size-1)]    
    # set initial consumer size spectrum 
    predators[ind_min_pred_size:120, 1] <- plank_pred_sizes[ind_min_pred_size:120]       
    # set initial detritivore spectrum  
    detritivores[ind_min_detritivore_size:120, 1] <- detritivore_sizes[ind_min_detritivore_size:120]  
    # set initial detritus biomass density (g.m^-3) 
    detritus[1] <- init_detritus
    
    if(use_init){
      # set up with the initial values from previous run
      #(phyto+zoo)plankton size spectrum  
      predators[1:(ind_min_pred_size-1), 1] <- plank_pred_sizes[1:(ind_min_pred_size-1)]    
      # set initial consumer size spectrum from previous run
      predators[ind_min_pred_size:numb_size_bins, 1] <- plank_pred_sizes[ind_min_pred_size:numb_size_bins]
      # set initial detritivore spectrum from previous run
      detritivores[ind_min_detritivore_size:numb_size_bins, 1] <- detritivore_sizes[ind_min_detritivore_size:numb_size_bins] 
      detritus[1] <- init_detritus
    }
    
    #intrinsic natural mortality
    other_mort_pred <- natural_mort*10^(-0.25*log10_size_bins)
    other_mort_det <- natural_mort*10^(-0.25*log10_size_bins)
    
    #senescence mortality rate to limit large fish from building up in the 
    #system
    #same function as in Law et al 2008, with chosen parameters gives similar 
    #M2 values as in Hall et al. 2006
    senes_mort_pred <- const_senescence_mort*10^(exp_senescence_mort*(log10_size_bins-size_senescence))
    senes_mort_det <- const_senescence_mort*10^(exp_senescence_mort*(log10_size_bins-size_senescence))
    
    #Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!)
    # from Benoit & Rochet 2004 
    # here fish_mort_pred and fish_mort_pred= fixed catchability term for predators and 
    # detritivores to be estimated along with ind_min_det and ind_min_fish_pred
    fishing_mort_pred[ind_min_fish_pred:numb_size_bins, 1] <- fish_mort_pred*effort[1]
    fishing_mort_det[ind_min_det:numb_size_bins, 1] <- fish_mort_detritivore*effort[1]
    
    #output fisheries catches per yr at size
    catch_pred[ind_min_fish_pred:numb_size_bins, 1] <- fishing_mort_pred[ind_min_fish_pred:numb_size_bins, 1]*predators[ind_min_fish_pred:numb_size_bins, 1]*10^log10_size_bins[ind_min_fish_pred:numb_size_bins] 
    #output fisheries catches per yr at size
    catch_det[ind_min_det:numb_size_bins, 1] <- fishing_mort_det[ind_min_det:numb_size_bins, 1]*detritivores[ind_min_det:numb_size_bins, 1]*10^log10_size_bins[ind_min_det:numb_size_bins] 
    
    #iteration over time, N [days]
    for(i in 1:(numb_time_steps)){
      
      #--------------------------------
      # Calculate Growth and Mortality
      #--------------------------------
      if(temp_effect == T){
        pel_tempeffect <- exp(c1-activation_energy/(boltzmann*(sea_surf_temp+273)))
        ben_tempeffect <- exp(c1-activation_energy/(boltzmann*(sea_floor_temp+273)))
      }
      if(temp_effect == F){
        pel_tempeffect <- 1
        ben_tempeffect <- 1
      }
      
      # feeding rates
      # yr-1
      #(f_pel)
      feed_rate_pel <- pel_tempeffect[i] * 
        as.vector(((hr_volume_search*10^(log10_size_bins*metabolic_req_pred)*pref_pelagic)*(predators[, i]*log_size_increase)%*%(constant_growth)) / 
                    (1+handling*(hr_volume_search*10^(log10_size_bins*metabolic_req_pred)*pref_pelagic) * 
                       (predators[, i]*log_size_increase)%*%(constant_growth))) 
      # yr-1
      #(f_ben)
      feed_rate_bent <- pel_tempeffect[i] * 
        as.vector(((hr_volume_search*10^(log10_size_bins*metabolic_req_pred)*pref_benthos)*(detritivores[, i]*log_size_increase)%*%(constant_growth)) / 
                    (1+handling*(hr_volume_search*10^(log10_size_bins*metabolic_req_pred)*pref_benthos) *
                       (detritivores[, i]*log_size_increase)%*%(constant_growth)))
      # yr-1
      #(f_det)
      feed_rate_det <- ben_tempeffect[i] *
        ((1/10^log10_size_bins)*hr_vol_filter_benthos*10^(log10_size_bins*metabolic_req_detritivore)*detritus[i]) / 
        (1+handling*(1/10^log10_size_bins)*hr_vol_filter_benthos*10^(log10_size_bins*metabolic_req_detritivore)*detritus[i])
      
      # Predator growth integral 
      # yr-1
      growth_pred[, i] <- (1-defecate_prop)*growth_pred*(feed_rate_pel)+(1-def_low)*growth_detritivore*(feed_rate_bent)
      
      # Reproduction
      # yr-1
      if(dynamic_reproduction){
        reprod_pred[, i] <- (1-defecate_prop)*(1-(growth_pred+energy_pred))*(feed_rate_pel) +
          (1-defecate_prop)*(1-(growth_detritivore+energy_detritivore))*feed_rate_bent
      }
      
      # Predator death integrals 
      #Satiation level of predator for pelagic prey
      sat_pel <- ifelse(feed_rate_pel > 0,
                        feed_rate_pel/((hr_volume_search*10^(log10_size_bins*metabolic_req_pred)*pref_pelagic) * 
                                 (predators[,i]*log_size_increase)%*%(constant_growth)),
                        0)
      # yr-1
      pred_mort_pred[, i] <- as.vector((pref_pelagic*hr_volume_search*met_req_log10_size_bins)*(predators[, i]*sat_pel*log_size_increase)%*%(constant_mortality))
      
      # yr-1
      tot_mort_pred[, i] <- pred_mort_pred[, i]+pel_tempeffect[i]*other_mort_pred+senes_mort_pred+fishing_mort_pred[, i]
      
      # Benthos growth integral
      # yr-1
      growth_det[, i] <- (1-def_low)*growth_detritus*feed_rate_det
      
      #reproduction
      # yr-1
      if(dynamic_reproduction){
        reprod_det[, i] <- (1-def_low)*(1-(growth_detritus+energy_detritivore))*(feed_rate_det)
      }
      
      # Benthos death integral
      #Satiation level of predator for benthic prey 
      sat_ben <- ifelse(feed_rate_bent > 0,
                        feed_rate_bent/((hr_volume_search*10^(log10_size_bins*metabolic_req_detritivore)*pref_benthos) * 
                                 (detritivores[, i]*log_size_increase)%*%(constant_growth)),
                        0)
      # yr-1
      pred_mort_det[, i] <- ifelse(sat_ben > 0,
                          as.vector((pref_benthos*hr_volume_search*met_req_log10_size_bins) *
                                      (predators[, i]*sat_ben*log_size_increase)%*%(constant_mortality)),
                          0)
      
      # yr-1
      tot_mort_det[, i] <- pred_mort_det[, i]+ben_tempeffect[i]*other_mort_det+senes_mort_det+fishing_mort_det[, i]
      
      #total biomass density eaten by pred (g.m-2.yr-1)
      eatenbypred <- 10^log10_size_bins*feed_rate_pel*predators[, i]+10^log10_size_bins*feed_rate_bent*predators[, i] 
      
      #detritus output (g.m-2.yr-1)
      # losses from detritivore scavenging/filtering only:
      output_w <- sum(10^log10_size_bins*feed_rate_det*detritivores[, i]*log_size_increase)   
      
      #total biomass density defecated by pred (g.m-2.yr-1)
      defbypred <- defecate_prop*(feed_rate_pel)*10^log10_size_bins*predators[, i]+ def_low*(feed_rate_bent)*10^log10_size_bins*predators[, i]
      
      #------------------------------------------------
      # Increment values of detritus, predators & detritivores	for next time step  
      #------------------------------------------------
      
      #Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of 
      #detritus pool and solve for detritus biomass density in next time step 
      if(!ERSEM_det_input){
        #considering pelagic faeces as input as well as dead bodies from both 
        #pelagic and benthic communities and phytodetritus (dying sinking
        #phytoplankton)
        if(detritus_coupling){
          # pelagic spectrum inputs (sinking dead bodies and faeces) - export 
          # ratio used for "sinking rate" + benthic spectrum inputs (dead stuff
          # already on/in seafloor)
          input_w <- (sinking_rate[i] * 
                        (sum(defbypred[ind_min_pred_size:numb_size_bins]*log_size_increase) +
                           sum(pel_tempeffect[i] * 
                                 other_mort_pred[1:numb_size_bins]*predators[1:numb_size_bins, i]*10^(log10_size_bins[1:numb_size_bins])*log_size_increase) + 
                           sum(pel_tempeffect[i] * 
                                 senes_mort_pred[1:numb_size_bins]*predators[1:numb_size_bins,i]*10^(log10_size_bins[1:numb_size_bins])*log_size_increase)) +
                        (sum(ben_tempeffect[i] * 
                               other_mort_det[1:numb_size_bins]*detritivores[1:numb_size_bins,i]*10^(log10_size_bins[1:numb_size_bins])*log_size_increase) + 
                           sum(ben_tempeffect[i] * 
                                 senes_mort_det[1:numb_size_bins]*detritivores[1:numb_size_bins,i]*10^(log10_size_bins[1:numb_size_bins])*log_size_increase)))
          # )
        }
        if(!detritus_coupling){
          input_w <- sum(ben_tempeffect[i] *
                           other_mort_det[1:numb_size_bins]*detritivores[1:numb_size_bins, i]*10^(log10_size_bins[1:numb_size_bins])*log_size_increase) + 
            sum(ben_tempeffect[i]*senes_mort_det[1:numb_size_bins]*detritivores[1:numb_size_bins, i]*10^(log10_size_bins[1:numb_size_bins])*log_size_increase)
        }
        
        # get burial rate from Dunne et al. 2007 equation 3
        burial <- input_w*(0.013 + 0.53*input_w^2/(7+input_w)^2)
        
        # losses from detritivory + burial rate (not including remineralisation
        # bc that goes to p.p. after sediment, we are using realised p.p. as
        # inputs to the model) 
        dW <- input_w-(output_w+burial) 
        #biomass density of detritus g.m-2
        detritus[i+1] <- detritus[i]+dW*timesteps_years
      }
      
      if(ERSEM_det_input){
        detritus[i+1] <- detritus[i]
      }
      
      #----------------------------------------------
      # Pelagic Predator Density (nos.m-2)- solve for time + timesteps_years using
      # implicit time Euler upwind finite difference (help from Ken Andersen 
      # and Richard Law)
      
      # Matrix setup for implicit differencing 
      Ai_u <- Bi_u <- Si_u <- array(0, c(numb_size_bins, 1))   
      
      Ai_u[idx] <- (1/log(10))*-growth_pred[idx-1, i]*timesteps_years/log_size_increase
      Bi_u[idx] <- 1+(1/log(10))*growth_pred[idx, i]*timesteps_years/log_size_increase +tot_mort_pred[idx, i]*timesteps_years
      Si_u[idx] <- predators[idx, i]
      
      # Boundary condition at upstream end 
      Ai_u[ind_min_pred_size] <- 0
      Bi_u[ind_min_pred_size] <- 1
      Si_u[ind_min_pred_size] <- predators[ind_min_pred_size, i]
      
      # Invert matrix
      #recruitment at smallest consumer mass
      #continuation of plankton hold constant  
      if(use_init){
        predators[1:ind_min_pred_size, i+1] <- ui0[i]*10^(slope_phy_zoo[i]*log10_size_bins)[1:(ind_min_pred_size)] 
      }
      if(!use_init){
        predators[1:ind_min_pred_size, i+1] <- ui0[i+1]*10^(slope_phy_zoo[i+1]*log10_size_bins)[1:(ind_min_pred_size)] 
      }
      
      # apply transfer efficiency of 10% *plankton density at same size
      # reproduction from energy allocation
      if(dynamic_reproduction){
        predators[ind_min_pred_size, i+1] <- predators[ind_min_pred_size, i] +
          (sum(reprod_pred[(ind_min_pred_size+1):numb_size_bins, i]*10^log10_size_bins[(ind_min_pred_size+1):numb_size_bins]*predators[(ind_min_pred_size+1):numb_size_bins, i]*log_size_increase) *
             timesteps_years)/(log_size_increase*10^log10_size_bins[ind_min_pred_size])-(timesteps_years/log_size_increase)*(1/log(10)) *
          (growth_pred[ind_min_pred_size, i])*predators[ind_min_pred_size, i]-timesteps_years*tot_mort_pred[ind_min_pred_size, i]*predators[ind_min_pred_size, i]
      }
      
      #main loop calculation
      for(j in (ind_min_pred_size+1):(numb_size_bins)){
        predators[j, i+1] <- (Si_u[j]-Ai_u[j]*predators[j-1, i+1])/Bi_u[j]
      }
      
      #----------------------------------------------
      # Benthic Detritivore Density (nos.m-2) 
      Ai_v <- Bi_v <- Si_v <- array(0, c(numb_size_bins, 1))   
      #shorthand for matrix referencing
      idx <- (ind_min_detritivore_size+1):numb_size_bins  
      
      Ai_v[idx] <- (1/log(10))*-growth_det[idx-1,i]*timesteps_years/log_size_increase 
      Bi_v[idx] <- 1+(1/log(10))*growth_det[idx, i]*timesteps_years/log_size_increase +tot_mort_det[idx, i]*timesteps_years
      Si_v[idx] <- detritivores[idx, i]
      
      #boundary condition at upstream end
      Ai_v[ind_min_detritivore_size] <- 0
      Bi_v[ind_min_detritivore_size] <- 1
      Si_v[ind_min_detritivore_size] <- detritivores[ind_min_detritivore_size, i]  
      
      #invert matrix
      #recruitment at smallest detritivore mass  
      #hold constant continuation of plankton with sinking rate multiplier 
      detritivores[1:ind_min_detritivore_size, i+1] <- detritivores[1:ind_min_detritivore_size, i]
      
      # apply a very low of transfer efficiency 1%* total biomass of detritus
      #divided by minimum size
      if(dynamic_reproduction){
        detritivores[ind_min_detritivore_size, i+1] <- detritivores[ind_min_detritivore_size, i] +
          sum(reprod_det[(ind_min_detritivore_size+1):numb_size_bins, i]*10^log10_size_bins[(ind_min_detritivore_size+1):numb_size_bins] * 
                detritivores[(ind_min_detritivore_size+1):numb_size_bins, i]*log_size_increase)*timesteps_years/(log_size_increase*10^log10_size_bins[ind_min_detritivore_size]) - 
          (timesteps_years/log_size_increase)*(1/log(10))*(growth_det[ind_min_detritivore_size, i])*detritivores[ind_min_detritivore_size, i] - 
          timesteps_years*tot_mort_det[ind_min_detritivore_size, i]*detritivores[ind_min_detritivore_size, i]
      }
      
      #loop calculation
      for(j in (ind_min_detritivore_size+1):(numb_size_bins)){ 
        detritivores[j, i+1] <- (Si_v[j]-Ai_v[j]*detritivores[j-1, i+1])/Bi_v[j]
      }		
      rm(j)
      
      # increment fishing 
      fishing_mort_pred[ind_min_fish_pred:numb_size_bins, i+1] <- fish_mort_pred*effort[i+1]
      fishing_mort_det[ind_min_det:numb_size_bins, i+1] <- fish_mort_detritivore*effort[i+1]
      
      #output fisheries catches per yr at size
      catch_pred[ind_min_fish_pred:numb_size_bins, i+1] <- fishing_mort_pred[ind_min_fish_pred:numb_size_bins, i+1]*predators[ind_min_fish_pred:numb_size_bins, i+1] *
        10^log10_size_bins[ind_min_fish_pred:numb_size_bins] 
      #output fisheries catches per yr at size
      catch_det[ind_min_det:numb_size_bins, i+1] <- fishing_mort_det[ind_min_det:numb_size_bins, i+1]*detritivores[ind_min_det:numb_size_bins, i+1] *
        10^log10_size_bins[ind_min_det:numb_size_bins] 
    }
    #end time iteration
    
    return(list(predators = predators[,],
                growth_pred = growth_pred[,],
                pred_mort_pred = pred_mort_pred[,],
                detritivores = detritivores[,], 
                growth_det = growth_det[,],
                pred_mort_det = pred_mort_det[,],
                catch_pred = catch_pred[,],
                catch_det = catch_det[,],
                detritus = detritus[],
                params = params))
  })  
  # end with(params)
}
#end size-based model function


