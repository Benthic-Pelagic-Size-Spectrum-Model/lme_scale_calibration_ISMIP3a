###############################################################################
# Library of functions developed originally by JLB on 21/12/2022
#
# Functions were updated by Denisse Fierro Arcos so they can be used with data
# produced by updated `01_getinputs_ISIMIP3a.R` script.
#
# Date of update: 2024-08-06

# Loading libraries
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggplot2)
library(lubridate)
library(zoo)
library(lhs)
library(patchwork)
library(purrr)
library(janitor)
library(parallel)
# library(optimParallel)
source("dbpm_model_functions_NEW_DFA.R")
# library(gridExtra)


# Reading LME-scale time-series inputs (climate and fishing) -----
get_lme_inputs <- function(forcing_file = NULL, gridded_forcing = NULL, 
                           fishing_effort_file, LMEnumber, yearly = F,
                           file_out = NULL){
  #This function reads pre-processed inputs as spatially averaged means to 
  #estimate catchability (and if needed other model parameters)
  #Inputs:
  #forcing_file (character) - Full path to forcing file. This must be 
  #non-gridded data
  #gridded_forcing (character) - Full path to folder containing gridded forcing
  #files
  #fishing_effort_file (character) - Full path to fishing effort file
  #LMEnumber (numeric) - Unique ID identifying an LME
  #yearly (boolean) - Default is FALSE. If set to TRUE, it will return yearly
  #means for all forcing variables
  #file_out (character) - Optional. If provided, the function does not return
  #anything, but gridded data with stable spinup is returned.
  #
  #Output:
  #lme_clim (data frame) - Forcing data to be used in model calibration
  
  # Climate fishing inputs available via THREDDS:
  # http://portal.sf.utas.edu.au/thredds/catalog/gem/fishmip/ISIMIP3a/
  # InputData/DBPM_lme_inputs/obsclim/025deg/catalog.html?dataset=
  # fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/025deg/
  # DBPM_LME_effort_catch_input.csv
  
  lme_fish <- fread(fishing_effort_file) |> 
    #Subset data for LME of interest
    filter(region == LMEnumber) |> 
    mutate(nom_active_relative = total_nom_active/max(total_nom_active),
           nom_active_area_m2_relative = total_nom_active_area_m2/
                           max(total_nom_active_area_m2))
  
  if(!is.null(forcing_file)){
    # Climate forcing inputs available via THREDDS:
    # http://portal.sf.utas.edu.au/thredds/catalog/gem/fishmip/ISIMIP3a/
    # InputData/DBPM_lme_inputs/obsclim/025deg/catalog.html?dataset=
    # fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/025deg/
    # DBPM_LME_climate_inputs_slope.csv
    
    lme_clim <- fread(forcing_file) |> 
      filter(region == LMEnumber)
    
    if(yearly){
      # option 1: keep monthly time steps but monthly value is the yearly mean
      
      # Calculating yearly average for all variables
      lme_clim <- lme_clim |> 
        mutate(year = year(t)) |> 
        group_by(year, region) |> 
        summarise(across(sst:expcbot, \(x) mean(x, na.rm = T))) |> 
        ungroup()
      
      # option 2: do as per global DBPML (builds on option 1): 
      # uses yearly average and extend the time series to monthly inputs
      # interpolate missing values using na.approx
      lme_clim <- lme_clim |> 
        #Yearly mean is copied 12 times (for each month in a year)
        slice(rep(1:n(), times = 12)) |> 
        #Adding month as column
        mutate(month = rep(1:12, each = nrow(lme_clim)),
               #Create a date column from year and month columns
               t = ym(paste(year, month, sep = "-"))) |> 
        select(!month) |> 
        arrange(t) |> 
        #Interpolate any missing values
        mutate(across(c(sst:lphy, expcbot), ~ na.approx(.x))) |> 
        mutate(year = as.numeric(year))
    }
    
    #Adding effort data
    lme_clim <- lme_clim |>  
      mutate(year = as.numeric(year(t))) |> 
      left_join(lme_fish, by = c("year", "region"))
  }
  
  else if(!is.null(gridded_forcing)){
    lme_clim <- list.files(gridded_forcing, 
                           paste0("LME_", LMEnumber, "_"), full.names = T) |> 
      fread() |> 
      #Extracting digits identifying region
      mutate(region = as.numeric(str_extract(region, "\\d{1,2}")))
   
    if(yearly){
      # option 1: keep monthly time steps but monthly value is the yearly mean
      # Calculating yearly average for all variables
      lme_clim <- lme_clim |> 
        mutate(year = year(t)) |> 
        group_by(year, region, lat, lon) |> 
        summarise(across(sst:expcbot, \(x) mean(x, na.rm = T))) |> 
        ungroup()
        
      # option 2: do as per global DBPM (builds on option 1): 
      # uses yearly average and extend the time series to monthly inputs
      # interpolate missing values using na.approx
      lme_clim <- lme_clim |> 
        #Yearly mean is copied 12 times (for each month in a year)
        slice(rep(1:n(), times = 12)) |> 
        #Adding month as column
        mutate(month = rep(1:12, each = nrow(lme_clim)),
               #Create a date column from year and month columns
               t = ym(paste(year, month, sep = "-"))) |> 
        select(!c(month)) |> 
        arrange(lat, lon, t) |> 
        #Interpolate any missing values
        mutate(across(c(sst:lphy, expcbot), ~ na.approx(.x)))
    }
    
    #Adding effort data
    lme_fish_sub <- lme_fish |>  
      select(!deptho) |>  
      distinct()
    
    lme_clim <- lme_clim |> 
      left_join(lme_fish_sub, by = c("year", "region")) |> 
      arrange(lat, lon)
  }
  
  #Add a spin-up to these inputs prior to 1841 - 100 yrs at first value
  if(!is.null(file_out)){
    if(!dir.exists(dirname(file_out))){
      dir.create(dirname(file_out), recursive = T)
    }
    lme_clim |> 
      fwrite(file_out)
  }
  
  return(lme_clim)
}


# Running model with time series ----
run_model <- function(vals = X, input = lme_input, withinput = T){
  #Inputs:
  #vals (named numeric vector) - Single column with named rows containing LHS
  #parameters
  #input (data frame) - Forcing data produced by `get_lme_inputs` function
  #withinput (boolean) - Default is TRUE. ????
  #
  #Output:
  #If withinput set to TRUE:
  #input (data frame) - ????
  #If withinput set to FALSE:
  #result_set (data frame) - ???
  
  #Ensuring values are numeric
  if(!is.numeric(vals)){
    stop("Values provided in the 'vals' parameter must be numeric.")
  }
  
  # set up parameters
  params <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin.consumer.u = -3,
                      xmin.consumer.v = -3, tmax = nrow(input)/12,
                      tstepspryr = 12, search_vol = vals["search.vol"],
                      fmort.u = vals["f.u"], fminx.u = vals["f.minu"], 
                      fmort.v = vals["f.v"], fminx.v = vals["f.minv"], 
                      depth = mean(input$depth), er = input$er,
                      pp = input$intercept, slope = input$slope, 
                      sst = input$sst, sft = input$sbt, use.init = FALSE,
                      effort = input$nom_active_area_m2_relative)
  
  # run model through time
  # TO DO IN SIZEMODEL CODE: make fishing function like one in model template
  result_set <- sizemodel(params)
  
  if(withinput){
    lims_ubio <- result_set$params$ref:result_set$params$Nx
    # JB:  changed inputs to m2 so no need to divide by depth here
    # added 1:params$Neq to same 2040 time steps instead of 2041
    time_steps <- 1:result_set$params$Neq
    
    input$TotalUbiomass <- apply(result_set$U[lims_ubio, time_steps] * 
                                   params$dx*10^params$x[lims_ubio], 
                                 2, sum)#*min(params$depth,100)
    
    lims_vbio <- result_set$params$ref.det:result_set$params$Nx
    input$TotalVbiomass <- apply(result_set$V[lims_vbio, time_steps] * 
                                   params$dx*10^params$x[lims_vbio],
                                 2, sum)#*min(params$depth,100)
    
    input$TotalW <- result_set$W[time_steps]#*min(params$depth,100)
    
    #sum catches (currently in grams per m3 per year, across size classes) 
    # keep as grams per m2, then be sure to convert observed from tonnes per m2
    # per year to g.^-m2.^-yr (for each month)
    
    input$TotalUcatch <- apply(result_set$Y_u[, time_steps]*params$dx, 2, sum)
    #*min(params$depth,100)
    input$TotalVcatch <- apply(result_set$Y_v[, time_steps]*params$dx, 2, sum)
    #*min(params$depth,100)
    input$Totalcatch <- input$TotalUcatch + input$TotalVcatch
  
   return(input)
    
  }else{
    return(result_set)
  }
}


# Comparing observed and predicted fish biomass ----
getError <- function(lhs_params = X, lme_forcings = lme_input, 
                     corr = F, figure_folder = NULL){
  #Inputs:
  #lhs_params (data frame) - Contains LHS parameters as columns and it must 
  #have a single row
  #lme_forcings (data frame) - Forcing data produced by `get_lme_inputs` 
  #function
  #corr (boolean) - Default is FALSE. If set to TRUE, it will calculate the
  #correlation between predicted and observed values
  #figure_folder (character) - Optional. If provided, it must be the full path
  #to the folder where figures comparing observed and predicted data will be 
  #stored
  #
  #Output:
  #rmse (numeric) - RMSE value between observed and predicted catch
  
  #Running model for specific LME
  result <- run_model(lhs_params, lme_forcings)
  
  #Aggregate data by year (mean to conserve units)
  error_calc <- result |> 
    filter(year > 1949) |> 
    group_by(year) |> 
    summarise(mean_total_catch_yr = mean(Totalcatch),
              mean_obs_catch_yr = mean(catch_tonnes_area_m2, na.rm = T)) |> 
    #Converting units from tonnes to g
    mutate(mean_obs_catch_yr = mean_obs_catch_yr*1e6,
           #calculate and output error 
           #convert from tonnes to grams (m^-2*yr^-1)
           squared_error = (mean_obs_catch_yr-mean_total_catch_yr)^2)
  
  #Calculate RMSE 
  sum_se <- sum(error_calc$squared_error, na.rm = T)
  count <- sum(!is.na(error_calc$squared_error))
  rmse <- sqrt(sum_se/count)
  
  if(corr){
    #Calculate correlation between observed and predicted catches
    corr_nas <- data.frame(cor = cor(error_calc$mean_obs_catch_yr, 
                                     error_calc$mean_total_catch_yr, 
                                     use = "complete.obs"),
                           #Get number of rows with NA values
                           catchNA = sum(is.na(error_calc$mean_total_catch_yr)),
                           region = unique(lme_forcings$region))
  }
  
  #If a path to save figures is provided, create figures and save 
  if(!is.null(figure_folder)){
    if(!dir.exists(figure_folder)){
      dir.create(figure_folder)
    }
    #Plotting predicted and observed catches over time
    p1 <- ggplot()+
      geom_line(data = error_calc, aes(x = year, y = mean_total_catch_yr))+
      geom_point(data = error_calc, aes(x = year, y = mean_obs_catch_yr))+
      scale_x_continuous(breaks = seq(min(error_calc$year), 
                                      max(error_calc$year), by = 10))+
      theme_classic()+ 
      theme(axis.text = element_text(colour = "grey20", size = 12),
            text = element_text(size = 15))+
      labs(x = "Year", y = bquote("Total catch (g*"~yr^-1*"*"*m^-2*")"))
    
    #Plotting predicted vs observed catches
    p2 <- ggplot()+
      geom_point(data = error_calc, 
                 aes(x = mean_total_catch_yr, y = mean_obs_catch_yr))+
      geom_abline(slope = 1, intercept = 0)+
      theme_classic()+
      theme(axis.text = element_text(colour = "grey20", size = 12),
            text = element_text(size = 15))+
      labs(x = "Predicted", y = "Observed")
    
    #Creating a single plot
    p3 <- p1+p2+
      plot_annotation(title = paste0("LME #", unique(lme_forcings$region)),
                      theme = theme(plot.title = element_text(size = 16, 
                                                              hjust = 0.5)))
    
    #Creating path to save figure
    f_out <- file.path(figure_folder,
                       paste0("dbpm_pred_obs_catches_LME_", 
                              unique(lme_forcings$region), ".png"))
    
    #Saving composite figure
    ggsave(filename = f_out, plot = p3, width = 15, height = 10)
  }
  
  if(!corr){
    #Return RMSE
    return(rmse)
  }
  if(corr){
    return(corr_nas)
  }
}


# Tuning LHS parameters -----
LHSsearch <- function(LMEnum = LME, num_iter = 1, search_vol = "estimated",
                      forcing_file = NULL, gridded_forcing = NULL, 
                      fishing_effort_file, corr = F, figure_folder = NULL,
                      best_val_folder = NULL){
  #Inputs:
  #LMEnum (numeric) - Unique ID identifying LME
  #num_iter (numeric) - Number of individual runs. Default is 1.
  #search_vol (???) - Default is "estimated". ???
  #forcing_file (character) - Full path to forcing file. This must be 
  #non-gridded data
  #gridded_forcing (character) - Full path to folder containing gridded forcing
  #files
  #fishing_effort_file (character) - Full path to fishing effort file
  #corr (boolean) - Default is FALSE. If set to TRUE, it will calculate the
  #correlation between predicted and observed values
  #figure_folder (character) - Optional. If provided, it must be the full path
  #to the folder where figures comparing observed and predicted data will be 
  #stored
  #best_val_folder (character) - Optional. If provided, it muste be the full
  #path to the folder where LHS search results will be saved
  #
  #Output:
  #bestvals (data frame) - Contains the values for LHS parameters that resulted
  #in the best performing model based on RMSE values
  
  #Making function reproducible
  set.seed(1234)
  
  #Construct a hypercube with random numbers
  #num_iter defines number of rows in hypercube
  #columns represent five specific parameters needed
  sim <- data.frame(randomLHS(num_iter, 5))
  colnames(sim) <- c("f.u", "f.v", "f.minu", "f.minv", "search.vol")
  
  #Adjust range of mi size params, others go from 0-1
  sim <- sim |> 
    mutate(f.minu = f.minu*2, 
           f.minv = f.minv*2,
           # adjust range of search vol, others go from 0-1
           search.vol = search.vol+0.001)
  
  if(is.numeric(search_vol)){
    sim <- sim |> 
      mutate(search.vol = search_vol)
  }
  
  # use below to select a constant value for search.vol
  if(!is.null(forcing_file)){
    lme_input <- get_lme_inputs(forcing_file = forcing_file, 
                                fishing_effort_file = fishing_effort_file, 
                                LMEnumber = LMEnum)
  }
  if(!is.null(gridded_forcing)){
    lme_input <- get_lme_inputs(gridded_forcing = gridded_forcing, 
                                fishing_effort_file = fishing_effort_file, 
                                LMEnumber = LMEnum)
  }
  
  # parallelise using 75% of cores available using mclapply
  no_cores <- round((detectCores()*.75), 0)
  sim$rmse <- mclapply(1:nrow(sim), 
                       FUN = function(i) getError(unlist(sim[i,]),
                                                  lme_forcings = lme_input, 
                                                  corr, figure_folder), 
                       mc.cores = no_cores) |> 
    unlist()
  
  # check this time param set with lowest error
  bestvals <- sim |> 
    filter(rmse == min(rmse, na.rm = T)) |> 
    mutate(region = LMEnum)
  
  #Print row with lowest RMSE
  print(bestvals)
  
  #If folder to save values is provided - Save results
  if(!is.null(best_val_folder)){
    #Ensure folder exists
    if(!dir.exists(best_val_folder)){
      dir.create(best_val_folder, recursive = T)
    }
    
    #File path to save output
    fout <- file.path(best_val_folder, 
                      paste0("best-fishing-parameters_LME_", LMEnum,
                             "_searchvol_", search_vol, "_numb-iter_", 
                             num_iter, ".csv"))
    #Save output
    bestvals |> 
      fwrite(fout)
  }
  
  return(bestvals)
}


# Correlation and calibration plots ----
corr_calib_plots <- function(fishing_params, forcing_file, 
                              fishing_effort_file, figure_folder = NULL){
  #Inputs:
  #fishing_params (named numeric vector) - Single column with named rows 
  #containing LHS parameters
  #forcing_file (character) - Full path to forcing file. This must be 
  #non-gridded data
  #fishing_effort_file (character) - Full path to fishing effort file
  #figure_folder (character) - Optional. If provided, it must be the full path
  #to the folder where figures comparing observed and predicted data will be 
  #stored
  #
  #Output:
  #corr_nas (data.frame) - Contains the correlation between predicted and
  #observed values
  
  #Get forcing data
  lme_input <- get_lme_inputs(forcing_file = forcing_file, 
                              fishing_effort_file = fishing_effort_file, 
                              LMEnumber = fishing_params["region"])
  
  #Calculate correlations with tuned fishing parameters and save plots
  corr_nas <- getError(lhs_params = fishing_params, lme_forcings = lme_input, 
                       corr = T, figure_folder)
  
  return(corr_nas)
}
  

# Calculating stable spinup -----------------------------------------------
stable_spinup <- function(df, name_col, base_out = NULL){
  # This function calculates input variables for a stable spinup between 1741 
  # and 1840
  #
  # Inputs:
  # df (data frame) - Contains input variables from 1841
  # name_col (character) - Name of the variable for which to calculate spinup
  # base_out (character) - Optional. If provided, the function does not return
  # anything, but gridded data with stable spinup is saved using the path to
  # the folder as a base for naming outputs
  #
  # Outputs:
  # df_out (data frame) - Contains the input variable original data plus the
  # spinup period
  
  #If variable of interest is depth, simply return values for first time step
  if(name_col == "depth"){
    df_out <- df |> 
      filter(t == min(t)) |> 
      select(all_of(c("lat", "lon", "t", name_col))) |> 
      pivot_wider(id_cols = lat:lon , names_from = t, 
                  values_from = name_col) |> 
      arrange(lat, lon)
  }else{
    #Creating a 100 years stable spinup before 1841
    dates <- seq(as.Date("1741-01-01"), as.Date("1840-12-01"), by = "month")
    
    #Select variables of interest
    sub_df <- df |>
      select(all_of(c("lat", "lon", "t", name_col))) |> 
      rename("value" = name_col) |> 
      #Ensure date has date format
      mutate(t = date(t))
    
    #Calculate mean per grid cell for the first year
    spinup <- sub_df |> 
      filter(year(t) == min(year(t))) |> 
      group_by(lat, lon) |> 
      summarise(value = mean(value)) |> 
      ungroup()
    
    #Repeat 1200 months (100 years)
    df_out <- spinup |> 
      slice(rep(1:n(), times = length(dates))) |> 
      mutate(t = rep(dates, each = nrow(spinup)), .before = value) |> 
      #Add original data
      bind_rows(sub_df) |> 
      arrange(t, lat, lon) |> 
      #Rearrange data
      pivot_wider(id_cols = lat:lon, names_from = t, values_from = value)
  }
  
  #If base_out is provided, save data
  if(!is.null(base_out)){
    if(!dir.exists(dirname(base_out))){
      dir.create(dirname(base_out), recursive = T)
    }
    f_out <- str_replace(base_out, ".csv", 
                         paste0("_", name_col, "_1741-2010.csv"))
    df_out |> 
      mutate(var_name = name_col, .after = "lon") |> 
      fwrite(f_out)
  }else{
    #If base_out is not given, then return value
    df_out <- df_out |> 
      select(!lat:lon) |>
      data.matrix()
    return(df_out)
  }
}


# Preparing gridded inputs ------------------------------------------------
# This function get gridded input variables for the region of interest and 
# calculates a stable spinup between 1741 and 1840 
gridded_stable_spinup <- function(meta, gridded_forcing, fishing_effort_file){
  #Inputs:
  #meta (named numeric vector) - Single column with named rows containing 
  #region ID and file paths to save outputs
  #gridded_forcing (character) - Full path to folder containing gridded forcing
  #files
  #fishing_effort_file (character) - Full path to fishing effort file
  #
  #Output:
  #If file_out is set to NULL:
  #lme_clim (data frame) - Forcing data to be used in model calibration from
  #1741 to 2010
  #
  #If file_out is not NULL:
  #Nothing is returned, but lme_clim (data frame) is saved to file_out
  
  # get gridded inputs and run through all grid cells one timestep at a time
  lme_inputs_grid <- get_lme_inputs(gridded_forcing = gridded_forcing, 
                                    fishing_effort_file = fishing_effort_file, 
                                    LMEnumber = meta["region"], 
                                    file_out = meta["grid"])
  
  # List of variables to be processed
  vars_int <- c("sst", "sbt", "er", "intercept", "slope", "depth", 
                "nom_active_area_m2_relative") 
  
  # Create base file path to be used to save variable
  if(!is.na(meta["grid"])){
    gridded_spinup <- vars_int |> 
      map(~stable_spinup(lme_inputs_grid, ., meta["grid"]))
  }else{
    #Otherwise return variables in a named list
    gridded_spinup <- vars_int |> 
      map(~stable_spinup(lme_inputs_grid, .))
    
    names(gridded_spinup) <- vars_int
    
    return(gridded_spinup) 
  }
}


# Merging gridded input ---------------------------------------------------
merging_gridded_inputs <- function(LME_path_full){
  #Inputs:
  #LME_path_full (character) - Full path to folder containing gridded data
  #
  #Output:
  #lme_inputs_grid (list) - Gridded forcing data to be used in model with named
  #items identifying the input variables contained in each list item
    
  #Getting name of files to be used as gridded inputs
  files_in <- list.files(LME_path_full, pattern = "clim_",
                         full.names = T) |> 
    str_subset("gridded.csv", negate = T) 
  
  gridded_files_no_depth <- str_subset(files_in, "depth", negate = T) |> 
    #Read files, but ignore coordinates
    map(~fread(., drop = c("lat", "lon"))) |>
    bind_rows() |> 
    #Group by variable name
    group_by(var_name)
  
  lme_inputs_grid <- gridded_files_no_depth |> 
    group_split(.keep = F) |> 
    #Turn to data matrix
    map(~data.matrix(.))
  
  #Add names to list elements
  names(lme_inputs_grid) <- group_keys(gridded_files_no_depth) |> 
    pull(var_name)
  
  lme_inputs_grid$depth <- str_subset(files_in, "depth") |> 
    fread(drop = c("lat", "lon", "var_name")) |> 
    data.matrix()
  
  return(lme_inputs_grid)
}


# Calculating fishing parameters for gridded model ------------------------
calc_grid_params <- function(meta, f.effort = NULL, start_cond = NULL){
  #Inputs:
  #meta (named character vector) - Single column with named rows: region 
  #(unique ID for LME), base_out (path to folder where outputs will be saved),
  #grid (full file path including filename to be used to save gridded inputs),
  #non_grid (full file path including filename to be used to save non-gridded 
  #inputs)
  #f.effort (named numeric vector) - Optional. If inputs are provided, fisheries
  #effort will be set using bestvals_LME. Otherwise, fisheries effort is set 
  #to 0.
  #start_cond (numeric vector) - Optional. Minimum and maximum years to be 
  #used in calculating steady_state for the model.
  #
  #Output:
  #This function does not return anything. Instead the output is saved in the 
  #base folder provided in the meta parameter
  
  LMEnumber <- meta["region"]
  #Ensuring values are numeric
  if(!is.numeric(f.effort)){
    stop("Values provided in the 'f.effort' parameter must be numeric.")
  }
  
  LME_path_full <- meta["base_out"]
  if(!dir.exists(LME_path_full)){
    dir.create(LME_path_full)
  }
  
  # get initial values from LME-scale results
  lme_input_init <- fread(meta["non_grid"])
  
  #Run initial model with best parameters
  initial_results <- run_model(vals = f.effort, input = lme_input_init, 
                               withinput = F)
  
  #Calculate starting conditions
  #Check if year range was given
  if(!is.null(start_cond)){
    ind <- which(lme_input_init$year >= min(start_cond) & 
                   lme_input_init$year <= max(start_cond))
  }else{
    ind <- which(lme_input_init$year >= 1861 & lme_input_init$year <= 1960)
  }
  
  #Add parameters for initial and end year to calculate starting conditions
  U.initial <- rowMeans(initial_results$U[, ind])
  V.initial <- rowMeans(initial_results$V[, ind])
  W.initial <- mean(initial_results$W[ind])
  
  #Getting name of files to be used as gridded inputs
  lme_inputs_grid <- merging_gridded_inputs(LME_path_full) 
  
  if(!is.null(f.effort)){
    f.u <- f.effort["f.u"]
    f.v <- f.effort["f.v"]
    f.minu <- f.effort["f.minu"]
    f.minv <- f.effort["f.minv"]
  }else{ 
    f.u <- f.v <- f.minu <- f.minv <- 0 
  }
  
  #Get total number of years
  tmax <- length(unique(year(ymd(colnames(lme_inputs_grid$er)))))
  
  # set up params for each month, across grid cells
  gridded_params <- sizeparam(equilibrium = FALSE, dx = 0.1, 
                              xmin.consumer.u = -3, xmin.consumer.v = -3,
                              tmax = tmax, tstepspryr = 12,
                              search_vol = f.effort["search.vol"], 
                              fmort.u = f.u, fminx.u = f.minu, fmort.v = f.v, 
                              fminx.v = f.minv, depth = lme_inputs_grid$depth, 
                              er = lme_inputs_grid$er, 
                              pp = lme_inputs_grid$intercept, 
                              slope = lme_inputs_grid$slope, 
                              sst = lme_inputs_grid$sst, 
                              sft = lme_inputs_grid$sbt, use.init = TRUE,
                              effort = lme_inputs_grid$nom_active_area_m2_relative, 
                              U.initial = U.initial, V.initial = V.initial, 
                              W.initial = W.initial, 
                              Ngrid = nrow(lme_inputs_grid$depth))
  
  # save inputs and parameters object needed for plotting 
  save(gridded_params, 
       file = file.path(LME_path_full, 
                        paste0("grid_inputs_params_LME_", LMEnumber, ".RData")))
}


# Run gridded model by LME ----
rungridbyLME <- function(meta){
  #Run model across space and time using gridded inputs for each LME
  #Inputs:
  #meta (named character vector) - Single column with named rows: region 
  #(unique ID for LME), base_out (path to folder where outputs will be saved),
  #grid (full file path including filename to be used to save gridded inputs),
  #non_grid (full file path including filename to be used to save non-gridded 
  #inputs)
  
  load(list.files(meta["base_out"], "grid_inputs_params_", full.names = T))
  
  lme_input_init <- fread(meta["non_grid"])
  grid_results <- vector("list", nrow(lme_input_init))
  # run model  for full time period across all grid cells
  grid_results <- gridded_sizemodel(gridded_params, ERSEM.det.input = F,
                                    temp_effect = T, eps = 1e-5, 
                                    output = "aggregated", use_init = TRUE)
  
  # removing the stable spinup section to match dimensions with the code
  ind <- which(ymd(colnames(gridded_params$er)) >= min(lme_input_init$t))
  
  grid_results$U <- grid_results$U[, , ind]
  grid_results$V <- grid_results$V[, , ind]
  grid_results$Y.u <- grid_results$Y.u[, , ind]
  grid_results$Y.v <- grid_results$Y.v[, , ind]
  
  # save results from run
  save(grid_results, 
       file = file.path(meta["base_out"], 
                        paste0("gridded_model_results_LME_", 
                               meta["region"], ".RData")))
}


# Running model with gridded inputs ----
run_model_timestep <- function(input = lme_inputs_igrid, 
                               vals = unlist(bestvals_LMEs[14, ]), U.initial,
                               V.initial, W.initial){
  
  # set up params for each month, across grid cells
  params <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin.consumer.u = -3,
                      xmin.consumer.v = -3, tmax = 1/12, tstepspryr = 12,
                      search_vol = 0.64, fmort.u = vals["f.u"], 
                      fminx.u = vals["f.minu"], fmort.v = vals["f.v"], 
                      fminx.v = vals["f.minv"],  depth = input["depth"],
                      er = input["er"], pp = input["intercept"], 
                      slope = input["slope"], sst = input["sst"], 
                      sft = input["sbt"], use.init = TRUE,
                      effort = input["total_nom_active_area_m2"], 
                      U.initial = U.initial, V.initial = V.initial,
                      W.initial = W.initial)
  
  # run model for 1 timestep
  results <- sizemodel(params, use.init = T, ERSEM.det.input = F) 
  
  #sum catches (currently in grams per m3 per year, across size classes) 
  # keep as grams per m2, then be sure to convert observed from tonnes per m2 
  # per year to g.^-m2.^-yr (for each month)
  
  return(cbind(U = results$U[, 2], V = results$V[, 2], Y.u = results$Y.u[, 2],
               Y.v = results$Y.v[, 2], GG.u = results$GG.u[, 1],
               GG.v = results$GG.v[, 1], PM.u = results$PM.u[, 1],
               PM.v = results$PM.v[, 1], W = results$W[2]))
  # return next time step value for that grid
  }


# Adding outputs to gridded_sizemodel ----
getGriddedOutputs <- function(input = lme_inputs_grid, results = grid_results,
                              params = params){
  # returns all outputs of the model 
  input$TotalUbiomass <- input$TotalVbiomass <- input$TotalUcatch <-
    input$TotalVcatch<- input$Totalcatch <- NA
  # Merging coordinates together to use as unique identifier
  input <- input |>
    unite("cell", lat, lon, remove = F)
  cells <- unique(input$cell)
  
  # 2:(params$Neq+1) changed to 1:params$Neq because the newly processed data
  # does subsets data based on time ranges. It does not include December data
  # for year prior
  for(igrid in 1:length(cells)){
    input[input$cell == cells[igrid], ]$TotalUbiomass <- 
      apply(results$U[igrid, params$ref:params$Nx, 1:(params$Neq)] * 
              params$dx*10^params$x[params$ref:params$Nx],
            2, sum, na.rm = T)
    
    input[input$cell == cells[igrid], ]$TotalVbiomass <- 
      apply(results$V[igrid, params$ref:params$Nx, 1:(params$Neq)] * 
              params$dx*10^params$x[params$ref:params$Nx],
            2, sum, na.rm = T)
    
    #sum catches (currently in grams per m3 per year, across size classes) 
    #keep as grams per m2, then be sure to convert observed from tonnes per m2 
    #per year to g.^-m2.^-yr (for each month)
    input[input$cell == cells[igrid], ]$TotalUcatch <- 
      apply(results$Y_u[igrid, params$ref:params$Nx, 1:(params$Neq)] *
              params$dx, 2, sum, na.rm = T)
    
    input[input$cell == cells[igrid], ]$TotalVcatch <- 
      apply(results$Y_v[igrid, params$ref:params$Nx, 1:(params$Neq)] *
              params$dx, 2, sum, na.rm = T)
    
    input[input$cell == cells[igrid], ]$Totalcatch <- 
      input[input$cell == cells[igrid], ]$TotalUcatch +
      input[input$cell == cells[igrid], ]$TotalVcatch
  }
  ## and then multiply outputs by depth to get per m2
  return(input)
}
  
# Set up and run optimisations in parallel ----
fastOptim <- function(LMEnum, vary, fishing_effort_file, forcing_file = NULL, 
                      gridded_forcing = NULL, errorFun = getError, ...,
                      spareCores = 1, libraries = c("optimParallel")){
  
  # get inputs for LME
  if(!is.null(forcing_file)){
    lme_input <- get_lme_inputs(forcing_file = forcing_file, 
                                fishing_effort_file = fishing_effort_file, 
                                LMEnumber = LMEnum)
  }
  if(!is.null(gridded_forcing)){
    lme_input <- get_lme_inputs(gridded_forcing = gridded_forcing, 
                                fishing_effort_file = fishing_effort_file, 
                                LMEnumber = LMEnum)
  }
  
  args <- list(...)
  if("corr" %in% names(args)){
    corr <- args$corr
  }else{
    corr <- F
  }
  if("figure_folder" %in% names(args)){
    figure_folder <- args$figure_folder
  }else{
    figure_folder <- NULL
  }
  
  # set up workers
  # keep some spare core
  noCores <- detectCores()-spareCores 
  if(noCores < 1){
    stop("You should allow at least one core for this operation.")
  }
  cl <- makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, varlist = c("cl", "libraries", "lme_input"),
                envir = environment())
  clusterEvalQ(cl, {
    for(item in 1:length(libraries)){
      library(libraries[item], character.only = T)
    }
  })
  
  clusterEvalQ(cl, source("LME_calibration.R"))
  
  optim_result <- optimParallel(par = vary, fn = errorFun,
                                method = "L-BFGS-B",
                                lower = rep(1e5, length(vary)),
                                upper = rep(2, length(vary)),
                                parallel = list(loginfo = TRUE,
                                                cl = cl,
                                                forward = TRUE), 
                                lme_forcings = lme_input, corr = corr, 
                                figure_folder = figure_folder)
  stopCluster(cl)
  
  return(optim_result$par)
}


