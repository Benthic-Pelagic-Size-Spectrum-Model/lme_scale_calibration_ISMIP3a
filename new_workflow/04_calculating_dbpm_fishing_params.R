
# Loading libraries -------------------------------------------------------
source("new_workflow/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)


# Loading DBPM climate and fishing inputs ---------------------------------
region_int <- "fao-58"
dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica/", 
                         "monthly_weighted", 
                         paste0("dbpm_clim-fish-inputs_", region_int, 
                                "_1841-2010.parquet")) |> 
  read_parquet()


# Searching best fishing parameters values for area of interest -----------
#Path to folder where results will be stored
results_folder <- paste0("new_workflow/outputs/best_fish_vals_", region_int)
#Number of iterations
no_iter <- 100
params_calibration <- LHSsearch(num_iter = no_iter, forcing_file = dbpm_inputs, 
                                gridded_forcing = NULL, 
                                best_val_folder = results_folder, 
                                best_param = F)

## Creating plots with fishing parameters calculated above --------------
# Calculate errors and correlations with tuned fishing parameters and save plot
params_corr <- params_calibration |> 
  split(params_calibration$rmse) |>
  map_df(\(x) getError(x, dbpm_inputs, corr = T))

#Adding correlation to fishing parameter data frame
params_calibration <- params_calibration |> 
  select(!region) |> 
  bind_cols(params_corr) |> 
  filter(cor >= 0) |> 
  arrange(desc(cor), rmse) |> 
  relocate(region, .before = fmort_u)

#Saving results
params_calibration |> 
  write_parquet(file.path(
    results_folder,
    paste0("best-fishing-parameters_", region_int, 
           "_searchvol_estimated_numb-iter_100.parquet")))

params_calibration |> 
  slice(1) |> 
  corr_calib_plots(dbpm_inputs, results_folder)

params_calibration |> 
  arrange(rmse) |>
  slice(1) |> 
  corr_calib_plots(dbpm_inputs, results_folder)


## Optimising underperforming regions --------------------------------------
# Since correlation is below 0.5 and the plots comparing estimates and obs do 
# not look like a great fit, we will calculate fishing parameters again 
no_iter <- 500
params_calibration_optim <- LHSsearch(num_iter = no_iter, seed = 42,
                                      forcing_file = dbpm_inputs, 
                                      gridded_forcing = NULL, 
                                      best_val_folder = results_folder, 
                                      best_param = F)

params_corr <- params_calibration_optim |> 
  split(params_calibration_optim$rmse) |>
  map_df(\(x) getError(x, dbpm_inputs, corr = T))

#Adding correlation to fishing parameter data frame
params_calibration_optim <- params_calibration_optim |> 
  select(!region) |> 
  bind_cols(params_corr) |> 
  filter(cor >= 0) |> 
  arrange(desc(cor), rmse) |> 
  relocate(region, .before = fmort_u)

#Saving results
params_calibration_optim |> 
  write_parquet(file.path(
    results_folder,
    paste0("best-fishing-parameters_", region_int, 
           "_searchvol_estimated_numb-iter_", no_iter, ".parquet")))

params_calibration_optim |> 
  slice(1) |> 
  corr_calib_plots(dbpm_inputs, results_folder)

params_calibration_optim |> 
  arrange(rmse) |> 
  slice(1) |> 
  corr_calib_plots(dbpm_inputs, results_folder)


# Getting DBPM parameters -------------------------------------------------
fishing_params <- params_calibration_optim |> 
  arrange(rmse) |> 
  slice(1)

params <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3, tstepspryr = 12)

# Saving parameters
params |> 
  #Ensuring up to 10 decimal places are saved in file
  write_json(paste0("new_workflow/outputs/dbpm_size_params_", region_int, 
                    ".json"), digits = 10)


# Loading DBPM parameters -------------------------------------------------
# If paramaters were already saved, they can be read instead of being 
# recalculated
params <- read_json(paste0("new_workflow/outputs/dbpm_size_params_",
                           region_int, ".json"), 
                    simplifyVector = T)





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
  #Adding 1 to ignore initialisation values
  U.initial <- rowMeans(initial_results$U[, ind+1])
  V.initial <- rowMeans(initial_results$V[, ind+1])
  W.initial <- mean(initial_results$W[ind+1])
  
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
  tmax <- length(unique(lme_inputs_grid$year))
  
  # set up params for each month, across grid cells
  gridded_params <- sizeparam(equilibrium = FALSE, dx = 0.1, 
                              xmin.consumer.u = -3, xmin.consumer.v = -3,
                              tmax = tmax, tstepspryr = 12,
                              search_vol = f.effort["search.vol"], 
                              fmort.u = f.u, fminx.u = f.minu, fmort.v = f.v, 
                              fminx.v = f.minv, use.init = TRUE,
                              pred_initial = U.initial, 
                              detritivore_initial = V.initial, 
                              detrititus_initial = W.initial, 
                              Ngrid = nrow(lme_inputs_grid$depth))
  
  # save inputs and parameters object needed for plotting 
  save(gridded_params, 
       file = file.path(LME_path_full, 
                        paste0("grid_inputs_params_LME_", LMEnumber, ".RData")))
}
