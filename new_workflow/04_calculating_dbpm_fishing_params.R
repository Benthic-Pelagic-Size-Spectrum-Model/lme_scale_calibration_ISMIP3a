
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
results_folder <- paste0("/g/data/vf71/la6889/dbpm_inputs/east_antarctica/", 
                         "best_fish_vals_", region_int)
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
  write_json(paste0("/g/data/vf71/la6889/dbpm_inputs/east_antarctica/",
                    "dbpm_size_params_", region_int, ".json"), digits = 10)


# Loading DBPM parameters -------------------------------------------------
# If parameters were already saved, they can be read instead of being 
# recalculated
params <- read_json(paste0("/g/data/vf71/la6889/dbpm_inputs/east_antarctica/",
                           "dbpm_size_params_", region_int, ".json"), 
                    simplifyVector = T)

# Loading data ------------------------------------------------------------
fish_params <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica",
                         "best_fish_vals_fao-58",
                         paste0("best-fishing-parameters_FAO-58_searchvol_",
                                "estimated_numb-iter_1000.parquet")) |> 
  read_parquet() |> 
  slice(1)

dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                         "monthly_weighted", 
                         "dbpm_clim-fish-inputs_fao-58_1841-2010.parquet") |> 
  read_parquet()


# Run non-spatial DBPM ----------------------------------------------------
# This step is necessary to get the initial conditions to be used in the gridded
# DBPM
params <- sizeparam(dbpm_inputs, fish_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3, tstepspryr = 12)

init_results <- run_model(fish_params, dbpm_inputs, withinput = F)


# Prepare fishing parameters for gridded DBPM -----------------------------
pred_initial <- rowMeans(init_results$predators)
detritivore_initial <- rowMeans(init_results$detritivores)
detrititus_initial <- mean(init_results$detritus)

gridded_params <- sizeparam(dbpm_inputs, fish_params, xmin_consumer_u = -3, 
                            xmin_consumer_v = -3, tstepspryr = 12, use_init = T, 
                            pred_initial = pred_initial, 
                            detritivore_initial = detritivore_initial, 
                            detrititus_initial = detrititus_initial,
                            gridded = T)

#Save for use in gridded DBPM (step 05)
gridded_params |> 
  write_json(file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                       "gridded_params/dbpm_gridded_size_params_fao-58.json"),
             digits = 10)

