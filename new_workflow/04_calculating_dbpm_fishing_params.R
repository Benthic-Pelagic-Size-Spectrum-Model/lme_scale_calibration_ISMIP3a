
# Loading libraries -------------------------------------------------------
library(dplyr)
library(arrow)
library(lhs)
library(jsonlite)
library(lubridate)
source("new_workflow/useful_functions.R")


# Loading DBPM climate and fishing inputs ---------------------------------
dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                         "monthly_weighted", 
                         "dbpm_clim-fish-inputs_fao-58_1841-2010.parquet") |> 
  read_parquet()


# Searching best fishing parameters values for area of interest -----------
no_iter <- 100
params_calibration <- LHSsearch(num_iter = no_iter, forcing_file = dbpm_inputs, 
                                gridded_forcing = NULL, 
                                best_val_folder = 
                                  file.path("new_workflow/outputs",
                                            "best_fish_vals"))

## Creating plots with fishing parameters calculated above --------------
#Path to folder where figures will be stored
figure_folder <- "new_workflow/outputs/best_fish_vals"

# Calculate errors and correlations with tuned fishing parameters and save plot
params_correlation  <- corr_calib_plots(params_calibration, dbpm_inputs, 
                                        figure_folder)

#Adding correlation to fishing parameter data frame
params_calibration <- params_calibration |> 
  left_join(params_correlation, by = join_by(region)) |> 
  relocate(region, .before = fmort_u)

params_calibration

## Optimising underperforming regions --------------------------------------
# Since correlation is below 0.5 and the plots comparing estimates and obs do 
# not look like a great fit, we will calculate fishing parameters again 
no_iter <- 1000
params_calibration <- LHSsearch(num_iter = no_iter, forcing_file = dbpm_inputs, 
                                gridded_forcing = NULL, 
                                best_val_folder = 
                                  file.path("new_workflow/outputs",
                                            "best_fish_vals"))







# Getting DBPM parameters -------------------------------------------------
params <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3, tstepspryr = 12)

# Saving parameters
params |> 
  #Ensuring up to 10 decimal places are saved in file
  write_json("new_workflow/outputs/dbpm_size_params.json", digits = 10)


# Loading DBPM parameters -------------------------------------------------
# If paramaters were already saved, they can be read instead of being 
# recalculated
params <- read_json("new_workflow/outputs/dbpm_size_params.json", 
                    simplifyVector = T)


# result_set <- sizemodel(params)



# attach(params)











gridded_params <- sizeparam(,
                            
                            equilibrium = FALSE, dx = 0.1, 
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


