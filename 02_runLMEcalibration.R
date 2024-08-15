###### Run LHS search for each LME to get best values for fishing parameters


# Loading libraries -------------------------------------------------------
library(dplyr)
library(pbapply)
library(parallel)
library(purrr)
library(data.table)
source("LME_calibration.R")

# Defining basic variables ------------------------------------------------
no_iter <- 100
# other option is to specify a value for search_vol 
search_vol <- "estimated" 

#Fishing effort file location
fishing_effort_file <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                                 "processed_forcings/lme_inputs/obsclim/025deg",
                                 "DBPM_LME_effort_catch_input.csv")
#Non-gridded data
forcing_file <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                          "processed_forcings/lme_inputs/obsclim/025deg", 
                          "DBPM_LME_climate_inputs_slope.csv")
#Gridded data
gridded_forcing <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                             "processed_forcings/lme_inputs_gridcell/obsclim",
                             "025deg")



# Searching best fishing parameters values for all LMEs -------------------
# Faster using pbsapply, LHSsearch uses mclapply to run in parallel, but here 
# it is run sequentially if cl is not specified.
lmes_best_params <- t(pbsapply(X = 1:66, LHSsearch, num_iter = no_iter, 
                               search_vol = search_vol, 
                               forcing_file = forcing_file, 
                               gridded_forcing = NULL, 
                               fishing_effort_file = fishing_effort_file, 
                               corr = F, figure_folder = NULL, 
                               best_val_folder = "Output/best_fish_vals"))

#File path for file where parameters will be stored
file_out <- file.path("Output", 
                      paste0("best-fishing-parameters_LMEs_searchvol_", 
                             search_vol, "_numb-iter_", no_iter, ".csv"))


## Creating plots with best values for fishing parameters --------------
#Calculate errors and correlations with tuned fishing parameters and save plots

#Load best fishing parameters
bestvals <- fread(file_out)

#Path to folder where figures will be stored
figure_folder <- file.path("Output", paste0("bestvals_LMEs_iter_", no_iter))

bestvals_fit <- pbapply(bestvals, 1, corr_calib_plots, forcing_file, 
                        fishing_effort_file, figure_folder = figure_folder, 
                        cl = (detectCores()-4))

#Adding correlation to fishing parameter data frame
bestvals_fit <- bestvals_fit |> 
  bind_rows() |> 
  right_join(bestvals, by = join_by(region)) |> 
  relocate(c(region, cor, catchNA), .after = rmse)

#Saving results
out_file <- file.path(figure_folder, 
                      paste0("best-fishing-parameters_LMEs_searchvol_", 
                             search_vol, "_numb-iter_", no_iter, ".csv"))

bestvals_fit |> 
  fwrite(out_file)


## Optimising underperforming regions --------------------------------------
bestvals_fit <- fread(out_file)

to_be_refined <- bestvals_fit |> 
  filter(cor < 0.5 | rmse > 0.5 | catchNA > 0)

f_out <- "Output/best_fish_vals/optimised_underperforming_LMEs"
refined_best_params <- t(pbsapply(X = to_be_refined$region, LHSsearch, 
                                  num_iter = 1000, search_vol = "estimated", 
                                  forcing_file = forcing_file, 
                                  gridded_forcing = NULL, 
                                  fishing_effort_file = fishing_effort_file, 
                                  corr = F, figure_folder = NULL, 
                                  best_val_folder = f_out))

#Load all files with refined fishing parameters for under-performing regions
refined_best_params <- list.files(f_out, full.names = T) |> 
  map(~fread(.)) |> 
  map_df(~bind_rows(.)) |> 
  arrange(region)

#Figure folder
figure_folder <- file.path(figure_folder, "optimised_underperforming_LMEs")

refined_vals_fit <- pbapply(refined_best_params, 1, corr_calib_plots,
                            forcing_file, fishing_effort_file, 
                            figure_folder = figure_folder, 
                            cl = (detectCores()-4))

#Adding correlation to fishing parameter data frame
refined_vals_fit <- refined_vals_fit |> 
  bind_rows() |> 
  right_join(refined_best_params, by = join_by(region)) |> 
  relocate(c(region, cor, catchNA), .after = rmse)

#Creating file name to save optimised parameters only
f_out <- file.path(figure_folder, 
                   paste0("optimised-fishing-parameters_LMEs_searchvol_", 
                          search_vol, "_numb-iter_1000.csv"))

#Saving optimised fishing parameters
refined_vals_fit |> 
  fwrite(f_out)

#File path for file where parameters will be stored
file_out <- file.path("Output", 
                      paste0("refined-fishing-parameters_LMEs_searchvol_", 
                             search_vol, "_numb-iter_", no_iter, "-1000.csv"))


## Merging best fishing parameters for all LMEs ---------------------------
# Merging optimised parameters with the regions that did not need to be 
# optimised from the "bestvals_fit" variable. Saving results in a single file
bestvals_fit |> 
  #Removing the regions that needed to be refined
  filter(!(cor < 0.5 | rmse > 0.5 | catchNA > 0)) |> 
  #Replacing with new refined values
  bind_rows(refined_vals_fit) |> 
  arrange(region) |> 
  fwrite(file_out)

