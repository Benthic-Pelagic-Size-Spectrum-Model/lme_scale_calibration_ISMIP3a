
# Loading libraries -------------------------------------------------------
library(arrow)
library(jsonlite)
library(ggplot2)
source("new_workflow/useful_functions.R")


# Loading data ------------------------------------------------------------
fish_params <- file.path("new_workflow/outputs/best_fish_vals_fao-58", 
                         paste0("best-fishing-parameters_FAO-58_searchvol_",
                                "estimated_numb-iter_1000.parquet")) |> 
  read_parquet()

dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                         "monthly_weighted", 
                         "dbpm_clim-fish-inputs_fao-58_1841-2010.parquet") |> 
  read_parquet() |> 
  filter(year == 1961)




# Run non-spatial DBPM ----------------------------------------------------
# This step is necessary to get the initial conditions to be used in the gridded
# DBPM
params <- sizeparam(dbpm_inputs, fish_params[1,], xmin_consumer_u = -3, 
                    xmin_consumer_v = -3, tstepspryr = 12)

init_results <- run_model(fish_params[1,], dbpm_inputs, withinput = F)


# Prepare fishing parameters for gridded DBPM -----------------------------
pred_initial <- rowMeans(init_results$predators)
detritivore_initial <- rowMeans(init_results$detritivores)
detrititus_initial <- mean(init_results$detritus)

gridded_params <- sizeparam(dbpm_inputs, fish_params[1,], xmin_consumer_u = -3, 
                            xmin_consumer_v = -3, tstepspryr = 12, use_init = T, 
                            pred_initial = pred_initial, 
                            detritivore_initial = detritivore_initial, 
                            detrititus_initial = detrititus_initial,
                            gridded = T)

#Save for use in gridded DBPM (step 06)
gridded_params |> 
  write_json("new_workflow/gridded_params_testing.json", digits = 10)














result_set <- sizemodel(params)



plotsizespectrum(result_set, params, timeaveraged = F)

  
dbpm_sau <- read_parquet("data/dbpm_clim-fish-inputs_fao-58_1841-2010.parquet")

no_iter <- 100
params_calibration <- LHSsearch(num_iter = no_iter, forcing_file = dbpm_sau, 
                                gridded_forcing = NULL, best_param = F,
                                best_val_folder = 
                                  file.path("new_workflow/outputs",
                                            "best_fish_vals_sau_fixed_code"))

params_calibration <- read_parquet("new_workflow/outputs/best_fish_vals_sau/best-fishing-parameters_FAO-58_searchvol_estimated_numb-iter_100.parquet")
params_sau <- sizeparam(dbpm_sau, params_calibration[1,], xmin_consumer_u = -3, 
                        xmin_consumer_v = -3, tstepspryr = 12)
result_sau <- sizemodel(params_sau)
plotsizespectrum(result_sau, params_sau, timeaveraged = T)




library(purrr)
test <- params_calibration |> 
  split(params_calibration$rmse) |>
  map_df(\(x) getError(x, dbpm_sau, corr = T))

params_calibration <- params_calibration |> 
  select(!region) |> 
  bind_cols(test) |> 
  arrange(desc(abs(cor)), rmse)

#8 (maybe), 38 (lowest RMSE)
test_params <- params_calibration[38,]
params_sau <- sizeparam(dbpm_sau, test_params, xmin_consumer_u = -3, 
                        xmin_consumer_v = -3, tstepspryr = 12)

result_sau <- sizemodel(params_sau)
corr_calib_plots(test_params, dbpm_sau, 
                 "new_workflow/outputs/best_fish_vals_sau")

plotsizespectrum(result_sau, params_sau, timeaveraged = T)



corr_calib_plots(fish_params, dbpm_sau, "data")




test <- dbpm_inputs |> 
  filter(year == 1960)
params_test <- sizeparam(test, fish_params, xmin_consumer_u = -3, 
                        xmin_consumer_v = -3, tstepspryr = 12)

sizemodel_test(params_test)


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







library(tidyverse)
sau <- read_csv("data/SAU FAO 58 v50-1.csv")
names(sau)

sau_prep <- sau |> 
  filter(year <= 2010) |> 
  group_by(year) |> 
  summarise(catch_tonnes = sum(tonnes, na.rm = T)) |> 
  mutate(region = 58, .before = catch_tonnes) |> 
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         catch_tonnes_area_m2 = catch_tonnes/area_m2) 

DBPM_effort_catch_input <- effort_data |> 
  full_join(sau_prep) |> 
  mutate(region = paste0("FAO ", region))

#Saving summarised catch and effort data
DBPM_effort_catch_input |> 
  write_parquet("data/dbpm_effort-catch-inputs_fao-58.parquet")

#Removing individual data frames
rm(effort_data, catch_data)

#Joining with climate inputs
forcing_file <- clim_forcing_file |> 
  full_join(DBPM_effort_catch_input, by = c("region", "year")) 

forcing_file |> 
  write_parquet("data/dbpm_clim-fish-inputs_fao-58_1841-2010.parquet")


ggplot()+
  geom_line(data = sau_prep, aes(year, catch_tonnes), color = "red")+
  geom_line(data = catch_data, aes(year, catch_tonnes), color = "blue")

