###### Run LHS search for each LME to get best values for fishing parameters

source("LME_calibration.R")
library(dplyr)
library(pbapply)
library(readr)
library(purrr)
library(data.table)

# faster using pbsapply, in the LHSsearch pbapply has cl=6 which uses cluster 
#to run in parallel, but here it is run sequentially if cl is not specified.
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

#Searching best values for fishing parameters for all LMEs
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

#Stacking files before saving
bestvals <- list.files("Output/best_fish_vals", full.names = T) |> 
  map(~fread(.)) |> 
  map_df(~bind_rows(.)) |> 
  arrange(region)
  
#Saving data frame
write_csv(bestvals, file_out)

# WARNINGS: 
# 2. some LME (e.g. LME 4) do not run because the model called by LHSsearch 
# does not produce biomass 
# increase to 64 # even worst 
# decrease to 0.064 # now we have biomass!

# 3. catches are always too small compared to observed data 
# F.mort estimated in LHSsearch can only go to 1, 
# so increase effort in get_lme_inputs() by 
# using relative effort (effort_m2/max(effort_m2), with the highest value 
# being 1)  
# relative effort for each LME - not working well if search_vol = 0.064 as 
# catches remain low  
# doing this also means that effort is equal across LMEs (from 1 to close to 0
# in each LME)
# relative effort across LMEs - same as above and do not use (Julia)   

# now try increasing search_vol again to 0.64
# relative effort for each LME - not working as well as the above
# relative effort across LMEs - working well (best option) but LME 4 and others 
# not working again. 

# now estimating search vol + relative effort for each LME + iter = 100 (but 
# will need to increase to 1000 at least)


### Creating plots with best values for fishing parameters -------------

#### Check other model performance indicators using the above estimates
#bestvals<-data.frame(readRDS("bestvals_LMEs.RDS")) # these bestvalues don't 
#give the CalibrationPlot 

#Load best fishing parameters
bestvals <- read_csv(file_out)

#Create empty data frame to store results
bestvals_fit <- data.frame()

#Getting correlation values for each LME
for(i in bestvals$region){
  #Get forcing data for each LME
  lme_input <- get_lme_inputs(forcing_file = forcing_file, 
                              fishing_effort_file = fishing_effort_file, 
                              LMEnumber = LMEnum)
  
  #Calculate errors and correlations with tuned fishing parameters and save
  #plots
  out_folder <- file.path("Output", paste0("bestvals_LMEs_iter_", no_iter))
  out <- getError(lhs_params = unlist(bestvals[i,]), lme_forcings = lme_input, 
                  corr = T, figure_folder = out_folder)
  #Adding results to data frame created above
  bestvals_fit <- bestvals_fit |> 
    bind_rows(out)
  
  #Saving results at every step
  bestvals_fit |> 
    write_csv("Output/bestvals_corr_na-vals.csv")
}

bestvals <- bestvals |> 
  left_join(bestvals_fit, by = join_by(region))


#### TO DO: Check other model performance indicators using the above estimates

# mean total consumer biomass density
# biomass vs sst relationship
# biomass vs phyto (+ zoo) biomass relationship
# slope and intercept of size spectrum
# P:B ratio
# mean size and TL of catch
# correlation with catch time series
# correlation with relative biomass time series (and RAM Legacy)
# OTHER : mse, rpi - hipsey et al metrics for total catch and for disaggregated
# benthic ad pelagic catch

# Use optimParallel to get better "bestvals" for LMES that do not have good 
#enough fit to data

# refine<-which(bestvals$cor<0.5|bestvals$rmse>0.5|bestvals$catchNA>0)
# tic()
# for (i in 1:dim(bestvals[refine,])[1]){
#   vals<-unlist(bestvals[refine,][i,1:5])
  # optim_result<-fastOptim(LMEnum=refine[i],vary=vals, 
  #                         fishing_effort_file = fishing_effort_file, 
  #                         forcing_file = forcing_file, gridded_forcing = NULL,
  #                         errorFun = getError, corr = T, 
  #                         figure_folder = NULL)


#   bestvals[refine,][i,1:5]<-unlist(optim_result$par)[1:5]
#   bestvals[refine,][i,6]<-unlist(optim_result$value)
#   print(i)
# }
# toc()

to_be_refined <- bestvals |> 
  filter(cor < 0.5 | rmse > 0.5 | catchNA > 0)

for(i in 1:nrow(to_be_refined)){
  result <- LHSsearch(LMEnum = to_be_refined$region[[i]], num_iter = 1000, 
                      search_vol = "estimated", forcing_file = forcing_file, 
                      gridded_forcing = NULL, 
                      fishing_effort_file = fishing_effort_file, corr = F, 
                      figure_folder = NULL, 
                      best_val_folder = "Output/optimised_fish_vals")
  
  #Print message 
  print(paste0("LHS refined search completed for LME #", 
               to_be_refined$region[[i]]))
}

#Load all files with refined fishing parameters for under-performing regions
refined <- list.files("Output/optimised_fish_vals", full.names = T) |> 
  map(~read_csv(.)) |> 
  map_df(~bind_rows(.)) 


#File path for file where parameters will be stored
file_out <- file.path("Output", 
                      paste0("refined-fishing-parameters_LMEs_searchvol_", 
                             search_vol, "_numb-iter_", no_iter, "-1000.csv"))

# then we need to put these together with the other "bestvals", and save 
# results in a single file
bestvals |> 
  filter(cor > 0.5 | rmse < 0.5 | catchNA < 0) |> 
  bind_rows(refined) |> 
  write_csv(file_out)



