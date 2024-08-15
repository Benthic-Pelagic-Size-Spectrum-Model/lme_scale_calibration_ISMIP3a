#' @param search_vol Numeric value for the search volume. Default is .64.
#' @param savePlots Boolean value. If true, plots related to the selected LME 
#' will
#' be saved in the same folder as the grid output. Default is TRUE.



###############################################################################
# Library of functions developed originally by Ryan/Julia/Cami in 2022-2023
#
# Functions were updated by Denisse Fierro Arcos 
# Date of update: 2024-08-08
#
###############################################################################

# Loading libraries
#LME in Southern Ocean is called Antarctica (ID merged 147) LME 61
#From /rd/gem/private/fishmip_inputs/ISIMIP3a/fishmip_regions/FAO-LME_masks

library(tictoc)
library(dplyr)
library(lubridate)
library(data.table)
library(parallel)
# library(pbapply)

# also calls dbpm_model_functions.R
source("LME_calibration.R") 
# source("Plotting_functions_DBPM.R")


# Defining base variables -------------------------------------------------
#Location of folder where outputs will be saved
LME_path <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/DBPM/obsclim"

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

#Data frame containing information to process input files
meta <- data.frame(region = 1:66) |> 
  mutate(base_out = file.path(LME_path, paste0("LME_", region, "_output")), 
         grid = file.path(base_out, 
                          "obsclim_historical_DBPM_LME_inputs_gridded.csv"), 
         non_grid = file.path(base_out, 
                              "obsclim_historical_DBPM_LME_inputs_non-gridded.csv"))


# Saving forcing data to disk ---------------------------------------------
lapply(1:nrow(meta), 
       FUN = function(i) gridded_stable_spinup(unlist(meta[i,]), 
                                               gridded_forcing, 
                                               fishing_effort_file))

lapply(1:nrow(meta), 
       FUN = function(i) get_lme_inputs(forcing_file = forcing_file, 
                                        fishing_effort_file = fishing_effort_file,
                                        LMEnumber = meta[i,]$region,
                                        file_out = meta[i,]$non_grid))

# get the latest best values based on iterations and search_vol
# no_iter <- 100
# search_vol <- "estimated" 
# version <- "refined-fishing-parameters"
f.effort <- list.files("Output/", "refined-fishing", full.names = T) |> 
  fread() |> 
  # run only LMEs with a good correlation/error and all catches 
  filter(cor > 0.5 & rmse < 0.5 & catchNA == 0) 

meta <- meta |> 
  filter(region %in% f.effort$region)

#Run model across space and time using gridded inputs for each LME
mclapply(1:nrow(meta), 
         FUN = function(i) calc_grid_params(meta = unlist(meta[i,]),
                                            f.effort = unlist(f.effort[i,]), 
                                            start_cond = NULL), 
         mc.cores = round((detectCores()*.75), 0))




#### run gridded model by LME ----
#Run model across space and time using gridded inputs for each LME
  




rungridbyLME <- function(meta){
  load(list.files(meta["base_out"], "grid_inputs_params_", full.names = T))
  
  lme_input_init <- fread(meta["non_grid"])
  grid_results <- vector("list", nrow(lme_input_init))
  # run model  for full time period across all grid cells
  grid_results <- gridded_sizemodel(gridded_params, ERSEM.det.input = F,
                                    temp.effect = T, eps = 1e-5, 
                                    output = "aggregated", use.init = TRUE,
                                    burnin.len)
  
  #### TEST 3 - OK working 
  # LME 1
  # random bestvalues - i.e. the ones that picked manually best approximate
  # catches (0.1,0.5,1,1) 
  # search vol = 0.64 as OK for this LME
  # Fmort = first spread effort then calculate Fmort and catches as discussed 
  # with Julia 
  # gravity model option 2 with iter = 1 
  
  # removing the stable spinup section to have matching dimensions with the code
  # WARNING  move to plotting function for now as need to figure out catch trend
  ind <- which(ymd(colnames(gridded_params$er)) >= min(lme_input_init$t))
  
  grid_results$U <- grid_results$U[, , ind]
  grid_results$V <- grid_results$V[, , ind]
  grid_results$Y.u <- grid_results$Y.u[, , ind]
  grid_results$Y.v <- grid_results$Y.v[, , ind]
  # moved to plotting function as needed there
  # gridded_params$Neq <- 2040
  
  # save results from run
  fwrite(grid_results, file.path(LME_path_full, "gridded_model_results.csv"))
}


## RUN all LMEs 

# tic()
# pbapply::pbsapply(X=LMEnumber,rungridbyLME,yearly = FALSE, f.effort = TRUE,
# vals = vals, cl = detectCores() - 5)
# toc() # trying this now BUT: WARNING: Session forced to suspend due to system
# upgrade, restart, maintenance, or other issue. Your session data was saved 
# however running computations may have been interrupted.
# 
# ## OR 
# tic()
# mclapply(LMEnumber, function(x) rungridbyLME(x, yearly = FALSE, 
# f.effort = TRUE, vals = vals), mc.cores = detectCores()-5)
# toc() # 1536.929 sec elapsed/25 min. No working with all LMEs ?!

## OR standard loop (possibly very slow)

tic()
for(i in 1:length(LMEnumber)){
  # # trial 
  i = 42 # 41 is LME 61 considering the selection above
  print(paste0("Now working on LME", LMEnumber[i]))

  calc_grid_params(LMEnumber = LMEnumber[i],
               yearly = FALSE, # for get_lme_inputs()
               f.effort = f.effort) # for rungridbyLME()
}

toc()

# STACK at LME 61 

#### plot gridded and averaged output by LME ----

LMEnumber <- LMEnumber[1:40]

### now produce plots - EMPTY when run in // 
# but OK when run individually and from inside the function. 
tic()
mclapply(LMEnumber, function(x) plotgridbyLME(x), mc.cores = detectCores()-5)
toc() 
# problem running this too... 

#### global map of gridded output ----

# this function is very similar to the beginning of plotgridbyLME() above but
# considers only gridded outputs
# wth the final aim of aggregating all outputs to plot global map 

#### MOVE TO FUNCTION CODE
##### POSSIBLY WORTH DOING ALL PLOTS CALCULATIONS HERE AS THE MOST EXPENSIVE
# ACTION IS TO LOAD THE OUTPUTS - SO BETTER DOING IT ONLY ONCE ... 
# COMBINED WITH PLOTTING FUNCTION ABOVE... 
getGriddedOutputs_decade <- function(LME_path, LMEnumber){
  # # trial 
  # LMEnumber = 12
  
  # load outputs 
  LME_path_full <- paste0(LME_path, "LME_", LMEnumber, "_output")
  full_file_name <- paste0(LME_path_full, "/grid_results.rds")
  full_file_name_output <- paste0(LME_path_full, 
                                  "/grid_results_toCheckMap.rds")
  
  if(file.exists(full_file_name)){
    print(paste0("Now working on LME", LMEnumber))
    tic()
    grid_results <- readRDS(paste0(LME_path_full, "/grid_results.rds"))
    toc() # 49 sec
    
    # # ### WARNING need to comment out to figure out trend in biomass 
    # # removing the stable spinup section to have matching dimensions with the 
    # code
    # grid_results$U <- grid_results$U[,,1201:3241]
    # grid_results$V <- grid_results$V[,,1201:3241]
    # grid_results$Y.u <- grid_results$Y.u[,,1201:3241]
    # grid_results$Y.v <- grid_results$Y.v[,,1201:3241]
  
    # load inputs and param object 
    load(paste0(LME_path_full, "/grid_results_inputs_params.RData"))
    gridded_params$Neq <- 2040
  
    ### WARNING need to check depth integration and Neq - do they match what
    # used for inputs and for runLMEcalibration? 
    tic()
    out <- getGriddedOutputs(input = lme_inputs_grid, results = grid_results,
                             params = gridded_params)
    toc() # 18 sec
    
    # averaged across decades for U and V
    tic()
    out2 <- out |> 
      select(lat, lon, LME, t, TotalVcatch, TotalUcatch, TotalVbiomass, 
             TotalUbiomass) |> 
      mutate(year = year(t), 
             # original data resolution is month - can do year than decade. 
             decade = floor(year/10)*10) |> 
      group_by(LME, decade, lat, lon) |> 
      summarize(Ubiomass_year = mean(TotalUbiomass, na.rm = T),
                Vbiomass_year = mean(TotalVbiomass, na.rm = T), 
                Ucatch_year = mean(TotalUcatch, na.rm = T),
                Vcatch_year = mean(TotalVcatch, na.rm = T)) |> 
      ungroup() |> 
      mutate(biomass_year = Ubiomass_year+Vbiomass_year,
             catch_year = Ucatch_year+Vcatch_year)
    toc() # 2 sec
    
    # write.csv(x = out2, file = file.path(full_file_name_output))
    fwrite(x = out2, file = file.path(full_file_name_output))
  }
}

# apply function in //
# for now run first LMEs that worked 
trial <- LMEnumber[1:40] 
tic()
a <- mclapply(trial, function(x) getGriddedOutputs_decade(LME_path, 
                                                          LMEnumber = x), 
              mc.cores = detectCores()-5)
toc() # 81 sec 6 LMEs.  

# get all files and put them together
file_list <- list.files(path = LME_path, pattern = "toCheckMap", 
                        recursive = TRUE, full.name = TRUE)
df_to_plot <- rbindlist(mclapply(file_list, function(x) fread(x), 
                                 mc.cores = detectCores()-5))

### all LMEs togetehr and plot 
plot_global_raster <- function(df_to_plot, variable_to_plot, decade_to_plot){
  # # trial 
  # variable_to_plot = "biomass_year"
  # decade_to_plot = 2010
  df <- df_to_plot |> 
    select(all_of(c("lat", "lon", variable_to_plot, "decade"))) |>  
    rename("variable" = variable_to_plot) |> 
    # unique() |> 
    filter(decade == decade_to_plot) |> 
    select(!decade) 
  
  p1 <- df |> 
    ggplot(aes(lon, lat, fill = variable))+
    geom_tile()

  name <- paste0(LME_path, "global_", variable_to_plot, decade_to_plot, "map",
                 ".pdf")
  
  ggsave(name, p1, device = "pdf")
}

plot_global_raster(df_to_plot = df_to_plot, variable_to_plot = "biomass_year", 
                   decade_to_plot = 2010)
plot_global_raster(df_to_plot = df_to_plot, variable_to_plot = "catch_year", 
                   decade_to_plot = 2010)

#### tesing ----

# Test 1- compare with results from DBPM, no fishing (checking code consistency)
# # 1.	Test 1: run yearly = TRUE, no fishing (effort = 0), search volume = 64. 
# # b.	Compare these plots with the 3a tcb netcdf file: 
# #   i.	Extract LME 14 from this file 
# # ii.	Produce plots
# 
# tcb <- nc_open("/rd/gem/private/fishmip_outputs/ISIMIP3a/ctrlclim/netcdf/
#dbpm_ctrlclim_nobasd_1deg_nat_default_tcb_global_monthly_1961_2010.nc")
# #tcb <- nc_open("/rd/gem/private/fishmip_outputs/ISIMIP3a/ctrlclim/netcdf/
#dbpm_ctrlclim_nobasd_1deg_nat_default_tcblog10_global_monthly_1961_2010.nc")
# 
# 
# lat_mask <- tcb$dim$lat$vals >= min(lat_ext) & tcb$dim$lat$vals <= 
# max(lat_ext)
# lon_mask <- tcb$dim$lon$vals >= min(lon_ext) & tcb$dim$lon$vals <= 
# max(lon_ext)
# time_mask <- tcb$dim$time$vals >= 1841 & tcb$dim$time$vals <= 2010


# library(ncdf4)
# tcb <- nc_open("/rd/gem/private/fishmip_outputs/ISIMIP3a/ctrlclim/netcdf/
# dbpm_ctrlclim_nobasd_1deg_nat_default_tcb_global_monthly_1961_2010.nc")
# lat_mask <- tcb$dim$lat$vals >= min(lat_ext) & tcb$dim$lat$vals <= 
# max(lat_ext)
# lon_mask <- tcb$dim$lon$vals >= min(lon_ext) & tcb$dim$lon$vals <= 
# max(lon_ext)
# time_mask <- tcb$dim$time$vals >= 1841 & tcb$dim$time$vals <= 2010
# 
# ncdData <- ncvar_get(tcb, "tcb", 
#                      start = c(which(lon_mask)[1], which(lat_mask)[1], 
                       #which(time_mask)[1]), 
#                      count = c(sum(lon_mask), sum(lat_mask), sum(time_mask)))
# dimnames(ncdData) <- list(tcb$var$tcb$dim[[1]]$vals[which(lon_mask)],
#                        tcb$var$tcb$dim[[2]]$vals[which(lat_mask)],
#                        tcb$var$tcb$dim[[3]]$vals[which(time_mask)])
# names(dimnames(ncdData)) <- c(tcb$var$tcb$dim[[1]]$name,
#                            tcb$var$tcb$dim[[2]]$name,
#                            tcb$var$tcb$dim[[3]]$name)
# 
# tcb_df <- reshape2::melt(ncdData)
# 
# tcb_decade_avg <- tcb_df |>
#   mutate(decade = as.integer(substr(time, 1, 3)) * 10) |> 
#   group_by(decade, lon, lat) |>  
#   summarize(avg_biomass = mean(value))
# 
# # values from test1
# out <- readRDS("rungridRes/test1.rds")
# biom_df <- out[,c(1,2,4,16,17)]
# biom_df <- biom_df |> mutate(totalB = TotalVbiomass + TotalUbiomass)
# biom_decade_avg <- biom_df |>
#   mutate(decade = as.integer(substr(t, 1, 3)) * 10) |> 
#   group_by(decade, lon, lat) |>  
#   summarize(avg_biomass = mean(totalB))
# 
# # different cells so cannot just substract 
# 
# p1 <- ggplot(tcb_decade_avg)+
#   geom_tile(aes(x = lon, y = lat, fill = avg_biomass)) +
#   geom_sf(data = world) +
#   coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
#   scale_fill_gradient2(low = "white", high = "red", 
#   name = "Average Biomass in g/m2") +
#   facet_wrap(~decade,ncol = 6) +
#   scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2],
#   by = 6)) +
#   scale_y_continuous(name = "Latitude") +
#   theme(legend.position = "bottom") +
#   ggtitle("tcb netcdf")
# 
# p2 <- ggplot(biom_decade_avg)+
#   geom_tile(aes(x = lon, y = lat, fill = avg_biomass)) +
#   geom_sf(data = world) +
#   coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
#   scale_fill_gradient2(low = "white", high = "red", 
#   name = "Average Biomass in g/m2") +
#   facet_wrap(~decade,ncol = 6) +
#   scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], 
#   by = 6)) +
#   scale_y_continuous(name = "Latitude") +
#   theme(legend.position = "bottom") +
#   ggtitle("test 1 runGrid")
# 
# pdf("test1res.pdf", onefile = T)
# p1
# p2
# dev.off()
# 
# #Test 2- model calibration/behaviour
# # 1. Spatial maps of TotalCatch and TotalBiomass (by decade)
# # 2. Community size spectra (U & V) - one line per grid cell - final decade?
# # 3. Plots of GG growth rates (see historyMatching.R for example)
# # 4. Compare time series to total catches at LME scale with obs catch
# # 5. Once checked, run for all LMEs
# 
# # TODO cami: compare model with empirical catches (see run lme calibration for
# plotting)
# 
