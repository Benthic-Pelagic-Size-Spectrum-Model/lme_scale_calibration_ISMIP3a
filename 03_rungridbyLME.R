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


#### plot gridded and averaged output by LME ---- (after gridded model was run)

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


