###############################################################################
# Library of supporting functions developed originally by Camilla Novaglio
# and Ryan Heneghan. 
# Functions were updated and split into smaller, more specific functions by
# Denisse Fierro Arcos. Functions were also moved to standalone script so they
# can be easily reused in other scripts.
#
# Date of update: 2024-07-02

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(data.table)
library(parallel)
library(dtplyr)
library(janitor)
library(terra)
library(lubridate)
source("dbpm_model_functions.R")

# Reorganising data extracted from GFDL-COBALT2-MOM6 model ----------------
#This function will transform the data so it appears longer (instead of wide)
#It will also add the correct area of the grid cell for the GFDL model
reorganise_gfdl <- function(gfdl_path, area_cell_df, region){
  #Inputs:
  #gfdl_path (character) - Full path to GFDL data 
  #area_cell_df (data frame) - Contains the area of the grid cell from the GFDL
  #model. It should contain three columns: lon, lat, and area_m2
  #region (character) - Region of interest
  #
  #Output:
  #df (data frame) - Contains reorganised GFDL data and adds correct grid cell
  #area
  
  if(str_detect(gfdl_path, "fao")){
    #Load data
    df <- fread(gfdl_path, header = F) |> 
      #Transpose data
      t() |> 
      #Use first row as column names
      row_to_names(1) |> 
      as.data.frame() |> 
      remove_empty("cols") |> 
      #Ensure data is numeric
      mutate(across(everything(), ~as.numeric(.x))) |> 
      #Reorganise data so values are in a single column
      pivot_longer(!c(lat,lon), values_to = "value", names_to = "date") |> 
      #Apply date format
      mutate(date = lubridate::as_date(date)) |> 
      #Add correct area of the grid cell
      left_join(area_cell_df, by = join_by(lat, lon)) |> 
      #Adding scenario to data
      mutate(scenario = str_extract(gfdl_path, 
                                    "(fao|lme)_inputs/(.*)/\\d{1,3}deg", 
                                    group = 2)) |> 
      #Remove any empty rows
      drop_na() |> 
      #Apply date format to date column
      mutate(month = lubridate::month(date, label = T), 
             year = year(date), region = region) |> 
      #Order by date
      arrange(date)
  }else{
    df <- fread(gfdl_path) |> 
      #Remove area column
      select(!area_m2) |> 
      remove_empty("cols") |> 
      #Reorganise data so values are in a single column
      pivot_longer(!c(lat,lon), values_to = "value", names_to = "date") |> 
      #Apply date format
      mutate(date = lubridate::my(date)) |> 
      #Add correct area of the grid cell
      left_join(area_cell_df, by = join_by(lat, lon)) |> 
      #Adding scenario to data
      mutate(scenario = str_extract(gfdl_path, 
                                    "(fao|lme)_inputs/(.*)/\\d{1,3}deg", 
                                    group = 2)) |> 
      #Remove any empty rows
      drop_na() |> 
      #Apply date format to date column
      mutate(month = lubridate::month(date, label = T), 
             year = year(date), region = region) |> 
      #Order by date
      arrange(date)
  }
  
  return(df)
}


# Calculating weighted mean by grid cell area -----------------------------
weighted_mean_by_area <- function(var_df){
  #This function calculates the weighted mean by grid cell area for different 
  #climate variables.
  #
  #Inputs:
  #var_df (data frame) - Data frame output from reorganise_gfdl() function
  #Outputs:
  #weighted_df (data frame) - Data frame including the weighted mean
  weighted_df <- var_df |> 
    group_by(date, month, year, scenario, region) |> 
    summarise(weighted_mean = weighted.mean(value, area_m2),
              total_area_km2 = sum(area_m2)*1e-6) |>
    ungroup()
  return(weighted_df)
}


# Adding spinup to GFDL data ----------------------------------------------
add_spinup_data <- function(spinup_df, weighted_gfdl, variable){
  #This function adds spin up information from ctrlclim data to GFDL variables. 
  #
  #Inputs:
  #spinup_df (data frame) - Containing spinup information from control data
  #weighted_df (data frame) - Data frame from weighted_mean_by_area() function
  #variable (character) - Name of GFDL weighted variable - to be added to data
  #frame
  #Outputs:
  #df_spinup (data frame) - Data frame combining GFDL weighted data and spinup
  #data
  
  df_spinup <- spinup_df |> 
    bind_rows(weighted_gfdl) |> 
    arrange(date) |> 
    # reorder columns 
    relocate(c("region", "date", "year", "month", "scenario", "weighted_mean",
               "total_area_km2")) |> 
    # add name of variable to data frame
    mutate(variable_name = variable)
  
  return(df_spinup)
}


# Calculating weighted means for GFDL data --------------------------------
calc_inputs <- function(path_ctrl, path_obs){
  #This function calculates fixed weighted.mean for depth and total area, 
  # and monthly weighted.mean for each climate variable. It also creates and 
  # adds spin up information from ctrlclim data to GFDL variables. 
  #
  #Inputs:
  #path_ctrl (character) - Full file path for control data
  #path_obs (character) - Full file path for observations data
  #region (character) - Region of interest
  #Outputs:
  #(list) - Containing weighted mean (by area of grid cell) for ctrlclim and
  #obsclim GFDL outputs
  
  #Getting file name for control and observations datasets
  file_name_obs <- basename(path_obs)
  file_name_ctrl <- basename(path_ctrl)
  
  #Extracting variable name from observations filename
  variable <- str_extract(file_name_obs, "clim_(.+)_\\d{2}arcmin",
                          group = 1) |> 
    #If unit is present, it will be disregarded
    str_split_i("_", i = 1) 
  
  #Getting name of region
  region <- str_extract(file_name_ctrl, "arcmin_(.*)_monthly", group = 1)
  
  #Loading grid in the resolution matching the obsclim and ctrlclim data
  area_frame <- list.files(file.path("/g/data/vf71/shared_resources",
                                     "grid_cell_area_ESMs/isimip3a"),
                           #Get resolution of data from file name
                           paste0(str_extract(file_name_ctrl, 
                                              "_\\d{2}arcmin_"), ".*csv$"),
                           full.names = T) |> 
    read_csv() |> 
    #Rename columns
    rename("lon" = "x", "lat" = "y", "area_m2" = "cellareao")
  
  if(variable == "deptho"){
    obsclim_df <- fread(path_obs) |> 
      select(!area_m2) |> 
      mutate(scenario = str_extract(path_obs, 
                                    "(fao|lme)_inputs/(.*)/\\d{1,3}deg", 
                                    group = 2)) |> 
      #Add correct area of the grid cell
      left_join(area_frame, by = join_by(lat, lon))
  }else{
    #Load ctrlclim data
    ctrlclim_df <- reorganise_gfdl(path_ctrl, area_frame, region)
    
    #Load obsclim data
    obsclim_df <- reorganise_gfdl(path_obs, area_frame, region)
  }
  
  # Calculate weighted mean for GFDL variables
  if(variable == "deptho"){
    # calculate fixed variables - mean depth and area of FAO/LME
    # shortcut if you need to extract ctrlclim values too for all variables
    weighted_mean_ctrl_final <- weighted_mean_obs_final <- obsclim_df |> 
      summarise(weighted_mean = weighted.mean(m, area_m2),
                total_area_km2 = sum(area_m2)*1e-6,
                variable_name = "deptho") |> 
      # Add name of region from file path
      mutate(region = region, .before = weighted_mean) |> 
      mutate(date = NA, year = NA, month = NA, 
             scenario = "obsclim", .after = region)
    weighted_mean_ctrl_final <- weighted_mean_ctrl_final |> 
      mutate(scenario = "ctrlclim")
  }else{
    # ctrlclim
    weighted_mean_ctrl <- weighted_mean_by_area(ctrlclim_df)
    # obsclim 
    weighted_mean_obs <- weighted_mean_by_area(obsclim_df)
    
    # SPINUP, based on ctrlclim and used for both ctrlclim and obsclim
    spinup <- weighted_mean_ctrl |>
      filter(year >= 1961, year <= 1980) |> 
      slice(rep(1:n(), times = 6)) |>
      mutate(year = rep(1841:1960, each = 12),
             date = my(paste(month, year, sep = "-")),
             scenario = "spinup") 
    
    # add spinup to control
    weighted_mean_ctrl_final <- add_spinup_data(spinup, weighted_mean_ctrl, 
                                                variable)
    
    # add spin up to obsclim
    weighted_mean_obs_final <- add_spinup_data(spinup, weighted_mean_obs, 
                                               variable)
  }
  return(list(weighted_mean_obs_final = weighted_mean_obs_final, 
              weighted_mean_ctrl_final = weighted_mean_ctrl_final))
}



# Flattening list of lists based on sublist names -------------------------
flatten_list_list <- function(list_list, name_pattern, column_values){
  # This function binds together data frames contained within a list of lists
  # based on the sublist names identified with the pattern given. 
  #
  # Inputs:
  # list_list (list) - A list with sublists containing data frames. The sublist
  # names must share a similar name to identify correct data frames to merge
  # name_pattern (character) - Pattern that helps identify data frames to be
  # merged
  #
  # Outputs:
  # df_flatten (data frame) - Flatten list of lists
  #obsclim GFDL outputs
  
  #Select data frame within sublist based on name pattern
  df_flatten <- map(list_list, ~.[str_detect(names(.), name_pattern)]) |> 
    #Bind all data frames by row
    map_df(~bind_rows(.)) |> 
    #Rearrange data
    pivot_wider(names_from = variable_name, 
                values_from = all_of(column_values)) |> 
    clean_names()
  
  if(sum(str_detect(df_flatten$region, "FAO")) == 0){
    depth <- df_flatten |> 
      drop_na(deptho) |> 
      select(deptho, region, lat, lon)
    
    df_flatten <- df_flatten |> 
      filter(is.na(deptho)) |> 
      select(!deptho) |> 
      left_join(depth, by = join_by(lat, lon, region))
  }
  
  return(df_flatten)
}


# Extracting data from GFDL outputs ---------------------------------------
calc_inputs_all <- function(file_path_ctrl, file_path_obs, region_choice,
                            out_path_ctrl, out_path_obs){
  # This function extracts depth and climate variable files and applies it to 
  # the region of interest (LME/FAO), and saves a csv with depth and the climate
  # variable as columns. Nothing is returned, instead results are saved to disk.
  # Inputs:
  # file_path_ctrl (character) - Full file path for control data
  # file_path_obs (character) - Full file path for observations data
  # region_choice (numeric) - A number identifying the FAO or LME of interest
  # out_path_ctrl (character) - Path to folder where ctrlclim data will be saved
  # out_path_obs (character) - Path to folder where obsclim data will be saved
  #
  # Outputs:
  # None - This function saves outputs in specified paths
  
  message("Processing region: ", region_choice)
  
  #Search ctrlclim and obsclim files available for the region of interest
  obs_files <- list.files(file_path_obs, 
                          pattern = paste0("(-|_)LME(-|_)", region_choice,
                                           "_"), full.names = T)
  ctrl_files <- list.files(file_path_ctrl, 
                           pattern = paste0("(-|_)LME(-|_)", region_choice, 
                                            "_"), full.names = T) 
  
  #Apply calc_inputs function to all files in list
  obs_ctrl_data <- map2(ctrl_files, obs_files, calc_inputs)
  
  #Divide list into ctrlclim and obsclim data frames 
  output_obs_all_variables <- flatten_list_list(obs_ctrl_data, "obs_final",
                                                "weighted_mean")
  output_ctrl_all_variables <- flatten_list_list(obs_ctrl_data, "ctrl_final",
                                                 "weighted_mean")
  
  #Check if output folder provided exist, if not, create them
  if(!dir.exists(out_path_ctrl)){
    dir.create(out_path_ctrl, recursive = T)
  }
  if(!dir.exists(out_path_obs)){
    dir.create(out_path_obs, recursive = T)
  }
  
  reg <- output_ctrl_all_variables |> distinct(region) |> pull()
  
  #Create output file name
  out_file_obsclim <- file.path(out_path_obs, paste0("obsclim_spinup_", reg, 
                                                     "_all_inputs.csv")) 
  out_file_ctrlclim <- file.path(out_path_ctrl, paste0("ctrlclim_spinup_", 
                                                       reg, "_all_inputs.csv")) 
  
  #Saving data frames containing all data: spinup + obsclim/ctrlclim
  fwrite(x = output_obs_all_variables, file = out_file_obsclim)
  fwrite(x = output_ctrl_all_variables, file = out_file_ctrlclim)
}


# Processing obsclim/ctrlclim GFDL data -----------------------------------
gridded_inputs <- function(ctrl_files, obs_files){
  # This function loads gridded data and reorganises it for further analyses.
  #
  #Inputs:
  #ctrl_files (character) - Full file path for control data
  #obs_files (character) - Full file path for observations data
  #region (character) - Region of interest
  #Outputs:
  #(list) - Containing reorganised data for ctrlclim and obsclim GFDL outputs
  
  #Getting file name for control and observations datasets
  file_name_obs <- basename(obs_files)
  file_name_ctrl <- basename(ctrl_files)
  
  #Extracting variable name from observations filename
  variable <- str_extract(file_name_obs, "clim_(.+)_\\d{2}arcmin",
                          group = 1) |> 
    #If unit is present, it will be disregarded
    str_split_i("_", i = 1) 
  
  #Getting name of region
  region <- str_extract(file_name_ctrl, "arcmin_(.*)_monthly", group = 1)
  
  #Loading grid in the resolution matching the obsclim and ctrlclim data
  area_frame <- list.files(file.path("/g/data/vf71/shared_resources",
                                     "grid_cell_area_ESMs/isimip3a"),
                           #Get resolution of data from file name
                           paste0(str_extract(file_name_ctrl, 
                                              "_\\d{2}arcmin_"), ".*csv$"),
                           full.names = T) |> 
    read_csv() |> 
    #Rename columns
    rename("lon" = "x", "lat" = "y", "area_m2" = "cellareao")
  
  
  if(variable == "deptho"){
    obsclim_df <- fread(obs_files) |> 
      select(!area_m2) |> 
      rename(value = m) |> 
      mutate(date = NA, .before = value) |> 
      #Add correct area of the grid cell
      left_join(area_frame, by = join_by(lat, lon)) |> 
      mutate(scenario = str_extract(obs_files, 
                                    "(fao|lme)_inputs/(.*)/\\d{1,3}deg", 
                                    group = 2),
             month = NA, year = NA, region = region, variable_name = "deptho")
    ctrlclim_df <- fread(ctrl_files) |> 
      select(!area_m2) |> 
      rename(value = m) |> 
      mutate(date = NA, .before = value) |> 
      #Add correct area of the grid cell
      left_join(area_frame, by = join_by(lat, lon)) |> 
      mutate(scenario = str_extract(ctrl_files, 
                                    "(fao|lme)_inputs/(.*)/\\d{1,3}deg", 
                                    group = 2),
             month = NA, year = NA, region = region, variable_name = "deptho")
  }else{
    #Load ctrlclim data
    ctrlclim_df <- reorganise_gfdl(ctrl_files, area_frame, region) |> 
      # add name of variable to data frame
      mutate(variable_name = variable)
    
    #Load obsclim data
    obsclim_df <- reorganise_gfdl(obs_files, area_frame, region) |> 
      # add name of variable to data frame
      mutate(variable_name = variable)
  }
  
  return(list(obsclim_df = obsclim_df, 
              ctrlclim_df = ctrlclim_df))
}


dbpm_calcs <- function(flatten_list){
  # This function calculates DBPM inputs from flatten lists containing all
  # input variables from GFDL
  #
  #Inputs:
  #flatten_list (data frame) - Data frame containing all DBPM input variables
  #from the gridded_inputs() function
  #obs_files (character) - Full file path for observations data
  #region (character) - Region of interest
  #Outputs:
  #df (data frame) - Containing reorganised gridded GFDL data
  
  #Check if a "depth" column is included in the data frame
  if(sum(grepl("depth", names(flatten_list))) == 0){
    #If "depth" does not exist, add it
    flatten_list <- flatten_list |> 
      mutate(depth = 200)
  }else{
    #Ensure "depth" column in correctly labelled
    names(flatten_list) <- str_replace(
      names(flatten_list), "deptho_m|deptho", "depth")
  }
  
  #Calculate DBPM inputs
  df <- flatten_list |> 
    mutate(sphy = phypico_vint,  lphy = phyc_vint - phypico_vint) |> 
    #Removing columns not needed
    select(-c(phyc_vint, phypico_vint)) |> 
    # name columns as in "dbpm_model_functions.R" script
    rename(t = date, sbt = tob, sst = tos, 
           expcbot = expc_bot) |> 
    #Calculate slope and intercept
    mutate(er = getExportRatio(sphy, lphy, sst, depth),
           er = ifelse(er < 0, 0, ifelse(er > 1, 1, er)),
           intercept = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                     mmid = 10^-10.184,
                                     mmax = 10^-5.25, depth, 
                                     output = "intercept"),
           slope = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                 mmid = 10^-10.184, mmax = 10^-5.25, depth, 
                                 output = "slope")) |> 
    #Reorganise columns following original script
    relocate("region", .before = t) |> 
    relocate(all_of(c("sst", "sbt", "er", "intercept", "slope", "sphy", "lphy",
                      "depth", "area_m2", "expcbot", "year", "month")), 
             .after = t) |> 
    #Keep up to 3 decimal places
    mutate(across(sst:slope, ~round(.x, 3)),
           across(sphy:lphy, ~round(.x, 5)))
  
  return(df)
}


# Extracting data from GFDL outputs in gridded form -----------------------
calc_inputs_gridded <- function(file_path_ctrl, file_path_obs, region_choice,
                                out_path_ctrl, out_path_obs){
  # This function extracts climate variable data for a region of interest 
  # (LME/FAO), and saves a csv with depth and the climate
  # variable as columns. Nothing is returned, instead results are saved to disk.
  # Inputs:
  # file_path_ctrl (character) - Full file path for control data
  # file_path_obs (character) - Full file path for observations data
  # region_choice (numeric) - A number identifying the FAO or LME of interest
  # out_path_ctrl (character) - Path to folder where ctrlclim data will be saved
  # out_path_obs (character) - Path to folder where obsclim data will be saved
  #
  # Outputs:
  # None - This function saves outputs in specified paths
  
  message("Processing region: ", region_choice)
  
  #Search ctrlclim and obsclim files available for the region of interest
  obs_files <- list.files(file_path_obs, 
                          pattern = paste0("(-|_)LME(-|_)", region_choice,
                                           "_"), full.names = T)
  ctrl_files <- list.files(file_path_ctrl, 
                           pattern = paste0("(-|_)LME(-|_)", region_choice, 
                                            "_"), full.names = T) 
  
  #Apply calc_inputs function to all files in list
  obs_ctrl_data <- map2(ctrl_files, obs_files, gridded_inputs)
  
  #Divide list into ctrlclim and obsclim data frames 
  output_obs_all_variables <- flatten_list_list(obs_ctrl_data, "obsclim",
                                                "value") |> 
    dbpm_calcs()
  
  output_ctrl_all_variables <- flatten_list_list(obs_ctrl_data, "ctrlclim",
                                                 "value") |> 
    dbpm_calcs()
  
  # Calculate spinup: 6 cycles of 20 years from ctrlclim data
  spinup <- output_ctrl_all_variables |> 
    arrange(lat, lon, t) |> 
    filter(year >= 1961 & year <= 1980)
  
  spinup <- spinup |> 
    #Repeat six times
    slice(rep(1:n(), times = 6)) |>
    #Relabel dates
    mutate(i = rep(1:6, each = nrow(spinup)),
           year = year-(120-(20*(i-1))),
           t = my(paste(month, year, sep = "-")),
           scenario = "spinup") |>
    select(!i)
  
  
  #Check if output folder provided exist, if not, create them
  if(!dir.exists(out_path_ctrl)){
    dir.create(out_path_ctrl, recursive = T)
  }
  if(!dir.exists(out_path_obs)){
    dir.create(out_path_obs, recursive = T)
  }
  
  #Getting name of region
  reg <- output_ctrl_all_variables |> distinct(region) |> pull()
  
  #Create output file name
  out_file_obsclim <- file.path(out_path_obs, paste0("obsclim_historical_", 
                                                     reg, "_all_inputs.csv")) 
  out_file_ctrlclim <- file.path(out_path_ctrl, paste0("ctrlclim_historical_", 
                                                       reg, "_all_inputs.csv"))
  out_file_spinup <- file.path(out_path_ctrl, paste0("ctrlclim_spinup_", 
                                                     reg, "_all_inputs.csv"))
  
  #Saving data frames containing all data: spinup + obsclim/ctrlclim
  fwrite(x = output_obs_all_variables, file = out_file_obsclim)
  fwrite(x = output_ctrl_all_variables, file = out_file_ctrlclim)
  fwrite(x = spinup, file = out_file_spinup)
}


# Transform spatRaster to data frame --------------------------------------
ras_to_df <- function(ras){
  # This function transforms a spatRaster object to a data frame. It uses the
  # variable name in the raster to label the data frame column containing
  # the grid cell values
  # Inputs:
  # ras (spatRaster) - This raster must have the variable name recorded
  #
  # Outputs:
  # ras_df (data frame) - Contains grid cell values of spatRaster as a data 
  # frame
  
  # Extracting raster layer name
  variable <- varnames(ras)
  
  #Keeping index for layer name
  names(ras) <- str_extract(names(ras), "\\d{1,3}")
  
  ras_df <- ras |> 
    #Keeping coordinate values
    as.data.frame(xy = T) |> 
    pivot_longer(!x:y, names_to = "index", values_to = variable) |> 
    #Changes indices from characters to integers
    mutate(index = as.integer(index)) |> 
    #Removes any special characters
    clean_names() |> 
    arrange(x, y)
  
  return(ras_df)
}



# Merge data frames if coordinates and time are equal ---------------------
merge_equal <- function(ras_df1, ras_df2){
  #If coordinates (labelled x, y) and index are the same across data frames, 
  #then add variable to first dataset
  # Inputs:
  # ras_df1 (data frame) - Main data frame to which data will be added
  # ras_df2 (data frame) - Secondary data frame from which data comes from
  #
  # Outputs:
  # ras_df (data frame) - Contains merged data frame if coordinates and time
  # are identical
  if(identical(ras_df1 |> select(x:index),
               ras_df2 |> select(x:index))){
    ras_df <- ras_df1 |> 
      bind_cols(ras_df2 |> select(!x:index))
    
  return(ras_df)
  }else{
    message("coordinates and time are not identical between data frames.")
  }
}


# Saving gridded ESM data as data frames ----------------------------------
getGCM <- function(folder_path, save_path, getdepth = T){
  # This function transforms gridded ESM data into into csv files and saves 
  # them to the specified path
  #
  # Inputs:
  # folder_path (character) - Full file path where gridded ESM data is located
  # save_path (character) - Full file path where transformed ESM data will be
  # saved
  # getdepth (boolean) - If TRUE, then depth variable is loaded and saved as 
  # data frame. Default is TRUE. 
  #
  # Outputs:
  # None - This function saves outputs in specified paths
  
  #Ensure save_path exist, otherwise create it
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = T)
  }
  
  #Get file paths for variables of interest
  files_int <- list.files(folder_path, full.names = TRUE) |> 
    str_subset("phypico-vint|phyc-vint|tob|tos")
  
  #Getting range of years included in dataset from file name
  yr_range <- str_extract(files_int[1], "\\d{4}_\\d{4}") |> 
    str_split("_", simplify = T)
  
  #Create sequence of dates for each month between year range
  #To be used to correct dates in data frame
  yr_range <- data.frame(t = seq(from = ym(paste0(yr_range[1], "-01")),
                                 to = ym(paste0(yr_range[2], "-12")),
                                 by = "month")) |> 
    rowid_to_column("index")
  
  #Get scenario
  scenario <- str_extract(folder_path, "global_inputs/(.*)/\\d{1,3}deg", 
                          group = 1)
  
  #Calculating large phytoplankton (lphy)
  phyc_ras <- rast(str_subset(files_int, "phyc-"))
  phypico_ras <- rast(str_subset(files_int, "phypico-"))
  lphy_ras <- phyc_ras-phypico_ras
  
  #Removing raster that is not needed
  rm(phyc_ras)
  
  #Updating variable name
  varnames(lphy_ras) <- "lphy"
  varnames(phypico_ras) <- "sphy"
  
  #Transforming to data frame
  lphy_df <- ras_to_df(lphy_ras)
  sphy_df <- ras_to_df(phypico_ras)
  
  #Removing rasters that are not needed
  rm(phypico_ras, lphy_ras)
  
  #Loading TOB and TOS
  to_zb_df <- rast(str_subset(files_int, "tob_")) |> 
    ras_to_df()
  to_zs_df <- rast(str_subset(files_int, "tos_")) |> 
    ras_to_df()
  
  #Merge datasets if they are equal
  pp <- merge_equal(lphy_df, sphy_df) |> 
    merge_equal(to_zb_df) |> 
    merge_equal(to_zs_df)
  
  #Removing data frames that are not needed
  rm(lphy_df, sphy_df, to_zb_df, to_zs_df)
  
  #Adding dates
  pp <- pp |> 
    left_join(yr_range, join_by(index)) |> 
    select(!index) |> 
    #Add scenario
    mutate(scenario = scenario) |> 
    #Renaming columns
    rename("lon" = "x", "lat" = "y", "sbt" = "tob", "sst" = "tos") |> 
    relocate(t, .after = lat)
  
  #Save outputs
  file_out <- basename(str_subset(files_int, "tob")) |> 
    str_replace("tob", "all_forcings") |>
    str_replace(".nc", ".csv")
  
  pp |> 
    fwrite(file.path(save_path, file_out))
  
  #Saving depth 
  if(getdepth == T){
    #Get path for depth data 
    depth_filename <- list.files(folder_path, pattern = "deptho", 
                                 full.names = T)
    #Load raster
    depth_df <- rast(depth_filename) |> 
      #Transform to data frame
      as.data.frame(xy = T, na.rm = F) |> 
      rename("lon" = "x", "lat" = "y", "depth" = "deptho") |> 
      #Add grid cell ID
      rowid_to_column("gridnum") |> 
      #Remove any grid cells without data
      drop_na() |> 
      #Add scenario
      mutate(scenario = scenario) |> 
      relocate(gridnum, .after = depth)
    
    #Save depth data frame
    depth_save_name <- file.path(save_path,
                                 str_replace(basename(depth_filename), 
                                             ".nc", ".csv"))
    depth_df |> 
      fwrite(depth_save_name)
  }
}


# Calculate spinup from gridded ESM data ----------------------------------
calc_input_spinup_gridcell <- function(base_path, save_path){
  # This function uses gridded ESM data, creates spinup data and saves 
  # them to the specified path
  #
  # Inputs:
  # base_path (character) - Full file path where gridded ESM data is located
  # save_path (character) - Full file path where spinup data will be saved
  #
  # Outputs:
  # None - This function saves outputs in specified paths
  
  #Ensure save_path exist, otherwise create it
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = T)
  }
  
  #Get path for file containing all model forcings
  file_path <- list.files(base_path, "all_forcings_", full.names = T)
  
  #Load file with all forcings
  ctrlclim <- fread(file_path)
  
  #Load file with depth
  gridnum_depth <- list.files(base_path, "deptho", full.names = T) |> 
    fread() |> 
    select(!scenario)
  
  #Merge files
  ctrlclim <- ctrlclim |>
    full_join(gridnum_depth, join_by(lon, lat)) |>
    arrange(t, gridnum)
  
  #Remove variable - not needed
  rm(gridnum_depth)
  
  #Calculate spin-up from data between 1961 and 1980
  spinup <- ctrlclim |>
    filter(t >= "1961-01-01", t <= "1980-12-01") |>
    mutate(year = year(t),
           month = month(t),
           scenario = "spinup") 
  
  #Remove variable - not needed
  rm(ctrlclim)
  
  #Calculating year for spinup
  year_df <- spinup |> 
    select(year, month) |> 
    slice(rep(1:n(), times = 6)) |> 
    #Relabel dates
    mutate(i = rep(1:6, each = nrow(spinup)),
           year = year-(120-(20*(i-1))),
           t = my(paste(month, year, sep = "-"))) |> 
    pull(t)
  
  spinup <- spinup |>
    select(!c(year, month)) |> 
    #Repeat six times
    slice(rep(1:n(), times = 6)) |> 
    #Relabel dates
    mutate(t = year_df) 
  
  rm(year_df)
  
  #Create file path to save spinup
  file_out <- file.path(save_path, str_replace(basename(file_path), 
                                               "all_forcings", "spinup"))
  #save spin up
  spinup |> 
    fwrite(file_out)
}
