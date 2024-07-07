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
    mutate(date = as_date(date)) |> 
    #Add correct area of the grid cell
    left_join(area_cell_df, by = join_by(lat, lon)) |> 
    #Adding scenario to data
    mutate(scenario = str_extract(gfdl_path, 
                                  "fao_inputs/(.*)/\\d{2,3}deg", group = 1)) |> 
    #Remove any empty rows
    drop_na() |> 
    #Apply date format to date column
    mutate(month = lubridate::month(date, label = T), 
           year = year(date), region = region) |> 
    #Order by date
    arrange(date)
  
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
    # rename weighted_mean column according to variable 
    mutate(variable_name = variable)
  
  return(df_spinup)
}


# Calculating weighted means for GFDL data --------------------------------
calc_inputs <- function(path_crtl, path_obs){
  #This function calculates fixed weighted.mean for depth and total area, 
  # and monthly weighted.mean for each climate variable. It also creates and 
  # adds spin up information from ctrlclim data to GFDL variables. 
  #
  #Inputs:
  #path_crtl (character) - Full file path for control data
  #path_obs (character) - Full file path for observations data
  #region (character) - Region of interest
  #Outputs:
  #(list) - Containing weighted mean (by area of grid cell) for ctrlclim and
  #obsclim GFDL outputs
  
  #Getting file name for control and observations datasets
  file_name_obs <- basename(path_obs)
  file_name_crtl <- basename(path_crtl)
  
  #Extracting variable name from observations filename
  variable <- str_extract(file_name_obs, "obsclim_(.+)_15arc",  group = 1) |> 
    #If unit is present, it will be disregarded
    str_split_i("_", i = 1) 
  
  #Getting name of region
  region <- str_extract(file_name_crtl, "arcmin_(.*)_monthly", 
                       group = 1)
  
  #Loading grid in the resolution matching the obsclim and ctrlclim data
  area_frame <- list.files(file.path("/rd/gem/private/shared_resources", 
                                     "grid_cell_area_ESMs/isimip3a"),
                           #Get resolution of data from file name
                           paste0(str_extract(file_name_crtl, 
                                              "_\\d{2}arcmin_"), ".*csv$"),
                           full.names = T) |> 
    read_csv() |> 
    #Rename columns
    rename("lon" = "x", "lat" = "y", "area_m2" = "cellareao")
  
  #Load ctrlclim data
  ctrlclim_df <- reorganise_gfdl(path_crtl, area_frame, region)
  
  #Load obsclim data
  obsclim_df <- reorganise_gfdl(path_obs, area_frame, region)
  
  # Calculate weighted mean for GFDL variables
  if(variable == "deptho_m"){
    # calculate fixed variables - mean depth and area of FAO 
    # shortcut if you need to extract ctrlclim values too for all variables
    weighted_mean_crtl_final <- weighted_mean_obs_final <- obsclim_df |> 
      summarise(deptho_m = weighted.mean(value, area_m2),
                total_area_km2 = sum(area_m2)*1e-6) |> 
      # Add name of region from file path
      mutate(region = region)
  }else{
    # ctrlclim
    weighted_mean_crtl <- weighted_mean_by_area(ctrlclim_df)
    # obsclim 
    weighted_mean_obs <- weighted_mean_by_area(obsclim_df)
    
    # SPINUP, based on ctrlclim and used for both ctrlclim and obsclim
    spinup <- weighted_mean_crtl |>
      filter(year >= 1961, year <= 1980) |> 
      slice(rep(1:n(), times = 6)) |>
      mutate(year = rep(1841:1960, each = 12),
             date = my(paste(month, year, sep = "-")),
             scenario = "spinup") 
    
    # add spinup to control
    weighted_mean_crtl_final <- add_spinup_data(spinup, weighted_mean_crtl, 
                                                variable)
    
    # add spin up to obsclim
    weighted_mean_obs_final <- add_spinup_data(spinup, weighted_mean_obs, 
                                               variable)
  }
  return(list(weighted_mean_obs_final = weighted_mean_obs_final, 
              weighted_mean_crtl_final = weighted_mean_crtl_final))
}



# Flattening list of lists based on sublist names -------------------------
flatten_list_list <- function(list_list, name_pattern){
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
    pivot_wider(names_from = variable_name, values_from = weighted_mean) |> 
    clean_names()
  return(df_flatten)
}


# Extracting data from GFDL outputs ---------------------------------------
calc_inputs_all <- function(file_path_crtl, file_path_obs, region_choice,
                            out_path_ctrl, out_path_obs){
  # This function extracts depth and climate variable files and applies it to 
  # the region of interest (LME/FAO), and saves a csv with depth and the climate
  # variable as columns. Nothing is returned, instead results are saved to disk.
  # Inputs:
  # file_path_crtl (character) - Full file path for control data
  # file_path_obs (character) - Full file path for observations data
  # region_choice (numeric) - A number identifying the FAO or LME of interest
  # out_path_ctrl (character) - Path to folder where ctrlclim data will be saved
  # out_path_obs (character) - Path to folder where obsclim data will be saved
  #
  # Outputs:
  # None - This function saves outputs in specified paths
  
  #Creating file name pattern to search relevant data files
  file_pattern <- paste0("FAO-LME-", region_choice, "_")
  
  #Search ctrlclim and obsclim files available for the region of interest
  obs_files <- list.files(file_path_obs, pattern = file_pattern, 
                          full.names = T)
  ctrl_files <- list.files(file_path_crtl, pattern = file_pattern, 
                           full.names = T) 
  
  #Apply calc_inputs function to all files in list
  obs_ctrl_data <- map2(ctrl_files, obs_files, calc_inputs)
  
  #Divide list into ctrlclim and obsclim data frames 
  output_obs_all_variables <- flatten_list_list(obs_ctrl_data, "obs_final")
  output_ctrl_all_variables <- flatten_list_list(obs_ctrl_data, "crtl_final")
  
  #Check if output folder provided exist, if not, create them
  if(!dir.exists(out_path_ctrl)){
    dir.create(out_path_ctrl)
  }
  if(!dir.exists(out_path_obs)){
    dir.create(out_path_obs)
  }
  
  #Create output file name
  out_file <- file.path(out_path_obs, paste0("obsclim_spinup", 
                                             file_pattern, "all_inputs.csv")) 
  this_destination_path_obs <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a", 
                                      "/processed_forcings/fao_inputs/obsclim/",
                                      "0.25deg/observed_FAO_", region_choice, ".csv")
  this_destination_path_ctrl <- paste0("/rd/gem/private/fishmip_inputs/", 
                                       "ISIMIP3a/processed_forcings/fao_inputs",
                                       "/ctrlclim/0.25deg/control_FAO_", 
                                       region_choice, ".csv")
  
  fwrite(x = output_obs_all_variables, 
         file = file.path(this_destination_path_obs))
  fwrite(x = output_crtl_all_variables, 
         file = file.path(this_destination_path_ctrl))
  
}



