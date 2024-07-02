
rm(list=ls())

# ------- STEP 1: GET GCM INPUTS FOR DYNAMIC BENTHIC-PELAGIC SIZE SPECTRUM MODEL

# Ryan ran this step for ISIMIP3b and provided the inputs. Inputs are provided 
# at the global 1 deg grid cell level. 
#### Ryan is running this step for ISIMIP3a, 3 scales: 
# 1 FAO aggregated - for calibration, using 0.25deg inputs only and focusing on
# obsclim
# 2 FAO grid cell level at 1 deg and 0.25 deg - for model running, need both 
# inputs and both scales 

## For testing, directory on Ryan's machine
#setwd("/Users/rheneghan/Desktop/Papers/DBPM_Attribution/fao_inputs")

#### ISIMIP3a scale 1 -----
# Use FAO inputs provided by Denisse in September 2023
# NOTE: we are using only 0.25 deg inputs as these are best when calculating 
#means across LMEs. 

# library(raster)
library(terra)
# library(stringr)
library(tidyverse)
library(tictoc)
library(data.table)
library(parallel)
library(dtplyr)

#### 1. Define main function: for each LME, and each variable, apply the 
# calc_inputs_LME function -----

# 1A. calc_inputs_FAO() (USING FOR FAO REGIONS, SEPTEMBER 2023)
# this function calculates FAO fixed weighted.mean depth and total area, 
# and monthly weighted.mean for each climate variable. 
# it then add spin up from control to the cliamte variables.  

calc_inputs_FAO <- function(file_name_obs, file_name_crtl, file_path_crtl, 
                            file_path_obs){
  
  # # trial, deleted by RFH to tidy code (Sep 2023) 
  
  # RFH trial
  # file_path_obs <- "./climate/obsclim/"
  #   file_name_obs <- "gfdl-mom6-cobalt2_obsclim_expc-bot_15arcmin_FAO-LME-21_monthly_1961_2010.csv"
  # file_path_crtl <- "./climate/ctrlclim/"
  # file_name_crtl <- "gfdl-mom6-cobalt2_ctrlclim_expc-bot_15arcmin_FAO-LME-21_monthly_1961_2010.csv"
  
  # extract variable name
  variable <- str_match(file_name_obs, 
                        "gfdl-mom6-cobalt2_obsclim_\\s*(.*?)\\s*_15arcmin")[2]
  
  ##############################################################################
  #### RESHAPE DATA (so FAO data conforms to LME data shape, RFH Sep. 2023) ####
  ##############################################################################
  
  # work on CONTROLCLIM first, OPEN AND RESHAPE DATA 
  # dataframe of areas for 0.25degree, global
  # area_frame <- data.frame("lon" = rep(seq(-179.875, 179.875, by = 0.25),
  #                                      by = 720), 
  #                          "lat" = rep(seq(-89.875, 89.875, by = 0.25), 
  #                                      each = 1440), 
                           # "area_m2" = 1e6*as.vector(t(as.matrix(area(raster(nrows = 720, 
                           #                                                   ncol = 1440))))))
  
  #CHECK: May need to replace this with new masks
  #Create spat raster object - 0.25 deg resolution
  area_frame <- rast(resolution = 0.25, nrows = 720, ncol = 1440) |> 
    #Calculate area of grid cell 
    cellSize() |> 
    #Transform to data frame
    as.data.frame(xy = T) |>
    #Rename columns
    rename("lon" = "x", "lat" = "y", "area_m2" = "area")
  
  fao_crtl_raw <- fread(file.path(file_path_crtl, file_name_crtl), 
                        header = FALSE)
  
  #CHECK: Need to verify what is being removed and replace hard coded row
  dat1 <- fao_crtl_raw[-3, ]
  dat_names <- unlist(as.vector(dat1[, 1]))
  all_dat <- t(as.matrix(dat1[, -1]))
  
  dates <- strsplit(dat_names[-c(1:2)], "-")
  years <- unlist(lapply(dates, '[[', 1))
  months <- lapply(dates, function(x){month.abb[as.numeric(x[2])]})
  col_names <- c("lat", "lon", paste(months, years, sep = "_"))
  
  fao_crtl <- as.data.frame(all_dat)
  colnames(fao_crtl) <- col_names
  
  # fao_crtl <- left_join(fao_crtl, area_frame, 
  #                       by = c("lat" = "lat", "lon" = "lon"))
  # Some boundary cells are NAs, remove (very small number)
  # fao_crtl <- fao_crtl |> na.omit() 
  
  #Adding area of grid cell
  fao_crtl <- fao_crtl |> 
    left_join(area_frame, by = join_by(lat, lon)) |> 
    #Removing empty cells
    drop_na()
  
  
  # then work on OBSERVED, OPEN AND RESHAPE DATA
  fao_raw <- fread(file.path(file_path_obs, file_name_obs), header = FALSE)
  
  #CHECK: Need to verify what is being removed and replace hard coded row
  dat1 <- fao_raw[-3,]
  dat_names <- unlist(as.vector(dat1[,1]))
  all_dat <- t(as.matrix(dat1[,-1]))
  
  dates <- strsplit(dat_names[-c(1:2)], "-")
  years <- unlist(lapply(dates, '[[', 1))
  months <- lapply(dates, function(x){month.abb[as.numeric(x[2])]})
  col_names <- c("lat", "lon", paste(months, years, sep = "_"))
  
  fao <- as.data.frame(all_dat)
  colnames(fao) <- col_names
  
  # fao <- left_join(fao, area_frame, by = c("lat" = "lat", "lon" = "lon"))
  # fao <- fao |> na.omit()
  
  #Adding area of grid cell
  fao <- fao |> 
    left_join(area_frame, by = join_by(lat, lon)) |> 
    #Removing any empty rows
    drop_na()
  
  ## Clear what we don't need
  rm("fao_crtl_raw", "fao_raw", "dat1", "dates", "all_dat", "years", "months", 
     "col_names")
  
  ##############################################################################
  ##############################################################################
  
  
  # 2 methods tested to calculate weighed mean
  # METHOD 1 based on dplyr and cell area
  # METHOD 2 based on raster and cos(lat) - NEEDS TO BE UPDATED
  # for lme 1, results are the same for depth (fixed variable) 
  # expc-bot_mol_m-2_s-1 (first of inputs) across methods (checked trends). 
  # We assume this is the case for all LME and inputs 
  # method 1 is faster, method 2 allow easy LME plotting 
  
  # METHOD 1
  if (variable == "deptho_m"){
    # calculate fixed variables - mean depth and area of FAO 
    weighted_mean_obs <- fao |> 
      summarise(deptho_m = weighted.mean(m, area_m2),
                area_km2 = sum(area_m2)/1e6) |> 
      ## Change here! From LME to FAO, RFH Sep. 2023
      mutate(FAO = str_extract(file_name_obs, "(?<=FAO-LME-).+(?=_monthly)")) 
    
    weighted_mean_obs_final <- weighted_mean_obs
    # shortcut if you need to extract ctrlclim values too for all variables
    weighted_mean_crtl_final <- weighted_mean_obs 
    
  }else{
    
    # CONTROL 
    weighted_mean_crtl <- fao_crtl |> 
      # gather(key = "Date",
      #        value = "value",
      #        -c("lat", "lon","area_m2"))|> 
      pivot_longer(-c("lat", "lon", "area_m2"), names_to = "Date", 
                   values_to = "value") |> 
      group_by(Date)   |> 
      summarise(weighted_mean_crtl = weighted.mean(value, area_m2),
                area_km2 = sum(area_m2)/1e6) |> 
      ungroup()  |> 
      ## Change here! From LME to FAO, RFH Sep. 2023
      mutate(FAO = str_extract(file_name_crtl, "(?<=FAO-LME-).+(?=_monthly)"), 
             Month = str_extract(Date, "[[:upper:]]+[[:lower:]]+"),
             Year = str_extract(Date, "\\d+"), 
             Date = lubridate::my(paste(Month, Year, sep = "." ))) |> 
      arrange(Date)
    
    # to plot METHOD 1 and compare with METHOD 2 -->>  
    # DELETED BY RFH (SEP 2023) TO TIDY CODE 
    
    # SPINUP, based on ctrlclim and used for both ctrlclim and observed
    spinup <- weighted_mean_crtl |>
      dplyr::select(-Date) |> 
      filter(Year >= 1961, Year <= 1980) |> 
      slice(rep(1:n(), times = 6)) |>
      mutate(Year = as.character(rep(1841:1960, each = 12)),
             Date = lubridate::my(paste(Month,Year, sep = "." ))) 
    
    # add spinup to control and check that's all OK 
    weighted_mean_crtl_final <- weighted_mean_crtl |> 
      full_join(spinup) |> 
      arrange(Date) |> 
      # reorder columns 
      relocate(any_of(c("FAO", "Date", "Year", "Month", "weighted_mean_crtl", 
                        "area_km2"))) |> 
      # rename weighted_mean column according to variable 
      rename(variable = "weighted_mean_crtl")
    
    # reorder columns 
    # weighted_mean_crtl_final <- weighted_mean_crtl_final[,c("FAO", "Date", 
    #                                                         "Year", "Month",
    #                                                         "weighted_mean_crtl",
    #                                                         "area_km2")]
    
    # rename weighted_mean column according to variable 
    # names(weighted_mean_crtl_final)[5] <- variable
    
    # OBSERVED 
    weighted_mean_obs<-fao |> 
      # gather(key = "Date",
      #        value = "value",
      #        -c("lat", "lon","area_m2")) |>
      pivot_longer(-c("lat", "lon", "area_m2"), names_to = "Date", 
                   values_to = "value") |> 
      group_by(Date) |> 
      summarise(weighted_mean_obs = weighted.mean(value, area_m2),
                area_km2 = sum(area_m2)/1e6) |> 
      ungroup() |> 
      ## Change here! From LME to FAO, RFH Sep. 2023
      mutate(FAO = str_extract(file_name_crtl, "(?<=FAO-LME-).+(?=_monthly)"), 
             Month = str_extract(Date, "[[:upper:]]+[[:lower:]]+"),
             Year = str_extract(Date, "\\d+"), 
             Date = lubridate::my(paste(Month, Year, sep = "." ))) |> 
      arrange(Date)
    
    # add spin up to observed and plot to check 
    spinup <- spinup |> 
      rename(weighted_mean_obs = weighted_mean_crtl)
    
    weighted_mean_obs_final <- weighted_mean_obs |> 
      full_join(spinup) |> 
      arrange(Date) |> 
      # reorder columns 
      relocate(any_of(c("FAO", "Date", "Year", "Month", "weighted_mean_obs", 
                        "area_km2"))) |> 
      # rename weighted_mean column according to variable 
      rename(variable = "weighted_mean_obs")
    
    # reorder columns 
    ## Change here! From LME to FAO, RFH Sep. 2023
    # weighted_mean_obs_final<-weighted_mean_obs_final[,c("FAO", "Date", "Year",
    #                                                     "Month",
    #                                                     "weighted_mean_obs", 
    #                                                     "area_km2")]
    
    # rename weighed_mean column according to variable 
    # names(weighted_mean_obs_final)[5] <- variable
    
  }
  
  # # METHOD 2, DELETED BY RFH (SEP 2023) TO TIDY CODE
  
  return(list(weighted_mean_obs_final = weighted_mean_obs_final, 
              weighted_mean_crtl_final = weighted_mean_crtl_final)) 
  
}

# 1B. calc_inputs_all_FAO()
# this function extracts depth and climate variable files
# applies calc_inputs_FAO() to each file 
# and saves a csv with depth and all climate variables as columns for one LME
calc_inputs_all_FAO<-function(this_FAO){
  
  # # trial 
  # this_LME = 2
  
  file_path_obs <- "/rd/gem/private/fishmip_inputs/ISIMIP3a/fao_inputs/obsclim/0.25deg"
  file_path_crtl <- "/rd/gem/private/fishmip_inputs/ISIMIP3a/fao_inputs/ctrlclim/0.25deg"
  
  this_FAO_new <- paste0("FAO-LME-", this_FAO, "_")
  
  # this adds the whole path
  fao_obs <- list.files(file_path_obs, pattern = this_FAO_new) #, full.names = TRUE) 
  fao_ctrl <- list.files(file_path_crtl, pattern = this_FAO_new) #, full.names = TRUE) 
  
  tic()
  output_obs <- list()
  output_crtl <- list()
  
  # loop across inputs: depth, phyc tob, expt, phypico, tos. 
  for(i in 1:length(fao_obs)){
    a <- calc_inputs_FAO(file_name_obs = fao_obs[[i]], 
                         file_name_crtl = fao_ctrl[[i]], 
                         file_path_crtl = file_path_crtl, 
                         file_path_obs = file_path_obs)
    output_obs[[i]] <- a$weighted_mean_obs_final
    output_crtl[[i]] <- a$weighted_mean_crtl_final
  } 
  toc() 
  
  # all inputs together for one FAO
  output_obs_all_variables <- Reduce(merge, output_obs)
  output_crtl_all_variables <- Reduce(merge, output_crtl)
  
  # write output files - temporary path - need to save on gem48! DONE
  # this_destination_path_obs <- paste0("/data/home/camillan/dbpm/Output/", 
  #                                     "observed_LME_", this_LME, ".csv")
  # this_destination_path_ctrl <- paste0("/data/home/camillan/dbpm/Output/",
  #                                      "control_LME_", this_LME, ".csv")
  this_destination_path_obs <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a", 
                                      "/processed_forcings/fao_inputs/obsclim/",
                                      "0.25deg/observed_FAO_", this_FAO, ".csv")
  this_destination_path_ctrl <- paste0("/rd/gem/private/fishmip_inputs/", 
                                       "ISIMIP3a/processed_forcings/fao_inputs",
                                       "/ctrlclim/0.25deg/control_FAO_", 
                                       this_FAO, ".csv")
  
  fwrite(x = output_obs_all_variables, 
         file = file.path(this_destination_path_obs))
  fwrite(x = output_crtl_all_variables, 
         file = file.path(this_destination_path_ctrl))
  
}

#### 2. apply the functions above to each FAO region -----

this_FAO = c(77, 31, 41, 87, 57, 58, 71, 81, 21, 51, 34, 27, 47, 48, 61, 67, 88) 

tic()
for (i in 1:length(this_FAO)){
  message("Processing #", i, " of ", length(this_FAO))
  calc_inputs_all_FAO(this_FAO[[i]])
}
toc() 

# 3. read in printed csv file for each FAO and merge info into a unique file ----

# WARNING not there right now as were printed in DBPM gem48 instance and then
# only the final product was moved to new folder I think! 
# DONE - now re-run and printed in re/gem
newly_written_files_observed <- file.path("/rd/gem/private/fishmip_inputs",
                                          "ISIMIP3a/processed_forcings", 
                                          "fao_inputs/obsclim/0.25deg") |> 
  list.files(pattern = "observed", full.names = TRUE)

# pick one randomly and check 
# map(newly_written_files_observed[[8]], fread)

# combine files 
combined_FAO_inputs <- rbindlist(mclapply(X = newly_written_files_observed, 
                                          FUN = fread, mc.cores = 40))
FAO_area <- combined_FAO_inputs |> 
  group_by(FAO) |> 
  dplyr::select(FAO, area_km2) |> 
  distinct()

# head(combined_LME_inputs)
# sort(unique(combined_LME_inputs$LME)) # WARNING LME 0 missing. 

#### 4. check calculation below in terms of sphy and sphy and adopt same ----
# variable names

# function below: 
# names(pp) <- c("lon", "lat", "t", "lphy", "sphy", "sbt", "sst")
# depth file: lat, lon, depth 
# the 2 files are saved directly withing the function as RData  
# the function is applied to different protocols
# for CMIP63a protocols are observed and control at 0.25 deg and 1 deg 
# resolution. However, you are only using observed at LME scale for calibration 
# - so run only this scenario. 

# checked with Julia 1/09/2022
# sphy = phypico-vint_mol_m-2
# lphy = phyc-vint_mol_m-2 - phypico-vint_mol_m-2 

combined_FAO_inputs <- combined_FAO_inputs |> 
  mutate(sphy = `phypico-vint`,  lphy = `phyc-vint` - `phypico-vint`) |> 
  dplyr::select(-c(`phyc-vint`,`phypico-vint`, area_km2))

#### 5. add effort and catches ------
# also tried _adminCountry - all the same
# effort_old <- read_csv("/rd/gem/private/users/yannickr/DKRZ_EffortFiles/
#effort_histsoc_1841_2010_revised.csv") 
effort <- fread(file.path("/rd/gem/private/users/yannickr/DKRZ_EffortFiles", 
                          "effort_isimip3a_histsoc_1841_2010.csv"))
effort_FAO <- effort |> 
  dplyr::filter(LME == 0)

# ### check LME 1 
# effort_old<-effort_old |> 
#   filter(LME == 1) |> 
#   group_by(LME, Year) |> 
#   summarise(tot = sum(NomActive)) |> 
#   mutate(type = "old")
# 
# effort_new<-effort |> 
#   filter(LME == 1) |> 
#   group_by(LME, Year) |> 
#   summarise(tot = sum(NomActive)) |> 
#   mutate(type = "new")
# 
# all<-effort_old |> 
#   full_join(effort_new)
# 
# ggplot(all,aes(x = Year, y = tot, group = type, color = type))+
#   geom_line()+
#   facet_wrap(~type)

# calculate climate inputs by Year as effort is by Year 
# no - skip as Julia would like monthly inputs

# WARNING - should this be annual mean or sum? - CHECKED with JULIA, do mean()
# CAMI THESE ARE RATES (value/m2) and what you are after is the mean annual RATE 
# mean(value/m2) - this is not annual absolute value sum(value) as you keep 
# thinking! STOP ASKING. 


DBPM_FAO_climate_inputs <- combined_FAO_inputs

# add LME total area and calculate effort/m2  
FAO_area<-combined_FAO_inputs |> 
  # this is FAO_area total area
  dplyr::select(FAO, area_m2) |> 
  # all depths are almost definitely > 200m in FAO regions excluding LMEs
  mutate(deptho_m = 200) |> 
  distinct()

# calculate sum of effort by LME/by total area of LME 
DBPM_FAO_effort_input <- effort_FAO |> 
  group_by(Year, fao_area) |> 
  summarize(NomActive = sum(NomActive), .groups = "drop") |> 
  ungroup() |> 
  full_join(FAO_area, by = c("fao_area" = "FAO")) |> 
  mutate(NomActive_area_m2 = NomActive/(area_km2*1e6))

# do the same with catches 
# catch<-read_csv("/rd/gem/private/users/yannickr/DKRZ_EffortFiles/
#calibration_catch_histsoc_1850_2004.csv")
catch <- file.path("/rd/gem/private/users/yannickr/DKRZ_EffortFiles", 
                   "catch-validation_isimip3a_histsoc_1850_2004.csv") |> 
  read_csv()
catch_FAO <- catch |>
  dplyr::filter(LME == 0)

DBPM_FAO_catch_input<-catch_FAO |> 
  mutate(catch_tonnes = Reported+IUU) |> 
  group_by(Year, fao_area) |> 
  # catch is in tonnes, checked in FishingEffort Rproject, 
  summarize(catch_tonnes = sum(catch_tonnes), .groups = "drop") |> 
  # also Reg advise to exclude discards 
  ungroup() |> 
  full_join(FAO_area, by = c("fao_area" = "FAO")) |> 
  mutate(catch_tonnes_area_m2 = catch_tonnes/(area_km2*1e6))

DBPM_FAO_effort_catch_input <- DBPM_FAO_effort_input |> 
  full_join(DBPM_FAO_catch_input)

head(DBPM_FAO_effort_catch_input)


#### 6. Plot to check ----

# WARNING - keep going with the checks

my_theme <- theme_bw()+
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9), 
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key.size = unit(0.1, "cm")) 


plot_df <- split(DBPM_FAO_effort_catch_input, 
                 DBPM_FAO_effort_catch_input$fao_area)

plot_df[[1]] 
head(plot_df[[2]])

Value = "NomActive"
toKeep = "3"

plot <- ggplot(data = plot_df[[3]], aes(x = Year, y = NomActive))+
  ggtitle(paste("FAO", toKeep, sep = " "))+
  annotate("rect", xmin = 1841, xmax = 1960, ymin = 0, ymax = Inf, 
           fill = "#b2e2e2", alpha = 0.4)+ # spin-up edf8fb
  annotate("rect", xmin = 1961, xmax = 2010, ymin = 0, ymax = Inf, 
           fill = "#238b45", alpha = 0.4)+ # projection 66c2a4
  geom_point(size = 1)+
  geom_line() +
  my_theme

pdf("Output/Effort_FAO1_check.pdf")
plot
dev.off()

#### 7. go to step 2 - batch create inputs.R ----
# in batch create input, the code uses getgridin_ISIMIP3b to calculate 
# intercept and slope the below are the steps you need from that function 
# (checked with Julia - line 14 to 22)

# NOTE/WARNING - in this repo GetPPIntSlope() is in 
# source("dbpm_model_functions.R"). 
# this code was taken from the dbpm_isimip_3a repo 
# (original source = dbpm_isimip_3b - I think as things changed...)
# source("input_funcs.R")
source("dbpm_model_functions.R")

# name columns as in code
# head(DBPM_LME_climate_inputs)
DBPM_FAO_climate_inputs_renamed <- DBPM_FAO_climate_inputs |>
  left_join(FAO_area, by = join_by("FAO")) |>
  mutate(deptho_m = 200, area_m2 = area_km2*1e6) |> 
  dplyr::select(FAO, Date, deptho_m, area_m2, `expc-bot`, tob, tos, sphy, 
                lphy) |> 
  rename(t = Date, depth = deptho_m, sbt = tob, sst = tos, 
         expcbot = `expc-bot`) |> 
  # reorder columns 
  relocate(any_of(c("FAO", "t", "lphy", "sphy", "sbt", "sst", "depth", 
                    "area_m2", "expcbot")))

# DBPM_FAO_climate_inputs_renamed <-DBPM_FAO_climate_inputs_renamed[,c("FAO", 
# "t", "lphy", "sphy", "sbt", "sst", "depth", "area_m2", "expcbot")]

### WARNING - change the below with summarise() and test the difference - 
# same problem as per gridded inputs?  
head(DBPM_FAO_climate_inputs_renamed)

DBPM_FAO_climate_inputs_slope <- DBPM_FAO_climate_inputs_renamed |>
  group_by(FAO, t, lphy, sphy, sbt, sst, depth, area_m2, expcbot) |> 
  summarise(er = getExportRatio(sphy, lphy, sst, depth),
            er = ifelse(er < 0, 0, ifelse(er > 1, 1, er)),
            intercept = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                      mmid = 10^-10.184,
                                      mmax = 10^-5.25, depth, 
                                      output = "intercept"),
            slope = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                  mmid = 10^-10.184, mmax = 10^-5.25, depth, 
                                  output = "slope"), .groups = "drop") |> 
  ungroup() |> 
  relocate("FAO", "t", "sst", "sbt", "er", "intercept", "slope", "sphy", 
           "lphy", "depth", "area_m2","expcbot")

head(DBPM_FAO_climate_inputs_slope)

#### 8. Save results ---- 
# when depth integration in GetPPIntSlope() is on
# fwrite(x = DBPM_LME_climate_inputs_slope,
#        file.path("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings
#/lme_inputs/obsclim/0.25deg/",
#                  "DBPM_LME_climate_inputs_slope_depthIntegrated.csv"))
fwrite(x = DBPM_FAO_climate_inputs_slope, 
       file.path("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings",
                 "fao_inputs/obsclim/0.25deg/DBPM_FAO_climate_inputs_slope.csv"))
fwrite(x = DBPM_FAO_effort_catch_input, 
       file.path("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings", 
                 "fao_inputs/obsclim/0.25deg/DBPM_FAO_effort_catch_input.csv"))

# print into public folder too
fwrite(x = DBPM_FAO_climate_inputs_slope, 
       file.path("/rd/gem/public/fishmip/ISIMIP3a/InputData/DBPM_fao_inputs", 
                 "obsclim/0.25deg/DBPM_FAO_climate_inputs_slope.csv"))
fwrite(x = DBPM_FAO_effort_catch_input,
       file.path("/rd/gem/public/fishmip/ISIMIP3a/InputData/DBPM_fao_inputs", 
                 "obsclim/0.25deg/DBPM_FAO_effort_catch_input.csv"))

##### END LME scale -----


#### ISIMIP3a scale 2 -----
# use FAO inputs provided by Denisse in September 2023

rm(list=ls())

#### 1. Define main function: for each FAO, and each variable, apply the ----
# calc_inputs_FAO function -----

# 1A. calc_inputs_FAO()
# see function description above - this is adjusted for gridcell level but names
# are kept - hence e.g. weighted_mean_obs
calc_inputs_FAO <- function(file_name_obs, file_name_crtl, file_path_crtl, 
                            file_path_obs, degree){
  
  # extract variable name
  if(degree == "1deg"){
  variable <- str_match(file_name_obs, 
                        "gfdl-mom6-cobalt2_obsclim_\\s*(.*?)\\s*_60arcmin")[2]
  
  #CHECK: May need to replace this with new masks
  # dataframe of areas for 1 degree, global
  area_frame <- data.frame("lon" = rep(seq(-179.5, 179.5, by = 1), by = 180), 
                           "lat" = rep(seq(-89.5, 89.5, by = 1), each = 360), 
                           "area_m2" = rast(nrows = 180, ncol = 360) |> 
                             cellSize() |> 
                             as.data.frame() |> 
                             pull(area))
                           # "area_m2" = 1e6*as.vector(t(as.matrix(area(raster(nrows = 180, 
                           #                                                   ncol = 360))))))
  }
  
  if(degree == "0.25deg"){
  variable <- str_match(file_name_obs, 
                        "gfdl-mom6-cobalt2_obsclim_\\s*(.*?)\\s*_15arcmin")[2]
  
  # dataframe of areas for 0.25 degree, global
  # area_frame <- data.frame("lon" = rep(seq(-179.875,179.875, by = 0.25), by = 720), 
  #                          "lat" = rep(seq(-89.875,89.875, by = 0.25), each = 1440), 
  #                          "area_m2" = 1e6*as.vector(t(as.matrix(area(raster(nrows = 720, 
  #                                                                            ncol = 1440))))))
  #CHECK: May need to replace this with new masks
  #Create spat raster object - 0.25 deg resolution
  area_frame <- rast(resolution = 0.25, nrows = 720, ncol = 1440) |> 
    #Calculate area of grid cell 
    cellSize() |> 
    #Transform to data frame
    as.data.frame(xy = T) |>
    #Rename columns
    rename("lon" = "x", "lat" = "y", "area_m2" = "area")
  
  }
 
  ##############################################################################
  #### RESHAPE DATA (so FAO data conforms to LME data shape, RFH Sep. 2023) ####
  ##############################################################################
  
  # work on CONTROLCLIM first, OPEN AND RESHAPE DATA 

  
  fao_crtl_raw <- fread(file.path(file_path_crtl, file_name_crtl), 
                        header = FALSE)
  
  dat1 <- fao_crtl_raw[-3, ]
  dat_names <- unlist(as.vector(dat1[ ,1]))
  all_dat <- t(as.matrix(dat1[ ,-1]))
  
  na_rows <- which(!complete.cases(all_dat))
  
  dates <- strsplit(dat_names[-c(1:2)], "-")
  years <- unlist(lapply(dates, "[[", 1))
  months <- lapply(dates, function(x){month.abb[as.numeric(x[2])]})
  col_names <- c("lat", "lon", paste(months, years, sep = "_"))
  
  fao_crtl <- as.data.frame(all_dat)
  colnames(fao_crtl) <- col_names
  fao_crtl <- fao_crtl |> 
    left_join(area_frame, by = join_by("lat", "lon")) |> 
    # A small handful of cells on edge of domain have nas in most fao regions.
    # Consistent across all variables
    drop_na()
  
  # fao_crtl <- left_join(fao_crtl, area_frame, 
  #                       by = c("lat" = "lat", "lon" = "lon"))
  # A small handful of cells on edge of domain have nas in most fao regions.
  # Consistent across all variables
  # fao_crtl <- fao_crtl |> na.omit()
  
  # then work on OBSERVED, OPEN AND RESHAPE DATA
  fao_raw <- fread(file.path(file_path_obs, file_name_obs), header = FALSE)
  
  dat1 <- fao_raw[-3,]
  dat_names <- unlist(as.vector(dat1[,1]))
  all_dat <- t(as.matrix(dat1[,-1]))
  
  dates <- strsplit(dat_names[-c(1:2)], "-")
  years <- unlist(lapply(dates, '[[', 1))
  months <- lapply(dates, function(x){month.abb[as.numeric(x[2])]})
  col_names <- c("lat", "lon", paste(months, years, sep = "_"))
  
  fao <- as.data.frame(all_dat)
  colnames(fao) <- col_names
  
  fao <- fao |> 
    left_join(area_frame, by = join_by("lat", "lon")) |> 
    # A small handful of cells on edge of domain have nas in most fao regions.
    # Consistent across all variables
    drop_na()
  
  # fao <- left_join(fao, area_frame, by = c("lat" = "lat", "lon" = "lon"))
  # fao <- fao |> na.omit()
  
  ## Clear what we don't need
  rm("fao_crtl_raw", "fao_raw", "dat1", "dates", "all_dat", "years", "months",
     "col_names")
  
  ##############################################################################
  ##############################################################################
  
  # NOTE - through the code "weighed_mean" is only because of code recycling 
  # from fao-aggregated level approach above 
  # here we are considering grid cells so there is no averadging and weighting  
  
  if (variable == "deptho_m"){
    
    # calculate fixed variables - mean depth and area of fao 
    weighted_mean_obs <- fao |> 
      mutate(FAO = str_extract(file_name_obs, "(?<=FAO-LME-).+(?=_monthly)"))
    
    #CHECK: Is this right? Both variables are the same?
    weighted_mean_obs_final <- weighted_mean_obs
    weighted_mean_crtl_final <- weighted_mean_obs 
    
  }else{
    
    # CONTROL 
    weighted_mean_crtl <- fao_crtl |> 
      # gather(key = "Date",
      #        value = "value",
      #        -c("lat", "lon","area_m2")) |>
      pivot_longer(-c("lat", "lon", "area_m2"), names_to = "Date", 
                   values_to = "value") |> 
      mutate(FAO = str_extract(file_name_crtl, "(?<=FAO-LME-).+(?=_monthly)"),
             Month = str_extract(Date, "[[:upper:]]+[[:lower:]]+"),
             Year = str_extract(Date, "\\d+"), 
             Date = lubridate::my(paste(Month, Year, sep = "."))) |> 
      arrange(Date) |> 
      dplyr::select(-area_m2) |> 
      # reorder columns 
      relocate(any_of(c("lat", "lon", "FAO", "Date", "Year", "Month", "value")))
    
    # reorder columns 
    # weighted_mean_crtl <- weighted_mean_crtl[, c("lat", "lon", "FAO", "Date",
    #                                              "Year", "Month", "value")]
    
    # rename weighed_mean column according to variable 
    names(weighted_mean_crtl)[ncol(weighted_mean_crtl)] <- variable
    
    # OBSERVED 
    weighted_mean_obs <- fao |> 
      # gather(key = "Date",
      #        value = "value",
      #        -c("lat", "lon","area_m2")) |> 
      pivot_longer(-c("lat", "lon", "area_m2"), names_to = "Date", 
                   values_to = "value") |> 
      mutate(FAO = str_extract(file_name_crtl, "(?<=FAO-LME-).+(?=_monthly)"),
             Month = str_extract(Date,"[[:upper:]]+[[:lower:]]+"),
             Year = str_extract(Date, "\\d+"), 
             Date = lubridate::my(paste(Month,Year, sep = "." ))) |> 
      arrange(Date) |> 
      dplyr::select(-area_m2) |> 
      # reorder columns 
      relocate(any_of(c("lat", "lon", "FAO", "Date", "Year", "Month", "value")))
    
    # reorder columns 
    # weighted_mean_obs <- weighted_mean_obs[, c("lat", "lon", "FAO", "Date", 
    #                                            "Year", "Month", "value")]
    
    # rename weighed_mean column according to variable 
    names(weighted_mean_obs)[ncol(weighted_mean_obs)] <- variable
    
  }
  
  return(list(weighted_mean_obs_final = weighted_mean_obs, 
              weighted_mean_crtl_final = weighted_mean_crtl)) 
  
}

## grid cell level FAO are too big to save with spin up, need to process first 
# then add spinup


# 1.B calc_inputs_all_FAO()
# see function description above - this is adjusted for gridcell level but names
# are kept - hence e.g. weighted_mean_obs
# see LME level 1 for explanation on how this function is applied. 

# NOTE/WARNING - in this repo GetPPIntSlope() is in 
# source("dbpm_model_functions.R"). 
# this code was taken from the dbpm_isimip_3a repo 
# (original source = dbpm_isimip_3b - I think as things changed...)
# source("input_funcs.R")
source("dbpm_model_functions.R")

calc_inputs_all_FAO<-function(this_FAO, degree){
  
  
  ## WARNING - you need to add the 0.25deg option later on - i.e. ifelse()
  file_path_obs<- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/fao_inputs/", 
                        "obsclim/", degree)
  file_path_crtl<- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/fao_inputs/",
                          "ctrlclim/", degree)
  
  if(degree == "1deg"){
    # dataframe of areas for 1degree, global
    area_frame <- data.frame("lon" = rep(seq(-179.5, 179.5, by = 1), by = 180), 
                             "lat" = rep(seq(-89.5, 89.5, by = 1), each = 360), 
                             "area_m2" = rast(nrows = 180, ncol = 360) |> 
                               #Calculate area of grid cell 
                               cellSize() |> 
                               #Transform to data frame
                               pull(area)
                             # "area_m2" = 1e6*as.vector(t(as.matrix(area(raster(nrows = 180,
                             #                                                   ncol = 360))))))
  }
  
  if(degree == "0.25deg"){
    # dataframe of areas for 0.25degree, global
    # area_frame <- data.frame("lon" = rep(seq(-179.875,179.875, by = 0.25), 
    #                                      by = 720),
    #                          "lat" = rep(seq(-89.875,89.875, by = 0.25), 
    #                                      each = 1440),
                             # "area_m2" = 1e6*as.vector(t(as.matrix(area(raster(nrows = 720, 
                             #                                                   ncol = 1440))))))
    #CHECK: May need to replace this with new masks
    #Create spat raster object - 0.25 deg resolution
    area_frame <- rast(resolution = 0.25, nrows = 720, ncol = 1440) |> 
      #Calculate area of grid cell 
      cellSize() |> 
      #Transform to data frame
      as.data.frame(xy = T) |>
      #Rename columns
      rename("lon" = "x", "lat" = "y", "area_m2" = "area")
    
  }
  
  this_FAO_new<-paste0("FAO-LME-", this_FAO, sep = "")
  
  fao_obs <- list.files(file_path_obs, pattern = this_FAO_new)
  fao_ctrl <- list.files(file_path_crtl, pattern = this_FAO_new)  
  
#  tic()
  output_obs<-list()
  output_crtl<-list()
  
  # loop across inputs: depth, phyc tob, expt, phypico, tos. 
  for(i in 1:length(fao_obs)){
    print(i)
    a<-calc_inputs_FAO(file_name_obs = fao_obs[[i]], 
                       file_name_crtl = fao_ctrl[[i]], 
                       file_path_crtl = file_path_crtl, 
                       file_path_obs = file_path_obs,
                       degree)
    output_obs[[i]] <- a$weighted_mean_obs_final
    output_crtl[[i]] <- a$weighted_mean_crtl_final
  }
 # toc() #
  
  rm("a")
  
  # all inputs together for one FAO
  
  # CTRLCLIM 
  # METHOD 1 
  output_crtl_all_variables<-Reduce(cbind,output_crtl)
  rm("output_crtl")
  
  # colnames(output_crtl_all_variables2)
  output_crtl_all_variables <- output_crtl_all_variables[, c(1:7, 14, 21, 28, 
                                                             35)]
  output_crtl_all_variables <- left_join(output_crtl_all_variables, area_frame,
                                       by = c("lat" = "lat", "lon" = "lon"))
  
  # If no depth variable
  if(!("depth" %in% colnames(output_crtl_all_variables))){ 
  output_crtl_all_variables<- output_crtl_all_variables |> 
    arrange(lat, lon, Date) |> 
    mutate(depth = 200)
  }
  #toc() # 4 sec
  
  
  output_crtl_all_variables <- output_crtl_all_variables |> 
    mutate(sphy = `phypico-vint`,
           lphy = `phyc-vint` - `phypico-vint`) |> 
    dplyr::select(-c(`phyc-vint`,`phypico-vint`))
  
  output_crtl_all_variables <- output_crtl_all_variables |> 
    dplyr::select(lat, lon, FAO, Date, Year, Month, depth, area_m2, `expc-bot`,
                  tob, tos, sphy, lphy) |> 
    rename(t = Date, sbt = tob, sst = tos, expcbot = `expc-bot`) |> 
    relocate(any_of(c("lat","lon","FAO", "t", "Year", "Month", "lphy", "sphy", 
                      "sbt", "sst", "depth", "area_m2", "expcbot")))
  
  # output_crtl_all_variables <- output_crtl_all_variables[,c("lat","lon","FAO",
  #                                                           "t", "Year", "Month", 
  #                                                           "lphy", "sphy", "sbt",
  #                                                           "sst", "depth", 
  #                                                           "area_m2", "expcbot")]
  
  # # dplyr method, WORKS FOR ME - RFH SEP 2023
  output_crtl_all_variables <- output_crtl_all_variables |>
    ungroup() |>
    mutate(er = getExportRatio(sphy, lphy, sst, depth),
           er = ifelse(er < 0, 0, ifelse(er > 1, 1, er)),
           er = round(er, 3),
           sst = round(sst, 3),
           sbt = round(sbt, 3),
           intercept = round(GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                           mmid = 10^-10.184, mmax = 10^-5.25, 
                                           depth, output = "intercept"), 3),
           slope = round(GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                       mmid = 10^-10.184, mmax = 10^-5.25, 
                                       depth, output = "slope"), 3),
           sphy = round(sphy, 5),
           lphy = round(lphy, 5)) |> 
    relocate("lat", "lon", "FAO", "t", "sst", "sbt", "er", "intercept", "slope",
             "sphy", "lphy", "depth", "area_m2", "expcbot")

  ### CHECK METHOD
# output_crtl_all_variables[1,]
 # lat     lon FAO          t    sst   sbt    er intercept  slope    sphy    lphy depth   area_m2
#  1 5.125 -47.875  31 1961-01-01 27.677 1.945 0.061    -0.509 -1.123 0.05052 0.01246   200 766261927
 # expcbot Year Month
#  1 3.638025e-09 1961   Jan
  # intercept = -0.509, round(GetPPIntSlope(0.05052,0.01246,mmin=10^-14.25,
  #                                         mmid=10^-10.184,mmax=10^-5.25,200,
  #                                         output="intercept"), 3)
  # slope = -1.123, round(GetPPIntSlope(0.05052,0.01246,mmin=10^-14.25,
  #                                     mmid=10^-10.184,mmax=10^-5.25,200,
  #                                     output="slope"), 3)
  # er = 0.06, round(getExportRatio(0.05052,0.01246,28.407,200), 3)
  
  ## Save output
  this_destination_path_ctrl <- paste0("/rd/gem/private/fishmip_inputs/", 
                                       "ISIMIP3a/processed_forcings/", 
                                       "fao_inputs_gridcell/ctrlclim/", degree,
                                       "/control_historical_FAO_", this_FAO, 
                                       "_", degree, ".csv")
  fwrite(x = output_crtl_all_variables, 
         file = file.path(this_destination_path_ctrl))
  rm("output_crtl_all_variables")
  
  
  ### CALCULATE SPINUP
  print("getting spinup")
  kk <- fread(this_destination_path_ctrl)
  
  spinup<- kk |> arrange(lat, lon, Year) |>
    dplyr::select(-t) |> 
    dplyr::filter(Year >= 1961, Year <=1980) 
  
  # 6 cycles of 20 years - here you are just repeating the spinup above for 6 times,
  # I couldn't figure out how to modify years with dplyr, so using a loop, RFH
  num_times <- 6
  
  for(i in 1:num_times){
    curr_spin <- spinup |> 
      dplyr::mutate(Year = Year - (120-(20*(i-1))))
    
    if(i == 1){
      all_spin <- curr_spin
    }else{
      all_spin <- rbind(all_spin, curr_spin)
    }
    rm("curr_spin")
  }
  
  rm("spinup")
  
  spinup <- all_spin
  
  rm("all_spin")
  
  spinup <- spinup |> 
    mutate(t = lubridate::my(paste(Month, Year, sep = "." ))) |>
    dplyr::select(-c(Month, Year))
  
  ### JUST SAVE SPINUP
  this_destination_path_ctrl <- paste0("/rd/gem/private/fishmip_inputs/",
                                       "ISIMIP3a/processed_forcings/",
                                       "fao_inputs_gridcell/ctrlclim/", degree, 
                                       "/control_spinup_FAO_", this_FAO, "_",
                                       degree, ".csv")
  fwrite(x = spinup, file = file.path(this_destination_path_ctrl))
  rm("spinup")
  
  ### ADD SPINUP TO CTRLCLIM AND SAVE
#  print("saving full ctrlclim")
#  kk <- kk |> dplyr::select(-c(Month, Year)) |> dplyr::mutate(t = as.Date(t))

#  kk <- rbind(spinup, kk)
  
#  this_destination_path_ctrl <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/
  # processed_forcings/fao_inputs_gridcell/ctrlclim/", degree, "/control_FAO_", 
  # this_FAO, "_", degree, ".csv")
#  fwrite(x = kk, file = file.path(this_destination_path_ctrl))
#  rm("kk")
  
  ## Delete ctrlclim file for current FAO that does not include spinup
#  to_delete <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/
  # fao_inputs_gridcell/ctrlclim/", degree, "/control_no_spinup_FAO_", this_FAO,
  # "_", degree, ".csv")
 # unlink(to_delete)
  
  # OBSCLIM
  
  # METHOD 1 - much faster but less reliable 
  # reduce from list to df and column bind all inputs from loop above 
  # except for depth and grid cell area which do not have date 
  output_obs_all_variables<-Reduce(cbind,output_obs) 
  rm("output_obs")
  
  # lat, lon, FAO, Date, Year, Month are repeated so need to remove these column 
  #repetitions 
  output_obs_all_variables <- output_obs_all_variables[,c(1:7, 14, 21, 28, 35)]
  # add grid cell area and depth (set depth to 200m for all grid cells)
  output_obs_all_variables <- output_obs_all_variables |> 
    left_join(area_frame, by = join_by("lat", "lon"))
  
  if(!("depth" %in% colnames(output_obs_all_variables))){ # If no depth variable
  output_obs_all_variables <- output_obs_all_variables |> 
    arrange(lat, lon, Date) |>
    mutate(depth = 200)
  }
  
  output_obs_all_variables<-output_obs_all_variables |> 
    mutate(sphy = `phypico-vint`, lphy = `phyc-vint` - `phypico-vint`) |> 
    dplyr::select(-c(`phyc-vint`,`phypico-vint`))
  
  ### steps 4 and 7 in LME scale 1 are now moved in this function for LME scale 
  # 2 as the final output is LME files and not one aggregated final file. 
  ### steps 3, 5, 6, 8 not needed here
  
  output_obs_all_variables <- output_obs_all_variables |> 
    dplyr::select(lat, lon, FAO, Date, Year, Month, area_m2, depth, `expc-bot`,
                  tob, tos, sphy, lphy) |> 
    rename(t = Date, sbt = tob, sst = tos, expcbot = `expc-bot`) |> 
    relocate(any_of(c("lat", "lon", "FAO", "t", "Year", "Month", "lphy", "sphy",
                      "sbt", "sst", "depth", "area_m2", "expcbot")))
  
  # output_obs_all_variables<-output_obs_all_variables[,c("lat","lon","FAO", "t", 
  #                                                       "Year","Month","lphy", 
  #                                                       "sphy", "sbt", "sst",
  #                                                       "depth", "area_m2",
  #                                                       "expcbot")]
  
  # # dplyr method, CHECKED - WORKS FOR ME, RFH SEP 2023
  output_obs_all_variables <- output_obs_all_variables |>
    ungroup() |>
    mutate(er = getExportRatio(sphy, lphy, sst, depth),
           er = ifelse(er < 0, 0, ifelse(er > 1, 1, er)),
           er = round(er, 3),
           sst = round(sst, 3),
           sbt = round(sbt, 3),
           intercept = round(GetPPIntSlope(sphy, lphy, mmin = 10^-14.25,
                                           mmid = 10^-10.184, mmax = 10^-5.25, 
                                           depth, output = "intercept"), 3),
           slope = round(GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                       mmid = 10^-10.184, mmax = 10^-5.25, 
                                       depth, output = "slope"), 3),
           sphy = round(sphy, 5),
           lphy = round(lphy, 5)) |> 
    relocate("lat", "lon", "FAO", "t", "sst", "sbt", "er", "intercept", 
             "slope", "sphy", "lphy", "depth", "area_m2", "expcbot")
  
  # write output files 
  this_destination_path_obs <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a",
                                      "/processed_forcings/fao_inputs_gridcell",
                                      "/obsclim/", degree, 
                                      "/observed_historical_FAO_", this_FAO, 
                                      "_", degree, ".csv")

  fwrite(x = output_obs_all_variables, 
         file = file.path(this_destination_path_obs))
  rm("output_obs_all_variables")
  
  ##### ADD SPINUP
 # kk <- data.table::fread(this_destination_path_obs)
  
#  kk <- kk |> dplyr::select(-c(Month, Year)) |> dplyr::mutate(t = as.Date(t))
  
#  kk <- rbind(spinup, kk)
  
#  this_destination_path_obs <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/
  # processed_forcings/fao_inputs_gridcell/obsclim/", degree, "/observed_FAO_", 
  # this_FAO, "_", degree, ".csv")
#  print("saving full obsclim")
#  fwrite(x = kk, file = file.path(this_destination_path_obs))
 # rm("kk", "spinup")
  
  ## Delete obsclim file for current FAO that does not include spinup
#  to_delete <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/
  # fao_inputs_gridcell/obsclim/", degree, "/observed_no_spinup_FAO_", this_FAO, 
  # "_", degree, ".csv")
#  unlink(to_delete)
  
}


#### 2. apply the functions above to each FAO -----
this_FAO = c(77, 31, 41, 87, 57, 58, 71, 81, 21, 51, 34, 27, 47, 48, 61, 67, 88) 

tic()
for (i in 1:length(this_FAO)){
  message("Processing #", i, " of ", length(this_FAO))
  message("0.25deg")
  calc_inputs_all_FAO(this_FAO[[i]], "0.25deg")
  message("1deg")
  calc_inputs_all_FAO(this_FAO[[i]], "1deg")
  
}
toc() 






### check outputs

### using code in dbpm_isimip_3a and function to clacualte inputs stored there 
# when rinning function above - OK
# LME 14 - output_obs_all_variables_slope (output_obs_all_variables_slope is the same for spin-ip)

# lat   lon LME   t            sst   sbt    er intercept  slope   sphy    lphy depth     area_m2       expcbot
# <dbl> <dbl> <chr> <date>     <dbl> <dbl> <dbl>     <dbl>  <dbl>  <dbl>   <dbl> <dbl>       <dbl>         <dbl>
# 1 -54.5 -64.5 14    1841-01-01  7.95  7.88 0.239    -0.772 -1.00  0.145  0.143    111. 7179721107. 0.000000202  
# 2 -54.5 -64.5 14    1841-02-01  8.43  8.41 0.245    -0.789 -0.991 0.108  0.120    111. 7179721107. 0.000000154  
# 3 -54.5 -64.5 14    1841-03-01  8.39  8.38 0.247    -0.862 -0.989 0.0878 0.0994   111. 7179721107. 0.000000106  
# 4 -54.5 -64.5 14    1841-04-01  7.84  7.86 0.232    -1.16  -1.01  0.0717 0.0650   111. 7179721107. 0.0000000528 
# 5 -54.5 -64.5 14    1841-05-01  7.16  7.19 0.195    -1.74  -1.05  0.0566 0.0305   111. 7179721107. 0.0000000181 
# 6 -54.5 -64.5 14    1841-06-01  6.25  6.27 0.159    -2.42  -1.11  0.0414 0.0124   111. 7179721107. 0.00000000733

# # when loading the printed file - OK
# check<-read.csv("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/
# lme_inputs_gridcell/obsclim/1deg/observed_LME_14.csv")
# head(check)
# lat   lon LME          t      sst      sbt        er  intercept      slope
# 1 -54.5 -64.5  14 1841-01-01 7.948583 7.880984 0.2385617 -0.7724442 -1.0013043
# 2 -54.5 -64.5  14 1841-02-01 8.433881 8.406563 0.2447559 -0.7890344 -0.9906357
# 3 -54.5 -64.5  14 1841-03-01 8.389428 8.384884 0.2465563 -0.8622664 -0.9890503
# 4 -54.5 -64.5  14 1841-04-01 7.841928 7.859843 0.2324489 -1.1557924 -1.0086839
# 5 -54.5 -64.5  14 1841-05-01 7.162827 7.194753 0.1952542 -1.7380429 -1.0544107
# 6 -54.5 -64.5  14 1841-06-01 6.245687 6.273799 0.1594469 -2.4159994 -1.1060848
# sphy       lphy    depth    area_m2      expcbot
# 1 0.14509845 0.14296427 111.0009 7179721107 2.022832e-07
# 2 0.10794740 0.12006476 111.0009 7179721107 1.541405e-07
# 3 0.08777226 0.09939922 111.0009 7179721107 1.063319e-07
# 4 0.07172815 0.06498957 111.0009 7179721107 5.279879e-08
# 5 0.05660887 0.03050872 111.0009 7179721107 1.811899e-08
# 6 0.04137433 0.01239686 111.0009 7179721107 7.334109e-09

# when calculating slope and intercept using the function - OK
# intercept = GetPPIntSlope(sphy,lphy,mmin=10^-14.25,mmid=10^-10.184,mmax=10^-5.25,
# depth,output="intercept")
# slope = GetPPIntSlope(sphy,lphy,mmin=10^-14.25,mmid=10^-10.184,mmax=10^-5.25,
# depth,output="slope"))
# GetPPIntSlope(0.14509845,0.14296427,mmin=10^-14.25,mmid=10^-10.184,mmax=10^-5.25,
# 111.0009,output="intercept") # -0.7724442
# GetPPIntSlope(0.14509845,0.14296427,mmin=10^-14.25,mmid=10^-10.184,mmax=10^-5.25,
# 111.0009,output="slope") # -1.001304

# when using code in this repo (copied from dbpm_isimip_3a) 
# and function to calculate inputs stored HERE - OK values are the same 
# lat   lon LME   t            sst   sbt    er intercept  slope   sphy    lphy depth     area_m2       expcbot
# <dbl> <dbl> <chr> <date>     <dbl> <dbl> <dbl>     <dbl>  <dbl>  <dbl>   <dbl> <dbl>       <dbl>         <dbl>
# 1 -54.5 -64.5 14    1841-01-01  7.95  7.88 0.239    -0.772 -1.00  0.145  0.143    111. 7179721107. 0.000000202  
# 2 -54.5 -64.5 14    1841-02-01  8.43  8.41 0.245    -0.789 -0.991 0.108  0.120    111. 7179721107. 0.000000154  
# 3 -54.5 -64.5 14    1841-03-01  8.39  8.38 0.247    -0.862 -0.989 0.0878 0.0994   111. 7179721107. 0.000000106  
# 4 -54.5 -64.5 14    1841-04-01  7.84  7.86 0.232    -1.16  -1.01  0.0717 0.0650   111. 7179721107. 0.0000000528 
# 5 -54.5 -64.5 14    1841-05-01  7.16  7.19 0.195    -1.74  -1.05  0.0566 0.0305   111. 7179721107. 0.0000000181 
# 6 -54.5 -64.5 14    1841-06-01  6.25  6.27 0.159    -2.42  -1.11  0.0414 0.0124   111. 7179721107. 0.00000000733

# check 1960 instead for spin-up time

# printed file using dbpm_isimip_3a
# head(check |> filter(t >= "1960-01-01", t < "1970-01-01"))
# lat   lon LME          t      sst      sbt        er  intercept      slope
# 1 -54.5 -64.5  14 1960-01-01 8.579283 8.525311 0.2226460 -0.8927402 -1.0135689
# 2 -54.5 -64.5  14 1960-02-01 8.909275 8.869553 0.2415954 -0.7875424 -0.9901062
# 3 -54.5 -64.5  14 1960-03-01 8.963905 8.961565 0.2330122 -0.9561290 -0.9991341
# 4 -54.5 -64.5  14 1960-04-01 8.595880 8.623443 0.2113011 -1.3052466 -1.0259838
# 5 -54.5 -64.5  14 1960-05-01 7.855768 7.881054 0.1818500 -1.8346303 -1.0653755
# 6 -54.5 -64.5  14 1960-06-01 7.518065 7.550831 0.1516842 -2.4610482 -1.1082789
# sphy       lphy    depth    area_m2      expcbot
# 1 0.14789762 0.12676890 111.0009 7179721107 1.727242e-07
# 2 0.10694323 0.11966546 111.0009 7179721107 1.590143e-07
# 3 0.09020357 0.09109531 111.0009 7179721107 1.015521e-07
# 4 0.07720137 0.05746743 111.0009 7179721107 5.236447e-08
# 5 0.05905584 0.02809974 111.0009 7179721107 2.056590e-08
# 6 0.03932662 0.01149323 111.0009 7179721107 7.078745e-09

# values caclaulted using the functions above - OK SAME AGAIN... 
# output_obs_all_variables_slope |> filter(t >= "1960-01-01", t < "1970-01-01")
# lat   lon LME   t            sst   sbt    er interc…¹  slope   sphy    lphy
# <dbl> <dbl> <chr> <date>     <dbl> <dbl> <dbl>    <dbl>  <dbl>  <dbl>   <dbl>
# 1 -54.5 -64.5 14    1960-01-01  8.58  8.53 0.223   -0.893 -1.01  0.148  0.127  
# 2 -54.5 -64.5 14    1960-02-01  8.91  8.87 0.242   -0.788 -0.990 0.107  0.120  
# 3 -54.5 -64.5 14    1960-03-01  8.96  8.96 0.233   -0.956 -0.999 0.0902 0.0911 
# 4 -54.5 -64.5 14    1960-04-01  8.60  8.62 0.211   -1.31  -1.03  0.0772 0.0575 
# 5 -54.5 -64.5 14    1960-05-01  7.86  7.88 0.182   -1.83  -1.07  0.0591 0.0281 
# 6 -54.5 -64.5 14    1960-06-01  7.52  7.55 0.152   -2.46  -1.11  0.0393 0.0115 


##### END LME grid cell scale -----



#### ISIMIP3a scale 3 -----
# get gridded GCM inputs for ISIMIP3b phase 1 protocol - from GFDL-ESM4, 
# IPSL-CM6A-LR
# monthly time steps, so need to set up time-varying plankton input into code 
#(as in Q_F, start with climatology, then apply dynamical forcing) 
# gridded values 1 by 1 degree lat/lon
# use parallel to do runs for a bunch of grids cell at the same time

# ------------------------------------------------------ 
# What inputs are needed from GCMS?
# Depends on method used to get the plankton size spectrum:

#  use Woodworth-Jefcoats 2013 GCB paper method:
# get small and large phytoplankton densities (if diazotroph density provided, 
#add to large phyto),
# to get slope and intercept, also get median size of consumer and minimum size 
#of phytoplankton everything else same as above

# --------------------------------------------------------

### CN prepare input files for ISIMIP3a ----
# download controlclim files from DKRZ on 02/09/2022
# download all file in this directory to ctrlclim 
# scp -r b381217@levante.dkrz.de:/work/bb0820/ISIMIP/ISIMIP3a/SecondaryInputData/
#climate/ocean/ctrlclim/global/monthly/historical/GFDL-MOM6-COBALT2/* ./
# move one degree files to 1deg folder (from inside 0.25deg folder): mv *onedeg*'
#../1deg
# obsclim file dowloaded earlier (no record kept) - obsclim are also in 
#SecondaryInputData folder (but from CESM2 not GFDL)
# move all files into a subfolder: mv * subfolder (in case you need) 

# Depth files download on 03/09/2022 
# ctrlclim: /work/bb0820/ISIMIP/ISIMIP3a/SecondaryInputData/climate/ocean/ctrlclim/
#global/fixed/historical/GFDL-MOM6-COBALT2
# files: gfdl-mom6-cobalt2_ctrlclim_deptho_15arcmin_global_fixed.nc; 
#gfdl-mom6-cobalt2_ctrlclim_deptho_onedeg_global_fixed.nc
# obsclim: /work/bb0820/ISIMIP/ISIMIP3a/InputData/climate/ocean/obsclim/global/
#fixed/historical/GFDL-MOM6-COBALT2
# files: gfdl-mom6-cobalt2_obsclim_deptho_15arcmin_global_fixed.nc; 
#gfdl-mom6-cobalt2_obsclim_deptho_onedeg_global_fixed.nc

# RECAP on files 
# obsclim 0.25deg = 39 files 
# obsclim 1deg = 39 files (78 in all)
# obsclime DKRZ = 76 files 

# crtlclim all folder = 74 files (2 of which are obsclim)
# crtlclim DKRZ = 74 files (2 of which are obsclim)

#### code starting ----

rm(list=ls())

library(RNetCDF) 
library(reshape2) 
library(abind) 
library(tictoc)

#### CN ISMIP63a adapted function -----

getGCM<-function(gcmPath = './inputs/', protocol, gcm = 'IPSL-CM6A-LR', 
                 savepath, getdepth = T, vers = 2){
  
  # # CN trial
  # gcmPath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/"
  # protocol = "0.25deg"
  # gcm = "obsclim"
  # savepath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/"
  # getdepth = T
  # vers = 3 # see meaning below
  
  # getGCM(gcmPath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/", 
  #        protocol = "0.25deg", 
  #        gcm = "obsclim", 
  #        savepath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/", 
  #        getdepth = T, 
  #        vers = 3)
  
  phypico_file = list.files(path = paste(gcmPath, gcm, '/', protocol, '/', 
                                         sep = ''), pattern = '*phypico-vint*', 
                            full.names = TRUE)
  phyc_file = list.files(path = paste(gcmPath, gcm, '/', protocol, '/', 
                                      sep = ''), pattern = '*phyc-vint*', 
                         full.names = TRUE)
  to_zb_file = list.files(path = paste(gcmPath, gcm, '/', protocol, '/', 
                                       sep = ''), pattern = '*_tob_*',
                          full.names = TRUE)
  to_zs_file = list.files(path = paste(gcmPath, gcm, '/', protocol, '/', 
                                       sep = ''), pattern = '*_tos_*', 
                          full.names = TRUE)
  
  # get large phy
  # checked with Julia 1/09/2022
  # sphy = phypico-vint_mol_m-2
  # lphy = phyc-vint_mol_m-2 - phypico-vint_mol_m-2 
  lphy <- var.get.nc(open.nc(phyc_file), 
                     'phyc-vint') - var.get.nc(open.nc(phypico_file), 
                                               'phypico-vint')
  
  t <- var.get.nc(open.nc(phypico_file), 'time')
  lon <- var.get.nc(open.nc(phypico_file), 'lon')
  lat <- var.get.nc(open.nc(phypico_file), 'lat')
  
  # Format lphy
  dimnames(lphy) <- list(lon=lon,lat=lat,t=t)
  pp <- melt(lphy)
  names(pp) <-  c("lon","lat","t","lphy")
  pp <- pp[!is.na(pp[,"lphy"]),]
  rm(list = ('lphy'))
  
  # sphy
  sphy <- var.get.nc(open.nc(phypico_file), 'phypico-vint')
  sphy <- as.vector(sphy)
  sphy <- sphy[!is.na(sphy)]
  pp$sphy <- sphy
  rm(list = ('sphy'))
  
  # bottom temperature
  to_zb <- var.get.nc(open.nc(to_zb_file), 'tob')
  
  to_zb <- as.vector(to_zb)
  to_zb <- to_zb[!is.na(to_zb)]
  pp$sbt <- to_zb
  rm(list = ('to_zb'))
  
  # Surface temperature
  to_zs <- var.get.nc(open.nc(to_zs_file), 'tos')
  
  to_zs <- as.vector(to_zs)
  to_zs <- to_zs[!is.na(to_zs)]
  pp$sst <- to_zs
  rm(list = ('to_zs'))
  
  # Standardise colnames
  names(pp) <- c("lon", "lat", "t", "lphy", "sphy", "sbt", "sst")
  
  if(getdepth == T){
    # get depth
    # NOTE: one depth file for gcm(e.g. obsclim)/protocol(e.g. 0.25deg) 
    #combination
    depth_nc <- open.nc(list.files(path = paste(gcmPath, gcm, '/', protocol, 
                                                '/', sep = ''), 
                                   pattern = 'deptho',
                                   full.names = TRUE))
    depth <- var.get.nc(depth_nc, 'deptho') # Depth in metres
    dimnames(depth) <- list(lon=var.get.nc(depth_nc, 'lon'),
                            lat=var.get.nc(depth_nc, 'lat'))
    depth <- melt(depth)
    names(depth) <- c("lon", "lat", "depth")
    depth$gridnum <- 1:length(depth[,1])
    
    # Remove land values (na and 0 depth)
    depth <- depth[!is.na(depth[,"depth"]),]
    depth <- depth[depth[,'depth'] != 0,]
    
    ## Save depth
    depth_save_name <- paste(savepath, gcm, '/', protocol, '/', gcm, "_", 
                             protocol, "_depth.RData", sep = '')
    save(depth, file = depth_save_name, version = vers)
  }
  
  ## Save processed forcings
  print(paste('Now saving forcings for', gcm, protocol, sep = ' '))
  pp_save_name <- paste(savepath, gcm, '/', protocol, '/', gcm, "_", 
                        protocol, ".RData", sep = '')
  save(pp, file = pp_save_name, version = vers) # version = the workspace 
  #format version to use. 
  # NULL specifies the current default format (3). 
  # Version 1 was the default from R 0.99.0 to R 1.3.1 and 
  # version 2 from R 1.4.0 to 3.5.0. 
  # Version 3 is supported from R 3.5.0.
  # "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/obsclim/
  #0.25deg/obsclim_0.25deg.RData"
  
  #remove any objects no longer needed 
  if(getdepth == T){
    rm(pp, depth)
  }else{rm(pp)}
}

#### apply getGCM() to all combo of protocols ----

tic()
getGCM(gcmPath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/", 
       protocol = "0.25deg", gcm = "obsclim", 
       savepath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/", 
       getdepth = T, vers = 3)
toc() # 8.533967 min 
getGCM(gcmPath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/", 
       protocol = "0.25deg", gcm = "ctrlclim", 
       savepath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/", 
       getdepth = T, vers = 3)
getGCM(gcmPath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/", protocol = "1deg", 
       gcm = "obsclim", 
       savepath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/", 
       getdepth = T, vers = 3)
getGCM(gcmPath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/", protocol = "1deg",
       gcm = "ctrlclim",
       savepath = "/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/",
       getdepth = T, vers = 3)

#### CN Calculate spin-up for CMIP63a climate forcings ----

### got here with code cleaning ..... but code below works and was used to 
#calculate the spinup 






## go to batch_run_create_inputs_ISIMIP3b and getgridin_ISIMIP3b.R 
## Option also slow - might need to create a spinup as scenarios file 

rm(list=ls())

# # calc spinup for ctrlclim
# load("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/ctrlclim/
#1deg/ctrlclim_1deg.RData")
# ctrlclim<-pp
#
# load("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/ctrlclim/
#1deg/ctrlclim_1deg_depth.RData")
# gridnum_depth<-select(depth, - depth)
# nrow(gridnum_depth) # 41934 lat/lon cell
#
# load("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/obsclim/
#1deg/obsclim_1deg.RData")
# obsclim<-pp

# assume they are the same...
# load("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/obsclim/
#1deg/obsclim_1deg_depth.RData")
# gridnum_depth<-select(depth, - depth)
# nrow(gridnum_depth) # 41934 lat/lon cell

# calculate the spinup only for crtl and run it as it was a scenario... 

calc_input_spinup_gridcell<-function(inputPath, protocol, subset){
  
  # # # # trial
  # inputPath<-"/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/"
  # protocol = "1deg"
  
  input_file_crtlclim <- file.path(paste0(inputPath, "ctrlclim", '/', protocol), 
                                   paste0("ctrlclim", '_', protocol, ".RData"))
  load(input_file_crtlclim)
  # ctrlclim<-lazy_dt(pp)
  ctrlclim<-pp
  # head(ctrlclim)
  rm(input_file_crtlclim, pp)
  
  # input_file_obsclim <- file.path(paste0(inputPath, "obsclim", '/', protocol), 
  #paste0("obsclim", '_', protocol, ".RData"))
  # load(input_file_obsclim)
  # obsclim<-lazy_dt(pp)
  # # obsclim<-pp
  # rm(input_file_obsclim, pp)
  
  gridnum_depth <- file.path(paste0(inputPath, "ctrlclim", '/', protocol), 
                             paste0("ctrlclim", '_', protocol,"_depth",
                                    ".RData"))
  load(gridnum_depth)
  # gridnum_depth<-lazy_dt(depth) |> select(- depth)
  gridnum_depth<-depth |> select(- depth)
  # head(gridnum_depth)
  
  # calculate timesteps, as t is confusing (inputs: monthly_1961_2010.nc)
  # and add gridnumber from depth file
  Date<-seq(as.Date("1961-01-01"), as.Date("2010-12-01"), by="month")
  t<-ctrlclim |> select(t) |> unique()
  time<-data.frame(t = t, Date = Date)
  # time<-lazy_dt(time)
  rm(Date)
  # head(time)
  
  # CONTROL
  tic()
  ctrlclim<-ctrlclim |>
    full_join(time) |>
    full_join(gridnum_depth) |>
    arrange(Date, gridnum)
  toc() # 21 sec # 5 with lazy
  
  # head(ctrlclim)
  
  # calculate spin-up
  tic()
  spinup<-ctrlclim |>
    filter(Date >= "1961-01-01", Date <="1980-12-01") |>
    slice(rep(1:n(), times = 6)) #|>
  # arrange(Date, gridnum)
  toc() # 13 sec # 0 with lazy
  
  # head(spinup)
  
  # calc new date and t for spin-up
  Date = seq(as.Date("1841-01-01"), as.Date("1960-12-01"), by="months")
  # new_t<-seq((720-length(Date)),720-1) # WARNING - need to fix
  new_t<-seq((time$t[1]-length(Date)),time$t[1]-1)
  Date_df<-data.frame(Date = Date, t = new_t)
  # Date_df<-lazy_dt(Date_df)
  rm(Date, new_t)
  # head(Date_df)
  
  tic()
  new_date<-gridnum_depth |>
    full_join(Date_df, by = character()) |>
    arrange(Date, gridnum)
  toc() # 30 sec # 0 with lazy
  # rm(new_date)
  
  # head(new_date)
  
  # replace new date to spin-up file
  # dtplyr
  tic()
  spinup<-spinup |>
    select(-Date, -t)
  toc()
  
  # head(spinup)
  
  # trying to figure out problem
  # new_date<-as.data.frame(new_date)
  Date_new<-new_date$Date
  t_new<-new_date$t
  rm(new_date)
  length(Date_new) # 60384960
  length(t_new) # 60384960
  nrow(spinup) # 60384960
  
  tic()
  spinup<-spinup |>
    mutate(Date = Date_new,
           t = t_new)
  toc()# 0.01
  # rm(Date_new, t_new)
  
  # head(spinup)
  
  # print final input_file_crtlclim file and remove from environment
  crtlclim_SpinUP <- file.path(paste0(inputPath, "spinup", '/', protocol), 
                               paste0("ctrlclim", '_', protocol,"_SpinUp", 
                                      ".RData"))
  
  tic()
  save(spinup, file = crtlclim_SpinUP, version = 3)
  toc() # 14 sec 
  
  # # check 
  # load("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/spinup/
  #1deg/ctrlclim_1deg_SpinUp.RData")
  # head(spinup)
  
  # tic()
  # fwrite(spinup, crtlclim_SpinUP)
  # toc() # 16
  
  rm(ctrlclim, crtlclim_SpinUP, spinup)
  
  # # add spinup to ctrlclim
  # tic()
  # ctrlclim<-ctrlclim |>
  #   full_join(spinup) |>
  #   arrange(Date, gridnum)
  # toc() # 20 sec # 1 with lazy
  # # rm(spinup)
  # 
  # tic()
  # ctrlclim<-as_tibble(ctrlclim)
  # toc() # 593.884 # 10 min
  # 
  # # print final input_file_crtlclim file and remove from environment
  # crtlclim_withSpinUP <- file.path(paste0(inputPath, "ctrlclim", '/',
  #protocol), paste0("ctrlclim", '_', protocol,"_withSpinUp", ".RData"))
  # save(ctrlclim, file = crtlclim_withSpinUP, version = 3)
  # # fwrite(x = ctrlclim, file = crtlclim_withSpinUP)
  # rm(ctrlclim)
  # 
  # # OBSERVED
  # # add spinup to obsclim file
  # tic()
  # obsclim<-obsclim |>
  #   full_join(time) |>
  #   full_join(gridnum_depth) |>
  #   arrange(Date, gridnum)
  # toc() # 21 sec
  # 
  # # add spinup to crtlclim
  # tic()
  # obsclim<-obsclim |>
  #   full_join(spinup)
  # toc() # 20 sec
  # 
  # # arrange file as per original
  # tic()
  # obsclim<-obsclim |>
  #   arrange(Date, gridnum)
  # toc() # 24 sec
  # 
  # # # Plot one gridcell to check - seems OK
  # # trial<-obsclim |>
  # #   filter(gridnum == 1) |>
  # #   mutate(Year = format(as.Date(Date), "%Y"))  |>
  # #   group_by(Year) |>
  # #   summarise(lphy = mean(lphy)) |>
  # #   ungroup() |>
  # #   mutate(Year = as.numeric(Year))
  # #
  # # plot<-ggplot(trial, aes(x = Year, y = lphy))+
  # #   geom_point()+
  # #   geom_line()+
  # #   annotate("rect", xmin=1961, xmax=1980, ymin=-Inf, ymax=Inf, 
  #fill = "#b2e2e2", alpha = 0.4)
  # #
  # # pdf("Output/plot1_lazy.pdf", height = 4, width = 6)
  # # plot
  # # dev.off()
  # 
  # tic()
  # obsclim<-as_tibble(obsclim)
  # toc() # 593.884 # 10 min
  # 
  # # print final input_file_crtlclim file and remove from environment
  # obsclim_withSpinUP <- file.path(paste0(inputPath, "obsclim", '/', protocol),
  #paste0("obsclim", '_', protocol,"_withSpinUp", ".RData"))
  # save(obsclim, file = obsclim_withSpinUP, version = 3)
  
  
  # rm(obsclim, spinup, time) # ANYTHING ELSE?
  
}

# tic()
# calc_input_spinup_gridcell(inputPath = "/rd/gem/private/fishmip_inputs/
#ISIMIP3a/processed_forcings/", protocol = "1deg", subset = NA)
# toc() # 235.699 # 4 min

# tic() # FILE TOO LARGE
# calc_input_spinup_gridcell(inputPath = "/rd/gem/private/fishmip_inputs/
#ISIMIP3a/processed_forcings/", protocol = "0.25deg")
# toc() # ERRROR Error: cannot allocate vector of size 10.7 Gb

# ONLY SPIN-UP for ctrlclim to run once as it's the same for both scenarios 

tic()
calc_input_spinup_gridcell(inputPath = paste0("/rd/gem/private/fishmip_inputs/", 
                                              "ISIMIP3a/processed_forcings/"),
                           protocol = "1deg")
toc() # 1 min but with lazy and in wired format if saves as Rdata - could try 
#as csv. # not lazy 1.5 min

tic()
calc_input_spinup_gridcell(inputPath = paste0("/rd/gem/private/fishmip_inputs/", 
                                              "ISIMIP3a/processed_forcings/"),
                           protocol = "0.25deg")
toc() # 27 min

# # OPTION 2 at gridcell level ...
# 
# # get gridnumber from depth file 
# load("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/ctrlclim/1deg/
# ctrlclim_1deg_depth.RData")
# gridnum_depth<-select(depth, - depth)
# 
# # # how many grid cell?
# # ncell<-ctrlclim |> filter(t == 720) # consider one time-step
# # nrow(ncell) # 670589 0.25deg deg # 41934 1deg - 16 times more (4*4)
# # nrow(gridnum_depth)
# 
# # prepare data - get gridcell ID 
# tic()
# ctrlclim<-ctrlclim |> 
#   full_join(gridnum_depth)
# toc() # 19 sec  
# 
# this_cell<-gridnum_depth$gridnum
# 
# # calculate Date and t of spinup 
# tic()
# trial<-ctrlclim |> filter(gridnum == this_cell[1])
# toc()# 1 sec
# 
# # inputs: monthly_1961_2010.nc (this shouod be done once)
# Date<-seq(as.Date("1961-01-01"), as.Date("2010-12-01"), by="month")
# t<-unique(trial$t)
# time<-data.frame(t = t, Date = Date)
# 
# # define function that caclautes spin up for each grid cell 
# calc_spinup_gridcell<-function(this_cell){
#   
#   # # trial  
#   this_cell_new = this_cell
#   # protocol = "1deg"
#   
#   trial<-ctrlclim |> filter(gridnum == this_cell_new)
#   
#   trial<-trial |> 
#     full_join(time) |>
#     mutate(Year = format(as.Date(Date), "%Y"),  
#            Month = format(as.Date(Date), "%m"))
#   
#   trial_spinup<-trial |> 
#     filter(Year >= 1961, Year <=1980) |>
#     slice(rep(1:n(), times = 6)) |> 
#     mutate(Year = as.character(rep(1841:1960, each = 12)),
#            Date = lubridate::my(paste(Month,Year, sep = "." ))) |> 
#     select(-t)
#   
#   # calcualte t 
#   # runs : 50 year 
#   # 50*12 = 600 months 
#   # spinup : 120 years 
#   # 120*12 = 1440 months 
#   # t1 should be 720 - 1440 = -720
#   
#   Date_spinup<-unique(trial_spinup$Date)
#   t_spinup<-seq((t[1]-length(Date_spinup)),t[1]-1)
#   
#   time_spinup<-data.frame(t = t_spinup, Date = Date_spinup)
#   trial_spinup<-trial_spinup |> 
#     full_join(time_spinup)
#   
#   # print spin up for each grid cell. 
#   this_destination_path <- paste0("/rd/gem/private/fishmip_inputs/ISIMIP3a/
# processed_forcings/", 
#                                   "ctrlclim/", 
#                                   protocol, 
#                                   "/TempSpinup/",
#                                   "gridnum_", 
#                                   this_cell, 
#                                   ".csv")
#   
#   fwrite(x = trial_spinup, file = file.path(this_destination_path))
#   
# }
# 
# # try function 
# tic()
# calc_spinup_gridcell(this_cell = this_cell[1], time, protocol = "1deg")
# toc()

# tic() # 9 min = 150 files. 30 min for 500 files (1 chunk below). 46 h for all
#files - NOT POSSIBLE 
# for(i in 1:length(this_cell)){
#   
#   # i = 1
#   message("Processing #", i, " of ", length(this_cell))
#   calc_spinup_gridcell(this_cell[[i]])
#   
# }
# toc()

# # prepare for loop in //
# # 30 min - not even 500 1st chunk ... 
# chunk_size <- 500 # chunk size for processing
# cell_chunk <- split(this_cell, ceiling(seq_along(this_cell)/chunk_size))
# length(cell_chunk) # 84 for deg1 
# 
# # trial 
# cell_chunk<-cell_chunk[1] # it shiould take less than 30 min 
# 
# protocol = "1deg"
# 
# tic() # strt 9:30 - it should take much less than 30 min ... 
# for(i in 1:length(cell_chunk)){
# 
#   # i = 1
#   file_chunk <- cell_chunk[[i]]
#   message("Processing chunk #", i, " of ", length(cell_chunk))
#   mclapply(X = file_chunk, FUN = calc_spinup_gridcell, mc.cores = 40)
# 
# }
# toc()





