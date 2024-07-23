# STEP 1: GET GCM INPUTS FOR DYNAMIC BENTHIC-PELAGIC SIZE SPECTRUM MODEL
#### Ryan is running this step for ISIMIP3a, 3 scales: 
# 1 FAO aggregated - for calibration, using 0.25deg inputs only and focusing on
# obsclim
# 2 FAO grid cell level at 1 deg and 0.25 deg - for model running, need both 
# inputs and both scales 

#### ISIMIP3a scale 1 -----
# Use FAO inputs provided by Denisse in September 2023
# NOTE: we are using only 0.25 deg inputs as these are best when calculating 
#means across LMEs. 


# Loading libraries -------------------------------------------------------
library(tidyverse)
library(data.table)
library(parallel)
library(dtplyr)
source("supporting_functions.R")
source("dbpm_model_functions.R")

# Apply calc_inputs_all() function to each FAO region (0.25 degree res) -----
file_path_obs <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs",
                           "obsclim/025deg")
file_path_ctrl <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs", 
                            "ctrlclim/025deg")

region_choice <- c(21, 27, 31, 34, 41, 47, 48, 51, 57, 58, 61, 67, 71, 77, 
                   81, 87, 88)

#Applying function to all chosen regions (gridded outputs)
#Define paths for gridded outputs
out_path_obs <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                          "processed_forcings/fao_inputs_gridcell/obsclim", 
                          "025deg")
out_path_ctrl <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                           "processed_forcings/fao_inputs_gridcell/ctrlclim", 
                           "025deg")

region_choice |> 
  map(~calc_inputs_gridded(file_path_ctrl, file_path_obs, ., out_path_ctrl, 
                           out_path_obs))

#Applying weighting function to all chosen regions
#Defining paths for weighted outputs
out_path_obs <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                          "processed_forcings/fao_inputs/obsclim/025deg")
out_path_ctrl <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                           "processed_forcings/fao_inputs/ctrlclim/025deg")
region_choice |> 
  map(~calc_inputs_all(file_path_ctrl, file_path_obs, ., out_path_ctrl, 
                       out_path_obs))

# Merging processed inputs into a single file -----------------------------
combined_FAO_inputs <- list.files(out_path_obs, full.names = TRUE) |> 
  #Note that the amount of cores available will depend on compute size chosen
  mclapply(FUN = fread, mc.cores = 28) |> 
  rbindlist() |> 
  #Original comment: all depths are almost definitely > 200m in FAO regions 
  mutate(deptho_m = 200,
         #Keep only the ID identifying the FAO region
         region = as.integer(str_remove(region, "FAO-LME-"))) |> 
  as_tibble()

#Get a list of FAO regions and their total area in km2 to calculate effort/m2
FAO_area <- combined_FAO_inputs |> 
  distinct(region, total_area_km2, deptho_m) |> 
  as_tibble()
  
#From original code: Variable names were standardised to use common terminology
# sphy = phypico-vint_mol_m-2
# lphy = phyc-vint_mol_m-2 - phypico-vint_mol_m-2 
combined_FAO_inputs <- combined_FAO_inputs |> 
  mutate(sphy = phypico_vint,  lphy = phyc_vint - phypico_vint) |> 
  #Removing columns not needed
  select(-c(phyc_vint, phypico_vint))

# Loading effort and catches data -----------------------------------------
effort_file_path <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/DKRZ_EffortFiles"

#Effort data
effort_FAO <- file.path(effort_file_path,
                        "effort_isimip3a_histsoc_1841_2010.csv") |> 
  fread() |> 
  filter(LME == 0) |> 
  as_tibble() |> 
  # calculate sum of effort by LME/by total area of LME 
  group_by(Year, fao_area) |> 
  summarize(total_nom_active = sum(NomActive, na.rm = T), 
            .groups = "drop") |> 
  ungroup() |>
  full_join(FAO_area, by = c("fao_area" = "region")) |> 
  mutate(total_nom_active_area_m2 = total_nom_active/(total_area_km2*1e6))

#Catches data
FAO_catch_input <- list.files(effort_file_path,
                                   "catch-validation_isimip3a_histsoc",
                                   full.names = T) |> 
  read_csv() |>
  filter(LME == 0) |> 
  # catch is in tonnes, checked in FishingEffort Rproject
  mutate(catch_tonnes = Reported+IUU) |> 
  group_by(Year, fao_area) |> 
  summarize(catch_tonnes = sum(catch_tonnes), .groups = "drop") |> 
  # also Reg advise to exclude discards 
  ungroup() |> 
  full_join(FAO_area, by = c("fao_area" = "region")) |> 
  mutate(catch_tonnes_area_m2 = catch_tonnes/(total_area_km2*1e6))

#Merging catches and effort data
DBPM_FAO_effort_catch_input <- effort_FAO |> 
  full_join(FAO_catch_input)

#Removing individual data frames
rm(effort_FAO, FAO_catch_input)


# Plotting fish and catch data --------------------------------------------
# Creating plots to ensure data makes sense - The original code was changed
# slightly to match the original saved image

#Split dataset as items in a list based on FAO area
plot_df <- DBPM_FAO_effort_catch_input |> 
  group_by(fao_area) |> 
  group_split() |> 
  #Select first FAO region
  first()

#Plotting data
plot_df |> 
  ggplot(aes(Year, total_nom_active))+
  ggtitle(paste("FAO region #", unique(plot_df$fao_area), sep = " "))+
  # spin-up edf8fb
  annotate("rect", xmin = 1841, xmax = 1960, ymin = 0, ymax = Inf, 
           fill = "#b2e2e2", alpha = 0.4)+ 
  # projection 66c2a4
  annotate("rect", xmin = 1961, xmax = 2010, ymin = 0, ymax = Inf, 
           fill = "#238b45", alpha = 0.4)+ 
  geom_point(size = 1)+
  geom_line()+
  theme_bw()+
  theme(labs(y = "Total nom active"),
        text = element_text(size = 11),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 11),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 9),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key.size = unit(0.1, "cm")) 

#Saving result that matches previous work
ggsave("Output/Effort_FAO1_check_DFA.pdf", device = "pdf")

#Removing variables not in use
rm(plot_df)


# Calculating intercept and slope -----------------------------------------
DBPM_FAO_climate_inputs_slope <- combined_FAO_inputs |>
  mutate(area_m2 = total_area_km2*1e6) |> 
  #Remove unused columns and reorder them
  select(region, date, tos, tob, sphy, lphy, deptho_m, area_m2, 
         expc_bot) |> 
  # name columns as in "dbpm_model_functions.R" script
  rename(FAO = region, t = date, depth = deptho_m, sbt = tob, sst = tos, 
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
  relocate(all_of(c("er", "intercept", "slope")), .before = sphy)


# Saving catch and effort, and inputs data --------------------------------
#Folder where outputs will be stored
folder_out <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                        "processed_forcings/fao_inputs/obsclim/025deg")

#Saving DBPM inputs 
DBPM_FAO_climate_inputs_slope |> 
  fwrite(file.path(folder_out, "DBPM_FAO_climate_inputs_slope.csv"))

#Saving catch and effort data 
DBPM_FAO_effort_catch_input |> 
  fwrite(file.path(folder_out, "DBPM_FAO_effort_catch_input.csv"))



# Apply calc_inputs_all() function to each FAO region (1 degree res) --------
file_path_obs <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs",
                           "obsclim/1deg")
file_path_ctrl <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs", 
                            "ctrlclim/1deg")
out_path_obs <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                          "processed_forcings/fao_inputs_gridcell/obsclim/1deg")
out_path_ctrl <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                           "processed_forcings/fao_inputs_gridcell/ctrlclim", 
                           "1deg")

#Applying function to all chosen regions
region_choice |> 
  map(~calc_inputs_gridded(file_path_ctrl, file_path_obs, ., out_path_ctrl, 
                           out_path_obs))



# ISIMIP3a scale 3 --------------------------------------------------------
#Applying getGCM() function to all experiments and resolutions
base_path <- "/g/data/vf71/fishmip_inputs/ISIMIP3a"

#0.25 degree datasets
getGCM(folder_path = file.path(base_path, "global_inputs/ctrlclim/025deg"),
       save_path = file.path(base_path, "processed_forcings/ctrlclim/025deg")) 

getGCM(folder_path = file.path(base_path, "global_inputs/obsclim/025deg"),
       save_path = file.path(base_path, "processed_forcings/obsclim/025deg"))

#1 degree datasets
getGCM(folder_path = file.path(base_path, "global_inputs/ctrlclim/1deg"),
       save_path = file.path(base_path, "processed_forcings/ctrlclim/1deg"))

getGCM(folder_path = file.path(base_path, "global_inputs/obsclim/1deg"),
       save_path = file.path(base_path, "processed_forcings/obsclim/1deg"))










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



