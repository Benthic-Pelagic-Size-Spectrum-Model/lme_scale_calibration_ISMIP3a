# STEP 1: GET GCM INPUTS FOR DYNAMIC BENTHIC-PELAGIC SIZE SPECTRUM MODEL

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



# Calculate spinup from gridded ctrlclim data -----------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/processed_forcings"

calc_input_spinup_gridcell(base_path = file.path(base_folder, "ctrlclim/1deg"), 
                           save_path = file.path(base_folder, "spinup/1deg"))

calc_input_spinup_gridcell(base_path = file.path(base_folder, "ctrlclim/025deg"), 
                           save_path = file.path(base_folder, "spinup/025deg"))


