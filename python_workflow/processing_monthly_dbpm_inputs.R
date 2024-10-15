# STEP 1: GET GCM INPUTS FOR DYNAMIC BENTHIC-PELAGIC SIZE SPECTRUM MODEL

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(data.table)
library(parallel)
library(dtplyr)
library(arrow)
library(janitor)
source("supporting_functions.R")
source("dbpm_model_functions.R")



# Define functions --------------------------------------------------------
dbpm_input_calc <- function(file_path, path_out){
  # Inputs:
  # - file_path (character) File path where monthly data is located
  # - path_out (character) File path where outputs should be stored as parquet
  # files
  
  #Outputs:
  # - None. This function saves results as parquet files in the path provided.
  
  df <- read_parquet(file_path) |> 
    clean_names() |> 
    rename(sphy = phypico_vint_mol_m_2) |> 
    mutate(lphy = phyc_vint_mol_m_2 - sphy) |> 
    select(!phyc_vint_mol_m_2) |> 
    rename(t = time, depth = depth_m, sbt = tob_deg_c, sst = tos_deg_c, 
           expcbot = expc_bot_mol_m_2_s_1, area_m2 = tot_area_m2) |> 
    mutate(er = getExportRatio(sphy, lphy, sst, depth),
           er = ifelse(er < 0, 0, ifelse(er > 1, 1, er)),
           intercept = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                     mmid = 10^-10.184, mmax = 10^-5.25, depth, 
                                     output = "intercept"),
           slope = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, mmid = 10^-10.184, 
                                 mmax = 10^-5.25, depth, output = "slope")) |> 
    relocate(all_of(c("er", "intercept", "slope")), .before = sphy)
  
  #Save results
  df |> 
    write_parquet(path_out)
}


# LMEs --------------------------------------------------------------------
#Calculating additional inputs from GFDL data
monthly_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs", 
                            "monthly_weighted_mean") |> 
  list.files(full.names = T)

for(f in monthly_inputs){
  f_out  <- str_replace(f, "all-inputs", "climate-inputs-slope")
  dbpm_input_calc(f, f_out)
}


#Getting information about depth and area. We can choose any experiment within
#the region of interest
depth_area <- read_parquet(monthly_inputs[1], 
                           col_select = c("depth_m", "tot_area_m2")) |> 
  distinct()

## Loading effort and catches data ----------------------------------------
effort_file_path <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/DKRZ_EffortFiles"

#Effort data
effort_LME <- file.path(effort_file_path,
                        "effort_isimip3a_histsoc_1841_2010.csv") |> 
  fread() |> 
  #Removing LME = 0 because there is no information about area
  filter(LME == 61) |> 
  # calculate sum of effort by LME/by total area of LME 
  group_by(Year, LME) |> 
  summarize(total_nom_active = sum(NomActive, na.rm = T), 
            .groups = "drop") |> 
  ungroup() |>
  rename(year = Year, region = LME) |> 
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         total_nom_active_area_m2 = total_nom_active/area_m2)

#Catches data
LME_catch_input <- file.path(effort_file_path,
                             "catch-validation_isimip3a_histsoc_1850_2004.csv") |> 
  fread() |>
  #Removing LME = 0 because there is no information about area
  filter(LME == 61) |> 
  # catch is in tonnes, checked in FishingEffort Rproject
  mutate(catch_tonnes = Reported+IUU) |> 
  group_by(Year, LME) |> 
  summarize(catch_tonnes = sum(catch_tonnes), .groups = "drop") |> 
  # also Reg advise to exclude discards 
  ungroup() |> 
  rename(year = Year, region = LME) |> 
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         catch_tonnes_area_m2 = catch_tonnes/area_m2)

#Merging catches and effort data
DBPM_LME_effort_catch_input <- effort_LME |> 
  full_join(LME_catch_input)

#Removing individual data frames
rm(effort_LME, LME_catch_input)


## Plotting fish and catch data -------------------------------------------
DBPM_LME_effort_catch_input |> 
  ggplot(aes(year, total_nom_active))+
  ggtitle("LME # 61")+
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
ggsave("Output/Effort_LME_61_check_DFA.pdf", device = "pdf")

#Removing variables not in use
rm(plot_df)

## Saving catch and effort, and inputs data -------------------------------
#Saving catch and effort data 
DBPM_LME_effort_catch_input |> 
  write_parquet(file.path("/g/data/vf71/la6889/dbpm_inputs/",
                          "DBPM_LME_61_effort_catch_input.parquet"))

