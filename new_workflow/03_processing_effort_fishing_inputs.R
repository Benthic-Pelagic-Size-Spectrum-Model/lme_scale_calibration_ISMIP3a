# Processing fishing and effort inputs for Dynamic Benthic-Pelagic Size 
# Spectrum Model (DBPM)

# Loading libraries -------------------------------------------------------
library(data.table)
library(dplyr)
library(arrow)
library(ggplot2)

# Loading DBPM climate inputs ---------------------------------------------
# We will load 'ctrlclim' data to get the mean depth and area of the region
# of interest
depth_area <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                        "monthly_weighted") |> 
  list.files("ctrlclim_dbpm_clim", full.names = T) |> 
  read_parquet(col_select = c("depth_m", "tot_area_m2")) |> 
  distinct()

## Loading effort and catches data ----------------------------------------
effort_data <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a/DKRZ_EffortFiles",
                         "effort_isimip3a_histsoc_1841_2010.csv") |> 
  #Keep only columns with relevant information
  fread(select = c("Year", "fao_area", "NomActive")) |> 
  #Selecting data for area of interest - FAO 58
  filter(fao_area == 58) |> 
  # calculate sum of effort by area
  group_by(Year, fao_area) |> 
  summarise(total_nom_active = sum(NomActive, na.rm = T), 
            .groups = "drop") |> 
  ungroup() |>
  rename(year = Year, region = fao_area) |> 
  # Adding depth and area information for the area of interest
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         total_nom_active_area_m2 = total_nom_active/area_m2,
         nom_active_relative = total_nom_active/max(total_nom_active),
         nom_active_area_m2_relative = total_nom_active_area_m2/
         max(total_nom_active_area_m2))

#Catches data
catch_data <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a/DKRZ_EffortFiles",
                        "catch-validation_isimip3a_histsoc_1850_2004.csv") |> 
  fread(select = c("Year", "fao_area", "Reported", "IUU")) |>
  #Selecting area of interest
  filter(fao_area == 58) |> 
  # catch is in tonnes. This was checked in "FishingEffort" project
  mutate(catch_tonnes = Reported+IUU) |> 
  group_by(Year, fao_area) |> 
  summarise(catch_tonnes = sum(catch_tonnes), .groups = "drop") |> 
  # also Reg advise to exclude discards 
  ungroup() |> 
  rename(year = Year, region = fao_area) |> 
  # Adding depth and area information for the area of interest
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         catch_tonnes_area_m2 = catch_tonnes/area_m2)

#Merging catches and effort data
DBPM_effort_catch_input <- effort_data |> 
  full_join(catch_data)

#Removing individual data frames
rm(effort_data, catch_data)

## Plotting fish and catch data -------------------------------------------
DBPM_effort_catch_input |> 
  ggplot(aes(year, total_nom_active))+
  annotate("rect", xmin = 1841, xmax = 1960, ymin = 0, ymax = Inf, 
           fill = "#b2e2e2", alpha = 0.4)+ 
  annotate("rect", xmin = 1961, xmax = 2010, ymin = 0, ymax = Inf, 
           fill = "#238b45", alpha = 0.4)+ 
  geom_point(size = 1)+
  geom_line()+
  scale_x_continuous(expand = c(.01, 0), breaks = seq(1850, 2010, 20))+
  scale_y_continuous(expand = c(.02, 0))+
  theme_bw()+
  labs(y = "Total nom active", title = "FAO Major Fishing Area # 58")+
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) 

#Saving result that matches previous work
ggsave("new_workflow/outputs/effort_fao-58.pdf", device = "pdf", dpi = 300)


## Saving catch and effort, and inputs data -------------------------------
DBPM_effort_catch_input |> 
  write_parquet(
    file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica",
              "monthly_weighted/dbpm_effort-catch-inputs_fao-58.parquet"))

