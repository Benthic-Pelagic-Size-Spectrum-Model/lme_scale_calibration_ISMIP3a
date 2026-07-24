# Creating a single file with all ocean and fishing inputs to run the 
# Dynamic Benthic-Pelagic Size Spectrum Model (DBPM) calibration runs

# Activate local R library
.libPaths("/g/data/vf71/la6889/R_personal_lib/")

# Loading libraries -------------------------------------------------------
library(dplyr)
library(tidyr)
library(arrow)
library(purrr)
library(janitor)
library(ggplot2)
library(stringr)
library(lubridate)

# Define base variables ---------------------------------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_lme_inputs"
fishing_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a"
fao_lme <- list.dirs(base_folder, recursive = FALSE, full.names = FALSE) |> 
  str_subset(pattern = "fao_lme-")


# Loading fishing datasets ------------------------------------------------
global_effort_data <- read_csv_arrow(
  file.path(fishing_folder, "DKRZ_EffortFiles",
            "yearly_effort_fao-lme_isimip3a_histsoc_1841_2010.csv"))

global_catch_watson <- read_csv_arrow(
  file.path(fishing_folder, "DKRZ_EffortFiles",
            "yearly_catch_fao-lme_isimip3a_histsoc_1869_2017.csv"))

global_ss_catches_summ <- read_csv_arrow(
  file.path(fishing_folder, "effort_catch_data", 
            "summary_size_spectrum_catches_fao-lme.csv"))

# Applying workflow to all regions
for(f in fao_lme){
  fao_lme_id <- as.numeric(str_extract(f, "[0-9]+"))
  #Creating path to ocean inputs
  forcing_folder <- file.path(base_folder, f, "monthly_weighted")
  
  # Loading DBPM climate inputs ---------------------------------------------
  # We will double the timesteps for the stable spinup period
  stable_spin <- list.files(
    forcing_folder, pattern = "^stable-spin_dbpm", full.names = TRUE) |> 
    read_parquet() |> 
    replicate(2, expr = _, simplify = FALSE) |> 
    bind_rows() |> 
    mutate(time = seq(as_date("1641-01-01"), as_date("1840-12-31"), 
                      by = "month"), 
           year = year(time), month = month(time, label = TRUE, abbr = FALSE))
    
  # We will load climate inputs to merge with catch and effort data before 
  # saving results
  clim_forcing_file <- list.files(
    forcing_folder, pattern = "obsclim|spinup", full.names = TRUE) |>
    str_subset("inputs_fao_lme") |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows(stable_spin) |>
    arrange(time) |> 
    clean_names() |> 
    mutate(time = as_date(time))
  
  ## Dynamic stable spinup period for the Arctic only -----------------------
  if(fao_lme_id == 64){
    spinup <- read_parquet(list.files(
      forcing_folder, pattern = "spinup", full.names = TRUE))
    
    dyn_spinup <- spinup |> 
      group_by(month) |> 
      summarise(across(where(is.double) & !c(year, time), 
                       ~ mean(.x, na.rm = TRUE))) |> 
      mutate(month = factor(month, levels = month.name, ordered = TRUE)) |> 
      arrange(month) |> 
      replicate(200, expr = _, simplify = FALSE) |> 
      bind_rows() |> 
      mutate(region = str_replace(str_to_upper(f), "-", " "), 
             scenario = "stable-spin", 
             time = seq(as_date("1641-01-01"), as_date("1840-12-31"), 
                        by = "month"), year = year(time), .before = month)
    
    clim_forcing_file <- list.files(forcing_folder, pattern = "obsclim",
                                    full.names = TRUE) |> 
      read_parquet() |> 
      bind_rows(dyn_spinup, spinup) |> 
      arrange(time) |> 
      clean_names() |> 
      mutate(time = as_date(time))
  }
    
  # Getting the mean depth and area of the region of interest
  depth_area <- clim_forcing_file |> 
    distinct(depth, area_m2)

  ## Loading effort data ----------------------------------------------------
  effort_data <- global_effort_data |> 
    #Selecting data for area of interest
    filter(region == fao_lme_id) 
  
  # Extend fishing effort (nom_active_area_m2_relative) starting in 1741
  # Repeat 1841 value for entire stable spin up period
  effort_stable_spin <- effort_data |> 
    filter(year == 1841) |> 
    pull(total_nom_active)
  
  effort_data <- tibble(year = seq(1641, 1840), region = fao_lme_id,
                        total_nom_active = effort_stable_spin) |> 
    bind_rows(effort_data) |> 
    # Adding depth and area information for the area of interest
    mutate(depth = depth_area$depth, 
           area_m2 = depth_area$area_m2,
           total_nom_active_area_m2 = total_nom_active/area_m2,
           nom_active_relative = total_nom_active/max(total_nom_active),
           nom_active_area_m2_relative = total_nom_active_area_m2/
             max(total_nom_active_area_m2))
  
  rm(effort_stable_spin)
  
  # Loading catches data ----------------------------------------------------
  #From Watson et al 2018
  catch_watson <- global_catch_watson |> 
    #Selecting area of interest
    filter(region == fao_lme_id & year <= 2010) |> 
    mutate(depth = depth_area$depth, 
           area_m2 = depth_area$area_m2,
           catch_tonnes_area_m2 = catch_tonnes/area_m2) |> 
    relocate(catch_tonnes, .before = catch_tonnes_area_m2)

  #From Pauly et al 2020
  if(fao_lme_id < 100){
    reg_name <- str_c("LME ", fao_lme_id)
  }else{
    reg_name <- str_c("FAO ", (fao_lme_id-100))
  }
  catch_pauly <- read.csv(
    list.files(file.path(fishing_folder, "SAU_catch_data"), 
               pattern = str_c(reg_name, " v50-1.csv"), full.names = TRUE)) |> 
    # Keep data up to 2010 and removing discards to match processing of Watson
    # data
    filter(year <= 2010 & catch_type != "Discards") |> 
    group_by(year) |> 
    #Calculate total tonnes caught per year
    summarise(catch_tonnes_pauly = sum(tonnes, na.rm = TRUE))
  
  # Load minimum and maximum fish sizes harvested 
  ss_catches_summ <- global_ss_catches_summ |> 
    filter(area == fao_lme_id & year <= 2010) |> 
    select(!c(region, area))
  
  catch_data <- catch_watson |> 
    full_join(catch_pauly, by = "year") |> 
    mutate(catch_pauly_tonnes_area_m2 = catch_tonnes_pauly/area_m2) |> 
    arrange(year) |> 
    filter(!if_all(c(catch_tonnes_area_m2, catch_pauly_tonnes_area_m2), 
                   is.na)) |> 
    rowwise() |>
    mutate(min_catch_density = min(
      catch_tonnes_area_m2, catch_pauly_tonnes_area_m2, na.rm = TRUE),
      max_catch_density = max(
        catch_tonnes_area_m2, catch_pauly_tonnes_area_m2, na.rm = TRUE)) |> 
    select(!c(region, depth, area_m2)) |> 
    full_join(ss_catches_summ, by = "year")
    
  rm(catch_pauly, catch_watson)

  # Merging catch and effort data -------------------------------------------
  DBPM_effort_catch_input <- effort_data |> 
    full_join(catch_data, by = "year") |> 
    mutate(region = reg_name,
           region_name = unique(ss_catches_summ$region_name)) |> 
    relocate(region_name, .after = region) |> 
    filter(year <= 2010)
  
  #Saving summarised catch and effort data
  DBPM_effort_catch_input |> 
    write_parquet(file.path(
      forcing_folder, paste0("dbpm_effort-catch-inputs_", f, ".parquet")))
  
  #Removing individual data frames
  rm(effort_data, catch_data, ss_catches_summ)
  
  #Joining with climate inputs
  forcing_file <- clim_forcing_file |> 
    select(!region) |> 
    full_join(DBPM_effort_catch_input) |> 
    relocate(region, region_name, .after = scenario)


  ## Plotting fish and catch data -------------------------------------------
  forcing_file |> 
    filter(year >= 1841) |> 
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
    labs(y = "Total nom active", 
         title = unique(DBPM_effort_catch_input$region))+
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 10), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) 
  
  folder_out <- file.path("/g/data/vf71/fishmip_outputs/ISIMIP3a",
                          "fao_lme_outputs", f)
  if(!dir.exists(folder_out)){
    dir.create(folder_out, recursive = TRUE)
  }
  
  #Saving results - only save once per FAO region
  fout <- file.path(folder_out, paste0("effort_", f, ".pdf"))
  if(!file.exists(fout)){
    ggsave(fout, device = "pdf", dpi = 300)
  }

  ## Saving catch and effort, and inputs data -------------------------------
  fout_forcing <- paste0(
    "dbpm_clim-fish-inputs_", f, "_", min(forcing_file$year), "-", 
    max(forcing_file$year), ".parquet")
  if(fao_lme_id == 64){
    fout_forcing <- paste0(
      "dbpm_dynamic_clim-fish-inputs_", f, "_", min(forcing_file$year), "-", 
      max(forcing_file$year), ".parquet")
  }
  
  forcing_file |> 
    write_parquet(file.path(forcing_folder, fout_forcing))
}

