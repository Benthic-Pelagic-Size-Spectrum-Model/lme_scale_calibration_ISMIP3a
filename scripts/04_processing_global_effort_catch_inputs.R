# Processing catch and effort inputs for Dynamic Benthic-Pelagic Size 
# Spectrum Model (DBPM)
# This step needs to be completed only once at a global scale

# Activate local R library
.libPaths("/g/data/vf71/la6889/R_personal_lib/")

# Loading libraries -------------------------------------------------------
library(arrow)
library(dplyr)
library(janitor)
library(tidyr)
library(stringr)



# library(purrr)
# library(ggplot2)
# library(lubridate)


# Defining base folder
fishing_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a"

# Processing global effort data -------------------------------------------
# Creating summaries of effort per year and region of interest
effort_data_global <- file.path(fishing_folder, "DKRZ_EffortFiles",
                                "effort_isimip3a_histsoc_1841_2010.csv") |> 
  read_csv_arrow(col_select = c("Year", "fao_area", "LME", "NomActive")) |> 
  clean_names() |>
  mutate(region = case_when(lme == 0 ~ fao_area + 100, .default = lme)) |> 
  # calculate sum of effort by area
  group_by(year, region) |> 
  summarise(total_nom_active = sum(nom_active, na.rm = TRUE)) |> 
  ungroup()

# Saving summarised data
effort_data_global |> 
  write_csv_arrow(
    file.path(fishing_folder, "DKRZ_EffortFiles",
              "yearly_effort_fao-lme_isimip3a_histsoc_1841_2010.csv")) 


# Processing global catch data (Watson) -----------------------------------
# This step needs to be completed only once as it is done at a global scale
catch_data_fishmip <- file.path(fishing_folder, "DKRZ_EffortFiles",
                                "catch_histsoc_1869_2017_EEZ_addFAO.csv") |>
  read_csv_arrow(col_select = c("Year", "fao_area", "LME",
                                "Reported", "IUU", "FGroup")) |> 
  clean_names() |> 
  mutate(region = case_when(lme == 0 ~ fao_area + 100, .default = lme))

catch_watson <- catch_data_fishmip |> 
  summarise(tot_reported = sum(reported, na.rm = TRUE),
            tot_iuu = sum(iuu, na.rm = TRUE), .by = c(year, region)) |> 
  # catch is in tonnes. This was checked in "FishingEffort" project
  mutate(catch_tonnes = tot_reported + tot_iuu) |> 
  # also Reg advise to exclude discards
  select(!starts_with("tot_")) 

# Saving summarised data
catch_watson |> 
  write_csv_arrow(
    file.path(fishing_folder, "DKRZ_EffortFiles",
              "yearly_catch_fao-lme_isimip3a_histsoc_1869_2017.csv")) 


# Processing size spectrum - catches --------------------------------------
# Size spectra data from Reg (original file name: SizesinLMET2.csv")
ss_catches <- read_csv_arrow(file.path(
  fishing_folder, "effort_catch_data", "size_spectrum_catches_fao-lme.csv")) |> 
  clean_names()

# Loading names of LMEs and FAO regions
lme_names <- read_csv_arrow(
  "/g/data/vf71/shared_resources/fao_lme_masks/fao-major_lme_keys.csv",
  col_select = c("fao_lme", "corrected_name"))

# Finding maximum and minimum weight classes of catches (per year and AOI)
# The 'Log10MidWt' column was selected based on the 'Plot Size Data 8.R' script 
# from Reg that produces the figure of size spectrum plots for each region
ss_catches_summ <- ss_catches |> 
  left_join(lme_names, by = c("area" = "fao_lme")) |> 
  group_by(year, area, corrected_name) |> 
  summarise(min_fished_weight_class = min(log10mid_wt, na.rm = TRUE),
            max_fished_weight_class = max(log10mid_wt, na.rm = TRUE), 
            .groups = "drop") |> 
  mutate(region = case_when(area < 100 ~ paste0("LME ", area),
                            .default = paste0("FAO ", area)), .after = year) |> 
  rename(region_name = corrected_name) 

# Per-SPECTRUM (U pelagic / V benthic) fished-size window --------------------
# Reg's size_spectrum_catches file is fish-only and taxon-less, so the min/max
# above collapse ALL taxa into ONE window and miss krill/invertebrates. For 
# DBPM's two-spectrum fishery we derive a per-spectrum window from the FGroup 
# catch (catch_histsoc, already used above), mapping each FGroup to a gram size
# range and to U vs V:
#   fish FGroups: W = 0.01 * L^3 (g), L = cm size-class bounds
#     <30cm -> L[4,30]  30-90cm -> L[30,90]  >=90cm -> L[90,200] 
#     < 90cm -> L[4,90]
#   inverts (fixed): krill[1,2] shrimp[3,60] lobsterscrab[100,4000] 
#     cephalopods[20,6000] demersalmollusc[2,500]
#   U = all fish + krill + cephalopods (water-column predators, incl. demersal 
#   fish)
#   V = shrimp, lobsterscrab, demersalmollusc (benthos)
# Each FGroup maps to a gram size RANGE [lo,hi]: fish W=0.01*L^3 with the 
# cm-class bounds, where the smallest classes' LOWER bound is the realistic 
# minimum FISHED length 10 cm (~10 g, e.g. anchovy/sardine gear onset), NOT 
# 4 cm (larvae). Window = [min lo, max hi] over FGroups holding >=0.5% of that 
# spectrum-year-region catch (drops trace bycatch). min = lower edge = smallest
# fished size (region-specific: small-pelagic LMEs ~10 g, toothfish LMEs stay 
# high, krill 1 g).
# NA where a spectrum was not fished that year. (U max later overridden by the 
# real WtMax below.)
uv_summ <- catch_data_fishmip |> 
  mutate(catch = reported + iuu,
         # shrimp 5 g, krill 1 g
         lo = case_when(f_group == "krill" ~ 1, f_group == "shrimp" ~ 5,
                        f_group == "lobsterscrab" ~ 100, 
                        f_group == "cephalopods" ~ 20,
                        # krill 1 g = fine mesh
                        f_group == "demersalmollusc" ~ 10,
                        # 5 cm ~ 1.25 g (graded pelagics)
                        str_detect(f_group, "<30cm")   ~ 0.01 * 5^3,
                        str_detect(f_group, "30-90cm") ~ 0.01 * 30^3,
                        str_detect(f_group, ">=90cm")  ~ 0.01 * 90^3,
                        str_detect(f_group, "<90cm")   ~ 0.01 * 10^3, 
                        .default = NA_real_),
         hi = case_when(f_group == "krill" ~ 2, f_group == "shrimp" ~ 60,
                        f_group == "lobsterscrab" ~ 4000, 
                        f_group == "cephalopods" ~ 6000,
                        f_group == "demersalmollusc" ~ 500,
                        str_detect(f_group, "<30cm")   ~ 0.01 * 30^3,
                        str_detect(f_group, "30-90cm") ~ 0.01 * 90^3,
                        str_detect(f_group, ">=90cm")  ~ 0.01 * 200^3,
                        str_detect(f_group, "<90cm")   ~ 0.01 * 90^3, 
                        .default = NA_real_),
         spectrum = if_else(f_group %in% c("shrimp", "lobsterscrab",
                                           "demersalmollusc"), 
                            "detritivores", "predators")) |>
  filter(catch > 0, !is.na(lo)) |>
  group_by(region, year, spectrum) |>
  mutate(frac = catch / sum(catch)) |>
  filter(frac >= 0.005) |>
  # lower/upper edges
  summarise(min_fished = log10(min(lo)), max_fished = log10(max(hi)), 
            .groups = "drop") |>
  pivot_wider(names_from = spectrum, values_from = c(min_fished, max_fished))

ss_catches_summ <- ss_catches_summ |>
  left_join(uv_summ, by = c("area" = "region", "year")) |>
  # HYBRID max: the pelagic U max uses Reg's real (WtMax-based) 
  # max_fished_weight_class rather than the coarse FGroup cm-class max; the 
  # benthic V max keeps the invert-group max (the fish-dominated 
  # max_fished_weight_class over-extends benthos). min (U, V) stays 
  # FGroup-derived.
  mutate(max_fished_predators = max_fished_weight_class)

# Saving summarised data
ss_catches_summ |> 
  write_csv_arrow(file.path(fishing_folder, "effort_catch_data", 
                            "summary_size_spectrum_catches_fao-lme.csv"))

