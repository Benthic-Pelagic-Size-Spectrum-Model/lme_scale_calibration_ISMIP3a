### AUTHOR: RYAN HENEGHAN, SEPTEMBER 2023
### This script takes forcings for obsclim and ctrlclim
### processed by Denisse, and converts them to a form suitable for
### "01_getinputsISIMIP3a.R"

forcings_path <- "/rd/gem/private/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/created_by_denisse_20230524/"
clim <- c("obsclim", "ctrlclim")
resolution <- c("15arcmin", "60arcmin")