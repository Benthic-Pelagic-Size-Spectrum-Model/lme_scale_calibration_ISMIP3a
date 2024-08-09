# calculate spinup which will then be used for the observed
library(dplyr)
library(lubridate)

# CHECK the Year dimension and that 12*ncol(fao_crtl) is the right approach
# OK step by step + plot of values below 

spinup <- weighted_mean_crtl |>
  select(-Date) |> 
  filter(year >= 1961, year <= 1980) 

# 6 cycles of 20 years - here you are just repeating the spinup above for 6 times
spinup <- spinup |> 
  slice(rep(1:n(), times = 6))

# change the date on the 6 cycles above to cover the whole spinup period
# you only need to change the years as months are the same for each year
# dimensions to consider are years, months and coordinates 
# years: 1841:1960
# months: 12
# coordinates: nrow(lme_crtl) which gives the number of unique lat/lon combinations fixed across year/month
# rep(1841:1960, each = 12*nrow(lme_crtl))
spinup <- spinup |> 
  mutate(year = as.character(rep(1841:1960, each = 12*nrow(fao_crtl))), 
         Date = my(paste(Month, year, sep = "-")))

# add spinup to control and check that's all OK 
weighted_mean_crtl_final <- weighted_mean_crtl |> 
  full_join(spinup) |> 
  arrange(Date)

# add spin up to observed and plot to check 
weighted_mean_obs_final <- weighted_mean_obs |> 
  full_join(spinup) |> 
  arrange(Date) 

spinup<- output_crtl_all_variables |>
  select(-t) |> 
  filter(Year >= 1961, Year <= 1980) 

# 6 cycles of 20 years - here you are just repeating the spinup above for 6 times
spinup <- spinup |> 
  slice(rep(1:n(), times = 6)) 

# change the date on the 6 cycles above to cover the whole spinup period
# you only need to change the years as months are the same for each year
# dimensions to consider are years, months and coordinates 
# years: 1841:1960
# months: 12
# coordinates: nrow(lme_crtl) which gives the number of unique lat/lon combinations fixed across year/month
# rep(1841:1960, each = nrow(spinup)/(6*20))) # 6 reps of 20 years
num_reps <- nrow(spinup)/(6*20)

spinup <- spinup |> 
  mutate(year = as.character(rep(1841:1960, each = num_reps)), 
         Date = lubridate::my(paste(Month, year, sep = "-"))) 

# add spinup to control and check that's all OK 
output_crtl_all_variables <- output_crtl_all_variables |> 
  full_join(spinup) |> 
  arrange(Date)