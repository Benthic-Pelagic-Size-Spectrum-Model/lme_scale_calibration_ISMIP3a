
# Loading libraries -------------------------------------------------------
library(dplyr)
library(arrow)
library(lhs)
library(jsonlite)
library(lubridate)
source("new_workflow/useful_functions.R")


# Loading DBPM climate and fishing inputs ---------------------------------
dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                         "monthly_weighted", 
                         "dbpm_clim-fish-inputs_fao-58_1841-2010.parquet") |> 
  read_parquet()


# Searching best fishing parameters values for LME ------------------------
test <- LHSsearch(num_iter = 20, 
                  forcing_file = dbpm_inputs, 
                  gridded_forcing = NULL, 
                  best_val_folder = "new_workflow/outputs/best_fish_vals")

test




result <- run_model(fishing_params, dbpm_inputs)
test <- getError(fishing_params, dbpm_inputs)

# Getting DBPM parameters -------------------------------------------------
params <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3, tstepspryr = 12)

# Saving parameters
params |> 
  #Ensuring up to 10 decimal places are saved in file
  write_json("new_workflow/outputs/dbpm_size_params.json", digits = 10)


# Loading DBPM parameters -------------------------------------------------
# If paramaters were already saved, they can be read instead of being 
# recalculated
params <- read_json("new_workflow/outputs/dbpm_size_params.json", 
                    simplifyVector = T)


# result_set <- sizemodel(params)



# attach(params)


no_iter <- 100

test <- LHSsearch(num_iter = no_iter, 
           forcing_file = dbpm_inputs, 
           gridded_forcing = NULL, 
           best_val_folder = "new_workflow/outputs/best_fish_vals")







#Line 357
initial_date <- min(as_date(dbpm_inputs$time))-month(1)
test_pred <- array(0, c(numb_size_bins, numb_time_steps+1)) |> 
  data.frame()
colnames(test_pred) <- c(initial_date, as_date(dbpm_inputs$time))
row.names(test_pred) <- as.character(log10_size_bins)

test_pred[which(log10_size_bins < 0), initial_date]
plank_pred_sizes




test_pred <- expand.grid(time = c(initial_date, as_date(dbpm_inputs$time)), 
                         log10_size = log10_size_bins) |> 
  mutate(vals = NA)

#Line 405
test_pred |> 
  mutate(vals = case_when(time == initial_date & log10_size < 0 ~ )


plank_pred_sizes[1:120]



log10_size_bins







params <- sizeparam(xmin.consumer.u = -3, xmin.consumer.v = -3, 
                    tmax = nrow(input)/12, tstepspryr = 12, 
                    search_vol =  0.39221252,#vals["search.vol"],
                    fmort.u = 0.79832977,#vals["f.u"], 
                    fminx.u = 1.9683917,#vals["f.minu"], 
                    fmort.v = 0.83945284,#vals["f.v"], 
                    fminx.v = 0.82696956,#vals["f.minv"], 
                    depth = mean(input$depth), er = input$export_ratio,
                    pp = input$intercept, slope = input$slope, 
                    sst = input$tos, sft = input$tob, 
                    effort = input$nom_active_area_m2_relative)






gridded_params <- sizeparam(,
                            
                            equilibrium = FALSE, dx = 0.1, 
                            xmin.consumer.u = -3, xmin.consumer.v = -3,
                            tmax = tmax, tstepspryr = 12,
                            search_vol = f.effort["search.vol"], 
                            fmort.u = f.u, fminx.u = f.minu, fmort.v = f.v, 
                            fminx.v = f.minv, depth = lme_inputs_grid$depth, 
                            er = lme_inputs_grid$er, 
                            pp = lme_inputs_grid$intercept, 
                            slope = lme_inputs_grid$slope, 
                            sst = lme_inputs_grid$sst, 
                            sft = lme_inputs_grid$sbt, use.init = TRUE,
                            effort = lme_inputs_grid$nom_active_area_m2_relative, 
                            U.initial = U.initial, V.initial = V.initial, 
                            W.initial = W.initial, 
                            Ngrid = nrow(lme_inputs_grid$depth))


