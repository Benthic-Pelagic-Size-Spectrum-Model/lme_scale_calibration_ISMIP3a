
# Loading libraries -------------------------------------------------------
library(arrow)
library(lhs)
library(dplyr)
library(jsonlite)
source("new_workflow/useful_functions.R")


# Loading DBPM climate and fishing inputs ---------------------------------
dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                         "monthly_weighted", 
                         "dbpm_clim-fish-inputs_fao-58_1841-2010.parquet") |> 
  read_parquet()


# Defining initial fishing parameters -------------------------------------
set.seed(1234)

#Number of rows to be included in fishing parameters data frame
num_iter <- 1

#Construct a hypercube with random numbers. Columns represent five specific 
#parameters needed to run DBPM
fishing_params <- data.frame(randomLHS(num_iter, 5))
#Renaming columns 
colnames(fishing_params) <- c("fmort_u", "fmort_v", "fminx_u", "fminx_v", 
                              "search_vol")

#Adjust range of mi size params, others go from 0-1
fishing_params <- fishing_params |> 
  mutate(fminx_u = fminx_u*2, 
         fminx_v = fminx_v*2,
         # adjust range of search vol, others go from 0-1
         search_vol = search_vol+0.001)


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


result_set <- sizemodel(params)




# attach(params)





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


