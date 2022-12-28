# JLB - 21/12/2022
# This function reads in LME-scale time-series inputs (climate and fishing) pre-processed 
# as spatially averaged means to estimate catchability (and if needed other model parameters)

library(devtools)
library(tidyverse)

# Function to run model
run_model<-function(input_filepath, output_filepath,model_filepath, LMEnumber=42){
 
  # source model script form Github - NOT WORKING, using local
  # source_url(url = "https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm/blob/master/size-based-models/dynamic_sizebased_model_functions.R")
  source(file = "dbpm_model_functions.R")
 
  # read climate forcing inputs from THREDDS
  lme_clim<-read_csv(file="http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_climate_inputs_slope.csv")
  lme_clim<-subset(lme_clim, LME %in% LMEnumber)
 
  # read climate fishing inputs from THREDDS
  lme_fish<-read_csv(file="http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_effort_catch_input.csv")
  lme_fish<-subset(lme_clim, LME %in% LMEnumber)
  
  #create finer inputs
  
  
  # set up params
  params <- sizeparam(equilibrium = FALSE
                      ,dx = 0.1
                      ,xmin.consumer.u = -3
                      ,xmin.consumer.v = -3
                      ,tmax = length(lme_clim$sst)/12
                      ,tstepspryr  =  12
                      ,fmort.u = 0.0
                      ,fminx.u = 0
                      ,fmort.v = 0.0
                      ,fminx.v = -1
                      ,depth = median(lme_clim$depth)
                      ,er = lme_clim$er
                      ,pp = lme_clim$intercept
                      ,slope = lme_clim$slope
                      ,sst = lme_clim$sst
                      ,sft = lme_clim$sbt
                      ,use.init = FALSE)      
  
  # run model through time
 
  result_set <- sizemodel(params) 
  
  # return outputs of the model to compare with data
  
}


# Error function
getError <-function(model,LMEnumber){
  
  
}


# Run Optim to estimate Catchability Param



