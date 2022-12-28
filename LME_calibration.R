# JLB - 21/12/2022
# This function reads in LME-scale time-series inputs (climate and fishing) pre-processed 
# as spatially averaged means to estimate catchability (and if needed other model parameters)

library(devtools)
library(tidyverse)
library(lubridate)

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
  lme_fish<-subset(lme_fish, LME %in% LMEnumber)
  
  #create monthly inputs for fishing
  lme_clim$Year<-year(lme_clim$t)
  lme_clim<-left_join(lme_clim,lme_fish,by="Year")
  # convert fishing effort and catch per yr (divide values by 12)
  lme_clim$NomActive
  lme_clim$NomActive_area_m2
  lme_clim$catch_tonnes
  lme_clim$catch_tonnes_area_m2
  
  # set up params
  params <- sizeparam(equilibrium = FALSE
                      ,dx = 0.1
                      ,xmin.consumer.u = -3
                      ,xmin.consumer.v = -3
                      ,tmax = length(lme_clim$sst)/12
                      ,tstepspryr  =  12
                      ,search_vol = 64.0
                      ,fmort.u = 0.5
                      ,fminx.u = 0
                      ,fmort.v = 0.5
                      ,fminx.v = -1
                      ,depth = mean(lme_clim$depth)
                      ,er = lme_clim$er
                      ,pp = lme_clim$intercept
                      ,slope = lme_clim$slope
                      ,sst = lme_clim$sst
                      ,sft = lme_clim$sbt
                      ,use.init = FALSE,effort = lme_clim$NomActive)      
  
  # run model through time
 
  result_set <- sizemodel(params) 
  
  # returns all outputs of the model 
  # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
  return(result_set)
}


# Error function
getError <-function(model=result_set,LMEnumber){
  
  # get outputs of model to compare with data
  # sum across all sizes fished and months to get total catch per year by functional group
  # converted to g ww per m^2, multiplying by depth assumption
  
  TotalUbiomass <- apply(result_set$U[params$ref:params$Nx,]*params$dx*10^params$x[params$ref:params$Nx],2,sum)*min(params$depth,100)
  TotalVbiomass <- apply(result_set$V[params$ref.det:params$Nx,]*params$dx*10^params$x[params$ref.det:params$Nx],2,sum)*min(params$depth,100)
  TotalW <- result_set$W[]*min(params$depth,100)
  
  #sum catches (in grams per m3 per month, across size classes) 
  
  TotalUcatch <- apply(result_set$Y.u[,]*params$dx,2,sum)*min(params$depth,100)
  TotalVcatch <- apply(result_set$Y.v[,]*params$dx,2,sum)*min(params$depth,100)
  
  
}


# Run Optim to estimate Catchability Param



