# JLB - 21/12/2022
# This function reads in LME-scale time-series inputs (climate and fishing) pre-processed 
# as spatially averaged means to estimate catchability (and if needed other model parameters)

library(devtools)
library(tidyverse)
library(lubridate)

# source model script form Github - NOT WORKING, using local
# source_url(url = "https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm/blob/master/size-based-models/dynamic_sizebased_model_functions.R")
source(file = "dbpm_model_functions.R")


# Select LME and get inputs ready

get_lme_inputs<-function(LMEnumber=42){

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

#TO DO HERE: need to add a spin-up to these inputs prior to 1841 - 100 yrs at first value
return (lme_clim)
}

# Function to run model
run_model<-function(vals = X,input=lme_clim, LMEnumber=42){
  
  f.u<-as.numeric(vals[1])
  f.v<-as.numeric(vals[2])
   
  # set up params
  params <- sizeparam(equilibrium = FALSE
                      ,dx = 0.1
                      ,xmin.consumer.u = -3
                      ,xmin.consumer.v = -3
                      ,tmax = length(input$sst)/12
                      ,tstepspryr  =  12
                      ,search_vol = 0.64
                      ,fmort.u = f.u
                      ,fminx.u = 1.0
                      ,fmort.v = f.v
                      ,fminx.v = 1.0
                      ,depth = mean(input$depth)
                      ,er = input$er
                      ,pp = input$intercept
                      ,slope = input$slope
                      ,sst = input$sst
                      ,sft = input$sbt
                      ,use.init = FALSE,effort = input$NomActive)      
  
  # run model through time
  
  #TO DO IN SIZEMODEL CODE: make fishing function like one in model template
  
  result_set <- sizemodel(params) 
  
  # returns all outputs of the model 
  # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
  
  input$TotalUbiomass <- apply(result_set$U[params$ref:params$Nx,]*params$dx*10^params$x[params$ref:params$Nx],2,sum)*min(params$depth,100)
  input$TotalVbiomass <- apply(result_set$V[params$ref.det:params$Nx,]*params$dx*10^params$x[params$ref.det:params$Nx],2,sum)*min(params$depth,100)
  input$TotalW <- result_set$W[]*min(params$depth,100)
  
  #sum catches (currently in grams per m3 per year, across size classes) 
  # convert to grams per m2, then tonnes per m2 per year (for each month)
  
  input$TotalUcatch <- apply(result_set$Y.u[,]*params$dx,2,sum)*min(params$depth,100)/1e6
  input$TotalVcatch <- apply(result_set$Y.v[,]*params$dx,2,sum)*min(params$depth,100)/1e6
  input$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
  
  return(input)
  
}


# # Error function
getError <-function(vals = X,input,LMEnumber){

result<-run_model(vals,input,LMEnumber)

## aggregate by year (mean to conserve units)
out <- result %>% group_by(Year) %>% filter(Year > "1949") %>% summarise(TotalCatchPerYr=mean(Totalcatch),ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))

## calculate and output error
out$squared_error <- (out$ObsCatchPerYr- out$TotalCatchPerYr)^2

rmse<-sqrt(sum(out$squared_error,na.rm=T)/sum(!is.na(out$squared_error)))

return(rmse)

}




######## Carry out LHS param search

# now could try again with lhs instead of the regular grid of parameters
library(lhs)
library(pbapply)

set.seed(1234)

# num "individual runs"
num_iter=100


X <- data.frame(randomLHS(num_iter, 2))
# rows are iterations, columns are specific parameters
colnames(X)<-c("f.u","f.v")

#if the first column is met_coef and you want it to be bound by say, 0.25 and 0.5, as above you do this:

lme_input<-get_lme_inputs(LMEnumber = 42)
X$rmse<-pbapply(X,1, getError, input=lme_input, LMEnumber = 42)

saveRDS(X,file = "lhs_res_LME42_b.RDS")
#X<-readRDS(file = "lhs_res_LME42.RDS")

# check this time param set with lowest error
findmin<-which(X$rmse==min(X$rmse))
bestvals<-X[findmin,c(1,2)]

# check run
out<-run_model(X[findmin,c(1:2)],input=lme_input)
## aggregate by year (mean to conserve units)
out <- out %>% group_by(Year) %>% filter(Year > "1949") %>% summarise(TotalCatchPerYr=mean(Totalcatch),ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))

p1<-ggplot() +
  geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
  geom_point(data = out, aes(x = Year, y = ObsCatchPerYr)) +
  theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
                          axis.text.y = element_text(colour="grey20", size=12),
                          text=element_text(size=16)) + 
  labs(x = 'Year',
       y = 'Tonnes per year per m2') 

p2<-ggplot() +
  geom_point(data = out, aes(x = TotalCatchPerYr, y = ObsCatchPerYr)) +
  geom_abline(slope=1,intercept=0) +
  theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
                          axis.text.y = element_text(colour="grey20", size=12),
                          text=element_text(size=16)) + 
  labs(x = 'Predicted',
       y = 'Observed') 



#### Set up function to run optimParallel across LMEs
# # each run would need to get files for LME then do the error, return the best set of params, and make a diagnostic plot per LME
# # Run Optim to estimate Catchability Param
# # use best set of LHS results as initial values
# # 
# #

library(optimParallel)
# set up workers
noCores <- parallel::detectCores() - 1 # keep some spare core
cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, varlist = "cl",envir=environment())
clusterEvalQ(cl, {
  library(optimParallel)
})

optim_result <- optim(par=as.numeric(bestvals),getError, input=lme_input, LMEnumber=42)

stopCluster(cl)

