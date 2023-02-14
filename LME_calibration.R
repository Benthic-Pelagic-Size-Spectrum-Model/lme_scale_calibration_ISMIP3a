# JLB - 21/12/2022
# This function reads in LME-scale time-series inputs (climate and fishing) pre-processed 
# as spatially averaged means to estimate catchability (and if needed other model parameters)

library(devtools)
library(tidyverse)
library(lubridate)
library(parallel)
library(tictoc)
library(lhs)
library(pbapply)
library(patchwork)
library(optimParallel)

# source model script form Github - NOT WORKING, using local
# source_url(url = "https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm/blob/master/size-based-models/dynamic_sizebased_model_functions.R")
source(file = "dbpm_model_functions.R")



# Select LME and get inputs ready

get_lme_inputs<-function(LMEnumber=14, gridded=F){

if (gridded!=T) {
# read climate forcing inputs from THREDDS
lme_clim<-read_csv(file="http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_climate_inputs_slope.csv")
lme_clim<-subset(lme_clim, LME %in% LMEnumber)
}
  
if (gridded==T) {
    # read climate forcing inputs from gem, here testing on LME14
    #filename =paste("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/lme_inputs_gridcell/obsclim/1deg/observed_LME_",LMEnumber,".csv",sep="") 
    filename="observed_LME_14.csv"
    lme_clim<-read_csv(file=filename)
          }
  
  
# read climate fishing inputs from THREDDS
lme_fish<-read_csv(file="http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_effort_catch_input.csv")
lme_fish<-subset(lme_fish, LME %in% LMEnumber)

#create monthly inputs for fishing
lme_clim$Year<-year(lme_clim$t)
lme_clim<-left_join(lme_clim,lme_fish,by="Year")
# convert fishing effort and catch per yr (divide values by 12)
lme_clim$NomActive_area_m2
lme_clim$catch_tonnes_area_m2

#TO DO HERE: need to add a spin-up to these inputs prior to 1841 - 100 yrs at first value
return (lme_clim)
}



# Function to run model
run_model<-function(vals = X,input=lme_input,withinput=T){
  
  f.u<-as.numeric(vals[1])
  f.v<-as.numeric(vals[2])
  f.minu<-as.numeric(vals[3])
  f.minv<-as.numeric(vals[4])
   
  # set up params
  params <- sizeparam(equilibrium = FALSE
                      ,dx = 0.1
                      ,xmin.consumer.u = -3
                      ,xmin.consumer.v = -3
                      ,tmax = length(input$sst)/12
                      ,tstepspryr  =  12
                      ,search_vol = 0.64
                      ,fmort.u = f.u
                      ,fminx.u = f.minu
                      ,fmort.v = f.v
                      ,fminx.v = f.minv
                      ,depth = mean(input$depth)
                      ,er = input$er
                      ,pp = input$intercept
                      ,slope = input$slope
                      ,sst = input$sst
                      ,sft = input$sbt
                      ,use.init = FALSE,effort = input$NomActive_area_m2)      
  
  # run model through time
  
  #TO DO IN SIZEMODEL CODE: make fishing function like one in model template
  
  result_set <- sizemodel(params) 
  
  # returns all outputs of the model 
  # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
  
  if (withinput==T){
  input$TotalUbiomass <- apply(result_set$U[params$ref:params$Nx,]*params$dx*10^params$x[params$ref:params$Nx],2,sum)*min(params$depth,100)
  input$TotalVbiomass <- apply(result_set$V[params$ref.det:params$Nx,]*params$dx*10^params$x[params$ref.det:params$Nx],2,sum)*min(params$depth,100)
  input$TotalW <- result_set$W[]*min(params$depth,100)
  
  #sum catches (currently in grams per m3 per year, across size classes) 
  # keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
  
  input$TotalUcatch <- apply(result_set$Y.u[,]*params$dx,2,sum)*min(params$depth,100)
  input$TotalVcatch <- apply(result_set$Y.v[,]*params$dx,2,sum)*min(params$depth,100)
  input$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
  
  return(input)
  }
  
  if (withinput==F) {
    return(result_set)
  }
  
}


# # Error function
getError <-function(vals = X,input=lme_input){

result<-run_model(vals, input)

## aggregate by year (mean to conserve units)
out <- result %>% group_by(Year) %>% filter(Year > "1949") %>% summarise(TotalCatchPerYr=mean(Totalcatch),ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))

## convert units

out$ObsCatchPerYr<-out$ObsCatchPerYr*1e6

## calculate and output error, convert observed from tonnes to grams (per m2 per yr)
out$squared_error <- (out$ObsCatchPerYr- out$TotalCatchPerYr)^2

rmse<-sqrt(sum(out$squared_error,na.rm=T)/sum(!is.na(out$squared_error)))

return(rmse)

}


######## Carry out LHS param search

# now could try again with lhs instead of the regular grid of parameters

LHSsearch<-function(X=LME,iter=1) {

  LMEnum=X
  set.seed(1234)
  # num "individual runs"
  num_iter=iter
  sim <- data.frame(randomLHS(num_iter, 4))
  # rows are iterations, columns are specific parameters
  colnames(sim)<-c("f.u","f.v","f.minu","f.minv")
  # adjust range of mi size params, others go form 0-1
  sim[,"f.minu"]<- sim[,"f.minu"]*2
  sim[,"f.minv"]<- sim[,"f.minv"]*2
  
  lme_input<-get_lme_inputs(LMEnumber=LMEnum)
  
  
   
# cl <- makeCluster(6)
# #setDefaultCluster(cl = cl)
# clusterExport(cl, varlist = c("LMEnum","lme_input","sim","iter"),envir=environment())
# clusterEvalQ(cl, {
#   source("LME_calibration.R")
#    })
# 

# in pbapply setting cl = 6 calls mcapply to set up cluster (so dont need above stuff)
sim$rmse<-pbapply(sim,1, getError, input=lme_input,cl=6)
  
#stopCluster(cl)

  #X<-readRDS(file = "lhs_res_LME42_b.RDS")
  # check this time param set with lowest error
  findmin<-which(sim$rmse==min(sim$rmse))
  bestvals<- sim[findmin,]
  print(bestvals[c(1:5)])
  # 
  # # check run
  # out<-run_model(bestvals[c(1:4)],input=lme_input)
  # ## aggregate by year (mean to conserve units)
  # out <- out %>% group_by(Year) %>% filter(Year > "1949") %>% summarise(TotalCatchPerYr=mean(Totalcatch),ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))
  # out$ObsCatchPerYr<-out$ObsCatchPerYr*1e6
  # 
  # p1<-ggplot() +
  #   geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
  #   geom_point(data = out, aes(x = Year, y = ObsCatchPerYr)) +
  #   theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
  #                           axis.text.y = element_text(colour="grey20", size=12),
  #                           text=element_text(size=16)) + 
  #   labs(x = 'Year',
  #        y = 'Grams per year per m2') 
  # 
  # p2<-ggplot() +
  #   geom_point(data = out, aes(x = TotalCatchPerYr, y = ObsCatchPerYr)) +
  #   geom_abline(slope=1,intercept=0) +
  #   theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
  #                           axis.text.y = element_text(colour="grey20", size=12),
  #                           text=element_text(size=16)) + 
  #   labs(x = 'Predicted',
  #        y = 'Observed',title=paste("LME ",LMEnum,sep=""))
  # 
  # p3<-p1+p2
  # 
  # ggsave(paste("dbpm_LME",LMEnum,".png",sep=""), width=15, height=10)
  # 
  # 
  return(bestvals)

}




# Function to run model for each LME with gridded inputs, after run_LME_calibration
run_model_timestep<-function(input=lme_inputs_igrid, vals = unlist(bestvals_LMEs[14,]),U.initial, V.initial, W.initial){
  
  f.u<-as.numeric(vals[1])
  f.v<-as.numeric(vals[2])
  f.minu<-as.numeric(vals[3])
  f.minv<-as.numeric(vals[4])
  
  # set up params for each month, across grid cells
  params <- sizeparam (equilibrium = FALSE
                      ,dx = 0.1
                      ,xmin.consumer.u = -3
                      ,xmin.consumer.v = -3
                      ,tmax = 1/12
                      ,tstepspryr  =  12
                      ,search_vol = 0.64
                      ,fmort.u = f.u
                      ,fminx.u = f.minu
                      ,fmort.v = f.v
                      ,fminx.v = f.minv
                      ,depth = input["depth"]
                      ,er = input["er"]
                      ,pp = input["intercept"]
                      ,slope = input["slope"]
                      ,sst = input["sst"]
                      ,sft = input["sbt"]
                      ,use.init = TRUE,effort = input["NomActive_area_m2"], U.initial =U.initial,V.initial = V.initial,W.initial = W.initial)      
  
  # run model for 1 timestep
  
  results <- sizemodel(params,use.init = T,ERSEM.det.input = F) 
  
  # input$TotalUbiomass <- apply(result_set$U[params$ref:params$Nx,]*params$dx*10^params$x[params$ref:params$Nx],2,sum)*min(params$depth,100)
  # input$TotalVbiomass <- apply(result_set$V[params$ref.det:params$Nx,]*params$dx*10^params$x[params$ref.det:params$Nx],2,sum)*min(params$depth,100)
  # input$TotalW <- result_set$W[]*min(params$depth,100)
  # 
  #sum catches (currently in grams per m3 per year, across size classes) 
  # keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
  
  # input$TotalUcatch <- apply(result_set$Y.u[,]*params$dx,2,sum)*min(params$depth,100)
  # input$TotalVcatch <- apply(result_set$Y.v[,]*params$dx,2,sum)*min(params$depth,100)
  # input$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
  # 
  return(cbind(U=results$U[,2],V=results$V[,2],Y.u=results$Y.u[,2],Y.v=results$Y.v[,2],GG.u=results$GG.u[,1],GG.v=results$GG.v[,1],PM.u=results$PM.u[,1],PM.v=results$PM.v[,1],W=results$W[2]))
  # return next time step value for that grid
  }
  
  


