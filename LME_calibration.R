# JLB - 21/12/2022
# This function reads in LME-scale time-series inputs (climate and fishing) pre-processed 
# as spatially averaged means to estimate catchability (and if needed other model parameters)

# rm(list=ls())
library(devtools)
library(tidyverse)
library(lubridate)
library(parallel)
library(tictoc)
library(lhs)
library(pbapply)
library(patchwork)
library(optimParallel)
# plots library
library(rnaturalearth)
library(sf)
library(gridExtra)

# source model script form Github - NOT WORKING, using local
# source_url(url = "https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm/blob/master/size-based-models/dynamic_sizebased_model_functions.R")
source(file = "dbpm_model_functions.R")

# Select LME and get inputs ready

get_lme_inputs<-function(LMEnumber=14, gridded=F,yearly=F){

  # # trial
  # LMEnumber=1
  # gridded=F
  # yearly=F

if (!gridded) {

 # read climate forcing inputs from THREDDS
 # lme_clim<-read_csv(file="http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_climate_inputs_slope.csv")
  
  # # trial 
  # LMEnumber = 1
  # gridded = F
  # yearly = F
  
  lme_clim<-read_csv(file="/rd/gem/public/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_climate_inputs_slope.csv")
 
  lme_clim<-subset(lme_clim, LME %in% LMEnumber)

  if (yearly){
    
    # option 1: keep monthly time steps but each month value is the yearly average 
    
    lme_clim<-lme_clim %>% 
      mutate(year = sub("\\-.*", "", t))
    
    time<-lme_clim %>% 
      select(t, year) %>% 
      distinct()
    
    lme_clim<-lme_clim %>%
      gather(key, value, -c(LME, t, year)) %>% 
      group_by(key, year, LME) %>% 
      summarise(mean = mean(value)) %>% 
      ungroup() %>% 
      spread(key, mean)

    lme_clim<-lme_clim %>% 
      full_join(time) %>% 
      arrange(t)
    
    # option 2: do as per global DBPML (builts on option 1): 
    # use yearly average for January each year
    # extend the time series to weekly inputs
    # interpolate missing values (from Jan year 1 to jan year 2 etc) using na.omit
    
    # for last year you should consider December instead as the run ends in Dec 2010
    # Because these monthly values are yearly avarage, January == December

    lme_clim<-lme_clim %>% # one value each year (yearly average) 
      mutate(month = month(t)) %>%  # consider only January
      filter(month %in% c(1)) %>% 
      select(-month)
    
    final_month<-lme_clim$t[nrow(lme_clim)]
    
    # lme_clim<-lme_clim %>% 
    #   mutate(t = as.Date(ifelse(t == final_month, as.Date("2010-12-01"), t)))

    # if above not working, do this instead:
    lme_clim<-lme_clim %>%
      mutate(t = as.character(t),
             t = ifelse(t == final_month, "2010-12-01", t),
             t = as.Date(t))
    
    # create a weekly time vector 
    extend<-data.frame(t = seq(lme_clim$t[1], lme_clim$t[nrow(lme_clim)], by="month")) %>% 
      mutate(LME = LMEnumber, 
             year = as.character(year(t)), 
             area_m2 = unique(lme_clim$area_m2), 
             depth = unique(lme_clim$depth))
    
    library(zoo) # na.approx as in dbpm global
    lme_clim<-lme_clim %>% 
      select(-t) %>% 
      full_join(extend) %>% 
      arrange(t) %>% 
      mutate(er = na.approx(er), 
             expcbot = na.approx(expcbot),
             intercept = na.approx(intercept),
             lphy = na.approx(lphy),
             sbt = na.approx(sbt),
             slope = na.approx(slope),
             sphy = na.approx(sphy),
             sst = na.approx(sst))
  
  }

}
  
if (gridded) {
  
    # read climate forcing inputs from gem, here testing on LME14
    # trial 
    # CN WARNING the 0.25deg inuts have not been produced yet (27/07/2023) as we were focusing on 1deg first  
    filename =paste("/rd/gem/private/fishmip_inputs/ISIMIP3a/processed_forcings/lme_inputs_gridcell/obsclim/1deg/observed_LME_",LMEnumber,".csv",sep="")
    # filename="observed_LME_14.csv"
    lme_clim<-read_csv(file=filename)

    # multiple grids mean multiple area_m2 ... 
    unique(lme_clim$area_m2)
    # and multiple depths
    unique(lme_clim$depth)
    
    if (yearly){

      # create a key - this needs to be the final lat/lon/t format
      key<-lme_clim %>% 
        select(lat,lon, depth, area_m2, t) %>% 
        distinct()
      
      # crate a second key - this is where you store grid-cell specific and not time-variant info
      key2<-key %>% 
        select(-t) %>% 
        distinct()
      
      # option 1: keep monthly time steps but each month value is the yearly average 
      
      lme_clim<-lme_clim %>% 
        mutate(year = sub("\\-.*", "", t))
      
      time<-lme_clim %>% 
        select(t, year) %>% 
        distinct()
      
      lme_clim<-lme_clim %>% 
        gather(key, value, -c(LME, t, year, lat, lon)) %>% 
        group_by(key, year, LME, lat, lon) %>% 
        summarise(mean = mean(value)) %>% 
        ungroup() %>% 
        spread(key, mean)
      
      lme_clim<-lme_clim %>% 
        full_join(time) %>% 
        arrange(lat,lon,t) # according to key above 
      
      # option 2: do as per global DBPML (builts on option 1): 
      # use yearly average for January each year
      # extend the time series to weekly inputs
      # interpolate missing values (from Jan year 1 to jan year 2 etc) using na.omit
      
      # for last year you should consider December instead as the run ends in Dec 2010
      # Because these monthly values are yearly average, January == December
      
      lme_clim<-lme_clim %>% # one value for Jan each year (yearly average) 
        mutate(month = month(t)) %>%
        filter(month %in% c(1)) %>% 
        select(-month)
      
      final_month<-lme_clim$t[nrow(lme_clim)]
      
      lme_clim<-lme_clim %>%
        mutate(t = as.character(t),
               t = ifelse(t == final_month, "2010-12-01", t),
               t = as.Date(t))
      # check 
      # lme_clim$t[nrow(lme_clim)]
      
      # create a third key - with year instead of t  
      key3<-key %>% 
        mutate(year = as.character(year(t))) %>% 
        select(-c(t)) %>%
        distinct()
      
      # create a weekly time vector 
      extend<-data.frame(t = seq(lme_clim$t[1], lme_clim$t[nrow(lme_clim)], by="month")) %>% 
        mutate(LME = LMEnumber, 
               year = as.character(year(t))) %>% 
        full_join(key3) %>% 
        arrange(lat,lon, t)
      
      library(zoo) # na.approx as in dbpm global
      lme_clim<-lme_clim %>%
        select(-t) %>% # Need to delete t here otherwise the join does not work properly
        full_join(extend) %>% 
        arrange(lat,lon, t) %>%
        group_by(lat, lon, LME, area_m2, depth) %>% 
        summarise(er = na.approx(er), 
               expcbot = na.approx(expcbot),
               intercept = na.approx(intercept),
               lphy = na.approx(lphy),
               sbt = na.approx(sbt),
               slope = na.approx(slope),
               sphy = na.approx(sphy),
               sst = na.approx(sst)) %>% 
        ungroup()
    
      # add t
      lme_clim<-lme_clim %>% 
        cbind("t" = extend$t)
      
    }
    
 }
  
  # read climate fishing inputs from THREDDS
  #lme_fish<-read_csv(file="http://portal.sf.utas.edu.au/thredds/fileServer/gem/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_effort_catch_input.csv")
  lme_fish<-read_csv(file="/rd/gem/public/fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/0.25deg/DBPM_LME_effort_catch_input.csv")
  
  # check dataset to se if you can calcualte relative effort here or need to do in 01_getgriddedinputs.R
  # head(lme_fish)
  
  # use relative effort - across all LMEs 
  #lme_fish$NomActive_area_m2_relative<-lme_fish$NomActive_area_m2/max(lme_fish$NomActive_area_m2)
  
  # subset data 
  lme_fish<-subset(lme_fish, LME %in% LMEnumber)
  
  # # or use relative effort - for each LME independently 
  lme_fish$NomActive_relative<-lme_fish$NomActive/max(lme_fish$NomActive) # same as below
  lme_fish$NomActive_area_m2_relative<-lme_fish$NomActive_area_m2/max(lme_fish$NomActive_area_m2)
  
  #create monthly inputs for fishing
  lme_clim$Year<-year(lme_clim$t)

  if(gridded == FALSE){
    # here we don't have multiple grids - so only 1 (weighted.mean) depth and 1 (weighted.mean) area 
    # which should be the same for lme_clim and lme_fish - NEED TO CHECK 
    lme_clim<-lme_clim %>% 
      left_join(lme_fish, by=c ("Year", "LME", "area_m2"))
    }else if(gridded == TRUE){
    # here we have multiple grids - so multiple depth and area
    # need to remove area and depth before merging
    # NomActive_area_m2 and catch_tonnes_area_m2 are already /m2 and hence work at the grid cell level too
    lme_fish<-lme_fish %>% 
      select(-c(area_m2, deptho_m)) %>% 
      unique()
    lme_clim<-lme_clim %>% 
      left_join(lme_fish,by=c ("Year", "LME")) 
    
  }
    
  # View(lme_clim)
  # convert fishing effort and catch per yr (divide values by 12)

  # ##### WARNING not sure this is OK Check With Julia.
  # # Catch and effort are yearly values and here repeated for each months. 
  # # Climate inputs are gm2 values (no time dimension, so??) repeated each month.
  # # also here done considering aggregated inputs and runs (gridded = F, Yearly = F) so
  # # need to consider the other options too
  # lme_clim$NomActive<-lme_clim$NomActive/12
  # lme_clim$NomActive_area_m2<-lme_clim$NomActive_area_m2/12
  # lme_clim$catch_tonnes<-lme_clim$catch_tonnes/12
  # lme_clim$catch_tonnes_area_m2<-lme_clim$catch_tonnes_area_m2/12

  # try increasing effort as catches are too low - this works. 
  # lme_clim$NomActive<-lme_clim$NomActive*1000
  
  # select(lme_clim, c(t, NomActive, NomActive_area_m2, NomActive_relative, NomActive_area_m2_relative))
  
  # if (yearly!=T) {
  #   # could use a smoother to get intrannual variation working better
  # }

  #TO DO HERE: need to add a spin-up to these inputs prior to 1841 - 100 yrs at first value
  return (lme_clim)
}



# Function to run model
run_model<-function(vals = X,input=lme_input,withinput=T){
  
  # # CN trial
  # vals = sim
  # input=lme_input
  # withinput=T
  
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
                      ,search_vol = vals[5]
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
                      ,use.init = FALSE
                      ,effort = input$NomActive_area_m2_relative)      
  
  # run model through time

  # TO DO IN SIZEMODEL CODE: make fishing function like one in model template
  
  result_set <- sizemodel(params) 
  
  # # CN CHECK 
  # NAs for size classes >90
  # dim(result_set$U)
  # result_set$U[,2040]
  # result_set$V[,10] # here is when it starts giving NAs
  
  # returns all outputs of the model 
  # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
  
  if (withinput==T){
  
  # JB:  changed inputs to  m2 so no need to divide by depth here
  # also added 1:params$Neq to same 2040 time steps instead of 2041
  input$TotalUbiomass <- apply(result_set$U[params$ref:params$Nx,1:params$Neq]*params$dx*10^params$x[params$ref:params$Nx],2,sum)#*min(params$depth,100)
  input$TotalUbiomass 
  
  input$TotalVbiomass <- apply(result_set$V[params$ref.det:params$Nx,1:params$Neq]*params$dx*10^params$x[params$ref.det:params$Nx],2,sum)#*min(params$depth,100)
  input$TotalW <- result_set$W[1:params$Neq]#*min(params$depth,100)
  
  #sum catches (currently in grams per m3 per year, across size classes) 
  # keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
  
  input$TotalUcatch <- apply(result_set$Y.u[,1:params$Neq]*params$dx,2,sum)#*min(params$depth,100)
  input$TotalVcatch <- apply(result_set$Y.v[,1:params$Neq]*params$dx,2,sum)#*min(params$depth,100)
  input$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
  
  return(input)
  }
  
  if (withinput==F) {
    return(result_set)
  }
  
}


# # Error function
getError <-function(vals = X,input=lme_input){

  # # trial 
  # vals = sim
  # input=lme_input
  
  result<-run_model(vals, input)

  ## aggregate by year (mean to conserve units)
  out <- result %>% 
    group_by(Year) %>% 
    filter(Year > "1949") %>% 
    summarise(TotalCatchPerYr=mean(Totalcatch),
              ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))

  ## convert units
  # CN from tonnes to g
  out$ObsCatchPerYr<-out$ObsCatchPerYr*1e6

  ## calculate and output error, convert observed from tonnes to grams (per m2 per yr)
  out$squared_error <- (out$ObsCatchPerYr- out$TotalCatchPerYr)^2

  rmse<-sqrt(sum(out$squared_error,na.rm=T)/sum(!is.na(out$squared_error)))

  return(rmse)

}


######## Carry out LHS param search

# now could try again with lhs instead of the regular grid of parameters

LHSsearch<-function(X=LME,iter=1) {

  # # trial
  # X = 1
  # iter = 10
  
  LMEnum=X
  set.seed(1234)
  # num "individual runs"
  num_iter=iter
  sim <- data.frame(randomLHS(num_iter,5))
  # rows are iterations, columns are specific parameters
  colnames(sim)<-c("f.u","f.v","f.minu","f.minv","search.vol")
  # adjust range of mi size params, others go form 0-1
  sim[,"f.minu"]<- sim[,"f.minu"]*2
  sim[,"f.minv"]<- sim[,"f.minv"]*2
  sim[,"search.vol"]<- sim[,"search.vol"] + 0.001
  
  lme_input<-get_lme_inputs(LMEnumber=LMEnum)
  
  # # check averadge and gridded inputs are the same in terms of effort 
  # lme_input[,c(2,15,16)]
  # lme_input_gridded<-get_lme_inputs(LMEnumber=LMEnum, gridded=T)
  # lme_input_gridded[,c(4,16,17)]
   
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
  print(bestvals[c(1:6)])
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
  


### this needs to be faster or maybe just part of outputs of gridded_sizemodel
getGriddedOutputs<-function(input = lme_inputs_grid,
                            results = grid_results,
                            params = params){
  # returns all outputs of the model 
  # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
  input$TotalUbiomass <-   input$TotalVbiomass <- input$TotalUcatch <- input$TotalVcatch<- input$Totalcatch <- NA  
  cells<-unique(input$cell)
  
  # remove depth adjustment as run inn /m2 now
  # depthadj<-ifelse(params$depth>100,100,params$depth)
  
  for(igrid in 1:length(cells)) {
    
    input[input$cell==cells[igrid],]$TotalUbiomass <- apply(results$U[igrid,params$ref:params$Nx,2:(params$Neq+1)]*params$dx*10^params$x[params$ref:params$Nx],
                                                            2,sum,na.rm=T)#*depthadj[igrid]
    input[input$cell==cells[igrid],]$TotalVbiomass <- apply(results$V[igrid,params$ref:params$Nx,2:(params$Neq+1)]*params$dx*10^params$x[params$ref:params$Nx],
                                                            2,sum,na.rm=T)#*depthadj[igrid]
    # input[input$t==time[itime]]$W <- results$W[,itime]*min(params$depth,100)
    #sum catches (currently in grams per m3 per year, across size classes) 
    #keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
    input[input$cell==cells[igrid],]$TotalUcatch <- apply(results$Y.u[igrid,params$ref:params$Nx,2:(params$Neq+1)]*params$dx,2,sum,na.rm=T)#*depthadj[igrid]
    input[input$cell==cells[igrid],]$TotalVcatch <- apply(results$Y.v[igrid,params$ref:params$Nx,2:(params$Neq+1)]*params$dx,2,sum,na.rm=T)#*depthadj[igrid]
    input[input$cell==cells[igrid],]$Totalcatch <- input[input$cell==cells[igrid],]$TotalUcatch +   input[input$cell==cells[igrid],]$TotalVcatch
  }
  
  ## and then multiply outputs by depth to get per m2
  
  return(input)
  
}
  


