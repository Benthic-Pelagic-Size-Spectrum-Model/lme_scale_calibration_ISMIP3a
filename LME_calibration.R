# JLB - 21/12/2022
# This function reads in LME-scale time-series inputs (climate and fishing) 
# pre-processed 
# as spatially averaged means to estimate catchability (and if needed other
# model parameters)


library(tidyverse)
library(lubridate)
library(zoo)
library(lhs)

# library(devtools)
# library(parallel)
# library(tictoc)
# library(pbapply)
# library(patchwork)
# library(optimParallel)
# plots library
# library(rnaturalearth)
# library(sf)
# library(gridExtra)


# source model script from Github - NOT WORKING, using local
source("dbpm_model_functions.R")


#Testing
fishing_effort_file <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                                 "processed_forcings/lme_inputs/obsclim/025deg",
                                 "DBPM_LME_effort_catch_input.csv")

#Non-gridded
forcing_file <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                          "processed_forcings/lme_inputs/obsclim/025deg", 
                          "DBPM_LME_climate_inputs_slope.csv")

#Gridded
gridded_forcing <- file.path("/g/data/vf71/fishmip_inputs/ISIMIP3a",
                          "processed_forcings/lme_inputs_gridcell/obsclim",
                          "1deg")
  

# Select LME and get inputs ready
get_lme_inputs <- function(forcing_file = NULL, gridded_forcing = NULL, 
                           fishing_effort_file, LMEnumber, yearly = F){
  #Inputs:
  #forcing_file (character) - Full path to forcing file. This must be 
  #non-gridded data
  #gridded_forcing (character) - Full path to folder containing gridded forcing
  #files
  #fishing_effort_file (character) - Full path to fishing effort file
  #LMEnumber (numeric) - Unique ID identifying an LME
  #yearly (boolean) - Default is FALSE. If set to TRUE, it will return yearly
  #means for all forcing variables
  #
  #Output:
  #lme_clim (data frame) - Forcing data to be used in model calibration
  
  # Climate fishing inputs available via THREDDS:
  # http://portal.sf.utas.edu.au/thredds/catalog/gem/fishmip/ISIMIP3a/
  # InputData/DBPM_lme_inputs/obsclim/025deg/catalog.html?dataset=
  # fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/025deg/
  # DBPM_LME_effort_catch_input.csv
  lme_fish <- read_csv(fishing_effort_file) |> 
    # subset data 
    filter(region == LMEnumber) |> 
    # use relative effort - for each LME independently - PREFERRED
    mutate(nom_active_relative = total_nom_active/max(total_nom_active, 
                                                      na.rm = T),
           nom_active_area_m2_relative = total_nom_active_area_m2/
             max(total_nom_active_area_m2, na.rm = T))
  
  if(!is.null(forcing_file)){
    # Climate forcing inputs available via THREDDS:
    # http://portal.sf.utas.edu.au/thredds/catalog/gem/fishmip/ISIMIP3a/
    # InputData/DBPM_lme_inputs/obsclim/025deg/catalog.html?dataset=
    # fishmip/ISIMIP3a/InputData/DBPM_lme_inputs/obsclim/025deg/
    # DBPM_LME_climate_inputs_slope.csv
  
    lme_clim <- read_csv(forcing_file) |> 
      filter(region == LMEnumber)

    if(yearly){
      # option 1: keep monthly time steps but monthly value is the yearly mean
      
      # Calculating yearly average for all variables
      lme_clim <- lme_clim |> 
        mutate(year = year(t)) |> 
        group_by(year, region) |> 
        summarise(across(sst:expcbot, \(x) mean(x, na.rm = T))) |> 
        ungroup()
      
      # option 2: do as per global DBPML (builds on option 1): 
      # uses yearly average and extend the time series to monthly inputs
      # interpolate missing values using na.approx
      lme_clim <- lme_clim |> 
        #Yearly mean is copied 12 times (for each month in a year)
        slice(rep(1:n(), times = 12)) |> 
        #Adding month as column
        mutate(month = rep(1:12, each = nrow(lme_clim)),
               #Create a date column from year and month columns
               t = ym(paste(year, month, sep = "-"))) |> 
        select(!month) |> 
        arrange(t) |> 
        #Interpolate any missing values
        mutate(across(c(sst:lphy, expcbot), ~ na.approx(.x))) |> 
        mutate(year = as.numeric(year))
    }
    
    #Adding effort data
    lme_clim <- lme_clim |>  
      mutate(year = as.numeric(year(t))) |> 
      left_join(lme_fish, by = c("year", "region"))
  }
  
  else if(!is.null(gridded_forcing)){
    lme_clim <- list.files(gridded_forcing, 
                           paste0("LME_", LMEnumber, "_"), full.names = T) |> 
      read_csv() |> 
      #Extracting digits identifying region
      mutate(region = as.numeric(str_extract(region, "\\d{1,2}")))
   
    if(yearly){
      # option 1: keep monthly time steps but monthly value is the yearly mean
      # Calculating yearly average for all variables
      lme_clim <- lme_clim |> 
        mutate(year = year(t)) |> 
        group_by(year, region, lat, lon) |> 
        summarise(across(sst:expcbot, \(x) mean(x, na.rm = T))) |> 
        ungroup()
        
      # option 2: do as per global DBPM (builds on option 1): 
      # uses yearly average and extend the time series to monthly inputs
      # interpolate missing values using na.approx
      lme_clim <- lme_clim |> 
        #Yearly mean is copied 12 times (for each month in a year)
        slice(rep(1:n(), times = 12)) |> 
        #Adding month as column
        mutate(month = rep(1:12, each = nrow(lme_clim)),
               #Create a date column from year and month columns
               t = ym(paste(year, month, sep = "-"))) |> 
        select(!c(month)) |> 
        arrange(lat, lon, t) |> 
        #Interpolate any missing values
        mutate(across(c(sst:lphy, expcbot), ~ na.approx(.x)))
    }
    
    #Adding effort data
    lme_fish_sub <- lme_fish |>  
      select(!deptho) |>  
      distinct()
    
    lme_clim <- lme_clim |> 
      left_join(lme_fish_sub, by = c("year", "region"))
  }
  
  # Original comments from line 260
  # View(lme_clim)
  # convert fishing effort and catch per yr (divide values by 12)

  # ##### WARNING not sure this is OK Check With Julia.
  # # Catch and effort are yearly values and here repeated for each month. 
  # # Climate inputs are gm2 values (no time dimension, so??) repeated each
  # month.
  # # also here done considering aggregated inputs and runs (gridded = F, 
  # Yearly = F) so
  # # need to consider the other options too
  # lme_clim$NomActive<-lme_clim$NomActive/12
  # lme_clim$NomActive_area_m2<-lme_clim$NomActive_area_m2/12
  # lme_clim$catch_tonnes<-lme_clim$catch_tonnes/12
  # lme_clim$catch_tonnes_area_m2<-lme_clim$catch_tonnes_area_m2/12

  # try increasing effort as catches are too low - this works. 
  # lme_clim$NomActive<-lme_clim$NomActive*1000
  
  # select(lme_clim, c(t, NomActive, NomActive_area_m2, NomActive_relative, 
  # NomActive_area_m2_relative))
  
  # if (yearly!=T) {
  #   # could use a smoother to get intrannual variation working better
  # }

  #TO DO HERE: need to add a spin-up to these inputs prior to 1841 - 100 yrs 
  #at first value
  return(lme_clim)
}


# Function to run model
run_model <- function(vals = X, input = lme_input, withinput = T){
  #Ensure all columns are numeric
  vals <- vals |> 
    mutate(across(everything(), ~as.numeric(.x)))
  
  # set up parameters
  params <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin.consumer.u = -3,
                      xmin.consumer.v = -3, tmax = nrow(input)/12,
                      tstepspryr = 12, search_vol = vals$search.vol,
                      fmort.u = vals$f.u, fminx.u = vals$f.minu, 
                      fmort.v = vals$f.v, fminx.v = vals$f.minv, 
                      depth = mean(input$depth, na.rm = F), er = input$er,
                      pp = input$intercept, slope = input$slope, 
                      sst = input$sst, sft = input$sbt, use.init = FALSE,
                      effort = input$NomActive_area_m2_relative)
  
  # run model through time
  # TO DO IN SIZEMODEL CODE: make fishing function like one in model template
  
  result_set <- sizemodel(params) 
  
  # # CN CHECK 
  # NAs for size classes >90
  # dim(result_set$U)
  # result_set$U[,2040]
  # result_set$V[,10] # here is when it starts giving NAs
  
  # returns all outputs of the model 
  # saveRDS(result_set,
  #filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
  
  if(withinput == T){
    lims_ubio <- result_set$params$ref:result_set$params$Nx
    # JB:  changed inputs to m2 so no need to divide by depth here
    # added 1:params$Neq to same 2040 time steps instead of 2041
    time_steps <- 1:params$Neq
    
    input$TotalUbiomass <- apply(result_set$U[lims_ubio, time_steps] * 
                                   params$dx*10^params$x[lims_ubio], 
                                 2, sum)#*min(params$depth,100)
    
    lims_vbio <- result_set$params$ref.det:result_set$params$Nx
    input$TotalVbiomass <- apply(result_set$V[lims_vbio, time_steps] * 
                                   params$dx*10^params$x[lims_vbio],
                                 2, sum)#*min(params$depth,100)
    
    input$TotalW <- result_set$W[time_steps]#*min(params$depth,100)
    
    #sum catches (currently in grams per m3 per year, across size classes) 
    # keep as grams per m2, then be sure to convert observed from tonnes per m2
    # per year to g.^-m2.^-yr (for each month)
    
    input$TotalUcatch <- apply(result_set$Y.u[, time_steps]*params$dx, 2, sum)
    #*min(params$depth,100)
    input$TotalVcatch <- apply(result_set$Y.v[, time_steps]*params$dx, 2, sum)
    #*min(params$depth,100)
    input$Totalcatch <- input$TotalUcatch + input$TotalVcatch
  
   return(input)
  }
  
  if(withinput == F){
    return(result_set)
  }
}


# Error function
getError <- function(vals = X, input = lme_input, lme = NULL){
  # # trial
  # vals = sim
  # input=lme_input
  # lme = 1

  # if(!is.null(lme)) input<-get_lme_inputs(LMEnumber=lme)
  
  result <- run_model(vals, input)

  ## aggregate by year (mean to conserve units)
  out <- result |> 
    group_by(Year) |> 
    filter(Year > "1949") |> 0
    summarise(TotalCatchPerYr = mean(Totalcatch),
              ObsCatchPerYr = mean(catch_tonnes_area_m2, na.rm = T))

  ## convert units
  # CN from tonnes to g
  out$ObsCatchPerYr <- out$ObsCatchPerYr*1e6

  ## calculate and output error, convert observed from tonnes to grams (per m2
  # per yr)
  out$squared_error <- (out$ObsCatchPerYr- out$TotalCatchPerYr)^2

  rmse <- sqrt(sum(out$squared_error, na.rm = T) / 
                 sum(!is.na(out$squared_error)))
  return(rmse)
}


######## Carry out LHS param search

# now could try again with lhs instead of the regular grid of parameters
LHSsearch <- function(X = LME, iter = 1, search_vol = "estimated"){
  # # trial
  # X = 1
  # iter = 100
  # search_vol=0.64 # not estimated 
  LMEnum = X
  set.seed(1234)
  # num "individual runs"
  num_iter = iter
  sim <- data.frame(randomLHS(num_iter, 5))
  # rows are iterations, columns are specific parameters
  colnames(sim) <- c("f.u", "f.v", "f.minu", "f.minv", "search.vol")
  # adjust range of mi size params, others go form 0-1
  sim[, "f.minu"]<- sim[, "f.minu"]*2
  sim[, "f.minv"]<- sim[, "f.minv"]*2
  
  # adjust range of search vol, others go form 0-1
  # runif(n=iter, min=0.064, max=1.0)
  sim[, "search.vol"] <- sim[, "search.vol"]+0.001 
  
  if(is.numeric(search_vol)){
    sim[, "search.vol"]<- search_vol
    # use below to select a constant value for search.vol
    # sim[,"search.vol"]<- 0.2
  }
  
  lme_input <- get_lme_inputs(LMEnumber = LMEnum)
  
  # # check average and gridded inputs are the same in terms of effort 
  # lme_input[,c(2,15,16)]
  # lme_input_gridded<-get_lme_inputs(LMEnumber=LMEnum, gridded=T)
  # lme_input_gridded[,c(4,16,17)]
   
  # cl <- makeCluster(6)
  # #setDefaultCluster(cl = cl)
  # clusterExport(cl, varlist = c("LMEnum","lme_input","sim","iter"),
  # envir=environment())
  # clusterEvalQ(cl, {
  #   source("LME_calibration.R")
  #    })
  # 

  # in pbapply setting cl = 6 calls mcapply to set up cluster (so dont need 
  # above stuff)
  sim$rmse <- pbapply(sim, 1, getError, input = lme_input, cl = 6)
  
  #stopCluster(cl)

  #X<-readRDS(file = "lhs_res_LME42_b.RDS")
  # check this time param set with lowest error
  findmin <- which(sim$rmse == min(sim$rmse, na.rm = T))
  bestvals <- sim[findmin, ]
  print(bestvals[c(1:6)])
  # 
  # # check run
  # out<-run_model(bestvals[c(1:4)],input=lme_input)
  # ## aggregate by year (mean to conserve units)
  # out <- out |> group_by(Year) |> filter(Year > "1949") |> 
  #summarise(TotalCatchPerYr=mean(Totalcatch),
  #ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))
  # out$ObsCatchPerYr<-out$ObsCatchPerYr*1e6
  # 
  # p1<-ggplot() +
  #   geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
  #   geom_point(data = out, aes(x = Year, y = ObsCatchPerYr)) +
  #   theme_classic() + theme(axis.text.x = element_text(colour="grey20", 
  #size=12, angle=90, hjust=.5, vjust=.5),
  #                           axis.text.y = element_text(colour="grey20", 
  #size=12),
  #                           text=element_text(size=16)) + 
  #   labs(x = 'Year',
  #        y = 'Grams per year per m2') 
  # 
  # p2<-ggplot() +
  #   geom_point(data = out, aes(x = TotalCatchPerYr, y = ObsCatchPerYr)) +
  #   geom_abline(slope=1,intercept=0) +
  #   theme_classic() + theme(axis.text.x = element_text(colour="grey20", 
  #size=12, angle=90, hjust=.5, vjust=.5),
  #                           axis.text.y = element_text(colour="grey20", 
  #size=12),
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


# Function to run model for each LME with gridded inputs, after 
# run_LME_calibration
run_model_timestep <- function(input = lme_inputs_igrid, 
                               vals = unlist(bestvals_LMEs[14, ]), U.initial,
                               V.initial, W.initial){
  f.u <- as.numeric(vals[1])
  f.v <- as.numeric(vals[2])
  f.minu <- as.numeric(vals[3])
  f.minv <- as.numeric(vals[4])
  
  # set up params for each month, across grid cells
  params <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin.consumer.u = -3,
                      xmin.consumer.v = -3, tmax = 1/12, tstepspryr = 12,
                      search_vol = 0.64, fmort.u = f.u, fminx.u = f.minu,
                      fmort.v = f.v, fminx.v = f.minv, depth = input["depth"],
                      er = input["er"], pp = input["intercept"], 
                      slope = input["slope"], sst = input["sst"], 
                      sft = input["sbt"], use.init = TRUE,
                      effort = input["NomActive_area_m2"], 
                      U.initial = U.initial, V.initial = V.initial,
                      W.initial = W.initial)
  
  # run model for 1 timestep
  
  results <- sizemodel(params, use.init = T, ERSEM.det.input = F) 
  
  # input$TotalUbiomass <- apply(result_set$U[params$ref:params$Nx,]*
  #params$dx*10^params$x[params$ref:params$Nx],2,sum)*min(params$depth,100)
  # input$TotalVbiomass <- apply(result_set$V[params$ref.det:params$Nx,]*
  #params$dx*10^params$x[params$ref.det:params$Nx],2,sum)*min(params$depth,100)
  # input$TotalW <- result_set$W[]*min(params$depth,100)
  # 
  #sum catches (currently in grams per m3 per year, across size classes) 
  # keep as grams per m2, then be sure to convert observed from tonnes per m2 
  # per year to g.^-m2.^-yr (for each month)
  
  # input$TotalUcatch <- apply(result_set$Y.u[,]*params$dx,2,sum)*
  #min(params$depth,100)
  # input$TotalVcatch <- apply(result_set$Y.v[,]*params$dx,2,sum)*
  #min(params$depth,100)
  # input$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
  # 
  return(cbind(U = results$U[, 2], V = results$V[, 2], Y.u = results$Y.u[, 2],
               Y.v = results$Y.v[, 2], GG.u = results$GG.u[, 1],
               GG.v = results$GG.v[, 1], PM.u = results$PM.u[, 1],
               PM.v = results$PM.v[, 1], W = results$W[2]))
  # return next time step value for that grid
  }


### this needs to be faster or maybe just part of outputs of gridded_sizemodel
getGriddedOutputs <- function(input = lme_inputs_grid, results = grid_results,
                              params = params){
  # returns all outputs of the model 
  # saveRDS(result_set,filename=
  #paste("dbpm_calibration_LMEnumber_catchability.rds"))
  input$TotalUbiomass <- input$TotalVbiomass <- input$TotalUcatch <-
    input$TotalVcatch<- input$Totalcatch <- NA
  cells <- unique(input$cell)
  
  # remove depth adjustment as run inn /m2 now
  # depthadj<-ifelse(params$depth>100,100,params$depth)
  
  for(igrid in 1:length(cells)){
    input[input$cell == cells[igrid], ]$TotalUbiomass <- 
      apply(results$U[igrid, params$ref:params$Nx, 2:(params$Neq+1)] * 
              params$dx*10^params$x[params$ref:params$Nx],
            2, sum, na.rm = T)#*depthadj[igrid]
    
    input[input$cell == cells[igrid], ]$TotalVbiomass <- 
      apply(results$V[igrid, params$ref:params$Nx, 2:(params$Neq+1)] * 
              params$dx*10^params$x[params$ref:params$Nx],
            2, sum, na.rm = T)#*depthadj[igrid]
    
    # input[input$t==time[itime]]$W <- results$W[,itime]*min(params$depth,100)
    
    #sum catches (currently in grams per m3 per year, across size classes) 
    #keep as grams per m2, then be sure to convert observed from tonnes per m2 
    #per year to g.^-m2.^-yr (for each month)
    input[input$cell == cells[igrid], ]$TotalUcatch <- 
      apply(results$Y.u[igrid, params$ref:params$Nx, 2:(params$Neq+1)] *
              params$dx, 2, sum, na.rm = T)#*depthadj[igrid]
    
    input[input$cell == cells[igrid], ]$TotalVcatch <- 
      apply(results$Y.v[igrid, params$ref:params$Nx, 2:(params$Neq+1)] *
              params$dx, 2, sum, na.rm = T)#*depthadj[igrid]
    
    input[input$cell == cells[igrid], ]$Totalcatch <- 
      input[input$cell == cells[igrid], ]$TotalUcatch +
      input[input$cell == cells[igrid], ]$TotalVcatch
  }
  ## and then multiply outputs by depth to get per m2
  return(input)
}
  
### fastOptim to set up and run optimisations in parallel
fastOptim <- function(lme, vary, errorFun = getError, spareCores = 1,
                      libraries = c("optimParallel")){
  # get input
  lme_input <- get_lme_inputs(LMEnumber = lme)
  
  # set up workers
  # keep some spare core
  noCores <- parallel::detectCores()-spareCores 
  if(noCores < 1){
    stop("You should allow at least one core for this operation.")
  }
  cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
  parallel::setDefaultCluster(cl = cl)
  parallel::clusterExport(cl, 
                          varlist = c("cl", "libraries", "lme_input"),
                          envir=environment())
  parallel::clusterEvalQ(cl, {
    for(item in 1:length(libraries)){
      library(libraries[item], character.only = T)
    }
  })
  parallel::clusterEvalQ(cl, source("LME_calibration.R"))
  
  optim_result <- optimParallel::optimParallel(par = vary, fn = errorFun,
                                               method = "L-BFGS-B",
                                               lower = rep(0, length(vary)),
                                               upper = rep(2, length(vary)),
                                               parallel = list(loginfo = TRUE,
                                                               forward = TRUE), 
                                               input = lme_input)
  parallel::stopCluster(cl)
  return(optim_result$par)
}





