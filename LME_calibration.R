###############################################################################
# Library of functions developed originally by JLB on 21/12/2022
#
# Functions were updated by Denisse Fierro Arcos so they can be used with data
# produced by updated `01_getinputs_ISIMIP3a.R` script.
#
# Date of update: 2024-08-06

# Loading libraries
library(tidyverse)
library(lubridate)
library(zoo)
library(lhs)
library(pbapply)
library(patchwork)
source("dbpm_model_functions.R")

# library(parallel)
# library(optimParallel)
# library(gridExtra)


# This function reads in LME-scale time-series inputs (climate and fishing) 
# pre-processed 
# as spatially averaged means to estimate catchability (and if needed other
# model parameters)
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
      left_join(lme_fish_sub, by = c("year", "region")) |> 
      arrange(lat, lon)
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
  #Inputs:
  #vals (named numeric vector) - Single column with named rows containing LHS
  #parameters
  #input (data frame) - Forcing data produced by `get_lme_inputs` function
  #withinput (boolean) - Default is TRUE. ????
  #
  #Output:
  #If withinput set to TRUE:
  #input (data frame) - ????
  #If withinput set to FALSE:
  #result_set (data frame) - ???
  
  #Ensuring values are numeric
  if(!is.numeric(vals)){
    stop("Values provided in the 'vals' parameter must be numeric.")
  }
  
  # set up parameters
  params <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin.consumer.u = -3,
                      xmin.consumer.v = -3, tmax = nrow(input)/12,
                      tstepspryr = 12, search_vol = vals["search.vol"],
                      fmort.u = vals["f.u"], fminx.u = vals["f.minu"], 
                      fmort.v = vals["f.v"], fminx.v = vals["f.minv"], 
                      depth = mean(input$depth, na.rm = F), er = input$er,
                      pp = input$intercept, slope = input$slope, 
                      sst = input$sst, sft = input$sbt, use.init = FALSE,
                      effort = input$nom_active_area_m2_relative)
  
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
  
  if(withinput){
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
    
  }else{
    return(result_set)
  }
}


# Error function
getError <- function(lhs_params = X, lme_forcings = lme_input, 
                     figure_folder = NULL){
  #Inputs:
  #lhs_params (data frame) - Contains LHS parameters as columns and it must 
  #have a single row
  #lme_forcings (data frame) - Forcing data produced by `get_lme_inputs` 
  #function
  #figure_folder (character) - Optional. If provided, it must be the full path
  #to the folder where figures comparing observed and predicted data will be 
  #stored.
  #
  #Output:
  #rmse (numeric) - RMSE value between observed and predicted catch
  
  #Running model for specific LME
  result <- run_model(lhs_params, lme_forcings)
  
  #Aggregate data by year (mean to conserve units)
  error_calc <- result |> 
    filter(year > 1949) |> 
    group_by(year) |> 
    summarise(mean_total_catch_yr = mean(Totalcatch, na.rm = T),
              mean_obs_catch_yr = mean(catch_tonnes_area_m2, na.rm = T)) |> 
    #Converting units from tonnes to g
    mutate(mean_obs_catch_yr = mean_obs_catch_yr*1e6,
           #calculate and output error 
           #convert from tonnes to grams (m^-2*yr^-1)
           squared_error = (mean_obs_catch_yr-mean_total_catch_yr)^2)
  
  #Calculate RMSE
  sum_se <- sum(error_calc$squared_error, na.rm = T)
  count <- sum(!is.na(error_calc$squared_error))
  rmse <- sqrt(sum_se/count)
  
  #If a path to save figures is provided, create figures and save 
  if(!is.null(figure_folder)){
    if(!dir.exists(figure_folder)){
      dir.create(figure_folder)
    }
    #Plotting predicted and observed catches over time
    p1 <- ggplot()+
      geom_line(data = error_calc, aes(x = year, y = mean_total_catch_yr))+
      geom_point(data = error_calc, aes(x = year, y = mean_obs_catch_yr))+
      scale_x_continuous(breaks = seq(min(error_calc$year), 
                                      max(error_calc$year), by = 10))+
      theme_classic()+ 
      theme(axis.text.x = element_text(colour = "grey20", size = 12),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 15))+
      labs(x = "Year", y = bquote("Total catch (g*"~yr^-1*"*"*m^-2*")"))
    
    #Plotting predicted vs observed catches
    p2 <- ggplot()+
      geom_point(data = error_calc, 
                 aes(x = mean_total_catch_yr, y = mean_obs_catch_yr))+
      geom_abline(slope = 1, intercept = 0)+
      theme_classic()+
      theme(axis.text.x = element_text(colour = "grey20", size = 12),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 15))+
      labs(x = "Predicted", y = "Observed")
    
    #Creating a single plot
    p3 <- p1+p2+
      plot_annotation(title = paste0("LME #", unique(lme_forcings$region)),
                      theme = theme(plot.title = element_text(size = 16, 
                                                              hjust = 0.5)))
    
    #Creating path to save figure
    f_out <- file.path(figure_folder,
                       paste0("dbpm_pred_obs_catches_LME_", 
                              unique(lme_forcings$region), ".png"))
    
    #Saving composite figure
    ggsave(filename = f_out, plot = p3, width = 15, height = 10)
  }
  
  #Return RMSE
  return(rmse)
}


######## Carry out LHS param search
# now could try again with lhs instead of the regular grid of parameters
LHSsearch <- function(LMEnum = LME, num_iter = 1, search_vol = "estimated",
                      forcing_file = NULL, gridded_forcing = NULL, 
                      fishing_effort_file, figure_folder = NULL){
  #Inputs:
  #LMEnum (numeric) - Unique ID identifying LME
  #num_iter (numeric) - Number of individual runs. Default is 1.
  #search_vol (???) - Default is "estimated". ???
  #forcing_file (character) - Full path to forcing file. This must be 
  #non-gridded data
  #gridded_forcing (character) - Full path to folder containing gridded forcing
  #files
  #fishing_effort_file (character) - Full path to fishing effort file
  #figure_folder (character) - Optional. If provided, it must be the full path
  #to the folder where figures comparing observed and predicted data will be 
  #stored.
  #Output:
  #bestvals (data frame) - Contains the values for LHS parameters that resulted
  #in the best performing model based on RMSE values
  
  #Making function reproducible
  set.seed(1234)
  
  #Construct a hypercube with random numbers
  #num_iter defines number of rows in hypercube
  #columns represent five specific parameters needed
  sim <- data.frame(randomLHS(num_iter, 5))
  colnames(sim) <- c("f.u", "f.v", "f.minu", "f.minv", "search.vol")
  
  #Adjust range of mi size params, others go from 0-1
  sim <- sim |> 
    mutate(f.minu = f.minu*2, 
           f.minv = f.minv*2,
           # adjust range of search vol, others go from 0-1
           search.vol = search.vol+0.001)
  # runif(n=iter, min=0.064, max=1.0))
  
  if(is.numeric(search_vol)){
    sim <- sim |> 
      mutate(search.vol = search_vol)
  }
  
  # use below to select a constant value for search.vol
  # sim[,"search.vol"]<- 0.2
  if(!is.null(forcing_file)){
    lme_input <- get_lme_inputs(forcing_file = forcing_file, 
                                fishing_effort_file = fishing_effort_file, 
                                LMEnumber = LMEnum)
  }
  if(!is.null(gridded_forcing)){
    lme_input <- get_lme_inputs(gridded_forcing = gridded_forcing, 
                                fishing_effort_file = fishing_effort_file, 
                                LMEnumber = LMEnum)
  }
  
  # # check averaged and gridded inputs are the same in terms of effort 
  # lme_input[,c(2,15,16)]
  # lme_input_gridded<-get_lme_inputs(LMEnumber=LMEnum, gridded=T)
  # lme_input_gridded[,c(4,16,17)]
  
  # in pbapply setting cl = 6 calls mcapply to set up cluster 
  sim$rmse <- pbapply(sim, 1, getError, lme_forcings = lme_input, figure_folder,
                      cl = 6)
  
  # check this time param set with lowest error
  bestvals <- sim |> 
    filter(rmse == min(rmse, na.rm = T)) |> 
    mutate(region = LMEnum)
  
  #Print row with lowest RMSE
  print(bestvals)
  
  return(bestvals)
}


# Function to run model for each LME with gridded inputs, after 
# run_LME_calibration
run_model_timestep <- function(input = lme_inputs_igrid, 
                               vals = unlist(bestvals_LMEs[14, ]), U.initial,
                               V.initial, W.initial){
  
  # set up params for each month, across grid cells
  params <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin.consumer.u = -3,
                      xmin.consumer.v = -3, tmax = 1/12, tstepspryr = 12,
                      search_vol = 0.64, fmort.u = vals["f.u"], 
                      fminx.u = vals["f.minu"], fmort.v = vals["f.v"], 
                      fminx.v = vals["f.minv"],  depth = input["depth"],
                      er = input["er"], pp = input["intercept"], 
                      slope = input["slope"], sst = input["sst"], 
                      sft = input["sbt"], use.init = TRUE,
                      effort = input["total_nom_active_area_m2"], 
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
  library(parallel)
  library(optimParallel)
  
  # get input
  lme_input <- get_lme_inputs(LMEnumber = lme)
  
  # set up workers
  # keep some spare core
  noCores <- detectCores()-spareCores 
  if(noCores < 1){
    stop("You should allow at least one core for this operation.")
  }
  cl <- makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, varlist = c("cl", "libraries", "lme_input"),
                envir=environment())
  clusterEvalQ(cl, {
    for(item in 1:length(libraries)){
      library(libraries[item], character.only = T)
    }
  })
  
  clusterEvalQ(cl, source("LME_calibration.R"))
  
  optim_result <- optimParallel(par = vary, fn = errorFun,
                                method = "L-BFGS-B",
                                lower = rep(0, length(vary)),
                                upper = rep(2, length(vary)),
                                parallel = list(loginfo = TRUE,
                                                forward = TRUE), 
                                input = lme_input)
  stopCluster(cl)
  return(optim_result$par)
}





