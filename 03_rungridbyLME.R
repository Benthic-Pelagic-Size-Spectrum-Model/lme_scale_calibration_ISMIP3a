#' @title rungridbyLME
#' @description run model across space and time using gridded inputs for each LME
#' @param LMEnumber The LME number you want to run
#' @param yearly Correspond to the yearly argument in the get_lme_inputs function 
#' after initialisation. Default is FALSE
#' @param f.effort Boolean value. If true, fisheries effort will be set using 
#' bestvals_LME.rds. If false, fisheries effort is set to 0. Default is TRUE.
#' @param search_vol Numeric value for the search volume. Default is .64.
#' @param savePlots Boolean value. If true, plots related to the selected LME will
#' be saved in the same folder as the grid output. Default is TRUE.
#' 
#' 

#### set environment ----

rm(list=ls())

library(tictoc)
library(raster)
library(tidyverse)
library(data.table)

select<-dplyr::select
summarise <-dplyr::summarise

source("LME_calibration.R") # also calls dbpm_model_functions.R
source("Plotting_functions_DBPM.R")

LME_path <- "/rd/gem/private/fishmip_outputs/ISIMIP3a/DBPM/obsclim/"

# get the latest best values based on iterations and search_vol
no_iter = 100
search_vol = "estimated" 
version = "refined_"
vals <- data.frame(readRDS(paste0("Output/", version,"bestvals_LMEs_cor_searchvol_", search_vol,"_iter_",no_iter,".RDS")))

# run only LMEs with a good correlation/error and all catches 
LMEnumber<-which(vals$cor>0.5 & vals$rmse<0.5 & vals$catchNA==0) 

#### run gridded model by LME ----

#### MOVE TO FUNCTION CODE
rungridbyLME <- function(LMEnumber = 14, 
                         yearly = FALSE, 
                         f.effort = TRUE, 
                         vals = vals){
  
  # # CN trial
  # LMEnumber = 61
  # yearly = FALSE
  # f.effort = TRUE
  # vals = vals
  

  LME_path_full <- paste0(LME_path, "/rd/gem/private/fishmip_outputs/ISIMIP3a/DBPM/obsclim/","LME_",LMEnumber,"_output")
  if(!file.exists(LME_path_full)) dir.create(LME_path_full)
  
  ## Tests ----
  # test 1 - check no fishing runs and compare to the matching netcdf? 
  # lme_input_grids -> Yearly = T, gridded_params -> f.u = 0 and f.v = o 
  # test 2 - checking the fishing component.
  # lme_inout_grids -> Yearly = F, gridded_params -> f.u, f.v, f.minu, f.minv as per values below. 
  
  # get initial values from LME-scale results
  lme_input_init <-get_lme_inputs(LMEnumber = LMEnumber, gridded = F, yearly = F)
  
  # # vals <- readRDS("bestvals_LMEs.RDS")
  # vals <- readRDS("Output/bestvals_LMEs_iter_10.RDS")
  # 
  # # run model using time-averaged inputs
  # initial_results<-run_model(vals=c(0.1,0.5,1,1), # vals[LMEnumber,], # WARNING - check: replace with fake best values that work
  #                            input=lme_input_init,
  #                            withinput = F)
  
  ## OR ### CHECK VERSION ABOVE IS THE SAME # moved outside function
  # vals <- data.frame(readRDS("Output/bestvals_LMEs_iter_10.RDS")) 
  
  initial_results<-run_model(vals=unlist(vals[LMEnumber,]), 
                             input=lme_input_init,
                             withinput = F)
  
  U.initial<-rowMeans(initial_results$U[,240:1440])
  V.initial<-rowMeans(initial_results$V[,240:1440])
  W.initial<-mean(initial_results$W[240:1440])
  
  # plot to check initial values
  # plotsizespectrum(initial_results,params=initial_results$params,
  #                  itime=240:1440,
  #                  timeaveraged = TRUE)
  
  # get gridded inputs and run through all grid cells one timestep at a time
  
  ### corrected for gridded = F, now need to correct for gridded = T - CHECKED 
  # also need to check yearly = TRUE or comment it to avoid use - thus far not used. 
  lme_inputs_grid<- get_lme_inputs(LMEnumber = LMEnumber, 
                                   gridded = T, 
                                   yearly = yearly)[,c("lat","lon", "LME", "t","sst",
                                                       "sbt","er","intercept","slope",
                                                       "depth","NomActive_area_m2_relative" )]

  time<-unique(lme_inputs_grid$t)
  grid_results<-vector("list", length(time))
  
  # #### WARNING - more chcks here, er in particular
  # # CN: check that aggregated inputs match gridded inputs - similar...
  # # in one I calculated values using weihgted inputs across grid cells, in the other I caclaulted values for each grid cell
  # unique(lme_input_init$Year)
  # a<-lme_input_init %>%
  #   group_by(Year) %>%
  #   summarise(sst = mean(sst), 
  #             sbt = mean(sbt), 
  #             er = mean(er),
  #             intercept = mean(intercept),
  #             slope = mean(slope), 
  #             depth = mean(depth), 
  #             NomActive_area_m2 = mean(NomActive_area_m2, na.rm = T)) %>%
  #   filter(Year >= 1950)
  # 
  # ## WARNING!!!!  - er is double here! NEED TO CHECK - checked and not sure why.... 
  # b<-lme_inputs_grid %>%
  #   mutate(Year = year(t)) %>%
  #   group_by(Year) %>%
  #   summarise(sst = mean(sst), 
  #             sbt = mean(sbt), 
  #             er = mean(er),
  #             intercept = mean(intercept),
  #             slope = mean(slope), 
  #             depth = mean(depth), 
  #             NomActive_area_m2 = mean(NomActive_area_m2, na.rm = T)) %>%
  #   filter(Year >= 1950)
  # 
  
  ##### reorganise outputs storage for later
  ntime<-length(time)
  # ngrid<-dim(subset(lme_inputs_grid,t==time[1]))[1]
  params<-initial_results$params
  
  ###################### TEST GRIDDED MODEL
  
  lme_inputs_grid$cell <- paste(lme_inputs_grid$lat,lme_inputs_grid$lon,sep="_")
  
  depth_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = depth)
  
  er_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = er)
  
  intercept_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = intercept)
  
  slope_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = slope)
  
  sst_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = sst)
  
  sbt_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = sbt)
  
  effort_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = NomActive_area_m2_relative)
  
  ## Creating a 100 years stable spinup before 1841
  spinFunc <- function(var_name){
    dates <- seq(as.Date("1741-01-01"), as.Date("1840-12-01"), by = "month")
    var_pre <- array(NA, dim = c(dim(var_name)[1],length(dates)),
                     dimnames = list(dimnames(var_name)[[1]],as.character(dates)))
    var_pre[, 1:ncol(var_pre)] <- apply(var_name[,-1][,1:12], 1, mean)
    var_name <- cbind(var_name[,1],var_pre, var_name[,-1])
  }
  
  er_grid <- spinFunc(er_grid)
  intercept_grid <- spinFunc(intercept_grid)
  slope_grid <- spinFunc(slope_grid)
  sst_grid <- spinFunc(sst_grid)
  sbt_grid <- spinFunc(sbt_grid)
  effort_grid <- spinFunc(effort_grid)
  
  # adjusting cells to be 3240 values per cells (100 more years)(not used but 
  # could be useful)
  # cell <- lme_inputs_grid$cell
  # for(iCell in unique(cell)){
  #   pos <- which(cell == iCell)
  #   end <- pos[length(pos)]
  #   cell <- c(cell[1:end],rep(iCell,1200),cell[(end+1):length(cell)])
  # }
  # cell <- cell[1:(length(cell)-2)] #the code above produces a NA with the last cell
  
  if(f.effort){
    f.u<-as.numeric(vals[LMEnumber,1])
    f.v<-as.numeric(vals[LMEnumber,2])
    f.minu<-as.numeric(vals[LMEnumber,3])
    f.minv<-as.numeric(vals[LMEnumber,4])
  } else { 
    f.u <- f.v <- 0 
    f.minu <- f.minv <- 0 
    }
  
  search_vol<-as.numeric(vals[LMEnumber,5])

  # Making values constant through time
  # er_grid[,3:dim(er_grid)[2]] <- er_grid[,2]
  # intercept_grid[,3:dim(intercept_grid)[2]] <- intercept_grid[,2]
  # slope_grid[,3:dim(slope_grid)[2]] <- slope_grid[,2]
  # sst_grid[,3:dim(sst_grid)[2]] <- sst_grid[,2]
  # sbt_grid[,3:dim(sbt_grid)[2]] <- sbt_grid[,2]
  
  ## WARNING search_vol to adjust according to runLMEcalibration
  # needs to match runLMEcalibration (LME_calibration.R/run_model())
  
  # set up params for each month, across grid cells
  gridded_params <- sizeparam(equilibrium = FALSE
                              ,dx = 0.1
                              ,xmin.consumer.u = -3
                              ,xmin.consumer.v = -3
                              ,tmax = dim(er_grid[,-1])[2]/12
                              ,tstepspryr  =  12
                              ,search_vol = search_vol
                              ,fmort.u = f.u
                              ,fminx.u = f.minu
                              ,fmort.v = f.v
                              ,fminx.v = f.minv
                              ,depth = data.matrix(depth_grid[,-1][,1])
                              ,er = data.matrix(er_grid[,-1])
                              ,pp = data.matrix(intercept_grid[,-1])
                              ,slope = data.matrix(slope_grid[,-1])
                              ,sst = data.matrix(sst_grid[,-1])
                              ,sft = data.matrix(sbt_grid[,-1])
                              ,use.init = TRUE
                              ,effort = data.matrix(effort_grid[,-1])
                              ,U.initial = U.initial
                              ,V.initial = V.initial
                              ,W.initial = W.initial
                              ,Ngrid=dim(depth_grid)[1])      
  
  # run model  for full time period across all grid cells
  tic()
  grid_results<-gridded_sizemodel(gridded_params,
                                  ERSEM.det.input=F,
                                  temp.effect=T,
                                  eps=1e-5,
                                  output="aggregated",
                                  use.init = TRUE,
                                  burnin.len)
  toc()# 65.50608 min 
  
  
  #### Arrivata qui with LME 61 - works OK
  

  #### TEST 3 - OK working 
  # LME 1
  # random bestvalues - i.e. the ones that picked manually best approximate catches (0.1,0.5,1,1) 
  # search vol = 0.64 as OK for this LME
  # Fmort = first spread effort then calculate Fmort and catches as discussed with Julia 
  # gravity model option 2 with iter = 1 
  
  # removing the stable spinup section to have matching dimensions with the code
  # WARNING  move to plotting function for now as need to figure out catch trend
  grid_results$U <- grid_results$U[,,1201:3241]
  grid_results$V <- grid_results$V[,,1201:3241]
  grid_results$Y.u <- grid_results$Y.u[,,1201:3241]
  grid_results$Y.v <- grid_results$Y.v[,,1201:3241]
  # moved to plotting function as needed there
  # gridded_params$Neq <- 2040
  
  # save results from run
  saveRDS(grid_results,paste0(LME_path_full,"/grid_results.rds"))
  # save inputs and params object needed for plotting 
  save(lme_input_init, lme_inputs_grid, gridded_params, file =paste0(LME_path_full,"/grid_results_inputs_params.RData"))
  
}

## TEST 

# tic()
# rungridbyLME(LMEnumber = 1, 
#              yearly = FALSE, # for get_lme_inputs()
#              f.effort = TRUE, # for rungridbyLME()
#              vals = vals)
# toc() 
# 
# plotgridbyLME(LMEnumber = 1)

## RUN all LMEs 

# tic()
# pbapply::pbsapply(X=LMEnumber,rungridbyLME,yearly = FALSE, f.effort = TRUE, vals = vals, cl = detectCores() - 5)
# toc() # trying this now BUT: WARNING: Session forced to suspend due to system upgrade, restart, maintenance, or other issue. Your session data was saved however running computations may have been interrupted.
# 
# ## OR 
# tic()
# mclapply(LMEnumber, function(x) rungridbyLME(x, yearly = FALSE, f.effort = TRUE, vals = vals), mc.cores = detectCores()-5)
# toc() # 1536.929 sec elapsed/25 min. No working with all LMEs ?!

## OR standard loop (possibly very slow)

tic()
for (i in 1:length(LMEnumber)){
  
  # # trial 
  i = 42 # 41 is LME 61 considering the selection above
  print(paste0("Now working on LME", LMEnumber[i]))
  
  rungridbyLME(LMEnumber = LMEnumber[i],
               yearly = FALSE, # for get_lme_inputs()
               f.effort = TRUE, # for rungridbyLME()
               vals = vals)
  
}
toc()

# STACK at LME 61 

#### plot gridded and averadged output by LME ----

LMEnumber<-LMEnumber[1:40]

### now produce plots - EMPTY when run in // but OK when run individually and from inside the function. 
tic()
mclapply(LMEnumber, function(x) plotgridbyLME(x), mc.cores = detectCores()-5)
toc() 
# problem running this too... 

#### global map of gridded output ----

# this function is very similar to the beginning of plotgridbyLME() above but considers only gridded outputs
# wth the final aim of aggregating all outputs to plot global map 

#### MOVE TO FUNCTION CODE
##### POSSIBLY WORTH DOING ALL PLOTS CALCUALTIONS HERE AS THE MOST EXPENSIVE ACTION IS TO LOAD THE OUTPUTS - SO BETTER DOING IT ONLY ONCE ... 
# COMBINED WITH PLOTTING FUNCTION ABOVE... 
getGriddedOutputs_decade<-function(LME_path, LMEnumber){
  
  # # trial 
  # LMEnumber = 12
  
  # load outputs 
  LME_path_full = paste0(LME_path, "LME_",LMEnumber,"_output")
  full_file_name<-paste0(LME_path_full,"/grid_results.rds")
  full_file_name_output<-paste0(LME_path_full,"/grid_results_toCheckMap.rds")
  
  if(file.exists(full_file_name)){
    
    print(paste0("Now working on LME", LMEnumber))
    
    tic()
    grid_results <- readRDS(paste0(LME_path_full,"/grid_results.rds"))
    toc() # 49 sec
    
    # # ### WARNING need to comment out to figure out trend in biomass 
    # # removing the stable spinup section to have matching dimensions with the code
    # grid_results$U <- grid_results$U[,,1201:3241]
    # grid_results$V <- grid_results$V[,,1201:3241]
    # grid_results$Y.u <- grid_results$Y.u[,,1201:3241]
    # grid_results$Y.v <- grid_results$Y.v[,,1201:3241]
  
    # load inputs and param object 
    load(paste0(LME_path_full,"/grid_results_inputs_params.RData"))
    gridded_params$Neq <- 2040
  
    ### WARNING need to check depth integration and Neq - do they match what used for inputs and for runLMEcalibration? 
    tic()
    out<-getGriddedOutputs(input=lme_inputs_grid,results=grid_results,params=gridded_params)
    toc() # 18 sec
    
    # averaged across decades for U and V
    tic()
    out2<-out %>% 
      select(lat,lon,LME,t,TotalVcatch, TotalUcatch, TotalVbiomass, TotalUbiomass) %>% 
      mutate(year = year(t), 
             decade = floor(year/10)*10) %>% # original data resolution is month - can do year than decade. 
      group_by(LME, decade, lat, lon) %>% 
      summarize(
        Ubiomass_year = mean(TotalUbiomass,na.rm=T),
        Vbiomass_year = mean(TotalVbiomass,na.rm=T), 
        Ucatch_year = mean(TotalUcatch,na.rm=T),
        Vcatch_year = mean(TotalVcatch,na.rm=T)) %>% 
      ungroup() %>% 
      mutate(biomass_year = Ubiomass_year+Vbiomass_year,
             catch_year = Ucatch_year+Vcatch_year)
    toc() # 2 sec
    
    # write.csv(x = out2, file = file.path(full_file_name_output))
    fwrite(x = out2, file = file.path(full_file_name_output))
    
  }
}

# apply function in //
trial<-LMEnumber # for now run first LMEs that worked 
tic()
a<-mclapply(trial, function(x) getGriddedOutputs_decade(LME_path,LMEnumber = x), mc.cores = detectCores()-5)
toc() # 81 sec 6 LMEs.  

# get all files and put them together
file_list<-list.files(path = LME_path, pattern = "toCheckMap", recursive = TRUE, full.name = TRUE)
df_to_plot<-rbindlist(mclapply(file_list, function(x) fread(x), mc.cores = detectCores()-5))

### all LMEs togetehr and plot 

plot_global_raster<-function(df_to_plot, variable_to_plot, decade_to_plot){
  
  # # trial 
  # variable_to_plot = "biomass_year"
  # decade_to_plot = 2010
  
  raster_to_plot<-df_to_plot %>% 
    select(lat,lon,eval(variable_to_plot),decade) %>%  
    unique() %>% 
    filter(decade == decade_to_plot) %>% 
    select(lon,lat,eval(variable_to_plot)) %>% 
    relocate(lon,lat,eval(variable_to_plot)) %>% # order is key for raster function below   
    rasterFromXYZ()

  name = paste0(LME_path, "global_", variable_to_plot, decade_to_plot, "map", ".pdf")
  pdf(name)
  plot(raster_to_plot)
  dev.off()
  
  name = paste0("/home/dbpm/lme_scale_calibration_ISMIP3a/Output/", "global_", variable_to_plot, "_", decade_to_plot, "_map", ".pdf") # easier to check but needs to be deleted
  pdf(name)
  plot(raster_to_plot)
  dev.off()
  
}

plot_global_raster(df_to_plot = df_to_plot, variable_to_plot = "biomass_year", decade_to_plot = 2010)
plot_global_raster(df_to_plot = df_to_plot, variable_to_plot = "catch_year", decade_to_plot = 2010)

#### tesing ----

# #Test 1- compare with results from DBPM, no fishing (checking code consistency)
# # 1.	Test 1: run yearly = TRUE, no fishing (effort = 0), search volume = 64. 
# # b.	Compare these plots with the 3a tcb netcdf file: 
# #   i.	Extract LME 14 from this file 
# # ii.	Produce plots
# 
# tcb <- nc_open("/rd/gem/private/fishmip_outputs/ISIMIP3a/ctrlclim/netcdf/dbpm_ctrlclim_nobasd_1deg_nat_default_tcb_global_monthly_1961_2010.nc")
# #tcb <- nc_open("/rd/gem/private/fishmip_outputs/ISIMIP3a/ctrlclim/netcdf/dbpm_ctrlclim_nobasd_1deg_nat_default_tcblog10_global_monthly_1961_2010.nc")
# 
# 
# lat_mask <- tcb$dim$lat$vals >= min(lat_ext) & tcb$dim$lat$vals <= max(lat_ext)
# lon_mask <- tcb$dim$lon$vals >= min(lon_ext) & tcb$dim$lon$vals <= max(lon_ext)
# time_mask <- tcb$dim$time$vals >= 1841 & tcb$dim$time$vals <= 2010


# library(ncdf4)
# tcb <- nc_open("/rd/gem/private/fishmip_outputs/ISIMIP3a/ctrlclim/netcdf/dbpm_ctrlclim_nobasd_1deg_nat_default_tcb_global_monthly_1961_2010.nc")
# lat_mask <- tcb$dim$lat$vals >= min(lat_ext) & tcb$dim$lat$vals <= max(lat_ext)
# lon_mask <- tcb$dim$lon$vals >= min(lon_ext) & tcb$dim$lon$vals <= max(lon_ext)
# time_mask <- tcb$dim$time$vals >= 1841 & tcb$dim$time$vals <= 2010
# 
# ncdData <- ncvar_get(tcb, "tcb", 
#                      start = c(which(lon_mask)[1], which(lat_mask)[1], which(time_mask)[1]), 
#                      count = c(sum(lon_mask), sum(lat_mask), sum(time_mask)))
# dimnames(ncdData) <- list(tcb$var$tcb$dim[[1]]$vals[which(lon_mask)],
#                        tcb$var$tcb$dim[[2]]$vals[which(lat_mask)],
#                        tcb$var$tcb$dim[[3]]$vals[which(time_mask)])
# names(dimnames(ncdData)) <- c(tcb$var$tcb$dim[[1]]$name,
#                            tcb$var$tcb$dim[[2]]$name,
#                            tcb$var$tcb$dim[[3]]$name)
# 
# tcb_df <- reshape2::melt(ncdData)
# 
# tcb_decade_avg <- tcb_df %>%
#   mutate(decade = as.integer(substr(time, 1, 3)) * 10) %>% 
#   group_by(decade, lon, lat) %>%  
#   summarize(avg_biomass = mean(value))
# 
# # values from test1
# out <- readRDS("rungridRes/test1.rds")
# biom_df <- out[,c(1,2,4,16,17)]
# biom_df <- biom_df %>% mutate(totalB = TotalVbiomass + TotalUbiomass)
# biom_decade_avg <- biom_df %>%
#   mutate(decade = as.integer(substr(t, 1, 3)) * 10) %>% 
#   group_by(decade, lon, lat) %>%  
#   summarize(avg_biomass = mean(totalB))
# 
# # different cells so cannot just substract 
# 
# p1 <- ggplot(tcb_decade_avg)+
#   geom_tile(aes(x = lon, y = lat, fill = avg_biomass)) +
#   geom_sf(data = world) +
#   coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
#   scale_fill_gradient2(low = "white", high = "red", name = "Average Biomass in g/m2") +
#   facet_wrap(~decade,ncol = 6) +
#   scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], by = 6)) +
#   scale_y_continuous(name = "Latitude") +
#   theme(legend.position = "bottom") +
#   ggtitle("tcb netcdf")
# 
# p2 <- ggplot(biom_decade_avg)+
#   geom_tile(aes(x = lon, y = lat, fill = avg_biomass)) +
#   geom_sf(data = world) +
#   coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
#   scale_fill_gradient2(low = "white", high = "red", name = "Average Biomass in g/m2") +
#   facet_wrap(~decade,ncol = 6) +
#   scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], by = 6)) +
#   scale_y_continuous(name = "Latitude") +
#   theme(legend.position = "bottom") +
#   ggtitle("test 1 runGrid")
# 
# pdf("test1res.pdf", onefile = T)
# p1
# p2
# dev.off()
# 
# #Test 2- model calibration/behaviour
# # 1. Spatial maps of TotalCatch and TotalBiomass (by decade)
# # 2. Community size spectra (U & V) - one line per grid cell - final decade?
# # 3. Plots of GG growth rates (see historyMatching.R for example)
# # 4. Compare time series to total catches at LME scale with obs catch
# # 5. Once checked, run for all LMEs
# 
# # TODO cami: compare model with empirical catches (see run lme calibration for plotting)
# 
