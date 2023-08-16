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

rm(list=ls())

rungridbyLME <- function(LMEnumber = 14, 
                         yearly = FALSE, 
                         f.effort = TRUE, 
                         # search_vol = 0.64,
                         savePlots = F){
  
  # # CN trial
  # LMEnumber = 14
  # yearly = FALSE
  # f.effort = TRUE
  # # search_vol = 0.64
  # savePlots = TRUE

  # setup
  source("LME_calibration.R")
  
  LME_path <- paste0("/rd/gem/private/fishmip_outputs/ISIMIP3a/DBPM/obsclim/LME_",LMEnumber,"_output")
  if(!file.exists(LME_path)) dir.create(LME_path)
  
  ## Tests ----
  # test 1 - check no fishing runs and compare to the matching netcdf? 
  # lme_input_grids -> Yearly = T, gridded_params -> f.u = 0 and f.v = o 
  # test 2 - checking the fishing component.
  # lme_inout_grids -> Yearly = F, gridded_params -> f.u, f.v, f.minu, f.minv as per values below. 
  
  # get initial values from LME-scale results
  lme_input_init <-get_lme_inputs(LMEnumber = LMEnumber, gridded = F, yearly = F)
  
  ### WARNING - update with latest bestvalue dataset!
  # vals <- readRDS("bestvals_LMEs.RDS")
  vals <- readRDS("Output/bestvals_LMEs_iter_10.RDS")
  
  ## CN WARNING - this gives NAs for all size classes >90 - ask Julia if this is OK
  ## CN CORRECTION: added a search_vol param in function argument because
  ## search_vol needs to be consistent with what we are using to run the final 
  ## model and was not (it was set to 64 inside the run_model() function) 
  # run model using time-averaged inputs
  initial_results<-run_model(vals=vals[LMEnumber,],
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
  
  ### WARNING - CN error here - CORRECTED due to wrong input resolution in function 
  ### WARNING LME.x and LME.y to correct - CORRECTED
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
  
  ## NOTE - here is where you change fishing effort 
  effort_grid<-lme_inputs_grid %>%
    pivot_wider(id_cols=cell,names_from = t, values_from = NomActive_area_m2_relative)
  
  # # check effort # very low 
  # effort_grid[1:10, 1:10]
  
  
  
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
    # CN correct 
    f.u<-as.numeric(vals[LMEnumber,1])
    f.v<-as.numeric(vals[LMEnumber,2])
    # f.u<-as.numeric(vals[1])
    # f.v<-as.numeric(vals[2])
  } else  f.u <- f.v <- 0
  
  # CN corrected by adding LMEnumber 
  f.minu<-as.numeric(vals[LMEnumber,3])
  f.minv<-as.numeric(vals[LMEnumber,4])
  
  # Making values constant through time
  # er_grid[,3:dim(er_grid)[2]] <- er_grid[,2]
  # intercept_grid[,3:dim(intercept_grid)[2]] <- intercept_grid[,2]
  # slope_grid[,3:dim(slope_grid)[2]] <- slope_grid[,2]
  # sst_grid[,3:dim(sst_grid)[2]] <- sst_grid[,2]
  # sbt_grid[,3:dim(sbt_grid)[2]] <- sbt_grid[,2]
  
  ## WARNING search_vol to adjust according to runLMEcalibration
  # now using 0.064 as in runLMEcalibration (LME_calibration.R/run_model())
  
  # set up params for each month, across grid cells
  gridded_params <- sizeparam(equilibrium = FALSE
                              ,dx = 0.1
                              ,xmin.consumer.u = -3
                              ,xmin.consumer.v = -3
                              ,tmax = dim(er_grid[,-1])[2]/12
                              ,tstepspryr  =  12
                              ,search_vol = 0.064
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
  grid_results<-gridded_sizemodel(gridded_params,
                                  ERSEM.det.input=F,
                                  # U_mat,
                                  # V_mat,
                                  # W_mat,
                                  temp.effect=T,
                                  eps=1e-5,
                                  output="aggregated",
                                  use.init = TRUE,
                                  burnin.len)
  
  # # Checks
  # U <- grid_results$U
  # dim(U)
  # U[1,,2040]
  # sum(is.na(U))
  # sum(any(U < 0))
  # dim(U) # this grid cell X size X time 
  # sum(U, na.rm = TRUE)
  # is.na(U[])
  # 
  # V <- grid_results$V
  # V[1,,2041]
  # sum(is.na(V))
  # sum(any(V < 0))
  
  saveRDS(grid_results,paste0(LME_path,"/grid_results.rds"))
  # grid_results <- readRDS(paste0(LME_path,"/grid_results.rds"))
  
  #removing the stable spinup section to have matching dimensions with the code
  grid_results$U <- grid_results$U[,,1201:3241]
  grid_results$V <- grid_results$V[,,1201:3241]
  grid_results$Y.u <- grid_results$Y.u[,,1201:3241]
  grid_results$Y.v <- grid_results$Y.v[,,1201:3241]
  gridded_params$Neq <- 2040
  
  
  
  
  
  
  
  

  
  ### WARNING need to check depth integration and Neq - do they match what used for inputs and for runLMEcalibration? 
  out<-getGriddedOutputs(input=lme_inputs_grid,results=grid_results,params=gridded_params)
  saveRDS(out,paste0(LME_path,"/out_results.rds"))
  # out <- readRDS(paste0(LME_path,"/out_results.rds"))
  
  ## WARNING - Save results here and move plots to another function?? 
  #### CHECK OUTPUTS!!
  
  cells<-unique(out$cell)
  out$cell<-as.factor(out$cell)
  
  colnames(out)
  head(out[,c("t","Totalcatch","TotalVcatch","TotalUcatch","TotalVbiomass","TotalUbiomass")])
  
  # download world map
  # world <- ne_download(category = "cultural", 
  #                      type = "admin_0_countries", 
  #                      scale = "large",
  #                      returnclass = "sf")
  # 
  # saveRDS(world,"worldMap.rds")
  world <- readRDS("worldMap.rds")
  
  # Plots ----
  ## Maps of biomass averaged across decades for U and V
  biom_df <- out[,c(1,2,4,16,17)]
  biom_df <- biom_df %>% mutate(totalB = TotalVbiomass + TotalUbiomass)
  
  # calculate the mean biomass for each decade
  biom_decade_avg <- biom_df %>%
    mutate(decade = as.integer(substr(t, 1, 3)) * 10) %>% 
    group_by(decade, lon, lat) %>%  
    summarize(avg_Ubiomass = mean(TotalUbiomass,na.rm=T),
              avg_Vbiomass = mean(TotalVbiomass,na.rm=T))
  
  # Axis
  lat_ext <- unique(biom_decade_avg)$lat[c(1,dim(biom_decade_avg)[1])]
  lon_ext <- unique(biom_decade_avg)$lon[c(1,dim(biom_decade_avg)[1])]
  
  # facet plot U
  p1 <- ggplot(biom_decade_avg)+
    geom_tile(aes(x = lon, y = lat, fill = avg_Ubiomass)) +
    geom_sf(data = world) +
    coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
    scale_fill_gradient2(low = "white", high = "red", name = "Average Biomass in g/m2") +
    facet_wrap(~decade,ncol = 6) +
    scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], by = 6)) +
    scale_y_continuous(name = "Latitude") +
    theme(legend.position = "bottom") +
    ggtitle("Maps of biomass averaged across decades for U")
  
  # facet plot V
  p2 <- ggplot(biom_decade_avg)+
    geom_tile(aes(x = lon, y = lat, fill = avg_Vbiomass)) +
    geom_sf(data = world) +
    coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
    scale_fill_gradient2(low = "white", high = "red", name = "Average Biomass in g/m2") +
    facet_wrap(~decade,ncol = 6) +
    scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], by = 6)) +
    scale_y_continuous(name = "Latitude") +
    theme(legend.position = "bottom") +
    ggtitle("Maps of biomass averaged across decades for V")
  
  ## Time series of biomass 1841-2020 U + V 
  # calculate the mean biomass across gridcell
  biom_grid_avg <- biom_df %>%
    group_by(t) %>%  
    summarize(avg_biomass = mean(totalB))  
  
  # time series of biomass
  p3 <- ggplot(biom_grid_avg)+
    geom_line(aes(x = t, y = avg_biomass)) +
    scale_y_continuous(name = "Average biomass in g/m2")+
    scale_x_date(name = "Time in year") +
    ggtitle("Average biomass through time")
  
  ## Size spectra U + V averaged per decade per longitude
  # Using grid_results$U and V. Dim are gridcell*size*time (monthly)
  
  totBiom <- grid_results$U + grid_results$V
  # averaging per decade. First decade is 108 month then the rest is 120
  
  
  
  
  
  
  
  
  ### CN WARNING should this be 2040 now that we saved up to 2040? BELOW TOO
  
  
  decade_start <- c(1,seq(109,2041,by = 120))
  spectra_decade_avg <- array(NA, dim = c(dim(totBiom)[1:2],length(decade_start)),
                              dimnames = list("gridCell" = 1:dim(totBiom)[1],
                                              "size" = grid_results$params$x,
                                              "decade" = seq(1840,2010,by = 10)))
  for(iTime in 1:(length(decade_start)-1)){
    t_start <- decade_start[iTime]
    t_end <- decade_start[iTime+1]
    tempBiom <- totBiom[,,t_start:t_end]
    avgBiom <- apply(tempBiom,c(1,2),mean)
    spectra_decade_avg[,,iTime] <- avgBiom
  }
  
  
  # average per longitude
  biom_df$cell <- paste(biom_df$lat,biom_df$lon)
  cell <- unique(biom_df$cell) # each grid cell is a combo of lat and long
  cell_lat <- substr(cell, 1,5)
  lat_list <- vector("list", length = length(unique(cell_lat)))
  names(lat_list) <- unique(cell_lat)
  for(iCell in unique(cell_lat)) lat_list[[iCell]] <- which(cell_lat == iCell)
  # this vector contains the id of grid cells having the same latitude
  
  spectra_grid_avg <- array(NA, dim = c(length(lat_list),dim(spectra_decade_avg)[2:3]),
                            dimnames = list("latitude" = names(lat_list),
                                            "size" = dimnames(spectra_decade_avg)[[2]],
                                            "decade" = dimnames(spectra_decade_avg)[[3]]))
  for (iCell in names(lat_list)) {
    tempBiom <- spectra_decade_avg[lat_list[[iCell]],,]
    if(length(dim(tempBiom)) == 3) avgBiom <- apply(tempBiom,c(2,3),mean) # if = 2 means only one lat value
    else avgBiom <- tempBiom
    spectra_grid_avg[iCell,,] <- avgBiom
  }
  
  # show sizes only between $ref and #Nx
  spectra_grid_avg <- spectra_grid_avg[,grid_results$params$ref:grid_results$params$Nx,]
  
  plot_dat <- reshape2::melt(spectra_grid_avg)
  
  p4 <- ggplot(plot_dat) +
    geom_line(aes(x = size, y = value, color = latitude), alpha = .5) +
    facet_wrap(~decade) +
    scale_y_continuous(trans = "log10", name = "Biomass in g") +
    scale_x_continuous(name = "Size in log10 g") +
    ggtitle("Size spectra averaged per decade per longitude (U+V)")
  
  ## plot growth rate GG.u + GG.v per decade per gridcell
  totGrowth <- grid_results$GG.u + grid_results$GG.v
  
  # averaging per decade. First decade is 108 month then the rest is 120
  decade_start <- c(1,seq(109,2041,by = 120))
  growth_decade_avg <- array(NA, dim = c(dim(totBiom)[1:2],length(decade_start)),
                             dimnames = list("gridCell" = 1:dim(totBiom)[1],
                                             "size" = grid_results$params$x,
                                             "decade" = seq(1840,2010,by = 10)))
  for(iTime in 1:(length(decade_start)-1)){
    t_start <- decade_start[iTime]
    t_end <- decade_start[iTime+1]
    tempGrowth <- totGrowth[,,t_start:t_end]
    avgGrowth <- apply(tempGrowth,c(1,2),mean)
    growth_decade_avg[,,iTime] <- avgGrowth
  }
  
  # average per longitude - using lat_list from previous plot
  growth_grid_avg <- array(NA, dim = c(length(lat_list),dim(growth_decade_avg)[2:3]),
                           dimnames = list("latitude" = names(lat_list),
                                           "size" = dimnames(growth_decade_avg)[[2]],
                                           "decade" = dimnames(growth_decade_avg)[[3]]))
  
  for (iCell in names(lat_list)) {
    tempBiom <- growth_decade_avg[lat_list[[iCell]],,]
    if(length(dim(tempBiom)) == 3) avgBiom <- apply(tempBiom,c(2,3),mean) # if = 2 means only one lat value
    else avgBiom <- tempBiom
    growth_grid_avg[iCell,,] <- avgBiom
  }
  
  # show sizes only between $ref and #Nx
  growth_grid_avg <- growth_grid_avg[,grid_results$params$ref:grid_results$params$Nx,]
  
  plot_dat <- reshape2::melt(growth_grid_avg)
  
  p5 <- ggplot(plot_dat) +
    geom_line(aes(x = size, y = value, color = latitude)) +
    facet_wrap(~decade) +
    scale_y_continuous(trans = "log10", name = "Relative growth rate per year") +
    scale_x_continuous(name = "Size in log10 g") +
    ggtitle("Growth rate averaged per decade per longitude (U+V)")
  
  ## plot total catch U and V from out
  catch_df <- out[,c(1,2,4,14,15)]
  catch_df <- catch_df %>% mutate(totalC = TotalVcatch + TotalUcatch)
  
  # calculate the mean catch for each decade
  catch_decade_avg <- catch_df %>%
    mutate(decade = as.integer(substr(t, 1, 3)) * 10) %>% 
    group_by(decade, lon, lat) %>%  
    summarize(avg_Ucatch = mean(TotalUcatch),
              avg_Vcatch = mean(TotalVcatch))  
  
  # plot map facets of average U catch per decade
  p6 <- ggplot(catch_decade_avg)+
    geom_tile(aes(x = lon, y = lat, fill = avg_Ucatch)) +
    geom_sf(data = world) +
    coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
    scale_fill_gradient2(low = "white", high = "red", name = "Average Catch in X") +
    facet_wrap(~decade,ncol = 6) +
    scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], by = 6)) +
    scale_y_continuous(name = "Latitude") +
    theme(legend.position = "bottom") +
    ggtitle("Maps of catches averaged across decades for U")
  
  # plot map facets of average V catch per decade
  p7 <- ggplot(catch_decade_avg)+
    geom_tile(aes(x = lon, y = lat, fill = avg_Vcatch)) +
    geom_sf(data = world) +
    coord_sf(xlim = lon_ext, ylim = lat_ext, expand = FALSE) +
    scale_fill_gradient2(low = "white", high = "red", name = "Average Catch in X") +
    facet_wrap(~decade,ncol = 6) +
    scale_x_continuous(name = "Longitude", breaks = seq(lon_ext[1],lon_ext[2], by = 6)) +
    scale_y_continuous(name = "Latitude") +
    theme(legend.position = "bottom") +
    ggtitle("Maps of catches averaged across decades for V")
  
  ## plot time series of catches U  + V with empirical data
  # calculate total catch across gridcells
  catch_grid_avg <- catch_df %>%
    group_by(t) %>%  
    summarize(avg_catch = mean(totalC))  
  # empirical data is in lme_input_init$catch_tonnes_area_m2,it's monthly, convert to yearly
  catch_grid_avg$empirical <- lme_input_init$catch_tonnes_area_m2
  
  p8 <- ggplot(catch_grid_avg) +
    geom_line(aes(x = t, y = avg_catch))+
    geom_point(aes(x = t, y = empirical)) +
    ggtitle("Time series of catches versus empirical data (U+V)")
  
  # Save the plots in a PDF file
  if(savePlots){
    
    # # CN change 
    # LME_plot<-list(p1,p2,p3,p4,p5,p6,p7,p8)
    # # pdf(paste0(LME_path,"/plots.pdf"), height = 8, width = 6)
    # tic()
    # pdf("plots.pdf", height = 8, width = 6)
    # marrangeGrob(grobs = LME_plot, nrow=2, ncol=1)
    # dev.off()
    # toc()
    
    name = ifelse(f.effort == FALSE, "/plots_no_fishing.pdf", "/plots_fishing.pdf")
    
    pdf(paste0(LME_path, name), height = 8, width = 6, onefile = T)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    dev.off()
    
    # print also in local for fast checking 
    name = ifelse(f.effort == FALSE, "_no_fishing.pdf", "_fishing.pdf")
    pdf(paste0("Output/LME_",LMEnumber, name), height = 8, width = 6, onefile = T)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    dev.off()
    
    
  }
}

## TESTS:
# test run for LME14 AND a different LME
library(tictoc)
tic()
rungridbyLME(LMEnumber = 14, 
             yearly = FALSE, # for get_lme_inputs()
             f.effort = FALSE, # for rungridbyLME()
             # search_vol = 0.64,# for sizeparam() but indicated as value there now - can change. 
             savePlots = TRUE)
toc() # 43.18363 min for LME 14 

tic()
rungridbyLME(LMEnumber = 14, 
             yearly = FALSE, # for get_lme_inputs()
             f.effort = TRUE, # for rungridbyLME()
             # search_vol = 0.64,# for sizeparam() but indicated as value there now - can change. 
             savePlots = TRUE)
toc() # 43.18363 min for LME 14 




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
