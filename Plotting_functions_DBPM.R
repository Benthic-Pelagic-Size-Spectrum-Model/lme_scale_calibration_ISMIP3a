
# CN 20/08/2023 - create a plotting file for DBPM based on plots prepared by Romain Forestier and JB


plotgridbyLME<-function(LMEnumber = 1){
  
  # # trial 
  # LMEnumber = 2
  
  # load outputs 
  LME_path = paste0("/rd/gem/private/fishmip_outputs/ISIMIP3a/DBPM/obsclim/LME_",LMEnumber,"_output")
  grid_results <- readRDS(paste0(LME_path,"/grid_results.rds"))
  
  # # ### WARNING need to comment out to figure out trend in biomass 
  # # removing the stable spinup section to have matching dimensions with the code
  # grid_results$U <- grid_results$U[,,1201:3241]
  # grid_results$V <- grid_results$V[,,1201:3241]
  # grid_results$Y.u <- grid_results$Y.u[,,1201:3241]
  # grid_results$Y.v <- grid_results$Y.v[,,1201:3241]

  # load inputs and param object 
  load(paste0(LME_path,"/grid_results_inputs_params.RData"))
  # lme_inputs_grid = lme_inputs_grid
  # gridded_params = gridded_params
  gridded_params$Neq <- 2040

  if(sum(gridded_params$effort)>0){f.effort = TRUE}else{f.effort = FALSE}
  
  ### WARNING need to check depth integration and Neq - do they match what used for inputs and for runLMEcalibration? 
  out<-getGriddedOutputs(input=lme_inputs_grid,results=grid_results,params=gridded_params)

  ## CN- is this necessary as you already saved gridded_output and have out here?? 
  # saveRDS(out,paste0(LME_path,"/out_results.rds"))

  cells<-unique(out$cell)
  out$cell<-as.factor(out$cell)

  # colnames(out)
  head(out[,c("t","Totalcatch","TotalVcatch","TotalUcatch","TotalVbiomass","TotalUbiomass")])

  # download world map
  # world <- ne_download(category = "cultural", 
  #                      type = "admin_0_countries", 
  #                      scale = "large",
  #                      returnclass = "sf")
  # 
  # saveRDS(world,"worldMap.rds")
  world <- readRDS("Tmp_data/worldMap.rds")

  ####### Maps of biomass ----
  # averaged across decades for U and V
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
  
  ## CN plot too big - selelct 
  biom_decade_avg<-biom_decade_avg %>% 
    filter(decade %in% c(1950, 2010))
  
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

  
  # WARNING - CN version - something wrong in trends so calculating this too - to finish   
  biom_decade_avg_CN<- biom_df %>%
    mutate(year = year(t)) %>% 
    group_by(year, lat, lon) %>% 
    summarize(avg_Ubiomass = mean(TotalUbiomass,na.rm=T),
              avg_Vbiomass = mean(TotalVbiomass,na.rm=T)) %>% 
    ungroup()
  
  ####### Biomass trends ---- 
  # calculate the mean biomass across gridcell
  biom_grid_avg <- biom_df %>%
    mutate(Year = year(t)) %>% 
    group_by(Year) %>%  
    summarize(avg_biomass = mean(totalB)) %>% 
    ungroup()

  # time series of biomass
  p3 <- ggplot(biom_grid_avg, aes(x = Year, y = avg_biomass))+
    geom_line() +
    scale_y_continuous(name = "Average biomass in g/m2")+
    scale_x_continuous(name = "Time in year") +
    ggtitle("Average biomass through time")

  ####### Size spectra ---- 
  # U + V averaged per decade per longitude
  # Using grid_results$U and V. Dim are gridcell*size*time (monthly)
  
  totBiom <- grid_results$U + grid_results$V
  # averaging per decade. First decade is 108 month then the rest is 120

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

  ####### Growth rate ---- 
  # GG.u + GG.v per decade per gridcell
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

  ####### Maps of catches ---- 
  catch_df <- out[,c(1,2,4,14,15)]
  catch_df <- catch_df %>% mutate(totalC = TotalVcatch + TotalUcatch)

  # calculate the mean catch for each decade
  catch_decade_avg <- catch_df %>%
    mutate(decade = as.integer(substr(t, 1, 3)) * 10) %>%
    group_by(decade, lon, lat) %>%
    summarize(avg_Ucatch = mean(TotalUcatch),
              avg_Vcatch = mean(TotalVcatch))

  ## CN plot too big - select 
  catch_decade_avg<-catch_decade_avg %>% 
    filter(decade %in% c(1950, 2010))
  
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

  ####### Catch trends ---- 
  # U + V with empirical data
  # calculate total catch across gridcells
  catch_grid_avg <- catch_df %>%
    mutate(Year = year(t)) %>% 
    filter(Year >=1950) %>% 
    group_by(Year) %>%  
    summarize(avg_catch = mean(totalC)) #WARNING weighted.mean here would be better 
  # empirical data is in lme_input_init$catch_tonnes_area_m2,it's monthly, convert to yearly

  catch_empirical<-lme_input_init %>% 
    mutate(Year = year(t)) %>% 
    group_by(Year) %>%  
    filter(Year >=1950) %>%
    summarize(empirical = mean(catch_tonnes_area_m2)) %>% 
    mutate(empirical = empirical*1e06) # CN WARNING this needs to be transformed in g m2 

  catch_grid_avg<-catch_grid_avg %>% 
    full_join(catch_empirical)

  p8 <- ggplot(catch_grid_avg) +
    geom_line(aes(x = Year, y = avg_catch))+
    geom_point(aes(x = Year, y = empirical)) +
    ggtitle("Time series of catches versus empirical data (U+V)")

  ####### Save plots ---- 

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
  # name = ifelse(f.effort == FALSE, "_no_fishing.pdf", "_fishing.pdf")
  # ## adding gravity option to plot name for testing 
  # pdf(paste0("Output/Gravity2_LME_",LMEnumber, name), height = 8, width = 6, onefile = T)
  # print(p1)
  # print(p2)
  # print(p3)
  # print(p4)
  # print(p5)
  # print(p6)
  # print(p7)
  # print(p8)
  # dev.off()
  
  # # increase the size of some plots as are not visible 
  # pdf(paste0("Output/Gravity2_BiomassU_LME_",LMEnumber, name))
  # p1
  # dev.off()
  # 
  # pdf(paste0("Output/Gravity2_BiomassV_LME_",LMEnumber, name))
  # p2
  # dev.off()
  # 
  # pdf(paste0("Output/Gravity2_CatchU_LME_",LMEnumber, name))
  # p6
  # dev.off()
  # 
  # pdf(paste0("Output/Gravity2_CatchV_LME_",LMEnumber, name))
  # p7
  # dev.off()
  
}

