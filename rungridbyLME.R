##### run model across space and time using gridded inputs for each LME

source("LME_calibration.R")

# get initial values from LME-scale results
lme_input_14<-get_lme_inputs(LMEnumber = 14,gridded = F, yearly = F)
vals <- readRDS("bestvals_LMEs.RDS")
# run model using time-averaged inputs
initial_results<-run_model(vals=vals[14,],input=lme_input_14,withinput = F)
U.initial<-rowMeans(initial_results$U[,240:1440])
V.initial<-rowMeans(initial_results$V[,240:1440])
W.initial<-mean(initial_results$W[240:1440])
# plot to check initial values
plotsizespectrum(initial_results,params=initial_results$params,
                 itime=240:1440,
                 timeaveraged = TRUE)

# WARNING 
# test 1 - check no fishing runs, are they matching the netcdf? Yearly = T, f.u = 0 and f.v = o 
# test 2 - checking the fishing component. yearly = F, f.u, f.v, f.minu, f.minv as per values below.     

# get gridded inputs and run through all grid cells one timestep at a time
lme_inputs_grid<-
  get_lme_inputs(LMEnumber = 14, gridded = T, yearly = T)[,c("lat","lon", "LME.x", "t","sst",
                                                             "sbt","er","intercept","slope",
                                                             "depth","NomActive_area_m2" )]
time<-unique(lme_inputs_grid$t)
grid_results<-vector("list", length(time))

##### reorganise outputs storage for later
ntime<-length(time)
ngrid<-dim(subset(lme_inputs_grid,t==time[1]))[1]
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
  pivot_wider(id_cols=cell,names_from = t, values_from = NomActive_area_m2)

f.u<-as.numeric(vals[1])
f.v<-as.numeric(vals[2])
f.minu<-as.numeric(vals[3])
f.minv<-as.numeric(vals[4])

# Making values constant through time

# er_grid[,3:dim(er_grid)[2]] <- er_grid[,2]
# intercept_grid[,3:dim(intercept_grid)[2]] <- intercept_grid[,2]
# slope_grid[,3:dim(slope_grid)[2]] <- slope_grid[,2]
# sst_grid[,3:dim(sst_grid)[2]] <- sst_grid[,2]
# sbt_grid[,3:dim(sbt_grid)[2]] <- sbt_grid[,2]

# set up params for each month, across grid cells
gridded_params <- sizeparam (equilibrium = FALSE
                             ,dx = 0.1
                             ,xmin.consumer.u = -3
                             ,xmin.consumer.v = -3
                             ,tmax = dim(er_grid[,-1])[2]/12
                             ,tstepspryr  =  12
                             ,search_vol = .64 # just for this test instead of .64
                             ,fmort.u = 0#f.u
                             ,fminx.u = f.minu
                             ,fmort.v = 0#f.v
                             ,fminx.v = f.minv
                             ,depth = data.matrix(depth_grid[,-1][,1])
                             ,er = data.matrix(er_grid[,-1])
                             ,pp = data.matrix(intercept_grid[,-1])
                             ,slope = data.matrix(slope_grid[,-1])
                             ,sst = data.matrix(sst_grid[,-1])
                             ,sft = data.matrix(sbt_grid[,-1])
                             ,use.init = TRUE
                             ,effort = data.matrix(effort_grid[,-1])
                             ,U.initial =U.initial
                             ,V.initial = V.initial
                             ,W.initial = W.initial
                             ,Ngrid=dim(depth_grid)[1])      

# run model  for full time period across all grid cells
grid_results<-gridded_sizemodel(gridded_params,
                                ERSEM.det.input=F,
                                U_mat,V_mat,
                                W_mat,
                                temp.effect=T,
                                eps=1e-5,
                                output="aggregated",
                                use.init = TRUE,
                                burnin.len)

# Checks
# U <- grid_results$U
# U[1,,2041]
# sum(is.na(U))
# sum(any(U < 0))
# 
# V <- grid_results$V
# V[1,,2041]
# sum(is.na(V))
# sum(any(V < 0))
saveRDS(grid_results,"rungridRes/test1_grid.rds")

out<-getGriddedOutputs(input=lme_inputs_grid,results=grid_results,params=gridded_params)
# save gridded object instead

saveRDS(out,"rungridRes/test2.rds")
out <- readRDS("rungridRes/test2.rds")
#### CHECK OUTPUTS!!

cells<-unique(out$cell)
out$cell<-as.factor(out$cell)

# ggplot(filter(out,cell==cells[1]), aes(x=t,y=TotalUbiomass)) + geom_line()
# ggplot(filter(out,cell==cells[1]), aes(x=t,y=TotalVbiomass)) + geom_line()
# ggplot(filter(out,cell==cells[1]), aes(x=t,y=Totalcatch)) + geom_line()

p1<-ggplot(out, aes(x=t,y=log10(TotalUbiomass),group=cell)) + 
  geom_line(aes(color=sst))+theme(legend.position = "none") + 
  scale_color_continuous()
p2<-ggplot(out, aes(x=t,y=log10(TotalVbiomass),group=cell)) + 
  geom_line(aes(color=sst))+theme(legend.position = "none")+ 
  scale_color_continuous()
p3<-ggplot(out, aes(x=t,y=log10(Totalcatch),group=cell)) + 
  geom_line(aes(color=sst))+theme(legend.position = "right")+ 
  scale_color_continuous()

p1 + p2 + p3


# save results
#saveRDS(grid_results,"grid_results.RDS")

# 1.	Test 1: run yearly = TRUE, no fishing (effort = 0), search volume = 64. 
# a.	Plot:
#   i.	maps of biomass for last decade (average across 1961-1970 etc … up to 1991-2010) for U + V, 
# ii.	time series of biomass 1841-2020 U + V (U consider only size greater than a certain size, 
# check 3a netcdf code but this should be included into the getGriddedOutputs) – mean across grid cells, 
# iii.	size spectra U + V average per decade per grid cell and all grid cells on the same plot.  
# iv. plot growth rate GG.u + GG.v per decade per gridcell
# b.	Compare these plots with the 3a tcb netcdf file: 
#   i.	Extract LME 14 from this file 
# ii.	Produce plots i-iii  
# TODO cami: compare model with empirical catches (see run lme calibration for plotting)
# use history matching and plotsizespectrum() for help

# setup data
library(tidyverse)
library(rnaturalearth)
library(sf)

biom_df <- out[,c(1,2,4,16,17)]
biom_df <- biom_df %>% mutate(totalB = TotalVbiomass + TotalUbiomass)

# calculate the mean biomass for each decade
df_decade_avg <- biom_df %>%
  mutate(decade = as.integer(substr(t, 1, 3)) * 10) %>% 
  group_by(decade, lon, lat) %>%  
  summarize(avg_biomass = mean(totalB))  

# download world map
world <- ne_download(category = "cultural", 
                     type = "admin_0_countries", 
                     scale = "large",
                     returnclass = "sf")

# plot map facets of average biomass per decade
p1 <- ggplot(df_decade_avg)+
  geom_tile(aes(x = lon, y = lat, fill = avg_biomass)) +
  geom_sf(data = world) +
  coord_sf(xlim = c(-68.5,-52.5), ylim = c(-54.5,-34.5), expand = FALSE) +
  scale_fill_gradient2(low = "white", high = "red", name = "Avg Biomass") +
  facet_wrap(~decade) +
  theme(legend.position = "bottom") +
  ggtitle("Map of average biomass per decade")

# calculate the mean biomass across gridcell
df_grid_avg <- biom_df %>%
  group_by(t) %>%  
  summarize(avg_biomass = mean(totalB))  

# plot of time series of biomass

p2 <- ggplot(df_grid_avg)+
  geom_line(aes(x = t, y = avg_biomass)) +
  ggtitle("Average biomass through time")

# size spectra
spectra_df <- out[,c(1,2,4,8,9)]

spectra_decade_avg <- spectra_df %>%
  mutate(decade = as.integer(substr(t, 1, 3)) * 10) %>% 
  group_by(decade, lon, lat) %>%  
  summarize(avg_int = mean(intercept), avg_slope = mean(slope)) 

# TODO average per longitude too using $cell

p3 <- ggplot(spectra_decade_avg) +
  geom_abline(aes(slope = avg_slope, intercept =  avg_int, color = lat)) +
  facet_wrap(~decade)

# Using grid_results from nom, assuming that arrays dims are gridcell*size*time
# dim is 132 181 2041
# time is monthly
# show sizes only between $ref and #Nx

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

plot_dat <- reshape2::melt(spectra_decade_avg)

p4 <- ggplot(plot_dat) +
  geom_line(aes(x = size, y = value, color = gridCell)) +
  facet_wrap(~decade) +
  scale_y_continuous(trans = "log10")

# growth

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

plot_dat <- reshape2::melt(growth_decade_avg)

p5 <- ggplot(plot_dat) +
  geom_line(aes(x = size, y = value, color = gridCell)) +
  facet_wrap(~decade) +
  scale_y_continuous(trans = "log10")

p <- list(p1,p2,p3,p4,p5)
# saving
library(gridExtra)

# Save the plots in a PDF file
pdf("rungridRes/grid_results_test1.pdf", onefile = T)
p1
p2
p3
p4
p5
dev.off()

# add p6 map of catch using TotalCatch


# TO DO: 
#Make a folder for each LME  output with the output files 
# and the following figures in PDFs:

#Test 1- compare with results from DBPM, no fishing (checking code consistency)

#Test 2- model calibration/behaviour
# 1. Spatial maps of TotalCatch and TotalBiomass (by decade)
# - use Camis raster code
# 2. Community size spectra (U & V) - one line per grid cell - final decade?
# 3. Plots of GG growth rates (see historyMatching.R for example)
# 4. Compare time series to total catches at LME scale with obs catch
# 5. Once checked, run for all LMEs

# TODO make it a function to automatically save a pdf of plots per LME