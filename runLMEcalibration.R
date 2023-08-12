###### Run LHS search for each lme to get an initial rough guess for fihsing parameters
# lmenum=66
# LMEs<-data.frame(LMEnum=1:lmenum,f.u=rep(0,lmenum),f.v=rep(0,lmenum),f.minu=rep(0,lmenum),f.minv=rep(0,lmenum),rmse=rep(0,lmenum))
# for (i in 1:lmenum) LMEs[i,c(2:6)]<-LHSsearch(LMEs[i,1],iter=100)
# saveRDS(LMEs,"bestvals_LMEs.RDS")
source("LME_calibration.R")
# faster using pbsapply, in the LHSsearch pbapply has cl=6 which uses cluster to run in parallel, but here it is run sequentially if cl is not specified.
lmenum=1 #66
no_iter = 100 # 100 
no_cores <- parallel::detectCores() - 1
lmes<-t(pbapply::pbsapply(X=1:lmenum,LHSsearch,iter=no_iter))
saveRDS(lmes,paste0("Output/bestvals_LMEs_iter_",no_iter,".RDS"))

############### Make plots

#### Check other model performance indicators using the above estimates
bestvals<-data.frame(readRDS(paste0("Output/bestvals_LMEs_iter_",no_iter,".RDS")))

### WARNING for some LME the above did not work... e.g. LME 4. 

# add column for correlation:
bestvals$cor<-rep(0,lmenum)

# pdf(paste0("CalibrationPlots_iter_",no_iter,".pdf"),height = 6, width = 8)

for (i in 1:66){
  
  # trial 
  i = 1
  
  lme_input<-get_lme_inputs(LMEnumber=i, gridded=F,yearly=F)
  
  out<-run_model(vals = unlist(bestvals[i,c(1:4)]),input=lme_input,withinput=T)
  
  ### CN this is copy-paste from getError(): 
  ## aggregate by year (mean to conserve units)
  out <- out %>% group_by(Year) %>% 
    filter(Year > "1949") %>% 
    summarise(TotalCatchPerYr=mean(Totalcatch),
              ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))
  # convert units
  # CN from tonnes to g
  out$ObsCatchPerYr<-out$ObsCatchPerYr*1e6
  
  ### CHECK OPTIONS: 
  
  #### 1 # raw effort (early values each month) 
  # Year TotalCatchPerYr ObsCatchPerYr
  # <dbl>           <dbl>         <dbl>
  #   1  1950        0.000993         0.220
  # 2  1951        0.00107          0.246
  # 3  1952        0.00101          0.276
  # 4  1953        0.00103          0.276
  # 5  1954        0.00108          0.273
  # 6  1955        0.00105          0.246
  # 7  1956        0.00103          0.278
  # 8  1957        0.00105          0.289
  # 9  1958        0.00119          0.299
  # 10  1959        0.00130          0.311
  
  #### 2 # effort/12 (early values/12 each month) 
  # Year TotalCatchPerYr ObsCatchPerYr
  # <dbl>           <dbl>         <dbl>
  #   1  1950     0.0000828        0.0183
  # 2  1951       0.0000895        0.0205
  # 3  1952       0.0000845        0.0230
  # 4  1953       0.0000861        0.0230
  # 5  1954       0.0000899        0.0227
  # 6  1955       0.0000872        0.0205
  # 7  1956       0.0000861        0.0232
  # 8  1957       0.0000873        0.0240
  # 9  1958       0.0000994        0.0249
  # 10  1959       0.000109         0.0259
  
  #### 3 # change back to old effort input file 
  ## checked - same values and trends for LME 1 across effort datasets. 
  # the trend in modelled catches for this version follows the trend in effort BUT 
  # the trend in modelled catches for the old version where LME calibration was working do not.  
  
  
  # bestvals$rmse[i]<-sqrt(sum(out$squared_error,na.rm=T)/sum(!is.na(out$squared_error)))
  
  bestvals$cor[i]<-cor(out$ObsCatchPerYr,out$TotalCatchPerYr,use="complete.obs")
  
  p1<-ggplot() +
    geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
    geom_point(data = out, aes(x = Year, y = ObsCatchPerYr)) +
    theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
                            axis.text.y = element_text(colour="grey20", size=12),
                            text=element_text(size=16)) + 
    labs(x = 'Year',
         y = 'Grams per m2 per year') 
  
  
  modelled <- ggplot() +
    geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
    theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
                            axis.text.y = element_text(colour="grey20", size=12),
                            text=element_text(size=16)) + 
    labs(x = 'Year',
         y = 'Grams per m2 per year') 
  
  p2<-ggplot() +
    geom_point(data = out, aes(x = TotalCatchPerYr, y = ObsCatchPerYr)) +
    geom_abline(slope=1,intercept=0) +
    theme_classic() + theme(axis.text.x = element_text(colour="grey20", size=12, angle=90, hjust=.5, vjust=.5),
                            axis.text.y = element_text(colour="grey20", size=12),
                            text=element_text(size=16)) + 
    labs(x = 'Predicted',
         y = 'Observed',title=paste("LME ",i,sep=""))
  
  p3<-p1+p2
  
  print(p3)
  
}

dev.off()

saveRDS(bestvals,paste0("Output/bestvals_LMEs_cor_",no_iter,".RDS"))

#### TO DO: Check other model performance indicators using the above estimates

# mean total consumer biomass density
# biomass vs sst relationship
# biomass vs phyto (+ zoo) biomass relationship
# slope and intercept of size spectrum
# P:B ratio
# mean size and TL of catch
# correlation with catch time series
# correlation with relative biomass time series (and RAM Legacy)
# OTHER : mse, rpi - hipsey et al metrics for total catch and for disaggregated benthic ad pelagic catch

#### TO DO: Refine estimates by running optimParallel across LMEs using "bestvals" as initial values.

# library(optimParallel)
#  set up workers
# noCores <- parallel::detectCores() - 1 # keep some spare core
# cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
# setDefaultCluster(cl = cl)
# clusterExport(cl, varlist = "cl",envir=environment())
# clusterEvalQ(cl, {
#   library(optimParallel)
#   source("LME_calibration.R")
# })
# 
# LMEnum=22
# lme_input<-get_lme_inputs(LMEnumber=LMEnum)
# optim_out<- optimParallel::optimParallel(par=as.numeric(bestvals[LMEnum,c(2:5)]),getError, input=lme_input, LMEnumber=LMEnum)
# stopCluster(cl)
# 
# bestvals[LMEnum,c(1:4)]<-optim_out$par
# 
# bestvals[LMEnum,5]<-optim_out$value
# 
# saveRDS(bestvals,"optim_vals_LME.RDS")
# 
# 
# ########## to save all error estimates:
# 
# # # set up parameter set to run initial simulations
# # set.seed(1234)
# # 
# # # num "individual runs"
# # num_iter=100
# # simset <- data.frame(randomLHS(num_iter, 4))
# # # rows are iterations, columns are specific parameters
# # colnames(simset )<-c("f.u","f.v","f.minu","f.minv")
# # # adjust range of mi size params, others go form 0-1
# # simset [,"f.minu"]<-simset [,"f.minu"]*2
# # simset [,"f.minv"]<-simset [,"f.minv"]*2
# # 
# # LMEnumber=c(1:66)
# # 
# # simset<-expand_grid(LMEnumber,simset)
# # 
# # simset$rmse<-0
# # 
# # #  set up workers for parallel runs
# # noCores <- parallel::detectCores() - 1 # keep some spare core
# # cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
# # setDefaultCluster(cl = cl)
# # clusterExport(cl, varlist = "cl",envir=environment())
# # clusterEvalQ(cl, {
# #   library(optimParallel)
# #   source("LME_calibration.R")
# # })
# # 
# # # run model and output error term for each row (in parallel)
# # simset$rmse<-pbapply(X=simset[,c(1:5)],1,getError,cl=cl)
# # 
# # stopCluster(cl)
# # 
# # #save the file
# # saveRDS(simset,"simset_LMEs.RDS")
# 
