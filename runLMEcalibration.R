###### Run LHS search for each lme to get an initial rough guess for fihsing parameters
# lmenum=66
# LMEs<-data.frame(LMEnum=1:lmenum,f.u=rep(0,lmenum),f.v=rep(0,lmenum),f.minu=rep(0,lmenum),f.minv=rep(0,lmenum),rmse=rep(0,lmenum))
# for (i in 1:lmenum) LMEs[i,c(2:6)]<-LHSsearch(LMEs[i,1],iter=100)
# saveRDS(LMEs,"bestvals_LMEs.RDS")

# faster using pbsapply, in the LHSsearch pbapply has cl=6 which uses cluster to run in parallel, but here it is run sequentially if cl is not specified.
lmenum=66
lmes<-t(pbsapply(X=1:lmenum,LHSsearch,iter=300,cl=2))
saveRDS(lmes,"bestvals_LMEs_iter300.RDS")

############### Make plots

#### Check other model performance indicators using the above estimates
bestvals<-data.frame(readRDS("bestvals_LMEs.RDS_iter300"))

# add column for correlation:
bestvals$cor<-rep(0,lmenum)

pdf("CalibrationPlots_iter100.pdf",height = 6, width = 8)

for (i in 1:66){
  
  lme_input<-get_lme_inputs(LMEnumber=i)
  
  out<-run_model(unlist(bestvals[i,c(1:4)]),input=lme_input)
  ## aggregate by year (mean to conserve units)
  out <- out %>% group_by(Year) %>% filter(Year > "1949") %>% summarise(TotalCatchPerYr=mean(Totalcatch),ObsCatchPerYr=mean(catch_tonnes_area_m2,na.rm=T))
  # convert units
  out$ObsCatchPerYr<-out$ObsCatchPerYr*1e6
  # out$squared_error <- (out$ObsCatchPerYr- out$TotalCatchPerYr)^2
  # 
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

saveRDS(bestvals,"bestvals_LMEs_cor.RDS")

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
