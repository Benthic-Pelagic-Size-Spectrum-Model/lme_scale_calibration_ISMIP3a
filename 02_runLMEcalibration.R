###### Run LHS search for each lme to get an initial rough guess for 
# fishing parameters
# lmenum=66
# LMEs<-data.frame(LMEnum=1:lmenum,f.u=rep(0,lmenum),f.v=rep(0,lmenum),
#f.minu=rep(0,lmenum),f.minv=rep(0,lmenum),rmse=rep(0,lmenum))
# for (i in 1:lmenum) LMEs[i,c(2:6)]<-LHSsearch(LMEs[i,1],iter=100)
# saveRDS(LMEs,"bestvals_LMEs.RDS")

source("LME_calibration.R")
library(tictoc)
library(parallel)
library(pbapply)

# faster using pbsapply, in the LHSsearch pbapply has cl=6 which uses cluster 
#to run in parallel, but here it is run sequentially if cl is not specified.
lmenum <- 66  
no_iter <- 100
# other option is to specify a value for search_vol 
search_vol <- "estimated" 
no_cores <- parallel::detectCores() - 1
tic()
lmes <- t(pbapply::pbsapply(X = 1:lmenum, LHSsearch, iter = no_iter, 
                            search_vol = search_vol)) 
toc()
saveRDS(lmes, paste0("Output/bestvals_LMEs_searchvol_", search_vol, "_iter_",
                     no_iter, ".RDS"))

# WARNINGS: 
# 2. some LME (e.g. LME 4) do not run because the model called by LHSsearch 
# does not produce biomass 
# increase to 64 # even worst 
# decrease to 0.064 # now we have biomass!

# 3. catches are always too small compared to observed data 
# F.mort estimated in LHSsearch can only go to 1, 
# so increase effort in get_lme_inputs() by 
# using relative effort (effort_m2/max(effort_m2), with the highest value 
# been 1)  
# relative effort for each LME - not working wellif search_vol = 0.064 as 
# catches remain low  
# doing this also means that effort is equal across LMEs (from 1 to close to 0
# in each LME)
# relative effort across LMEs - same as above and do not use (Julia)   

# now try increasing search_vol again to 0.64
# relative effort for each LME - not working as well as the above
# relative effort across LMEs - working well (best option) but LME 4 and others 
# not working again. 

# now estimating search vol + relative effort for each LME + iter = 100 (but 
# will need to increase to 1000 at least)

############### Make plots

#### Check other model performance indicators using the above estimates
#bestvals<-data.frame(readRDS("bestvals_LMEs.RDS")) # these bestvalues don't 
#give the CalibrationPlot 
bestvals <- data.frame(readRDS(paste0("Output/bestvals_LMEs_searchvol_", 
                                      search_vol, "_iter_", no_iter, ".RDS")))

# add column for correlation:
bestvals$cor <- rep(0, lmenum)
# add column for NA catches:
bestvals$catchNA <- rep(0, lmenum)

pdf(paste0("Output/CalibrationPlots_searchvol_", search_vol, "_iter_", 
           no_iter, ".pdf"), height = 6, width = 8)

for(i in 1:66){
  # # trial 
  # i = 1

  ### has the function above worked for all LMEs?
  # run only if the function above produced bestvalues
  if(length(unlist(bestvals[i, c(1:5)])) > 0){ 
  lme_input <- get_lme_inputs(LMEnumber = i, gridded = F, yearly = F)
  # try with highest possible values of Fmort to see if catch are in line with
  # observed
  out <- run_model(vals = unlist(bestvals[i, 1:5]), input = lme_input, 
                   withinput = T)
   
  ### CN this is copy-paste from getError(): 
  ## aggregate by year (mean to conserve units)
  out <- out |> 
    group_by(Year) |> 
    filter(Year > "1949") |> 
    summarise(TotalCatchPerYr = mean(Totalcatch),
              ObsCatchPerYr = mean(catch_tonnes_area_m2, na.rm=T))
  # convert units
  out$ObsCatchPerYr <- out$ObsCatchPerYr*1e6
  
  # bestvals$rmse[i]<-sqrt(sum(out$squared_error,na.rm=T)/
  #sum(!is.na(out$squared_error)))
  
  bestvals$cor[i] <- cor(out$ObsCatchPerYr, out$TotalCatchPerYr,
                         use = "complete.obs")
  bestvals$catchNA[i] <- sum(is.na(out$TotalCatchPerYr))
  
  p1 <- ggplot() +
    geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
    geom_point(data = out, aes(x = Year, y = ObsCatchPerYr)) +
    theme_classic() + 
    theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90,
                                     hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(colour = "grey20", size = 12),
          text = element_text(size = 16)) + 
    labs(x = "Year", y = "Grams per m2 per year") 
  
  modelled <- ggplot() +
    geom_line(data = out, aes(x = Year, y = TotalCatchPerYr)) +
    theme_classic() +
    theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90,
                                     hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(colour = "grey20", size = 12),
          text = element_text(size = 16)) + 
    labs(x = "Year", y = "Grams per m2 per year") 
  
  p2 <- ggplot() +
    geom_point(data = out, aes(x = TotalCatchPerYr, y = ObsCatchPerYr)) +
    geom_abline(slope = 1, intercept = 0) +
    theme_classic() + 
    theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90,
                                     hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(colour = "grey20", size = 12),
          text = element_text(size = 16)) + 
    labs(x = "Predicted", y = "Observed", title = paste0("LME ", i))
  
  p3 <- p1+p2
  
  print(p3)
  }
  # end of if(unlist ...)
}


dev.off()
# saveRDS(bestvals,paste0("Output/bestvals_LMEs_cor_searchvol_", 
# search_vol,"_iter_",no_iter,".RDS"))

#### TO DO: Check other model performance indicators using the above estimates

# mean total consumer biomass density
# biomass vs sst relationship
# biomass vs phyto (+ zoo) biomass relationship
# slope and intercept of size spectrum
# P:B ratio
# mean size and TL of catch
# correlation with catch time series
# correlation with relative biomass time series (and RAM Legacy)
# OTHER : mse, rpi - hipsey et al metrics for total catch and for disaggregated
# benthic ad pelagic catch

#### TO DO: Refine estimates by running optimParallel across LMEs using 
# "bestvals" as initial values.

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
# optim_out<- optimParallel::optimParallel(par=as.numeric(bestvals[LMEnum,
# c(2:5)]),getError, input=lme_input, LMEnumber=LMEnum)
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

# Use optimParallel to get better "bestvals" for LMES that do not have good 
#enough fit to data

# refine<-which(bestvals$cor<0.5|bestvals$rmse>0.5|bestvals$catchNA>0)
# tic()
# for (i in 1:dim(bestvals[refine,])[1]){
#   vals<-unlist(bestvals[refine,][i,1:5])
#   optim_result<-fastOptim(lme=refine[i],vary=vals)
#   bestvals[refine,][i,1:5]<-unlist(optim_result$par)[1:5]
#   bestvals[refine,][i,6]<-unlist(optim_result$value)
#   print(i)
# }
# toc()

refine <- which(bestvals$cor < 0.5 | bestvals$rmse > 0.5 | bestvals$catchNA > 0)

tic()
for(i in 1:dim(bestvals[refine, ])[1]){
  vals <- unlist(bestvals[refine, ][i, 1:5])
  result <- LHSsearch(X = refine[i], iter = 1000)
  bestvals[refine, ][i, 1:6] <- result
  print(i)
}

toc()

saveRDS(bestvals, paste0("Output/refined_bestvals_LMEs_cor_searchvol_", 
                         search_vol, "_iter_", no_iter, ".RDS"))

## then need to put these back with  the other "bestvals".

# test_lme4<-LHSsearch(4,iter=100)
# vals<-unlist(bestvals[3,1:5])
# newvals<-fastOptim(lme=2,unlist(vals),getError)

