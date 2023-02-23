##### run model across space and time using gridded inputs for each LME

source("LME_calibration.R")

# get initial values from LME-scale results
lme_input_14<-get_lme_inputs(LMEnumber = 14,gridded = F)
vals <- readRDS("bestvals_LMEs.RDS")
# run model using time-averaged inputs
initial_results<-run_model(vals=vals[14,],input=lme_input_14,withinput = F)
U.initial<-rowMeans(initial_results$U[,240:1440])
V.initial<-rowMeans(initial_results$V[,240:1440])
W.initial<-mean(initial_results$W[240:1440])
# plot to check initial values
plotsizespectrum(initial_results,params=initial_results$params,itime=240:1440,timeaveraged = TRUE)

# get gridded inputs and run through all grid cells one timestep at a time

lme_inputs_grid<-
  get_lme_inputs(LMEnumber = 14,gridded = T)[,c("lat","lon", "LME.x", "t","sst",
                                                "sbt","er","intercept","slope",
                                                "depth","NomActive_area_m2" )]
time<-unique(lme_inputs_grid$t)
grid_results<-vector("list", length(time))

##### reorganise outputs storage for later
ntime<-length(time)
ngrid<-dim(subset(lme_inputs_grid,t==time[1]))[1]
params<-initial_results$params

# U=V=Y.u=Y.v=GG.u=GG.v=PM.u=PM.v=array(0,c(ngrid,params$Nx,ntime))
# W=array(0,c(ngrid,params$Nx,ntime))
# 
# results<-list(U=U,V=V,Y.u=Y.u,Y.v=Y.v,PM.u=PM.u,PM.v=PM.v,GG.u=GG.u, GG.v=GG.v)
# 


# ####################### RUN IT
# 
# for (itime in 1:length(time)){
#   input_igrid<-
#     subset(lme_inputs_grid,t==time[itime])[,c("depth", "er","intercept","slope",
#                                               "sst","sbt","NomActive_area_m2")]
#   grid_results[[itime]]<-pbapply(X=as.matrix(input_igrid),c(1),run_model_timestep,
#                                  vals = unlist(vals), U.initial=U.initial, 
#                                  V.initial=V.initial, W.initial=W.initial)
#   #spread effort
#   bio<-rep(0,length=132)
#   x<-params$x
#   dx<-params$dx
#   Nx<-params$Nx
#   fmin.u<-params$Fref.u
#   fmin.v<-params$Fref.v
#   
#   for (i in 1:132) 
#     bio[i]<-sum(grid_results[[itime]][[i]]$U[fmin.u:Nx,2]*10^x[fmin.u:Nx]*dx) + 
#     sum(grid_results[[itime]][[i]]$V[fmin.v:Nx,2]*10^x[fmin.v:Nx]*dx)
#   
#   prop_bio<-bio/sum(bio)
#   input_igrid$NomActive_area_m2<-prop_bio*input_igrid$NomActive_area_m2
#   #rerun
#   grid_results[[itime]]<-pbapply(X=as.matrix(input_igrid),c(1),run_model_timestep,vals = unlist(vals), U.initial=U.initial, V.initial=V.initial, W.initial=W.initial)
# }
# 
# saveRDS(grid_results,"LME_14.rds")
# 
# grid_results<-readRDS("LME_14.rds")
# 
# 
# 
# grid_results<-readRDS("LME_14.rds")
# 
# for (itime in 1:ntime){
#   for (igrid in 1:ngrid){
# results$U[igrid,,itime]<-grid_results[[itime]][[igrid]]$U[,2]
# results$V[igrid,,itime]<-grid_results[[itime]][[igrid]]$V[,2]
# results$Y.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$Y.u[,1]
# results$Y.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$Y.v[,2]
# results$PM.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$PM.u[,1]
# results$PM.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$PM.v[,1]
# results$GG.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$GG.u[,1]
# results$GG.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$GG.v[,1]
# results$W[igrid,itime]<-grid_results[[itime]][[igrid]]$W[2]
# 
#   }
# }
# 
# 
# ### get outputs into lat/lon summarise biomass, catches etc.
# 
# agg_outputs<-function(input=lme_inputs_grid,results=results,params=params){
#   # returns all outputs of the model 
#   # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
#     for (itime in 1: ntime){
#     input[input$t==time[itime],]$TotalUbiomass <- apply(results$U[,params$ref:params$Nx,itime]*params$dx*10^params$x[params$ref:params$Nx],1,sum)*min(params$depth,100)
#     input[input$t==time[itime],]$TotalVbiomass <- apply(results$V[,params$ref:params$Nx,itime]*params$dx*10^params$x[params$ref:params$Nx],1,sum)*min(params$depth,100)
#    # input[input$t==time[itime]]$W <- results$W[,itime]*min(params$depth,100)
#   #sum catches (currently in grams per m3 per year, across size classes) 
#   #keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
#     input[input$t==time[itime]]$TotalUcatch <- apply(results$Y.u[,params$ref:params$Nx,itime]*params$dx,2,sum)*min(params$depth,100)
#     input[input$t==time[itime]]$TotalVcatch <- apply(results$Y.v[,params$ref:params$Nx,itime]*params$dx,2,sum)*min(params$depth,100)
#     input[input$t==time[itime]]$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
#     
#     return(input)
#   }
#   
# }
#   
# 
# output<-agg_outputs(input=lme_inputs_grid,results=results)
#   
#   

# ####################### RUN IT
# ## Not working
# for (itime in 1:length(time)){
#   input_igrid<-
#     subset(lme_inputs_grid,t==time[itime])[,c("depth", "er","intercept","slope",
#                                               "sst","sbt","NomActive_area_m2")]
#   grid_results[[itime]]<-pbapply(X=as.matrix(input_igrid),c(1),run_model_timestep,
#                                  vals = unlist(vals), U.initial=U.initial, 
#                                  V.initial=V.initial, W.initial=W.initial)
#   #spread effort
#   bio<-rep(0,length=132)
#   x<-params$x
#   dx<-params$dx
#   Nx<-params$Nx
#   fmin.u<-params$Fref.u
#   fmin.v<-params$Fref.v
#   
#   for (i in 1:132) 
#     bio[i]<-sum(grid_results[[itime]][[i]]$U[fmin.u:Nx,2]*10^x[fmin.u:Nx]*dx) + 
#     sum(grid_results[[itime]][[i]]$V[fmin.v:Nx,2]*10^x[fmin.v:Nx]*dx)
#   
#   prop_bio<-bio/sum(bio)
#   input_igrid$NomActive_area_m2<-prop_bio*input_igrid$NomActive_area_m2
#   #rerun
#   grid_results[[itime]]<-pbapply(X=as.matrix(input_igrid),c(1),run_model_timestep,vals = unlist(vals), U.initial=U.initial, V.initial=V.initial, W.initial=W.initial)
# }
# 
# saveRDS(grid_results,"LME_14.rds")
# 
# grid_results<-readRDS("LME_14.rds")
# 
# #not working
# for (itime in 1:ntime){
#   for (igrid in 1:ngrid){
#     results$U[igrid,,itime]<-grid_results[[itime]][[igrid]]$U[,2]
#     results$V[igrid,,itime]<-grid_results[[itime]][[igrid]]$V[,2]
#     results$Y.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$Y.u[,1]
#     results$Y.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$Y.v[,2]
#     results$PM.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$PM.u[,1]
#     results$PM.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$PM.v[,1]
#     results$GG.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$GG.u[,1]
#     results$GG.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$GG.v[,1]
#     results$W[igrid,itime]<-grid_results[[itime]][[igrid]]$W[2]
#     
#   }
# }
# 
# ### get outputs into lat/lon summarise biomass, catches etc.
# 
# agg_outputs<-function(input=lme_inputs_grid,results=results,params=params){
#   # returns all outputs of the model 
#   # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
#   for (itime in 1: ntime){
#     input[input$t==time[itime],]$TotalUbiomass <- apply(results$U[,params$ref:params$Nx,itime]*params$dx*10^params$x[params$ref:params$Nx],1,sum)*min(params$depth,100)
#     input[input$t==time[itime],]$TotalVbiomass <- apply(results$V[,params$ref:params$Nx,itime]*params$dx*10^params$x[params$ref:params$Nx],1,sum)*min(params$depth,100)
#     # input[input$t==time[itime]]$W <- results$W[,itime]*min(params$depth,100)
#     #sum catches (currently in grams per m3 per year, across size classes) 
#     #keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
#     input[input$t==time[itime]]$TotalUcatch <- apply(results$Y.u[,params$ref:params$Nx,itime]*params$dx,2,sum)*min(params$depth,100)
#     input[input$t==time[itime]]$TotalVcatch <- apply(results$Y.v[,params$ref:params$Nx,itime]*params$dx,2,sum)*min(params$depth,100)
#     input[input$t==time[itime]]$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
#     
#     return(input)
#   }
#   
# }
# 
# # not working
# output<-agg_outputs(input=lme_inputs_grid,results=results)

###compare with catch data

### TO DO: test against code without fishing

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

er_grid[,3:dim(er_grid)[2]] <- er_grid[,2]
intercept_grid[,3:dim(intercept_grid)[2]] <- intercept_grid[,2]
slope_grid[,3:dim(slope_grid)[2]] <- slope_grid[,2]
sst_grid[,3:dim(sst_grid)[2]] <- sst_grid[,2]
sbt_grid[,3:dim(sbt_grid)[2]] <- sbt_grid[,2]

# set up params for each month, across grid cells
gridded_params <- sizeparam (equilibrium = FALSE
                             ,dx = 0.1
                             ,xmin.consumer.u = -3
                             ,xmin.consumer.v = -3
                             ,tmax = dim(er_grid[,-1])[2]/12
                             ,tstepspryr  =  12
                             ,search_vol = 0.64
                             # ,fmort.u = f.u
                             ,fmort.u = 0
                             ,fminx.u = f.minu
                             # ,fmort.v = f.v
                             ,fmort.v = 0
                             ,fminx.v = f.minv
                             ,depth = data.matrix(depth_grid[,-1][,1])
                             ,er = data.matrix(er_grid[,-1])
                             ,pp = data.matrix(intercept_grid[,-1])
                             ,slope = data.matrix(slope_grid[,-1])
                             ,sst = data.matrix(sst_grid[,-1])
                             ,sft = data.matrix(sbt_grid[,-1])
                             ,use.init = TRUE,effort = data.matrix(effort_grid[,-1])
                             ,U.initial =U.initial
                             ,V.initial = V.initial
                             ,W.initial = W.initial
                             ,Ngrid=dim(depth_grid)[1])      

# run model  for full time period across all grid cells
grid_results<-gridded_sizemodel(gridded_params,ERSEM.det.input=F,U_mat,V_mat,W_mat,temp.effect=T,eps=1e-5,output="aggregated",
                                use.init = TRUE, burnin.len)

out<-getGriddedOutputs(input=lme_inputs_grid,results=grid_results,params=gridded_params)

#### CHECK OUTPUTS!!

cells<-unique(out$cell)

ggplot(filter(out,cell==cells[1]), aes(x=t,y=TotalUbiomass)) + geom_line()
ggplot(filter(out,cell==cells[1]), aes(x=t,y=TotalVbiomass)) + geom_line()
ggplot(filter(out,cell==cells[1]), aes(x=t,y=Totalcatch)) + geom_line()
