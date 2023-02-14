##### run model across space and time using gridded inputs for each LME

library(tidyverse)

# get initial values from LME-scale results
lme_input_14<-get_lme_inputs(LMEnumber = 14,gridded = F)
vals <- readRDS("~/Dropbox/DBPM_ISIMIP_3a/lme_scale_calibration_ISMIP3a/bestvals_LMEs.RDS")
initial_results<-run_model(vals=unlist(vals),input=lme_input_14,withinput = F)
U.initial<-rowMeans(initial_results$U[,1:12])
V.initial<-rowMeans(initial_results$V[,1:12])
W.initial<-mean(initial_results$W[1:12])

# get gridded inputs and run through all grid cells one timestep at a time

lme_inputs_grid<-get_lme_inputs(LMEnumber = 14,gridded = T)[,c("lat","lon", "LME.x", "t","sst","sbt","er","intercept","slope","depth","NomActive_area_m2" )]
time<-unique(lme_inputs_grid$t)
grid_results<-vector("list", length(time))

##### reorganise outputs storage for later
ntime<-length(time)
ngrid<-dim(subset(lme_inputs_grid,t==time[1]))[1]
params<-initial_results$params

U=V=Y.u=Y.v=GG.u=GG.v=PM.u=PM.v=array(0,c(ngrid,params$Nx,ntime))
W=array(0,c(ngrid,params$Nx,ntime))

results<-list(U=U,V=V,Y.u=Y.u,Y.v=Y.v,PM.u=PM.u,PM.v=PM.v,GG.u=GG.u, GG.v=GG.v)

####################### RUN IT

for (itime in 1:length(time)){
  input_igrid<-subset(lme_inputs_grid,t==time[itime])[,c("depth", "er","intercept","slope","sst","sbt","NomActive_area_m2")]
  grid_results[[itime]]<-pbapply(X=as.matrix(input_igrid),c(1),run_model_timestep,vals = unlist(vals), U.initial=U.initial, V.initial=V.initial, W.initial=W.initial)
  #spread effort
  bio<-rep(0,length=132)
  x<-params$x
  dx<-params$dx
  Nx<-params$Nx
  fmin.u<-params$Fref.u
  fmin.v<-params$Fref.v
  for (i in 1:132) bio[i]<-sum(grid_results[[itime]][[i]]$U[fmin.u:Nx,2]*10^x[fmin.u:Nx]*dx) + sum(grid_results[[itime]][[i]]$V[fmin.v:Nx,2]*10^x[fmin.v:Nx]*dx)
  prop_bio<-bio/sum(bio)
  input_igrid$NomActive_area_m2<-prop_bio*input_igrid$NomActive_area_m2
  #rerun
  grid_results[[itime]]<-pbapply(X=as.matrix(input_igrid),c(1),run_model_timestep,vals = unlist(vals), U.initial=U.initial, V.initial=V.initial, W.initial=W.initial)
}

saveRDS(grid_results,"LME_14.rds")

grid_results<-readRDS("LME_14.rds")



grid_results<-readRDS("LME_14.rds")

for (itime in 1:ntime){
  for (igrid in 1:ngrid){
results$U[igrid,,itime]<-grid_results[[itime]][[igrid]]$U[,2]
results$V[igrid,,itime]<-grid_results[[itime]][[igrid]]$V[,2]
results$Y.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$Y.u[,1]
results$Y.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$Y.v[,2]
results$PM.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$PM.u[,1]
results$PM.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$PM.v[,1]
results$GG.u[igrid,,itime]<-grid_results[[itime]][[igrid]]$GG.u[,1]
results$GG.v[igrid,,itime]<-grid_results[[itime]][[igrid]]$GG.v[,1]
results$W[igrid,itime]<-grid_results[[itime]][[igrid]]$W[2]

  }
}


### get outputs into lat/lon summarise biomass, catches etc.

agg_outputs<-function(input=lme_inputs_grid,results=results,params=params){
  # returns all outputs of the model 
  # saveRDS(result_set,filename=paste("dbpm_calibration_LMEnumber_catchability.rds"))
    for (itime in 1: ntime){
    input[input$t==time[itime],]$TotalUbiomass <- apply(results$U[,params$ref:params$Nx,itime]*params$dx*10^params$x[params$ref:params$Nx],1,sum)*min(params$depth,100)
    input[input$t==time[itime],]$TotalVbiomass <- apply(results$V[,params$ref:params$Nx,itime]*params$dx*10^params$x[params$ref:params$Nx],1,sum)*min(params$depth,100)
   # input[input$t==time[itime]]$W <- results$W[,itime]*min(params$depth,100)
  #sum catches (currently in grams per m3 per year, across size classes) 
  #keep as grams per m2, then be sure to convert observed from tonnes per m2 per year to g.^-m2.^-yr (for each month)
    input[input$t==time[itime]]$TotalUcatch <- apply(results$Y.u[,params$ref:params$Nx,itime]*params$dx,2,sum)*min(params$depth,100)
    input[input$t==time[itime]]$TotalVcatch <- apply(results$Y.v[,params$ref:params$Nx,itime]*params$dx,2,sum)*min(params$depth,100)
    input[input$t==time[itime]]$Totalcatch <- input$TotalUcatch +   input$TotalVcatch
    
    return(input)
  }
  
}
  

output<-agg_outputs(input=lme_inputs_grid,results=results)
  
  

###compare with catch data

### TO DO: test against code without fishing
