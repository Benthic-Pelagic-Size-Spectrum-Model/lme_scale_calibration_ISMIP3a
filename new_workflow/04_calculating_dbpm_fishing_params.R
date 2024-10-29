
# Loading libraries -------------------------------------------------------
library(arrow)
library(lhs)
library(dplyr)
library(jsonlite)
source("new_workflow/useful_functions.R")


# Loading DBPM climate and fishing inputs ---------------------------------
dbpm_inputs <- file.path("/g/data/vf71/la6889/dbpm_inputs/east_antarctica", 
                         "monthly_weighted", 
                         "dbpm_clim-fish-inputs_fao-58_1841-2010.parquet") |> 
  read_parquet()


# Defining initial fishing parameters -------------------------------------
set.seed(1234)

#Number of rows to be included in fishing parameters data frame
num_iter <- 100

#Construct a hypercube with random numbers. Columns represent five specific 
#parameters needed to run DBPM
fishing_params <- data.frame(randomLHS(num_iter, 5))
#Renaming columns 
colnames(fishing_params) <- c("fmort_u", "fmort_v", "fminx_u", "fminx_v", 
                              "search_vol")

#Adjust range of mi size params, others go from 0-1
fishing_params <- fishing_params |> 
  mutate(fminx_u = fminx_u*2, 
         fminx_v = fminx_v*2,
         # adjust range of search vol, others go from 0-1
         search_vol = search_vol+0.001)


# Getting DBPM parameters -------------------------------------------------
params <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3, tstepspryr = 12)

# Saving parameters
params |> 
  write_json("new_workflow/outputs/dbpm_size_params.json")


# Loading DBPM parameters -------------------------------------------------
# If paramaters were already saved, they can be read instead of being 
# recalculated
params <- read_json("new_workflow/outputs/dbpm_size_params.json", 
                    simplifyVector = T)










params <- sizeparam(xmin.consumer.u = -3, xmin.consumer.v = -3, 
                    tmax = nrow(input)/12, tstepspryr = 12, 
                    search_vol =  0.39221252,#vals["search.vol"],
                    fmort.u = 0.79832977,#vals["f.u"], 
                    fminx.u = 1.9683917,#vals["f.minu"], 
                    fmort.v = 0.83945284,#vals["f.v"], 
                    fminx.v = 0.82696956,#vals["f.minv"], 
                    depth = mean(input$depth), er = input$export_ratio,
                    pp = input$intercept, slope = input$slope, 
                    sst = input$tos, sft = input$tob, 
                    effort = input$nom_active_area_m2_relative)



# sizemodel <- function(params, 
ERSEM.det.input = F
temp_effect = T
eps = 1e-5
output="aggregated"
use_init = F
    
# Input parameters to vary:
# time series of intercept of plankton size spectrum (estimated from GCM, 
# biogeophysical model output or satellite data).		
ui0 <- 10^params$pp          
# time series of slope of plankton size spectrum (estimated from GCM, 
#biogeophysical model output or satellite data).  	
r.plank <- params$r.plank
# time series of temperature in water column
sst <- params$sst
#time series of temperature near seabed
sft <- params$sft
#time series of export ratio (read in sizeparam)
sinking.rate <- params$sinking.rate


#--------------------------------------------------------------------------
# Initialising matrices
#--------------------------------------------------------------------------

#q1 is a square matrix holding the log(predatorsize/preysize) for all 
#combinations of sizes
x <- params$x
y <- x

q1  <- matrix(NA, length(x), length(y))
for(i in 1:length(y)){
  q1[, i] <- y[i] - x
}

#q2 is the reverse matrix holding the log(preysize/predatorsize) for all 
#combinations of sizes
q2 <- matrix(-q1, length(x), length(y))	

#matrix for recording the two size spectra 
V <- U <- array(0, c(length(x), params$Neq+1))

#vector to hold detritus biomass density (g.m-2)
W <- array(0,params$Neq+1)

#matrix for keeping track of growth and reproduction from ingested food:
R_v <- R_u <- GG_v <- GG_u <- array(0, c(length(x), params$Neq+1)) 

#matrix for keeping track of predation mortality
PM_v <- PM_u <- array(0, c(length(x), params$Neq+1))   

#matrix for keeping track of catches  
Y_v <- Y_u <- array(0, c(length(x), params$Neq+1))  

#matrix for keeping track of  total mortality (Z)
Z_v <- Z_u <- array(0, c(length(x), params$Neq+1))

#matrix for keeping track of senescence mortality and other (intrinsic) 
#mortality
SM_v <- SM_u <- OM_v <- OM_u <- array(0, length(x))

#empty vector to hold fishing mortality rates at each size class at time
Fvec_v <- Fvec_u <- array(0, c(length(x), params$Neq+1))

#lookup tables for terms in the integrals which remain constant over time
gphi  <- gphi_f(q1, params$q0, params$sd.q)
mphi  <- mphi_f(q2, params$q0, params$sd.q, params$alpha.u)

#lookup table for components of 10^(alpha.u*x)
expax <- expax_f(x, params$alpha.u)

#--------------------------------------------------------------------------
# Numerical integration
#--------------------------------------------------------------------------

# set up with the initial values from param
#(phyto+zoo)plankton size spectrum  
U[1:(params$ref-1), 1] <- params$U.init[1:(params$ref-1)]
# set initial consumer size spectrum 
U[params$ref:120, 1] <- params$U.init[params$ref:120]       
# set initial detritivore spectrum  
V[params$ref.det:120, 1] <- params$V.init[params$ref.det:120]  
# set initial detritus biomass density (g.m^-3) 
W[1] <- params$W.init

if(use_init == TRUE){
  # set up with the initial values from previous run
  #(phyto+zoo)plankton size spectrum  
  U[1:(ref-1), 1] <- U.init[1:(ref-1)]  
  # set initial consumer size spectrum from previous run
  U[ref:length(x), 1] <- U.init[ref:length(x)]
  # set initial detritivore spectrum from previous run
  V[ref.det:length(x), 1] <- V.init[ref.det:length(x)] 
  W[1] <- W.init
}

#intrinsic natural mortality
OM_u <- params$mu0*10^(-0.25*x)
OM_v <- params$mu0*10^(-0.25*x)

#senescence mortality rate to limit large fish from building up in the 
#system
#same function as in Law et al 2008, with chosen parameters gives similar 
#M2 values as in Hall et al. 2006
SM_u <- params$k.sm*10^(params$p.s*(x-params$xs))
SM_v <- params$k.sm*10^(params$p.s*(x-params$xs))

#Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!)

#Fvec[Fref:Nx] = 0.09*(x[Fref:Nx]/ log10(exp(1)) ) + 0.04 
# from Benoit & Rochet 2004 

# here Fmort.u and Fmort.u= fixed catchability term for U and V to be
#estimated along with Fref.v and Fref.u
Fvec_u[params$Fref.u:params$Nx, 1] <- params$Fmort.u*params$effort[1]
Fvec_v[params$Fref.v:params$Nx, 1] <- params$Fmort.v*params$effort[1]

#output fisheries catches per yr at size
Y_u[Fref.u:Nx, 1] <- Fvec_u[Fref.u:Nx, 1]*U[Fref.u:Nx, 1]*10^x[Fref.u:Nx] 
#output fisheries catches per yr at size
Y_v[Fref.v:Nx, 1] <- Fvec_v[Fref.v:Nx, 1]*V[Fref.v:Nx, 1]*10^x[Fref.v:Nx] 

#iteration over time, N [days]
# Initial progress bar
#pb = txtProgressBar(min = 0, max = Neq, initial = 1, style = 3) 

for(i in 1:(Neq)){
  #  setTxtProgressBar(pb, i) # Update progress bar
  # if(W[i]=="NaN"|W[i]<0)
  # {
  #   #browser()
  #   U[,i]<-"NaN"
  #   return(list(U=U[,],GG_u=GG_u[,],PM_u=PM_u[,],V=V[,],GG_v=GG_v[,],
  #PM_v=PM_v[,],Y_u=Y_u[,],Y_v=Y_v[,],W=W[], params=params))
  # }
  
  #  ONLY IF RUNNING TO EQUILIBRIUM:
  #  below is to skip unecessary timesteps if model has reached equilibirum
  # if(i>100 & equilibrium==T )
  # {
  #   # browser()
  #   if(max(abs((log(U[-c(1:91),(i-1)]))-(log(U[-c(1:91),(i-2)]))))<eps
  #      &max(abs((log(V[-c(1:91),(i-1)]))-(log(V[-c(1:91),(i-2)]))))<eps
  #      &(abs((log(W[(i-1)]))-(log(W[(i-2)]))))<eps)
  #   {
  #     U[,Neq]<-U[,i-1]
  #     PM_v[,Neq]<-PM_v[,i-1]
  #     GG_v[,Neq]<-GG_v[,i-1]
  #     V[,Neq]<-V[,i-1]
  #     PM_u[,Neq]<-PM_u[,i-1]
  #     GG_u[,Neq]<-GG_u[,i-1]
  #     Y_u[,Neq]<-Y_u[,i-1]
  #     Y_v[,Neq]<-Y_v[,i-1]
  #     W[Neq]<-W[i-1]
  #     return(list(U=U[,],GG_u=GG_u[,],PM_u=PM_u[,],V=V[,],GG_v=GG_v[,],
  #PM_v=PM_v[,],Y_u=Y_u[,],Y_v=Y_v[,],W=W[], params=params))
  #   }
  # }
  # 
  
  
  #--------------------------------
  # Calculate Growth and Mortality
  #--------------------------------
  if(temp_effect == T){
    pel_tempeffect <- exp(c1-E/(Boltzmann*(sst+273)))
    ben_tempeffect <- exp(c1-E/(Boltzmann*(sft+273)))
  }
  if(temp_effect == F){
    pel_tempeffect <- 1
    ben_tempeffect <- 1
  }
  
  # feeding rates
  # yr-1
  f_pel <- pel_tempeffect[i] * 
    as.vector(((A.u*10^(x*alpha.u)*pref.pel)*(U[, i]*dx)%*%(gphi)) / 
                (1+handling*(A.u*10^(x*alpha.u)*pref.pel) * 
                   (U[, i]*dx)%*%(gphi))) 
  # yr-1
  f_ben <- pel_tempeffect[i] * 
    as.vector(((A.u*10^(x*alpha.u)*pref.ben)*(V[, i]*dx)%*%(gphi)) / 
                (1+handling*(A.u*10^(x*alpha.u)*pref.ben) *
                   (V[, i]*dx)%*%(gphi)))
  # yr-1
  f_det <- ben_tempeffect[i] *
    ((1/10^x)*A.v*10^(x*alpha.v)*W[i]) / 
    (1+handling*(1/10^x)*A.v*10^(x*alpha.v)*W[i])
  
  # Predator growth integral 
  # yr-1
  GG_u[, i] <- (1-def.high)*K.u*(f_pel)+(1-def.low)*K.v*(f_ben)
  
  # Reproduction
  # yr-1
  if(repro.on == 1){
    R_u[, i] <- (1-def.high)*(1-(K.u+AM.u))*(f_pel) +
      (1-def.high)*(1-(K.v+AM.v))*f_ben
  }
  
  # Predator death integrals 
  #Satiation level of predator for pelagic prey
  sat_pel <- ifelse(f_pel > 0,
                    f_pel/((A.u*10^(x*alpha.u)*pref.pel) * 
                             (U[,i]*dx)%*%(gphi)),
                    0)
  # yr-1
  PM_u[, i] <- as.vector((pref.pel*A.u*expax)*(U[, i]*sat_pel*dx)%*%(mphi))
  # yr-1
  #PM_u[,i]<-as.vector((1-f_pel)*(A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(mphi))
  
  # yr-1
  Z_u[, i] <- PM_u[, i]+pel_tempeffect[i]*OM_u+SM_u+Fvec_u[, i]
  
  # Benthos growth integral
  # yr-1
  GG_v[, i] <- (1-def.low)*K.d*f_det
  
  #reproduction
  # yr-1
  if(repro.on == 1){
    R_v[, i] <- (1-def.low)*(1-(K.d+AM.v))*(f_det)
  }
  
  # Benthos death integral
  #Satiation level of predator for benthic prey 
  sat_ben <- ifelse(f_ben > 0,
                    f_ben/((A.u*10^(x*alpha.v)*pref.ben) * 
                             (V[, i]*dx)%*%(gphi)),
                    0)
  # yr-1
  PM_v[, i] <- ifelse(sat_ben > 0,
                      as.vector((pref.ben*A.u*expax) *
                                  (U[, i]*sat_ben*dx)%*%(mphi)),
                      0)
  # yr-1
  #PM_v[,i]<-as.vector((1-f_ben)*(A.u*10^(x*alpha.u)*pref.ben)*(U[,i]*dx)%*%(mphi))
  
  # yr-1
  Z_v[, i] <- PM_v[, i]+ben_tempeffect[i]*OM_v+SM_v+Fvec_v[, i]
  
  #total biomass density eaten by pred (g.m-2.yr-1)
  eatenbypred <- 10^x*f_pel*U[, i]+10^x*f_ben*U[, i] 
  
  #detritus output (g.m-2.yr-1)
  # losses from detritivore scavenging/filtering only:
  output_w <- sum(10^x*f_det*V[, i]*dx)   
  
  #total biomass density defecated by pred (g.m-2.yr-1)
  defbypred <- def.high*(f_pel)*10^x*U[, i]+ def.low*(f_ben)*10^x*U[, i]
  
  #------------------------------------------------
  # Increment values of W,U & V	for next time step  
  #------------------------------------------------
  
  #Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of 
  #detritus pool and solve for detritus biomass density in next time step 
  if(ERSEM.det.input == F){
    #considering pelagic faeces as input as well as dead bodies from both 
    #pelagic and benthic communities and phytodetritus (dying sinking
    #phytoplankton)
    if(det.coupling == 1.0){
      # pelagic spectrum inputs (sinking dead bodies and faeces) - export 
      # ratio used for "sinking rate" + benthic spectrum inputs (dead stuff
      # already on/in seafloor)
      input_w <- (sinking.rate[i] * 
                    (sum(defbypred[ref:Nx]*dx) +
                       sum(pel_tempeffect[i] * 
                             OM_u[1:Nx]*U[1:Nx, i]*10^(x[1:Nx])*dx) + 
                       sum(pel_tempeffect[i] * 
                             SM_u[1:Nx]*U[1:Nx,i]*10^(x[1:Nx])*dx)) +
                    (sum(ben_tempeffect[i] * 
                           OM_v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx) + 
                       sum(ben_tempeffect[i] * 
                             SM_v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx)))
      # )
    }
    if(det.coupling == 0.0){
      input_w <- sum(ben_tempeffect[i] *
                       OM_v[1:Nx]*V[1:Nx, i]*10^(x[1:Nx])*dx) + 
        sum(ben_tempeffect[i]*SM_v[1:Nx]*V[1:Nx, i]*10^(x[1:Nx])*dx)
    }
    
    # get burial rate from Dunne et al. 2007 equation 3
    burial <- input_w*(0.013 + 0.53*input_w^2/(7+input_w)^2)
    
    # change in detritus biomass density (g.m-2.yr-1)                     
    # this one assumes immeidate additional losses to sediment?
    # dW<-input_w - (output_w + W[i])     
    # losses due to detritivory only:
    # dW<-input_w - (output_w)     
    
    # losses from detritivory + burial rate (not including remineralisation
    # bc that goes to p.p. after sediment, we are using realised p.p. as
    # inputs to the model) 
    dW <- input_w-(output_w+burial) 
    #biomass density of detritus g.m-2
    W[i+1] <- W[i]+dW*delta_t
  }
  
  if(ERSEM.det.input == T){
    W[i+1] <- W[i]
  }
  
  #----------------------------------------------
  # Pelagic Predator Density (nos.m-2)- solve for time + delta_t using
  # implicit time Euler upwind finite difference (help from Ken Andersen 
  # and Richard Law)
  
  # Matrix setup for implict differencing 
  Ai_u <- Bi_u <- Si_u <- array(0, c(length(x), 1))   
  
  Ai_u[idx] <- (1/log(10))*-GG_u[idx-1, i]*delta_t/dx
  Bi_u[idx] <- 1+(1/log(10))*GG_u[idx, i]*delta_t/dx +Z_u[idx, i]*delta_t
  Si_u[idx] <- U[idx, i]
  
  # Boundary condition at upstream end 
  Ai_u[ref] <- 0
  Bi_u[ref] <- 1
  Si_u[ref] <- U[ref, i]
  
  # Invert matrix
  #recruitment at smallest consumer mass
  #continuation of plankton hold constant  
  # U[1:ref,i+1]<-U[1:ref,i] 
  
  if(use_init==T){
    U[1:ref, i+1] <- ui0[i]*10^(r.plank[i]*x)[1:(ref)] 
  }
  if(use_init == F){
    U[1:ref, i+1] <- ui0[i+1]*10^(r.plank[i+1]*x)[1:(ref)] 
  }
  
  # apply transfer efficency of 10% *plankton density at same size  
  # if (repro.on==0)  U[ref,i+1]<-0.1*u.init.f(x,ui0,r.plank)[ref]       
  # reproduction from energy allocation
  if(repro.on == 1){
    U[ref, i+1] <- U[ref, i] +
      (sum(R_u[(ref+1):Nx, i]*10^x[(ref+1):Nx]*U[(ref+1):Nx, i]*dx) *
         delta_t)/(dx*10^x[ref])-(delta_t/dx)*(1/log(10)) *
      (GG_u[ref, i])*U[ref, i]-delta_t*Z_u[ref, i]*U[ref, i]
  }
  
  #main loop calculation
  for(j in (ref+1):(Nx)){
    #U[j,i+1]<-U[j,i]-(delta_t/dx)*(1/log(10))*(GG_u[j,i])*U[j,i] + 
    #(delta_t/dx)*(1/log(10))*GG_u[j-1,i]*U[j-1,i] -delta_t*Z_u[j,i]*U[j,i]
    
    U[j, i+1] <- (Si_u[j]-Ai_u[j]*U[j-1, i+1])/Bi_u[j]
  }
  
  #----------------------------------------------
  # Benthic Detritivore Density (nos.m-2) 
  Ai_v <- Bi_v <- Si_v <- array(0, c(length(x), 1))   
  #shorthand for matrix referencing
  idx <- (ref.det+1):Nx  
  
  Ai_v[idx] <- (1/log(10))*-GG_v[idx-1,i]*delta_t/dx 
  Bi_v[idx] <- 1+(1/log(10))*GG_v[idx, i]*delta_t/dx +Z_v[idx, i]*delta_t
  Si_v[idx] <- V[idx, i]
  
  #boundary condition at upstream end
  Ai_v[ref.det] <- 0
  Bi_v[ref.det] <- 1
  Si_v[ref.det] <- V[ref.det, i]  
  
  #invert matrix
  #recruitment at smallest detritivore mass  
  #hold constant continution of plankton with sinking rate multiplier 
  V[1:ref.det, i+1] <- V[1:ref.det, i]
  
  # apply a very low of transfer effiency 1%* total biomass of detritus
  #divided by minimum size
  # if (repro.on==0)  V[ref.det,i+1]<-(0.01*W[i+1])/(10^x[ref]) 
  if(repro.on == 1){
    V[ref.det, i+1] <- V[ref.det, i] +
      sum(R_v[(ref.det+1):Nx, i]*10^x[(ref.det+1):Nx] * 
            V[(ref.det+1):Nx, i]*dx)*delta_t/(dx*10^x[ref.det]) - 
      (delta_t/dx)*(1/log(10))*(GG_v[ref.det, i])*V[ref.det, i] - 
      delta_t*Z_v[ref.det, i]*V[ref.det, i]
  }
  
  #loop calculation
  for(j in (ref.det+1):(Nx)){ 
    #V[j,i+1]<-V[j,i]-(delta_t/dx)*(1/log(10))*(GG_v[j,i])*V[j,i] + 
    #(delta_t/dx)*(1/log(10))*GG_v[j-1,i]*V[j-1,i]-delta_t*Z_v[j,i]*V[j,i] 
    V[j, i+1] <- (Si_v[j]-Ai_v[j]*V[j-1, i+1])/Bi_v[j]
  }		
  rm(j)
  
  # increment fishing 
  Fvec_u[Fref.u:Nx, i+1] <- Fmort.u*effort[i+1]
  Fvec_v[Fref.v:Nx, i+1] <- Fmort.v*effort[i+1]
  
  #output fisheries catches per yr at size
  Y_u[Fref.u:Nx, i+1] <- Fvec_u[Fref.u:Nx, i+1]*U[Fref.u:Nx, i+1] *
    10^x[Fref.u:Nx] 
  #output fisheries catches per yr at size
  Y_v[Fref.v:Nx, i+1] <- Fvec_v[Fref.v:Nx, i+1]*V[Fref.v:Nx, i+1] *
    10^x[Fref.v:Nx] 
}
#end time iteration

return(list(U = U[,],
            GG_u = GG_u[,],
            PM_u = PM_u[,],
            V = V[,], 
            GG_v = GG_v[,],
            PM_v = PM_v[,],
            Y_u = Y_u[,],
            Y_v = Y_v[,],
            W = W[],
            params=params))
# return(list(U=U[,Neq+1],GG_u=GG_u[,Neq],PM_u=PM_u[,Neq],V=V[,Neq+1],
#GG_v=GG_v[,Neq],PM_v=PM_v[,Neq],Y=Y[,Neq],W=W[Neq+1], params=params))
  })  
  # end with(params)
}
#end size-based model function








# 
# 
# 
# # Tuning LHS parameters -----
# LHSsearch <- function(LMEnum, forcing_file, num_iter = 1, 
#                       search_vol = "estimated",
#                       #gridded_forcing = NULL, fishing_effort_file, 
#                       corr = F, figure_folder = NULL,
#                       best_val_folder = NULL){
#   #Inputs:
#   #- LMEnum (numeric) - Unique ID identifying LME
#   #- num_iter (numeric) - Number of individual runs. Default is 1.
#   #- search_vol (???) - Default is "estimated". ???
#   #- forcing_file (character) - Full path to forcing file (non-gridded) that must
#   # contain both climate and fishing (catch and effort) data
#   #- gridded_forcing (character) - Full path to folder containing gridded forcing
#   # files
#   #- fishing_effort_file (character) - Full path to fishing effort file
#   #- corr (boolean) - Default is FALSE. If set to TRUE, it will calculate the
#   # correlation between predicted and observed values
#   #- figure_folder (character) - Optional. If provided, it must be the full path
#   # to the folder where figures comparing observed and predicted data will be 
#   # stored
#   #- best_val_folder (character) - Optional. If provided, it muste be the full
#   # path to the folder where LHS search results will be saved
#   #
#   #Output:
#   #bestvals (data frame) - Contains the values for LHS parameters that resulted
#   #in the best performing model based on RMSE values
#   
#   
#   # use below to select a constant value for search.vol
#   # if(!is.null(forcing_file)){
#   lme_input <- read_parquet(forcing_file)
#   # }
#   # if(!is.null(gridded_forcing)){
#   #   lme_input <- get_lme_inputs(gridded_forcing = gridded_forcing, 
#   #                               fishing_effort_file = fishing_effort_file, 
#   #                               LMEnumber = LMEnum)
#   # }
#   
#   # parallelise using 75% of cores available using mclapply
#   # no_cores <- round((detectCores()*.75), 0)
#   sim$rmse <- mclapply(1:nrow(sim), 
#                        FUN = function(i) getError(unlist(sim[i,]),
#                                                   lme_forcings = lme_input, 
#                                                   corr, figure_folder), 
#                        mc.cores = no_cores) |> 
#     unlist()
#   
#   # check this time param set with lowest error
#   bestvals <- sim |> 
#     filter(rmse == min(rmse, na.rm = T)) |> 
#     mutate(region = LMEnum)
#   
#   #Print row with lowest RMSE
#   print(bestvals)
#   
#   #If folder to save values is provided - Save results
#   if(!is.null(best_val_folder)){
#     #Ensure folder exists
#     if(!dir.exists(best_val_folder)){
#       dir.create(best_val_folder, recursive = T)
#     }
#     
#     #File path to save output
#     fout <- file.path(best_val_folder, 
#                       paste0("best-fishing-parameters_LME_", LMEnum,
#                              "_searchvol_", search_vol, "_numb-iter_", 
#                              num_iter, ".csv"))
#     #Save output
#     bestvals |> 
#       fwrite(fout)
#   }
#   
#   return(bestvals)
# }








gridded_params <- sizeparam(,
                            
                            equilibrium = FALSE, dx = 0.1, 
                            xmin.consumer.u = -3, xmin.consumer.v = -3,
                            tmax = tmax, tstepspryr = 12,
                            search_vol = f.effort["search.vol"], 
                            fmort.u = f.u, fminx.u = f.minu, fmort.v = f.v, 
                            fminx.v = f.minv, depth = lme_inputs_grid$depth, 
                            er = lme_inputs_grid$er, 
                            pp = lme_inputs_grid$intercept, 
                            slope = lme_inputs_grid$slope, 
                            sst = lme_inputs_grid$sst, 
                            sft = lme_inputs_grid$sbt, use.init = TRUE,
                            effort = lme_inputs_grid$nom_active_area_m2_relative, 
                            U.initial = U.initial, V.initial = V.initial, 
                            W.initial = W.initial, 
                            Ngrid = nrow(lme_inputs_grid$depth))


