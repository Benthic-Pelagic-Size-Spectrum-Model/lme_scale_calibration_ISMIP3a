#### model to run across grid cells, spreading out effort, uses 3D arrays ------
gridded_sizemodel <- function(params, ERSEM.det.input = F, U_mat, V_mat, W_mat, 
                              temp.effect = T, eps = 1e-5, 
                              output = "aggregated", use.init = F, burnin.len){
  # # trial
  # params = gridded_params
  # ERSEM.det.input=F
  # temp.effect=T
  # eps=1e-5
  # output="aggregated"
  # use.init = TRUE
  # # burnin.len
  
  # WARNING CN - commented for testing 
  with(params, {
    #--------------------------------------------------------------------------
    # Model for a dynamical ecosystem comprised of: two functionally distinct 
    # size spectra (predators and detritivores), size structured primary 
    # producers and an unstructured detritus resource pool. 
    # time implicit upwind difference discretization algorithm (from Ken 
    # Andersen)
    # (Numerical Recipes and Richards notes, see "C://...standardised test 
    # implicit.txt" for basic example)
    # U represents abundance density (m-2) of "fish" and V abundance density
    # (m-2) of "detritivores"
    # W is the pool of detritus (expressed in biomass density, g.m-2) - not 
    # size-based
    # Fish feed on other fish and benthic detritivores with size preference 
    # function
    # Benthic detritivores feed on detritus produced from pelagic spectrum 
    # Senescence Mortality also included.  Options to include dynamic 
    # reproduction and predator handling time (but not currently used).
    #
    # Code modified for global fishing mortality rate application. 
    # JLB 17/02/2014
    # Code modified to include temperature scaling on senescence and detrital 
    # flux. RFH 18/06/2020
    # -------------------------------------------------------------------------
    # Input parameters to vary:
    
    # # trial
    # pp = params$pp
    # r.plank = params$r.plank
    # sst = params$sst
    # sft = params$sft
    # sinking.rate = params$sinking.rate
    # x = params$x
    # Ngrid = params$Ngrid
    # Neq= params$Neq
    # q0 = params$q0
    # alpha.u = params$alpha.u
    # sd.q = params$sd.q
    # U.init = params$U.init
    # V.init = params$V.init
    # ref = params$ref
    # ref.det = params$ref.det
    # W.init = params$W.init
    # mu0 = params$mu0
    # k.sm = params$k.sm
    # p.s =params$p.s
    # xs = params$xs
    # # FISHING
    # Fmort.u = params$Fmort.u
    # effort = params$effort
    # Fmort.v = params$Fmort.v
    # Fref.v = params$Fref.v
    # Fref.u = params$Fref.u
    #
    # Nx = params$Nx
    # c1 = params$c1
    # E = params$E
    # Boltzmann = params$Boltzmann
    # A.u = params$A.u
    # A.v = params$A.v
    # alpha.u = params$alpha.u
    # alpha.v = params$alpha.v
    # pref.pel = params$pref.pel
    # pref.ben = params$pref.ben
    # dx = params$dx
    # handling = params$handling
    # repro.on = params$repro.on
    # def.high = params$def.high
    # K.u = params$K.u
    # K.v = params$K.v
    # AM.u = params$AM.u
    # AM.v = params$AM.v
    # def.low = params$def.low
    # def.high = params$def.high
    # depth = params$depth
    # K.d = params$K.d
    # K.u = params$K.u
    # K.v = params$K.v
    # K.sm =params$k.sm
    # det.coupling = params$det.coupling
    # delta_t = params$delta_t
    # idx = params$idx

    # gridded time series of intercept of plankton size spectrum (estimated 
    # from GCM, biogeophysical model output or satellite data).
    ui0 <- 10^pp          
    # gridded time series of slope of plankton size spectrum (estimated from
    # GCM, biogeophysical model output or satellite data).
    r.plank <- r.plank   
    
    # gridded time series of temperature in water column
    sst <- sst            
    #gridded time series of temperature near seabed
    sft <- sft		      
    # gridded time series of export ratio ( read in sizeparam)
    sinking.rate <- sinking.rate  	      
    
    
    
    #--------------------------------------------------------------------------
    # Functions called within sizemodel()
    #--------------------------------------------------------------------------
    
    
    # Function to build a lookup table for diet preference of all combinations 
    # of predator and prey body size: diet preference (in the predator spectrum
    # only)
    
    phi.f = function(q){
      q = q
      q0 = q0
      sd.q = sd.q
      
      phi = ifelse(q > 0, 
                   exp(-(q-q0)*(q-q0)/(2*sd.q*sd.q))/(sd.q*sqrt(2.0*pi)),
                   0) 
      
      return(phi)
      
      # or normalise feeding kernel to sum to 1:
      # return(phi/sum(phi))
      # this is commented out because in previous work normalising did not 
      # produce realistic growth rates
      
      
    }
    
    #Functions to build lookup tables for components of integration which remain
    #constant
    #growth
    gphi.f = function(q1){
      return(10^(-q1)*phi.f(q1))
      }	
    #mortality
    mphi.f = function(q2){
      return(10^(alpha.u*q2)*phi.f(q2))
      }	
    
    #Function to build lookup table for components of 10^(alpha*x)
    expax.f = function(x, alpha.u){
      return(10^(alpha.u*x)) 
      }
    
    #Function to compute convolution products: can be for growth or mortality
    #depending on phi
    convolution.f = function(phi, u) {
      res = matrix(NA, length(x), 1)
      for(i in 1:length(res)){
        res[i] = 0.0
        res[i] = sum(phi[,i] * u * dx)
      }
      return(res)
    }
    
    #--------------------------------
    # Growth and Mortality Equations
    #--------------------------------
    
    #predators
    #death.u<-function(u,A,expax,mphi){return(pref.pel*A*expax*(convolution.f(mphi,u)))}
    #use faster matrix method instead of convolution loop function
    death.u <- function(u, A, expax, mphi){
      return((pref.pel*A*expax)*(u*dx)%*%(mphi))
      } 
    #growth.u<-function(K0,K1,A,expax,gphi,u,v){return((pref.pel*K0*A*expax*(convolution.f(gphi,u)))+
    # (pref.ben*K1*A*expax*(convolution.f(gphi,v))))}
    growth.u <- function(K0, K1, A, expax, gphi, u, v){
      return((pref.pel*K0*A*expax)*(u*dx)%*%(gphi) + 
               (pref.ben*K1*A*expax)*(v*dx)%*%(gphi))
      }
    #detritivores
    #death.v<-function(u,A,expax,mphi){return((pref.ben)*A*expax*(convolution.f(mphi,u)))}
    death.v <- function(u, A, expax, mphi){
      return((pref.ben*A*expax)*(u*dx)%*%(mphi))
      }
    growth.v <- function(K1, A, w){
      return((1/10^x)*(K1*A*10^(x*0.75)*w))
      }
    #detritus output (g.m-3.yr-1)
    out.w <- function(A, v, w){
      return(sum((A*(10^(x*0.75))*w*v)*dx))}
    
    #--------------------------------------------------------------------------
    # Initialising matrices
    #--------------------------------------------------------------------------
    
    #q1 is a square matrix holding the log(predatorsize/preysize) for all 
    #combinations of sizes
    
    y = x
    
    q1  = matrix(NA, length(x), length(y))
    for(i in 1:length(y)){
      q1[,i] = y[i] - x
      }
    
    #q2 is the reverse matrix holding the log(preysize/predatorsize) for all 
    #combinations of sizes
    q2 = matrix(-q1, length(x), length(y))	
    
    #matrix for recording the two size spectra 
    V = U = array(0, c(Ngrid, length(x), Neq+1))
    
    #vector to hold detrtitus biomass density (g.m-2)
    W = array(0, c(Ngrid, Neq+1))
    
    #matrix for keeping track of growth and reproduction from ingested food:
    R.v = R.u = GG.v = GG.u = array(0, c(Ngrid, length(x), Neq+1))
    
    #matrix for keeping track of predation mortality
    PM.v = PM.u = array(0, c(Ngrid, length(x), Neq+1))   
    
    #matrix for keeping track of catches  
    Y.v = Y.u = array(0, c(Ngrid, length(x), Neq+1)) 
    
    #matrix for keeping track of total mortality (Z)
    Z.v = Z.u = array(0, c(Ngrid, length(x), Neq+1))
    
    #matrix for keeping track of senescence mortality and other (intrinsic) 
    #mortality
    SM.v = SM.u = OM.v = OM.u = array(0, length(x))
    
    #empty vector to hold fishing mortality rates at each size class at time
    Fvec.v = Fvec.u = array(0, c(Ngrid, length(x), Neq+1))
    
    
    #lookup tables for terms in the integrals which remain constant over time
    gphi  = gphi.f(q1)
    mphi  = mphi.f(q2)
    
    #lookup table for components of 10^(alpha.u*x)
    expax = expax.f(x, alpha.u)
    
    
    #--------------------------------------------------------------------------
    # Numerical integration
    #--------------------------------------------------------------------------
    
    # set up with the initial values from param - same for all grid cells
    for (j in 1:Ngrid){
      #(phyto+zoo)plankton size spectrum  
      U[1:Ngrid, 1:(ref-1), 1] <- U.init[1:(ref-1)]
      # set initial consumer size spectrum 
      U[1:Ngrid, ref:120, 1] <- U.init[ref:120]           
      # set initial detritivore spectrum  
      V[1:Ngrid, ref.det:120, 1] <- V.init[ref.det:120]  
      # set initial detritus biomass density (g.m^-3) 
      W[1:Ngrid, 1] <- W.init
    }
    
    if(use.init == TRUE){
     for (j in 1:Ngrid){
       # set up with the initial values from previous run
       #(phyto+zoo)plankton size spectrum 
       U[j, 1:(ref-1), 1] <- U.init[1:(ref-1)]    
       # set initial consumer size spectrum from previous run
       U[j, ref:length(x), 1] <- U.init[ref:length(x)]
       # set initial detritivore spectrum from previous run
       V[j, ref.det:length(x), 1] <- V.init[ref.det:length(x)]  
       W[j, 1] <- W.init
     }
    }
    
    #intrinsic natural mortality
    OM.u <- mu0*10^(-0.25*x)
    OM.v <- mu0*10^(-0.25*x)
    
    #senescence mortality rate to limit large fish from building up in the 
    #system
    #same function as in Law et al 2008, with chosen parameters gives similar 
    #M2 values as in Hall et al. 2006
    SM.u=k.sm*10^(p.s*(x-xs))
    SM.v=k.sm*10^(p.s*(x-xs))
    
    #Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!)
    
    # from Benoit & Rochet 2004 
    #Fvec[Fref:Nx] = 0.09*(x[Fref:Nx]/ log10(exp(1)) ) + 0.04 
    
    
    # here Fmort.u and Fmort.u= fixed catchability term for U and V to be 
    #estimated along with Fref.v and Fref.u
    
    # # CN check effort matrix 
    # dim(effort) # grid cell X time 
    # dim(Fvec.u) # grid cell X size X time 
    # # HERE selecting first time step - 
    # # does this mean starting values only and the time iteration is below? 
    
    Fvec.u[, Fref.u:Nx, 1] = Fmort.u*effort[,1]
    Fvec.v[, Fref.v:Nx, 1] = Fmort.v*effort[,1]
    
    # # CN check above 
    # Fvec.u[1,131,1] # OK this is catchability * effort 
    
    #output fisheries catches per yr at size
    
    # # CN check catch matrix 
    # dim(Y.u) # grid cell X size X time # here again selecting the first time step 
    # sum(Y.u[,Fref.u:Nx,1]) # 4408.997
    # sum(Y.v[,Fref.v:Nx,1]) # 94.88772 # so there is initla effort and initial catch 
    
    Y.u[, Fref.u:Nx, 1] <-
      Fvec.u[, Fref.u:Nx, 1]*U[, Fref.u:Nx, 1]*10^x[Fref.u:Nx] 
    #output fisheries catches per yr at size
    Y.v[, Fref.v:Nx, 1] <- 
      Fvec.v[, Fref.v:Nx, 1]*V[, Fref.v:Nx, 1]*10^x[Fref.v:Nx] 
    
    # # CN check - U and catches on different unit? 
    # U[1,131,1] # 0.875756
    # Y.u[1,131,1] # initial value 0.0005576951
    # Y.u[1,131,3]
    # Fvec.u[1,131,1] # initial value 6.368156e-05
    # Fvec.u[1,131,2] # 0 as not calcualted for next time step 

    #iteration over time, N [days]
    
    # Initial progress bar
    pb = txtProgressBar(min = 0, max = Neq, initial = 1, style = 3) 
    
    # time
    for (i in 1:(Neq)){ 
      if (i < Neq){
        # get proportion of total fishable biomass for each grid cell
        #output rates of fisheries catches per yr at size
        prop.b <- (apply(U[, Fref.u:Nx, i]*10^x[Fref.u:Nx]*dx, 1, sum) + 
                     apply(V[, Fref.v:Nx, i]*10^x[Fref.v:Nx]*dx, 1, sum)) / 
          sum(apply(U[, Fref.u:Nx, i]*10^x[Fref.u:Nx]*dx, 1, sum) + 
                apply(V[, Fref.v:Nx, i]*10^x[Fref.v:Nx]*dx, 1, sum))
        #check sums to 1
        #sum(prop.b)
        
        # redistribute total effort across gridcells according to proportion of 
        # biomass in that gridcell
        effort[, i+1] = gravitymodel(effort[, i+1], prop.b, depth = depth,
                                     iter = 1)
        # for option 2 iter = 1
        
      } # end gravity model 
      
      # grid cells
      for(j in 1:(Ngrid)) { 
        # Update progress bar
        setTxtProgressBar(pb, i) 
      
      # if(W[i]=="NaN"|W[i]<0)
      # {
      #   #browser()
      #   U[,i]<-"NaN"
      #   return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],
        #PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))
      # }
      
      #  ONLY IF RUNNING TO EQUILIBRIUM:
      #  below is to skip unecessary timesteps if model has reached equilibrium
      # if(i>100 & equilibrium==T )
      # {
      #   # browser()
      #   if(max(abs((log(U[-c(1:91),(i-1)]))-(log(U[-c(1:91),(i-2)]))))<eps
      #      &max(abs((log(V[-c(1:91),(i-1)]))-(log(V[-c(1:91),(i-2)]))))<eps
      #      &(abs((log(W[(i-1)]))-(log(W[(i-2)]))))<eps)
      #   {
      #     U[,Neq]<-U[,i-1]
      #     PM.v[,Neq]<-PM.v[,i-1]
      #     GG.v[,Neq]<-GG.v[,i-1]
      #     V[,Neq]<-V[,i-1]
      #     PM.u[,Neq]<-PM.u[,i-1]
      #     GG.u[,Neq]<-GG.u[,i-1]
      #     Y.u[,Neq]<-Y.u[,i-1]
      #     Y.v[,Neq]<-Y.v[,i-1]
      #     W[Neq]<-W[i-1]
      #     return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],
        #PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))
      #   }
      # }
      # 
      
        #--------------------------------
        # Calculate Growth and Mortality
        #--------------------------------
        if(temp.effect == T){
          pel.Tempeffect = exp(c1 - E/(Boltzmann*(sst[j, i]+ 273)))
          ben.Tempeffect = exp(c1 - E/(Boltzmann*(sft[j, i]+ 273)))
        }            
        
        if(temp.effect == F){
          pel.Tempeffect = 1
          ben.Tempeffect = 1
        }
        
        # feeding rates
        # yr-1
        f.pel <- pel.Tempeffect * 
          as.vector(((A.u*10^(x*alpha.u)*pref.pel[j])*(U[j,,i]*dx)%*%(gphi)) / 
                      (1+handling*(A.u*10^(x*alpha.u)*pref.pel[j]) *
                         (U[j, , i]*dx)%*%(gphi))) 
        # yr-1
        f.ben <- pel.Tempeffect * 
          as.vector(((A.u*10^(x*alpha.u)*pref.ben[j])*(V[j,,i]*dx)%*%(gphi)) / 
                      (1+handling*(A.u*10^(x*alpha.u)*pref.ben[j]) *
                         (V[j, , i]*dx)%*%(gphi))) 
        # yr-1 
        f.det <- ben.Tempeffect * 
          ((1/10^x)*A.v*10^(x*alpha.v)*W[j, i]) / 
          (1+handling*(1/10^x)*A.v*10^(x*alpha.v)*W[j, i]) 
        
        # Predator growth integral 
        #yr-1
        GG.u[j, , i] <- (1-def.high)*K.u*(f.pel)+(1-def.low)*K.v*(f.ben)
        
        # Reproduction
        #yr-1
        if(repro.on == 1){
          R.u[j, , i] <- (1-def.high)*(1-(K.u+AM.u))*(f.pel) + 
            (1-def.high)*(1-(K.v+AM.v))*f.ben
        }
      
          
        # Predator death integrals 
        #Satiation level of predator for pelagic prey
        sat.pel <- ifelse(f.pel > 0, f.pel / 
                            ((A.u*10^(x*alpha.u)*pref.pel[j]) * 
                               (U[j, , i]*dx)%*%(gphi)), 0)
        
        #yr-1 
        PM.u[j, , i] <- as.vector((pref.pel[j]*A.u*expax) * 
                                    (U[j, , i]*sat.pel*dx)%*%(mphi))
        
        #yr-1 
        #PM.u[,i]<-as.vector((1-f.pel)*(A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(mphi))
      
        #summing mortality terms - note no temperature effect on senescence 
        #mortality
        #yr-1
        Z.u[j, , i] <- PM.u[j, , i] + pel.Tempeffect*OM.u + SM.u + 
          Fvec.u[j, , i]
        
        # Benthos growth integral
        #yr-1
        GG.v[j, , i] <- (1-def.low)*K.d*f.det
        
        #reproduction
        #yr-1
        if(repro.on == 1){
          R.v[j, , i] <- (1-def.low)*(1-(K.d+AM.v))*(f.det)
        }
      
      # Benthos death integral
      #Satiation level of predator for benthic prey  
      sat.ben <- ifelse(f.ben > 0, 
                        f.ben/((A.u*10^(x*alpha.v)*pref.ben[j]) * 
                                 (V[j, , i]*dx)%*%(gphi)),
                        0)
      #yr-1
      PM.v[j, , i] <- ifelse(sat.ben > 0,
                             as.vector((pref.ben[j]*A.u*expax) * 
                                         (U[j, , i]*sat.ben*dx)%*%(mphi)),
                             0)
      
      #yr-1
      #PM.v[,i]<-as.vector((1-f.ben)*(A.u*10^(x*alpha.u)*pref.ben)*(U[,i]*dx)%*%(mphi))
      
      #yr-1
      Z.v[j, , i] <- PM.v[j, , i]+ben.Tempeffect*OM.v+SM.v+Fvec.v[j, , i]
      
      #total biomass density eaten by pred (g.m-2.yr-1)
      eatenbypred <- 10^x*f.pel*U[j, , i]+10^x*f.ben*U[j, , i] 
      
      #detritus output (g.m-2.yr-1)
      # losses from detritivore scavenging/filtering only:
      output.w <- sum(10^x*f.det*V[j, , i]*dx)   
      
      #total biomass density defecated by pred (g.m-2.yr-1)
      defbypred <- def.high*(f.pel)*10^x*U[j, , i] + 
        def.low*(f.ben)*10^x*U[j, , i]
      
      #------------------------------------------------
      # Increment values of W,U & V	for next time step  
      #------------------------------------------------
      
      #Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of 
      # detritus pool and solve for detritus biomass density in next time step 
      if(ERSEM.det.input == F){
        #considering pelagic faeces as input as well as dead bodies from both 
        #pelagic and benthic communities and phytodetritus (dying sinking 
        #phytoplankton)
        if(det.coupling == 1.0){
          # pelagic spectrum inputs (sinking dead bodies and faeces) - export 
          # ratio used for "sinking rate"
          #  + benthic spectrum inputs (dead stuff - already on/in seafloor)
          input.w <- (sinking.rate[j, i] * 
                        (sum(defbypred[ref:Nx]*dx) + 
                           sum(pel.Tempeffect*OM.u[1:Nx]*U[j, 1:Nx, i] *
                                 10^(x[1:Nx])*dx) + 
                           sum(pel.Tempeffect*SM.u[1:Nx]*U[j, 1:Nx, i] * 
                                 10^(x[1:Nx])*dx)) +
                        (sum(ben.Tempeffect*OM.v[1:Nx]*V[j, 1:Nx, i] * 
                               10^(x[1:Nx])*dx) + 
                           sum(ben.Tempeffect*SM.v[1:Nx]*V[j, 1:Nx, i] * 
                                 10^(x[1:Nx])*dx)))
          # )
        }
        
        if(det.coupling == 0.0){
          input.w <- sum(ben.Tempeffect*OM.v[1:Nx]*V[j, 1:Nx, i] * 
                           10^(x[1:Nx])*dx) + 
            sum(ben.Tempeffect*SM.v[1:Nx]*V[j, 1:Nx, i]*10^(x[1:Nx])*dx)
        }
        
        # get burial rate from Dunne et al. 2007 equation 3
        burial <- input.w*(0.013+0.53*input.w^2/(7+input.w)^2)
        
        # change in detritus biomass density (g.m-2.yr-1)                     
        
        # this one assumes immediate additional losses to sediment?
        # dW<-input.w - (output.w + W[i])     
        
        # losses due to detritivory only:
        # dW<-input.w - (output.w)     
        
        # losses from detritivory + burial rate (not including remineralisation 
        # bc that goes to p.p. after sediment, we are using realised p.p. as 
        # inputs to the model) 
        dW <- input.w-(output.w+burial)
        
        #biomass density of detritus g.m-2
        W[j, i+1] = W[j, i] + dW*delta_t  
      }
      
      if(ERSEM.det.input == T) {
         W[j, i+1] = W[j, i+1]
       }
      
      #----------------------------------------------
      # Pelagic Predator Density (nos.m-2) - solve for time + delta_t using 
      # implicit time Euler upwind finite difference (help from Ken Andersen 
      # and Richard Law)
      
      # Matrix setup for implicit differencing 
      Ai.u <- Bi.u <- Si.u <- array(0, c(length(x), 1))   
      Ai.u[idx] <- (1/log(10))*-GG.u[j, idx-1, i]*delta_t/dx
      Bi.u[idx] <- 1+(1/log(10))*GG.u[j, idx, i] * 
        delta_t/dx+Z.u[j, idx, i]*delta_t
      Si.u[idx] <- U[j, idx, i]
      
      # Boundary condition at upstream end 
      Ai.u[ref] <- 0
      Bi.u[ref] <- 1
      Si.u[ref] <- U[j, ref, i]
      
      # Invert matrix
      #recruitment at smallest consumer mass
      #continuation of plankton hold constant  
      # U[1:ref,i+1]<-U[1:ref,i] 
      
      if(use.init == T){
        U[j, 1:ref, i+1] <- ui0[j, i]*10^(r.plank[j, i]*x)[1:(ref)] 
      }
      
      if(use.init == F){
        U[j, 1:ref, i+1] <- ui0[j, i]*10^(r.plank[j, i]*x)[1:(ref)] 
      }
      
      # apply transfer efficiency of 10% *plankton density at same size  
      # if (repro.on==0)  U[ref,i+1]<-0.1*u.init.f(x,ui0,r.plank)[ref]       
      # reproduction from energy allocation
      if(repro.on == 1){
        U[j, ref, i+1] <- U[j, ref, i] +
          (sum(R.u[j, (ref+1):Nx, i]*10^x[(ref+1):Nx] * 
                 U[j, (ref+1):Nx, i]*dx)*delta_t)/(dx*10^x[ref]) - 
          (delta_t/dx)*(1/log(10))*(GG.u[j, ref, i])*U[j, ref, i] -
          delta_t*Z.u[j, ref, i]*U[j, ref, i]
      }
      
      #main loop calculation
      for(sz in (ref+1):(Nx)){
        #U[j,i+1]<-U[j,i]-(delta_t/dx)*(1/log(10))*(GG.u[j,i])*U[j,i] + 
        #(delta_t/dx)*(1/log(10))*GG.u[j-1,i]*U[j-1,i] -delta_t*Z.u[j,i]*U[j,i]
        U[j, sz, i+1] <- (Si.u[sz]-Ai.u[sz]*U[j, sz-1, i+1])/Bi.u[sz]
      }
      
      #----------------------------------------------
      # Benthic Detritivore Density (nos.m-2) 
      Ai.v <- Bi.v <- Si.v <- array(0, c(length(x), 1))   
      #shorthand for matrix referencing
      idx = (ref.det+1):Nx  
      
      Ai.v[idx] <- (1/log(10))*-GG.v[j, idx-1, i]*delta_t/dx 
      Bi.v[idx] <- 1+(1/log(10))*GG.v[j, idx, i]*delta_t/dx +
        Z.v[j, idx, i]*delta_t
      Si.v[idx] <- V[j, idx, i]
      
      #boundary condition at upstream end
      Ai.v[ref.det] <- 0
      Bi.v[ref.det] <- 1
      Si.v[ref.det] <- V[j, ref.det, i]  
      
      #invert matrix
      #recruitment at smallest detritivore mass
      #hold constant continution of plankton with sinking rate multiplier 
      V[j, 1:ref.det, i+1] <- V[j, 1:ref.det, i]
      
      # apply a very low of transfer efficiency 1%* total biomass of detritus 
      #divided by minimum size
      # if (repro.on==0)  V[ref.det,i+1]<-(0.01*W[i+1])/(10^x[ref]) 
      if(repro.on == 1){
        V[j, ref.det, i+1] <- V[j, ref.det, i] +
          sum(R.v[j, (ref.det+1):Nx, i]*10^x[(ref.det+1):Nx] * 
                V[j, (ref.det+1):Nx, i]*dx)*delta_t/(dx*10^x[ref.det]) - 
          (delta_t/dx)*(1/log(10))*(GG.v[j, ref.det, i])*V[j, ref.det, i] -
          delta_t*Z.v[j, ref.det, i]*V[j, ref.det, i]
      }
      
      #loop calculation
      for(sz in (ref.det+1):(Nx)){ 
        #V[j,i+1]<-V[j,i]-(delta_t/dx)*(1/log(10))*(GG.v[j,i])*V[j,i] + 
        #(delta_t/dx)*(1/log(10))*GG.v[j-1,i]*V[j-1,i]-delta_t*Z.v[j,i]*V[j,i] 
        V[j, sz, i+1] <- (Si.v[sz]-Ai.v[sz]*V[j, sz-1, i+1])/Bi.v[sz]
      }		
     
      # update Fmort 
      if (i < Neq){
        Fvec.u[j, Fref.u:Nx, i+1] = Fmort.u*effort[j, i+1]
        Fvec.v[j, Fref.v:Nx, i+1] = Fmort.v*effort[j, i+1]
      }

      #output rates of fisheries catches per yr at size
      Y.u[j, Fref.u:Nx, i+1] <- Fvec.u[j, Fref.u:Nx, i+1] *
        U[j, Fref.u:Nx, i+1]*10^x[Fref.u:Nx]
      #output rates of fisheries catches per yr at size
      Y.v[j, Fref.v:Nx, i+1] <- Fvec.v[j, Fref.v:Nx, i+1] *
        V[j, Fref.v:Nx, i+1]*10^x[Fref.v:Nx]
      
      ### ORIGINAL
      #output rates of fisheries catches per yr at size
      # Y.u[,Fref.u:Nx,i+1]<-Fvec.u[j,Fref.u:Nx,i]*U[j,Fref.u:Nx,i+1]*10^x[Fref.u:Nx]
      #output rates of fisheries catches per yr at size
      # Y.v[,Fref.v:Nx,i+1]<-Fvec.v[j,Fref.v:Nx,i]*V[j,Fref.v:Nx,i+1]*10^x[Fref.v:Nx]
    
      }#end across grid   
        
    }#end time iteration

    return(list(U = U[,,],
                GG.u = GG.u[,,],
                PM.u = PM.u[,,],
                V = V[,,],
                GG.v = GG.v[,,],
                PM.v = PM.v[,,],
                Y.u = Y.u[,,],
                Y.v = Y.v[,,],
                W = W[,], 
                params = params, 
                effort = effort[,]))
  
    # return(list(U=U[,Neq+1],GG.u=GG.u[,Neq],PM.u=PM.u[,Neq],V=V[,Neq+1],
      #GG.v=GG.v[,Neq],PM.v=PM.v[,Neq],Y=Y[,Neq],W=W[Neq+1], params=params))
    
  }) # end with(params)
  
}#end gridded size-based model function



#### model to run per grid cell or averaged over an area ------
sizemodel <- function(params, ERSEM.det.input = F, U_mat, V_mat, W_mat, 
                      temp.effect = T, eps = 1e-5, output="aggregated",
                      use.init = F, burnin.len){
  with(params,{
    #--------------------------------------------------------------------------
    # Model for a dynamical ecosystem comprised of: two functionally distinct 
    # size spectra (predators and detritivores), size structured primary 
    # producers and an unstructured detritus resource pool. 
    # time implicit upwind difference discretization algorithm (from Ken 
    # Andersen)
    # (Numerical Recipes and Richards notes, see "C://...standardised test 
    # implicit.txt" for basic example)
    # U represents abundance density (m-2) of "fish" and V abundance density 
    # (m-2) of "detritivores"
    # W is the pool of detritus (expressed in biomass density, g.m-2) - not 
    # size-based
    # Fish feed on other fish and benthic detritivores with size preference 
    # function 
    # Benthic detritivores feed on detritus produced from pelagic spectrum 
    # Senescence Mortality also included.  Options to include dynamic 
    # reproduction and predator handling time (but not currently used).
    #
    # Code modified for global fishing mortality rate application. 
    # JLB 17/02/2014
    # Code modified to include temperature scaling on senescence and detrital 
    # flux. RFH 18/06/2020
    # -------------------------------------------------------------------------
    
    # Input parameters to vary:
    # time series of intercept of plankton size spectrum (estimated from GCM, 
    # biogeophysical model output or satellite data).		
    ui0 <- 10^pp          
    # time series of slope of plankton size spectrum (estimated from GCM, 
    #biogeophysical model output or satellite data).  	
    r.plank <- r.plank
    # time series of temperature in water column
    sst <- sst
    #time series of temperature near seabed
    sft <- sft
    #time series of export ratio (read in sizeparam)
    sinking.rate <- sinking.rate
    
    #--------------------------------------------------------------------------
    # Functions called within sizemodel()
    #--------------------------------------------------------------------------
    
    #Function to build a lookup table for diet preference of all combinations of 
    #predator and prey body size: diet preference (in predator spectrum only)
    phi.f = function(q){
      q = q
      q0 = q0
      sd.q = sd.q
      phi = ifelse(q > 0,
                   exp(-(q-q0)*(q-q0)/(2*sd.q*sd.q))/(sd.q*sqrt(2.0*pi)),
                   0) 
      
      return(phi)
      
      # or normalise feeding kernel to sum to 1:
      # return(phi/sum(phi))
      # this is commented out because in previous work normalising did not 
      # produce realistic growth rates
    }
    
    #Functions to build lookup tables for components of integration which 
    #remain constant
    #growth
    gphi.f = function(q1){
      return(10^(-q1)*phi.f(q1))
    }
    
    #mortality
    mphi.f = function(q2){
      return(10^(alpha.u*q2)*phi.f(q2))
    }
    
    #Function to build lookup table for components of 10^(alpha*x)
    expax.f = function(x, alpha.u){
      return(10^(alpha.u*x))
    }
    
    #Function to compute convolution products: can be for growth or mortality 
    #depending on phi
    convolution.f = function(phi, u){
      res = matrix(NA, length(x), 1)
      for(i in 1:length(res)){
        res[i] = 0.0
        res[i] = sum(phi[, i]*u*dx)
      }
      return(res)
    }
    
    #--------------------------------
    # Growth and Mortality Equations
    #--------------------------------
    
    #predators
    #death.u<-function(u,A,expax,mphi){return(pref.pel*A*expax*
    #(convolution.f(mphi,u)))}
    #use faster matrix method instead of convolution loop function
    death.u <- function(u, A, expax, mphi){
      return((pref.pel*A*expax)*(u*dx)%*%(mphi))
    }
    
    #growth.u<-function(K0,K1,A,expax,gphi,u,v){return((pref.pel*K0*A*expax*
    #(convolution.f(gphi,u)))+ (pref.ben*K1*A*expax*(convolution.f(gphi,v))))}
    growth.u <- function(K0, K1, A, expax, gphi, u, v){
      return((pref.pel*K0*A*expax)*(u*dx)%*%(gphi)+
               (pref.ben*K1*A*expax)*(v*dx)%*%(gphi))
    }
    
    #detritivores
    #death.v<-function(u,A,expax,mphi){return((pref.ben)*A*expax*
    #(convolution.f(mphi,u)))}
    death.v <- function(u, A, expax, mphi){
      return((pref.ben*A*expax)*(u*dx)%*%(mphi))
    }
    growth.v <- function(K1, A, w){
      return((1/10^x)*(K1*A*10^(x*0.75)*w))
    }
    
    #detritus output (g.m-3.yr-1)
    out.w <- function(A, v, w){
      return(sum((A*(10^(x*0.75))*w*v)*dx))
      }     
    
    #--------------------------------------------------------------------------
    # Initialising matrices
    #--------------------------------------------------------------------------
    
    #q1 is a square matrix holding the log(predatorsize/preysize) for all 
    #combinations of sizes
    y = x
    
    q1  = matrix(NA, length(x), length(y))
    for(i in 1:length(y)){
      q1[, i] = y[i] - x
      }
    
    #q2 is the reverse matrix holding the log(preysize/predatorsize) for all 
    #combinations of sizes
    q2 = matrix(-q1, length(x), length(y))	
    
    #matrix for recording the two size spectra 
    V = U = array(0, c(length(x), Neq+1))
    
    #vector to hold detrtitus biomass density (g.m-2)
    W = array(0,Neq+1)
    
    #matrix for keeping track of growth and reproduction from ingested food:
    R.v = R.u = GG.v = GG.u = array(0, c(length(x), Neq+1)) 
    
    #matrix for keeping track of predation mortality
    PM.v = PM.u = array(0, c(length(x), Neq+1))   
    
    #matrix for keeping track of catches  
    Y.v = Y.u = array(0, c(length(x), Neq+1))  
    
    #matrix for keeping track of  total mortality (Z)
    Z.v = Z.u = array(0, c(length(x), Neq+1))
    
    #matrix for keeping track of senescence mortality and other (intrinsic) 
    #mortality
    SM.v = SM.u = OM.v = OM.u = array(0, length(x))
    
    #empty vector to hold fishing mortality rates at each size class at time
    Fvec.v = Fvec.u = array(0, c(length(x), Neq+1))
    
    #lookup tables for terms in the integrals which remain constant over time
    gphi  = gphi.f(q1)
    mphi  = mphi.f(q2)
    
    #lookup table for components of 10^(alpha.u*x)
    expax = expax.f(x, alpha.u)
    
    #--------------------------------------------------------------------------
    # Numerical integration
    #--------------------------------------------------------------------------
    
    # set up with the initial values from param
    #(phyto+zoo)plankton size spectrum  
    U[1:(ref-1), 1] <- U.init[1:(ref-1)]    
    # set initial consumer size spectrum 
    U[ref:120, 1] <- U.init[ref:120]       
    # set initial detritivore spectrum  
    V[ref.det:120, 1] <- V.init[ref.det:120]  
    # set initial detritus biomass density (g.m^-3) 
    W[1] <- W.init
    
    if(use.init == TRUE){
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
    OM.u <- mu0*10^(-0.25*x)
    OM.v <- mu0*10^(-0.25*x)
    
    #senescence mortality rate to limit large fish from building up in the 
    #system
    #same function as in Law et al 2008, with chosen parameters gives similar 
    #M2 values as in Hall et al. 2006
    SM.u = k.sm*10^(p.s*(x-xs))
    SM.v = k.sm*10^(p.s*(x-xs))
    
    #Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!)
    
    #Fvec[Fref:Nx] = 0.09*(x[Fref:Nx]/ log10(exp(1)) ) + 0.04 
    # from Benoit & Rochet 2004 
     
    # here Fmort.u and Fmort.u= fixed catchability term for U and V to be
    #estimated along with Fref.v and Fref.u
    Fvec.u[Fref.u:Nx, 1] = Fmort.u*effort[1]
    Fvec.v[Fref.v:Nx, 1] = Fmort.v*effort[1]
    
    #output fisheries catches per yr at size
    Y.u[Fref.u:Nx, 1] <- Fvec.u[Fref.u:Nx, 1]*U[Fref.u:Nx, 1]*10^x[Fref.u:Nx] 
    #output fisheries catches per yr at size
    Y.v[Fref.v:Nx, 1] <- Fvec.v[Fref.v:Nx, 1]*V[Fref.v:Nx, 1]*10^x[Fref.v:Nx] 
    
    #iteration over time, N [days]
    # Initial progress bar
    #pb = txtProgressBar(min = 0, max = Neq, initial = 1, style = 3) 
    
    for(i in 1:(Neq)){
      #  setTxtProgressBar(pb, i) # Update progress bar
      # if(W[i]=="NaN"|W[i]<0)
      # {
      #   #browser()
      #   U[,i]<-"NaN"
      #   return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],
      #PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))
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
       #     PM.v[,Neq]<-PM.v[,i-1]
       #     GG.v[,Neq]<-GG.v[,i-1]
       #     V[,Neq]<-V[,i-1]
       #     PM.u[,Neq]<-PM.u[,i-1]
       #     GG.u[,Neq]<-GG.u[,i-1]
       #     Y.u[,Neq]<-Y.u[,i-1]
       #     Y.v[,Neq]<-Y.v[,i-1]
       #     W[Neq]<-W[i-1]
       #     return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],
      #PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))
       #   }
       # }
       # 

      
      #--------------------------------
      # Calculate Growth and Mortality
      #--------------------------------
      if(temp.effect == T){
        pel.Tempeffect = exp(c1-E/(Boltzmann*(sst+273)))
        ben.Tempeffect = exp(c1-E/(Boltzmann*(sft+273)))
      }
      if(temp.effect == F){
        pel.Tempeffect = 1
        ben.Tempeffect = 1
      }
      
      # feeding rates
      # yr-1
      f.pel <- pel.Tempeffect[i] * 
        as.vector(((A.u*10^(x*alpha.u)*pref.pel)*(U[, i]*dx)%*%(gphi)) / 
                    (1+handling*(A.u*10^(x*alpha.u)*pref.pel) * 
                       (U[, i]*dx)%*%(gphi))) 
      # yr-1
      f.ben <- pel.Tempeffect[i] * 
        as.vector(((A.u*10^(x*alpha.u)*pref.ben)*(V[, i]*dx)%*%(gphi)) / 
                    (1+handling*(A.u*10^(x*alpha.u)*pref.ben) *
                       (V[, i]*dx)%*%(gphi)))
      # yr-1
      f.det <- ben.Tempeffect[i] *
        ((1/10^x)*A.v*10^(x*alpha.v)*W[i]) / 
        (1+handling*(1/10^x)*A.v*10^(x*alpha.v)*W[i])
      
      # Predator growth integral 
      # yr-1
      GG.u[, i] <- (1-def.high)*K.u*(f.pel)+(1-def.low)*K.v*(f.ben)
      
      # Reproduction
      # yr-1
      if(repro.on == 1){
        R.u[, i] <- (1-def.high)*(1-(K.u+AM.u))*(f.pel) +
          (1-def.high)*(1-(K.v+AM.v))*f.ben
      }
      
      # Predator death integrals 
      #Satiation level of predator for pelagic prey
      sat.pel <- ifelse(f.pel > 0,
                        f.pel/((A.u*10^(x*alpha.u)*pref.pel) * 
                                 (U[,i]*dx)%*%(gphi)),
                        0)
      # yr-1
      PM.u[, i] <- as.vector((pref.pel*A.u*expax)*(U[, i]*sat.pel*dx)%*%(mphi))
      # yr-1
      #PM.u[,i]<-as.vector((1-f.pel)*(A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(mphi))
      
      # yr-1
      Z.u[, i] <- PM.u[, i]+pel.Tempeffect[i]*OM.u+SM.u+Fvec.u[, i]
      
      # Benthos growth integral
      # yr-1
      GG.v[, i] <- (1-def.low)*K.d*f.det
      
      #reproduction
      # yr-1
      if(repro.on == 1){
        R.v[, i] <- (1-def.low)*(1-(K.d+AM.v))*(f.det)
      }
      
      # Benthos death integral
      #Satiation level of predator for benthic prey 
      sat.ben <- ifelse(f.ben > 0,
                        f.ben/((A.u*10^(x*alpha.v)*pref.ben) * 
                                 (V[, i]*dx)%*%(gphi)),
                        0)
      # yr-1
      PM.v[, i] <- ifelse(sat.ben > 0,
                          as.vector((pref.ben*A.u*expax) *
                                      (U[, i]*sat.ben*dx)%*%(mphi)),
                          0)
      # yr-1
      #PM.v[,i]<-as.vector((1-f.ben)*(A.u*10^(x*alpha.u)*pref.ben)*(U[,i]*dx)%*%(mphi))
      
      # yr-1
      Z.v[, i] <- PM.v[, i]+ben.Tempeffect[i]*OM.v+SM.v+Fvec.v[, i]
      
      #total biomass density eaten by pred (g.m-2.yr-1)
      eatenbypred <- 10^x*f.pel*U[, i]+10^x*f.ben*U[, i] 
      
      #detritus output (g.m-2.yr-1)
      # losses from detritivore scavenging/filtering only:
      output.w <- sum(10^x*f.det*V[, i]*dx)   
      
      #total biomass density defecated by pred (g.m-2.yr-1)
      defbypred <- def.high*(f.pel)*10^x*U[, i]+ def.low*(f.ben)*10^x*U[, i]
      
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
          input.w <- (sinking.rate[i] * 
                        (sum(defbypred[ref:Nx]*dx) +
                           sum(pel.Tempeffect[i] * 
                                 OM.u[1:Nx]*U[1:Nx, i]*10^(x[1:Nx])*dx) + 
                           sum(pel.Tempeffect[i] * 
                                 SM.u[1:Nx]*U[1:Nx,i]*10^(x[1:Nx])*dx)) +
                        (sum(ben.Tempeffect[i] * 
                               OM.v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx) + 
                           sum(ben.Tempeffect[i] * 
                                 SM.v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx)))
          # )
        }
        if(det.coupling == 0.0){
          input.w <- sum(ben.Tempeffect[i] *
                           OM.v[1:Nx]*V[1:Nx, i]*10^(x[1:Nx])*dx) + 
            sum(ben.Tempeffect[i]*SM.v[1:Nx]*V[1:Nx, i]*10^(x[1:Nx])*dx)
        }
        
        # get burial rate from Dunne et al. 2007 equation 3
        burial <- input.w*(0.013 + 0.53*input.w^2/(7+input.w)^2)
        
        # change in detritus biomass density (g.m-2.yr-1)                     
        # this one assumes immeidate additional losses to sediment?
        # dW<-input.w - (output.w + W[i])     
        # losses due to detritivory only:
        # dW<-input.w - (output.w)     
        
        # losses from detritivory + burial rate (not including remineralisation
        # bc that goes to p.p. after sediment, we are using realised p.p. as
        # inputs to the model) 
        dW <- input.w-(output.w+burial) 
        #biomass density of detritus g.m-2
        W[i+1] = W[i]+dW*delta_t
      }
      
      if(ERSEM.det.input == T){
        W[i+1] = W[i]
      }
      
      #----------------------------------------------
      # Pelagic Predator Density (nos.m-2)- solve for time + delta_t using
      # implicit time Euler upwind finite difference (help from Ken Andersen 
      # and Richard Law)
      
      # Matrix setup for implict differencing 
      Ai.u <- Bi.u <- Si.u <- array(0, c(length(x), 1))   
      
      Ai.u[idx] <- (1/log(10))*-GG.u[idx-1, i]*delta_t/dx
      Bi.u[idx] <- 1+(1/log(10))*GG.u[idx, i]*delta_t/dx +Z.u[idx, i]*delta_t
      Si.u[idx] <- U[idx, i]
      
      # Boundary condition at upstream end 
      Ai.u[ref] <- 0
      Bi.u[ref] <- 1
      Si.u[ref] <- U[ref, i]
      
      # Invert matrix
      #recruitment at smallest consumer mass
      #continuation of plankton hold constant  
      # U[1:ref,i+1]<-U[1:ref,i] 
      
      if(use.init==T){
        U[1:ref, i+1] <- ui0[i]*10^(r.plank[i]*x)[1:(ref)] 
      }
      if(use.init == F){
        U[1:ref, i+1] <- ui0[i+1]*10^(r.plank[i+1]*x)[1:(ref)] 
      }
      
      # apply transfer efficency of 10% *plankton density at same size  
      # if (repro.on==0)  U[ref,i+1]<-0.1*u.init.f(x,ui0,r.plank)[ref]       
      # reproduction from energy allocation
      if(repro.on == 1){
        U[ref, i+1] <- U[ref, i] +
          (sum(R.u[(ref+1):Nx, i]*10^x[(ref+1):Nx]*U[(ref+1):Nx, i]*dx) *
             delta_t)/(dx*10^x[ref])-(delta_t/dx)*(1/log(10)) *
          (GG.u[ref, i])*U[ref, i]-delta_t*Z.u[ref, i]*U[ref, i]
      }
      
      #main loop calculation
      for(j in (ref+1):(Nx)){
        #U[j,i+1]<-U[j,i]-(delta_t/dx)*(1/log(10))*(GG.u[j,i])*U[j,i] + 
        #(delta_t/dx)*(1/log(10))*GG.u[j-1,i]*U[j-1,i] -delta_t*Z.u[j,i]*U[j,i]
        
        U[j, i+1] <- (Si.u[j]-Ai.u[j]*U[j-1, i+1])/Bi.u[j]
      }
      
      #----------------------------------------------
      # Benthic Detritivore Density (nos.m-2) 
      Ai.v <- Bi.v <- Si.v <- array(0, c(length(x), 1))   
      #shorthand for matrix referencing
      idx = (ref.det+1):Nx  
      
      Ai.v[idx] <- (1/log(10))*-GG.v[idx-1,i]*delta_t/dx 
      Bi.v[idx] <- 1+(1/log(10))*GG.v[idx, i]*delta_t/dx +Z.v[idx, i]*delta_t
      Si.v[idx] <- V[idx, i]
      
      #boundary condition at upstream end
      Ai.v[ref.det] <- 0
      Bi.v[ref.det] <- 1
      Si.v[ref.det] <- V[ref.det, i]  
      
      #invert matrix
      #recruitment at smallest detritivore mass  
      #hold constant continution of plankton with sinking rate multiplier 
      V[1:ref.det, i+1] <- V[1:ref.det, i]
      
      # apply a very low of transfer effiency 1%* total biomass of detritus
      #divided by minimum size
      # if (repro.on==0)  V[ref.det,i+1]<-(0.01*W[i+1])/(10^x[ref]) 
      if(repro.on == 1){
        V[ref.det, i+1] <- V[ref.det, i] +
          sum(R.v[(ref.det+1):Nx, i]*10^x[(ref.det+1):Nx] * 
                V[(ref.det+1):Nx, i]*dx)*delta_t/(dx*10^x[ref.det]) - 
          (delta_t/dx)*(1/log(10))*(GG.v[ref.det, i])*V[ref.det, i] - 
          delta_t*Z.v[ref.det, i]*V[ref.det, i]
      }
      
      #loop calculation
      for(j in (ref.det+1):(Nx)){ 
        #V[j,i+1]<-V[j,i]-(delta_t/dx)*(1/log(10))*(GG.v[j,i])*V[j,i] + 
        #(delta_t/dx)*(1/log(10))*GG.v[j-1,i]*V[j-1,i]-delta_t*Z.v[j,i]*V[j,i] 
        V[j, i+1] <- (Si.v[j]-Ai.v[j]*V[j-1, i+1])/Bi.v[j]
      }		
      rm(j)
      
      # increment fishing 
      Fvec.u[Fref.u:Nx, i+1] = Fmort.u*effort[i+1]
      Fvec.v[Fref.v:Nx, i+1] = Fmort.v*effort[i+1]
      
      #output fisheries catches per yr at size
      Y.u[Fref.u:Nx, i+1] <- Fvec.u[Fref.u:Nx, i+1]*U[Fref.u:Nx, i+1] *
        10^x[Fref.u:Nx] 
      #output fisheries catches per yr at size
      Y.v[Fref.v:Nx, i+1] <- Fvec.v[Fref.v:Nx, i+1]*V[Fref.v:Nx, i+1] *
        10^x[Fref.v:Nx] 
    }
    #end time iteration
    
    return(list(U = U[,],
                GG.u = GG.u[,],
                PM.u = PM.u[,],
                V = V[,], 
                GG.v = GG.v[,],
                PM.v = PM.v[,],
                Y.u = Y.u[,],
                Y.v = Y.v[,],
                W = W[],
                params=params))
    # return(list(U=U[,Neq+1],GG.u=GG.u[,Neq],PM.u=PM.u[,Neq],V=V[,Neq+1],
    #GG.v=GG.v[,Neq],PM.v=PM.v[,Neq],Y=Y[,Neq],W=W[Neq+1], params=params))
  })  
  # end with(params)
}
#end size-based model function


sizeparam <- function(equilibrium = F, dx = 0.1, xmin = -12, xmax = 6,
                      xmin.consumer.u = -7, xmin.consumer.v = -7, tmax = 100,
                      tstepspryr = 48, fmort.u = 0.0, fminx.u = 1, 
                      fmort.v = 0.0, fminx.v = 1, er = 0.5, pp = -3, 
                      slope = -1, lat = NA, lon = NA, depth = 500, sst = 20,
                      sft = 20, effort = 1.0, search_vol = 640, use.init = F, 
                      U.initial = NA, V.initial = NA, W.initial = NA, 
                      Ngrid = NA){
  #----------------------------------------------------------------------------
  # FUNCTION TO GET Parameters of model
  #----------------------------------------------------------------------------
  
  #----------------------------------------------------------------------------
  # Parameters:
  #----------------------------------------------------------------------------
  param = list()
  
  # store grid info: lat, lon and depth 
  param$lat = lat
  param$lon = lon
  param$depth = depth
  
  # plankton parameters   
  # intercept of phyto-zooplankton spectrum
  param$pp = pp
  # slope of phyto-zooplankton spectrum
  param$r.plank = slope
  
  # export ratio
  param$er = er
  
  # temperature parameters 
  # sea-surface temperature - degrees celsius
  param$sst = sst
  # near sea-floor temperature - degrees celsius
  param$sft = sft
  
  # fishing parameters 
  # fishing mortality rate for predators 
  param$Fmort.u = fmort.u
  # fishing mortality rate for detritivores
  param$Fmort.v = fmort.v
  # minimum log10 body size fished for predators
  param$min.fishing.size.u = fminx.u 
  # minimum log10 body size fished for detritivores
  param$min.fishing.size.v = fminx.v 

  # get effort and rescale 
  param$effort = effort
  
  # get effort and rescale so that max is set to 1
  #param$effort = effort/max(effort)
  
  # benthic-pelagic coupling parameters
  # set predator coupling to benthos, depth dependent - 0.75 above 500 m, 0.5
  # between 500-1800 and 0 below 	1800m (suggestions of values from Clive
  # Trueman based on stable isotope work, and proportion of biomass, 	Rockall 
  # Trough studies)
  
  param$pref.ben = 0.8*exp(-1/250*depth)
  
  # preference for pelagic prey 
  param$pref.pel = 1-param$pref.ben 
  # detritus coupling on? 1= yes, 0 = no
  param$det.coupling = 1.0  
  # fraction of sinking detritus reaching the seafloor (from export ratio input)
  param$sinking.rate = param$er 
  # fraction of sinking detritus reaching the seafloor (SAME AS pref.ben)
  # param$sinking.rate = (0.8*exp(-1/250*depth))*0.01 
  
  # feeding and energy budget parameters
  # Mean log10 predator prey mass ratio  100:1.
  param$q0 = 2.0
  # 0.6 to 1.0 for lognormal prey preference function. 
  param$sd.q = 1.0
  # originally 640, but if 64 then using Quest-fish default of 64 hourly rate
  # volume searched constant m3.yr-1 for fish. need to check this value, its
  # quite large.
  param$A.u = search_vol    
  # hourly rate volume filtered constant m3.yr-1 for benthos. this value yields 
  # believable growth curve. Approximatelty 10 times less than A.u
  param$A.v = param$A.u*0.1    
  #exponent for metabolic requirements plus swimming for predators (Ware et al 
  # 1978)
  param$alpha.u = 0.82  
  #exponent =0.75 for whole organism basal (sedentary) metabolic rate 
  #(see growth.v) from Peters (1983) and Brown et al. (2004) for detritivores
  param$alpha.v = 0.75  
  # fraction of ingested food that is defecated (Peters,1983)
  param$def.high = 0.3  
  # low = low quality (K) food, high = high quality (K) food
  param$def.low = 0.5   
  #net growth conversion efficiency for organisms in the "predator" spectrum
  #from Ware (1978)
  param$K.u = 0.3       
  #net growth conversion efficiency for organisms in the "detritivore" spectrum
  param$K.v = 0.2
  #net growth conversion efficency for detritus
  param$K.d = param$K.v
  #fraction of energy required for maintenance & activity etc.
  param$AM.u = 0.5
  param$AM.v = 0.7
  # if handling time is > 0 Type II functional response, if = 0 linear (no 
  # predator satiation) - 5.7e-7
  param$handling = 0
  # dynamic reproduction on=1,off=0  
  param$repro.on = 1
  # constant used in Jennings et al. 2008 Proc B to standardize 
  #metabolism-temperature effects
  # for Boltzmann equation. Derived from Simon's fit to Andy Clarke's data
  param$c1 = 25.22 
  # activation energy, eV
  param$E = 0.63
  # Boltzmann's constant
  param$Boltzmann = 8.62*10^-5
  
  # "other" mortality parameters
  # residual natural mortality
  param$mu0 = 0.2
  # size at sensenscence
  param$xs = 3
  # exponent for senescence mortality 
  param$p.s = 0.3
  # constant for senescence mortality 
  param$k.sm = 0.2
  
  #----------------------------------------------------------------------------
  # Parameters for numerical integration (size & time discretisation)
  #----------------------------------------------------------------------------
  # size increment after discretization for integration (log body weight)
  param$dx =  dx	  
  # minimum log10 body size of plankton
  param$xmin = xmin  
  # minimum log10 body size in dynamics predators
  param$x1 = xmin.consumer.u
  # minimum log10 body size in dynamics benthic detritivores
  param$x1.det = xmin.consumer.v
  # maximum log10 body size of predators
  param$xmax = xmax
  # maximum log10 body size before senescence kicks in (departure form 
  # linearity)
  param$xmax2 = 4.0
  
  #Vector with size bins
  param$x = seq(param$xmin, param$xmax, param$dx)
  param$Nx = length(param$x)
  
  #position in x of x1
  param$ref = ((param$x1-param$xmin)/param$dx)+1
  
  #position in x of x1.det
  param$ref.det = ((param$x1.det-param$xmin)/param$dx)+1 
  
  #position in F vector corresponding to smallest size fished in U
  param$Fref.u = ((param$min.fishing.size.u-param$xmin)/param$dx)+1
  #position in F vector corresponding to smallest size fished in V
  param$Fref.v = ((param$min.fishing.size.v-param$xmin)/param$dx)+1
  
  #short hand for matrix indexing
  param$idx = 2:param$Nx
  
  # integration parameters 
  # number of years to run model
  param$tmax = tmax
  # discretisation of year 
  param$delta_t = (1/tstepspryr)
  # number of time bins 
  param$Neq = param$tmax/param$delta_t	
  # number of grid cells
  param$Ngrid = Ngrid	
  
  # Set initial values for size spectra & detritus
  # param$U.init=read.table("pp_0.006F_0pel_0.5ben_0.5sink_0.5_U.txt")
  # param$V.init=read.table("pp_0.006F_0pel_0.5ben_0.5sink_0.5_V.txt")
  # param$W.init=scan("pp_0.006F_0pel_0.5ben_0.5sink_0.5_Det.txt",quiet=TRUE)
  
  # (phyto+zoo)plankton + pelagic predator size spectrum
  param$U.init <- 10^param$pp[1]*10^(param$r.plank[1]*param$x)  
  # set initial detritivore spectrum  
  param$V.init <- param$sinking.rate[1]*10^param$pp[1] * 
    10^(param$r.plank[1]*param$x)
  # abritrary initial value for detritus
  param$W.init <- 0.00001
  
  if(use.init == TRUE){
    param$U.init <- U.initial
    param$V.init <- V.initial
    param$W.init <- W.initial
  }
  
  param$equilibrium = equilibrium
  
  return(param)
}
#end sizeparam function


plotsizespectrum <- function(modeloutput, params, itime = params$Neq, 
                             timeaveraged = F){
  with(params, {
    ###############################PLOT######################################
    # plot changes in the two size spectra over time
    ui0 <- 10^pp[Neq]
    
    U <- modeloutput$U[, itime]
    V <- modeloutput$V[, itime]
    W <- modeloutput$W[itime]
    
    if(timeaveraged == T){
      U <- rowMeans(U)
      V <- rowMeans(V)
      W <- mean(W)
    }
    
    maxy = max(log10(U))
    miny = -20
    
    # par(mfrow=c(1,1))
    plot(x[ref:Nx], log10(U[ref:Nx]), type = "l", col = "blue", cex = 1.6,
         ylab = c("Log abundance density [1/m3]"), 
         xlab = c("log Body mass [g]"), 
         xlim = c(x1.det, xmax), ylim = c(miny, maxy))
    
    points(x[ref.det:Nx], log10(V[ref.det:Nx]), type = "l", col = "red", 
           cex = 1.6, ylab = c(""), xlab = c(""))
    
    # text(x1.det+abs(0.05*x1.det),miny+abs(1.3*miny),
    #paste("day = ",itime,sep=""),pos=4)
    
    text(x1.det+abs(0.05*x1.det), miny+abs(0.3*miny),
         paste("pel.pref=", round(pref.pel, 2), sep = ""), pos = 4)
    
    text(x1.det+abs(0.05*x1.det), miny+abs(0.2*miny), 
         paste("ben.pref=", round(pref.ben, 2), sep = ""), pos = 4)
    
    #text(x1.det+abs(0.05*x1.det),miny+abs(0.1*miny),
    #paste("PP=",round(ui0,3),sep=""),pos=4)
    
    legend(xmax-0.48*(xmax-x1.det), maxy-0.05*maxy, 
           c("Predators", "Detritivores"), col = c("blue","red"), lwd = 2.0)
  })
  # end with(params) 
}
# end plot function


# Method to calculate intercept from appendix of Barnes et al. 2010 JPR -----
GetPPIntSlope <- function(sphy, lphy, mmin = 10^-14.25, mmid = 10^-10.184,
                          mmax = 10^-5.25, depth, output = "slope"){
  # need to convert sphy and lphy from mol C / m^3 to g C / m^3
  sphy = (sphy*12.0107)
  lphy = (lphy*12.0107)
  
  #  if it's depth integrated units are /m^-2 and need to divide my depth if 
  # using depth integrated inputs   
  ## CN in ageement with JB - remove depth integration of inputs  
  # sphy= sphy/min(depth,100)
  # lphy= lphy/min(depth,100)
   
  #mmin=10^-14.25 gww, approx 0.2 ESD
  #division between small and large phytoplankton in GFDL-Topaz
  #mmid=10^-10.184 gww, approax 5 ESD
  #mmax=10^-5.25 gww, approx 200 ESD
  
  
  #from Appendix of Barnes 2010 JPR, used in Woodworth-Jefcoats et al. 2013 GCB
  #the scaling of biomass with body mass can be described as B=aM^b
  
  # the exponent b (also equivalent to the slope b in a log B vs log M 
  # relationship) can be assumed:
  # 0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
  # most studies seem to suggest N~M^-1, so can assume that and test 
  # sensitivity of our results to this assumption. 
  
  #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
  #in log10 gww
  midsmall <- log10((mmin+mmid)/2) 
  #in log10 gww
  midlarge <- log10((mmid+mmax)/2)
  
  #convert to log10 (gww/size class median size) for log10 abundance
  small <- log10((sphy*10)/10^midsmall)
  #convert to log10 (gww/size class median size) for log10 abundance
  large <- log10((lphy*10)/10^midlarge)
  
  b <- (small-large)/(midsmall-midlarge)
  
  #a is really log10(a), same a when small, midsmall are used
  a <- large - b*midlarge  
  # log10(a) will equal to the log10 density of particles at mass = 1 g and at
  # log10 (mass)=0
  # pp<-c(a,b) 
  # a could be used directly to replace 10^pp in sizemodel()
  if(output =="slope"){
    return(b)
  }
  if(output =="intercept"){
    return(a)
  }
}


getExportRatio <- function(sphy, lphy, sst, depth){
  ptotal <- sphy+lphy
  plarge <- lphy/ptotal
  psmall <- sphy/ptotal
  er <- (exp(-0.032*sst)*((0.14*psmall)+(0.74*(plarge))) + 
           (0.0228*(plarge)*(depth*0.004)))/(1+(depth*0.004)) 
  return(er)
}

# Added 1 May, 2013
# Methods to calculate:

# 1. primary producer size distribution (as described by MB10, 50, 90)
# 2. primary prodcuer slope and intercept (from Barnes et al. and Simon's notes)
# 3.detritus input dependent on primary production and depth - proportion of 
# sinking phytoplankton for detritus input
# 4. minimum consumer size
# 5. input for preferred predator-prey mass ratio function - beta and sigma

# These are calculated from measurements of temperature, primary production 
GetPPSpectrum <- function(pp = 0.77, sst = 20, chla = NA){
  # List of parameters
  # sst= temperature (Celsius) - INPUT, in degrees C
  # P= primary production Note: needs to be in g C m-2 d-1  - INPUT
  # chla = chl a biomass Note: mg Chl m-3, if used instead of pp input !=NA
  
  # OUTPUTS:
  # a= intercept phytoplankton size spectrum (orginal units are (x) cell mass 
  #pg C and (y) log10 pg C m-3 converted to g vs. number density per m^3)
  # b= slope phytoplankton size spectrum 
  # mb50= as Barnes et al
  
  # /* phytoplankton size spectrum parameters */
  log10pp <- log10(pp*1000)
  
  # predict log 10 (pg C) median phytoplankton cell size
  mb50log <- (-0.081*sst)+(0.875*log10pp)-0.661
  
  # predict phytoplankton ss intercept from PP as g C m-2 d-1
  # spectrum units are log10 pg C m-3 (y) v cell mass log10 pg C (x)
  mbint <- (0.554*log10pp)+8.087
  
  # predict phytoplankton ss slope from PP as g C m-2 d-1 
  # spectrum units are log10 pg C m-3 (y) v cell mass log10 pg C (x)
  mbslope <- (0.182*log10pp)-1.715
  
  #fixed log predator=prey mass ratio to calculate minimum consumer size in 
  #the model
  beta <- 2.0
  
  # need to convert log 10 pg C to g C (10^-12) then g C to g wet weight 
  #(10^-1; 1 g C = 10 g ww)
  mb50a <- 10^(mb50log)*10^-12*10^-1
  
  # mb50 now in log 10 (grams wet weight)
  mb50 <- log10(mb50a)
  
  # limits to prediction
  if(any(mbslope > -1.05)){
    mbslope = -1.05
  }  
  if(any(mbslope < -1.45)){
    mbslope = -1.45
  }
  
  return(list(
    # OUTPUTS: 
    # a is intercept - log 10 (pg C m^-3) - convert to log10 (g ww m^3) -> 
    #log10 (density m^3) 
    a = log10((10^mbint)*10^-13),
    # pp slope
    b = mbslope,
    # log 10 (median pp size, g ww),
    med.x = mb50,
    # estimate min consumer size from PPMR and median phytoplankton size 
    #(see also Woodworth et al 2012)
    min.cons.x = log10(10^beta*10^mb50)
  ))
}


optimizeYield <- function(params, com, pp = pp, sst = sst, sft = sft, 
                          U_mat = U, V_mat = V, sizeindex = (1:params$Nx)){
  #----------------------------------------------------------------------------
  #
  # Function to Optimise yield , based on a single parameter, (F-even, a single
  # F across all sizes) 
  # 
  #----------------------------------------------------------------------------
  
  ## ------------------------------------------------------------------
  ##  Do a single-parameter optimization of  yield
  ## ------------------------------------------------------------------
  
  # Helper function to return the quantity to be optimized based on the
  # fishing mortality:
  totalYield <- function(newF, params, com, sizeindex){
    print(newF)
    params$Fmort = newF
    params$tstepdays = 1.0
    #change to 10 years( not running all the way to equilib within this)
    params$tmaxyears = 10   
    # total number of steps in a period
    params$Neq = as.integer(365.0*params$tmaxyears/params$tstepdays)	
    # time step (rates per year)  
    params$delta_t = params$tmaxyears/params$Neq
    params$siz_time = 1:params$Neq
    
    params$U.init[, 1] = com$U
    params$V.init[, 1] = com$V
    params$W.init[1] = com$W 
    
    com <- sizemodel(params = params, pp, sst, sft, U_mat, V_mat, 
                     temp.effect = T) 
    
    # mean yield through time for specific size range (or entire fished range)    
    yield <- com$Y[sizeindex]
    totalYield <- sum(yield)*1e-6
    -totalYield
  }
  
  result <- optimize(f = totalYield, interval = c(0, 1.5), maximum = F,
                     params = params, com = com)
  
  # result <- nlm(f=totalYield,p=0.8,params=params, com=com)
  # result <- genoud(fn=totalYield, nvars=1,Domains=cbind(0,2), max=FALSE,
  #params=params, com=com)
  # result<-GenSA(par=0.1, lower=0, upper=3, fn=totalYield,params=params,
  #com=com,control=list(maxit=100))
  
  return(list(result, com))
}

# # get proportion of total fishable biomass for each grid cell
# #output rates of fisheries catches per yr at size
# prop.b<-(apply(U[,Fref.u:Nx,i]*10^x[Fref.u:Nx]*dx,1,sum) + 
# apply(V[,Fref.v:Nx,i]*10^x[Fref.v:Nx]*dx,1,sum))/
# sum(apply(U[,Fref.u:Nx,i]*10^x[Fref.u:Nx]*dx,1,sum) +
# apply(V[,Fref.v:Nx,i]*10^x[Fref.v:Nx]*dx,1,sum))
# #check sums to 1
# #sum(prop.b)
# dim(effort)
# plot(colSums(effort_grid[,-1]))
# depth

gravitymodel <- function(effort = effort[, 3000], prop.b, depth, iter){
  # redistribute total effort across gridcells according to proportion of
  # biomass in that gridcell
  # using graivity model, Walters & Bonfil, 1999, Gelchu & Pauly 2007
  # ideal free distribution - Blanchard et al 2008
  # new_effort= prop.b*effort
  eff <- as.vector(effort)
  d <- unlist(as.vector(depth))
  
  for(j in 1:iter){
    a <- eff
    # suit = prop.b*(1-d/max(d))*(1-a/0.001)
    # # CN option 2 from below 
    suit <- prop.b*(1-d/max(d))
    # # CN option 3 from below 
    # suit = prop.b
    # rescale:
    rel.suit <- suit/sum(suit)
    neweffort <- eff+rel.suit*eff
    mult <- sum(eff)/sum(neweffort)
    #gradient drives individuals to best locations at equilibrium (assumed to 
    #reached after 10000 iterations)
    eff <- neweffort*mult
  }
  
 # plot(prop.b,effort_grid[,3200],cex=(1-d/max(d)))
 # plot(prop.b,eff,cex=(1-d/max(d)))
 # 
 #  rmax<-max(rel.suit)
 #  r<-ifelse(suit >= rmax,rmax,s1)
 #  newden<-ifelse(s1 >= rmax,density[iter,] + rmax*density[iter,],0)
 #  mult <- totalpop[i,n]/sum(newden * area)
 #  density[iter,]<-newden*mult
  
  # make depth dependent as well
  # d<-unlist(as.vector(depth_grid[,2]))
  # suit = prop.b*(1-d/max(d))
  # rescale:
  # rel.suit = suit/sum(suit)
  # effort[,i+1] = rel.suit*effort[,i+1]
  # plot(d,effort[,i+1])
  # 
  #for IFD would calculate "rmax" then work out which grid cells are < rmax
  
  return(eff)
}
