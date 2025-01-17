
#Method to calculate intercept from appendix of Barnes et al. 2010 JPR


# Added 1 May, 2013
# Methods to calculate:

# 1. primary producer size distribution (as described by MB10, 50, 90)
# 2. primary prodcuer slope and intercept (from Barnes et al. and Simon's notes)
# 3.detritus input dependent on primary production and depth - proportion of sinking phytoplankton for detritus input
# 4. minimum consumer size
# 5. input for preferred predator-prey mass ratio function - beta and sigma

# These are calculated from measurements of temperature, primary production 

GetPPSpectrum<-function(pp = 0.77, sst = 20, chla=NA) {
  
  # List of parameters
  # sst= temperature (celsius) - INPUT, in degrees C
  # P= primary production Note: needs to be in g C m-2 d-1  - INPUT
  # chla = chl a biomass Note: mg Chl m-3, if used instead of pp input !=NA
  
  
  # OUTPUTS:
  # a= intercept phytoplankton size spectrum (orginal units are (x) cell mass pg C and (y) log10 pg C m-3 converted to g vs. number density per m^3)
  # b= slope phytoplankton size spectrum 
  # mb50= as Barnes et al
  
  # /* phytoplankton size spectrum parameters */
  
  log10pp=log10(pp*1000)
  
  # predict log 10 (pg C) median phytoplankton cell size
  
  mb50log=(-0.081*sst)+(0.875*log10pp)-0.661
  
  # predict phytoplankton ss intercept from PP as g C m-2 d-1
  # spectrum units are log10 pg C m-3 (y) v cell mass log10 pg C (x)
  
  mbint=(0.554*log10pp)+8.087
  # predict phytoplankton ss slope from PP as g C m-2 d-1 
  # spectrum units are log10 pg C m-3 (y) v cell mass log10 pg C (x)
  
  mbslope=(0.182*log10pp)-1.715
  
  
  # fixed log predator=prey mass ratio to calculate minimum consumer size in the model
  
  beta = 2.0
  
  # need to convert log 10 pg C to g C (10^-12) then g C to g wet weight (10^-1; 1 g C = 10 g ww)
  
  mb50a = 10^(mb50log)*10^-12*10^-1
  
  # mb50 now in log 10 (grams wet weight)
  
  mb50 = log10(mb50a)
  
  if (any(mbslope > -1.05)){mbslope = -1.05}  # limits to prediction
  if (any(mbslope < -1.45)){mbslope = -1.45}
  
  
  return(list(
    
    # OUTPUTS: 
    # a is intercept - log 10 (pg C m^-3) - convert to log10 (g ww m^3) -> log10 (density m^3) 
    a=log10((10^mbint)*10^-13),
    # pp slope
    b = mbslope,
    # log 10 (median pp size, g ww),
    med.x = mb50,
    # estimate min consumer size from PPMR and median phytoplankton size ( see also Woodworth et al 2012)
    min.cons.x = log10(10^beta*10^mb50)
    
  ))
  
}




optimizeYield <- function(params,com,pp=pp, sst=sst, sft=sft, U_mat=U, V_mat=V,sizeindex=(1:params$Nx)) {
  
  
  #------------------------------------------------------------------------------------- 
  #
  # Function to Optimise yield , based on a single parameter, (F-even, a single F across all sizes) 
  # 
  #-------------------------------------------------------------------------------------
  
  
  ## ------------------------------------------------------------------
  ##  Do a single-parameter optimization of  yield
  ## ------------------------------------------------------------------
  
  
  # Helper function to return the quantity to be optimized based on the
  # fishing mortality :
  #
  
  
  
  totalYield <- function(newF,params, com,sizeindex) {
    
    print(newF)
    
    params$Fmort = newF
    
    
    params$tstepdays = 1.0
    params$tmaxyears = 10   #change to 10 years( not running all the way to equilib within this)
    params$Neq=as.integer(365.0*params$tmaxyears/params$tstepdays)	# total number of steps in a period
    params$delta_t=params$tmaxyears/params$Neq               		# time step (rates per year)  
    params$siz_time=1:params$Neq
    
    
    params$U.init[,1]=com$U
    params$V.init[,1]=com$V
    params$W.init[1]=com$W 
    
    com <- sizemodel(params=params, pp, sst, sft, U_mat, V_mat, temp.effect=T) 
    
    
    # mean yield through time for specific size range ( or entire fished range)    
    
    yield = com$Y[sizeindex]
    
    
    totalYield <- sum(yield)*1e-6
    
    
    -totalYield
  }
  
  
  result <- optimize(f=totalYield,interval=c(0,1.5),maximum=F,params=params, com=com)
  
  # result <- nlm(f=totalYield,p=0.8,params=params, com=com)
  
  # result <- genoud(fn=totalYield, nvars=1,Domains=cbind(0,2), max=FALSE,params=params, com=com)
  
  # result<-GenSA(par=0.1, lower=0, upper=3, fn=totalYield,params=params,com=com,control=list(maxit=100))
  
  return(list(result,com))
  
}

