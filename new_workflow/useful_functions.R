###############################################################################
# Supporting DBPM functions
# Functions have been adapted from previous DBPM work
# 
# Edited by: Denisse Fierro Arcos
# Date of update: 2024-10-16


# Loading libraries -------------------------------------------------------
library(arrow)
library(janitor)
library(dplyr)


# Calculate export ratio from tabular data --------------------------------
getExportRatio <- function(sphy, lphy, sst, depth){
  #Inputs:
  # - sphy (numeric vector) Small phytoplankton (phypico-vint variable in GFDL).
  # Must be in units of mol*m^-2.
  # - lphy (numeric vector) Large phytoplankton. This is calculated by 
  # substracting small phytoplankton (phypico-vint in GFDL) from total 
  # phytoplankton (phyc-vint in GFDL). Must be in units of mol*m^-2.
  # - sst (numeric vector) Sea surface temperature (tos variable in GFDL). Must
  # be in units of degrees Celsius.
  # - depth (numeric vector) Depth of the water column (deptho variable in 
  # GFDL). Must be in meters.
  #
  # Outputs:
  # er (numeric vector) Export ratio
  
  # Total phytoplankton
  ptotal <- sphy+lphy
  
  # Proportion of each phytoplankton group
  plarge <- lphy/ptotal
  psmall <- sphy/ptotal
  
  #Export ratio
  er <- (exp(-0.032*sst)*((0.14*psmall)+(0.74*plarge))+
           (0.0228*plarge*depth*0.004))/(1+(depth*0.004))
  
  return(er)
}


# Calculating intercept from tabular data ---------------------------------
# Method to calculate intercept from appendix of Barnes et al. 2010 JPR
GetPPIntSlope <- function(sphy, lphy, mmin = 10^-14.25, mmid = 10^-10.184, 
                          mmax = 10^-5.25, depth, output = "slope"){
  #Inputs:
  # - sphy (numeric vector) Small phytoplankton
  # - lphy (numeric vector) Large phytoplankton
  # - depth (numeric vector) Depth
  # - mmin (numeric)  Default is 10**(-14.25). ????
  # - mmid (numeric)  Default is 10**(-10.184). ????
  # - mmax (numeric)  Default is 10**(-5.25). ????
  # - output (character) Default is 'both'. Select what outputs should be 
  # returned. Choose from 'both', 'slope', or 'intercept'
  # 
  # Outputs:
  # - (numeric vector) - Depends on value of 'output' parameter. 
    
  
  # Converting sphy and lphy from mol C / m^3 to g C / m^3
  sphy <- sphy*12.0107
  lphy <- lphy*12.0107
  
  # from Appendix of Barnes 2010 JPR, used in Woodworth-Jefcoats et al. 2013 GCB
  # the scaling of biomass with body mass can be described as B=aM^b
  # the exponent b (also equivalent to the slope b in a log B vs log M 
  # relationship) can be assumed:
  # 0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
  # most studies seem to suggest N~M^-1, so can assume that and test sensitivity
  # of our results to this assumption. 
  
  #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
  #Units in log10 gww
  midsmall <- log10((mmin+mmid)/2)
  midlarge <- log10((mmid+mmax)/2)
  
  #convert to log10 (gww/size class median size) for log10 abundance
  small <- log10((sphy*10)/10^midsmall)       
  large <- log10((lphy*10)/10^midlarge)
  
  b <- (small-large)/(midsmall-midlarge)
  
  #a is really log10(a), same a when small, midsmall are used
  a <- large-(b*midlarge)  
  
  # a could be used directly to replace 10^pp in sizemodel()
  if(output =="slope"){
    return(b)
  }else if(output =="intercept"){
    return(a)
  }
}



# Calculate DBPM from GFDL outputs ----------------------------------------
dbpm_input_calc <- function(file_path, path_out){
  # Inputs:
  # - file_path (character) File path where monthly data is located
  # - path_out (character) File path where outputs should be stored as parquet
  # files
  
  #Outputs:
  # - None. This function saves results as parquet files in the path provided.
  
  df <- read_parquet(file_path) |> 
    clean_names() |> 
    rename(sphy = phypico_vint_mol_m_2) |> 
    mutate(lphy = phyc_vint_mol_m_2 - sphy) |> 
    select(!phyc_vint_mol_m_2) |> 
    rename(t = time, depth = depth_m, sbt = tob_deg_c, sst = tos_deg_c, 
           expcbot = expc_bot_mol_m_2_s_1, area_m2 = tot_area_m2) |> 
    mutate(er = getExportRatio(sphy, lphy, sst, depth),
           er = ifelse(er < 0, 0, ifelse(er > 1, 1, er)),
           intercept = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, 
                                     mmid = 10^-10.184, mmax = 10^-5.25, depth, 
                                     output = "intercept"),
           slope = GetPPIntSlope(sphy, lphy, mmin = 10^-14.25, mmid = 10^-10.184, 
                                 mmax = 10^-5.25, depth, output = "slope")) |> 
    relocate(all_of(c("er", "intercept", "slope")), .before = sphy)
  
  #Save results
  df |> 
    write_parquet(path_out)
}
