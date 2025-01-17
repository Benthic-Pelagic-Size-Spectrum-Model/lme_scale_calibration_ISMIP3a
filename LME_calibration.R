###############################################################################
# Library of functions developed originally by JLB on 21/12/2022
#
# Functions were updated by Denisse Fierro Arcos so they can be used with data
# produced by updated `01_getinputs_ISIMIP3a.R` script.
#
# Date of update: 2024-08-06

# Loading libraries
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggplot2)
library(lubridate)
library(zoo)
library(lhs)
library(patchwork)
library(purrr)
library(janitor)
library(parallel)
# library(optimParallel)
source("dbpm_model_functions_NEW_DFA.R")
# library(gridExtra)


# Adding outputs to gridded_sizemodel ----
getGriddedOutputs <- function(input = lme_inputs_grid, results = grid_results,
                              params = params){
  # returns all outputs of the model 
  input$TotalUbiomass <- input$TotalVbiomass <- input$TotalUcatch <-
    input$TotalVcatch<- input$Totalcatch <- NA
  # Merging coordinates together to use as unique identifier
  input <- input |>
    unite("cell", lat, lon, remove = F)
  cells <- unique(input$cell)
  
  # 2:(params$Neq+1) changed to 1:params$Neq because the newly processed data
  # does subsets data based on time ranges. It does not include December data
  # for year prior
  for(igrid in 1:length(cells)){
    input[input$cell == cells[igrid], ]$TotalUbiomass <- 
      apply(results$U[igrid, params$ref:params$Nx, 1:(params$Neq)] * 
              params$dx*10^params$x[params$ref:params$Nx],
            2, sum, na.rm = T)
    
    input[input$cell == cells[igrid], ]$TotalVbiomass <- 
      apply(results$V[igrid, params$ref:params$Nx, 1:(params$Neq)] * 
              params$dx*10^params$x[params$ref:params$Nx],
            2, sum, na.rm = T)
    
    #sum catches (currently in grams per m3 per year, across size classes) 
    #keep as grams per m2, then be sure to convert observed from tonnes per m2 
    #per year to g.^-m2.^-yr (for each month)
    input[input$cell == cells[igrid], ]$TotalUcatch <- 
      apply(results$Y_u[igrid, params$ref:params$Nx, 1:(params$Neq)] *
              params$dx, 2, sum, na.rm = T)
    
    input[input$cell == cells[igrid], ]$TotalVcatch <- 
      apply(results$Y_v[igrid, params$ref:params$Nx, 1:(params$Neq)] *
              params$dx, 2, sum, na.rm = T)
    
    input[input$cell == cells[igrid], ]$Totalcatch <- 
      input[input$cell == cells[igrid], ]$TotalUcatch +
      input[input$cell == cells[igrid], ]$TotalVcatch
  }
  ## and then multiply outputs by depth to get per m2
  return(input)
}
  
# Set up and run optimisations in parallel ----
fastOptim <- function(LMEnum, vary, forcing_file = NULL, 
                      gridded_forcing = NULL, errorFun = getError, ...,
                      spareCores = 1, libraries = c("optimParallel")){
  
  # get inputs for LME
  if(!is.null(forcing_file)){
    lme_input <- forcing_file
  }
  if(!is.null(gridded_forcing)){
    lme_input <- gridded_forcing
  }
  
  args <- list(...)
  if("corr" %in% names(args)){
    corr <- args$corr
  }else{
    corr <- F
  }
  if("figure_folder" %in% names(args)){
    figure_folder <- args$figure_folder
  }else{
    figure_folder <- NULL
  }
  
  # set up workers
  # keep some spare core
  noCores <- detectCores()-spareCores 
  if(noCores < 1){
    stop("You should allow at least one core for this operation.")
  }
  cl <- makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, varlist = c("cl", "libraries", "lme_input"),
                envir = environment())
  clusterEvalQ(cl, {
    for(item in 1:length(libraries)){
      library(libraries[item], character.only = T)
    }
  })
  
  clusterEvalQ(cl, source("LME_calibration.R"))
  
  optim_result <- optimParallel(par = vary, fn = errorFun,
                                method = "L-BFGS-B",
                                lower = rep(1e5, length(vary)),
                                upper = rep(2, length(vary)),
                                parallel = list(loginfo = TRUE,
                                                cl = cl,
                                                forward = TRUE), 
                                lme_forcings = lme_input, corr = corr, 
                                figure_folder = figure_folder)
  stopCluster(cl)
  
  return(optim_result$par)
}


