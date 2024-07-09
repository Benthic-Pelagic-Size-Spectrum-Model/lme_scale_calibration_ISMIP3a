# Dynamic Benthic Pelagic Model (DBPM) calibration - ISIMIP3A protocol
This repository contains all code necessary to process inputs used by DBPM. Following protocol ISIMIP3A, this simulation uses inputs from GFDL-MOM6-COBALT2 (horizontal resolution: $0.25^{\circ}$).  
  
Two additional DBPM experiments forced by ACCESS-OM2-025 (horizontal resolution: $0.25^{\circ}$) and ACCESS-OM2-01 (horizontal resolution: $0.10^{\circ}$).  
  
**STEPS:**  
  
Step 1. LME-scale calibration: `runLMEcalibration.RDS`  
- this script estimates fishing mortality parameters (catchability and selectivities for each functional group)  
- need to check and adjust search volume parameter  
- creates `CalibrationPlots.pdf`, use this to visually inspect fits to catch data  
  
Step 2. For each LME run model with and without fishing  
- Initial tests carried out with LME 14  
- Takes parameter estimates from Step 1 to create initial values by re-running the LME scale inputs  
- Get inputs for `gridded_sizemodel()` and runs for LME  
- Produces diagnostic plots for each LME




