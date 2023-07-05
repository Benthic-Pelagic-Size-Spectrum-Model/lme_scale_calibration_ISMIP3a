# lme_scale_calibration_ISMIP3a

STEPS:

Step 1. LME-scale calibration:
"runLMEcalibration.RDS"
-this code estimates fihsing mortality paramters ( catchability and selectivities for each functional group)
- need to check and adjust search volume parameter
- Creates "CalibrationPlots.pdf" , use this to visually inspect fits to ctach data

Step 2. For each LME run model with and without fishing
- Initial Tests carried with LME 14
- Takes parameter estimates from Step 1 to create initial values by re-running the LME scale inputs
- Get inputs for gridded_sizemodel() and runs for LME
- Produces diagnostic plots for each LME
