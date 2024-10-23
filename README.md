# Dynamic Benthic Pelagic Model (DBPM) calibration - ISIMIP3A protocol
This repository contains all code necessary to process inputs used by DBPM. This repository has been redesigned to use both Python and R as part of the model workflow. Following protocol ISIMIP3A, this simulation uses inputs from GFDL-MOM6-COBALT2 at two horizontal resolutions: $0.25^{\circ}$ (original) and $1^{\circ}$ (coarsen).  
  
## Step 1. Processing DBPM climate inputs at a global scale
- Script [`01_processing_dbpm_global_inputs.ipynb`](new_workflow/01_processing_dbpm_global_inputs.ipynb) processes environmental data needed to force the DBPM model at a global scale. GFDL-MOM6-COBALT2 output files are transformed from `netCDF` to analysis ready `zarr` files. Files for `spinup` period are also created here.  
  
## Step 2. Processing DBPM climate inputs at a regional scale
- Script [`02_processing_dbpm_regional_inputs.ipynb`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a/blob/new_features/new_workflow/02_processing_dbpm_regional_inputs.ipynb) uses `zarr` files produced in the previous step to extract data for an area of interest. In this notebook, we used the [FAO Major Fishing Area 58: Indian Ocean, Antarctic And Southern](https://www.fao.org/fishery/en/area/fao:58/en).  

## Step 3. Processing DBPM fishing inputs at a regional scale
- Script [`03_processing_effort_fishing_inputs.R
`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a/blob/new_features/new_workflow/03_processing_effort_fishing_inputs.R) processes fishing catch and effort data for the area of interest. It also creates a single file including fishing and climate data, which has all variables needed to run DBPM within the boundaries of the area of interest.  


# Running this repository
The scripts in this repository were developed in NCI's Gadi, so the easiest way to run these script is to clone this repository to Gadi. However, before you can do this, you will need an NCI account, which are only available for researchers with an email address from an Australian institution. Further down in this document, we include information about how to create an NCI account if you do not have one already. Remember, you must have an email address for an Australian institution to create an NCI account.  
You can also run these scripts in your own computer or a different server, but you will need need access to the forcing data (i.e., GFDL-MOM6-COBALT2 outputs and fishing data) to run them. We include information about how to access these data.  
  
## Getting an NCI account
1. [Create an NCI user account](https://access-hive.org.au/getting_started/first_steps/#create-an-nci-user-account)  
  * You should use your Australian institution’s email account when signing up  
  * When it asks for a project to join:  
    * If possible, contact the NCI scheme manager at your institution to find out what NCI project you should use to sign up for your NCI account. This account will provide you with computing power.    
2. [Join relevant NCI projects](https://access-hive.org.au/getting_started/first_steps/#join-relevant-nci-projects)
  * Request to join the following NCI projects:  
    * vf71 - for access to GFDL-MOM6-COBALT2 outputs in analysis ready data format 
    * hh5 – for the Python conda environment   
  * Note that it can take one business day or longer to get approved for projects  
3. [Verify that you can log into NCI’s Gadi](https://access-hive.org.au/getting_started/first_steps/#login-to-gadi)  
  * Note that it usually takes more than 30 minutes for your account to be created  
  * You are also welcome to follow along with the [instructions to set up ssh keys](https://access-hive.org.au/getting_started/first_steps/#automatic-login), but this is optional.  

## Accessing forcing data 
The environmental data comes from GFDL-MOM6-COBALT2, which is available at two horizontal resolutions: $0.25^{\circ}$ (original model outputs) and ($1^{\circ}$, coarsen from original outputs).  The original GFDL-MOM6-COBALT2 outputs can be downloaded from the [Inter-Sectoral Impact Model Intercomparison Project (ISIMIP) Data Portal](https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/ocean/gfdl-mom6-cobalt2/) as `netCDF` files. However, you can also access GFDL-MOM6-COBALT2 outputs as `zarr` files from project `vf71` at the [National Computational Infrastructure (NCI)](https://nci.org.au/). 
  
The fishing data is available from ????







## Step 2. LME-scale calibration  
- Script [`02_processing_dbpm_regional_inputs.ipynb`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a/blob/main/02_runLMEcalibration.R) estimates fishing mortality parameters (catchability and selectivities for each functional group)    
- Checks and adjusts search volume parameter  
- Creates and saves calibration plots as in PDF format  
- Plots can be used to visually inspect fit of predicted catch to observed catch data  
  
## Step 3. Running DBPM model (with and without fishing)
- Script [`03_processing_effort_fishing_inputs.R
`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/lme_scale_calibration_ISMIP3a/blob/main/03_rungridbyLME.R) takes parameters estimated in **step 2** to create initial values by re-running the LME scale inputs  
- Gets inputs for `gridded_sizemodel()` and runs for LME  
- Produces diagnostic plots for each LME  

# Running this repository
Although you could run all scripts included here in your local machine, you will need access to the environmental data to run them. The environmental data comes from one earth system model available at two horizontal resolutions: $0.25^{\circ}$ (original model outputs) and ($1^{\circ}$, coarsen from original outputs).  

Outputs from GFDL-MOM6-COBALT2 are originally available at the [Inter-Sectoral Impact Model Intercomparison Project (ISIMIP)](https://www.isimip.org/) as netCDF files. We extracted data for each LME and transformed files to `zarr` and `parquet` formats, which are now available under project `vf71` at the [National Computational Infrastructure (NCI)](https://nci.org.au/). 
  
The scripts in this repository were developed in NCI's Gadi, so the easiest way to run these script is to clone this repository to Gadi. However, to do this you will need an NCI account, which are only available for researchers with an email address from an Australian institution. You can still run these scripts locally, but you will need to download the environmental data from the [ISIMIP Data Portal](https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/ocean/gfdl-mom6-cobalt2/) or NCI.    
  
In the next section, we include information about how to create an NCI account if you do not have one already. Remember, you must have an email address for an Australian institution to create an NCI account.  
  
## Getting an NCI account
1. [Create an NCI user account](https://access-hive.org.au/getting_started/first_steps/#create-an-nci-user-account)  
  * You should use your Australian institution’s email account when signing up  
  * When it asks for a project to join:  
    * If possible, contact the NCI scheme manager at your institution to find out what NCI project you should use to sign up for your NCI account. This account will provide you with computing power.    
2. [Join relevant NCI projects](https://access-hive.org.au/getting_started/first_steps/#join-relevant-nci-projects)
  * Request to join the following NCI projects:  
    * vf71 - for access to GFDL-MOM6-COBALT2 outputs in analysis ready data format 
    * hh5 – for the Python conda environment   
  * Note that it can take one business day or longer to get approved for projects  
3. [Verify that you can log into NCI’s Gadi](https://access-hive.org.au/getting_started/first_steps/#login-to-gadi)  
  * Note that it usually takes more than 30 minutes for your account to be created  
  * You are also welcome to follow along with the [instructions to set up ssh keys](https://access-hive.org.au/getting_started/first_steps/#automatic-login), but this is optional.  


