# three_lvl_cSMART_sim_study
This repository contains R code and supporting files to simulate longitudinal 
data from a three-level clustered SMART, generate inference results, and 
create final output tables for publication.

- 00_Notes/.  *# Files which help facilitate the flow of the simulation*
  - Driver_Three_Level_Simulations.xlsx *# A driver file which controls the overhead of the simulation study*
  - Consistency_Check/. *# Results of sanity-checking whether data-generative code produces data with specified distributions*
  - Pre_R_MC_Parameters/. *# Monte Carlo estimates of pre-response conditional parameters*
  - Sim_Parms_From_ASIC/. *# Parameters from ASIC estimates, which will be cleaned and put into the driver file*
  - Model_Fitting_Wrapper.R *# Wrapper functions which facilitate implementing model described in https://arxiv.org/abs/2503.08987 for prototypical cSMARTs. Made available to ease implementation of method for external researchers. Not used in simulation code*
- 01_Programs/. *# R code which implements the simulations*
  - Process_ASIC_Parms.Rmd *# Processes raw ASIC estimates and makes slight adjustments to jibe with data-generative model*
  - Pre_R_Cond_Parm_MC.Rmd *# Runs simulations of pre-response data to estimate pre-response conditional means/variances*
  - Run_MC_Analysis.Rmd *# Simulate and fit data as dictated in driver file. Conduct inference*
  - Create_Tables.Rmd *# Create camera-ready output*
  - Create_Tables.html *# Shows further simulation results and MC confidence intervals*
  - Check_Distribution_Consistency.Rmd *# Check whether data-generative code is producing data with specified distributions*
  - Functions/ *# Helper functions*
- 02_Data/. *# Raw output from each simulation iteration*
- 03_Analysis/. *# Processed simulation output*
- 04_Output/. *# Camera-ready output*

*Driver_Three_Level_Simulations.xlsx* contains all the data-generative parameters for all simulation runs. 
The "Dictionary" sheet contains a description of each necessary field. Each row of the "main" sheet represents
a unique data-generative environment to simulate. To simulate/fit data according to this model: (1) Specify the pre-response
variance structure in the "var_parm_settings_pre_r" sheet in the driver file; (2) Run Pre_R_Cond_Parm_MC.Rmd,
which will calculate/store pre-response conditional parameters; (3) Specify the overall mean/variance structure, as
well as the model fitting specifications, in the "main" sheet of the driver file (which will require specifying distributions in
the "cond_param_settings" and "var_parm_settings_pre/post_r" sheets). Run the Run_MC_Analysis.Rmd file. These steps will
create inference output in the 03_Analysis/ folder, which Create_Tables.Rmd can turn into camera-ready results.




These simulations run in R, and require the following packages:
readxl, openxlsx, dplyr, MASS, stats, doRNG, doParallel, xtable, writexl, ggplot2, patchwork.
This code is provided for academic use.
