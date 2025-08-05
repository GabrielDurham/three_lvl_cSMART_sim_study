################
################
################
### PREAMBLE ###
################
################
################


#### PURPOSE: This file contains folder definitions and library calls

#### DATE CREATED:  13 OCT 2023
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         18 FEB 2025 (GJD): Added ggplot2, patchwork, and figures_path
####                22 FEB 2025 (GJD): Suppressed start-up messages

# Call libraries
suppressPackageStartupMessages({
library(readxl)
library(openxlsx)
library(dplyr)
library(MASS)
library(stats)
library(doRNG)
library(doParallel)
library(xtable)
library(writexl)
library(ggplot2)
library(patchwork)
})


# Establish paths
notes_path <- "../../00_Notes/"
pre_r_param_path <- "../../00_Notes/Pre_R_MC_Parameters/"
ASIC_parms_path <- "../../00_Notes/Sim_Parms_From_ASIC/"
helper_funcs_path <- "../Functions/"
model_fit_func_path <- "../Functions/Model_Fitting/"
data_path <- "../../02_Data/"
analysis_path <- "../../03_Analysis/"
output_path <- "../../04_Output/"
figures_path <- "../../04_Output/Figures/"