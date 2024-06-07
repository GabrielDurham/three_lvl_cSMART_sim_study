################
################
################
### PREAMBLE ###
################
################
################


#### PURPOSE: This file contains folder definitions and library calls

#### DATE CREATED: 13 OCT 2023
#### PROGRAMMER: GABRIEL DURHAM (GJD)
#### EDITS: 

# Call libraries
library(readxl)
library(openxlsx)
library(dplyr)
library(MASS)
library(stats)
library(doRNG)
library(doParallel)
library(xtable)
# suppressWarnings(library(readxl))
# suppressWarnings(library(openxlsx))
# suppressWarnings(library(dplyr))
# suppressWarnings(library(MASS))
# suppressWarnings(library(stats))
# suppressWarnings(library(doRNG))

# Establish paths
notes_path <- "../00_Notes/"
pre_r_param_path <- "../00_Notes/Pre_R_MC_Parameters/"
helper_funcs_path <- "Functions/"
model_fit_func_path <- "Functions/Model_Fitting/"
data_path <- "../02_Data/"
analysis_path <- "../03_Analysis/"
output_path <- "../04_Output/"