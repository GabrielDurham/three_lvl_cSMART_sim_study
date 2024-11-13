################################
################################
################################
### GENERAL DATA MAINTENANCE ###
################################
################################
################################

#### PURPOSE: This file contains code that serves to support general data
####          flow and maintenance for the simulation study.


#### DATE CREATED: 23 OCT 2023
#### PROGRAMMER: GABRIEL DURHAM (GJD)
#### EDITS: 14 JAN 2024 (GJD) - Reworked Process_Driver_Row_Main() function
####        25 JAN 2024 (GJD) - Fixed Get_Export_File_Name() function (added suffix)
####                            Fixed Partition_Vector() function so it wouldn't create
####                              too few partitions
####        15 MAY 2024 (GJD) - Added alternate fit settings to Read_Driver_File()
####        24 MAY 2024 (GJD) - Modified Import_Pre_R_MC_Results() to combine data
####                              data from multiple runs if specified in driver file
####        27 MAY 2024 (GJD) - Added Summarize_Emp_Data()
####        03 JUN 2024 (GJD) - Added outcome_prefix input to Summarize_Emp_Data()
####        10 JUL 2024 (GJD) - Added functionality to Import_Pre_R_MC_Results()
####                              to allow for centering errors
####        01 AUG 2024 (GJD) - Dropped unnecessary columns in cp_settings driver
####                              sheet
####        10 SEP 2024 (GJD) - Fixed issue in processing driver file row where it didn't
####                              read NAs correctly

#Read in Driver File
### path = File location for driver file
### driver_sheet = Sheet name for simulation driver
### cp_sheet = Sheet name for conditional parameter settings (default "cond_param_settings")
### var_parm_pre_r_sheet = Sheet name for pre-response variance settings
### var_parm_post_r_sheet = Sheet name for post-response variance settings
### fit_settings_sheet = Sheet name for model fit settings
### alt_fit_settings_sheet = List of vectors c(alt-fit type, 
###                                            corresponding sheet name)
#### Outputs list with simulation driver and conditional parameter settings
Read_Driver_File <- function(path, driver_sheet="main", 
                             cp_sheet="cond_param_settings",
                             var_parm_pre_r_sheet="var_parm_settings_pre_r",
                             var_parm_post_r_sheet="var_parm_settings_post_r",
                             fit_settings_sheet="fit_settings",
                             alt_fit_settings_types_sheets=NULL){
  Output <- NULL
  Output[["driver"]] <- as.data.frame(read_xlsx(path=path, sheet=driver_sheet, 
                                                col_names=TRUE))
  Output[["cp_settings"]] <- as.data.frame(read_xlsx(path=path, sheet=cp_sheet, 
                                                     col_names=TRUE))
  #Don't store unnecessary cp_settings columns
  cps_cols_to_drop <- c("Notes")
  Output[["cp_settings"]] <- 
    Output[["cp_settings"]][ , !names(Output[["cp_settings"]]) %in% cps_cols_to_drop]
  Output[["var_parms_pre_r"]] <- as.data.frame(read_xlsx(path=path, 
                                                         sheet=var_parm_pre_r_sheet, 
                                                         col_names=TRUE))
  Output[["var_parms_post_r"]] <- as.data.frame(read_xlsx(path=path, 
                                                          sheet=var_parm_post_r_sheet, 
                                                          col_names=TRUE))
  Output[["fit_settings"]] <- as.data.frame(read_xlsx(path=path, 
                                                      sheet=fit_settings_sheet, 
                                                      col_names=TRUE))
  if (!is.null(alt_fit_settings_types_sheets)) {
    Output[["alt_fit_settings"]] <- NULL
    for (type_sheet in alt_fit_settings_types_sheets) {
      type <- type_sheet[1]
      sheet_name <- type_sheet[2]
      Output[["alt_fit_settings"]][[type]] <- as.data.frame(read_xlsx(path=path, 
                                                                      sheet=sheet_name, 
                                                                      col_names=TRUE))
    }
  }
  
  
  return(Output)
}



# Get File Name for Exporting Data
## Scans path folder for existing files to ensure no overwriting. Suffix is
## YYYYMMDD date format
### path = Path where past results are stored
### prefix = Prefix of file names (i.e., file name prior to date suffix)
### file_suffix = File type "should be ".filetype"
Get_Export_File_Name <- function(path, prefix, file_suffix=".xlsx"){
  current_date <- strsplit(toString(Sys.Date()), split='-')[[1]]
  date_suffix <- paste0(current_date[1], current_date[2], current_date[3])
  
  file_names <- list.files(path, pattern = paste0("^", prefix), all.files = TRUE)
  
  conflicting_file_names <- file_names[grepl(paste0("^", prefix, date_suffix), file_names)]
  
  n_existing_files <- length(conflicting_file_names)
  
  if (n_existing_files==0) {Output <- paste0(prefix, date_suffix, file_suffix)
  } else {Output <- paste0(prefix, date_suffix, letters[n_existing_files], file_suffix)}
  return(Output)
}
# Import Pre-Response Monte Carlo Results
## Imports most recent version of results
### path = Path where past results are stored
### prefix = Prefix of file names (i.e., file name prior to date suffix)
### pre_r_marg_parms = Driver sheet of pre-response marginal parameters
###                     Sheet will specify stored run for given settings, if none
###                     then we assume it's in the most recent run
Import_Pre_R_MC_Results <- function(path, prefix, pre_r_marg_parms=NULL, center=FALSE) {
  file_names <- list.files(path, pattern = paste0("^", prefix), all.files = TRUE)
  recent_file <- sort(file_names)[length(file_names)]
  
  # Check to make sure there are any settings in the most recent file and/or
  # if there are any specified stored runs
  if (!is.null(pre_r_marg_parms)) {
    any_most_recent <- sum(is.na(pre_r_marg_parms$stored_run))>0
    any_stored_runs <- sum(!is.na(pre_r_marg_parms$stored_run))>0
  } else {
    any_most_recent <- TRUE
    any_stored_runs <- FALSE
  }
  # If so, pull them
  if (any_most_recent) {
    temp_output_r <- read.xlsx(xlsxFile=paste0(path=path, recent_file), 
                               sheet="Responders")
    temp_output_nr <- read.xlsx(xlsxFile=paste0(path=path, recent_file), 
                                sheet="Nonresponders")
  } else {
    temp_output_r <- data.frame()
    temp_output_nr <- data.frame()
  }
  # If there are any stored runs, pull them
  if (any_stored_runs) {
    stored_settings <- pre_r_marg_parms[!is.na(pre_r_marg_parms$stored_run), ]
    for (row in rownames(stored_settings)) {
      setting <- stored_settings[row, "cond_param_setting_pre_r"]
      stored_data_file_name <- stored_settings[row, "stored_run"]
      # Access full MC data from the run
      stored_data_r <- read.xlsx(xlsxFile=paste0(path=path, stored_data_file_name), 
                                 sheet="Responders")
      stored_data_nr <- read.xlsx(xlsxFile=paste0(path=path, stored_data_file_name), 
                                  sheet="Nonresponders")
      
      # Add the runs for the given setting
      temp_output_r <- 
        rbind(temp_output_r, 
              stored_data_r[stored_data_r$cond_param_setting_pre_r==setting, ])
      temp_output_nr <- 
        rbind(temp_output_nr, 
              stored_data_nr[stored_data_r$cond_param_setting_pre_r==setting, ])
    }
  }
  # Re-order and re-label rows
  temp_output_r <- temp_output_r[order(temp_output_r$cond_param_setting_pre_r, 
                                       temp_output_r$n_i), ]
  rownames(temp_output_r) <- (1:nrow(temp_output_r))
  temp_output_nr <- temp_output_nr[order(temp_output_nr$cond_param_setting_pre_r, 
                                         temp_output_nr$n_i), ]
  rownames(temp_output_nr) <- (1:nrow(temp_output_nr))
  
  # Center MC errors so they average to zero
  if (center) {
    # Pull actual p_r to center, otherwise take empirical rate
    if (!is.null(pre_r_marg_parms)) {
      p_r <- merge(x=temp_output_r, 
                   y=pre_r_marg_parms[,c("cond_param_setting_pre_r", "p_r")],
                   by="cond_param_setting_pre_r", all.x=TRUE)$p_r
    } else {
      n_r <- temp_output_r$n_obs
      n_nr <- temp_output_nr$n_obs
      p_r <- n_r/(n_r+n_nr)
    }
    for (var in c("e_0", "e_1", "y_0", "y_1")) {
      mean_lab <- paste0("mean_", var)
      temp_output_r[[paste0("old_", mean_lab)]] <- temp_output_r[[mean_lab]]
      temp_output_nr[[paste0("old_", mean_lab)]] <- temp_output_nr[[mean_lab]]
      
      center_var <- p_r*temp_output_r[[mean_lab]] + (1-p_r)*temp_output_nr[[mean_lab]]
      
      temp_output_r[[mean_lab]] <- temp_output_r[[mean_lab]]-center_var
      temp_output_nr[[mean_lab]] <- temp_output_nr[[mean_lab]]-center_var
    }
  }
  
  
  Output <- NULL
  Output[["R"]] <- temp_output_r
  Output[["NR"]] <-temp_output_nr
  return(Output)
}


# Partition a Vector into K Parts - Used for parallelization
## Used for splitting up driver rows for parallel computing
### vec = List to partition
### k = Number of parts
Partition_Vector <- function(vec, k){
  n_items <- length(vec)
  part_size <- ceiling((n_items-k)/k)
  
  partition_groups <- (1:k)
  for (part in (1:k)) {
    partition_groups <- c(partition_groups, rep(x=part, times=part_size))
  }
  partition_groups <- sort(partition_groups[(1:n_items)])
  Output <- split(x=vec, f=partition_groups)
  return(Output)
}




# Process a Driver Row - Main Simulation Driver
### driver_row = Row of driver file to process
#### Returns a list with simulation parameters
Process_Driver_Row_Main <- function(driver_row){
  Output <- NULL
  # Make Output a list of all columns of driver_row
  for (col in colnames(driver_row)) {Output[[col]] <- driver_row[[col]][1]}
  # Overwrite specific list elements
  Output[["cp_setting"]] <- driver_row$cond_parm_setting[1]
  #Output[["n_clusters"]] <- driver_row$n_clusters[1]
  Output[["cluster_sizes"]] <- eval(parse(text=driver_row[1,"cluster_sizes"]))
  
  # Default cluster size randomization probabilities are balanced
  if (is.na(driver_row$p_cluster_size)[1]) {
    Output[["p_cluster_size"]] <- rep(1/length(Output[["cluster_sizes"]]),
                                      length(Output[["cluster_sizes"]]))
  } else {
    Output[["p_cluster_size"]] <- eval(parse(text=driver_row[1,"p_cluster_size"]))
  }
  
  # Default DTR assignment probabilities are balanced
  if (is.na(driver_row$p_dtr)[1]) {
    if (driver_row$SMART_structure=="prototypical") {
      Output[["p_dtr"]] <- rep(0.25, 4)
    }
  } else {
    Output[["p_dtr"]] <- eval(parse(text=driver_row[1,"p_dtr"]))
  }
  Output[["min_dtr_obs"]] <- eval(parse(text=driver_row[1,"min_dtr_obs"]))
  Output[["min_path_obs"]] <- eval(parse(text=driver_row[1,"min_path_obs"]))
  
  Output[["covars"]] <- c()
  if (driver_row$n_covar>0) {
    for (i in (1:driver_row$n_covar)) {
      Output[["covars"]] <- c(Output[["covars"]], paste0("x_", i))
    }
  }
  # Default to no alternate models
  if (is.na(driver_row$n_alt_fit)) {Output[["n_alt_fit"]] <- 0}
  
  return(Output)
}


# Define SMART Pathways
### SMART_structure = Structure of SMART (e.g., "prototypical")
#### Returns data frame with embedded DTRs and their response/nonresponse pathway
Define_SMART_Pathways <- function(SMART_structure){
  temp_output <- data.frame(dtr=character(), path_r=numeric(), path_nr=numeric())
  if (SMART_structure=="prototypical") {
    temp_output <- rbind(temp_output, 
                         data.frame(dtr="d_1_1", path_r=1, path_nr=2, 
                                    a_1=1, a_2r=NA, a_2nr=1))
    temp_output <- rbind(temp_output, 
                         data.frame(dtr="d_1_m1", path_r=1, path_nr=3, 
                                    a_1=1, a_2r=NA, a_2nr=-1))
    temp_output <- rbind(temp_output, 
                         data.frame(dtr="d_m1_1", path_r=4, path_nr=5, 
                                    a_1=-1, a_2r=NA, a_2nr=1))
    temp_output <- rbind(temp_output, 
                         data.frame(dtr="d_m1_m1", path_r=4, path_nr=6, 
                                    a_1=-1, a_2r=NA, a_2nr=-1))
  }
  
  Output <- temp_output
  return(Output)
}

# Calculate Weighted Variance
### x = Numerical vector containing values for which weighted variance 
###     to be computed
### w = Numerical vector of weights
weighted.var <- function(x, w) {
  n <- length(x)
  x_bar_w <- weighted.mean(x=x, w=w)
  num <- sum(w*((x-x_bar_w)^2))
  denom <- ((n-1)/n)*sum(w)
  Output <- num/denom
  return(Output)
}

# Calculate Weighted Covariance
### x = Numerical vector containing values for which weighted covariance 
###     to be computed
### y = Numerical vector containing values for which weighted covariance 
###     to be computed
### w = Numerical vector of weights
weighted.cov <- function(x, y, w) {
  n <- length(x)
  x_bar_w <- weighted.mean(x=x, w=w)
  y_bar_w <- weighted.mean(x=y, w=w)
  num <- sum(w*((x-x_bar_w)*(y-y_bar_w)))
  denom <- ((n-1)/n)*sum(w)
  Output <- num/denom
  return(Output)
}

# Calculate Empirical Mean/Variance for a Single Setting
### in_wide_data = Data with row per cluster, outcomes organized as 
###                paste0("y_p_" person_id, "_t_", time)
###                person_id should be 1,..., n_i
###                times should be c(0,1,...)
### outcome_prefix = Prefix of outcome. E.g., "y" means
###                  outcomes organized as
###                  paste0("y_p_" person_id, "_t_", time) 
###                  "e" means outcomes organized as
###                  paste0("e_p_" person_id, "_t_", time)
### cluster_size = Number of units in each cluster
### times = Vector of times
### weight_col = Column name of weights, if NULL then will not use weights
Summarize_Emp_Data <- function(in_wide_data, cluster_size, 
                               outcome_prefix="y",
                               times=c(0,1,2),
                               weight_col=NULL) {
  people <- (1:cluster_size)
  T <- length(times)
  # Assign weights
  if (is.null(weight_col)) {weights <- rep(1, nrow(in_wide_data))
  } else {weights <- in_wide_data[[weight_col]]}
  
  # Get total distribution
  temp_var <- matrix(NA, nrow=T*cluster_size, ncol=T*cluster_size)
  temp_mean <- matrix(NA, nrow=T*cluster_size, ncol=1)
  for (p_1 in people) {
    for (t_1 in times) {
      row <- (p_1-1)*length(times) + match(t_1, times)
      row_y_label <- paste0(outcome_prefix, "_p_", p_1, "_t_", t_1)
      
      temp_mean[row, 1] <- weighted.mean(x=in_wide_data[[row_y_label]], 
                                         w=weights)
      temp_var[row, row] <- weighted.var(x=in_wide_data[[row_y_label]],
                                         w=weights)
      for (p_2 in people[people>=p_1]) {
        #Only loop through person/time combos we haven't seen yet
        times_2 <- times
        if (p_1==p_2) {if (t_1<max(times)) {times_2 <- times[times>t_1]
        } else {times_2 <- NA}
        }
        
        if (!anyNA(times_2)) {for (t_2 in times_2) {
          col <- (match(p_2, people)-1)*length(times) + match(t_2, times)
          col_y_label <- paste0(outcome_prefix, "_p_", p_2, "_t_", t_2)
          # Calculate covariance
          covar <- weighted.cov(x=in_wide_data[[row_y_label]], 
                                y=in_wide_data[[col_y_label]],
                                w=weights)
          temp_var[row, col] <- covar
          temp_var[col, row] <- covar
        }}
      }
    }
  }
  emp_mean_vec <- temp_mean
  emp_var_mat <- temp_var
  
  # Summarize Distribution
  avg_mean_parms <- NULL
  for (t in times) {
    total_obs_t <- c()
    total_wgts <- c()
    for (person in (1:cluster_size)) {
      total_obs_t <- c(total_obs_t, 
                       in_wide_data[[paste0(outcome_prefix, "_p_", person, 
                                            "_t_", t)]])
      total_wgts <- c(total_wgts, weights)
    }
    avg_mean_parms[[paste0(outcome_prefix, "_", t)]] <- 
      weighted.mean(x=total_obs_t, w=total_wgts)
  }
  # Collect all diagonal blocks
  d_blocks <- vector("list", length=cluster_size)
  for (i in (1:cluster_size)) {
    start_index <- (i-1)*length(times) + 1
    end_index <- start_index + length(times) - 1
    d_blocks[[i]] <- emp_var_mat[(start_index:end_index), (start_index:end_index)]
  }
  # Collect all off-diagonal blocks
  od_blocks <- vector("list", length=(cluster_size^2-cluster_size)/2)
  index <- 1
  for (i in (1:(cluster_size-1))) {for (ii in ((i+1):cluster_size)) {
    start_index_row <- (i-1)*length(times) + 1
    end_index_row <- start_index_row + length(times) - 1
    
    start_index_col <- (ii-1)*length(times) + 1
    end_index_col <- start_index_col + length(times) - 1
    
    od_blocks[[index]] <- 
      emp_var_mat[(start_index_row:end_index_row), (start_index_col:end_index_col)]
    od_blocks[[index+1]] <- 
      emp_var_mat[(start_index_col:end_index_col), (start_index_row:end_index_row)]
    index <- index+2
  }}
  # Average elements of matrices
  # Diagonal blocks
  avg_d_block <- matrix(NA, nrow=length(times), ncol=length(times))
  for (row in (1:nrow(avg_d_block))) {for (col in (1:ncol(avg_d_block))) {
    element_vec <- c()
    for (i in (1:length(d_blocks))) {
      element_vec <- c(element_vec, d_blocks[[i]][row, col])
    }
    avg_d_block[row, col] <- mean(element_vec)
  }}
  # Off-Diagonal Blocks
  avg_od_block <- matrix(NA, nrow=length(times), ncol=length(times))
  for (row in (1:nrow(avg_od_block))) {for (col in (1:ncol(avg_od_block))) {
    element_vec <- c()
    for (i in (1:length(od_blocks))) {
      element_vec <- c(element_vec, od_blocks[[i]][row, col])
    }
    avg_od_block[row, col] <- mean(element_vec)
  }}
  
  # Re-Organize as Parameter Names
  emp_var_parms <- NULL
  for (t_1 in times) {for (t_2 in times[times>=t_1]) {
    if (t_1==t_2) {
      emp_var_parms[[paste0("sigma2_", t_1)]] <- avg_d_block[(t_1+1), (t_1+1)]
      emp_var_parms[[paste0("rho_", t_1)]] <- avg_od_block[(t_1+1), (t_1+1)]
    } else {
      emp_var_parms[[paste0("phi_", t_1, t_2)]] <- avg_d_block[(t_1+1), (t_2+1)]
      emp_var_parms[[paste0("rho_", t_1, t_2)]] <- avg_od_block[(t_1+1), (t_2+1)]
    }
  }}
  
  
  Output <- list(mean=list(vec=emp_mean_vec,
                           parms=avg_mean_parms),
                 var=list(mat_full=emp_var_mat,
                          mat_means=list(diag=avg_d_block,
                                         off_diag=avg_od_block),
                          parms=emp_var_parms))
  return(Output)
}