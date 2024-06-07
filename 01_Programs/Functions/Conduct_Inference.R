#########################
#########################
#########################
### CONDUCT_INFERENCE ###
#########################
#########################
#########################


#### PURPOSE: This file contains functions for conducting inference on output from
####          driver runs

#### DATE CREATED: 21 MAY 2024
#### PROGRAMMER: GABRIEL DURHAM (GJD)
#### EDITS: 



# Contrast row definitions #

# Create a Contrast Definition Row for a Mean Difference at Single Time
### driver_parms = Output of Process_Driver_Row_Main() function
### t = Time we're comparing mean outcomes at
### crit_t = Critical t value (time of second decision). Default is 1
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### a_i = i^th decision in DTR 1
### a_ip = i^th decision in DTR 2
Create_Individual_Contrast_Row_Single_t <- function(driver_parms, t, crit_t=1, 
                                                    model_fit_params, a_1, a_1p, 
                                                    a_2r=NULL, a_2rp=NULL,
                                                    a_2nr=NULL, a_2nrp=NULL) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter
  t_diff <- t-crit_t
  t_ind <- ifelse(t>crit_t, 1, 0)
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    temp_output[1,"gamma_2"] <- ((1-t_ind)*t + t_ind*crit_t)*(a_1-a_1p)
    temp_output[1,"gamma_4"] <- t_ind*t_diff*(a_1-a_1p)
    temp_output[1,"gamma_5"] <- t_ind*t_diff*(a_2nr-a_2nrp)
    temp_output[1,"gamma_6"] <- t_ind*t_diff*((a_1*a_2nr) - (a_1p*a_2nrp))
  }
  
  Output <- temp_output
  return(Output)
}

# Create a Contrast Definition Row for AUC Comparisons
### driver_parms = Output of Process_Driver_Row_Main() function
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### t_0 = First time point
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### a_i = i^th decision in DTR 1
### a_ip = i^th decision in DTR 2
Create_Individual_Contrast_Row_AUC <- function(driver_parms,
                                               model_fit_params,
                                               t_0, crit_t, t_max,
                                               a_1, a_1p, a_2r=NULL, a_2rp=NULL,
                                               a_2nr=NULL, a_2nrp=NULL, scale=TRUE) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter
  
  x_1 <- (1/2)*(crit_t^2 - t_0^2) + crit_t*(t_max-crit_t)
  x_2 <- (1/2)*(t_max^2 - crit_t^2) - crit_t*(t_max-crit_t)
  scale_factor <- ifelse(scale, (1/(t_max-t_0)), 1)
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    temp_output[1,"gamma_2"] <- x_1*(a_1-a_1p)*scale_factor
    temp_output[1,"gamma_4"] <- x_2*(a_1-a_1p)*scale_factor
    temp_output[1,"gamma_5"] <- x_2*(a_2nr-a_2nrp)*scale_factor
    temp_output[1,"gamma_6"] <- x_2*((a_1*a_2nr) - (a_1p*a_2nrp))*scale_factor
  }
  
  Output <- temp_output
  return(Output)
}


# Create a Contrast Definition Row for Main Effect Comparisons
### driver_parms = Output of Process_Driver_Row_Main() function
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### test_type = Main Effect to Test (A1, A2, Int)
Create_Individual_Contrast_Row_ME <- function(driver_parms,
                                              model_fit_params, crit_t, 
                                              t_max, test_type) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    if (test_type=="A1") {
      temp_output[1,"gamma_2"] <- (crit_t)*2
      temp_output[1,"gamma_4"] <- (t_max-crit_t)*2
    } else if (test_type=="A2") {
      temp_output[1,"gamma_5"] <- (t_max-crit_t)*2
    } else if (test_type=="Int") {
      temp_output[1,"gamma_6"] <- (t_max-crit_t)*2
    }
  }
  Output <- temp_output
  return(Output)
}




# Create a Contrast Matrix for Input into hypothesis_testing() Function
### driver_parms = Output of Process_Driver_Row_Main() function
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### test_type = The outcome to test
###             Potential options: "single_t" = Single time point 
###                                   (Must input measured_t)
###                                "AUC" = Area Under Curve
###                                   (Must input t_0, t_max)
###                                "A1" = Main Effect A1
###                                "A2" = Main Effect A2
###                                "Int" = A1/A2 interaction
### t_0 = First time point
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### comps = Which comparisons to run - Defaults to all 6 pairwise comparisons
###         Must be formatted as a list of two dataframes (with same number of rows)
###           Two columns in each dataframe (Correspond to (A1,A2))
###           We then run comparisons of DTR_1 vs DTR_2 with the rows of 
###           comps[[1]] defining DTR_1 and the rows of comps[[2]] defining DTR_2.
###           E.g. If comps[[1]][i,]=[1,1], comps[[2]][i,]=[-1,1], then the ith
###           comparison we'd make is (1,1) vs (-1,1) (Y(1,1)-Y(-1,1))
#### Returns list:
####    [["contrast_mat"]] = Contrast matrix for direct input to hypothesis_testing()
####    [["labels"]] = Dataframe of labels to append to hypothesis_testing() output
Create_Contrast_Mat <- function(driver_parms, model_fit_params, 
                                test_type, crit_t,
                                t_0=NULL, t_max=NULL,
                                measured_t=NULL, comps=NULL) {
  # Default to all comparisons
  if (is.null(comps)) {
    comps <- list(data.frame(), data.frame())
    SMART_structure_df <- Define_SMART_Pathways(driver_parms[["SMART_structure"]])
    for (i in (1:(nrow(SMART_structure_df)-1))) {
      for (ii in ((i+1):nrow(SMART_structure_df))) {
        comps[[1]] <- rbind(comps[[1]], data.frame(a_1=SMART_structure_df[i,"a_1"],
                                                   a_2r=SMART_structure_df[i,"a_2r"],
                                                   a_2nr=SMART_structure_df[i,"a_2nr"]))
        comps[[2]] <- rbind(comps[[2]], data.frame(a_1=SMART_structure_df[ii,"a_1"],
                                                   a_2r=SMART_structure_df[ii,"a_2r"],
                                                   a_2nr=SMART_structure_df[ii,"a_2nr"]))
      }
    }
  }
  temp_out <- data.frame()
  if (test_type=="single_t") {
    for (row in (1:nrow(comps[[1]]))) {
      temp_out <- rbind(temp_out,
                        Create_Individual_Contrast_Row_Single_t(driver_parms=driver_parms,
                                                                model_fit_params=model_fit_params, 
                                                                t=measured_t, 
                                                                crit_t=crit_t,
                                                                a_1=comps[[1]][row, "a_1"], 
                                                                a_1p=comps[[2]][row, "a_1"],
                                                                a_2r=comps[[1]][row, "a_2r"], 
                                                                a_2rp=comps[[2]][row, "a_2r"],
                                                                a_2nr=comps[[1]][row, "a_2nr"], 
                                                                a_2nrp=comps[[2]][row, "a_2nr"]))
    }
    comp_type <- paste0("E[Y_", measured_t, "]")
  } else if (test_type=="AUC") {
    for (row in (1:nrow(comps[[1]]))) {
      temp_out <- rbind(temp_out,
                        Create_Individual_Contrast_Row_AUC(driver_parms=driver_parms,
                                                           model_fit_params=model_fit_params,
                                                           t_0=t_0, 
                                                           crit_t=crit_t, 
                                                           t_max=t_max,
                                                           a_1=comps[[1]][row, "a_1"], 
                                                           a_1p=comps[[2]][row, "a_1"],
                                                           a_2r=comps[[1]][row, "a_2r"], 
                                                           a_2rp=comps[[2]][row, "a_2r"],
                                                           a_2nr=comps[[1]][row, "a_2nr"], 
                                                           a_2nrp=comps[[2]][row, "a_2nr"]))
    }
    comp_type <- "AUC"
  } else if (test_type %in% c("A1", "A2", "Int")) {
    temp_out <- Create_Individual_Contrast_Row_Proto_ME(driver_parms=driver_parms,
                                                        model_fit_params=model_fit_params,
                                                        crit_t=crit_t, 
                                                        t_max=t_max,
                                                        test_type=test_type)
    comp_type <- test_type
  }
  temp_out$zero <- 0
  if (test_type %in% c("single_t", "AUC")) {
    if (driver_parms[["SMART_structure"]]=="prototypical") {
      dtr_1_labels <- paste0("d_", 
                             ifelse(comps[[1]][, "a_1"]==1, 1, "m1"),
                             "_",
                             ifelse(comps[[1]][, "a_2nr"]==1, 1, "m1"))
      dtr_2_labels <- paste0("d_", 
                             ifelse(comps[[2]][, "a_1"]==1, 1, "m1"),
                             "_",
                             ifelse(comps[[2]][, "a_2nr"]==1, 1, "m1"))
      label_out <- data.frame(type=comp_type,
                              comparison=paste0("(", comps[[1]][, "a_1"], ",", 
                                                comps[[1]][, "a_2nr"], 
                                                ") vs (",
                                                comps[[2]][, "a_1"], ",", 
                                                comps[[2]][, "a_2nr"], ")"),
                              dtr_1=dtr_1_labels, 
                              dtr_2=dtr_2_labels)
    }
  } else if (test_type %in% c("A1", "A2", "Int")) {
    label_out <- data.frame(type=comp_type,
                            comparison=c("Main Effect"),
                            dtr_1=".",
                            dtr_2=".")
  }
  
  Output <- list(contrast_mat=temp_out, labels=label_out)
  return(Output)
}



# Define All Contrasts for Analysis
### driver_parms = Output of Process_Driver_Row_Main() function
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### test_types = The outcomes to test
###              Potential options: "single_t" = End of Study
###                                "AUC" = Area Under Curve
###                                   (Must input t_0, t_max)
###                                "A1" = Main Effect A1
###                                "A2" = Main Effect A2
###                                "Int" = A1/A2 interaction
### t_0 = First time point
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### comps = Which comparisons to run - Defaults to all 6 pairwise comparisons
###         Must be formatted as a list of two dataframes (with same number of rows)
###           Two columns in each dataframe (Correspond to (A1,A2))
###           We then run comparisons of DTR_1 vs DTR_2 with the rows of 
###           comps[[1]] defining DTR_1 and the rows of comps[[2]] defining DTR_2.
###           E.g. If comps[[1]][i,]=[1,1], comps[[2]][i,]=[-1,1], then the ith
###           comparison we'd make is (1,1) vs (-1,1) (Y(1,1)-Y(-1,1))
#### Returns list of lists (outer list indexed by outcome):
####    [["contrast_mat"]] = Contrast matrix for direct input to hypothesis_testing()
####    [["labels"]] = Dataframe of labels to append to hypothesis_testing() output
Define_All_Contrasts <- function(driver_parms, model_fit_params, test_types,
                                 t_0, crit_t, t_max, comps=NULL) {
  temp_output <- NULL
  for (test_type in test_types) {
    temp_output[[test_type]] <- Create_Contrast_Mat(driver_parms=driver_parms, 
                                                    model_fit_params=model_fit_params, 
                                                    test_type=test_type, 
                                                    crit_t=1,
                                                    t_0=t_0, 
                                                    t_max=2,
                                                    measured_t=t_max, 
                                                    comps=comps)
  }
  Output <- temp_output
  return(Output)
}





# Conduct All Hypothesis Tests
### driver_parms = Output of Process_Driver_Row_Main() function
### model_output = Output of solve_SMART_Multilayer() function
### contrast_output = Output of Define_All_Contrasts() function
### alpha = Level of test
### use_t = Boolean indicating whether to use t-distribution
Conduct_Hypothesis_Tests <- function(driver_parms, model_output, 
                                     contrast_output, alpha, use_t) {
  q <- ifelse(driver_parms[["SMART_structure"]]=="prototypical", 7, NA)
  temp_output <- NULL
  for (test_type in names(contrast_output)) {
    contrasts <- contrast_output[[test_type]][["contrast_mat"]]
    temp_output[[test_type]] <- hypothesis_testing(result=model_output,
                                                   aimed_test=contrasts,
                                                   q=q,
                                                   alpha=alpha,
                                                   use_t=use_t)
    contrast_labels <- contrast_output[[test_type]][["labels"]]
    temp_output[[test_type]][,colnames(contrast_labels)] <- contrast_labels
  }
  Output <- temp_output
  return(Output)
}



# Pull Inference Parameters for All Fits
### driver_parms = Output of Process_Driver_Row_Main() function
### fit_settings_df = Dataframe of fit settings
Pull_Inference_Parms_All_Fits <- function(driver_parms, fit_settings_df) {
  temp_output <- NULL
  if (driver_parms[["n_fit"]]>0) {
    for (i in (1:driver_parms[["n_fit"]])) {
      fit_setting <- driver_parms[[paste0("fit_setting_", i)]]
      fit_setting_parms <- 
        fit_settings_df[fit_settings_df$fit_setting==fit_setting, ]
      temp_output[[paste0("fit_", i)]] <- list(
        alpha=fit_setting_parms$alpha[1],
        use_t=(fit_setting_parms$use_t[1]==1)
      )
    }
  }
  Output <- temp_output
  return(Output)
}



# Conduct Raw Inference for a Simulation Run (Single Row of Main Driver)
### simulation_output = Output of Execute_Driver_Row() function
### driver_parms = Output of Process_Driver_Row_Main() function
### inference_parms = Output of Pull_Inference_Parms_All_Fits() function
### alt_fit_settings_dfs = List of dataframe of driver sheets of alternate model 
###                        fit settings (indexed by alternate fit types)
### test_types = The outcomes to test
###              Potential options: "single_t" = End of Study
###                                "AUC" = Area Under Curve
###                                   (Must input t_0, t_max)
###                                "A1" = Main Effect A1
###                                "A2" = Main Effect A2
###                                "Int" = A1/A2 interaction
### t_0 = First time point
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### comps = Which comparisons to run - Defaults to all 6 pairwise comparisons
###         Must be formatted as a list of two dataframes (with same number of rows)
###           Two columns in each dataframe (Correspond to (A1,A2))
###           We then run comparisons of DTR_1 vs DTR_2 with the rows of 
###           comps[[1]] defining DTR_1 and the rows of comps[[2]] defining DTR_2.
###           E.g. If comps[[1]][i,]=[1,1], comps[[2]][i,]=[-1,1], then the ith
###           comparison we'd make is (1,1) vs (-1,1) (Y(1,1)-Y(-1,1))
Conduct_Inference_for_Sim_Run <- function(simulation_output, driver_parms, 
                                          inference_parms, 
                                          alt_fit_settings_dfs=NULL,
                                          test_types=c("single_t", "AUC"),
                                          t_0=0, crit_t=1, t_max=2, comps=NULL) {
  temp_output <- NULL
  for (fit in names(simulation_output)) {
    #Determine if an alternate fit
    alt_fit <- strsplit(x=fit, split="fit_")[[1]][1]=="alt_"
    fit_id <- strsplit(x=fit, split="fit_")[[1]][2]
    
    if (!alt_fit) {
      # Want a different dataframe of fits for each outcome type
      for (test_type in test_types) {
        temp_output[[fit]][[test_type]] <- data.frame()
      }
      
      fit_inf_parms <- inference_parms[[fit]]
      for (iter in (1:length(simulation_output[[fit]]))) {
        #Construct contrasts and run hypothesis testing
        iter_data <- simulation_output[[fit]][[iter]]
        contrast_output <- Define_All_Contrasts(driver_parms=driver_parms, 
                                                model_fit_params=iter_data[["summary_paras"]], 
                                                test_types=test_types,
                                                t_0=t_0, 
                                                crit_t=crit_t, 
                                                t_max=t_max, 
                                                comps=comps)
        
        iter_inference <- Conduct_Hypothesis_Tests(driver_parms=driver_parms,
                                                   model_output=iter_data, 
                                                   contrast_output=contrast_output, 
                                                   alpha=fit_inf_parms[["alpha"]], 
                                                   use_t=fit_inf_parms[["use_t"]])
        # Store output
        for (test_type in test_types) {
          iter_inference[[test_type]]$iter <- iter
          temp_output[[fit]][[test_type]] <- rbind(temp_output[[fit]][[test_type]],
                                                   iter_inference[[test_type]])
        }
      }
    } else {
      if (driver_parms[[paste0("alt_fit_type_", fit_id)]]=="static") {
        temp_output[[fit]] <- data.frame()
        if (driver_parms[["SMART_structure"]]=="prototypical") {
          for (iter in (1:length(simulation_output[[fit]]))) {
            iter_inference <- simulation_output[[fit]][[iter]][["comparison"]]
            # Create DTR labels
            labels <- data.frame()
            for (i in (1:nrow(iter_inference))) {
              dtr_1_raw <- strsplit(iter_inference$Comparison[i], "AI")[[1]][2]
              dtr_1_label <- ifelse(dtr_1_raw=="(1,1)-", "d_1_1",
                                    ifelse(dtr_1_raw=="(1,-1)-", "d_1_m1",
                                           ifelse(dtr_1_raw=="(-1,1)-", "d_m1_1",
                                                  ifelse(dtr_1_raw=="(-1,-1)-", "d_m1_m1", NA))))
              
              dtr_2_raw <- strsplit(iter_inference$Comparison[i], "AI")[[1]][3]
              dtr_2_label <- ifelse(dtr_2_raw=="(1,1)", "d_1_1",
                                    ifelse(dtr_2_raw=="(1,-1)", "d_1_m1",
                                           ifelse(dtr_2_raw=="(-1,1)", "d_m1_1",
                                                  ifelse(dtr_2_raw=="(-1,-1)", "d_m1_m1", NA))))
              labels <- rbind(labels, data.frame(dtr_1=dtr_1_label,
                                                 dtr_2=dtr_2_label))
            }
            iter_inference$iter <- iter
            iter_inference[,colnames(labels)] <- labels
            temp_output[[fit]] <- rbind(temp_output[[fit]],
                                        iter_inference)
          }
        }
      }
    }
  }
  Output <- temp_output
  return(Output)
}




# Derive True Differences for Static Model Fits
### driver_parms = Output of Process_Driver_Row_Main() function
### raw_fit_inference = Element of Conduct_Inference_for_Sim_Run() output corresponding
###                     to fit in question
### outcome_type = How static data was rolled up ("EOS" or "Total")
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
#### Returns copy of raw_sim_inference with true values attached
Derive_True_Diffs_Alt_Static <- function(driver_parms, raw_fit_inference, outcome_type, 
                                         cp_settings, pre_r_marg_parms) {
  pathway_structure <- Define_SMART_Pathways(SMART_structure=driver_parms[["SMART_structure"]])
  # Grab response probabilities and calculate marginal means
  # Should be same pre-r settings for all paths with same first treatment
  # Will grab paths for responders
  dtr_means <- data.frame()
  pathway_structure$p_r <- NA
  for (row in rownames(pathway_structure)) {
    r_path <- pathway_structure[row, "path_r"]
    nr_path <- pathway_structure[row, "path_nr"]
    var_str <- cp_settings[cp_settings$pathway==r_path, "pre_r_var_str"]
    pathway_structure[row, "p_r"] <- 
      pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==var_str, "p_r"]
    # Store like this to help with debugging
    p_r <- pathway_structure[row, "p_r"]
    
    mean_cols <- colnames(cp_settings)[grepl("mean_", names(cp_settings))]
    means_r <- cp_settings[cp_settings$pathway==r_path, mean_cols]
    means_nr <- cp_settings[cp_settings$pathway==nr_path, mean_cols]
    
    row_means <- p_r*means_r + (1-p_r)*means_nr
    row_means$dtr <- pathway_structure[row, "dtr"]
    dtr_means <- rbind(dtr_means, row_means[,c("dtr", mean_cols)])
  }
  # Isolate the actual mean of interest
  if (outcome_type=="EOS") {dtr_means$primary_mean <- dtr_means$mean_2
  } else if (outcome_type=="Total") {
    dtr_means$primary_mean <- dtr_means$mean_0 + dtr_means$mean_1 + dtr_means$mean_2
  }
  # Merge on primary means
  dtr_comps <- raw_fit_inference[,c("dtr_1", "dtr_2")]
  dtr_comps$row <- rownames(raw_fit_inference)
  dtr_comps <- merge(dtr_comps, dtr_means[,c("dtr", "primary_mean")], 
                     by.x="dtr_1", by.y="dtr")
  colnames(dtr_comps)[which(names(dtr_comps) == "primary_mean")] <- "primary_mean_1"
  
  dtr_comps <- merge(dtr_comps, dtr_means[,c("dtr", "primary_mean")], 
                     by.x="dtr_2", by.y="dtr")
  colnames(dtr_comps)[which(names(dtr_comps) == "primary_mean")] <- "primary_mean_2"
  dtr_comps$true_value <- dtr_comps$primary_mean_1-dtr_comps$primary_mean_2
  
  temp_output <- merge(raw_fit_inference, dtr_comps[,c("row", "true_value")],
                       by.x="row.names", by.y="row")
  rownames(temp_output) <- as.numeric(temp_output$Row.names)
  Output <- 
    temp_output[order(as.numeric(rownames(temp_output))),-which(names(temp_output) %in% c("Row.names"))]
  return(Output)
  
}





# Derive True Values for Estimands (And Bias/Coverage)
### driver_parms = Output of Process_Driver_Row_Main() function
### raw_sim_inference = Output of Conduct_Inference_for_Sim_Run() function
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### alt_fit_settings_dfs = List of dataframe of driver sheets of alternate model 
###                        fit settings (indexed by alternate fit types)
#### Returns copy of raw_sim_inference with true values attached
Derive_True_Diffs <- function(driver_parms, raw_sim_inference,
                              cp_settings, pre_r_marg_parms, 
                              alt_fit_settings_dfs) {
  temp_output <- raw_sim_inference
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    # Calculate_Mean_Parms_Prototypical() from "Generate_All_Simulation_Parameters.R"
    # code
    mean_parms <- Calculate_Mean_Parms_Prototypical(cp_settings=cp_settings,
                                                    pre_r_marg_parms=pre_r_marg_parms)
  }
  
  for (fit in names(raw_sim_inference)) {
    #Determine if an alternate fit
    alt_fit <- strsplit(x=fit, split="fit_")[[1]][1]=="alt_"
    fit_id <- strsplit(x=fit, split="fit_")[[1]][2]
    if (!alt_fit) {
      for (test_type in names(raw_sim_inference[[fit]])) {
        inference_results <- raw_sim_inference[[fit]][[test_type]]
        true_diff_eqs <- inference_results$Hypothesis_Estimand
        for (var in names(mean_parms)) {
          # Create a character version of the equation defined by Hypothesis_Estimand
          # by replacing variable names with the actual value
          true_diff_eqs <- gsub(pattern=var, 
                                replacement=as.character(mean_parms[[var]]), 
                                x=true_diff_eqs)
        }
        # Add a column to inference results for true values by evaluating
        # the true difference equation
        temp_output[[fit]][[test_type]]$true_value <- 
          sapply(true_diff_eqs, function(x) eval(parse(text=x)))
        
        # Calculate bias and coverage
        temp_output[[fit]][[test_type]]$bias <- 
          temp_output[[fit]][[test_type]]$Estimate - temp_output[[fit]][[test_type]]$true_value
        
        temp_output[[fit]][[test_type]]$coverage <- 
          ifelse((temp_output[[fit]][[test_type]]$true_value > temp_output[[fit]][[test_type]]$CI.Lower)&
                   (temp_output[[fit]][[test_type]]$true_value < temp_output[[fit]][[test_type]]$CI.Higher),
                 1, 0)
        
        temp_output[[fit]][[test_type]]$ci_width <- 
          temp_output[[fit]][[test_type]]$CI.Higher - temp_output[[fit]][[test_type]]$CI.Lower
      }
    } else {
      if (driver_parms[[paste0("alt_fit_type_", fit_id)]]=="static") {
        alt_fit_setting <- driver_parms[[paste0("alt_fit_setting_", fit_id)]]
        alt_fit_parms <- 
          alt_fit_settings_dfs[["static"]][alt_fit_settings_dfs[["static"]]$alt_fit_setting==alt_fit_setting, ]
        outcome_type <- alt_fit_parms$outcome_type[1]
        
        temp_output[[fit]] <- 
          Derive_True_Diffs_Alt_Static(driver_parms=driver_parms, 
                                       raw_fit_inference=raw_sim_inference[[fit]], 
                                       outcome_type=outcome_type, 
                                       cp_settings=cp_settings, 
                                       pre_r_marg_parms=pre_r_marg_parms)
        
        # Calculate bias and coverage
        temp_output[[fit]]$bias <- 
          temp_output[[fit]]$Estimate - temp_output[[fit]]$true_value
        
        temp_output[[fit]]$coverage <- 
          ifelse((temp_output[[fit]]$true_value > temp_output[[fit]]$CI.Lower)&
                   (temp_output[[fit]]$true_value < temp_output[[fit]]$CI.Higher),
                 1, 0)
        
        temp_output[[fit]]$ci_width <- 
          temp_output[[fit]]$CI.Higher - temp_output[[fit]]$CI.Lower
        
      }
    }
  }
  Output <- temp_output
  return(Output)
}








# Conduct Inference for an Entire Driver Run
### path = Path where simulation results are stored
### test_types = Primary test types to run
###              Acceptable: single_t (EOS), AUC, A1/A2/Int Main effects
Conduct_Inference_for_Driver_Run <- function(path, test_types) {
  full_driver <- readRDS(paste0(path, "/Driver"))
  sims_run <- list.files(path, pattern = paste0("^", "sim"), all.files = TRUE) 
  temp_output <- NULL
  for (sim in sims_run) {
    sim_id <- strsplit(sim, split="_")[[1]][2]
    # Pull setting and parameter information for row run
    driver_row <- full_driver[["driver"]][full_driver[["driver"]]$sim_label==sim_id,]
    driver_parms <- Process_Driver_Row_Main(driver_row=driver_row)
    cp_settings_df <- full_driver[["cp_settings"]]
    cp_settings <- 
      cp_settings_df[cp_settings_df$cond_parm_setting==driver_parms[["cond_parm_setting"]],]
    raw_sim_output <- readRDS(paste0(path, "/", sim))
    inference_parms <- 
      Pull_Inference_Parms_All_Fits(driver_parms=driver_parms, 
                                    fit_settings_df=full_driver[["fit_settings"]])
    # Conduct raw inference
    raw_sim_inference <- Conduct_Inference_for_Sim_Run(simulation_output=raw_sim_output, 
                                                       driver_parms=driver_parms,
                                                       inference_parms=inference_parms)
    # Calculate and merge on the true values
    temp_output[[sim]] <- 
      Derive_True_Diffs(driver_parms=driver_parms, 
                        raw_sim_inference=raw_sim_inference,
                        cp_settings=cp_settings,
                        pre_r_marg_parms=full_driver[["var_parms_pre_r"]],
                        alt_fit_settings_dfs=full_driver[["alt_fit_settings"]])
    
  }
  Output <- temp_output
  return(Output)
}

