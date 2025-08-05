###############################
###############################
###############################
### MODEL_FITTING_WRAPPER.R ###
###############################
###############################
###############################


#### PROGRAMMER: GABRIEL DURHAM (GJD)
#### DATE CREATED: 05 FEB 2024
#### EDIT HISTORY: 12 JUL 2025 (GJD) - Edited to remove extraneous functions
####                  and add more documentation
####
#### PURPOSE: This code contains functions that help ease model fitting.
####          The Fit_cSMART_Data_Prototypical() function is a final wrapper
####          which facilitates implementation of the piecewise linear marginal 
####          mean model outlined in: https://arxiv.org/abs/2503.08987



# Create a Matrix Indicating Consistency with DTR
### in_data = Input dataset, should have columns t (time) and a_1 (-1/1) and 
###           a_2(NA/-1/1). NAs in a_2 should correspond with responders
### cluster_var = Name of cluster variable
Create_Consistency_Mat <- function(in_data, cluster_var) {
  consistency_mat <- data.frame(d_1_1=numeric(), d_1_m1=numeric(),
                                d_m1_1=numeric(), d_m1_m1=numeric())
  clusters <- unique(in_data[[cluster_var]])
  for (id in clusters) {
    cluster_data <- in_data[in_data[[cluster_var]]==id,]
    a_1 <- cluster_data$a_1[1]
    a_2 <- cluster_data$a_2[1]
    
    # Code Consistent DTRs
    c_d_1_1 <- ifelse(((a_1==1)&(a_2==1)), 1, 0)
    c_d_1_m1 <- ifelse(((a_1==1)&(a_2==-1)), 1, 0)
    c_d_m1_1 <- ifelse(((a_1==-1)&(a_2==1)), 1, 0)
    c_d_m1_m1 <- ifelse(((a_1==-1)&(a_2==-1)), 1, 0)
    
    consistency_mat <- rbind(consistency_mat, 
                             data.frame(d_1_1=c_d_1_1,
                                        d_1_m1=c_d_1_m1,
                                        d_m1_1=c_d_m1_1,
                                        d_m1_m1=c_d_m1_m1))
  }
  rownames(consistency_mat) <- clusters
  
  Output <- consistency_mat
  return(Output)
}


# Create a Consistency Matrix for a Prototypical SMART
### in_data = Data from SMART
### cluster_var = Cluster ID variable name
### a_1_var = First treatment variable name
### a_2_var = Second treatment variable name (NA for responders)
Create_Consistency_Mat_Prototypical <- function(in_data, cluster_var, 
                                                a_1_var, a_2_var) {
  temp_output <- data.frame(d_1_1=numeric(), d_1_m1=numeric(),
                            d_m1_1=numeric(), d_m1_m1=numeric())
  clusters <- unique(in_data[[cluster_var]])
  for (id in clusters) {
    cluster_data <- in_data[in_data[[cluster_var]]==id,]
    a_1 <- cluster_data[1,a_1_var]
    a_2 <- ifelse(is.na(cluster_data[1,a_2_var]), 0, cluster_data[1,a_2_var])
    r <- ifelse(is.na(cluster_data[1,a_2_var]), 1, 0)
    
    # Code Consistent DTRs
    c_d_1_1 <- ifelse(((a_1==1)&(r==1))|((a_1==1)&(a_2==1)), 1, 0)
    c_d_1_m1 <- ifelse(((a_1==1)&(r==1))|((a_1==1)&(a_2==-1)), 1, 0)
    c_d_m1_1 <- ifelse(((a_1==-1)&(r==1))|((a_1==-1)&(a_2==1)), 1, 0)
    c_d_m1_m1 <- ifelse(((a_1==-1)&(r==1))|((a_1==-1)&(a_2==-1)), 1, 0)
    
    temp_output <- rbind(temp_output,
                         data.frame(d_1_1=c_d_1_1,
                                    d_1_m1=c_d_1_m1,
                                    d_m1_1=c_d_m1_1,
                                    d_m1_m1=c_d_m1_m1))
  }
  Output <- temp_output
  return(Output)
}













# Create a Contrast Matrix for Input into hypothesis_testing() Function
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
Create_Contrast_Mat <- function(model_fit_params, test_type, crit_t,
                                t_0=NULL, t_max=NULL,
                                measured_t=NULL, comps=NULL) {
  # Default to all comparisons
  if (is.null(comps)) {
    comps <- list(
      as.data.frame(matrix(
        c(1,1,
          1,1,
          1,1,
          1,-1,
          1,-1,
          -1,1), 
        nrow=6, byrow=TRUE)),
      as.data.frame(matrix(
        c(1,-1,
          -1,1,
          -1,-1,
          -1,1,
          -1,-1,
          -1,-1), 
        nrow=6, byrow=TRUE))
    )
  }
  temp_out <- data.frame()
  if (test_type=="single_t") {
    for (row in (1:nrow(comps[[1]]))) {
      temp_out <- rbind(temp_out,
                        Create_Individual_Contrast_Row_Single_t(model_fit_params=model_fit_params, 
                                                                t=measured_t, 
                                                                crit_t=crit_t,
                                                                a_1=comps[[1]][row, 1], 
                                                                a_1p=comps[[2]][row, 1],
                                                                a_2=comps[[1]][row, 2], 
                                                                a_2p=comps[[2]][row, 2]))
    }
    comp_type <- paste0("E[Y_", measured_t, "]")
  } else if (test_type=="AUC") {
    for (row in (1:nrow(comps[[1]]))) {
      temp_out <- rbind(temp_out,
                        Create_Individual_Contrast_Row_AUC(model_fit_params=model_fit_params,
                                                           t_0=t_0, 
                                                           crit_t=crit_t, 
                                                           t_max=t_max,
                                                           a_1=comps[[1]][row, 1], 
                                                           a_1p=comps[[2]][row, 1],
                                                           a_2=comps[[1]][row, 2], 
                                                           a_2p=comps[[2]][row, 2]))
    }
    comp_type <- "AUC"
  } else if (test_type %in% c("A1", "A2", "Int")) {
    temp_out <- Create_Individual_Contrast_Row_Proto_ME(model_fit_params=model_fit_params,
                                                        crit_t=crit_t, 
                                                        t_max=t_max,
                                                        test_type=test_type)
    comp_type <- test_type
  } else if (test_type=="second_stg_slope") {
    for (row in (1:nrow(comps[[1]]))) {
      temp_out <- 
        rbind(temp_out,
              Create_Individual_Contrast_Row_Slope2(model_fit_params=model_fit_params, 
                                                    a_1=comps[[1]][row, 1], 
                                                    a_1p=comps[[2]][row, 1],
                                                    a_2=comps[[1]][row, 2], 
                                                    a_2p=comps[[2]][row, 2]))
    }
    comp_type <- "Second Stage Slope"
  }
  temp_out$zero <- 0
  if (test_type %in% c("single_t", "AUC", "second_stg_slope")) {
    label_out <- data.frame(type=comp_type,
                            comparison=paste0("(", comps[[1]][, 1], ",", 
                                              comps[[1]][, 2], 
                                              ") vs (",
                                              comps[[2]][, 1], ",", 
                                              comps[[2]][, 2], ")"))
  } else if (test_type %in% c("A1", "A2", "Int")) {
    label_out <- data.frame(type=comp_type,
                            comparison=c("Main Effect"))
  }
  
  Output <- list(contrast_mat=temp_out, labels=label_out)
  return(Output)
}

# Create a Contrast Definition Row for AUC Comparisons
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### t_0 = First time point
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### a_i = i^th decision in DTR 1
### a_ip = i^th decision in DTR 2
Create_Individual_Contrast_Row_AUC <- function(model_fit_params,
                                               t_0, crit_t, t_max,
                                               a_1, a_1p, 
                                               a_2=NULL, a_2p=NULL, scale=TRUE) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter
  
  x_1 <- (1/2)*(crit_t^2 - t_0^2) + crit_t*(t_max-crit_t)
  x_2 <- (1/2)*(t_max^2 - crit_t^2) - crit_t*(t_max-crit_t)
  scale_factor <- ifelse(scale, (1/(t_max-t_0)), 1)
  
  temp_output[1,"gamma_2"] <- x_1*(a_1-a_1p)*scale_factor
  temp_output[1,"gamma_4"] <- x_2*(a_1-a_1p)*scale_factor
  temp_output[1,"gamma_5"] <- x_2*(a_2-a_2p)*scale_factor
  temp_output[1,"gamma_6"] <- x_2*((a_1*a_2) - (a_1p*a_2p))*scale_factor
  
  Output <- temp_output
  return(Output)
}



# Create a Contrast Definition Row for Main Effect Comparisons
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### crit_t = Critical t value (time of second decision). Default is 1
### t_max = Maximum time point
### test_type = Main Effect to Test (A1, A2, Int)
Create_Individual_Contrast_Row_Proto_ME <- function(model_fit_params, crit_t, 
                                                    t_max, test_type) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter
  
  if (test_type=="A1") {
    temp_output[1,"gamma_2"] <- (crit_t)*2
    temp_output[1,"gamma_4"] <- (t_max-crit_t)*2
  } else if (test_type=="A2") {
    temp_output[1,"gamma_5"] <- (t_max-crit_t)*2
  } else if (test_type=="Int") {
    temp_output[1,"gamma_6"] <- (t_max-crit_t)*2
  }
  
  Output <- temp_output
  return(Output)
}

# Create Design Matrices for Marginal Mean Model of Prototypical SMART
## Partition into coefficients of 1, a_1, a_2 and a_1*a_2
### in_data = Data to analyze
### covar_columns = Names of columns to include as covariates
### time_column = Name of time column
### crit_t = Critical t value (time of second decision). Default is 1
Create_Partial_Design_Mats_Prototypical <- function(in_data, covar_columns=NULL, 
                                                    time_column, crit_t) {
  N <- nrow(in_data)
  in_data$t <- in_data[[time_column]]
  # Create some time placeholder variables 
  t_diff <- in_data$t-crit_t
  t_ind <- ifelse(in_data$t>crit_t, 1, 0)
  
  gamma_0 <- integer(N)+1
  gamma_1 <- (1-t_ind)*in_data$t + (t_ind*crit_t)
  #gamma_2 <- ((1-t_ind)*in_data$t + (t_ind*crit_t))*in_data$a_1
  gamma_2 <- ((1-t_ind)*in_data$t + (t_ind*crit_t))
  gamma_3 <- t_ind*t_diff
  #gamma_4 <- t_ind*t_diff*in_data$a_1
  gamma_4 <- t_ind*t_diff
  #gamma_5 <- t_ind*t_diff*in_data$a_2
  gamma_5 <- t_ind*t_diff
  #gamma_6 <- t_ind*t_diff*in_data$a_1*in_data$a_2
  gamma_6 <- t_ind*t_diff

  
  Output <- NULL
  if (!is.null(covar_columns)) {
    covars <- as.matrix(in_data[,covar_columns])
    Output[["ones"]] <- as.matrix(cbind(gamma_0, gamma_1, gamma_3, covars))
  } else {Output[["ones"]] <- as.matrix(cbind(gamma_0, gamma_1, gamma_3))}
  Output[["a_1"]] <- as.matrix(cbind(gamma_2, gamma_4))
  Output[["a_2"]] <- as.matrix(gamma_5, ncol=1)
  colnames(Output[["a_2"]]) <- c("gamma_5")
  Output[["a_1*a_2"]] <- as.matrix(gamma_6, ncol=1)
  colnames(Output[["a_1*a_2"]]) <- c("gamma_6")
  
  return(Output)
}



# Reformat Summary Parameters and Variance Estimator Labels
### model_fit_output = Output of solve_SMART_Multilayer() function
Reformat_Output <- function(model_fit_output) {
  temp_output_paras <- model_fit_output[["summary_paras"]]
  temp_output_var <- model_fit_output[["var_estimator"]]
  for (row in rownames(temp_output_paras)) {
    broken_string <- strsplit(temp_output_paras[row, "Parameter"], split=" * ")
    temp_output_paras[row, "Parameter"] <- broken_string[[1]][length(broken_string[[1]])]
  }
  for (i in (1:nrow(temp_output_var))) {
    rowname <- rownames(temp_output_var)[i]
    colname <- colnames(temp_output_var)[i]
    broken_string_row <- strsplit(rowname, split=" * ")
    rownames(temp_output_var)[i] <- 
      broken_string_row[[1]][length(broken_string_row[[1]])]
    broken_string_col <- strsplit(colname, split=" * ")
    colnames(temp_output_var)[i] <- 
      broken_string_col[[1]][length(broken_string_col[[1]])]
  }
  params_in_new_order <- order(temp_output_paras$Parameter)
  Output <- list(summary_paras=temp_output_paras[params_in_new_order,],
                 var_estimator=temp_output_var[params_in_new_order, params_in_new_order])
  rownames(Output[["summary_paras"]]) <- 1:nrow(Output[["summary_paras"]])
  return(Output)
}


# Create a A Contrast Definition Row for a Mean Difference at Single Time
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### t = Time we're comparing mean outcomes at
### crit_t = Critical t value (time of second decision). Default is 1
### a_i = i^th decision in DTR 1
### a_ip = i^th decision in DTR 2
Create_Individual_Contrast_Row_Single_t <- function(model_fit_params,
                                                    t, crit_t,  
                                                    a_1, a_1p, 
                                                    a_2=NULL, a_2p=NULL) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter
  t_diff <- t-crit_t
  t_ind <- ifelse(t>crit_t, 1, 0)
  temp_output[1,"gamma_2"] <- ((1-t_ind)*t + t_ind*crit_t)*(a_1-a_1p)
  temp_output[1,"gamma_4"] <- t_ind*t_diff*(a_1-a_1p)
  temp_output[1,"gamma_5"] <- t_ind*t_diff*(a_2-a_2p)
  temp_output[1,"gamma_6"] <- t_ind*t_diff*((a_1*a_2) - (a_1p*a_2p))
  
  Output <- temp_output
  return(Output)
}





# Create a A Contrast Definition Row for a Mean Difference at Single Time
### model_fit_params = summary_paras output of solve_SMART_Multilayer() function
### t = Time we're comparing mean outcomes at
### crit_t = Critical t value (time of second decision). Default is 1
### a_i = i^th decision in DTR 1
### a_ip = i^th decision in DTR 2
Create_Individual_Contrast_Row_Slope2 <- function(model_fit_params,
                                                  t=NULL, crit_t=NULL,  
                                                  a_1, a_1p, 
                                                  a_2=NULL, a_2p=NULL) {
  temp_output <- as.data.frame(matrix(0, nrow=1, ncol=nrow(model_fit_params)))
  colnames(temp_output) <- model_fit_params$Parameter

  temp_output[1,"gamma_4"] <- a_1-a_1p
  temp_output[1,"gamma_5"] <- a_2-a_2p
  temp_output[1,"gamma_6"] <- (a_1*a_2) - (a_1p*a_2p)
  
  Output <- temp_output
  return(Output)
}


# Fit Prototypical Clustered SMART Data and Run Static Outcome Tests
### in_data = Input data
### outcome_var = Outcome variable to study
### cluster_var = Name of variable indicating cluster ID
### crit_t = Critical time point (second decision)
### covars = List of covariate variable names
### time_var = Name of time variable in in_data
### hypothesis_tests = List of hypothesis tests to run
###   Eligible types: "second_stg_slope", "single_t", "AUC", "A1", "A2", "Int"
### measured_t = Time point to measure outcome differences 
###   (for running "single_t" test)
### a_1_var = First-stage treatment variable name
### a_2_var = Second-stage treatment variable name (must be NA for responders)
### diag_str = Diagonal block variance structure
###   Eligible input: "Independence", "AR1","Exchangeable","Unstructured"
### off_diag_str = Off-diagonal block variance structure
###   Eligible input: "Independence", "Exchangeable", "Unstructured"
### var_homo_across_time = Boolean indicating whether to model variance as constant
###   across time
### var_homo_across_AI = Boolean indicating whether to model variance as constant
###   across embedded adaptive intervention
### diagonal_homo_across_AI = Boolean indicating whether to model within-person
###   correlation as constant across embedded adaptive intervention
### off_diagonal_homo_across_AI = Boolean indicating whether to model between-person
###   correlation as constant across embedded adaptive intervention
### diagonal_ICC_lower_thresh = Lower censoring bound for within-person correlation 
###   components. If an estimated within-person correlation component falls 
###   below diagonal_ICC_lower_thresh, diagonal_ICC_lower_thresh is used in its place.
###   Setting diagonal_ICC_lower_thresh to -1 means no censoring will take place
###   for within-person correlation.
### off_diagonal_ICC_lower_thresh = Lower censoring bound for between-person correlation 
###   components.
###   Setting diagonal_ICC_lower_thresh and off_diagonal_ICC_lower_thresh to 0 
###   effectively implements the "Enforcing Nonnegative Correlation" finite sample
###   adjustment of Pan et al.
### dof_adjustment = 0/1 flag indicating whether to use Pan et al.'s degree of freedom
###   finite sample adjustment
### use_t = 0/1 flag indicating whether to use Pan et al.'s "Student's t" finite
###   sample adjustment
### bias_correction = 0/1 flag indicating whether to use Pan et al.'s "Bias Correction"
###   finite sample adjustment, which corrects for bias in estimating Var[thetahat].
###   All finite sample adjustments are described, in detail, in Pan et al.
###   https://arxiv.org/abs/2405.00185
Fit_cSMART_Data_Prototypical <- function(in_data, outcome_var, cluster_var,
                                         crit_t, covars=NULL, time_var="t", 
                                         hypothesis_tests=NULL,
                                         measured_t=NULL, a_1_var="A1", a_2_var="A2",
                                         diag_str="Independence", 
                                         off_diag_str="Independence", 
                                         var_homo_across_time=FALSE,
                                         var_homo_across_AI=FALSE,
                                         diagonal_homo_across_AI=FALSE,
                                         off_diagonal_homo_across_AI=FALSE,
                                         diagonal_ICC_lower_thresh=0,
                                         off_diagonal_ICC_lower_thresh=0,
                                         dof_adjustment=0, use_t=0, bias_correction=0) {
  # Record time variables
  t_0 <- min(in_data[[time_var]])
  t_max <- max(in_data[[time_var]])
  if (is.null(measured_t)) {measured_t <- t_max}
  
  in_data <- as.data.frame(in_data)
  partial_d_mats <- Create_Partial_Design_Mats_Prototypical(in_data=in_data, 
                                                            covar_columns=covars, 
                                                            time_column=time_var, 
                                                            crit_t=crit_t)
  
  ind_mat <- data.frame(
    "ones"=integer(4)+1,
    "a_1"=c(1,1,-1,-1),
    "a_2"=c(1,-1,1,-1),
    "a_1*a_2"=c(1,-1,-1,1)
  )
  rownames(ind_mat) <- c("d_1_1", "d_1_m1", "d_m1_1", "d_m1_m1")
  consistency_mat <- Create_Consistency_Mat_Prototypical(in_data=in_data, 
                                                         cluster_var=cluster_var, 
                                                         a_1_var=a_1_var, 
                                                         a_2_var=a_2_var)
  # Assumes all individual-level units are measured at t_0
  cluster_sizes <- table(in_data[in_data[[time_var]]==t_0, cluster_var])
  
  temp_out <- solve_SMART_Multilayer(N=length(unique(in_data[[cluster_var]])),
                                     M=as.vector(cluster_sizes),
                                     max_T=length(unique(in_data[[time_var]])),
                                     D=4,
                                     Ind=ind_mat,
                                     Y=in_data[[outcome_var]],
                                     X=partial_d_mats,
                                     within=consistency_mat,
                                     weight=1/rowSums(consistency_mat),
                                     var_homo_across_time=var_homo_across_time,
                                     var_homo_across_AI=var_homo_across_AI,
                                     diagonal_structure=diag_str,
                                     diagonal_homo_across_AI=diagonal_homo_across_AI,
                                     diagonal_ICC_lower_thresh=diagonal_ICC_lower_thresh,
                                     off_diagonal_structure=off_diag_str,
                                     off_diagonal_homo_across_AI=off_diagonal_homo_across_AI,
                                     off_diagonal_ICC_lower_thresh=off_diagonal_ICC_lower_thresh,
                                     max_iter=100,
                                     dof_adjustment=dof_adjustment,
                                     use_t=use_t,
                                     bias_correction=bias_correction,
                                     verbose=0)
  
  reformatted_components <- Reformat_Output(model_fit_output=temp_out)
  temp_out[["summary_paras"]] <- reformatted_components[["summary_paras"]]
  temp_out[["var_estimator"]] <- reformatted_components[["var_estimator"]]
  hypothesis_test_output <- data.frame()
  if (!is.null(hypothesis_tests)) {
    for (test in hypothesis_tests) {
      contrast_input <- Create_Contrast_Mat(model_fit_params=temp_out[["summary_paras"]], 
                                            test_type=test, 
                                            crit_t=crit_t,
                                            measured_t=measured_t,
                                            t_0=t_0,
                                            t_max=t_max)
      hypothesis_test_results <- hypothesis_testing(result=temp_out,
                                                    aimed_test=
                                                      contrast_input[["contrast_mat"]],
                                                    alpha=0.05,
                                                    use_t=TRUE)
      
      hypothesis_test_output <- rbind(hypothesis_test_output,
                                      cbind(contrast_input[["labels"]],
                                            hypothesis_test_results))
    }
  }
  
  
  
  
  
  
  Output <- list(model_fit=temp_out, 
                 hypothesis_tests=hypothesis_test_output)
  
  return(Output)
}
