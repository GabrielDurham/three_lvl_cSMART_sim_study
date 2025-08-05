###########################
###########################
###########################
### SIMULATE SMART DATA ###
###########################
###########################
###########################

#### PURPOSE: This file contains code to support the Sim_SMART_Data() function,
####          which simulates data from a single SMART to driver specifications


#### DATE CREATED:  14 JAN 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         22 MAY 2024 (GJD)
####                    - Corrected issue in Randomize_DTR_n_i() which prevented
####                      simulation runs with only one cluster size
####                    - Corrected issue in Add_All_Covariates() which prevented
####                      simulation runs with no covariates
####                04 JUN 2024 (GJD) - Changed [["var_parms"]] label of 
####                      Obtain_All_Sim_Params() output to [["target_var"]]
####                07 JUL 2024 (GJD) - Cleaned up code to incorporate new
####                      Derive_Sim_Parms() function (replaced 
####                      Obtain_All_Sim_Params())
####                24 FEB 2025 (GJD) - Cleaned up Randomize_DTR_n_i() and created
####                      Randomize_DTR_Proto() to allow for complete randomization
####                      based DTR assignment



# Randomize DTR Assignments for a Prototypical SMART
### driver_parms = Output of Process_Driver_Row_Main() function
Randomize_DTR_Proto <- function(driver_parms) {
  dtr_list <- c("d_1_1", "d_1_m1", "d_m1_1", "d_m1_m1")
  if (driver_parms[["ra_structure"]]=="bernoulli") {
    dtr_probs <- driver_parms[["p_dtr"]]
    # Repeat until the minimum number of clusters randomized to a DTR is met
    below_min_dtr_rep <- TRUE
    while (below_min_dtr_rep) {
      dtr_assignments <- sample(x=dtr_list, 
                                size=driver_parms[["n_clusters"]], 
                                replace=TRUE, 
                                prob=dtr_probs)
      dtr_representations <- c()
      for (dtr in dtr_list) {
        dtr_representations <- c(dtr_representations, sum(dtr_assignments==dtr))
      }
      #Check whether the minimum number of clusters randomized to a DTR is met
      below_min_dtr_rep <- (min(dtr_representations)<driver_parms[["min_dtr_obs"]])
    }
  }
  if (driver_parms[["ra_structure"]]=="complete_ra") {
    dtrs_to_draw <- c()
    # Make a list of ~N DTR assignments, balanced between each DTR according to
    # p_dtr
    for (d in (1:length(dtr_list))) {
      dtr <- dtr_list[d]
      n_clusters_per_dtr <- 
        ceil(driver_parms[["n_clusters"]]*driver_parms[["p_dtr"]][d])
      dtrs_to_draw <- c(dtrs_to_draw, rep(dtr, n_clusters_per_dtr))
    }
    # Sample from the list without replacement to get roughly even DTR assignments
    dtr_assignments <- 
      sample(dtrs_to_draw, driver_parms[["n_clusters"]], replace=FALSE)
  }
  
  Output <- dtr_assignments
  return(Output)
}
  




# Randomize the DTR Assignments and Sample Sizes
### driver_parms = Output of Process_Driver_Row_Main() function
Randomize_DTR_n_i <- function(driver_parms) {
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    dtr_assignments <- Randomize_DTR_Proto(driver_parms=driver_parms)
  }
  
  size_probs <- driver_parms[["p_cluster_size"]]
  # If only one cluster size then just produce the vector without sampling
  if (length(driver_parms[["cluster_sizes"]])==1) {
    cluster_sizes <- rep(x=driver_parms[["cluster_sizes"]][1], 
                         times=driver_parms[["n_clusters"]])
  } else { 
    cluster_sizes <- sample(x=driver_parms[["cluster_sizes"]],
                            size=driver_parms[["n_clusters"]], 
                            replace=TRUE, 
                            prob=size_probs)
  }
  Output <- data.frame(cluster_id=(1:driver_parms[["n_clusters"]]),
                       n=cluster_sizes,
                       dtr=dtr_assignments)
  return(Output)
}

# Simulate Pre-Response Marginal Deviations and Response for a Single Cluster
### sim_parms = Output of Obtain_All_Sim_Params() function
### dtr = Embedded DTR of cluster
### n_i = Cluster size
Simulate_Y_0_1_R_One_Cluster <- function(driver_parms, sim_parms, dtr, n_i) {
  cluster_parms <- sim_parms[[paste0("n_", n_i)]][[dtr]]
  Sigma_0_1 <- cluster_parms[["sim_parms"]][["Sigma_01"]]
  e_0_1 <- mvrnorm(n = 1, mu=integer(2*n_i), Sigma=Sigma_0_1)
  
  e_0 <- e_0_1[(1:n_i)]
  e_1 <- e_0_1[((n_i+1):(2*n_i))]
  
  mean_helper_parms <- cluster_parms[["sim_parms"]][["mean_helper_parms"]]
  P_1 <- cluster_parms[["sim_parms"]][["P"]][["R"]][1]
  y_0 <- mean_helper_parms[["t_0"]] + e_0
  y_1 <- mean_helper_parms[["t_1"]] + P_1*y_0 + e_1
  #Grab centered outcomes and response
  y_0_c <- e_0
  y_1_c <- P_1*e_0 + e_1
  response <- 
    cluster_parms[["response_function"]](list(Y_0=y_0_c, Y_1=y_1_c))
  
  Output <- list(Y=list(Y_0=y_0, Y_1=y_1), 
                 e=list(e_0=e_0, e_1=e_1),
                 R=response)
  return(Output)
}


# Simulate a Set of Clusters and their Pre-Response Outcomes
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_parms = Output of Obtain_All_Sim_Params() function
Sim_Pre_R_Data <- function(driver_parms, sim_parms) {
  below_min_path_met <- TRUE
  path_structure <- Define_SMART_Pathways(SMART_structure=driver_parms[["SMART_structure"]])
  # Track whether each path is represented
  all_paths <- sort(unique(c(path_structure$path_r,path_structure$path_nr)))
  
  # Repeat until the minimum number of clusters randomized to a path is met
  while (below_min_path_met) {
    path_rep_df <- data.frame(path=all_paths, n_members=integer(length(all_paths)))
    cluster_assignments <- Randomize_DTR_n_i(driver_parms=driver_parms)
    cluster_assignments$r <- NA
    cluster_assignments$a_1 <- NA
    cluster_assignments$a_2 <- NA
    cluster_assignments$a_2r <- NA
    cluster_assignments$a_2nr <- NA
    pre_r_outcomes <- data.frame(cluster_id=numeric(), person_id=numeric(), 
                                 y_0=numeric(), y_1=numeric(), e_0=numeric(), 
                                 e_1=numeric(), r=numeric())
    for (row in rownames(cluster_assignments)) {
      cluster_id <- cluster_assignments[row, "cluster_id"]
      n_i <- cluster_assignments[row, "n"]
      d_i <- cluster_assignments[row, "dtr"]
      cluster_data <- Simulate_Y_0_1_R_One_Cluster(sim_parms=sim_parms, 
                                                   dtr=d_i, 
                                                   n_i=n_i)
      
      cluster_y <- cluster_data[["Y"]]
      cluster_errs <- cluster_data[["e"]]
      cluster_r <- cluster_data[["R"]]
      
      cluster_assignments[row, "r"] <- cluster_r
      dtr_structure <- sim_parms[[paste0("n_", n_i)]][[d_i]][["dtr"]]
      cluster_assignments[row, "a_1"] <- dtr_structure[["a_1"]]
      cluster_assignments[row, "a_2r"] <- dtr_structure[["a_2r"]]
      cluster_assignments[row, "a_2nr"] <- dtr_structure[["a_2nr"]]
      
      
      # Identify cluster path
      cluster_path <- ifelse(cluster_r==1,
                             path_structure[path_structure$dtr==d_i,"path_r"],
                             path_structure[path_structure$dtr==d_i,"path_nr"])
      # Update path representation tracker
      path_rep_df[path_rep_df$path==cluster_path, "n_members"] <- 
        path_rep_df[path_rep_df$path==cluster_path, "n_members"] + 1
      
      # Update pre_r_outcomes dataframe
      pre_r_outcomes <- rbind(pre_r_outcomes,
                              data.frame(cluster_id=rep(cluster_id, n_i),
                                         person_id=(1:n_i),
                                         y_0=cluster_y[["Y_0"]],
                                         y_1=cluster_y[["Y_1"]],
                                         e_0=cluster_errs[["e_0"]],
                                         e_1=cluster_errs[["e_1"]],
                                         r=rep(cluster_r, n_i)))
    }
    below_min_path_met <- (min(path_rep_df$n_members)<driver_parms[["min_path_obs"]])
  }
  cluster_assignments$a_2 <- ifelse(cluster_assignments$r==1, 
                                    cluster_assignments$a_2r, 
                                    cluster_assignments$a_2nr)
  
  Output <- list(cluster_data=cluster_assignments, outcomes=pre_r_outcomes)
  return(Output)
}


# Simulate Post-Response Outcomes (w/o Covs) for a Single Cluster of a Prototypical SMART
### cluster_id = ID of cluster whose outcomes to simulate
### sim_parms = Output of Obtain_All_Sim_Params() function
### pre_r_data = Output of Sim_Pre_R_Data() function
Sim_Y2_One_Cluster_Proto <- function(cluster_id, sim_parms, pre_r_data) {
  #Load in cluster level information
  cluster_lvl_data <- 
    pre_r_data[["cluster_data"]][pre_r_data[["cluster_data"]]$cluster_id==cluster_id,]
  d_i <- cluster_lvl_data$dtr
  n_i <- cluster_lvl_data$n
  # a_1 <- cluster_lvl_data$a_1
  # a_2nr <- cluster_lvl_data$a_2nr
  r_i <- cluster_lvl_data$r
  r_lab <- ifelse(r_i==1, "R", "NR")
  
  # Load in cluster level simulation parameters
  cluster_sim_parms <- sim_parms[[paste0("n_", n_i)]][[d_i]][["sim_parms"]]
  P <- cluster_sim_parms[["P"]][[r_lab]]
  Sigma_2 <- cluster_sim_parms[["Sigma_2"]][[r_lab]]
  # p_r <- sim_parms[[paste0("n_", n_i)]][[d_i]][["dist"]][["other_parms"]][["marg"]][["p_r]]
  mean_helper_parms <- cluster_sim_parms[["mean_helper_parms"]]
  
  # Generate Errors
  e_2 <- mvrnorm(n = 1, mu=integer(n_i), Sigma=Sigma_2)
  
  #Pull relevant outcome information
  outcomes <- 
    pre_r_data[["outcomes"]][pre_r_data[["outcomes"]]$cluster_id==cluster_id,]
  
  # Created the helper parms to avoid some of this repeated calculation
  Output <- mean_helper_parms[["t_2"]][[r_lab]] + 
    # (1 - P[2] - P[4])*gamma[["0"]] + 
    # (1 - P[4])*(gamma[["1"]] + gamma[["2"]]*a_1) + 
    # (gamma[["3"]] + gamma[["4"]]*a_1) + 
    P[2]*outcomes$y_0 + 
    P[3]*mean(outcomes$e_0) + 
    P[4]*outcomes$y_1 + 
    P[5]*mean(outcomes$e_1) + 
    # ( (gamma[["5"]] + gamma[["6"]]*a_1)*a_2nr )*( (1-r_i)/(1-p_r) ) + 
    # (lambda[["1"]] + lambda[["2"]]*a_1)*(r_i-p_r) + 
    e_2
  return(Output)
}



# Simulate Post-Response Outcomes
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_parms = Output of Obtain_All_Sim_Params() function
### pre_r_data = Output of Sim_Pre_R_Data() function
#### Returns copy of pre_r_data[["outcomes"]] with simulated y_2
Sim_Post_R_Outcomes <- function(driver_parms, sim_parms, pre_r_data) {
  temp_output <- pre_r_data[["outcomes"]]
  temp_output$y_2 <- NA
  for (id in unique(temp_output$cluster_id)) {
    if (driver_parms[["SMART_structure"]]=="prototypical") {
      cluster_y_2 <- Sim_Y2_One_Cluster_Proto(cluster_id=id, 
                                              sim_parms=sim_parms, 
                                              pre_r_data=pre_r_data)
    }
    
    temp_output[temp_output$cluster_id==id, "y_2"] <- cluster_y_2
  }
  Output <- temp_output
  return(Output)
}


# Turn Wide Data To Long (I.e., Make Single Outcome Column, Add Time Column)
### wide_output = Person-level data with outcomes y_0, y_1, y_2
Structure_Wide_Output <- function(wide_output){
  outcomes <- c("y_0", "y_1", "y_2")
  temp_list <- NULL
  for (time in c(0,1,2)) {
    time_lab <- paste0("t_", time)
    y_lab <- paste0("y_", time)
    temp_list[[time_lab]] <- 
      wide_output[,!(names(wide_output) %in% outcomes)]
    temp_list[[time_lab]]$t <- time
    temp_list[[time_lab]]$y <- wide_output[[y_lab]]
  }
  stacked_data <- rbind(temp_list[["t_0"]], temp_list[["t_1"]], temp_list[["t_2"]])
  Output <- stacked_data[order(stacked_data$cluster_id,
                               stacked_data$person_id,
                               stacked_data$t),]
  rownames(Output) <- 1:nrow(Output)
  return(Output)
}




# Simulate a Single Covariate
### in_data = Person-level dataset to add covariate to (must have "cluster_id" column)
### covar_level = Level of covariate
###           0: Individual
###           1: Cluster
### covar_dist = Distribution of covariate (accepts "uniform", "normal")
### new_name = Desired name of new covariate (default "new_covar")
#### Returns copy of in_data with new_x
Add_Single_Covariate <- function(in_data, covar_level, covar_dist, new_name="new_x"){
  temp_data <- in_data
  temp_data$temp_new_x <- NA
  clusters <- unique(temp_data$cluster_id)
  #Simulate individual-level covariate
  if (covar_level==0) {
    if (covar_dist=="uniform") {
      temp_data$temp_new_x <- runif(n=nrow(temp_data), min=-sqrt(12)/2, max=sqrt(12)/2)
    } else if (covar_dist=="normal") {
      temp_data$temp_new_x <- rnorm(n=nrow(temp_data), mean=0, sd=1)
    }
  }
  #Simulate cluster level covariate
  else if (covar_level==1) {
    if (covar_dist=="uniform") {
      for (id in clusters) {
        temp_data[temp_data$cluster_id==id, "temp_new_x"] <- 
          runif(n=1, min=-sqrt(12)/2, max=sqrt(12)/2)
      }
    } else if (covar_dist=="normal") {
      for (id in clusters) {
        temp_data[temp_data$cluster_id==id, "temp_new_x"] <- rnorm(n=1, mean=0, sd=1)
      }
    }
  }
  colnames(temp_data)[which(names(temp_data) == "temp_new_x")] <- new_name
  Output <- temp_data
  return(Output)
}


# Simulate All Person/Cluster Level Covariates
### in_data = Person-level dataset to add covariate to (must have "cluster_id" column)
### driver_parms = Output of Process_Driver_Row_Main() function
#### Returns copy of in_data with simulated covariates and a column for covariate effect
####  (i.e., X^T%*%eta)
Add_All_Covariates <- function(in_data, driver_parms) {
  temp_data <- in_data
  # Need to specify no covariate effect if no covariates
  if (driver_parms[["n_covar"]]==0) {temp_data$covar_effect <- 0
  } else {
    temp_data$covar_effect <- 0
    for (i in (1:driver_parms[["n_covar"]])) {
      cov_label <- paste0("x_", i)
      coef_label <- paste0("eta_", i)
      
      cov_level <- driver_parms[[paste0(cov_label, "_cluster_lvl")]]
      cov_dist <- driver_parms[[paste0(cov_label, "_dist")]]
      temp_data <- Add_Single_Covariate(in_data=temp_data,
                                        covar_level=cov_level,
                                        covar_dist=cov_dist,
                                        new_name=paste0("x_", i))
      # Capture covariate effects
      temp_data$covar_effect <- temp_data$covar_effect + 
        temp_data[[cov_label]]*driver_parms[[coef_label]]
    }
  }
  Output <- temp_data
  return(Output)
}


# Simulate Single SMART According to Driver Specifications
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_parms = Output of Obtain_All_Sim_Params() function
Sim_SMART_Data <- function(driver_parms, sim_parms) {
  # Randomize cluster size and DTR assignments and pre-response outcomes
  pre_response_data <- Sim_Pre_R_Data(driver_parms=driver_parms, 
                                      sim_parms=sim_parms)
  # Add post-response outcomes without covariates
  outcomes_wo_covs <- Sim_Post_R_Outcomes(driver_parms=driver_parms, 
                                          sim_parms=sim_parms, 
                                          pre_r_data=pre_response_data)
  # Add simulated covariates
  raw_outcomes <- Add_All_Covariates(in_data=outcomes_wo_covs,
                                     driver_parms=driver_parms)
  # Construct output dataset
  wide_outcomes <- data.frame(cluster_id=raw_outcomes$cluster_id,
                              person_id=raw_outcomes$person_id,
                              y_0=raw_outcomes$y_0 + raw_outcomes$covar_effect,
                              y_1=raw_outcomes$y_1 + raw_outcomes$covar_effect,
                              y_2=raw_outcomes$y_2 + raw_outcomes$covar_effect)
  wide_outcomes[,driver_parms[["covars"]]] <- raw_outcomes[,driver_parms[["covars"]]]
  
  person_data <- Structure_Wide_Output(wide_output=wide_outcomes)
  
  # Merge treatment assignment data
  treatment_assignments <- 
    pre_response_data[["cluster_data"]][, c("cluster_id", "a_1", "r", "a_2")]
  Output <- merge(x=person_data, 
                  y=treatment_assignments, 
                  by="cluster_id")
  # Reorder columns just because I think it looks nice
  cols_in_order <- c("cluster_id", "person_id", "t", "a_1", "r", "a_2", 
                     "y", driver_parms[["covars"]])
  Output <- Output[, cols_in_order]
  
  return(Output)
}