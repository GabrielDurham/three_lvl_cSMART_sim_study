############################
############################
############################
### ALT_HELPER_FUNCTIONS ###
############################
############################
############################

#### PURPOSE: This file contains code to support alternate simulation runs that
####          depend on multiple time points.


#### DATE CREATED:  04 JUN 2025
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         04 JUN 2025 (GJD)
####                06 JUN 2025 (GJD) - Changed "person" to "person_id" in
####                  simulated dataset to match main version.


# Build a Variance Matrix
## Build a variance matrix with AR1 within-person correlation and exchangeable
## between-person correlation
### marg_var = Marginal variances (i.e., diagonal entries)
### rho = AR1 correlation
### icc = Between-person correlation
Build_Variance_Matrix <- function(marg_var, rho, icc, n){
  S <- diag(rep(marg_var, n))
  t <- length(marg_var)
  corr_inner <- diag(rep(1,t))
  corr_outer <- diag(rep(icc, t))
  for (i in c(1:(t-1))) {for (ii in c((i+1):t)) {
    corr_inner[i,ii] <- rho^abs(i-ii)
    corr_inner[ii,i] <- rho^abs(i-ii)
    
    corr_outer[i,ii] <- icc
    corr_outer[ii,i] <- icc
  }}
  if (n==1) {temp_corr <- corr_inner
  } else {
    new_row <- corr_inner
    for (i in (2:n)) {new_row <- cbind(new_row, corr_outer)}
    temp_corr <- new_row
    begin_row <- corr_outer
    for (row in (2:n)) {
      new_row <- cbind(begin_row, corr_inner)
      if (row<n) {for (col in ((row+1):n)) {new_row <- cbind(new_row, corr_outer)}}
      begin_row <- cbind(begin_row, corr_outer)
    }
    temp_corr <- rbind(temp_corr, new_row)
  }
  corr <- temp_corr
  Output <- (S^0.5)%*%corr%*%(S^0.5)
  return(Output)
}

# Build Distribution Matrices for a Given DTR
## Build the mean vector and variance matrices for a given DTR. Pre-R mean and
## variance, post-R conditional means and variances. Also attaches some helper
## info like marginal response probability and a scaling factor for response.
### cp_settings = Relevent rows of alternate cp_settings sheet
### dtr = Embedded dtr in question
### rho = AR1 correlation
### icc = Between-person correlation
### n = Cluster size
Build_Distribution_Mats_DTR <- function(cp_settings, dtr, rho, icc, n) {
  dtr_row <- cp_settings[cp_settings$dtr==dtr, ]
  mean_pre_r <- rep(eval(parse(text=dtr_row[1,"mean_pre_r"])), n)
  mean_post_r_R <- rep(eval(parse(text=dtr_row[1,"mean_post_r_R"])), n)
  mean_post_r_NR <- rep(eval(parse(text=dtr_row[1,"mean_post_r_NR"])), n)
  
  var_pre_r <- eval(parse(text=dtr_row[1,"variance_pre_r"]))
  var_post_r_R <- eval(parse(text=dtr_row[1,"variance_post_r_R"]))
  var_post_r_NR <- eval(parse(text=dtr_row[1,"variance_post_r_NR"]))
  
  var_matrix_pre_R <- Build_Variance_Matrix(marg_var=var_pre_r,
                                            rho=rho,
                                            icc=icc,
                                            n=n)
  var_matrix_post_R <- 
    list(R=Build_Variance_Matrix(marg_var=var_post_r_R, rho=rho, icc=icc, n=n),
         NR=Build_Variance_Matrix(marg_var=var_post_r_NR, rho=rho, icc=icc, n=n))
  # Add scaling factor for response
  scaling_factor_1 <- 
    sqrt((1/n)*(var_pre_r[length(var_pre_r)] + (n-1)*icc*var_pre_r[length(var_pre_r)]))
  p_r <- dtr_row[1, "p_r"][[1]]
  
  #Grab marginal means for inference
  t <- 0
  marg_means <- NULL
  for (mean in eval(parse(text=dtr_row[1,"mean_pre_r"]))) {
    marg_means[[paste0("mean_", t)]] <- mean
    t <- t+1
  }
  post_r_means_R <- eval(parse(text=dtr_row[1,"mean_post_r_R"]))
  post_r_means_NR <- eval(parse(text=dtr_row[1,"mean_post_r_NR"]))
  for (i in (2:length(post_r_means_R))) {
    marg_means[[paste0("mean_", t)]] <- 
      p_r*post_r_means_R[i] + (1-p_r)*post_r_means_NR[i]
    t <- t+1
  }
  
  Output <- 
    list(pre_r=list(mean=mean_pre_r, var=var_matrix_pre_R,
                    scaling_factor_1=scaling_factor_1, p_r=p_r),
         post_r=list(
           R=list(mean=mean_post_r_R, var=var_matrix_post_R[["R"]]),
           NR=list(mean=mean_post_r_NR, var=var_matrix_post_R[["NR"]])
         ),
         dist = list(marg=list(mean=marg_means))
         )
  return(Output)
}

# Build All Distribution Matrices
## Construct distribution information for all DTRs/n
### driver_parms = Output of Process_Driver_Row_Main() function
### cp_settings = Relevent rows of alternate cp_settings sheet
Build_All_Distribution_Mats <- function(driver_parms, cp_settings) {
  temp_output <- NULL
  for (n in driver_parms[["cluster_sizes"]]) {
    n_lab <- paste0("n_", n)
    temp_output[[n_lab]] <- NULL
    for (dtr in cp_settings$dtr) {
      rho <- cp_settings[cp_settings$dtr==dtr, "target_rho"][[1]]
      icc <- cp_settings[cp_settings$dtr==dtr, "target_icc"][[1]]
      
      temp_output[[n_lab]][[dtr]] <- 
        Build_Distribution_Mats_DTR(cp_settings=cp_settings,
                                    dtr=dtr,
                                    rho=rho,
                                    icc=icc,
                                    n=n)
    }
  }
  Output <- temp_output
  return(Output)
}

# Response Function 1 - Alternate
## Slightly modified version of response function 1 for the alternate runs
### marg_p_r = Probability of response
### outcome_devs = Outcome deviations (centered outcomes - by true mean)
### scaling_factor = Scaling factor (sqrt(1/n (sigma_1 + (n-1)sigma_1*icc)))
Response_Function_1_Alt <- function(marg_p_r, outcome_devs, scaling_factor){
  beta_input <- pnorm(mean(outcome_devs)/scaling_factor)
  p_response <- qbeta(beta_input, shape1=(marg_p_r/(1-marg_p_r)), shape2=1, 
                      ncp = 0, log = FALSE)
  Output <- rbinom(n=1, size=1, prob=p_response)
  return(Output)
}





# Simulate SMART Data-Alternate
## Simulate SMART data from more than three time points
### driver_parms = Output of Process_Driver_Row_Main() function
### dist_data = Output from Build_All_Distribution_Mats() function
#### Requires t_0 = 0, and crit_t being the midpoint of t_0 and t_T (with time
####    points equally spaced between the two.)
Sim_SMART_Data_Alt <- function(driver_parms, dist_data) {
  cluster_data <- Randomize_DTR_n_i(driver_parms)
  SMART_structure <- 
    Define_SMART_Pathways(SMART_structure=driver_parms[["SMART_structure"]])
  # Simulate Pre-Response Outcomes
  temp_output <- data.frame()
  for (cluster in unique(cluster_data$cluster_id)) {
    cluster_row <- cluster_data[cluster_data$cluster_id==cluster, ]
    n_i <- cluster_row[1, "n"]
    d_i <- cluster_row[1, "dtr"]
    
    a_1 <- SMART_structure[SMART_structure$dtr==d_i, "a_1"]
    a_2r <- SMART_structure[SMART_structure$dtr==d_i, "a_2r"]
    a_2nr <- SMART_structure[SMART_structure$dtr==d_i, "a_2nr"]

    dist_data_pre_r <- dist_data[[paste0("n_", n_i)]][[d_i]][["pre_r"]]
    pre_r_outcomes <- 
      mvrnorm(n=1, mu=dist_data_pre_r[["mean"]], Sigma=dist_data_pre_r[["var"]])
    
    crit_t <- driver_parms[["crit_t"]]
    
    cluster_pre_r_data <- data.frame(
      person_id = rep(1:n_i, each=(crit_t+1)),
      t = rep(0:(crit_t), n_i),
      y = pre_r_outcomes
    )
    
    # Simulate Response
    outcome_devs <- 
      cluster_pre_r_data[cluster_pre_r_data$t==crit_t, "y"] - dist_data_pre_r[["mean"]][crit_t]
    
    response <- 
      Response_Function_1_Alt(marg_p_r=dist_data_pre_r[["p_r"]], 
                              outcome_devs=outcome_devs, 
                              scaling_factor=dist_data_pre_r[["scaling_factor_1"]])
    
    if (response==1) {
      dist_data_post_r <- dist_data[[paste0("n_", n_i)]][[d_i]][["post_r"]][["R"]]
      a_2 <- a_2r
    } else {
      dist_data_post_r <- dist_data[[paste0("n_", n_i)]][[d_i]][["post_r"]][["NR"]]
      a_2 <- a_2nr
    }
    
    post_r_outcomes <- 
      mvrnorm(n=1, mu=dist_data_post_r[["mean"]], Sigma=dist_data_post_r[["var"]])
    #Drop the first one - represents crit_t
    post_r_outcomes <- 
      post_r_outcomes[which(((1:length(post_r_outcomes)) %% (crit_t+1)) != 1)]
    cluster_post_r_data <- data.frame(
      person_id = rep(1:n_i, each=(crit_t)),
      t = rep((crit_t+1):(2*crit_t), n_i),
      y = post_r_outcomes
    )
    temp_cluster_data <- rbind(cluster_pre_r_data, cluster_post_r_data)
    temp_cluster_data <- 
      temp_cluster_data[order(temp_cluster_data$person_id, temp_cluster_data$t), ]
    
    temp_cluster_data$cluster_id <- cluster
    temp_cluster_data$r <- response
    temp_cluster_data$a_1 <- a_1
    temp_cluster_data$a_2 <- a_2
    
    temp_output <- rbind(temp_output, temp_cluster_data)
  }
  Output <- temp_output
  rownames(Output) <- (1:nrow(Output))
  
  return(Output)
}