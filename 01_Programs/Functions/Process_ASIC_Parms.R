############################
############################
############################
### PROCESS_ASIC_PARMS.R ###
############################
############################
############################

#### PURPOSE: This file contains code to support processing the parameters derived
####          from ASIC which will be used to guide the simulation settings.


#### DATE CREATED:  09 SEP 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         10 SEP 2024 - Fixed uniform shifting of DTR means when
####                  adjusting for ESs. Added sample size calculator




# Check Variance Validity (Post-Response)
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### var_parm_settings_pre_r = Dataframe of marginal pre-response variance parameter settings
### var_parm_settings_post_r = Dataframe of post-response variance parameter settings
### pre_r_cond_parms = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC)
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### cluster_sizes = Size of cluster to consider
### r_model = Response model
Check_Variance_Validity <- function(cp_settings, var_parm_settings_pre_r,
                                    var_parm_settings_post_r, pre_r_cond_parms,
                                    cluster_sizes=c(2,3), r_model=1) {
  var_parm_settings_pre_r$r_model <- r_model
  embedded_dtrs <- 
    unique(Define_SMART_Pathways(SMART_structure=cp_settings[1, "SMART_structure"])$dtr)
  temp_sim_parms <- NULL
  for (n_i in cluster_sizes) {
    temp_sim_parms[[paste0("n_", n_i)]] <- NULL
    pre_r_cond_parms_n_i <- list(
      R=pre_r_cond_parms[["R"]][pre_r_cond_parms[["R"]]$n_i==n_i, ],
      NR=pre_r_cond_parms[["NR"]][pre_r_cond_parms[["NR"]]$n_i==n_i, ]
    )
    for (dtr in embedded_dtrs) {
      temp_sim_parms[[paste0("n_", n_i)]][[dtr]] <- 
        Derive_DTR_Data_n(cp_settings=cp_settings, 
                          dtr=dtr,
                          pre_r_marg_parms=var_parm_settings_pre_r,
                          pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                          post_r_cond_parms=var_parm_settings_post_r)
    }
  }
  
  sigma_2_mat_dets <- data.frame()
  for (n_i in cluster_sizes) {for (dtr in embedded_dtrs) {
    sim_parms_n_d <- temp_sim_parms[[paste0("n_", n_i)]][[dtr]]
    for (r in names(sim_parms_n_d[["sim_parms"]][["Sigma_2"]])) {
      sigma_2_mat_dets <- 
        rbind(sigma_2_mat_dets,
              data.frame(n_i=n_i,
                         dtr=dtr,
                         response=r,
                         det=det(sim_parms_n_d[["sim_parms"]][["Sigma_2"]][[r]])))
    }
  }}
  
  Output <- sigma_2_mat_dets
  return(Output)
}


# Modify Variance Components
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### var_parm_settings_pre_r = Dataframe of marginal pre-response variance parameter settings
### var_parm_settings_post_r = Dataframe of post-response variance parameter settings
### pre_r_cond_parms = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC)
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### var_validity = Output of Check_Variance_Validity() function
### verbose = Boolean indicating whether to print output
Modify_Var_Comps <- function(cp_settings, var_parm_settings_pre_r,
                             var_parm_settings_post_r, pre_r_cond_parms,
                             var_validity, verbose=FALSE) {
  post_r_var_comp_strs <- unique(cp_settings$post_r_var_str)
  post_r_var_comps <- 
    var_parm_settings_post_r[var_parm_settings_post_r$cond_param_setting_post_r 
                             %in% post_r_var_comp_strs, ]
  SMART_structure <- cp_settings[1, "SMART_structure"]
  SMART_structure_df <- Define_SMART_Pathways(SMART_structure=SMART_structure)
  var_spec <- cp_settings[1, "var_spec"]
  
  new_cp_settings <- cp_settings
  #new_post_r_var <- var_parm_settings_post_r
  dtrs_modified <- c()
  for (row in rownames(var_validity)) {
    if (var_validity[row,"det"]<=0) {
      dtr <- var_validity[row,"dtr"]
      if (!(dtr %in% dtrs_modified)) {
        if (verbose) {message("Modifying DTR: ", dtr)}
        a_1 <- SMART_structure_df[SMART_structure_df$dtr==dtr, "a_1"]
        if (SMART_structure=="prototypical") {
          new_dtr <- 
            SMART_structure_df[SMART_structure_df$a_1==a_1&
                                 SMART_structure_df$dtr!=dtr, "dtr"]
          if (var_spec %in% c("marg_r", "marg_nr")) {
            old_setting <- cp_settings[cp_settings$dtr==dtr, "post_r_var_str"]
            new_setting <- cp_settings[cp_settings$dtr==new_dtr, "post_r_var_str"]
          }
        }
        new_cp_settings[new_cp_settings$dtr==dtr, "post_r_var_str"] <- 
          new_setting
        #new_post_r_var <- 
        #  new_post_r_var[new_post_r_var$cond_param_setting_post_r!=old_setting, ]
        dtrs_modified <- c(dtrs_modified, dtr)
      }
    }
  }
  Output <- new_cp_settings
}


# Expand CP Settings for Different Types of Comparisons
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### type = Type of comparison (n, eff_size)
### eff_sizes = Effect sizes (if type=="eff_size")
### var_parms_post_r = Dataframe of post-response variance parameters
### type = Type of comparison (n, eff_size)
### eff_sizes = Effect sizes (if type=="eff_size")
### comp = Comparison to make
### shift_a1_uniformly = Whether or not to shift the response and marginal
###   means uniformly (if adjusting to different effect sizes). E.g., if we 
###   shift the mean for d_1_1, should we also shift for d_1_m1. Not doing this
###   may cause inconsistencies in the d_1_m1 structure
###   Doing so prevents creating invalid variance structures in the new settings
### notes_base = Base text for Notes column
Expand_CP_Settings <- function(cp_settings, type, 
                               eff_sizes=c(0.01, 0.2, 0.5, 0.8, 1.2), 
                               var_parms_post_r=NULL,
                               comp=c("d_1_1", "d_m1_m1"), 
                               shift_a1_uniformly=TRUE,
                               notes_base="") {
  cp_settings_m <- cp_settings[!is.na(cp_settings$dtr), ]
  cp_settings_m <- cp_settings_m[cp_settings_m$dtr!="", ]
  cp_settings_c <- cp_settings[!is.na(cp_settings$pathway), ]
  # Adjust sample size - can be done in main driver sheet
  if (type=="n") {
    temp_output <- cp_settings
    temp_output$Notes <- notes_base
    # Adjust effect size of EOS outcome
  } else if (type=="eff_size") {if (cp_settings$var_spec[1]=="marg_r") {
    var_set_1 <- cp_settings[cp_settings$dtr==comp[1], "post_r_var_str"]
    var_set_2 <- cp_settings[cp_settings$dtr==comp[2], "post_r_var_str"]
    var_1 <- 
      var_parms_post_r[var_parms_post_r$cond_param_setting_post_r==var_set_1, "sigma2_2"]
    var_2 <- 
      var_parms_post_r[var_parms_post_r$cond_param_setting_post_r==var_set_2, "sigma2_2"]
    sigma_approx <- sqrt(mean(var_1, var_2))
    
    Y_1 <- cp_settings_m[cp_settings_m$dtr==comp[1],"mean_2"]
    Y_2 <- cp_settings_m[cp_settings_m$dtr==comp[2],"mean_2"]
    y_mean <- mean(c(Y_1, Y_2))
    
    temp_output <- data.frame()
    for (delta in eff_sizes) {
      temp_cp_settings <- cp_settings
      temp_cp_settings$Notes <- paste0(notes_base, " ES: ", delta)
      temp_cp_settings$cond_parm_setting <- 
        paste0(temp_cp_settings$cond_parm_setting, "_ES_", delta)
      SMART_str_df <- 
        Define_SMART_Pathways(SMART_structure=cp_settings[1, "SMART_structure"])
      r_path_1 <- SMART_str_df[SMART_str_df$dtr==comp[1], "path_r"]
      r_path_2 <- SMART_str_df[SMART_str_df$dtr==comp[2], "path_r"]
      shift_factor <- delta*sigma_approx/2
      
      # Shift means for DTR 1 (marginal and r-conditional)
      temp_cp_settings[(!is.na(temp_cp_settings$dtr))&
                         temp_cp_settings$dtr==comp[1], "mean_2"] <- 
        y_mean + shift_factor
      temp_cp_settings[(!is.na(temp_cp_settings$pathway))&
                         temp_cp_settings$pathway==r_path_1, "mean_2"] <- 
        cp_settings_c[cp_settings_c$pathway==r_path_1, "mean_2"] + 
        ((y_mean-Y_1) + shift_factor)
      # Shift means for DTR 2 (marginal and r-conditional)
      temp_cp_settings[(!is.na(temp_cp_settings$dtr))&
                         temp_cp_settings$dtr==comp[2], "mean_2"] <- 
        y_mean - shift_factor
      temp_cp_settings[(!is.na(temp_cp_settings$pathway))&
                         temp_cp_settings$pathway==r_path_2, "mean_2"] <- 
        cp_settings_c[cp_settings_c$pathway==r_path_2, "mean_2"] + 
        ((y_mean-Y_2) - shift_factor)
      # Shifting the other marginal means (with the same A1) by the same factor
      # keeps the R/NR conditional mean difference constant, which prevents 
      # the possibility that we create an invalid variance.
      if (shift_a1_uniformly) {for (c in c(1,2)) {
        a1 <- SMART_str_df[SMART_str_df$dtr==comp[c], "a_1"]
        for (d in SMART_str_df[SMART_str_df$a_1==a1, "dtr"]) {
          if (d!=comp[c]) {
            new_shift_factor <- 
              ifelse(c==1, 
                     (y_mean-Y_1) + shift_factor,
                     (y_mean-Y_2) - shift_factor)
            temp_cp_settings[(!is.na(temp_cp_settings$dtr))&
                               temp_cp_settings$dtr==d, "mean_2"] <- 
              cp_settings_m[cp_settings_m$dtr==d,"mean_2"] + new_shift_factor
          }
        }
      }}
      
      
      temp_output <- rbind(temp_output, temp_cp_settings)
    }
  }}
  Output <- temp_output
  return(Output)
}







# Create a Processed Version of Variance Settings
### cp_settings_df = Dataframe of conditional parameter settings
### cp_setting_base = Base conditional parameter setting to process
### var_parm_settings_pre_r = Dataframe of marginal pre-response variance parameter settings
### var_parm_settings_post_r = Dataframe of post-response variance parameter settings
### pre_r_cond_parms = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC)
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### type = Type of comparison (n, eff_size)
### eff_sizes = Effect sizes (if type=="eff_size")
### comp = Comparison to make
### notes_base = Base text for Notes column
### cluster_sizes = Cluster sizes to consider
Process_Variance_Settings <- function(cp_settings_df, cp_setting_base, 
                                      var_parm_settings_pre_r, var_parm_settings_post_r,
                                      pre_r_cond_parms, type, eff_sizes, cluster_sizes,
                                      comp=c("d_1_1", "d_m1_m1"), notes_base="") {
  base_settings <- 
    cp_settings_df[cp_settings_df$cond_parm_setting==cp_setting_base, ]
  base_validity <- 
    Check_Variance_Validity(cp_settings=base_settings,
                            var_parm_settings_pre_r=var_parm_settings_pre_r,
                            var_parm_settings_post_r=var_parm_settings_post_r, 
                            pre_r_cond_parms=pre_r_cond_parms,
                            cluster_sizes=cluster_sizes)
  new_parms <- Modify_Var_Comps(cp_settings=base_settings,
                                var_parm_settings_pre_r=var_parm_settings_pre_r,
                                var_parm_settings_post_r=var_parm_settings_post_r, 
                                pre_r_cond_parms=pre_r_cond_parms,
                                var_validity=base_validity, 
                                verbose=TRUE)
  temp_out <- Expand_CP_Settings(cp_settings=new_parms, 
                                 type=type, 
                                 eff_sizes=eff_sizes, 
                                 var_parms_post_r=var_parm_settings_post_r,
                                 shift_a1_uniformly = TRUE,
                                 comp=comp, 
                                 notes_base=notes_base)
  temp_out_valid <- TRUE
  for (note in unique(temp_out$Notes)) {
    temp_validity <- 
      Check_Variance_Validity(cp_settings=temp_out[temp_out$Notes==note, ],
                              var_parm_settings_pre_r=var_parm_settings_pre_r,
                              var_parm_settings_post_r=var_parm_settings_post_r, 
                              pre_r_cond_parms=pre_r_mc_parms,
                              cluster_sizes=cluster_sizes)
    if (min(temp_validity$det)<=0) {temp_out_valid <- FALSE}
  }
  if (temp_out_valid) {print("New output consistent with data generative model.")
  } else {print("New output inconsistent with data generative model.")}

  Output <- temp_out
  return(Output)
}



# Get the Required Sample Size from a Static SMART Based off Tim NeCamp's Work
### cp_settings_df = Dataframe of conditional parameter settings
### cp_setting_base = Base conditional parameter setting to process
### var_parm_settings_pre_r = Dataframe of marginal pre-response variance parameter settings
### var_parm_settings_post_r = Dataframe of post-response variance parameter settings
### pre_r_cond_parms = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC)
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### comp = Comparison to make
### cluster_sizes = Cluster sizes to consider
### r_model = Response model used
### power = Desired power
### level = Level of test
### verbose = Whether to print out variance comparisons
Grab_Sample_Size <- function(cp_settings_df, cp_setting_base, var_parm_settings_pre_r,
                             var_parm_settings_post_r, pre_r_cond_parms,
                             cluster_sizes, eff_sizes,
                             comp=c("d_1_1", "d_m1_m1"), r_model=1,
                             power=0.8, level=0.05, verbose=FALSE) {
  cp_settings <- 
    cp_settings_df[cp_settings_df$cond_parm_setting==cp_setting_base, ]
  var_parm_settings_pre_r$r_model <- r_model
  temp_sim_parms <- NULL
  for (n_i in cluster_sizes) {
    temp_sim_parms[[paste0("n_", n_i)]] <- NULL
    pre_r_cond_parms_n_i <- list(
      R=pre_r_cond_parms[["R"]][pre_r_cond_parms[["R"]]$n_i==n_i, ],
      NR=pre_r_cond_parms[["NR"]][pre_r_cond_parms[["NR"]]$n_i==n_i, ]
    )
    for (dtr in comp) {
      temp_sim_parms[[paste0("n_", n_i)]][[dtr]] <- 
        Derive_DTR_Data_n(cp_settings=cp_settings, 
                          dtr=dtr,
                          pre_r_marg_parms=var_parm_settings_pre_r,
                          pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                          post_r_cond_parms=var_parm_settings_post_r)
    }
  }
  # Grab variance components
  rho <- c()
  parms_n1 <- temp_sim_parms[[paste0("n_", cluster_sizes[1])]]
  for (dtr in comp) {
    marg_comps <- 
      parms_n1[[dtr]][["dist"]][["marg"]]
    nr_comps <- 
      parms_n1[[dtr]][["dist"]][["cond"]][["NR"]]
    rho <- c(rho, marg_comps[["var"]][["rho_2"]])
    var_m <- marg_comps[["var"]][["sigma2_2"]]
    var_nr <- nr_comps[["var"]][["sigma2_2"]]
    mean_m <- marg_comps[["mean"]][["mean_2"]]
    mean_nr <- nr_comps[["mean"]][["mean_2"]]
    if (verbose) {
      print("***************")
      print(paste0("*** ", dtr, " ***"))
      print("***************")
      message("Marginal Variance: ", round(var_m,2))
      message("E[(Y2-mu2)^2|R=0]: ", round(var_nr + (mean_m + mean_nr)^2, 2))
    }
  }
  # Check identical rho values
  if (rho[1]!=rho[2]) {print("Warning: Non-identical rho values. 
                             Using comp[1] distribution")}
  rho <- rho[1]/parms_n1[[dtr]][["dist"]][["marg"]][["var"]][["sigma2_2"]]
  p_1 <- 
    parms_n1[[comp[1]]][["dist"]][["other_parms"]][["marg"]][["p_r"]]
  p_2 <- 
    parms_n1[[comp[2]]][["dist"]][["other_parms"]][["marg"]][["p_r"]]
  z_beta <- qnorm(p=power)
  z_alpha2 <- qnorm(p=1-(level/2))
  z_term <- 
    for (m in cluster_sizes) {for (delta in eff_sizes) {
      sample_size <- 
        ( 4*((z_beta+z_alpha2)^2)/(m*(delta^2))  ) *
        (1 + (m-1)*rho ) * 
        (1 + 0.5*(2-p_1-p_2))
      print("***************")
      print(paste0("(n_i = ", m,
                   ", ES = ", delta,
                   ") Required Sample Size: ",
                   ceiling(sample_size)))
    }}
  print("***************")
  print(paste0("(rho=", round(rho, 3), ")"))
}