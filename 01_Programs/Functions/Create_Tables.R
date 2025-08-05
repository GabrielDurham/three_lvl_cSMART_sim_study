#####################
#####################
#####################
### CREATE_TABLES ###
#####################
#####################
#####################


#### PURPOSE: This file contains functions for creating tables from inference output

#### DATE CREATED:  10 SEP 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         24 FEB 2025 (GJD) - Added Create_Caption() and 
####                  Compare_Across_Iterations()
#### EDITS:         09 JUN 2025 (GJD) - Added Get_MC_CI(), Get MC_CI_RMSE_Bootstrap().
####                  - Tweaked Summarize_Multiple_Fits() to save all columns from
####                      Summarize_Single_Fit_Data()
####                




# Get Monte Carlo Confidence Interval for RMSE Using Bootstrap
### rel_column = A vector of estimates (each representing a quantity from a different
###       MC iteration)
Get_MC_CI_RMSE_Bootstrap <- function(bias_col, B=500, n_digits=3) {
  set.seed(3318)
  rmse_boot <- replicate(B, {
    sample_est <- sample(bias_col, replace = TRUE)
    sqrt(mean(sample_est^2))
  })
  lower_bound <- 
    format(round(quantile(rmse_boot, 0.025), n_digits), nsmall=n_digits)
  upper_bound <- 
    format(round(quantile(rmse_boot, 0.975), n_digits), nsmall=n_digits)
  Output <- paste0("(", lower_bound, ", ", upper_bound, ")")
  return(Output)
}


# Get Monte Carlo Confidence Interval
### rel_column = A vector of estimates (each representing a quantity from a different
###       MC iteration)
Get_MC_CI <- function(rel_column, n_digits=3) {
  std_err <- sd(rel_column)/sqrt(length(rel_column))
  estimate <- mean(rel_column)
  lower_bound <- 
    format(round(estimate - qnorm(0.975)*std_err, n_digits), nsmall=n_digits)
  upper_bound <- 
    format(round(estimate + qnorm(0.975)*std_err, n_digits), nsmall=n_digits)
  Output <- paste0("(", lower_bound, ", ", upper_bound, ")")
  return(Output)
}


# Process Output from a Single Fit
### fit_data = Single fit output from Conduct_Inference_for_Driver_Run() output
### filter_statement = Text statement filtering rows of fit_data to relevant rows
Summarize_Single_Fit_Data <- function(fit_data, filter_statement, n_digit_ci=3) {
  relevant_data <- 
    eval(parse(text=paste0("fit_data %>% filter(", filter_statement, ")")))
  relevant_data$reject <- 
    ifelse(relevant_data$CI.Lower>0|relevant_data$CI.Higher<0, 1, 0)
  relevant_data$reject_correct <- 
    ifelse(relevant_data$true_value>0,
           ifelse(relevant_data$CI.Lower>0, 1, 0),
           ifelse(relevant_data$CI.Higher<0, 1, 0))
  #relevant_data$relative_bias <- relevant_data$bias/relevant_data$true_value
  relevant_data$relative_bias <- relevant_data$bias/relevant_data$raw_2
  
  Output <- data.frame(bias=mean(relevant_data$bias),
                       relative_bias=mean(relevant_data$relative_bias),
                       sd=sd(relevant_data$Estimate),
                       RMSE=sqrt(mean(relevant_data$bias^2)),
                       coverage=mean(relevant_data$coverage),
                       power=mean(relevant_data$reject),
                       ci_width=mean(relevant_data$ci_width),
                       bias_CI=Get_MC_CI(relevant_data$bias, 
                                         n_digits=n_digit_ci),
                       relative_bias_CI=Get_MC_CI(relevant_data$relative_bias, 
                                                  n_digits=n_digit_ci),
                       RMSE_CI=Get_MC_CI_RMSE_Bootstrap(bias_col=relevant_data$bias,
                                                        n_digits=n_digit_ci),
                       coverage_CI=Get_MC_CI(relevant_data$coverage, 
                                             n_digits=n_digit_ci))
  return(Output)
}


# Summarize Multiple Fit Outputs
### fit_list = List of fits, indexed by description of fits
###            fit_list[["ex"]] should be a list with elements "fit_data" and
###            [["filter_statement"]] where
###             fit_list[["ex"]][["fit_data"]] = Single fit output from 
###               Conduct_Inference_for_Driver_Run() output
###             fit_list[["ex"]][["filter_statement"]] = Text statement filtering
###               rows of fit_list[["ex"]][["fit_data"]] to relevant rows
Summarize_Multiple_Fits <- function(fit_list, n_digit_ci=3) {
  temp_output <- data.frame()
  for (rowname in names(fit_list)) {
    row_output <- 
      Summarize_Single_Fit_Data(fit_data=fit_list[[rowname]][["fit_data"]],
                                filter_statement=fit_list[[rowname]][["filter_statement"]],
                                n_digit_ci=n_digit_ci)
    new_row <- data.frame(Description=rowname,
                          Bias=row_output$bias,
                          Relative_Bias=row_output$relative_bias,
                          SD=row_output$sd,
                          RMSE=row_output$RMSE,
                          Coverage=row_output$coverage,
                          Power=row_output$power,
                          CI_Width=row_output$ci_width)
    vars_in_new_row <- 
      c("bias", "relative_bias", "sd", "RMSE", "coverage", "power", "ci_width")
    for (var in colnames(row_output)) {if (!(var %in% vars_in_new_row)) {
      new_row[[var]] <- row_output[[var]]
    }}
    temp_output <- rbind(temp_output, new_row)
  }
  Output <- temp_output
  return(Output)
}



# Round and Format Numeric Columns of a Table
### df = Dataframe to round/format
### n_digits = Number of digits to round/format
Round_and_Format_DF <- function(df, n_digits=3) {
  temp_output <- df
  for (col in colnames(df)) {
    if (is.numeric(df[[col]])) {
      temp_output[[col]] <- format(round(temp_output[[col]], n_digits), 
                                   nsmall=n_digits)
    }
  }
  Output <- temp_output
  return(Output)
}


# Create a Caption Object for Easy LaTeX input
### table = Table to create caption for
### caption = Desired table caption
### label = Label for table
Create_Caption <- function(table, caption=NULL, label=NULL) {
  temp_output <- list()
  temp_output$pos <- list()
  temp_output$pos[[1]] <- c(nrow(table))
  out_label <- ""
  if (!is.null(caption)) {
    out_label <- paste0(out_label, "\\caption{", caption, "}")
  }
  if (!is.null(label)) {
    out_label <- paste0(out_label, "\\label{", label, "}\n")
  }
  temp_output$command <- c(out_label)
  
  Output <- temp_output
  return(Output)
}

# Compare a Statistic Across Simulations
## Usually only relevant for comparing different fits
### fit_list_1 = Inference output for the "base" fit
### fit_list_2 = Inference output for the comparison fit
### summary_stat = Which statistic to calculate
Compare_Across_Iterations <- function(fit_list_1, fit_list_2, summary_stat="SE_Ratio") {
  fit_data_1 <- fit_list_1[["fit_data"]]
  fit_data_2 <- fit_list_2[["fit_data"]]
  relevant_data_1 <- 
    eval(parse(text=paste0("fit_data_1 %>% filter(", fit_list_1[["filter_statement"]], ")")))
  relevant_data_2 <- 
    eval(parse(text=paste0("fit_data_2 %>% filter(", fit_list_2[["filter_statement"]], ")")))
  if (summary_stat=="SE_Ratio") {
    se_list <- merge(relevant_data_1[,c("iter", "Std.Err")], relevant_data_2[,c("iter", "Std.Err")], by="iter")
    Output <- mean(se_list$Std.Err.y/se_list$Std.Err.x)
  }
  return(Output)
}
