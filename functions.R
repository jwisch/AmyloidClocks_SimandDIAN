get_conversionTime <- function(df, id_name, value_name, PET_pos_threshold, time_col = "TimefromBaseline"){
  subject_col <- sym(id_name)
  val_col <- sym(value_name)
  time_col <- sym(time_col)
  
  # Step 1: Identify subjects with both Apos = 0 and Apos = 1
  valid_subjects <- df %>%
    group_by(!!subject_col) %>%
    filter(first(Apos) == 0 & any(Apos == 1)) %>%  # First Apos must be 0 and there must be a 1 later
    pull(!!subject_col) %>%
    unique()
  
  # Step 2: Identify the first instance of Apos changing from 0 to 1 in valid subjects
  transition_subjects <- df %>%
    filter(!!subject_col %in% valid_subjects) %>%
    group_by(!!subject_col) %>%
    filter(Apos == 1 & lag(Apos, default = 0) == 0) %>%
    pull(!!subject_col)
  
  # Step 3: For each SUBJECT_LABEL, interpolate TimefromBaseline where TESTVALUE == 0.479
  interpolated_times <- df %>%
    filter(!!subject_col %in% transition_subjects) %>%
    group_by(!!subject_col) %>%
    summarize(
      interpolated_time = approx(
        x = !!val_col, 
        y = !!time_col, 
        xout = PET_pos_threshold, 
        rule = 2)$y  # 'rule = 2' ensures extrapolation if needed
    )
  # Merge the interpolated times back into the original dataframe
  converters <- df %>%
    left_join(interpolated_times, by = id_name)
  
  
}


###############################################################################
# Function to perform linear regression and return slope (rate of change per year)
# calc_rate_of_change function
calc_rate_of_change <- function(sub_data, val_col, time_col) {
  
  # Ensure the value and time columns are numeric
  sub_data[[val_col]] <- as.numeric(sub_data[[val_col]])
  sub_data[[time_col]] <- as.numeric(sub_data[[time_col]])
  
  # Remove rows with NA values in the relevant columns
  sub_data <- sub_data %>%
    filter(!is.na(!!sym(val_col)) & !is.na(!!sym(time_col)))
  
  # Create the formula dynamically using the column names
  formula <- as.formula(paste(val_col, "~", time_col))
  
  # Perform linear regression using the dynamically created formula
  if (nrow(sub_data) > 1) {  # Ensure there is enough data for regression
    model <- lm(formula, data = sub_data)
    return(coef(model)[[time_col]])  # Return the slope coefficient
  } else {
    return(NA)  # Return NA if there isn't enough data
  }
}
compute_zscore <- function(df, col_name) {
  col_sym <- ensym(col_name)
  col_str <- as_string(col_sym)
  new_col <- paste0(col_str, "_Z")
  
  # Convert the column to numeric
  df[[col_str]] <- as.numeric(df[[col_str]])
  
  # Calculate mean and SD for Mutation == 0 and no duplicate IDs
  mean_val <- mean(df[df$Mutation == 0 & !duplicated(df$newid18),][[col_str]], na.rm = TRUE)
  sd_val <- sd(df[df$Mutation == 0 & !duplicated(df$newid18),][[col_str]], na.rm = TRUE)
  
  # Add the new z-score column to df
  df[[new_col]] <- (df[[col_str]] - mean_val) / sd_val
  
  return(df)
}


get_RofC_df <- function(df, id_name, value_name, time_name, flip_sign = FALSE) {
  
  # Define column symbols
  subject_col <- sym(id_name)
  val_col <- value_name
  time_col <- time_name
  
  # Apply the function to each subject and get the rate of change
  RofC <- df %>%
    group_by(!!subject_col) %>%
    group_modify(~ tibble(rate_of_change = calc_rate_of_change(.x, val_col, time_col)))
  
  if (flip_sign) {
    RofC$rate_of_change <- RofC$rate_of_change * -1
  }
  
  # Compute Z-scores and remove outliers
  RofC$Z <- (RofC$rate_of_change - mean(RofC$rate_of_change, na.rm = TRUE)) / 
    sd(RofC$rate_of_change, na.rm = TRUE)
  RofC <- data.frame(RofC)
  RofC <- RofC[RofC$Z < 5 & RofC$Z > -5 & !is.na(RofC$Z), ]
  RofC <- RofC[, !names(RofC) %in% c("Z")]
  # Fit clustering model

  fit <- Mclust(RofC[!is.na(RofC$rate_of_change), ]$rate_of_change, G = 2, 
                model = "V")
  
  # Assign classifications
  RofC$classification <- NA
  RofC[!is.na(RofC$rate_of_change), ]$classification <- fit$classification
  
  # --- Reclassification of the left-hand tail (as in your original code) ---
  row_index <- which.min(abs(RofC$rate_of_change - median(RofC$rate_of_change, na.rm = TRUE)))
  RofC[RofC$rate_of_change < median(RofC$rate_of_change), ]$classification <- 
    RofC[row_index, ]$classification
  
  # --- Ensure class 1 is lower distribution and class 2 is upper distribution ---
  class_means <- tapply(RofC$rate_of_change, RofC$classification, mean, na.rm = TRUE)
  
  if (length(class_means) == 2 && class_means[1] > class_means[2]) {
    # If class 1 is actually the higher group, swap labels
    RofC$classification <- ifelse(RofC$classification == 1, 2, 1)
  }
  
  RofC$classification <- factor(RofC$classification, levels = c(1, 2))
  
  return(RofC)
}


get_relAccumFlag <- function(df, RofC_col_name, rel_accum) {
  # ensure column name is handled correctly
  col_sym <- rlang::sym(RofC_col_name)
  
  # drop NAs and create new column
  df <- df[!is.na(df[[RofC_col_name]]), ]
  df$ReliableAccumulator <- ifelse(df[[RofC_col_name]] > rel_accum, 1, 0)
  
  return(df)
}


# Function to add suffix to column names except 'newid18'
rename_cols <- function(df, suffix) {
  colnames(df)[colnames(df) != "newid18"] <- paste0(colnames(df)[colnames(df) != "newid18"], "_", suffix)
  return(df)
}



get_ReliableAccumThreshold <- function(DF, Zscore = 1.96){
  # Identify all rate_of_change columns
  rate_cols <- grep("^rate_of_change_", names(DF), value = TRUE)
  
  # Create an empty list to store thresholds
  thresholds <- list()
  
  # Loop through each rate_of_change column
  for (rate_col in rate_cols) {
    # Identify the corresponding classification column
    class_col <- gsub("rate_of_change", "classification", rate_col)
    
    # Compute the threshold if classification column exists
    if (class_col %in% names(DF)) {
      thresholds[[rate_col]] <- mean(DF[DF[[class_col]] == 1, rate_col], na.rm = TRUE) +
        Zscore * sd(DF[DF[[class_col]] == 1, rate_col], na.rm = TRUE)
    }
  }
  
  # Convert to a named vector for easy access
  thresholds <- unlist(thresholds)
  
  return(thresholds)}


# Function to force the numbering scheme for the GMM
# Apply the function to all columns in the dataframe that contain only 1's, 2's, and NA's
# Function to force the numbering scheme and keep the factor levels
force_order_factor <- function(col) {
  # Convert to character first if it's a factor to facilitate manipulation
  col <- as.character(col)
  
  # Calculate the frequency of 1's and 2's (ignoring NAs)
  freq_1 <- sum(col == "1", na.rm = TRUE)
  freq_2 <- sum(col == "2", na.rm = TRUE)
  
  # If 2 is more frequent than 1, swap them
  if (freq_2 > freq_1) {
    col <- ifelse(col == "1", "2", ifelse(col == "2", "1", col))
  }
  
  # Convert back to factor, ensuring levels are set correctly
  col <- factor(col, levels = c("1", "2"))
  
  return(col)
}

#Plotting rate of change vs. baseline level
get_BollackPlot <- function(dataframes_list_index, xval){
  xval <- ensym(xval)
  p <- ggplot(dataframes_list[[dataframes_list_index]], aes(x = !!xval, y = rate_of_change, 
                                                            colour = factor(ReliableAccumulator), fill = factor(ReliableAccumulator))) + 
    geom_point() + theme_bw() + theme(legend.position = "bottom") +
    geom_smooth(aes(group = 1), method = "gam", colour = "black")
  return(p)
  
}



get_BaselineValforAccum <- function(DF, DF_BASELINE, thresholds, param){
  param_sym <- sym(param)
  tmp <- data.frame("newid18" = DF[, "newid18"],
                    "rate_of_change" = DF[, c(paste0("rate_of_change_", param))])
  tmp$ReliableAccumulator <- ifelse(as.numeric(tmp$rate_of_change) > thresholds[[paste0("rate_of_change_", param)]], 1, 0)
  df_baseline <- DF_BASELINE %>%
    filter(!is.na(!!param_sym)) %>%#,  # Remove rows where Alamar_AB42_AB40 is NA
    #abs(!!param_sym - mean(!!param_sym, na.rm = TRUE)) <= 3 * sd(!!param_sym, na.rm = TRUE)) %>% #remove outliers
    group_by(newid18) %>%
    slice(1) %>%  # Select the first row per group
    ungroup()
  tmp <- merge(tmp, df_baseline, by = "newid18", all = FALSE)
  tmp <- tmp[, c("newid18", "rate_of_change", "ReliableAccumulator", param)]
  
  tmp <- tmp %>%
    filter(abs(rate_of_change - mean(rate_of_change, na.rm = TRUE)) <= 3 * sd(rate_of_change, na.rm = TRUE))
  
  BaselineThresholdforReliableAccum <- cutpointr(tmp, !!param_sym, ReliableAccumulator, 
                                                 method = maximize_metric, metric = youden, na.rm = TRUE) #2156
  return(list(tmp, BaselineThresholdforReliableAccum))
  
  
}

#Identifying th eminimum baseline value to predict future accumulation
get_optimal_threshold <- function(df, df_RofC, CL_col, classification_col) {
  # Convert column names to symbols for tidy evaluation
  CL_sym <- sym(CL_col)
  class_sym <- sym(classification_col)
  
  # Remove duplicates based on 'newid18' and keep only relevant columns
  tmp <- df %>%
    distinct(newid18, .keep_all = TRUE) %>%
    select(newid18, !!CL_sym)
  
  # Merge with classification data from df_RofC
  tmp <- tmp %>%
    inner_join(df_RofC %>% select(newid18, !!class_sym), by = "newid18") %>%
    filter(!is.na(!!class_sym))
  
  # Compute the optimal threshold using cutpointr
  optimal_threshold <- cutpointr(tmp, x = !!CL_sym, class = !!class_sym, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
  result <- data.frame("Biomarker" = CL_col,
                       "optimal_cutpoint" = optimal_threshold$optimal_cutpoint,
                       "sensitivity" = optimal_threshold$sensitivity,
                       "specificity" = optimal_threshold$specificity)
  
  return(result)
}

get_groupedScatter <- function(df, yVal, ymin, ymax){
  yVal_sym <- ensym(yVal)
  
  # Calculate mean and SD for Mutation == 0
  stats <- df %>%
    filter(Mutation == 0) %>%
    summarise(
      mean_y = mean(!!yVal_sym, na.rm = TRUE),
      sd_y = sd(!!yVal_sym, na.rm = TRUE)
    )
  
  mean_val <- stats$mean_y
  sd_val <- stats$sd_y
  
  # Add z-scored column based on Mutation == 0
  df <- df %>%
    mutate(z_score = (!!yVal_sym - mean_val) / sd_val)
  
  # Plot
  p <- ggplot(df, aes(x = DIAN_EYO, y = z_score, colour = factor(Mutation), fill = factor(Mutation))) +
    geom_point() +
    geom_smooth(method = "gam") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab(paste0(as_label(yVal_sym), " Z")) + ylim(c(ymin, ymax))
  
  return(p)
}


get_RofC_scatter_and_ROC <- function(DF, predictor, outcome, predictor_label, outcome_label){
  predictor <- ensym(predictor)  # Capture column name as symbol
  outcome <- ensym(outcome)      # Capture column name as symbol
  cutpointr_result <- cutpointr(
    DF, 
    x = !!predictor, 
    class = !!outcome, 
    method = maximize_metric, 
    metric = youden, 
    na.rm = TRUE, 
    boot_runs = 1000
  )
  
  
  p2 <- ggplot(data.frame(cutpointr_result$roc_curve), aes(x = fpr, y = tpr)) +
    geom_line() +
    geom_line(data = data.frame(cutpointr_result$roc_curve), aes(x = fpr, y = tpr), colour = "#007B82") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Add the diagonal line (random classifier)
    labs(x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    annotate("text", x = 0.5, y = 0.3, hjust = 0, label = paste0("Sensitivity = ", round(cutpointr_result$sensitivity, 3))) +
    annotate("text", x = 0.5, y = 0.2, hjust = 0, label = paste0("Specificity = ", round(cutpointr_result$specificity, 3))) +
    geom_hline(yintercept = cutpointr_result$sensitivity, linetype = "dashed", colour = "#fc8d59") +
    geom_vline(xintercept = (1 - cutpointr_result$specificity), linetype = "dashed", colour = "#fc8d59")
  
  cp_df <- data.frame(cutpointr_result$roc_curve)
  cp_df$Group <- predictor_label
  return(list(p2, cp_df, cutpointr_result))
}



generate_EYOcomparison_plots <- function(DF, IDname, TimefromBaselineName, ValueofInterest,
                                         cutpointObj, TimefromAposName) {
  
  # Capture column names as symbols
  ID_sym <- ensym(IDname)
  TimeBase_sym <- ensym(TimefromBaselineName)
  Val_sym <- ensym(ValueofInterest)
  TimeApos_sym <- ensym(TimefromAposName)
  
  # Create working copy
  DF <- DF %>%
    group_by(!!ID_sym) %>%
    mutate(TimefromBaseline_PET = as.numeric(!!TimeApos_sym - first(!!TimeApos_sym))) %>%
    ungroup()
  
  # Subset for MAE calculation
  DF_mae <- DF %>%
    filter(!is.na(TimefromBaseline_PET),
           !is.na(!!TimeBase_sym),
           !!TimeBase_sym != 0,
           !!Val_sym > cutpointObj$optimal_cutpoint)
  
  MAE_val <- MLmetrics::MAE(DF_mae$TimefromBaseline_PET, DF_mae[[as_string(TimeBase_sym)]])
  # Model: TimefromApos ~ DIAN_EYO
  mod1 <- lmer(formula = as.formula(paste0(as_string(TimeApos_sym), " ~ DIAN_EYO + (1 | ", as_string(ID_sym), ")")),
               data = DF %>% filter(!!Val_sym > cutpointObj$optimal_cutpoint, DIAN_EYO > -25))
  r2_mod1 <- r.squaredGLMM(mod1)
  fix1 <- fixef(mod1)
  intercept1 <- round(fix1[1], 3)
  slope1 <- round(fix1["DIAN_EYO"], 3)
  
  # Model: ValueofInterest ~ DIAN_EYO
  mod2 <- lmer(formula = as.formula(paste0(as_string(Val_sym), " ~ DIAN_EYO + (1 | ", as_string(ID_sym), ")")),
               data = DF %>% filter(!!Val_sym > cutpointObj$optimal_cutpoint, DIAN_EYO > -25))
  r2_mod2 <- r.squaredGLMM(mod2)
  fix2 <- fixef(mod2)
  intercept2 <- round(fix2[1], 3)
  slope2 <- round(fix2["DIAN_EYO"], 3)
  
  #Annotate coordinates
  x_pos_p1 <- min(DF$DIAN_EYO, na.rm = TRUE) + 8
  y_pos_p1 <- 1.15 * min(DF[[as_string(Val_sym)]], na.rm = TRUE) 
  
  x_pos_p2 <- min(DF$DIAN_EYO, na.rm = TRUE) + 8
  y_pos_p2 <- 0.85 * max(DF[[as_string(TimeApos_sym)]], na.rm = TRUE)
  
  x_pos_p3 <- min(DF[[as_string(TimeBase_sym)]], na.rm = TRUE) + 1
  y_pos_p3 <- 0.85 * max(DF$TimefromBaseline_PET, na.rm = TRUE)
  
  # Plot 1: Value of interest vs DIAN_EYO
  p1 <- ggplot(DF %>% filter(!!Val_sym > cutpointObj$optimal_cutpoint),
               aes(x = DIAN_EYO, y = !!Val_sym, group = !!ID_sym)) +
    geom_point(alpha = 0.8) + geom_line(alpha = 0.3) +
    geom_smooth(method = "lm", aes(group = 1)) + theme_bw() +
    annotate("text", x = -22, y = y_pos_p1,
             label = paste0("y = ", slope2, "x + ", intercept2,
                            "\nR2 = ", round(r2_mod2[1], 3)), hjust = 0) +
    xlim(-25, 20)
  
  # Plot 2: TimefromApos vs DIAN_EYO
  p2 <- ggplot(DF %>% filter(!!Val_sym > cutpointObj$optimal_cutpoint),
               aes(x = DIAN_EYO, y = !!TimeApos_sym, group = !!ID_sym)) +
    geom_point(alpha = 0.8) + geom_line(alpha = 0.3) +
    geom_abline(intercept = 10, slope = 1, linetype = "dashed", colour = "red") +
    geom_smooth(method = "lm", aes(group = 1)) + theme_bw() +
    annotate("text", x = -22, y = -2.5, hjust = 0,
             label = paste0("y = ", slope1, "x + ", intercept1,
                            "\nR2 = ", round(r2_mod1[1], 3))) +
    xlim(-25, 20) + ylim(-5, 20)
  
  # Plot 3: Estimated vs actual time from baseline
  p3 <- ggplot(DF %>% filter(!!Val_sym > cutpointObj$optimal_cutpoint),
               aes(x = !!TimeBase_sym, y = TimefromBaseline_PET, group = !!ID_sym)) +
    geom_point(alpha = 0.8) + geom_line(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
    theme_bw() +
    annotate("text", x = x_pos_p3, y = y_pos_p3, hjust = 0,
             label = paste0("MAE = ", round(MAE_val, 2), " years")) +
    xlab("Actual Time from Baseline") + ylab("Estimated Time from Baseline (PET)")
  
  # Display side-by-side
  grid.arrange(p1, p2, p3, ncol = 3)
  
  # Return models and MAE
  return(list(model_value = mod2,
              model_time = mod1,
              MAE = MAE_val,
              r2_value = r2_mod2,
              r2_time = r2_mod1))
}


Zscore_to_Time <- function(xout_vec, result_df) {
  df <- data.frame(
    Zscore = xout_vec,
    
    Estimate = approx(y = result_df$Time_to_Positivity, 
                      x = result_df$Estimate, 
                      xout = xout_vec, rule = 2)$y,
    
    CI_Lower = approx(y = result_df$Time_to_Positivity, 
                      x = result_df$CI_Upper, 
                      xout = xout_vec, rule = 2)$y,
    
    CI_Upper = approx(y = result_df$Time_to_Positivity, 
                      x = result_df$CI_Lower, 
                      xout = xout_vec, rule = 2)$y
  )
  df$TimeRange <- df$CI_Upper - df$CI_Lower
  df$group <- rep(unique(result_df$group), length(xout_vec))
  
  return(df)
}



Estimate_to_CI <- function(estimate_vec, result_df) {
  df <- data.frame(
    Estimate = estimate_vec,
    
    CI_Lower = approx(
      x = result_df$Estimate,
      y = result_df$CI_Lower,
      xout = estimate_vec,
      rule = 2
    )$y,
    
    CI_Upper = approx(
      x = result_df$Estimate,
      y = result_df$CI_Upper,
      xout = estimate_vec,
      rule = 2
    )$y
  )
  
  df$CI_Range <- df$CI_Upper - df$CI_Lower
  df$group <- rep(unique(result_df$group), length(estimate_vec))
  return(df)
}


plot_biomarker_groups <- function(df, biomarker_col, time_col, result_obj, cp_obj) {
  # Ensure we handle strings as column names
  biomarker <- rlang::sym(biomarker_col)
  time_var  <- rlang::sym(time_col)
  
  # Subset dataframe based on cutpoint
  plot_df <- df[df[[biomarker_col]] > cp_obj$optimal_cutpoint, 
                c("newid18", "VISITAGEc", biomarker_col, time_col)]
  
  # Add group assignment
  plot_df$group <- ifelse(plot_df$newid18 %in% result_obj$Anegalways$newid18, "NegativeAlways",
                          ifelse(plot_df$newid18 %in% result_obj$Aposalways$newid18, "PositiveAlways",
                                 ifelse(plot_df$newid18 %in% result_obj$Converters$newid18, "Converters", "ruhRoh")))
  
  # Full plot
  p1 <- ggplot(plot_df[!plot_df$group == "ruhRoh",], 
               aes(x = !!time_var, y = !!biomarker, group = newid18, colour = group)) +
    geom_point(alpha = 0.5) + 
    geom_line(alpha = 0.6) + 
    theme_bw() + 
    # ylim(c(1.1, 5.5)) + 
    theme(legend.position = "bottom")
  
  # Converters-only plot
  p2 <- ggplot(plot_df[plot_df$group == "NegativeAlways",], 
               aes(x = !!time_var, y = !!biomarker, group = newid18, colour = group)) +
    geom_point(alpha = 0.5) + 
    geom_line(alpha = 0.6) + 
    theme_bw() + 
    # ylim(c(1.1, 5.5)) + 
    theme(legend.position = "none") +
    ylab("Negative Always")
  p3 <- ggplot(plot_df[plot_df$group == "Converters",], 
               aes(x = !!time_var, y = !!biomarker, group = newid18, colour = group)) +
    geom_point(alpha = 0.5) + 
    geom_line(alpha = 0.6) + 
    theme_bw() + 
    # ylim(c(1.1, 5.5)) + 
    theme(legend.position = "none") +
    ylab("Converters")
  p4 <- ggplot(plot_df[plot_df$group == "PositiveAlways",], 
               aes(x = !!time_var, y = !!biomarker, group = newid18, colour = group)) +
    geom_point(alpha = 0.5) + 
    geom_line(alpha = 0.6) + 
    theme_bw() + 
    # ylim(c(1.1, 5.5)) + 
    theme(legend.position = "none") +
    ylab("Positive Always")
  
  # Arrange plots side by side
  layout_matrix <- rbind(c(1, 1, 2), c(1, 1, 3), c(1, 1, 4))
  grid.arrange(p1, p2, p3, p4, layout_matrix = layout_matrix)
}


split_status_groups <- function(df, subject_col, status_col, age_col, atime_col, val_col) {
  subject_col <- enquo(subject_col)
  status_col  <- enquo(status_col)
  age_col     <- enquo(age_col)
  atime_col   <- enquo(atime_col)
  val_col     <- enquo(val_col)
  
  # Always negative: all 0, keep row with max age
  df_Anegalways <- df %>%
    group_by(!!subject_col) %>%
    filter(all(!!status_col == 0)) %>%
    filter(!!age_col == max(!!age_col)) %>%
    ungroup()
  
  # Always positive: all 1, keep row with min age
  df_Aposalways <- df %>%
    group_by(!!subject_col) %>%
    filter(all(!!status_col == 1)) %>%
    filter(!!age_col == min(!!age_col)) %>%
    ungroup()
  
  # Converters: 0 → 1
  df_converters <- df %>%
    group_by(!!subject_col) %>%
    filter(any(!!status_col == 0) & any(!!status_col == 1)) %>%
    summarise(
      last_neg_age    = max(if_else(!!status_col == 0, !!age_col, -Inf), na.rm = TRUE),
      first_pos_age   = min(if_else(!!status_col == 1, !!age_col, Inf), na.rm = TRUE),
      last_neg_atime  = max(if_else(!!status_col == 0, !!atime_col, -Inf), na.rm = TRUE),
      first_pos_atime = min(if_else(!!status_col == 1, !!atime_col, Inf), na.rm = TRUE),
      last_neg_val    = max(if_else(!!status_col == 0, !!val_col, -Inf), na.rm = TRUE),
      first_pos_val   = min(if_else(!!status_col == 1, !!val_col, Inf), na.rm = TRUE),
      midpoint_age    = mean(c(last_neg_age, first_pos_age)),
      midpoint_atime  = mean(c(last_neg_atime, first_pos_atime)),
      midpoint_val    = mean(c(last_neg_val, first_pos_val)),
      .groups = "drop"
    ) %>%
    transmute(
      !!quo_name(subject_col) := !!subject_col,
      !!quo_name(age_col)     := midpoint_age,
      !!quo_name(atime_col)   := midpoint_atime,
      !!quo_name(val_col)     := midpoint_val,
      !!quo_name(status_col)  := 1  # converters end up positive
    )
  
  # Return as a list
  list(
    Anegalways = df_Anegalways,
    Aposalways = df_Aposalways,
    Converters = df_converters
  )
}



# estimate_simulation_params <- function(
#     df,
#     id_col,
#     age_col,
#     biomarker_col,
#     minMeaningfulZ,
#     AposThresh
# ) {
#   # Non-standard evaluation to allow column names as strings
#   id_col <- rlang::sym(id_col)
#   age_col <- rlang::sym(age_col)
#   biomarker_col <- rlang::sym(biomarker_col)
#   
#   # Filter for IDs with more than 1 visit (i.e., longitudinal)
#   df_long <- df %>%
#     group_by(!!id_col) %>%
#     filter(n() > 1) %>%
#     arrange(!!id_col, !!age_col) %>%
#     ungroup()
#   
#   # 1. n_ids
#   n_ids <- n_distinct(df_long[[rlang::as_string(id_col)]])
#   
#   # 2. Average interval and interval noise
#   interval_summary <- df_long %>%
#     group_by(!!id_col) %>%
#     summarize(
#       diffs = diff(sort(!!age_col)),
#       .groups = "drop"
#     ) %>%
#     unnest(cols = c(diffs))
#   
#   avg_interval <- mean(interval_summary$diffs, na.rm = TRUE)
#   interval_noise <- sd(interval_summary$diffs, na.rm = TRUE)
#   
#   # 3. Estimate per-subject slopes (annualized change)
#   indiv_slopes <- df_long %>%
#     group_by(!!id_col) %>%
#     summarize(
#       slope = ifelse(n() > 1,
#                      coef(lm(!!biomarker_col ~ !!age_col))[2],
#                      NA_real_),
#       .groups = "drop"
#     ) %>%
#     drop_na()
# 
#   # 4. Focus on individuals whose biomarker starts above minMeaningfulZ
#   #    (so we’re estimating “annual_rate” in meaningful accumulation range)
#   ids_above_thresh <- df_long %>%
#     group_by(!!id_col) %>%
#     summarize(first_val = first(!!biomarker_col), .groups = "drop") %>%
#     filter(first_val >= minMeaningfulZ) %>%
#     pull(!!id_col)
#   
#   slopes_above <- indiv_slopes %>%
#     filter(!!id_col %in% ids_above_thresh)
#   
#   annual_rate <- mean(slopes_above$slope, na.rm = TRUE)
#   var_annual_rate <- sd(slopes_above$slope, na.rm = TRUE)
#   
#   # Return as named list
#   list(
#     n_ids = n_ids,
#     avg_interval = avg_interval,
#     interval_noise = interval_noise,
#     annual_rate = annual_rate,
#     var_annual_rate = var_annual_rate
#   )
# }

estimate_simulation_params <- function(
    df,
    id_col,
    age_col,
    biomarker_col,
    minMeaningfulZ,
    AposThresh
) {
  # Non-standard evaluation
  id_col <- rlang::sym(id_col)
  age_col <- rlang::sym(age_col)
  biomarker_col <- rlang::sym(biomarker_col)
  
  # --- 1. Filter for IDs with >1 visit ---
  df_long <- df %>%
    group_by(!!id_col) %>%
    filter(n() > 1) %>%
    arrange(!!id_col, !!age_col) %>%
    ungroup()
  
  # --- 2. n_ids ---
  n_ids <- n_distinct(df_long[[rlang::as_string(id_col)]])
  
  # --- 3. Interval summary ---
  interval_summary <- df_long %>%
    group_by(!!id_col) %>%
    summarize(
      diffs = diff(sort(!!age_col)),
      .groups = "drop"
    ) %>%
    unnest(cols = c(diffs))
  avg_interval <- mean(interval_summary$diffs, na.rm = TRUE)
  interval_noise <- sd(interval_summary$diffs, na.rm = TRUE)
  
  # --- 4. Compute per-interval slopes safely ---
  slopes_df <- df_long %>%
    group_by(!!id_col) %>%
    arrange(!!age_col, .by_group = TRUE) %>%
    mutate(
      delta_bio = lead(!!biomarker_col) - !!biomarker_col,
      delta_age = lead(!!age_col) - !!age_col,
      slope = delta_bio / delta_age,
      Z_mid = (!!biomarker_col + lead(!!biomarker_col)) / 2
    ) %>%
    ungroup() %>%
    # remove non-finite slopes (Inf, NaN, etc.)
    filter(is.finite(slope), is.finite(Z_mid))
  
  # --- 5. Compute region-specific slope summaries ---
  neg_prop_summary <- slopes_df %>%
    mutate(region = case_when(
      Z_mid < minMeaningfulZ ~ "below_minMeaningfulZ",
      Z_mid < AposThresh ~ "between_min_and_Apos",
      TRUE ~ "above_AposThresh"
    )) %>%
    group_by(region) %>%
    summarize(
      n = n(),
      n_neg = sum(slope < 0, na.rm = TRUE),
      prop_neg = n_neg / n,
      mean_slope = ifelse(any(is.finite(slope)), mean(slope, na.rm = TRUE), NA_real_),
      sd_slope = ifelse(any(is.finite(slope)), sd(slope, na.rm = TRUE), NA_real_),
      .groups = "drop"
    )
  
  # --- 6. Estimate per-subject annual slopes ---
  indiv_slopes <- df_long %>%
    group_by(!!id_col) %>%
    summarize(
      slope = {
        df_sub <- tibble(age = !!age_col, bio = !!biomarker_col) %>%
          filter(!is.na(age), !is.na(bio))
        if (nrow(df_sub) > 1) coef(lm(bio ~ age, data = df_sub))[2] else NA_real_
      },
      start_val = first(!!biomarker_col),
      .groups = "drop"
    ) %>%
    filter(is.finite(slope))
  
  ids_above_thresh <- df_long %>%
    group_by(!!id_col) %>%
    summarize(first_val = first(!!biomarker_col), .groups = "drop") %>%
    filter(first_val >= minMeaningfulZ) %>%
    pull(!!id_col)
  
  slopes_above <- indiv_slopes %>%
    filter(!!id_col %in% ids_above_thresh)
  
  annual_rate <- mean(slopes_above$slope, na.rm = TRUE)
  var_annual_rate <- sd(slopes_above$slope, na.rm = TRUE)
  
  # --- 7. Return results ---
  list(
    n_ids = n_ids,
    avg_interval = avg_interval,
    interval_noise = interval_noise,
    annual_rate = annual_rate,
    var_annual_rate = var_annual_rate,
    neg_slope_summary = neg_prop_summary
  )
}

adjust_crossing_times <- function(df, mean_within_below = 0.44, sd_within_below = 0.042){
  df_adjusted <- df %>%
    group_by(ID) %>%
    group_modify(~{
      sub <- .x
      if(!any(sub$ever_cross)) return(sub)
      
      first_cross <- unique(sub$first_cross)
      
      if(is.na(first_cross) || first_cross >= 3) return(sub)
  
  new_cross <- runif(1, 3, 10)
  shift <- new_cross - first_cross
  
  sub <- sub %>%
    mutate(TimefromBaseline = TimefromBaseline + shift)
  
  pre_cross <- sub$TimefromBaseline < new_cross
  
  n_pre <- sum(pre_cross)
  if (n_pre > 0) {
    sub$Z[pre_cross] <- abs(rnorm(n_pre, mean_within_below, sd_within_below))
  }
  
  sub
    }) %>% ungroup()
return(df_adjusted)}


