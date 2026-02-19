library(cutpointr)
library(dplyr)
library(ggplot2)
library(rlang)

run_cutpointr_analysis <- function(df, biomarker, outcome, y1, y2, 
                                   predictor_for_x = NULL, XLAB = NULL, YLAB = NULL) {
  # --- Step 1: Run bootstrapped cutpointr ---
  cp_boot <- cutpointr(
    df,
    x = !!sym(biomarker),
    class = !!sym(outcome),
    method = maximize_metric,
    metric = youden,
    na.rm = TRUE,
    boot_runs = 1000,
    boot_stratify = TRUE
  )
  
  # --- Step 2: Extract bootstrapped summary safely ---
  boot_summary <- summary(cp_boot)$boot[[1]]
  
  # --- Step 3: Remove infinite values before averaging ---
  if ("optimal_cutpoint" %in% names(boot_summary)) {
    boot_summary <- boot_summary %>%
      dplyr::filter(is.finite(optimal_cutpoint))
  }
  
  # --- Step 4: Compute mean and CI excluding Inf values ---
  threshold_mean <- mean(boot_summary$optimal_cutpoint, na.rm = TRUE)
  threshold_CI <- quantile(boot_summary$optimal_cutpoint, probs = c(0.025, 0.975), na.rm = TRUE)
  

  # --- Step 3: Run non-bootstrapped cutpointr for counts ---
  cp_noboot <- cutpointr(
    df,
    x = !!sym(biomarker),
    class = !!sym(outcome),
    method = maximize_metric,
    metric = youden,
    na.rm = TRUE
  )
  
  summary_df <- as.data.frame(cp_noboot$roc_curve)

  
  annot <- summary_df[summary_df$x.sorted == cp_noboot$optimal_cutpoint, ] %>%
    mutate(
      precision = tp / (tp + fp),
      recall = tp / (tp + fn),
      TP = tp,
      TN = tn,
      FP = fp,
      FN = fn
    ) %>%
    select(precision, recall, TP, TN, FP, FN)
  auc_val <- round(cp_noboot$AUC, 3) 
  # --- Step 4: Build annotation text ---
  annotation_text <- sprintf(
    "Threshold = %.2f [%s, %s]\nAUC = %.3f\nPrecision = %.3f\nRecall = %.3f",
    boot_summary[boot_summary$Variable == "optimal_cutpoint", "Median"],
    ifelse(is.na(boot_summary[boot_summary$Variable == "optimal_cutpoint", "5%"]), "-", sprintf("%.3f", boot_summary[boot_summary$Variable == "optimal_cutpoint", "5%"])),
    ifelse(is.na(boot_summary[boot_summary$Variable == "optimal_cutpoint", "5%"]), "-", sprintf("%.3f", boot_summary[boot_summary$Variable == "optimal_cutpoint", "95%"])),
    auc_val,
    annot$precision,
    annot$recall
  )
  

  
  # --- Step 5: ROC curve plot ---
  roc_df <- data.frame(cp_noboot$roc_curve)

  
  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "#007B82", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    coord_equal() +
    labs(
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    annotate(
      "label",
      x = 0.25, y = 0.2,
      label = annotation_text,
      hjust = 0, size = 4.2,
      color = "black",
      fill = "white", alpha = 0.8
    ) +
    theme_minimal(base_size = 14)
  
  # --- Step 6: Compute proportions ---
  annot <- annot %>%
    mutate(
      total = TP + TN + FP + FN,
      TP_perc = round(TP / total * 100, 1),
      FP_perc = round(FP / total * 100, 1),
      TN_perc = round(TN / total * 100, 1),
      FN_perc = round(FN / total * 100, 1)
    )
  
  # --- Step 7: Classification plot ---
  Apos_thresh <- boot_summary[boot_summary$Variable == "optimal_cutpoint", "Mean"]
  
  # Allow user to specify which column should be the x-axis (default = biomarker)
  if (is.null(predictor_for_x)) predictor_for_x <- biomarker
  
  class_plot <- ggplot(df, aes(x = !!sym(predictor_for_x), y = !!sym(biomarker))) +
    geom_point() +
    geom_vline(xintercept = Apos_thresh, linetype = "dashed") +
    theme_bw() +
    geom_hline(yintercept = Apos_thresh, linetype = "dashed") +
    annotate("label", x = 2.5, y = y2, hjust = 1,
             label = paste0("False Positive\n", annot$FP, " (", annot$FP_perc, "%)")) +
    annotate("label", x = 2.5, y = y1, hjust = 1,
             label = paste0("True Negative\n", annot$TN, " (", annot$TN_perc, "%)")) +
    annotate("label", x = 15, y = y2, hjust = 1,
             label = paste0("True Positive\n", annot$TP, " (", annot$TP_perc, "%)")) +
    annotate("label", x = 15, y = y1, hjust = 1,
             label = paste0("False Negative\n", annot$FN, " (", annot$FN_perc, "%)")) +
    xlab(XLAB) + ylab(YLAB)
  
  # --- Step 8: Return results ---
  return(list(
    summary = annot,
    threshold_mean = Apos_thresh,
    roc_plot = roc_plot,
    classification_plot = class_plot
  ))
}


plot_lmer_by_class <- function(df, yvar, class_col, xlab = "Time from baseline (years)", ylab = "Biomarker (Z)") {
  
  y_sym <- sym(yvar)
  class_sym <- sym(class_col)
  
  plot_full <- ggplot(df,
         aes(x = TimefromBaseline, y = !!y_sym, group = newid18)) +
    geom_point() +
    geom_line(alpha = 0.6) +
    theme_bw() +
    labs(x = xlab, y = ylab) +
    theme(plot.title = element_text(hjust = 0.5))
  
  make_plot <- function(class_value) {
    # Fit mixed model
    mod <- lmer(as.formula(paste(yvar, "~ TimefromBaseline + (1 | newid18)")),
                data = df %>% filter(!!class_sym == class_value))
    
    # Predictions over time range
    pred_df <- data.frame(TimefromBaseline = seq(
      min(df$TimefromBaseline, na.rm = TRUE),
      max(df$TimefromBaseline, na.rm = TRUE),
      length.out = 100
    ))
    
    emm <- emmeans(mod, ~ TimefromBaseline, at = list(TimefromBaseline = pred_df$TimefromBaseline))
    emm_sum <- summary(emm)
    
    pred_df$fit <- emm_sum$emmean
    pred_df$lower <- emm_sum$emmean - 1.96 * emm_sum$SE
    pred_df$upper <- emm_sum$emmean + 1.96 * emm_sum$SE
    
    # Fixed effect estimate
    time_estimate <- fixef(mod)["TimefromBaseline"]
    
    # Plot
    ggplot(df %>% filter(!!class_sym == class_value),
           aes(x = TimefromBaseline, y = !!y_sym, group = newid18)) +
      geom_point() +
      geom_line(alpha = 0.6) +
      geom_line(data = pred_df, aes(y = fit, x = TimefromBaseline),
                inherit.aes = FALSE, color = "#FFD700", size = 2, linetype = "dashed") +
      annotate("text",
               x = max(df$TimefromBaseline, na.rm = TRUE),
               y = max(df[[yvar]], na.rm = TRUE),
               label = paste0(round(time_estimate, 3), " Z/year"),
               hjust = 1, vjust = 1, color = "#007B82", size = 5) +
      theme_bw() +
      labs(x = xlab, y = ylab) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  plots <- list(
    low = make_plot(0),
    high = make_plot(1),
    full = plot_full
  )
  
  return(plots)
}


define_baseline_rel_to_RelAccum <- function(df, col_name, cp_obj, thresh_col_name = "ReliableAccumThresh") {
  
  col_sym <- sym(col_name)
  thresh_sym <- sym(thresh_col_name)
  # Get the earliest visit value per newid18
  earliest_vals <- df %>%
    filter(!is.na(!!col_sym)) %>%
    group_by(newid18) %>%
    slice_min(VISITAGEc, with_ties = FALSE) %>%
    ungroup() %>%
    select(newid18, val_at_min_age = !!col_sym)
  
  # Merge back and define ReliableAccumThresh
  df <- df %>%
    left_join(earliest_vals, by = "newid18") %>%
    mutate(
      !!thresh_sym := case_when(
        is.na(val_at_min_age) ~ NA_real_,
        val_at_min_age < cp_obj$optimal_cutpoint ~ 0,
        val_at_min_age >= cp_obj$optimal_cutpoint ~ 1
      )
    ) %>%
    select(-val_at_min_age)
  
  return(df)
}


combine_RofC_and_df <- function(df, biomarker_col, RofC_df, RofC_col, rel_accum) {
  
  # Ensure biomarker_col is a symbol for tidy evaluation
  biomarker_sym <- sym(biomarker_col)
  RofC_sym <- sym(RofC_col)
  
  # Get earliest visit info per subject
  df_earliest <- df %>%
    filter(!is.na(!!biomarker_sym)) %>%
    group_by(newid18) %>%
    slice_min(VISITAGEc, with_ties = FALSE) %>%  # keep only the row with minimum VISITAGEc
    ungroup() %>%
    select(newid18, VISITAGEc, SEX, RACE, apoe, Mutation, fam_mutation, CL_Z, !!biomarker_sym)
  # Merge with RofC dataframe using dplyr::select for dynamic column
  merged_df <- merge(
    RofC_df %>%
      select(newid18, !!RofC_sym, classification),
    df_earliest,
    by = "newid18",
    all = FALSE
  )
  
  
  merged_df$ReliableAccumulator <- ifelse(merged_df[[RofC_sym]] > rel_accum, 1, 0)
  return(merged_df)
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

plot_rate_of_change <- function(df, xvar, cutoff, bins = 30, adjust = 2,
                                xi, omega, alpha) {
  
  # Compute max histogram count to place annotation dynamically
  ymax <- ggplot_build(
    ggplot(df, aes(x = !!sym(xvar))) +
      geom_histogram(bins = bins)
  )$data[[1]]$count %>% max(na.rm = TRUE)
  
  ggplot(df, aes(x = !!sym(xvar), fill = classification)) +
    # Histogram
    geom_histogram(alpha = 0.3, position = "identity", bins = bins) +
    # Vertical line
    geom_vline(xintercept = cutoff, linetype = "dashed", colour = "#FFD700", size = 2) +
    
    # Dynamic annotation near top of plot
    annotate(
      geom = "label",
      x = 1.1 * cutoff,
      y = 0.8 * ymax,  # 90% of histogram height
      hjust = 0,
      vjust = 1,
      label = paste0(round(cutoff, 2), " Z/year\nxi = ",
                     xi, "\nomega = ", omega, "\nalpha = ", alpha),
      fill = "#FFD700",      # light yellow background
      color = "black",       # text color
      label.size = 0.3,      # border thickness
      label.r = unit(0.15, "lines")  # slight rounding of corners
    ) +
    
    xlim(c(-2.5, 2.5)) +
    xlab("Rate of Change (Z/year)") +
    # Smoothed density scaled to histogram height
    geom_density(aes(y = after_stat(density * n * diff(range(x)) / bins)),
                 alpha = 0, adjust = adjust) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) + 
    scale_fill_manual(values = c("#007B82", "grey67"))
}


get_RofC_vs_baseline_scatter <- function(df, RofC_col, biomarker_col, thresh_obj, 
                                         colour_col = "ReliableAccumulator",
                                         XLAB = NULL, YLAB = NULL,
                                         rel_accum_thresh) {
  
  # Convert to symbols for tidy evaluation
  RofC_sym <- sym(RofC_col)
  biomarker_sym <- sym(biomarker_col)
  colour_sym <- sym(colour_col)
  
  # Determine the xintercept value for the vertical line
  if (is.numeric(thresh_obj)) {
    vline_value <- thresh_obj
  } else if (!is.null(thresh_obj$threshold_mean)) {
    vline_value <- thresh_obj$threshold_mean
  } else {
    stop("thresh_obj must be numeric or have a $threshold_mean element")
  }
  
  p <- ggplot(df, aes(x = !!biomarker_sym, y = !!RofC_sym, colour = factor(!!colour_sym))) +
    geom_point() +
    theme_bw() +
    scale_colour_manual(values = c("#007B82", "grey67")) +
    geom_hline(yintercept = rel_accum_thresh, linetype = "solid") + 
    geom_smooth(method = "gam", aes(group = 1), colour = "#FFD700") +
    theme(legend.position = "none")
  
  p2 <- p +     geom_vline(xintercept = vline_value, linetype = "dashed") 
  
  # Apply axis labels if provided
  if (!is.null(XLAB)) p <- p + xlab(XLAB)
  if (!is.null(YLAB)) p <- p + ylab(YLAB)
  
  return(list(withoutThresh = p, withThresh = p2))
}




bootstrap_get_Time_to_Positivity_amended <- function(df, PET_pos_threshold, id_name, time_name, value_name,
         num_bootstraps = 1000, bootstrap_percent = 0.8, degree = 3, printIter = TRUE) {
  df_res <- list()
  
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  # List to store bootstrapped results
  df_bs <- vector("list", num_bootstraps)
  
  # Perform bootstrap resampling
  for (i in 947:num_bootstraps) {
    if(printIter == TRUE){
      print(i)
    }
    # Sample IDs with replacement
    sampled_ids <- sample(seq(from = 1, to = length(unique(df[[ id_name]]))), 
                          size = floor(bootstrap_percent * length(unique(df[[id_name]]))), replace = TRUE)
    ID_vec <- unique(df[[id_name]])
    ID_vec <- ID_vec[sampled_ids]
    # Subset the data for the sampled IDs
    df <- data.frame(df)
    sampled_data <- df[df[[id_name]] %in% (ID_vec),] 
    
    # Apply get_Time_to_Positivity on the sampled data
    result <- get_Time_to_Positivity(sampled_data,id_name, time_name, value_name, PET_pos_threshold, degree)
    
    # Store the actual_predicted_val from this iteration
    df_bs[[i]] <- data.frame("val_bs" = result$actual_predicted_val, 
                             "time_bs" = result$Time_to_Positivity)
    interpolated_val <- vector("list", length(Time_Window))
    
    for(j in 1:length(Time_Window)){
      
      interpolated_val[[j]] <- approx( result$Time_to_Positivity, result$actual_predicted_val, as.numeric(Time_Window[j]))$y
    }
    df_res[[i]] <- data.frame("Time_Window" = Time_Window)
    
    df_res[[i]]$interpolated_val <- unlist(interpolated_val)
    
  }
  
  # Convert results to a matrix or dataframe for easier processing
  bootstrap_matrix <- do.call(rbind, df_res)
  mean_result <- data.frame(bootstrap_matrix) %>%
    group_by(Time_Window) %>%
    summarize(interpolated_val_mean = median(interpolated_val, na.rm = TRUE))
  
  sd_result <- data.frame(bootstrap_matrix) %>%
    group_by(Time_Window) %>%
    summarize(interpolated_val_sd = sd(interpolated_val, na.rm = TRUE))
  
  
  ci_calc <- merge(mean_result, sd_result, by = "Time_Window")
  
  # Combine into a dataframe of confidence intervals
  ci_df <- data.frame(
    Time_to_Positivity = ci_calc$Time_Window,
    Estimate =ci_calc$interpolated_val_mean,
    CI_Lower = ci_calc$interpolated_val_mean - 1.96 * ci_calc$interpolated_val_sd,
    CI_Upper = ci_calc$interpolated_val_mean + 1.96 * ci_calc$interpolated_val_sd
  )
  
  # Only compute adjustment if there's a crossing over the threshold
  # Only compute adjustment if there's a crossing over the threshold
  tryCatch({
    
    valid_idx <- which(!is.na(ci_df$Estimate) & !is.na(ci_df$Time_Window))
    
    if (length(valid_idx) < 2) {
      stop("Not enough valid points for interpolation")
    }
    
    est <- ci_df$Estimate[valid_idx]
    time <- ci_df$Time_Window[valid_idx]
    
    if (any(est < PET_pos_threshold) & any(est > PET_pos_threshold)) {
      
      adj <- approx(
        x = est,
        y = time,
        xout = PET_pos_threshold,
        ties = mean
      )$y
      
      if (length(adj) == 1 && !is.na(adj)) {
        ci_df$Time_to_Positivity <- ci_df$Time_Window - adj
      } else {
        stop("Interpolation returned empty or invalid adjustment")
      }
      
    } else {
      # No crossing: keep original Time_Window
      ci_df$Time_to_Positivity <- ci_df$Time_Window
    }
    
  }, error = function(e) {
    
    warning(
      "Failed to calculate Time_to_Positivity; returning unadjusted time. Reason: ",
      conditionMessage(e)
    )
    
    ci_df$Time_to_Positivity <- ci_df$Time_Window
  })
  
  return(ci_df[, c("Time_to_Positivity", "Estimate", "CI_Lower", "CI_Upper")])
  
}





