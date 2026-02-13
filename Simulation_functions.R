#This function generates synthetic data based on my simplest conception of the problem
#Note that within the baseline_group assignment You can set the probability of the
#starting value being below the meaningful pathological level, within the key conversion window,
#within the range of accumualtion, or above the meaningful pathological level

generate_synthetic_data <- function(
    n_ids = 100,
    n_visits = 5,
    avg_interval = 1,            # average years between visits
    interval_noise = 0.2,        # variability in visit spacing
    minMeaningfulZ = 20,
    AposThresh = 18,
    maxMeaningfulZ = 100,
    annual_rate = 5,             # population-average annual change
    var_annual_rate = 1,         # between-individual variability in annual rate
    # test-retest variability parameters
    mean_within_above = 0.08,    # 8% above threshold
    sd_within_above   = 0.07,
    mean_within_below = 0.044,   # 4.4% below threshold
    sd_within_below   = 0.042
) {
  
  df_list <- list()
  
  for (id in 1:n_ids) {
    
    # ---- baseline CL assignment ----
    baseline_group <- sample(c("below", "conversionWindow", "between", "above"), 
                             size = 1, 
                             prob = c(0.3, 0.5, 0.2, 0.1))
    
    if (baseline_group == "below") {
      baseline_Z <- runif(1, min = 0, max = max(1e-6, minMeaningfulZ - 0.05))
    } else if (baseline_group == "conversionWindow") {
      baseline_Z <- runif(1, min = minMeaningfulZ, max = AposThresh + 0.05)
    } else if (baseline_group == "between") {
      baseline_Z <- runif(1, min = AposThresh + 0.05, max = maxMeaningfulZ - 0.05)
    } else {
      baseline_Z <- runif(1, min = maxMeaningfulZ, max = maxMeaningfulZ + 3)
    }

    # ---- individual-specific slope ----
    indiv_rate <- rnorm(1, mean = annual_rate, sd = var_annual_rate)
    
    # ---- visit times ----
    intervals <- abs(rnorm(n_visits - 1, mean = avg_interval, sd = interval_noise))
    timepoints <- c(0, cumsum(intervals))
    
    # ---- CL trajectory ----
    Z_values <- numeric(length(timepoints))
    Z_values[1] <- baseline_Z
    
    for (v in 2:length(timepoints)) {
      years <- timepoints[v] - timepoints[v - 1]
      
      if (Z_values[v - 1] < minMeaningfulZ) {
        noise_frac <- rnorm(1, mean = mean_within_below, sd = sd_within_below)
      } else {
        noise_frac <- rnorm(1, mean = mean_within_above, sd = sd_within_above)
      }
      
      # force non-negative
      noise_frac <- pmax(noise_frac, 0)
      
      # safe sd
      noise_sd <- abs(noise_frac) * max(Z_values[v - 1], 1e-6)
      
      noise <- rnorm(1, mean = 0, sd = noise_sd)
      
      if (Z_values[v - 1] < minMeaningfulZ) {
        Z_values[v] <- Z_values[v - 1] + noise
      } else {
        Z_values[v] <- Z_values[v - 1] + years * indiv_rate + noise
      }
    }
    
    df_list[[id]] <- data.frame(
      ID = id,
      TimefromBaseline = timepoints,
      Z = Z_values,
      indiv_rate = indiv_rate
    )
  }
  
  df <- do.call(rbind, df_list)
  return(df)
}



robust_bootstrap_TTP <- function(df, PET_pos_threshold, id_name, time_name, value_name,
                                 num_bootstraps = 1000, bootstrap_percent = 0.8,
                                 degree = 3, printIter = TRUE) {
  
  df_res <- list()
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  df_bs <- vector("list", num_bootstraps)
  
  successful_runs <- 0  # track how many bootstrap iterations produce valid results
  
  for (i in 1:num_bootstraps) {
    if (printIter) cat("Bootstrap iteration:", i, "\n")
    
    sampled_ids <- sample(seq_len(length(unique(df[[id_name]]))),
                          size = floor(bootstrap_percent * length(unique(df[[id_name]]))),
                          replace = TRUE)
    ID_vec <- unique(df[[id_name]])[sampled_ids]
    sampled_data <- df[df[[id_name]] %in% ID_vec, ]
    
    # Attempt model fit
    result <- tryCatch({
      get_Time_to_Positivity(sampled_data, id_name, time_name, value_name, PET_pos_threshold, degree)
    }, error = function(e) NULL)
    
    # Skip if model failed or returned invalid data
    if (is.null(result) || any(is.na(result$actual_predicted_val)) || any(is.na(result$Time_to_Positivity))) {
      next
    }
    
    # Mark successful bootstrap
    successful_runs <- successful_runs + 1
    
    # Interpolate results over the common Time_Window
    interpolated_val <- sapply(Time_Window, function(tw) {
      tryCatch({
        valid_idx <- which(!is.na(result$Time_to_Positivity) & !is.na(result$actual_predicted_val))
        if (length(valid_idx) < 2) {
          return(NA)
        } else {
          return(approx(result$Time_to_Positivity[valid_idx],
                        result$actual_predicted_val[valid_idx],
                        xout = tw)$y)
        }
      }, error = function(e) {
        return(NA)
      })
    })
    
    df_res[[i]] <- data.frame(Time_Window = Time_Window,
                              interpolated_val = interpolated_val)
  }
  
  # Check success rate
  success_rate <- successful_runs / num_bootstraps
  message(sprintf("Successful bootstrap rate: %.1f%%", success_rate * 100))
  
  if (success_rate < 0.9) {
    warning(sprintf("Only %.1f%% of bootstrap fits succeeded — skipping CI estimation.", success_rate * 100))
    return(data.frame(Time_to_Positivity = NA, Estimate = NA, CI_Lower = NA, CI_Upper = NA))
  }
  
  # Combine successful results
  bootstrap_matrix <- do.call(rbind, df_res)
  
  mean_result <- bootstrap_matrix %>%
    dplyr::group_by(Time_Window) %>%
    dplyr::summarize(interpolated_val_mean = median(interpolated_val, na.rm = TRUE), .groups = "drop")
  
  sd_result <- bootstrap_matrix %>%
    dplyr::group_by(Time_Window) %>%
    dplyr::summarize(interpolated_val_sd = sd(interpolated_val, na.rm = TRUE), .groups = "drop")
  
  ci_calc <- merge(mean_result, sd_result, by = "Time_Window")
  
  ci_df <- data.frame(
    Time_to_Positivity = ci_calc$Time_Window,
    Estimate = ci_calc$interpolated_val_mean,
    CI_Lower = ci_calc$interpolated_val_mean - 1.96 * ci_calc$interpolated_val_sd,
    CI_Upper = ci_calc$interpolated_val_mean + 1.96 * ci_calc$interpolated_val_sd
  )
  
  # Adjust to align PET_pos_threshold = 0
  tryCatch({
    if (any(ci_df$Estimate < PET_pos_threshold) & any(ci_df$Estimate > PET_pos_threshold)) {
      adjustment <- approx(ci_df$Estimate, ci_df$Time_to_Positivity, xout = PET_pos_threshold)$y
      ci_df$Time_to_Positivity <- ci_df$Time_to_Positivity - as.numeric(adjustment)
    }
  }, error = function(e) {
    warning("Failed to calculate Time_to_Positivity: ", conditionMessage(e))
  })
  
  return(ci_df[, c("Time_to_Positivity", "Estimate", "CI_Lower", "CI_Upper")])
}


generate_synthetic_data_opt2 <- function(
    n_ids = 100,
    n_visits = 5,
    avg_interval = 1,            # average years between visits
    interval_noise = 0.2,        # variability in visit spacing
    minMeaningfulZ = 20,
    AposThresh = 18,
    maxMeaningfulZ = 100,
    annual_rate = 5,             # population-average annual change
    var_annual_rate = 1,         # between-individual variability in annual rate
    # test-retest variability parameters
    mean_within_above = 0.08,    # 8% above threshold
    sd_within_above   = 0.07,
    mean_within_below = 0.044,   # 4.4% below threshold
    sd_within_below   = 0.042,
    plateau_decay = 0.15         # fraction of slope retained above plateau
) {
  
  df_list <- list()
  
  for (id in 1:n_ids) {
    
    # ---- baseline CL assignment ----
    baseline_group <- sample(c("below", "conversionWindow", "between", "above"), 
                             size = 1, 
                             prob = c(0.12, 0.75, 0.10, 0.03))
    
    if (baseline_group == "below") {
      baseline_Z <- runif(1, min = 0, max = max(1e-6, minMeaningfulZ - 0.05))
    } else if (baseline_group == "conversionWindow") {
      baseline_Z <- runif(1, min = minMeaningfulZ, max = AposThresh + 0.05)
    } else if (baseline_group == "between") {
      baseline_Z <- runif(1, min = AposThresh + 0.05, max = maxMeaningfulZ - 0.05)
    } else {
      baseline_Z <- runif(1, min = maxMeaningfulZ, max = maxMeaningfulZ + 3)
    }
    
    # ---- individual-specific slope ----
    indiv_rate <- rnorm(1, mean = annual_rate, sd = var_annual_rate)
    
    # ---- visit times ----
    intervals <- abs(rnorm(n_visits - 1, mean = avg_interval, sd = interval_noise))
    timepoints <- c(0, cumsum(intervals))
    
    # ---- CL trajectory ----
    Z_values <- numeric(length(timepoints))
    Z_values[1] <- baseline_Z
    
    for (v in 2:length(timepoints)) {
      years <- timepoints[v] - timepoints[v - 1]
      prevZ <- Z_values[v - 1]
      
      # --- within-subject noise logic ---
      if (prevZ < minMeaningfulZ) {
        noise_frac <- rnorm(1, mean = mean_within_below, sd = sd_within_below)
      } else if (prevZ < maxMeaningfulZ) {
        noise_frac <- rnorm(1, mean = mean_within_above, sd = sd_within_above)
      } else {
        # above plateau → return to below-threshold noise
        noise_frac <- rnorm(1, mean = mean_within_below, sd = sd_within_below)
      }
      noise_frac <- pmax(noise_frac, 0)
      noise_sd <- abs(noise_frac) * max(prevZ, 1e-6)
      noise <- rnorm(1, mean = 0, sd = noise_sd)
      
      # --- trajectory dynamics ---
      if (prevZ < minMeaningfulZ) {
        # below threshold: flat, only noise
        dZ <- noise
      } else if (prevZ < maxMeaningfulZ) {
        # normal growth region
        dZ <- years * indiv_rate + noise
      } else {
        # above plateau → flatten trajectory
        # keep small residual slope so it doesn’t freeze completely
        dZ <- years * indiv_rate * plateau_decay + noise
      }
      
      # constrain not to overshoot too high
      Z_values[v] <- prevZ + dZ
      if (Z_values[v] > maxMeaningfulZ + 5) {
        Z_values[v] <- maxMeaningfulZ + rnorm(1, 0, 0.5)  # gentle bound
      }
    }
    
    df_list[[id]] <- data.frame(
      ID = id,
      TimefromBaseline = timepoints,
      Z = Z_values,
      indiv_rate = indiv_rate
    )
  }
  
  df <- do.call(rbind, df_list)
  return(df)
}


estimate_simulation_params <- function(df, id_col, age_col, biomarker_col, minMeaningfulZ, AposThresh) {
  
  df <- df %>%
    dplyr::select(all_of(c(id_col, age_col, biomarker_col))) %>%
    dplyr::rename(ID = !!id_col, Age = !!age_col, Z = !!biomarker_col) %>%
    dplyr::filter(!is.na(Z), !is.na(Age))  # remove missing
  
  # --- Basic summaries ---
  n_ids <- n_distinct(df$ID)
  
  # average interval and noise
  intervals_df <- df %>%
    group_by(ID) %>%
    summarize(
      interval_mean = ifelse(n() > 1, mean(diff(sort(Age))), NA_real_),
      interval_sd   = ifelse(n() > 2, sd(diff(sort(Age))), NA_real_)
    )
  
  avg_interval   <- mean(intervals_df$interval_mean, na.rm = TRUE)
  interval_noise <- mean(intervals_df$interval_sd, na.rm = TRUE)
  
  # --- Slopes: only for IDs with ≥2 valid points ---
  slopes <- df %>%
    group_by(ID) %>%
    filter(n() > 1) %>%
    summarize(
      slope = {
        sub <- cur_data()
        if (sum(!is.na(sub$Z)) >= 2) {
          coef(lm(Z ~ Age, data = sub))[2]
        } else {
          NA_real_
        }
      },
      meanZ = mean(Z, na.rm = TRUE)
    )
  
  # --- Limit to above-threshold region ---
  slopes_above <- slopes %>% filter(meanZ >= minMeaningfulZ)
  annual_rate  <- mean(slopes_above$slope, na.rm = TRUE)
  var_annual_rate <- sd(slopes_above$slope, na.rm = TRUE)
  
  tibble(
    n_ids = n_ids,
    avg_interval = avg_interval,
    interval_noise = interval_noise,
    annual_rate = annual_rate,
    var_annual_rate = var_annual_rate
  )
}


# --- robust synthetic data generator ---
generate_synthetic_data_with_noise <- function(
    n_ids = 100,
    n_visits = 5,
    avg_interval = 1,
    interval_noise = 0.2,
    minMeaningfulZ = 20,
    AposThresh = 18,
    maxMeaningfulZ = 100,
    annual_rate = 5,
    var_annual_rate = 1,
    mean_within_above = 0.08,
    sd_within_above   = 0.07,
    mean_within_below = 0.044,
    sd_within_below   = 0.042,
    plateau_decay = 0.15,
    prob_neg = 0.5,
    xi = NULL,    # skew-normal location
    omega = NULL, # skew-normal scale
    alpha = NULL,  # skew-normal shape
    seed = 6432
) {
  # --- Dependencies ---
  if (!requireNamespace("sn", quietly = TRUE)) {
    stop("Package 'sn' is required. Please install it with install.packages('sn').")
  }
  set.seed(seed)
  # --- Sanity checks ---
  if (AposThresh <= minMeaningfulZ) warning("AposThresh <= minMeaningfulZ")
  if (maxMeaningfulZ <= AposThresh) warning("maxMeaningfulZ <= AposThresh")

  # epsilon to prevent zero-width intervals
  eps <- 1e-6
  
  # --- Helper: safe runif ---
  safe_runif <- function(n, minv, maxv) {
    if (is.na(minv) || is.na(maxv)) return(NA_real_)
    if (minv == maxv) return(rep(minv, n))
    if (minv > maxv) {
      tmp <- minv
      minv <- maxv
      maxv <- tmp
    }
    runif(n, minv, maxv)
  }
  
    
  df_list <- vector("list", n_ids)
  
  for (id in seq_len(n_ids)) {
    # --- baseline assignment ---
    baseline_group <- sample(c("below", "conversionWindow", "between", "above"),
                             size = 1, prob = c(0.12, 0.75, 0.10, 0.03))
    
    baseline_Z <- switch(
      baseline_group,
      below = safe_runif(1, 0, max(eps, minMeaningfulZ - 0.05)),
      conversionWindow = safe_runif(1, minMeaningfulZ, AposThresh + 0.05),
      between = safe_runif(1, AposThresh + 0.05, maxMeaningfulZ - 0.05),
      above = safe_runif(1, maxMeaningfulZ, maxMeaningfulZ + 3)
    )
    
    # --- draw indiv_rate ---
    if (!is.null(xi) && !is.null(omega) && !is.null(alpha)) {
      indiv_rate <- sn::rsn(1, xi = xi, omega = omega, alpha = alpha)
    } else {
      indiv_rate <- rnorm(1, mean = annual_rate, sd = var_annual_rate)
    }
    
    # --- timepoints ---
    intervals <- abs(rnorm(n_visits - 1, mean = avg_interval, sd = interval_noise))
    timepoints <- c(0, cumsum(intervals))
    
    Z_values <- numeric(length(timepoints))
    Z_true <- numeric(length(timepoints))
    noise_vec <- numeric(length(timepoints))
    Z_values[1] <- Z_true[1] <- baseline_Z
    
    for (v in 2:length(timepoints)) {
      prevZ_true <- Z_true[v - 1]
      prevZ_obs <- Z_values[v - 1]
      years <- timepoints[v] - timepoints[v - 1]
      
      # --- noise ---
      if (prevZ_obs < minMeaningfulZ) {
        noise_frac <- rnorm(1, mean_within_below, sd_within_below)
      } else if (prevZ_obs < maxMeaningfulZ) {
        noise_frac <- rnorm(1, mean_within_above, sd_within_above)
      } else {
        noise_frac <- rnorm(1, mean_within_below, sd_within_below)
      }
      noise_frac <- pmax(noise_frac, 0)
      noise_sd <- abs(noise_frac) * max(prevZ_obs, 1e-6)
      noise <- rnorm(1, 0, noise_sd)
      noise <- ifelse(noise < 0 & runif(1) > prob_neg, 0, noise)
      
      # --- deterministic "true" change ---
      if (prevZ_true < minMeaningfulZ) {
        dZ_true <- 0
      } else if (prevZ_true < maxMeaningfulZ) {
        dZ_true <- years * indiv_rate
      } else {
        dZ_true <- years * indiv_rate * plateau_decay
      }
      
      Z_true[v] <- min(maxMeaningfulZ + 5, prevZ_true + dZ_true)
      Z_values[v] <- min(maxMeaningfulZ + 5, prevZ_true + dZ_true + noise)
      noise_vec[v] <- noise
    }
    
    df_list[[id]] <- data.frame(
      ID = rep(id, length(timepoints)),
      TimefromBaseline = timepoints,
      Z_noisy = Z_values,
      Z_true = Z_true,
      indiv_rate = rep(indiv_rate, length(timepoints)),
      noise = noise_vec
    )
  } # <-- end for loop over IDs
  
  # --- combine all IDs ---
  df <- do.call(rbind, df_list)
  
  # --- identify converters ---
  cross_info <- df %>%
    dplyr::group_by(ID) %>%
    dplyr::summarize(
      start_below = first(Z_noisy) < AposThresh,
      ever_cross = any(Z_noisy > AposThresh),
      first_cross = ifelse(ever_cross, min(TimefromBaseline[Z_noisy > AposThresh]), NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::filter(start_below & ever_cross)
  
  converters <- cross_info %>% dplyr::filter(ever_cross & !is.na(first_cross))
  
  return(list(df = df, converters = converters))
} # <-- end function




# generate_synthetic_data_with_noise <- function(
#     n_ids = 100,
#     n_visits = 5,
#     avg_interval = 1,
#     interval_noise = 0.2,
#     minMeaningfulZ = 20,
#     AposThresh = 18,
#     maxMeaningfulZ = 100,
#     annual_rate = 5,
#     var_annual_rate = 1,
#     mean_within_above = 0.08,
#     sd_within_above   = 0.07,
#     mean_within_below = 0.044,
#     sd_within_below   = 0.042,
#     plateau_decay = 0.15,
#     prob_neg = 0.5,
#     xi = NULL,    # skew-normal location
#     omega = NULL, # skew-normal scale
#     alpha = NULL, # skew-normal shape
#     seed = 6432
# ) {
#   # --- Dependencies ---
#   if (!requireNamespace("sn", quietly = TRUE)) {
#     stop("Package 'sn' is required. Please install it with install.packages('sn').")
#   }
#   set.seed(seed)
#   
#   # --- Sanity checks ---
#   if (AposThresh <= minMeaningfulZ)
#     warning("AposThresh <= minMeaningfulZ — ranges will be adjusted to avoid invalid sampling.")
#   if (maxMeaningfulZ <= AposThresh)
#     warning("maxMeaningfulZ <= AposThresh — ranges will be adjusted to avoid invalid sampling.")
#   
#   # epsilon to prevent zero-width intervals
#   eps <- 1e-6
#   
#   # --- Helper: safe runif ---
#   safe_runif <- function(n, minv, maxv) {
#     if (is.na(minv) || is.na(maxv)) return(NA_real_)
#     if (minv == maxv) return(rep(minv, n))
#     if (minv > maxv) {
#       tmp <- minv
#       minv <- maxv
#       maxv <- tmp
#     }
#     runif(n, minv, maxv)
#   }
#   
#   df_list <- vector("list", n_ids)
#   
#   for (id in seq_len(n_ids)) {
#     # --- baseline assignment ---
#     baseline_group <- sample(c("below", "conversionWindow", "between", "above"),
#                              size = 1, prob = c(0.12, 0.75, 0.10, 0.03))
#     
#     baseline_Z <- switch(
#       baseline_group,
#       below = safe_runif(1, 0, max(eps, minMeaningfulZ - 0.05)),
#       conversionWindow = safe_runif(1, minMeaningfulZ, AposThresh + 0.05),
#       between = safe_runif(1, AposThresh + 0.05, maxMeaningfulZ - 0.05),
#       above = safe_runif(1, maxMeaningfulZ, maxMeaningfulZ + 3)
#     )
#     
#     # --- draw indiv_rate ---
#     if (!is.null(xi) && !is.null(omega) && !is.null(alpha)) {
#       indiv_rate <- sn::rsn(1, xi = xi, omega = omega, alpha = alpha)
#     } else {
#       indiv_rate <- rnorm(1, mean = annual_rate, sd = var_annual_rate)
#     }
#     
#     # --- timepoints ---
#     intervals <- abs(rnorm(n_visits - 1, mean = avg_interval, sd = interval_noise))
#     timepoints <- c(0, cumsum(intervals))
#     
#     Z_values <- numeric(length(timepoints))
#     Z_true <- numeric(length(timepoints))
#     noise_vec <- numeric(length(timepoints))
#     Z_values[1] <- Z_true[1] <- baseline_Z
#     
#     for (v in 2:length(timepoints)) {
#       prevZ_true <- Z_true[v - 1]
#       prevZ_obs <- Z_values[v - 1]
#       years <- timepoints[v] - timepoints[v - 1]
#       
#       # --- noise ---
#       if (is.na(prevZ_obs)) {
#         noise_frac <- mean_within_below
#       } else if (prevZ_obs < minMeaningfulZ) {
#         noise_frac <- rnorm(1, mean_within_below, sd_within_below)
#       } else if (prevZ_obs < maxMeaningfulZ) {
#         noise_frac <- rnorm(1, mean_within_above, sd_within_above)
#       } else {
#         noise_frac <- rnorm(1, mean_within_below, sd_within_below)
#       }
#       
#       noise_frac <- pmax(noise_frac, 0)
#       noise_sd <- abs(noise_frac) * max(prevZ_obs, eps)
#       noise <- rnorm(1, 0, noise_sd)
#       noise <- ifelse(noise < 0 & runif(1) > prob_neg, 0, noise)
#       
#       # --- deterministic "true" change ---
#       if (is.na(prevZ_true) || prevZ_true < minMeaningfulZ) {
#         dZ_true <- 0
#       } else if (prevZ_true < maxMeaningfulZ) {
#         dZ_true <- years * indiv_rate
#       } else {
#         dZ_true <- years * indiv_rate * plateau_decay
#       }
#       
#       Z_true[v] <- min(maxMeaningfulZ + 5, prevZ_true + dZ_true)
#       Z_values[v] <- min(maxMeaningfulZ + 5, prevZ_true + dZ_true + noise)
#       noise_vec[v] <- noise
#     }
#     
#     df_list[[id]] <- data.frame(
#       ID = rep(id, length(timepoints)),
#       TimefromBaseline = timepoints,
#       Z_noisy = Z_values,
#       Z_true = Z_true,
#       indiv_rate = rep(indiv_rate, length(timepoints)),
#       noise = noise_vec
#     )
#   }
#   
#   df <- do.call(rbind, df_list)
#   
#   cross_info <- df %>%
#     dplyr::group_by(ID) %>%
#     dplyr::summarize(
#       start_below = first(Z_noisy) < AposThresh,
#       ever_cross = any(Z_noisy > AposThresh, na.rm = TRUE),
#       first_cross = ifelse(ever_cross, min(TimefromBaseline[Z_noisy > AposThresh], na.rm = TRUE), NA_real_),
#       .groups = "drop"
#     ) %>%
#     dplyr::filter(start_below & ever_cross)
#   
#   converters <- cross_info %>% dplyr::filter(ever_cross & !is.na(first_cross))
#   
#   return(list(df = df, converters = converters))
# }
# 

# generate_synthetic_data_with_noise <- function(
#     n_ids = 100,
#     n_visits = 5,
#     avg_interval = 1,
#     interval_noise = 0.2,
#     minMeaningfulZ = 20,
#     AposThresh = 18,
#     maxMeaningfulZ = 100,
#     mean_slope_below = 0.044,
#     sd_slope_below = 0.042,
#     mean_slope_between = 0.08,
#     sd_slope_between = 0.07,
#     mean_slope_above = 0.05,
#     sd_slope_above = 0.03,
#     plateau_decay = 0.15,
#     prob_neg_below = 0.2,
#     prob_neg_between = 0.25,
#     prob_neg_above = 0.1
# ) {
#   
#   # enforce ordering
#   if (AposThresh <= minMeaningfulZ) warning("Problematic Structure: AposThresh <= minMeaningfulZ")
#   if (maxMeaningfulZ <= AposThresh) warning("Problematic Structure: maxMeaningfulZ <= AposThresh")
#   
#   df_list <- vector("list", n_ids)
#   
#   for (id in seq_len(n_ids)) {
#     # Determine baseline Z group
#     baseline_group <- sample(c("below", "conversionWindow", "between", "above"),
#                              size = 1,
#                              prob = c(0.12, 0.75, 0.10, 0.03))
#     
#     baseline_Z <- switch(baseline_group,
#                          below = runif(1, 0, max(1e-6, minMeaningfulZ - 0.05)),
#                          conversionWindow = runif(1, minMeaningfulZ, AposThresh + 0.05),
#                          between = runif(1, AposThresh + 0.05, maxMeaningfulZ - 0.05),
#                          above = runif(1, maxMeaningfulZ, maxMeaningfulZ + 3)
#     )
#     
#     # Select individual slope for each region
#     indiv_rate_below   <- rnorm(1, mean = mean_slope_below,   sd = sd_slope_below)
#     indiv_rate_between <- rnorm(1, mean = mean_slope_between, sd = sd_slope_between)
#     indiv_rate_above   <- rnorm(1, mean = mean_slope_above,   sd = sd_slope_above)
#     
#     intervals <- abs(rnorm(n_visits - 1, mean = avg_interval, sd = interval_noise))
#     timepoints <- c(0, cumsum(intervals))
#     
#     Z_values <- numeric(length(timepoints))
#     noise_vec <- numeric(length(timepoints))
#     Z_values[1] <- baseline_Z
#     
#     for (v in 2:length(timepoints)) {
#       prevZ <- Z_values[v - 1]
#       years <- timepoints[v] - timepoints[v - 1]
#       
#       # Assign rate and negative probability by region
#       if (prevZ < minMeaningfulZ) {
#         slope_mean <- indiv_rate_below
#         slope_sd   <- sd_slope_below
#         prob_neg   <- prob_neg_below
#       } else if (prevZ < AposThresh) {
#         slope_mean <- indiv_rate_between
#         slope_sd   <- sd_slope_between
#         prob_neg   <- prob_neg_between
#       } else {
#         slope_mean <- indiv_rate_above
#         slope_sd   <- sd_slope_above
#         prob_neg   <- prob_neg_above
#       }
#       
#       # Generate mostly positive noise, correct magnitude
#       raw <- abs(rnorm(1, mean = slope_mean, sd = slope_sd)) - slope_sd * sqrt(2/pi)
#       sign <- ifelse(runif(1) < prob_neg, -1, 1)
#       noise <- sign * raw
#       
#       # True trajectory increment
#       if (prevZ < minMeaningfulZ) {
#         dZ <- noise
#       } else if (prevZ < maxMeaningfulZ) {
#         dZ <- years * noise
#       } else {
#         dZ <- years * noise * plateau_decay
#       }
#       
#       Z_values[v] <- min(maxMeaningfulZ + 5, prevZ + dZ)
#       noise_vec[v] <- noise
#     }
#     
#     df_list[[id]] <- data.frame(
#       ID = id,
#       TimefromBaseline = timepoints,
#       Z = Z_values,
#       noise = noise_vec,
#       indiv_rate_below = indiv_rate_below,
#       indiv_rate_between = indiv_rate_between,
#       indiv_rate_above = indiv_rate_above
#     )
#   }
#   
#   df <- do.call(rbind, df_list)
#   
#   # --- identify converters ---
#   cross_info <- df %>%
#     dplyr::group_by(ID) %>%
#     dplyr::summarize(
#       start_below = first(Z) < AposThresh,
#       ever_cross = any(Z > AposThresh),
#       first_cross = ifelse(ever_cross, min(TimefromBaseline[Z > AposThresh]), NA_real_),
#       .groups = "drop"
#     ) %>%
#     dplyr::filter(start_below & ever_cross)
#   
#   converters <- cross_info %>% dplyr::filter(ever_cross & !is.na(first_cross))
#   
#   return(list(df = df, converters = converters))
# }
