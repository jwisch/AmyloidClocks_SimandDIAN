bs_DistSpecplot <- function(df, PET_pos_threshold, id_name, time_name, value_name,
                         num_bootstraps = 1000, bootstrap_percent = 0.8, degree = 3,  printIter = TRUE) {
  df_res <- list()
  
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  # List to store bootstrapped results
  df_bs <- vector("list", num_bootstraps)
  
  # Perform bootstrap resampling
  for (i in 1:num_bootstraps) {
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
  return(bootstrap_matrix)}


robust_bs_DistSpecplot <- function(
    df,
    PET_pos_threshold,
    id_name,
    time_name,
    value_name,
    num_bootstraps = 1000,
    bootstrap_percent = 0.8,
    degree = 3,
    printIter = TRUE
) {
  
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  unique_ids <- unique(df[[id_name]])
  n_ids <- length(unique_ids)
  
  df_res <- vector("list", num_bootstraps)
  
  # --- SAFE ITERATION ---
  safe_iteration <- function(i) {
    tryCatch({
      
      # Sample IDs
      sampled_idx <- sample(
        seq_len(n_ids),
        size = floor(bootstrap_percent * n_ids),
        replace = TRUE
      )
      sampled_ids <- unique_ids[sampled_idx]
      sampled_data <- df[df[[id_name]] %in% sampled_ids, ]
      
      # Model
      result <- get_Time_to_Positivity(
        sampled_data,
        id_name,
        time_name,
        value_name,
        PET_pos_threshold,
        degree
      )
      
      # Validate
      if (is.null(result) ||
          length(result$Time_to_Positivity) < 2 ||
          length(result$actual_predicted_val) < 2 ||
          all(is.na(result$Time_to_Positivity)) ||
          all(is.na(result$actual_predicted_val))) {
        return(NULL)
      }
      
      # Interpolation
      interpolated_val <- suppressWarnings(
        sapply(Time_Window, function(t) {
          approx(
            result$Time_to_Positivity,
            result$actual_predicted_val,
            xout = t,
            rule = 1,
            ties = mean
          )$y
        })
      )
      
      if (all(is.na(interpolated_val))) return(NULL)
      
      data.frame(
        Time_Window = Time_Window,
        interpolated_val = interpolated_val
      )
      
    }, error = function(e) {
      if (printIter) message("Iteration ", i, " failed: ", e$message)
      return(NULL)
    })
  }
  
  # --- RUN LOOP (this WILL go to num_bootstraps) ---
  for (i in seq_len(num_bootstraps)) {
    if (printIter) print(i)
    
    # try() catches anything tryCatch misses (important)
    res_i <- try(safe_iteration(i), silent = TRUE)
    
    if (inherits(res_i, "try-error") || is.null(res_i)) {
      df_res[[i]] <- NULL
    } else {
      df_res[[i]] <- res_i
    }
  }
  
  # --- Diagnostics ---
  num_failed <- sum(sapply(df_res, is.null))
  message("Completed all iterations.")
  message("Failed iterations: ", num_failed, " / ", num_bootstraps)
  
  df_res_clean <- df_res[!sapply(df_res, is.null)]
  
  if (length(df_res_clean) == 0) {
    warning("All iterations failed.")
    return(data.frame(
      Time_to_Positivity = numeric(0),
      Estimate = numeric(0),
      CI_Lower = numeric(0),
      CI_Upper = numeric(0)
    ))
  }
  
  # --- Aggregate ---
  bootstrap_matrix <- do.call(rbind, df_res_clean)
  
  
  return(bootstrap_matrix)
}


plot_hist_density_at_time <- function(data, rel_accum, 
                                      time_vals = c(0, 2),
                                      binwidth = 0.01,
                                      color_vec = c("#007B82","#00BFC4"),
                                      title = "",
                                      min_n = 100,
                                      scale_break = FALSE,
                                      scale_break1 = 0,
                                      scale_break2 = 0,
                                      scale_break_scale = 0.5) {
  
  
  
  # Subset once
  data_sub <- data[
    data$Time_Window %in% time_vals &
      data$interpolated_val > rel_accum,
  ]
  
  # Count rows per time window
  counts <- table(data_sub$Time_Window)
  
  # Keep only time windows meeting minimum N
  valid_times <- names(counts[counts >= min_n])
  
  if (length(valid_times) == 0) {
    stop("No time windows meet the minimum sample size requirement.")
  }
  
  # Filter data
  data_sub <- data_sub[data_sub$Time_Window %in% valid_times, ]
  
  # Preserve ordering from original time_vals
  valid_times <- intersect(time_vals, as.numeric(valid_times))
  data_sub$Time_Window <- factor(data_sub$Time_Window,
                                 levels = valid_times)
  
  # Match color vector length
  if (length(color_vec) < length(valid_times)) {
    stop("Not enough colors supplied for the number of valid time windows.")
  }
  color_vec <- color_vec[seq_along(valid_times)]
  
  # Name colors to match levels
  names(color_vec) <- valid_times
  
  p <- ggplot(data_sub, aes(x = interpolated_val, fill = Time_Window)) +
    
    # Histogram layer
    # geom_histogram(aes(y = after_stat(density)),
    #                binwidth = binwidth,
    #                position = "identity",
    #                alpha = 0.3,
    #                colour = NA) +
    
    # Density layer
    geom_density(
      aes(y = after_stat(density), color = Time_Window),
      linewidth = 1,
      adjust = 2,
      fill = NA  # ensures no shaded area
    ) +
    
    scale_fill_manual(values = color_vec) +
    scale_color_manual(values = color_vec) +
    
    theme_bw() +
    labs(title = title,
         x = "Z",
         y = NULL,
         fill = "Time Window",
         color = "Time Window") +
    
    theme(
      panel.grid = element_blank(),        # remove grid lines
      panel.border = element_blank(),      # remove full rectangle
      axis.ticks.y = element_blank(),      # remove y-axis ticks
      axis.text.y = element_blank(),       # remove y-axis labels
      axis.title.y = element_blank(),      # remove y-axis title
      axis.line.x = element_line(color = "black"),  # only bottom line
      legend.position = "bottom"
    )
  # Optional broken axis
  if (scale_break) {
    p <- p + ggbreak::scale_y_break(
      c(scale_break1, scale_break2),
      scales = scale_break_scale
    )
  }
  
  return(p)
}
