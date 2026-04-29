files <- list.files(
  path = "./RDS_Simulation_v02_fullBS/",
  pattern = "df_result",
  full.names = TRUE  # set to FALSE if you want just the file names
)

files


df_03 <- readRDS("./RDS_Simulation_v02_fullBS/df_result_Xi0Omegapoint5Alpha2_mean_within_above_0.03.RDS")
df_07 <- readRDS("./RDS_Simulation_v02_fullBS/df_result_Xi0Omegapoint5Alpha2_mean_within_above_0.07.RDS")
df_10 <- readRDS("./RDS_Simulation_v02_fullBS/df_result_Xi0Omegapoint5Alpha2_mean_within_above_0.1.RDS")
df_15 <- readRDS("./RDS_Simulation_v02_fullBS/df_result_Xi0Omegapoint5Alpha2_mean_within_above_0.15.RDS")
df_20 <- readRDS("./RDS_Simulation_v02_fullBS/df_result_Xi0Omegapoint5Alpha2_mean_within_above_0.2.RDS")

get_interpolate <- function(df){
  # Use data.table by group
  interp_df <- setDT(df)[, {
    interp <- approx(x = Time_to_Positivity,
                     y = Estimate,
                     xout = ttp_seq,
                     rule = 2)
    .(Time_to_Positivity = interp$x,
      interpolated_val   = interp$y)
  }, by = run]
  
  colnames(interp_df)[2] <- "Time_Window"
  return(interp_df)
  
}

ttp_seq  <- seq(from = -2, to = 10, by = 2)

interp_df_03 <- get_interpolate(df_03)
interp_df_07 <- get_interpolate(df_07)
interp_df_10 <- get_interpolate(df_10)
interp_df_15 <- get_interpolate(df_15)
interp_df_20 <- get_interpolate(df_20)



plot_hist_density_at_time <- function(data, rel_accum, 
                                      time_vals = c(0, 2),
                                      binwidth = 0.01,
                                      color_vec = c("#007B82","#00BFC4"),
                                      title = "",
                                      min_n = 0,
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
  # 
  # if (length(valid_times) == 0) {
  #   stop("No time windows meet the minimum sample size requirement.")
  # }
  
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
  # # Optional broken axis
  # if (scale_break) {
  #   p <- p + ggbreak::scale_y_break(
  #     c(scale_break1, scale_break2),
  #     scales = scale_break_scale
  #   )
  # }
  # 
  return(p)
}

p1 <- plot_hist_density_at_time(data = interp_df_03, rel_accum = -100, 
                          time_vals = c(-2, 0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ) + xlim(c(0, 7))

p2 <- plot_hist_density_at_time(data = interp_df_07, rel_accum = -100, 
                          time_vals = c(-2,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ) + xlim(c(0, 7))

p3 <- plot_hist_density_at_time(data = interp_df_10, rel_accum = -100, 
                          time_vals = c(-2,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) )  + xlim(c(0, 7))

p4 <- plot_hist_density_at_time(data = interp_df_15, rel_accum = -100, 
                          time_vals = c(-2,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ) + xlim(c(0, 7))

p5 <- plot_hist_density_at_time(data = interp_df_20, rel_accum = -100, 
                          time_vals = c(-2,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) )  + xlim(c(0, 7))

lemon::grid_arrange_shared_legend(p1, p2, p3, p4, p5, nrow = 5, ncol = 1)
graph2ppt(file = "./Figures/DistributionSeperation_Simulation_withinIndividualError.pptx", width = 3.5, height = 7)



# Convert to data.table for convenience
dt <- as.data.table(interp_df_03)

# List of unique time windows
time_windows <- sort(unique(dt$Time_Window))

# Compute density estimates per Time_Window
dens_list <- lapply(time_windows, function(tw) {
  d <- density(dt[Time_Window == tw, interpolated_val])
  list(time_window = tw, x = d$x, y = d$y)
})

x0 <- 2

prob_y <- sapply(dens_list, function(d) {
  # Find density at x0 using linear interpolation
  approx(d$x, d$y, xout = x0)$y
})

names(prob_y) <- sapply(dens_list, function(d) d$time_window)
prob_y

probabilities <- prob_y / sum(prob_y, na.rm = TRUE)
probabilities
