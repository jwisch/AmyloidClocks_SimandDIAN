# path to folder

dir_path <- "./Simulations20260424"

files <- list.files(
  path = dir_path,
  pattern = "^df_result_.*\\.RDS$",
  full.names = TRUE
)

grid <- c(-2 ,0, 2, 4, 6)

results_df <- map_dfr(files, function(f) {
  
  fname <- basename(f)
  
  matches <- str_match(
    fname,
    "df_result_Xi(-?[0-9]+\\.?[0-9]*)Omega(-?[0-9]+\\.?[0-9]*)Alpha(-?[0-9]+\\.?[0-9]*)_mean_within_above_(-?[0-9]+\\.?[0-9]*)\\.RDS"
  )
  
  if (is.na(matches[1])) {
    message("Skipping: ", fname)
    return(NULL)
  }
  
  Xi    <- as.numeric(matches[2])
  Omega <- as.numeric(matches[3])
  Alpha <- as.numeric(matches[4])
  error <- as.numeric(matches[5])
  
  dt <- as.data.table(readRDS(f))
  
  # ensure clean ordering (important for approx)
  setorder(dt, run, Time_to_Positivity)
  
  # interpolate per run
  interp_dt <- dt[, {
    approx_out <- approx(
      x = Time_to_Positivity,
      y = Estimate,
      xout = grid,
      rule = 2   # allows extrapolation at edges
    )
    
    .(Time_to_Positivity = approx_out$x,
      Estimate = approx_out$y)
  }, by = run]
  
  # attach metadata
  interp_dt[, `:=`(
    Xi = Xi,
    Omega = Omega,
    Alpha = Alpha,
    error = error
  )]
  
  interp_dt
})

plot_hist_density_at_time <- function(data, 
                                      time_vals = c(0, 2),
                                      binwidth = 0.01,
                                      color_vec = c("#007B82","#00BFC4"),
                                      title = "",
                                      scale_break = FALSE,
                                      scale_break1 = 0,
                                      scale_break2 = 0,
                                      scale_break_scale = 0.5) {
  
  
  
  # Subset once
  data_sub <- data[
    data$Time_to_Positivity %in% time_vals,
  ]
  
  data_sub$Time_to_Positivity <- factor(data_sub$Time_to_Positivity,
                                        levels = time_vals)
  
  color_vec <- color_vec[seq_along(time_vals)]
  
  # Name colors to match levels
  names(color_vec) <- time_vals
  
  p <- ggplot(data_sub, aes(x = Estimate, fill = Time_to_Positivity)) +
    
    
    # Density layer
    geom_density(
      aes(y = after_stat(density), color = Time_to_Positivity),
      linewidth = 1,
      adjust = 2,
      fill = NA  # ensures no shaded area
    ) + xlim(c(0, 10)) +
    
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

p1 <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0 & results_df$Omega == 1.5 & results_df$Alpha == -1 & results_df$error == 0,], 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "Xi = -, Omega = 1.5, Alpha = -1, Error = 0")

p2 <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0 & results_df$Omega == 1.5 & results_df$Alpha == -1 & results_df$error == 0.03,], 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.03")

p3 <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0 & results_df$Omega == 1.5 & results_df$Alpha == -1 & results_df$error == 0.07,], 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.07")



p4 <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0 & results_df$Omega == 1.5 & results_df$Alpha == -1 & results_df$error == 0.10,], 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.10")
p5 <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0 & results_df$Omega == 1.5 & results_df$Alpha == -1 & results_df$error == 0.15,], 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.15")
p6 <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0 & results_df$Omega == 1.5 & results_df$Alpha == -1 & results_df$error == 0.20,], 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.20")

lemon::grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, nrow = 1, ncol = 6)










p1a <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0.5 & results_df$Omega == 1 & results_df$Alpha == 1 & results_df$error == 0,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0")

p2a <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0.5 & results_df$Omega == 1 & results_df$Alpha == 1 & results_df$error == 0.03,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.03")

p3a <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0.5 & results_df$Omega == 1 & results_df$Alpha == 1 & results_df$error == 0.07,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.07")



p4a <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0.5 & results_df$Omega == 1 & results_df$Alpha == 1 & results_df$error == 0.10,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.10")
p5a <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0.5 & results_df$Omega == 1 & results_df$Alpha == 1 & results_df$error == 0.15,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.15")
p6a <- plot_hist_density_at_time(data = results_df[results_df$Xi == 0.5 & results_df$Omega == 1 & results_df$Alpha == 1 & results_df$error == 0.20,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.20")

lemon::grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, nrow = 1, ncol = 6)









p1b <- plot_hist_density_at_time(data = results_df[results_df$Xi == 1 & results_df$Omega == 0.5 & results_df$Alpha == 2 & results_df$error == 0,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0")

p2b <- plot_hist_density_at_time(data = results_df[results_df$Xi == 1 & results_df$Omega == 0.5 & results_df$Alpha == 2 & results_df$error == 0.03,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.03")

p3b <- plot_hist_density_at_time(data = results_df[results_df$Xi == 1 & results_df$Omega == 0.5 & results_df$Alpha == 2 & results_df$error == 0.07,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.07")



p4b <- plot_hist_density_at_time(data = results_df[results_df$Xi == 1 & results_df$Omega == 0.5 & results_df$Alpha == 2 & results_df$error == 0.10,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.10")
p5b <- plot_hist_density_at_time(data = results_df[results_df$Xi == 1 & results_df$Omega == 0.5 & results_df$Alpha == 2 & results_df$error == 0.15,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.15")
p6b <- plot_hist_density_at_time(data = results_df[results_df$Xi == 1 & results_df$Omega == 0.5 & results_df$Alpha == 2 & results_df$error == 0.20,], 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "grey10",  # darkest
                                   "grey30",  # 2nd darkest (your anchor)
                                   "grey50",
                                   "grey70",
                                   "grey90"   # lightest
                                 ),  # very light endpoint
                                 "Xi = -, Omega = 1.5, Alpha = -1, Error = 0.20")

lemon::grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, nrow = 1, ncol = 6)

lemon::grid_arrange_shared_legend(p1 + ggtitle(""), p2  + ggtitle(""), p3 + ggtitle(""), p4 + ggtitle(""), p5 + ggtitle(""), p6 + ggtitle(""),
                                  p1a + ggtitle(""), p2a + ggtitle(""), p3a + ggtitle(""), p4a + ggtitle(""), p5a + ggtitle(""), p6a + ggtitle(""),
                                  p1b + ggtitle(""), p2b + ggtitle(""), p3b + ggtitle(""), p4b + ggtitle(""), p5b + ggtitle(""), p6b + ggtitle(""),
                                  nrow = 3, ncol = 6)
graph2ppt(file = "./Figures/DistSpecificity_SyntheticData.pptx", width = 10, height = 7)
