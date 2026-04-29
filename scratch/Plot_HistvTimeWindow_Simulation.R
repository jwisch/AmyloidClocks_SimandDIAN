library(data.table)
library(ggplot2)
library(export)
# Extract all file names containing "ResultingCurves"
files <- list.files(
  path = "./RDS_Simulation_HistogramVarying/",
  pattern = "ResultingCurves",
  full.names = TRUE  # set to FALSE if you want just the file names
)

files
df <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_0_omega_0.5_xi_-1_ResultingCurves.RDS")
df1 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_0_omega_0.5_xi_0_ResultingCurves.RDS")
df2 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_0_omega_0.5_xi_1_ResultingCurves.RDS")

df3<- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_0.5_xi_-1_ResultingCurves.RDS")
df4 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_0.5_xi_0_ResultingCurves.RDS")
df5 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_0.5_xi_1_ResultingCurves.RDS")

df6 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_0_omega_2_xi_-1_ResultingCurves.RDS")
df7 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_0_omega_2_xi_0_ResultingCurves.RDS")
df8 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_0_omega_2_xi_1_ResultingCurves.RDS")

df9<- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_2_xi_-1_ResultingCurves.RDS")
df10 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_2_xi_0_ResultingCurves.RDS")
df11 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_2_xi_1_ResultingCurves.RDS")

ttp_seq  <- seq(from = -2, to = 10, by = 2)




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

interp_df <- get_interpolate(df)
interp_df1 <- get_interpolate(df1)
interp_df2 <- get_interpolate(df2)
interp_df3 <- get_interpolate(df3)
interp_df4 <- get_interpolate(df4)
interp_df5 <- get_interpolate(df5)
interp_df6 <- get_interpolate(df6)
interp_df7 <- get_interpolate(df7)
interp_df8 <- get_interpolate(df8)
interp_df9 <- get_interpolate(df9)
interp_df10 <- get_interpolate(df10)
interp_df11 <- get_interpolate(df11)



p1 <- plot_hist_density_at_time(data = interp_df, rel_accum = -100, 
                                time_vals = c(0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                   "#A93E36",  
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ) +
  ylim(c(0, 5)) + xlim(c(0, 10))

p2 <- plot_hist_density_at_time(data = interp_df1, rel_accum = -100, 
                          time_vals = c(-2 ,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  # darkest
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ) + ylim(c(0, 5)) + xlim(c(0, 10))
p3 <- plot_hist_density_at_time(data = interp_df2, rel_accum = -100, 
                          time_vals = c(-2 ,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  # darkest
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ) + ylim(c(0, 5)) + xlim(c(0, 10))



p4 <- plot_hist_density_at_time(data = interp_df3, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#A93E36",  # darkest
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ) +
  ylim(c(0, 5)) + xlim(c(0, 10))

p5 <- plot_hist_density_at_time(data = interp_df4, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#A93E36",  # darkest
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ) + ylim(c(0, 5)) + xlim(c(0, 10))
p6 <- plot_hist_density_at_time(data = interp_df5, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#A93E36",  # darkest
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ) + ylim(c(0, 5)) + xlim(c(0, 10))


p7 <- plot_hist_density_at_time(data = interp_df6, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) ) +
  ylim(c(0, 5)) + xlim(c(0, 10))

p8 <- plot_hist_density_at_time(data = interp_df7, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) ) + ylim(c(0, 5)) + xlim(c(0, 10))
p9 <- plot_hist_density_at_time(data = interp_df8, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) ) + ylim(c(0, 5)) + xlim(c(0, 10))




p10 <- plot_hist_density_at_time(data = interp_df9, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) ) +
  ylim(c(0, 5)) + xlim(c(0, 10))

p11 <- plot_hist_density_at_time(data = interp_df10, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) ) + ylim(c(0, 5)) + xlim(c(0, 10))
p12 <- plot_hist_density_at_time(data = interp_df11, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) ) + ylim(c(0, 5)) + xlim(c(0, 10))


gridExtra::grid.arrange(p1, p7, p2, p8, p3, p9, 
             p4, p10, p5, p11, p6, p12, nrow = 2, ncol = 6)
graph2ppt(file = "./Figures/DistributionSeperation_Simulation.pptx", width = 8, height = 7)

df1 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_-2_omega_2_xi_-1_ResultingCurves.RDS")
df2 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_-1_omega_1.5_xi_0_ResultingCurves.RDS")
df3<- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_1_omega_0.5_xi_0.5_ResultingCurves.RDS")
df4 <- readRDS("./RDS_Simulation_HistogramVarying/PiB_alpha_2_omega_0.5_xi_1_ResultingCurves.RDS")

interp_df1 <- get_interpolate(df1)
interp_df2 <- get_interpolate(df2)
interp_df3 <- get_interpolate(df3)
interp_df4 <- get_interpolate(df4)

p1 <- plot_hist_density_at_time(data = interp_df1, rel_accum = -100, 
                          time_vals = c(-2 ,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#8E4CC7",  # darkest
                            "#C77CFF",  # 2nd darkest (your anchor)
                            "#D79BFF",
                            "#E8C2FF",
                            "#F5E6FF"   # lightest
                          ) ) 

p2 <- plot_hist_density_at_time(data = interp_df2, rel_accum = -100, 
                          time_vals = c(-2 ,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#002F31", #(darkest)
                            "#007B82",# (second darkest — your anchor)
                            "#33A6AC",
                            "#7FCCCF",
                            "#CFEDEE"  # lightest
                          ) ) + xlim(c(0, 15))

p3 <- plot_hist_density_at_time(data = interp_df3, rel_accum = -100, 
                          time_vals = c(-2 ,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#004D16", #(darkest)
                            "#00BA38",# (second darkest — anchor)
                            "#33C95E",
                            "#80E09C",
                            "#CCF3D7"# (lightest)
                          ) ) + xlim(c(0, 15))

p4 <- plot_hist_density_at_time(data = interp_df4, rel_accum = -100, 
                          time_vals = c(-2 ,0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  # darkest
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ) + xlim(c(0, 15))

grid.arrange(p1, p2, p3, p4, nrow = 4)
graph2ppt(file = "./Figures/BetweenIndividualVariability_Individuallevelperfromance.pptx", width = 4, height = 7)

p1 <- plot_hist_density_at_time(data = interp_df1, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#8E4CC7",  # darkest
                                  "#C77CFF",  # 2nd darkest (your anchor)
                                  "#D79BFF",
                                  "#E8C2FF",
                                  "#F5E6FF"   # lightest
                                ) )  + theme(legend.position = "none")

p2 <- plot_hist_density_at_time(data = interp_df2, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#002F31", #(darkest)
                                  "#007B82",# (second darkest — your anchor)
                                  "#33A6AC",
                                  "#7FCCCF",
                                  "#CFEDEE"  # lightest
                                ) ) + xlim(c(0, 15)) + theme(legend.position = "none")

p3 <- plot_hist_density_at_time(data = interp_df3, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#004D16", #(darkest)
                                  "#00BA38",# (second darkest — anchor)
                                  "#33C95E",
                                  "#80E09C",
                                  "#CCF3D7"# (lightest)
                                ) ) + xlim(c(0, 15)) + theme(legend.position = "none")

p4 <- plot_hist_density_at_time(data = interp_df4, rel_accum = -100, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#A93E36",  # darkest
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ) + xlim(c(0, 15)) + theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, nrow = 4)
graph2ppt(file = "./Figures/BetweenIndividualVariability_Individuallevelperfromance_noLegend.pptx", width = 4, height = 7)

