#This script generates Figures 2A and 2B together.
##Note that when using grid.arrange to stack figures 2A and 2B, for some reason
##The breaks that were designed for FIgure 2B go away. As such, 
##We generate plot 2B twice - once with the breaks where it is saved off on its own
##And then a second time when it is stacked with FIgure 2A cartoon
##We then use powerpoint to replace FIgure 2B that does not include breaks with 
##The figure 2B that does include breaks


library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggbreak)
library(data.table)
library(export)
library(ggplotify)
# Folder containing your RDS files
folder <- "./RDS_Simulation_HistogramVarying/" #Comes from GenerateFig2Data.R

# List all RDS files in the folder
all_files <- list.files(folder, pattern = "\\.RDS$", full.names = TRUE)

all_files <- all_files[all_files %like% "MAE"]

# Keep only MAE RDS files
files <- list.files(
  "./RDS_Simulation_HistogramVarying",
  pattern = "MAE\\.RDS$",
  full.names = TRUE
)

df_all <- map_dfr(files, function(f) {
  
  fname <- basename(f)
  
  # Extract parameters
  matches <- str_match(
    fname,
    "alpha_([-0-9.]+)_omega_([-0-9.]+)_xi_([-0-9.]+)_MAE\\.RDS$"
  )
  
  alpha_val <- as.numeric(matches[,2])
  omega_val <- as.numeric(matches[,3])
  xi_val    <- as.numeric(matches[,4])
  
  # Read numeric vector
  mae_vec <- readRDS(f)
  
  # Ensure numeric
  mae_vec <- as.numeric(mae_vec)
  
  # Return long dataframe
  tibble(
    MAE = mae_vec,
    alpha = alpha_val,
    omega = omega_val,
    xi = xi_val
  )
})


# Convert to data.table
dt <- as.data.table(df_all)

# Compute summaries by alpha, omega, xi
dt_summary <- data.frame(dt[
  ,
  .(
    median_MAE = median(MAE, na.rm = TRUE),
    lower_2.5  = quantile(MAE, 0.025, na.rm = TRUE),
    upper_97.5 = quantile(MAE, 0.975, na.rm = TRUE)
  ),
  by = .(alpha, omega, xi)
])


p_MAE <- ggplot(dt_summary[dt_summary$xi %in% c(-1, 0, 1) &
                    dt_summary$alpha %in% c(-2, -1, 0, 1, 2),], 
       aes(x = factor(alpha), 
           y = median_MAE, 
           ymin = lower_2.5, 
           ymax = upper_97.5, 
           group = factor(omega), 
           colour = factor(omega))) +
  geom_point(size = 1.5, position = position_dodge(width = 0.3)) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.3)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Mean Average Error (Years)") + 
  xlab("Alpha") + 
  facet_wrap(~xi, nrow = 1, ncol = 5) +
  geom_hline(yintercept = 2.2, linetype = "dashed") +
  ggtitle("PiB") +
  scale_color_discrete(name = "Omega") +
  # This sets the lower and upper ranges explicitly
  scale_y_cut(breaks = c(4, 40), which = 3, scales = c(2, 0.5, 0.2))
p_MAE
graph2ppt(file = "./Figures/Simulation_MAE.pptx", width = 8, height = 6.5)


simulate_skewnorm_ggplot_2 <- function(
    xi_values = c(-1, 0, 1),
    omega_values = c(0.5, 1, 1.5, 2),
    alpha_values = c(-2, 0, 2),
    n = 5000
) {
  
  # --- Create parameter grid and simulate ---
  param_grid <- expand.grid(
    xi = xi_values,
    omega = omega_values,
    alpha = alpha_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  sim_data <- param_grid %>%
    rowwise() %>%
    mutate(values = list(rsn(n, xi = xi, omega = omega, alpha = alpha))) %>%
    unnest(cols = c(values)) %>%
    ungroup()
  
  # Factors for facetting + legend cleanliness
  sim_data <- sim_data %>%
    mutate(
      xi_f    = factor(xi),
      alpha_f = factor(alpha),
      omega_f = factor(omega)
    )
  
  # Final Plot:  
  # • xi → x-axis  
  # • alpha → facet rows  
  # • omega → color (overlay)  
  p <- ggplot(sim_data, aes(x = values, color = omega_f)) +
    geom_density(size = 1) +
    facet_grid(alpha_f ~ xi_f,
               labeller = labeller(
                 alpha_f = function(a) paste0("α = ", a),
                 xi_f    = function(x) paste0("ξ = ", x)
               )) +
    theme_bw(base_size = 14) +
    labs(
      x = "Skew-Normal Value",
      y = "Density",
      color = "ω (Scale)",
      title = "Skew-Normal Densities by α (skewness), ξ (location), ω (scale)"
    ) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      strip.text   = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    scale_x_continuous(limits = c(-5, 5))
  
  return(p)
}

p <- simulate_skewnorm_ggplot_2(
  xi_values = c(-1, 0, 1),
  omega_values = c(0.5, 1, 1.5, 2),
  alpha_values = c(-2, 0, 2)
)

grid.arrange(p, as.grob(p_MAE), nrow = 2, ncol = 1)

graph2ppt(file = "./Figures/Simulation_MAE_paired.pptx", width = 8, height = 8.5)

