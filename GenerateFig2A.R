
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggbreak)
library(data.table)
library(export)
library(tidyr)
library(ggplotify)


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

p
graph2ppt(file = "./Figures/Simulation_MAE_paired.pptx", width = 8, height = 8.5)



