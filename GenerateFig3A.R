library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggbreak)
library(patchwork)
dt_summary <- readRDS( "./Data/medianMAEspecificSimParamsnotest-rest.RDS") #BootstarppedMAEComparison_histogram
# alpha = -2 / omega = 2 / xi = -1
# alpha = -1 / omega = 1.5 / xi = 0
# alpha = 1 / omega = 1 / xi = 0.5
# alpha = 2 / omega = 0.5 / xi = 1
colnames(dt_summary) <- c("Alpha", "Omega", "Xi", "median", "CI_low", "CI_high")
# path to folder
dir_path <- "./Simulations20260424"

files <- list.files(
  path = dir_path,
  pattern = "^MAE_.*\\.RDS$",
  full.names = TRUE
)

results_df <- map_dfr(files, function(f) {
  
  fname <- basename(f)
  
  matches <- str_match(
    fname,
    "MAE_Xi(-?[0-9]+\\.?[0-9]*)Omega(-?[0-9]+\\.?[0-9]*)Alpha(-?[0-9]+\\.?[0-9]*)_mean_within_above_(-?[0-9]+\\.?[0-9]*)\\.RDS"
  )
  
  # sanity check (helps debugging if something still fails)
  if (any(is.na(matches))) {
    warning(paste("Failed to parse:", fname))
    return(NULL)
  }
  
  Xi     <- as.numeric(matches[2])
  Omega  <- as.numeric(matches[3])
  Alpha  <- as.numeric(matches[4])
  error  <- as.numeric(matches[5])
  
  varying <- readRDS(f)
  
  med <- median(varying, na.rm = TRUE)
  s   <- sd(varying, na.rm = TRUE)
  
  tibble(
    Xi = Xi,
    Omega = Omega,
    Alpha = Alpha,
    error = error,
    median = med,
    CI_low = med - 1.96 * s,
    CI_high = med + 1.96 * s
  )
})

results_df_updated <- results_df %>%
  left_join(
    dt_summary %>%
      select(Xi, Omega, Alpha,
             median_s = median,
             CI_low_s = CI_low,
             CI_high_s = CI_high),
    by = c("Xi", "Omega", "Alpha")
  ) %>%
  mutate(
    median  = ifelse(error == 0 & !is.na(median_s),  median_s,  median),
    CI_low  = ifelse(error == 0 & !is.na(CI_low_s),  CI_low_s,  CI_low),
    CI_high = ifelse(error == 0 & !is.na(CI_high_s), CI_high_s, CI_high)
  ) %>%
  select(-median_s, -CI_low_s, -CI_high_s)




get_plot <- function(Xi, Omega, Alpha, color1 = "grey70", color2 = "black"){
  p_main <- ggplot(results_df_updated[results_df_updated$Xi == Xi & results_df_updated$Omega == Omega & results_df_updated$Alpha == Alpha,], aes(x = error)) +
    
    geom_hline(yintercept = 2, linetype = "dashed") +
    
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                fill = color1, alpha = 0.4) +
    
    geom_line(aes(y = median),
              linewidth = 1.2, color = color2) +
    geom_vline(xintercept = 0.1, linetype = "dotted") +

    
    theme_classic() +
    theme(
      legend.position = "none"
    ) +
    labs(
      x = "Intra-Individual Error",
      y = "MAE (years)"
    ) + ylim(c(-30, 30))
  dt_summary$error <- 0
  p_inset <- ggplot(dt_summary[dt_summary$Xi == Xi & dt_summary$Omega == Omega & dt_summary$Alpha == Alpha,] , aes(x = error)) +
    
    geom_hline(yintercept = 2, linetype = "dashed") +
    
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                  fill = color1, alpha = 0.4) +
    
    geom_point(aes(y = median),
               linewidth = 1.1, color = color2) +
    
    theme_classic(base_size = 9) +
    theme(
      plot.background = element_rect(color = color2, linewidth = 0.5),
      axis.title = element_blank()
    )+ ylim(c(0, 6)) + scale_x_continuous(breaks = c(0))
  
  final_plot <- p_main +
    inset_element(
      p_inset,
      left = 0.78, bottom = 0.68,
      right = 0.98, top = 0.98
    )
  final_plot
  
  
  
}
p1 <- get_plot(Xi = 0, Omega = 1.5, Alpha = -1, color1 = "#00BFC4", color2 = "#007B82")
p2 <- get_plot(Xi = 0.5, Omega = 1, Alpha = 0, color1 = "#00BFC4", color2 = "#007B82")
p3 <- get_plot(Xi = 1, Omega = 0.5, Alpha = 2, color1 = "#00BFC4", color2 = "#007B82")

(p1 | p2 | p3)
graph2ppt(file = "./Figures/MAE_IntraindividualError3A.pptx", width =9, height = 4)
