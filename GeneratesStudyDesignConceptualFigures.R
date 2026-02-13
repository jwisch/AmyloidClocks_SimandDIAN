###############################################################################
##The purpose of this script is to generate figures that were used for 
##conceptual figure 1. This script generates synthetic data that was not
##used in any analysis. These figures were extracted and arranged in a powerpoint
##slide. Labelling was added post hoc. 
###############################################################################


# --- libraries ---
library(mgcv)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(grid)
library(export)
library(Budgeon)
library(gridExtra)
source("./Simulation_functions.R")
source("./functions.R")

source("./params_pib.R")
n_ids <- 300
num_of_scans_for_model <- 2

tmp <- generate_synthetic_data_with_noise(n_ids = n_ids, 
                                          n_visits = 8,
                                          avg_interval = avg_interval,
                                          interval_noise = interval_noise,
                                          minMeaningfulZ = minMeaningfulZ,
                                          AposThresh = AposThresh,
                                          maxMeaningfulZ = 25,
                                          annual_rate = annual_rate,
                                          var_annual_rate = var_annual_rate,
                                          mean_within_above = mean_within_above,
                                          sd_within_above   = sd_within_above,
                                          mean_within_below = mean_within_below,
                                          sd_within_below   = sd_within_below,
                                          plateau_decay = 0.15,
                                          prob_neg = 0.19,
                                          xi = xi,
                                          omega = omega,
                                          alpha = alpha)

df <- tmp[[1]]
converters <- tmp[[2]] #identify converters

df <- merge(df, converters, by = "ID", all.x = TRUE)
df <- df[, c("ID", "TimefromBaseline", "Z_noisy", "Z_true", "ever_cross", "first_cross")]
df[is.na(df$ever_cross),]$ever_cross <- FALSE
colnames(df)[3] <- "Z"

#Applying a manual shift to make it so not everyone converts within the first 5 years
df <- adjust_crossing_times(df,
                            mean_within_below = mean_within_below,
                            sd_within_below = sd_within_below)
df <- df[df$TimefromBaseline >= 0 & df$TimefromBaseline <= 16,]



##randomly sampling consecutive visits from the originally simulated data
#I will use the excluded data to get observed conversion dates outside the bounds
#of the data that I was training
df_sample <- df %>%
  group_by(ID) %>%
  arrange(TimefromBaseline, .by_group = TRUE) %>%  # ensures chronological order
  group_modify(~ {
    n <- nrow(.x)
    if (n < num_of_scans_for_model) {
      return(.x)  # not enough rows, return all
    } else {
      start_idx <- sample(1:(n - (num_of_scans_for_model - 1)), 1)
      return(.x[start_idx:(start_idx + (num_of_scans_for_model - 1)), ])
    }
  }) %>%
  ungroup()




df_result <- get_Time_to_Positivity(data.frame(df_sample[df_sample$Z > minMeaningfulZ,]), AposThresh,
                                    id_name = "ID", time_name = "TimefromBaseline",
                                    value_name = "Z")

df_result_true <- get_Time_to_Positivity(data.frame(df_sample[df_sample$Z_true > minMeaningfulZ,]), AposThresh,
                                         id_name = "ID", time_name = "TimefromBaseline",
                                         value_name = "Z_true")


colnames(df_result) <- c("ID", "Estimate", "Time_to_Positivity")
colnames(df_result_true) <- c("ID", "Estimate_true", "Time_to_Positivity_true")


#This is the non-modeled data
df_excluded <- anti_join(df, df_sample)
df_excluded_baseline <- df_excluded[df_excluded$Z > minMeaningfulZ,]
df_excluded_baseline <- df_excluded_baseline[!duplicated(df_excluded_baseline$ID),]

df_excluded_baseline$timefromApos <- approx(y = df_result$Time_to_Positivity, 
                                            x = df_result$Estimate, xout = df_excluded_baseline$Z)$y


ggplot(df, aes(x = TimefromBaseline, y = Z, group = ID)) + geom_line(colour = "black") + theme_bw() +
 # geom_line(data = df_sample, aes(x = TimefromBaseline, y = Z, group = ID), colour = "black") +
  xlab("Time from Baseline")

ggplot(df, aes(x = TimefromBaseline, y = Z, group = ID)) + geom_line(colour = "grey80") + theme_bw() +
  geom_line(data = df_sample, aes(x = TimefromBaseline, y = Z, group = ID), colour = "black") +
  xlab("Time from Baseline")


ggplot(df_excluded, aes(x = TimefromBaseline, y = Z, group = ID)) + geom_line(colour = "black") + theme_bw() +
  geom_line(data = df_sample, aes(x = TimefromBaseline, y = Z, group = ID), colour = "white", size = 2) +
  xlab("Time from Baseline")



individual_coeffs <- get_individ_lm(df_sample, "Z", "TimefromBaseline", "ID")

individual_coeffs <- individual_coeffs[individual_coeffs$Slope >= 0,]

eval_at_midpoint = eval_lmm(df_sample[df_sample[["ID"]] %in% individual_coeffs[["ID"]],], id_name = "ID", time_name = "TimefromBaseline", individual_coeffs) 

df_tau <- get_integration_estimates(eval_at_midpoint, id_name = "ID", degree = 3)

hor_adj <- get_horizontal_adjustment(df_tau, 1.62)

df_tau$Time_to_Positivity <- df_tau$tau - hor_adj


plot1_df <- merge(df_sample[!duplicated(df_sample$ID),], individual_coeffs, by = "ID", all = FALSE)

p1 <- ggplot(plot1_df, aes(x = Z, y = Slope)) + geom_point() + theme_bw() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), colour = "#007B82") + xlab("Baseline Z") + ylab("Annualized Rate of Change") + ylim(c(-0.1, 3.6))

p2 <- ggplot(df_tau, aes(x = mu_tau,  y = tau)) + geom_point() + theme_bw() +
  geom_smooth(method = "gam", se = FALSE, colour = "#007B82") +
 xlab("Baseline Z") + ylab("Integral of (1/ARC)dARC")

p3 <- ggplot(df_tau, aes(y = mu_tau, x = tau)) + geom_point() + theme_bw() +
  geom_smooth(method = "gam", se = FALSE, colour = "#007B82") + 
  xlab("Disease Progression Time") + ylab("Baseline Z")

grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
