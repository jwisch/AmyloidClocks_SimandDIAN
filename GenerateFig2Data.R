##This script lets you run the bootstrapping to generate the data that goes into the
#MAE and RMSE plots that compare the relative importance of the skew-normal parameters


xi_vec <- c(-1, -0.5, 0, 0.5, 1)
omega_vec <- c(0.5, 1, 1.5, 2)
alpha_vec <- c(-1, -0.5, 0.5, 1)


# --- libraries ---
library(mgcv)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(grid)
library(export)
source("./Simulation_functions.R")
source("./functions.R")
NBoots <- 1000 #Number of bootstraps

mean_within_above <- 0.07 #fixed within individual variance

n_ids <- 500
num_of_scans_for_model <- 2
minMeaningfulZ <- 0.9022
AposThresh <- 2.600
annual_rate <- 0.601
var_annual_rate <- 0.934
# mean_within_above <- 0.08 # PIB: 8% ± 7% from Tolboom 2009 paper
# # 10.3% ± 1.25% for pTau217 alzpath from brum 2023
# # 16.7% ± 1.8% for pTau181 gothenberg simoa from brum 2023
sd_within_above <- 0.03 * mean_within_above
mean_within_below <- mean_within_above / 2  # PIB: 4.4% ± 4.2% from Tolboom paper
sd_within_below <- 0.03 * mean_within_below

avg_interval <- 2.02
interval_noise <- 1.06




for(z in 1:length(xi_vec)){
  xi <- xi_vec[z]
  for(y in 1:length(omega_vec)){
    omega <- omega_vec[y]
    for(x in 1:length(alpha_vec)){
      alpha <- alpha_vec[x]

  seed <- 43266432

#Generating synthetic data
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
                                          prob_neg = 0.5,
                                          xi = xi,
                                          omega = omega,
                                          alpha = alpha,
                                          seed = seed)

df <- tmp[[1]]
converters <- tmp[[2]] #identify converters

df <- merge(df, converters, by = "ID", all.x = TRUE)
colnames(df)[3] <- "Z"
df <- df[, c("ID", "TimefromBaseline", "Z", "Z_true", "ever_cross", "first_cross")]
df[is.na(df$ever_cross),]$ever_cross <- FALSE


#Applying a manual shift to make it so not everyone converts within the first 5 years
df <- adjust_crossing_times(df,
                                mean_within_below = mean_within_below,
                                sd_within_below = sd_within_below)

##randomly sampling consecutive visits from the originally simulated data
#I will use the excluded data to get observed conversion dates outside the bounds
#of the data that I was training

MAE_list <- list()
R2_list <- list()
RMSE_list <- list()
df_result_list <- list()
df_result_true_list <- list()

df_interp_list <- list()
df <- df[df$TimefromBaseline < 16,]

for(j in 1:NBoots){
  print(j)
  set.seed((seed + j))
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

# df_result_true <- get_Time_to_Positivity(data.frame(df_sample[df_sample$Z_true > minMeaningfulZ,]), AposThresh,
#                                     id_name = "ID", time_name = "TimefromBaseline",
#                                     value_name = "Z_true")


colnames(df_result) <- c("ID", "Estimate", "Time_to_Positivity")
# colnames(df_result_true) <- c("ID", "Estimate_true", "Time_to_Positivity_true")


#This is the non-modeled data
df_excluded <- anti_join(df, df_sample)
df_excluded_baseline <- df_excluded[df_excluded$Z > minMeaningfulZ,]
df_excluded_baseline <- df_excluded_baseline[!duplicated(df_excluded_baseline$ID),]

df_excluded_baseline$timefromApos <- tryCatch({
  valid_idx <- which(!is.na(df_result$Time_to_Positivity) & !is.na(df_result$Estimate))
  if (length(valid_idx) < 2) {
    rep(NA, length(df_excluded_baseline$Z))
  } else {
    approx(
      x = df_result$Estimate[valid_idx],
      y = df_result$Time_to_Positivity[valid_idx],
      xout = df_excluded_baseline$Z
    )$y
  }
}, error = function(e) {
  rep(NA, length(df_excluded_baseline$Z))
})



# Function to linearly interpolate observed conversion from A- to A+
df_crossings <- df_excluded %>%
  group_by(ID) %>%
  # look at each row and the next row
  mutate(
    Z_next = lead(Z),
    Time_next = lead(TimefromBaseline),
    Z_next_true = lead(Z_true)
  ) %>%
  # keep only rows where a crossing occurs between current and next
  filter(!is.na(Z_next) & ((Z < AposThresh & Z_next >= AposThresh) | (Z > AposThresh & Z_next <= AposThresh))) %>%
  # compute linear interpolation
  mutate(
    Time_Crossing = TimefromBaseline +
      (AposThresh - Z) * (Time_next - TimefromBaseline) / (Z_next - Z)
  ) %>%
  mutate(
    Time_Crossing_true = TimefromBaseline +
      (AposThresh - Z_true) * (Time_next - TimefromBaseline) / (Z_next_true - Z_true)
  ) %>%
  select(ID, Time_Crossing, Time_Crossing_true) %>%
  ungroup()



df_excluded_baseline$Estimated_Time_Crossing <- df_excluded_baseline$TimefromBaseline - df_excluded_baseline$timefromApos
tmp <- df_excluded_baseline[!is.na(df_excluded_baseline$Estimated_Time_Crossing) & df_excluded_baseline$Estimated_Time_Crossing > 0,]

merged_df <- merge(df_crossings, tmp, by = "ID")

if(nrow(merged_df) > 10){
  
  df$TimeRound <- round(df$TimefromBaseline, 1)
  df_sample$TimeRound <- round(df_sample$TimefromBaseline, 1)
  
  
  setDT(df)[, in_sample := 0L]
  setDT(df)[setDT(df_sample), on = .(ID, TimeRound), in_sample := 1L]
  
  
  MAE_list[[j]] <- MLmetrics::MAE(merged_df$Time_Crossing, merged_df$Estimated_Time_Crossing)
  R2_list[[j]] <- MLmetrics::R2_Score(merged_df$Time_Crossing, merged_df$Estimated_Time_Crossing)
  RMSE_list[[j]] <- MLmetrics::RMSE(merged_df$Time_Crossing, merged_df$Estimated_Time_Crossing)
  
  df_result$run <- j
  
  df_result_list[[j]] <- df_result
  # Define target time grid
  time_grid <- seq(from = -5, to = 20, by = 0.5)
  
  # Interpolate within each run
  df_interp <- setDT(df_result)[, {
    # Ensure sorted order
    setorder(.SD, Time_to_Positivity)
    
    # Apply approx() only if enough points
    if (.N >= 2) {
      interp <- approx(x = Time_to_Positivity, y = Estimate, xout = time_grid, rule = 1, ties = "ordered")
      data.table(Time_to_Positivity = interp$x, Estimate = interp$y)
    } else {
      # Not enough points — return NA for all requested times
      data.table(Time_to_Positivity = time_grid, Estimate = NA_real_)
    }
  }, by = run]
  
  df_interp_list[[j]] <- df_interp
  
}else{
  j <- j - 1
  print("REDO")
}
MAE <- unlist(MAE_list)
R2 <- unlist(R2_list)
RMSE <- unlist(RMSE_list)
df_result <- rbindlist(df_result_list)
# df_result_true <- rbindlist(df_result_true_list)



saveRDS(df_result, paste0("./RDS_Simulation_HistogramVarying/PiB_alpha_", alpha, "_omega_", omega, "_xi_", xi, "_ResultingCurves.RDS"))
saveRDS(df_interp, paste0("./RDS_Simulation_HistogramVarying/PiB_alpha_", alpha, "_omega_", omega, "_xi_", xi, "_InterpolatedCurves.RDS"))

saveRDS(MAE, paste0("./RDS_Simulation_HistogramVarying/PiB_alpha_", alpha, "_omega_", omega, "_xi_", xi, "_MAE.RDS"))

saveRDS(RMSE, paste0("./RDS_Simulation_HistogramVarying/PiB_alpha_", alpha, "_omega_", omega, "_xi_", xi, "RMSE.RDS"))

saveRDS(R2, paste0("./RDS_Simulation_HistogramVarying/PiB_alpha_", alpha, "_omega_", omega, "_xi_", xi, "R2.RDS"))



}
    }
  }
}


    


