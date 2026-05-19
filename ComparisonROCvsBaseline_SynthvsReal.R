library(mgcv)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(grid)
library(export)
library(boot)
source("./Simulation_functions.R")
source("./functions.R")

spearman_boot <- function(data, indices) {
  
  d <- data[indices, ]
  
  cor(
    d$CL,
    d$rate_of_change,
    method = "spearman",
    use = "complete.obs"
  )
}


ADNI <- readRDS(".././Alamar_AmyloidClocks/Data/ADNI_matched.RDS")#not available to post on git because of data sharing restrictions
DIAN <- read.csv("../Alamar_AmyloidClocks/Data/cleaned_df_withBatchNorming_20260420.csv") #not available to post on git because of data sharing restrictions

#Generating synthetic options

Time_Window <- seq(from = -5, to = 20, by = 0.5)

# --- libraries ---

NBoots <- 1000 #Number of bootstraps

mean_within_above <- 0 #fixed within individual variance

n_ids <- 300
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

xi_vec <- c(-1, 0, 1)
omega_vec <- c(0.5, 1, 1.5, 2)
alpha_vec <- c(-2, 0, 2)
combos <- data.frame(expand.grid("xi" = xi_vec,
                     "omega" = omega_vec,
                     "alpha" = alpha_vec))



# --- Pre-allocate storage ---
# We expect nrow(combos) * length(omega_vec) * length(alpha_vec) * NBoots rows
total_combos <- nrow(combos) * length(omega_vec) * length(alpha_vec)
results_list <- vector("list", total_combos)
counter <- 1

for(z in 1:nrow(combos)){
  xi <- combos$xi[z]
  for(y in 1:length(omega_vec)){
    omega <- omega_vec[y]
    for(x in 1:length(alpha_vec)){
      alpha <- alpha_vec[x]
      
      seed <- 43266432
      
      # Generating synthetic data
      tmp <- generate_synthetic_data_with_noise(
        n_ids = n_ids, 
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
        seed = seed
      )
      
      df <- tmp[[1]]
      converters <- tmp[[2]] 
      
      df <- merge(df, converters, by = "ID", all.x = TRUE)
      colnames(df)[3] <- "Z"
      df <- df[, c("ID", "TimefromBaseline", "Z", "Z_true", "ever_cross", "first_cross")]
      df[is.na(df$ever_cross),]$ever_cross <- FALSE
      
      # Applying manual shift
      # (Assuming adjust_crossing_times is defined in your environment)
      df <- adjust_crossing_times(df,
                                  mean_within_below = mean_within_below,
                                  sd_within_below = sd_within_below)
      
      # Vector to store 1000 rhos for this specific xi/omega/alpha combo
      boot_rhos <- numeric(NBoots)
      
      for(j in 1:NBoots){
        set.seed((seed + j))
        
        # 1. Generate the sample (2 visits per ID)
        df_sample <- df %>%
          group_by(ID) %>%
          arrange(TimefromBaseline, .by_group = TRUE) %>%
          group_modify(~ {
            n <- nrow(.x)
            if (n < num_of_scans_for_model) {
              return(.x)
            } else {
              start_idx <- sample(1:(n - (num_of_scans_for_model - 1)), 1)
              return(.x[start_idx:(start_idx + (num_of_scans_for_model - 1)), ])
            }
          }) %>%
          ungroup()
        
        # 2. Extract baseline and annualized rate of change
        # Rate = (Z_final - Z_initial) / (Time_final - Time_initial)
        sample_stats <- df_sample %>%
          group_by(ID) %>%
          summarize(
            baseline_val = first(Z),
            annualized_rate = (last(Z) - first(Z)) / (last(TimefromBaseline) - first(TimefromBaseline)),
            .groups = "drop"
          )
        
        # 3. Calculate Spearman correlation and store
        boot_rhos[j] <- cor(sample_stats$baseline_val, 
                            sample_stats$annualized_rate, 
                            method = "spearman", 
                            use = "complete.obs")
      }
      
      # Store the full distribution for this combo
      results_list[[counter]] <- data.frame(
        xi = xi,
        omega = omega,
        alpha = alpha,
        boot_id = 1:NBoots,
        spearman_rho = boot_rhos
      )
      
      message(sprintf("Completed Combo %d/%d: xi=%0.1f, om=%0.1f, al=%0.1f", 
                      counter, total_combos, xi, omega, alpha))
      counter <- counter + 1
    }
  }
}

# --- Final Step: Combine all results into one master data frame ---
final_correlations <- bind_rows(results_list)

saveRDS(final_correlations, "./Data/ROCtoBaselineCorrelations_for_Simulations")

# You can now examine the distribution, for example:
ggplot(final_correlations, aes(x = spearman_rho, fill = factor(xi))) +
  geom_density(alpha = 0.5) +
  facet_grid(omega ~ alpha) +
  theme_minimal()



DIAN <- DIAN[!is.na(DIAN$VISITAGEc),]
DIAN <- DIAN %>%
  group_by(newid18) %>%
  mutate(TimefromBaseline = VISITAGEc - min(VISITAGEc, na.rm = TRUE)) %>%
  ungroup()

library(boot)

get_spearman_boot_ci <- function(data,
                                 id_col,
                                 biomarker_col,
                                 time_col,
                                 R = 1000,
                                 ci_type = "bca") {
  
  # Generate rate-of-change dataframe
  roc_df <- get_RofC_df(
    data[!is.na(data[[biomarker_col]]), ],
    id_col,
    biomarker_col,
    time_col
  )
  
  roc_df <- as.data.frame(roc_df)
  
  # Baseline rows
  baseline_df <- data[
    data[[time_col]] == 0,
    c(id_col, biomarker_col, time_col),
    drop = FALSE
  ]
  
  baseline_df <- as.data.frame(baseline_df)
  
  # Remove duplicated column names if present
  names(roc_df) <- make.unique(names(roc_df))
  names(baseline_df) <- make.unique(names(baseline_df))
  
  # Verify merge column exists
  if(!(id_col %in% names(roc_df))) {
    stop(paste(id_col, "not found in roc_df"))
  }
  
  if(!(id_col %in% names(baseline_df))) {
    stop(paste(id_col, "not found in baseline_df"))
  }
  
  # Merge
  roc_df <- merge(
    roc_df,
    baseline_df,
    by = id_col,
    all = FALSE
  )
  
  # Bootstrap function
  spearman_boot <- function(data, indices) {
    
    d <- data[indices, ]
    
    cor(
      d[[biomarker_col]],
      d$rate_of_change,
      method = "spearman",
      use = "complete.obs"
    )
  }
  
  # Run bootstrap
  boot_obj <- boot(
    data = roc_df,
    statistic = spearman_boot,
    R = R
  )
  
  # Confidence interval
  ci <- boot.ci(boot_obj, type = ci_type)
  
  if(ci_type == "bca") {
    lower <- ci$bca[4]
    upper <- ci$bca[5]
  } else if(ci_type == "perc") {
    lower <- ci$percent[4]
    upper <- ci$percent[5]
  } else if(ci_type == "basic") {
    lower <- ci$basic[4]
    upper <- ci$basic[5]
  } else if(ci_type == "norm") {
    lower <- ci$normal[2]
    upper <- ci$normal[3]
  } else {
    stop("Unsupported ci_type")
  }
  
  list(
    estimate = as.numeric(boot_obj$t0),
    lower_ci = as.numeric(lower),
    upper_ci = as.numeric(upper),
    boot_object = boot_obj
  )
}


result_CL <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "CL",
  time_col = "TimefromBaseline"
)

result_LUMIPULSE_CSF_pTau <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "LUMIPULSE_CSF_pTau",
  time_col = "TimefromBaseline"
)

result_plasmapTau217_Nico <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "plasmapTau217_Nico",
  time_col = "TimefromBaseline"
)

result_plasmapTau181_Nico <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "plasmapTau181_Nico",
  time_col = "TimefromBaseline"
)


result_plasmapTau181_jucker <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "plasmapTau181_jucker",
  time_col = "TimefromBaseline"
)

result_CSFpT181_Nico_corrected <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "CSFpT181_Nico_corrected",
  time_col = "TimefromBaseline"
)

result_CSFpT217_Nico_corrected <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "CSFpT217_Nico_corrected",
  time_col = "TimefromBaseline"
)

result_Alamar_pTau181 <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "Alamar_pTau181",
  time_col = "TimefromBaseline"
)

result_Alamar_pTau217 <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "Alamar_pTau217",
  time_col = "TimefromBaseline"
)

result_AlamarCSF_pTau217 <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "AlamarCSF_pTau217",
  time_col = "TimefromBaseline"
)

result_AlamarCSF_pTau181 <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "AlamarCSF_pTau181",
  time_col = "TimefromBaseline"
)

result_BD.pTau.217 <- get_spearman_boot_ci(
  data = DIAN,
  id_col = "newid18",
  biomarker_col = "BD.pTau.217",
  time_col = "TimefromBaseline"
)

result_ADNI_CL <- get_spearman_boot_ci(
  data = data.frame(ADNI),
  id_col = "PTID",
  biomarker_col = "CENTILOIDS.combat",
  time_col = "TimefromBaseline_PET"
)

result_ADNI_plasma <- get_spearman_boot_ci(
  data = data.frame(ADNI),
  id_col = "PTID",
  biomarker_col = "pT217_AB42_F",
  time_col = "TimefromBaseline_plasma"
)

realdata_results <- rbind(data.frame("name" = "CL_ADNI", "estimate" = result_ADNI_CL$estimate, "upper_ci" = result_ADNI_CL$upper_ci, 
                                     "lower_ci" = result_ADNI_CL$lower_ci),
                          data.frame("name" = "plasmapTau217_ADNI", "estimate" = result_ADNI_plasma$estimate, "upper_ci" = result_ADNI_plasma$upper_ci, 
                                     "lower_ci" = result_ADNI_plasma$lower_ci),
                          data.frame("name" = "Alamar_pTau181", "estimate" = result_Alamar_pTau181$estimate, "upper_ci" = result_Alamar_pTau181$upper_ci, 
                                     "lower_ci" = result_Alamar_pTau181$lower_ci),
                          data.frame("name" = "Alamar_pTau217", "estimate" = result_Alamar_pTau217$estimate, "upper_ci" = result_Alamar_pTau217$upper_ci, 
                                     "lower_ci" = result_Alamar_pTau217$lower_ci),
                          data.frame("name" = "AlamarCSF_pTau181", "estimate" = result_AlamarCSF_pTau181$estimate, "upper_ci" = result_AlamarCSF_pTau181$upper_ci, 
                                     "lower_ci" = result_AlamarCSF_pTau181$lower_ci),
                          data.frame("name" = "BD.pTau.217", "estimate" = result_BD.pTau.217$estimate, "upper_ci" = result_BD.pTau.217$upper_ci, 
                                     "lower_ci" = result_BD.pTau.217$lower_ci),
                          data.frame("name" = "CL_DIAN", "estimate" = result_CL$estimate, "upper_ci" = result_CL$upper_ci, "lower_ci" = result_CL$lower_ci),
                          data.frame("name" = "CSFpT181_Nico_corrected", "estimate" = result_CSFpT181_Nico_corrected$estimate, 
                                     "upper_ci" = result_CSFpT181_Nico_corrected$upper_ci, 
                                     "lower_ci" = result_CSFpT181_Nico_corrected$lower_ci),
                          data.frame("name" = "CSFpT217_Nico_corrected", "estimate" = result_CSFpT217_Nico_corrected$estimate, 
                                     "upper_ci" = result_CSFpT217_Nico_corrected$upper_ci, 
                                     "lower_ci" = result_CSFpT217_Nico_corrected$lower_ci),
                          data.frame("name" = "LUMIPULSE_CSF_pTau", "estimate" = result_LUMIPULSE_CSF_pTau$estimate, 
                                     "upper_ci" = result_LUMIPULSE_CSF_pTau$upper_ci, 
                                     "lower_ci" = result_LUMIPULSE_CSF_pTau$lower_ci),
                          data.frame("name" = "plasmapTau181_jucker", "estimate" = result_plasmapTau181_jucker$estimate, 
                                     "upper_ci" = result_plasmapTau181_jucker$upper_ci, 
                                     "lower_ci" = result_plasmapTau181_jucker$lower_ci),
                          data.frame("name" = "plasmapTau181_Nico", "estimate" = result_plasmapTau181_Nico$estimate, 
                                     "upper_ci" = result_plasmapTau181_Nico$upper_ci, 
                                     "lower_ci" = result_plasmapTau181_Nico$lower_ci),
                          data.frame("name" = "plasmapTau217_Nico", "estimate" = result_plasmapTau217_Nico$estimate, 
                                     "upper_ci" = result_plasmapTau217_Nico$upper_ci, 
                                     "lower_ci" = result_plasmapTau217_Nico$lower_ci))

syntheticdata_results <- setDT(final_correlations)[, .(mean(spearman_rho), sd(spearman_rho)), by = c("xi", "omega", "alpha")]
colnames(syntheticdata_results)[4:5] <- c("estimate", "SD")
syntheticdata_results$lower_ci <- syntheticdata_results$estimate - 1.96 * syntheticdata_results$SD
syntheticdata_results$upper_ci <- syntheticdata_results$estimate + 1.96 * syntheticdata_results$SD

names <- rownames(realdata_results)
realdata_results <- data.frame(realdata_results)
realdata_results$names <- names
syntheticdata_results$name <- paste0("xi", syntheticdata_results$xi, "omega",
                                      syntheticdata_results$omega, "alpha",
                                      syntheticdata_results$alpha)
realdata_results$group <- "researchData"
syntheticdata_results$group <- "syntheticData"

forplotting <- rbind(data.frame(realdata_results[, c("name", "estimate", "lower_ci", "upper_ci", "group")]),
                     data.frame(syntheticdata_results[, c("name", "estimate", "lower_ci", "upper_ci", "group")]))
saveRDS(forplotting, "./Data/ROCtoBaselineCorrelations_for_All.RDS")


ggplot(forplotting, aes(x = reorder(name, estimate), y = estimate, ymin = lower_ci, ymax = upper_ci, colour = group)) +
  geom_point() + geom_errorbar() + theme_bw() + ylab("Spearman Correlation with 95% bootstrapped confidence intervals") +
  xlab("Biomarker") + theme(legend.position = "bottom") + coord_flip() + scale_colour_manual(values = c("#E07529", "#7F7991"))
graph2ppt(file = "./Figures/SuppFig_CorrelationsBaselineandARC.pptx", width = 6, height = 8)
