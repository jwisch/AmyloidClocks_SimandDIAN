library(Budgeon)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cutpointr)
library(mclust)
library(lme4)
library(MuMIn)
library(rlang)
library(forestplot)
library(tableone)
library(export)
library(data.table)
library(GGally)
library(ggbreak)

source("./functions.R")
source("./functions_DistributionSpecificity.R")

source("./CleaningandPrepofRealData.R")



# #Thisi s how to filter out the subthreshold values.
# # bootstrap_matrix$interpolated_val[bootstrap_matrix$interpolated_val < rel_accum] <- NA
# 
# # result_PET_CSFpT217_Nico_Z <- bs_DistSpecplot(data.frame(df), cp_PET_CSFpT217_Nico_obj[[3]]$optimal_cutpoint,
# #                                                              "newid18", "TimefromBaseline",
# #                                                              "CSFpT217_Nico_Z", num_bootstraps = 1000)
# # result_PET_plasmapTau217_Nico_Z <- bs_DistSpecplot(data.frame(df), cp_PET_plasmapTau217_Nico_obj[[3]]$optimal_cutpoint,
# #                                                                   "newid18", "TimefromBaseline",
# #                                                                   "plasmapTau217_Nico_Z", num_bootstraps = 1000)
# # result_PET_Alamar_pTau217_Z <- bs_DistSpecplot(data.frame(df), cp_PET_Alamar_pTau217_obj[[3]]$optimal_cutpoint,
# #                                                               "newid18", "TimefromBaseline",
# #                                                               "Alamar_pTau217_Z", num_bootstraps = 1000)
# result_PET_AlamarCSF_pTau217_Z <- robust_bs_DistSpecplot(data.frame(df), cp_PET_AlamarCSF_pTau217_obj[[3]]$optimal_cutpoint,
#                                                                 "newid18", "TimefromBaseline",
#                                                                 "AlamarCSF_pTau217_Z", num_bootstraps = 1000)
# # result_pib_Z <- bs_DistSpecplot(data.frame(df), (18 - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE),
# #                                                "newid18", "TimefromBaseline",
# #                                                "CL_Z", num_bootstraps = 1000)
# # # saveRDS(result_PET_Alamar_pTau217_Z, "./DistSpecBootstraps/result_PET_Alamar_pTau217_Z_batchCorrected_20260420.RDS")
# # saveRDS(result_pib_Z, "./DistSpecBootstraps/result_pib_Z_20260420.RDS")
# 
# # saveRDS(result_PET_CSFpT217_Nico_Z, "./DistSpecBootstraps/result_PET_CSFpT217_Nico_Z_20260420.RDS")
# # saveRDS(result_PET_plasmapTau217_Nico_Z, "./DistSpecBootstraps/result_PET_plasmapTau217_Nico_Z_20260420.RDS")
# saveRDS(result_PET_AlamarCSF_pTau217_Z, "./DistSpecBootstraps/result_PET_AlamarCSF_pTau217_Z_batchCorrected_20260420.RDS")
# 
# # 
# # 
# # result_PET_CSFpT181_Nico_Z <- bs_DistSpecplot(data.frame(df), cp_PET_CSFpT181_Nico_obj[[3]]$optimal_cutpoint,
# #                                                              "newid18", "TimefromBaseline",
# #                                                              "CSFpT181_Nico_Z", num_bootstraps = 1000)
# # result_PET_plasmapTau181_Nico_Z <- bs_DistSpecplot(data.frame(df), cp_PET_plasmapTau181_Nico_obj[[3]]$optimal_cutpoint,
# #                                                                   "newid18", "TimefromBaseline",
# #                                                                   "plasmapTau181_Nico_Z", num_bootstraps = 1000)
# # result_PET_plasmapTau181_Jucker_Z <- bs_DistSpecplot(data.frame(df), cp_PET_plasmapTau181_Jucker_obj[[3]]$optimal_cutpoint,
# #                                                                     "newid18", "TimefromBaseline",
# #                                                                     "plasmapTau181_jucker_Z", num_bootstraps = 1000)
# # result_PET_Alamar_pTau181_Z <- bs_DistSpecplot(data.frame(df), cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint,
# #                                                             "newid18", "TimefromBaseline",
# #                                                               "Alamar_pTau181_Z", num_bootstraps = 1000)
# result_PET_AlamarCSF_pTau181_Z <- robust_bs_DistSpecplot(data.frame(df), cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint,
#                                                             "newid18", "TimefromBaseline",
#                                                               "AlamarCSF_pTau181_Z", num_bootstraps = 1000)
# # result_PET_LumipulseCSF_pTau181_Z <- bs_DistSpecplot(data.frame(df), cp_PET_LumipulseCSF_pTau181_obj[[3]]$optimal_cutpoint,
# #                                                             "newid18", "TimefromBaseline",
# #                                                               "LUMIPULSE_CSF_pTau_Z", num_bootstraps = 1000)
# # 
# # saveRDS(result_PET_CSFpT181_Nico_Z, "./DistSpecBootstraps/result_PET_CSFpT181_Nico_Z_20260420.RDS")
# # saveRDS(result_PET_plasmapTau181_Nico_Z, "./DistSpecBootstraps/result_PET_plasmapTau181_Nico_Z_20260420.RDS")
# # saveRDS(result_PET_plasmapTau181_Jucker_Z, "./DistSpecBootstraps/result_PET_plasmapTau181_Jucker_Z_20260420.RDS")
# # saveRDS(result_PET_Alamar_pTau181_Z, "./DistSpecBootstraps/result_PET_Alamar_pTau181_Z_batchCorrected_20260420.RDS")
# saveRDS(result_PET_AlamarCSF_pTau181_Z, "./DistSpecBootstraps/result_PET_AlamarCSF_pTau181_Z_batchCorrected_20260420.RDS")
# saveRDS( result_PET_LumipulseCSF_pTau181_Z, "./DistSpecBootstraps/result_PET_LumipulseCSF_pTau181_Z_batchCorrected_20260420.RDS")




# result_PET_BD.pTau.217_Z <- bs_DistSpecplot(data.frame(df), cp_PET_BD.pTau.217_obj[[3]]$optimal_cutpoint,
#                                                              "newid18", "TimefromBaseline",
#                                                              "BD.pTau.217_Z", num_bootstraps = 1000)
# saveRDS( result_PET_BD.pTau.217_Z, "./DistSpecBootstrapsresult_PET_BD.pTau.217_Z_20260420.RDS")
# 






result_PET_Alamar_pTau217_Z <- readRDS("./DistSpecBootstraps/result_PET_Alamar_pTau217_Z_batchCorrected_20260420.RDS")
result_pib_Z <- readRDS("./DistSpecBootstraps/result_pib_Z_20260420.RDS")
result_PET_CSFpT217_Nico_Z <-  readRDS("./DistSpecBootstraps/result_PET_CSFpT217_Nico_Z_20260420.RDS")
result_PET_plasmapTau217_Nico_Z <- readRDS("./DistSpecBootstraps/result_PET_plasmapTau217_Nico_Z_20260420.RDS")
result_PET_AlamarCSF_pTau217_Z <- readRDS("./DistSpecBootstraps/result_PET_AlamarCSF_pTau217_Z_batchCorrected_20260420.RDS")
result_PET_CSFpT181_Nico_Z <- readRDS("./DistSpecBootstraps/result_PET_CSFpT181_Nico_Z_20260420.RDS")
result_PET_plasmapTau181_Nico_Z <- readRDS("./DistSpecBootstraps/result_PET_plasmapTau181_Nico_Z_20260420.RDS")
result_PET_plasmapTau181_Jucker_Z <- readRDS("./DistSpecBootstraps/result_PET_plasmapTau181_Jucker_Z_20260420.RDS")
result_PET_Alamar_pTau181_Z <- readRDS("./DistSpecBootstraps/result_PET_Alamar_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_AlamarCSF_pTau181_Z <- readRDS("./DistSpecBootstraps/result_PET_AlamarCSF_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_LumipulseCSF_pTau181_Z <- readRDS("./DistSpecBootstraps/result_PET_LumipulseCSF_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_BD.pTau.217_Z <-  readRDS("./DistSpecBootstrapsresult_PET_BD.pTau.217_Z_20260420.RDS")




p1 <- plot_hist_density_at_time(data = result_pib_Z, rel_accum = rel_accum_pib_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c("#003B44",   # darker lead-in
                                                "#007B82",   # 2nd (your first required color)
                                                "#00BFC4",   # 3rd (your second required color)
                                                "#6FE3E6",   # lighter continuation
                                                "#CFF7F8") ,  # very light endpoint
                                "F1. Amyloid PET",
                                scale_break = TRUE,
                                scale_break1 = 7,
                                scale_break2 = 45,
                                scale_break_scale = 0.2)

p2 <- plot_hist_density_at_time(data = result_PET_BD.pTau.217_Z, rel_accum = rel_accum_BD.pTau.217_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "F2. CSF BD pT181 [Alamar]")

p3 <- plot_hist_density_at_time(data = result_PET_CSFpT181_Nico_Z, rel_accum = rel_accum_CSFpT181_Nico_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "F3. CSF pT181 [IP-MS]")

p4 <- plot_hist_density_at_time(data = result_PET_CSFpT217_Nico_Z, rel_accum = rel_accum_CSFpT217_Nico_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "F4. CSF pT217 [IP-MS]")

p5 <- plot_hist_density_at_time(data = result_PET_AlamarCSF_pTau217_Z, rel_accum = rel_accum_AlamarCSFpT217_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "F5. CSF pT217 [Alamar]")

p6 <- plot_hist_density_at_time(data = result_PET_plasmapTau217_Nico_Z, rel_accum = rel_accum_plasmapT217_Nico_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#A93E36",  # darkest
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ,  # very light endpoint
                                "F6. Plasma pT217 [IP-MS]")

p7 <- plot_hist_density_at_time(data = result_PET_LumipulseCSF_pTau181_Z, rel_accum = rel_accum_LumipulseCSFpT181_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "F7. CSF pT181 [Lumipulse]")

p8 <- plot_hist_density_at_time(data = result_PET_AlamarCSF_pTau181_Z, rel_accum = rel_accum_AlamarCSFpT181_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "grey10",  # darkest
                                  "grey30",  # 2nd darkest (your anchor)
                                  "grey50",
                                  "grey70",
                                  "grey90"   # lightest
                                ),  # very light endpoint
                                "F8. CSF pT181 [Alamar]")

p9 <- plot_hist_density_at_time(data = result_PET_plasmapTau181_Jucker_Z, rel_accum = rel_accum_plasmapT181_Jucker_Z, 
                                time_vals = c(-2 ,0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c(
                                  "#A93E36",  # darkest
                                  "#F8766D",  # 2nd darkest (your anchor)
                                  "#FB9A93",
                                  "#FDC1BD",
                                  "#FFE5E3"   # lightest
                                ) ,  # very light endpoint
                                "F9. Plasma pT181 [Simoa]")

p10 <- plot_hist_density_at_time(data = result_PET_Alamar_pTau181_Z, rel_accum = rel_accum_AlamarpT181_Z, 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "#A93E36",  # darkest
                                   "#F8766D",  # 2nd darkest (your anchor)
                                   "#FB9A93",
                                   "#FDC1BD",
                                   "#FFE5E3"   # lightest
                                 ) ,  # very light endpoint
                                 "F10. Plasma pT181 [Alamar]")

p11 <- plot_hist_density_at_time(data = result_PET_plasmapTau181_Nico_Z, rel_accum = rel_accum_plasmapT181_Nico_Z, 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "#A93E36",  # darkest
                                   "#F8766D",  # 2nd darkest (your anchor)
                                   "#FB9A93",
                                   "#FDC1BD",
                                   "#FFE5E3"   # lightest
                                 ) ,  # very light endpoint
                                 "F11. Plasma pT181 [IP-MS]")

p12 <- plot_hist_density_at_time(data = result_PET_Alamar_pTau217_Z, rel_accum = rel_accum_AlamarpT217_Z, 
                                 time_vals = c(-2 ,0, 2, 4, 6),
                                 binwidth = 0.01,
                                 color_vec =   c(
                                   "#A93E36",  # darkest
                                   "#F8766D",  # 2nd darkest (your anchor)
                                   "#FB9A93",
                                   "#FDC1BD",
                                   "#FFE5E3"   # lightest
                                 ) ,  # very light endpoint
                                 "F12. Plasma pT217 [Alamar]")


lemon::grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4, ncol = 3)
graph2ppt(file = "./Figures/BrianPlot_all12.pptx", width = 4, height = 7)



# Convert to data.table for convenience
tmp <- result_PET_plasmapTau181_Jucker_Z[result_PET_plasmapTau181_Jucker_Z$Time_Window %in% c(4, 6),]
tmp <- tmp[!is.na(tmp$interpolated_val),]
dt <- as.data.table(tmp)

# List of unique time windows
time_windows <- sort(unique(dt$Time_Window))
dt <- dt[!is.na(dt),]
# Compute density estimates per Time_Window
dens_list <- lapply(time_windows, function(tw) {
  d <- density(dt[Time_Window == tw, interpolated_val])
  list(time_window = tw, x = d$x, y = d$y)
})

x0 <- 2.7

prob_y <- sapply(dens_list, function(d) {
  # Find density at x0 using linear interpolation
  approx(d$x, d$y, xout = x0)$y
})

names(prob_y) <- sapply(dens_list, function(d) d$time_window)
prob_y

probabilities <- prob_y / sum(prob_y, na.rm = TRUE)
probabilities
