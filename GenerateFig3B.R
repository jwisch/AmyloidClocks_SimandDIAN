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
library(tidyr)
library(lme4)
library(lmerTest)
library(sn)
library(plotly)
library(forcats)
library(emmeans) # for confidence intervals
library(export)

source("./functions.R")
source("./Simulation_functions.R")
source("./Fig1and2Funcs.R")
source("./ADNI_prep_for_analysis.R")

source("./CleaningandPrepofRealData.R")


####### READ IN RESULTS FROM GENERATE_TIMETOPOSITIVITY_REALDATA
result_PET_CSFpT217_Nico_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_CSFpT217_Nico_Z_20260420.RDS")
result_PET_Alamar_pTau217_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_Alamar_pTau217_Z_batchCorrected_20260420.RDS")
result_pib_Z <- readRDS( "../Alamar_AmyloidClocks/BootstrapResults/result_pib_Z_20260420.RDS")
result_PET_AlamarCSF_pTau217_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_AlamarCSF_pTau217_Z_batchCorrected_20260420.RDS")
result_PET_plasmapTau217_Nico_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_plasmapTau217_Nico_Z_20260420.RDS")

result_PET_CSFpT181_Nico_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_CSFpT181_Nico_Z_20260420.RDS")
result_PET_plasmapTau181_Nico_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_plasmapTau181_Nico_Z_20260420.RDS")
result_PET_plasmapTau181_Jucker_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_plasmapTau181_Jucker_Z_20260420.RDS")
result_PET_Alamar_pTau181_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_Alamar_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_AlamarCSF_pTau181_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_AlamarCSF_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_LumipulseCSF_pTau181_Z <- readRDS( "../Alamar_AmyloidClocks/BootstrapResults/result_PET_LumipulseCSF_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_BD.pTau.181_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_BD.pTau.181_Z_20260420.RDS")
result_PET_BD.pTau.217_Z <- readRDS("../Alamar_AmyloidClocks/BootstrapResults/result_PET_BD.pTau.217_Z_20260420.RDS")



df$TimefromApos_CSFpT217_Z <- approx(y = result_PET_CSFpT217_Nico_Z$Time_to_Positivity, 
                                     x = result_PET_CSFpT217_Nico_Z$Estimate, xout = df$CSFpT217_Nico_Z)$y
df$TimefromApos_AlamarpT217_Z <- approx(y = result_PET_Alamar_pTau217_Z$Time_to_Positivity, 
                                        x = result_PET_Alamar_pTau217_Z$Estimate, xout = df$Alamar_pTau217_Z)$y
df$TimefromApos_Z <- approx(y = result_pib_Z$Time_to_Positivity, x = result_pib_Z$Estimate, xout = df$CL_Z)$y
df$TimefromApos_AlamarCSFpT217_Z <- approx(y = result_PET_AlamarCSF_pTau217_Z$Time_to_Positivity, 
                                           x = result_PET_AlamarCSF_pTau217_Z$Estimate, xout = df$AlamarCSF_pTau217_Z)$y

df$TimefromApos_CSFpT181_Z <- approx(y = result_PET_CSFpT181_Nico_Z$Time_to_Positivity, 
                                     x = result_PET_CSFpT181_Nico_Z$Estimate, xout = df$CSFpT181_Nico_Z)$y
df$TimefromApos_plasmapT181_Nico_Z <- approx(y = result_PET_plasmapTau181_Nico_Z$Time_to_Positivity, 
                                             x = result_PET_plasmapTau181_Nico_Z$Estimate, xout = df$plasmapTau181_Nico_Z)$y
df$TimefromApos_plasmapT217_Nico_Z <- approx(y = result_PET_plasmapTau217_Nico_Z$Time_to_Positivity, 
                                             x = result_PET_plasmapTau217_Nico_Z$Estimate, xout = df$plasmapTau217_Nico_Z)$y

df$TimefromApos_plasmapT181_Jucker_Z <- approx(y = result_PET_plasmapTau181_Jucker_Z$Time_to_Positivity, 
                                               x = result_PET_plasmapTau181_Jucker_Z$Estimate, xout = df$plasmapTau181_jucker_Z)$y
df$TimefromApos_AlamarpT181_Z <- approx(y = result_PET_Alamar_pTau181_Z$Time_to_Positivity, 
                                        x = result_PET_Alamar_pTau181_Z$Estimate, xout = df$Alamar_pTau181_Z)$y
df$TimefromApos_AlamarCSFpT181_Z <- approx(y = result_PET_AlamarCSF_pTau181_Z$Time_to_Positivity, 
                                           x = result_PET_AlamarCSF_pTau181_Z$Estimate, xout = df$AlamarCSF_pTau181_Z)$y
df$TimefromApos_LumipulseCSFpT181_Z <- approx(y = result_PET_LumipulseCSF_pTau181_Z$Time_to_Positivity, 
                                              x = result_PET_LumipulseCSF_pTau181_Z$Estimate, xout = df$LUMIPULSE_CSF_pTau_Z)$y
df$TimefromApos_BD.pTau.181_Z <- approx(y = result_PET_BD.pTau.181_Z$Time_to_Positivity, 
                                        x = result_PET_BD.pTau.181_Z$Estimate, xout = df$BD.pTau.181_Z)$y
df$TimefromApos_BD.pTau.217_Z <- approx(y = result_PET_BD.pTau.217_Z$Time_to_Positivity, 
                                        x = result_PET_BD.pTau.217_Z$Estimate, xout = df$BD.pTau.217_Z)$y
#now shifting to plotting these results

df_for_positivity <- df[, c("newid18", "visit", "TimefromBaseline", "DIAN_EYO", "VISITAGEc", 
                            "CSFpT181_Nico_Z", "CSFpT217_Nico_Z", "LUMIPULSE_CSF_pTau_Z",                      
                            "plasmapTau217_Nico_Z", "plasmapTau181_Nico_Z", "plasmapTau181_jucker_Z",                     
                            "Alamar_pTau181_Z", "Alamar_pTau217_Z", "AlamarCSF_pTau181_Z",                        
                            "AlamarCSF_pTau217_Z", "BD.pTau.217_Z" ,                             
                            "CL_Z",                   
                            "TimefromApos_AlamarpT217_Z","TimefromApos_Z","TimefromApos_AlamarCSFpT217_Z",             
                            "TimefromApos_CSFpT181_Z","TimefromApos_CSFpT217_Z",
                            "TimefromApos_plasmapT181_Nico_Z","TimefromApos_plasmapT217_Nico_Z","TimefromApos_plasmapT181_Jucker_Z",         
                            "TimefromApos_AlamarpT181_Z","TimefromApos_AlamarCSFpT181_Z","TimefromApos_LumipulseCSFpT181_Z",          
                            "TimefromApos_BD.pTau.217_Z")]

df_for_positivity$Apos_by_CSFpT181_Nico_Z <- ifelse(df_for_positivity$CSFpT181_Nico_Z > cp_PET_CSFpT181_Nico_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_CSFpT217_Nico_Z <- ifelse(df_for_positivity$CSFpT217_Nico_Z > cp_PET_CSFpT217_Nico_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_LUMIPULSE_CSF_pTau_Z <- ifelse(df_for_positivity$LUMIPULSE_CSF_pTau_Z > cp_PET_LumipulseCSF_pTau181_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_plasmapTau217_Nico_Z <- ifelse(df_for_positivity$plasmapTau217_Nico_Z > cp_PET_plasmapTau217_Nico_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_plasmapTau181_Nico_Z <- ifelse(df_for_positivity$plasmapTau181_Nico_Z > cp_PET_plasmapTau181_Nico_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_plasmapTau181_jucker_Z <- ifelse(df_for_positivity$plasmapTau181_jucker_Z > cp_PET_plasmapTau181_Jucker_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_AlamarpT181_Z <- ifelse(df_for_positivity$Alamar_pTau181_Z > cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_AlamarpT217_Z <- ifelse(df_for_positivity$Alamar_pTau217_Z > cp_PET_Alamar_pTau217_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_AlamarCSF_pT181_Z <- ifelse(df_for_positivity$AlamarCSF_pTau181_Z > cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_AlamarCSF_pT217_Z <- ifelse(df_for_positivity$AlamarCSF_pTau217_Z > cp_PET_AlamarCSF_pTau217_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_BD.pTau.217_Z <- ifelse(df_for_positivity$BD.pTau.217_Z > cp_PET_BD.pTau.217_obj[[3]]$optimal_cutpoint, 1, 0)
df_for_positivity$Apos_by_CL_Z <- ifelse(df_for_positivity$CL_Z > Apos_thresh_CL_Z, 1, 0)
saveRDS(df_for_positivity, "./Data/df_for_positivity_DIAN.RDS")

p_CSF_217 <- get_TimefromAposPlot(result_PET_CSFpT217_Nico_Z, cp_CSFpT217_Nico_Z$optimal_cutpoint, 
                                  "CSF pTau217 - Nico Z")
p_Alamar_217 <- get_TimefromAposPlot(result_PET_Alamar_pTau217_Z, cp_AlamarpT217_Z$optimal_cutpoint, 
                                     "Plasma pTau217 - Alamar Z")
p_plasma_217 <- get_TimefromAposPlot(result_PET_plasmapTau217_Nico_Z, cp_plasmapT217_Nico_Z$optimal_cutpoint, 
                                     "Plasma pTau217 - Nico Z")

p_pib  <- get_TimefromAposPlot(result_pib_Z, cp_pib_Z$optimal_cutpoint, 
                               "Amyloid PET - PiB Z")
p_AlamarCSF_217 <- get_TimefromAposPlot(result_PET_AlamarCSF_pTau217_Z, cp_AlamarCSFpT217_Z$optimal_cutpoint, 
                                        "CSF pTau217 - Alamar Z")
p_CSF_181 <- get_TimefromAposPlot(result_PET_CSFpT181_Nico_Z, cp_CSFpT181_Nico_Z$optimal_cutpoint, 
                                  "CSF pTau181 - Nico Z")
p_plasma_181 <- get_TimefromAposPlot(result_PET_plasmapTau181_Nico_Z, cp_plasmapT181_Nico_Z$optimal_cutpoint, 
                                     "Plasma pTau181 - Nico Z")
p_plasma_181_jucker <- get_TimefromAposPlot(result_PET_plasmapTau181_Jucker_Z, cp_plasmapT181_Jucker_Z$optimal_cutpoint, 
                                            "Plasma pTau181 - Jucker Z")
p_Alamar_181 <- get_TimefromAposPlot(result_PET_Alamar_pTau181_Z, cp_AlamarpT181_Z$optimal_cutpoint, 
                                     "Plasma pTau181 - Alamar Z")
p_AlamarCSF_181 <- get_TimefromAposPlot(result_PET_AlamarCSF_pTau181_Z, cp_AlamarCSFpT181_Z$optimal_cutpoint, 
                                        "CSF pTau181 - Alamar Z")
p_LumipulseCSF_181 <- get_TimefromAposPlot(result_PET_LumipulseCSF_pTau181_Z, cp_LumipulseCSFpT181_Z$optimal_cutpoint, 
                                           "CSF pTau181 - Lumipulse Z")
p_BD.pTau.217_Z <- get_TimefromAposPlot(result_PET_BD.pTau.217_Z, cp_BD.pTau.217_Z$optimal_cutpoint, 
                                        "CSF BD pTau217")
p_BD.pTau.181_Z <- get_TimefromAposPlot(result_PET_BD.pTau.181_Z, cp_BD.pTau.181_Z$optimal_cutpoint, 
                                        "CSF BD pTau181")


layout_matrix <- rbind(c(1, 1, 1, 1, 1, 1, 1, 1, 1), c(2, 2, 3, 3, 4,4, 5, 5, 5), c(6,6,6, 7,7,7, 8,8,8), c(9,9,9,10,10, 10,11,11, 11),
                       c(12, 12, 12, 12, 13, 13, 13, 13, 13))
grid.arrange(p_pib,
             p_CSF_217, p_AlamarCSF_217, p_Alamar_217, p_plasma_217,
             p_CSF_181, p_AlamarCSF_181, p_Alamar_181,
             p_LumipulseCSF_181, p_plasma_181, p_plasma_181_jucker,
             p_BD.pTau.181_Z, p_BD.pTau.217_Z,
             layout_matrix = layout_matrix)
graph2ppt(file = "./Figures/BootstrappedAmyloidTimeResults_20260420.pptx", width = 10, height = 9)

##################################################################################

##Getting observed conversion


p_error_CSFpT217_Nico <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                  value_col = "CSFpT217_Nico_Z", threshold = cp_PET_CSFpT217_Nico_obj[[3]]$optimal_cutpoint)
p_error_Alamar_pTau217 <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                   value_col = "Alamar_pTau217_Z", threshold = cp_PET_Alamar_pTau217_obj[[3]]$optimal_cutpoint)
p_error_pib <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", value_col = "CL_Z")
p_error_AlamarCSF_pTau217 <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                      value_col = "AlamarCSF_pTau217_Z", threshold = cp_PET_AlamarCSF_pTau217_obj[[3]]$optimal_cutpoint)
p_error_plasmapTau217_Nico <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                       value_col = "plasmapTau217_Nico_Z", threshold = cp_PET_plasmapTau217_Nico_obj[[3]]$optimal_cutpoint)

p_error_CSFpT181_Nico <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                  value_col = "CSFpT181_Nico_Z", threshold = cp_PET_CSFpT181_Nico_obj[[3]]$optimal_cutpoint)

p_error_plasmapTau181_Nico <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                       value_col = "plasmapTau181_Nico_Z", threshold = cp_PET_plasmapTau181_Nico_obj[[3]]$optimal_cutpoint)

p_error_plasmapTau181_jucker <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                         value_col = "plasmapTau181_jucker_Z", threshold = cp_PET_plasmapTau181_Jucker_obj[[3]]$optimal_cutpoint)
p_error_Alamar_pTau181 <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                   value_col = "Alamar_pTau181_Z", threshold = cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint)
p_error_AlamarCSF_pTau181 <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                      value_col = "AlamarCSF_pTau181_Z", threshold = cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint)
p_error_LumipulseCSF_pTau181 <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                         value_col = "LUMIPULSE_CSF_pTau_Z", threshold = cp_PET_LumipulseCSF_pTau181_obj[[3]]$optimal_cutpoint)

p_error_BD.pTau.181_Z <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                  value_col = "BD.pTau.181_Z", threshold = cp_PET_BD.pTau.181_obj[[3]]$optimal_cutpoint)

p_error_BD.pTau.217_Z <- plot_threshold_crossings(df, id_col = "newid18", age_col = "VISITAGEc", 
                                                  value_col = "BD.pTau.217_Z", threshold = cp_PET_BD.pTau.217_obj[[3]]$optimal_cutpoint)


# lemon::grid_arrange_shared_legend(p_error_pib[[1]],
#                                   p_error_CSFpT217_Nico[[1]],  p_error_AlamarCSF_pTau217[[1]], p_error_Alamar_pTau217[[1]], p_error_plasmapTau217_Nico[[1]],
#                                   p_error_CSFpT181_Nico[[1]], p_error_AlamarCSF_pTau181[[1]], p_error_Alamar_pTau181[[1]],
#                                   p_error_LumipulseCSF_pTau181[[1]],p_error_plasmapTau181_Nico[[1]], p_error_plasmapTau181_jucker[[1]],
#                                   p_error_BD.pTau.217_Z[[1]],
#                                   nrow = 3, ncol = 4)
# 
# graph2ppt(file = "./Figures/ObservedvsPredictedAgeatConversion.pptx", width = 10, height = 9)

results <- data.frame("Biomarker" = c("PiB", "CSFpT217_IPMS", "CSFpT217_Alamar",
                                      "PlasmapT217_Alamar", "PlasmapT217_IPMS",
                                      "CSFpT181_IPMS", "CSFpT181_Alamar", "PlasmapT181_Alamar",
                                      "CSFpT181_Lumipulse", "PlasmapT181_IPMS", "PlasmapT181_Simoa",
                                      "CSFBDpT181", "ADNI PET", "ADNI pTau217/AB42"),
                      "MAE" = c(p_error_pib[[2]],
                                p_error_CSFpT217_Nico[[2]],  p_error_AlamarCSF_pTau217[[2]], p_error_Alamar_pTau217[[2]], p_error_plasmapTau217_Nico[[2]],
                                p_error_CSFpT181_Nico[[2]], p_error_AlamarCSF_pTau181[[2]], p_error_Alamar_pTau181[[2]],
                                p_error_LumipulseCSF_pTau181[[2]],p_error_plasmapTau181_Nico[[2]], p_error_plasmapTau181_jucker[[2]],
                                p_error_BD.pTau.217_Z[[2]], p_error_PET[[2]], p_error_pT217[[2]]),
                      "Modality" = c("PET", "CSF", "CSF", "Plasma", "Plasma", 
                                     "CSF", "CSF", "Plasma", "CSF", "Plasma", "Plasma", "CSF", "PET", "Plasma"))



results %>%
  mutate(Biomarker = fct_reorder(Biomarker, MAE, .fun = mean)) %>%
  ggplot(aes(x = Biomarker, y = MAE, shape = Modality, colour = Modality)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_manual(values = c("black", "#007B82", "#F8766D")) +
  scale_shape_manual(values = c(19, 15,17)) + xlab("") + ylab("Mean Average Error (Years)")
graph2ppt(file = "./Figures/Evaluation_of_Real_Models_PointsvsMAE.pptx", width = 4.5, height = 5)

