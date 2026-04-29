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

source("./CleaningandPrepofRealData.R")


result_PET_BD.pTau.181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_BD.pTau.181_obj[[3]]$optimal_cutpoint,
                                                               "newid18", "TimefromBaseline",
                                                               "BD.pTau.181_Z", num_bootstraps = 1000)
result_PET_BD.pTau.217_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_BD.pTau.217_obj[[3]]$optimal_cutpoint,
                                                             "newid18", "TimefromBaseline",
                                                             "BD.pTau.217_Z", num_bootstraps = 1000)
saveRDS(result_PET_BD.pTau.181_Z, "./BootstrapResults/result_PET_BD.pTau.181_Z_20260420.RDS")
saveRDS(result_PET_BD.pTau.217_Z, "./BootstrapResults/result_PET_BD.pTau.217_Z_20260420.RDS")

 result_PET_CSFpT217_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_CSFpT217_Nico_obj[[3]]$optimal_cutpoint,
                                                              "newid18", "TimefromBaseline",
                                                              "CSFpT217_Nico_Z", num_bootstraps = 1000)
 result_PET_plasmapTau217_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_plasmapTau217_Nico_obj[[3]]$optimal_cutpoint,
                                                                   "newid18", "TimefromBaseline",
                                                                   "plasmapTau217_Nico_Z", num_bootstraps = 1000)
 result_PET_Alamar_pTau217_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_Alamar_pTau217_obj[[3]]$optimal_cutpoint,
                                                               "newid18", "TimefromBaseline",
                                                               "Alamar_pTau217_Z", num_bootstraps = 1000)
result_PET_AlamarCSF_pTau217_Z <- robust_bootstrap_get_Time_to_Positivity(data.frame(df[!is.na(df$AlamarCSF_pTau217_Z),]), cp_PET_AlamarCSF_pTau217_obj[[3]]$optimal_cutpoint,
                                                                          "newid18", "TimefromBaseline",
                                                                          "AlamarCSF_pTau217_Z", num_bootstraps = 1000)
 result_pib_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (18 - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE),
                                                "newid18", "TimefromBaseline",
                                                "CL_Z", num_bootstraps = 1000)
 #
 saveRDS(result_PET_CSFpT217_Nico_Z, "./BootstrapResults/result_PET_CSFpT217_Nico_Z_20260420.RDS")
saveRDS(result_PET_plasmapTau217_Nico_Z, "./BootstrapResults/result_PET_plasmapTau217_Nico_Z_20260420.RDS")
saveRDS(result_PET_Alamar_pTau217_Z, "./BootstrapResults/result_PET_Alamar_pTau217_Z_batchCorrected_20260420.RDS")
saveRDS(result_pib_Z, "./BootstrapResults/result_pib_Z_20260420.RDS")
saveRDS(result_PET_AlamarCSF_pTau217_Z, "./BootstrapResults/result_PET_AlamarCSF_pTau217_Z_batchCorrected_20260420.RDS")
 #
 #
 #
 result_PET_CSFpT181_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_CSFpT181_Nico_obj[[3]]$optimal_cutpoint ,
                                                              "newid18", "TimefromBaseline",
                                                              "CSFpT181_Nico_Z", num_bootstraps = 1000)
 result_PET_plasmapTau181_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_plasmapTau181_Nico_obj[[3]]$optimal_cutpoint,
                                                                   "newid18", "TimefromBaseline",
                                                                   "plasmapTau181_Nico_Z", num_bootstraps = 1000)
 result_PET_plasmapTau181_Jucker_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_plasmapTau181_Jucker_obj[[3]]$optimal_cutpoint,
                                                                     "newid18", "TimefromBaseline",
                                                                     "plasmapTau181_jucker_Z", num_bootstraps = 1000)
 result_PET_Alamar_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint,
                                                             "newid18", "TimefromBaseline",
                                                               "Alamar_pTau181_Z", num_bootstraps = 1000)
result_PET_AlamarCSF_pTau181_Z <- robust_bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint,
                                                                          "newid18", "TimefromBaseline",
                                                                          "AlamarCSF_pTau181_Z", num_bootstraps = 1000)
 result_PET_LumipulseCSF_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_LumipulseCSF_pTau181_obj[[3]]$optimal_cutpoint,
                                                             "newid18", "TimefromBaseline",
                                                               "LUMIPULSE_CSF_pTau_Z", num_bootstraps = 1000)
 #
 saveRDS(result_PET_CSFpT181_Nico_Z, "./BootstrapResults/result_PET_CSFpT181_Nico_Z_20260420.RDS")
 saveRDS(result_PET_plasmapTau181_Nico_Z, "./BootstrapResults/result_PET_plasmapTau181_Nico_Z_20260420.RDS")
 saveRDS(result_PET_plasmapTau181_Jucker_Z, "./BootstrapResults/result_PET_plasmapTau181_Jucker_Z_20260420.RDS")
 saveRDS(result_PET_Alamar_pTau181_Z, "./BootstrapResults/result_PET_Alamar_pTau181_Z_batchCorrected_20260420.RDS")
saveRDS(result_PET_AlamarCSF_pTau181_Z, "./BootstrapResults/result_PET_AlamarCSF_pTau181_Z_batchCorrected_20260420.RDS")
 saveRDS( result_PET_LumipulseCSF_pTau181_Z, "./BootstrapResults/result_PET_LumipulseCSF_pTau181_Z_batchCorrected_20260420.RDS")
