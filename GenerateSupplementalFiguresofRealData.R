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

library(emmeans) # for confidence intervals
 library(export)

source("./functions.R")
source("./Simulation_functions.R")
source("./Fig1and2Funcs.R")

source("./CleaningandPrepofRealData.R")

#################################################################################
##MAKING SUPP FIGURE COMPONENTS
#################################################################################
#PANELS A & B
plots_CL <- plot_lmer_by_class(df = df, yvar = "CL_Z",
                               class_col = "ReliableAccumThresh", xlab = "Years from baseline",
                               ylab = "Cortical Amyloid (Z)")

plots_ReliableAccumThresh_Alamar_pTau217_Z <- plot_lmer_by_class(df = df,
                                                                 yvar = "Alamar_pTau217_Z", class_col = "ReliableAccumThresh_Alamar_pTau217_Z",
                                                                 xlab = "Years from baseline", ylab = "Alamar Plasma pTau217 (Z)")

plots_ReliableAccumThresh_Alamar_pTau181_Z <- plot_lmer_by_class(df = df,
                                                                 yvar = "Alamar_pTau181_Z", class_col = "ReliableAccumThresh_Alamar_pTau181_Z",
                                                                 xlab = "Years from baseline", ylab = "Alamar Plasma pTau181 (Z)")

plots_ReliableAccumThresh_AlamarCSF_pTau217_Z <- plot_lmer_by_class(df = df,
                                                                 yvar = "AlamarCSF_pTau217_Z", class_col = "ReliableAccumThresh_AlamarCSF_pTau217_Z",
                                                                 xlab = "Years from baseline", ylab = "Alamar CSF pTau217 (Z)")

plots_ReliableAccumThresh_AlamarCSF_pTau181_Z <- plot_lmer_by_class(df = df,
                                                                 yvar = "AlamarCSF_pTau181_Z", class_col = "ReliableAccumThresh_AlamarCSF_pTau181_Z",
                                                                 xlab = "Years from baseline", ylab = "Alamar CSF pTau181 (Z)")
plots_ReliableAccumThresh_LumipulseCSF_pTau181_Z <- plot_lmer_by_class(df = df,
                                                                    yvar = "LUMIPULSE_CSF_pTau_Z", class_col = "ReliableAccumThresh_LumipulseCSF_pTau181_Z",
                                                                    xlab = "Years from baseline", ylab = "Lumipulse CSF pTau181 (Z)")

plots_ReliableAccumThresh_plasmapTau217_Nico_Z <- plot_lmer_by_class(df = df,
                                                                     yvar = "plasmapTau217_Nico_Z", class_col = "ReliableAccumThresh_plasmapTau217_Nico_Z",
                                                                     xlab = "Years from baseline", ylab = "Plasma pTau217- Nico (Z)")

plots_ReliableAccumThresh_plasmapTau181_Nico_Z <- plot_lmer_by_class(df = df,
   yvar = "plasmapTau181_Nico_Z", class_col = "ReliableAccumThresh_plasmapTau181_Nico_Z",
   xlab = "Years from baseline", ylab = "Plasma pTau181- Nico (Z)")

plots_ReliableAccumThresh_plasmapTau181_Jucker_Z <- plot_lmer_by_class(df = df,
                                                                       yvar = "plasmapTau181_jucker_Z", class_col = "ReliableAccumThresh_plasmapTau181_Jucker_Z",
                                                                       xlab = "Years from baseline", ylab = "Plasma pTau181- Jucker (Z)")

plots_ReliableAccumThresh_CSFpT217_Nico_Z <- plot_lmer_by_class(df = df,
                                                                yvar = "CSFpT217_Nico_Z", class_col = "ReliableAccumThresh_CSFpT217_Nico_Z",
                                                                xlab = "Years from baseline", ylab = "CSF pTau217 (Z)")

plots_ReliableAccumThresh_CSFpT181_Nico_Z <- plot_lmer_by_class(df = df,
   yvar = "CSFpT181_Nico_Z", class_col = "ReliableAccumThresh_CSFpT181_Nico_Z",
   xlab = "Years from baseline", ylab = "CSF pTau181 (Z)")

plots_ReliableAccumThresh_CSFpT181_Lumipulse <- plot_lmer_by_class(df = df,
                                                                yvar = "LUMIPULSE_CSF_pTau_Z", class_col = "ReliableAccumThresh_LumipulseCSF_pTau181_Z",
                                                                xlab = "Years from baseline", ylab = "CSF pTau181 (Z)")

#panel c
plots_hist_CSFpT217 <- plot_rate_of_change(RofC_CSFpT217_Nico_Z, "RofC_CSFpT217_Nico_Z", rel_accum_CSFpT217_Nico_Z,
                    xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "CSFpT217_Nico",]$xi), 3), 
                    omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "CSFpT217_Nico",]$omega), 3), 
                    alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "CSFpT217_Nico",]$alpha), 3))
plots_hist_CSFpT181 <- plot_rate_of_change(RofC_CSFpT181_Nico_Z, "RofC_CSFpT181_Nico_Z", rel_accum_CSFpT181_Nico_Z,
                                           xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "CSFpT181_Nico",]$xi), 3), 
                                           omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "CSFpT181_Nico",]$omega), 3), 
                                           alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "CSFpT181_Nico",]$alpha), 3))
plots_hist_AlamarCSFpT217 <- plot_rate_of_change(RofC_AlamarCSFpT217_Z, "RofC_AlamarCSFpT217_Z", rel_accum_AlamarCSFpT217_Z,
                                           xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarCSFpT217",]$xi), 3), 
                                           omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarCSFpT217",]$omega), 3), 
                                           alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarCSFpT217",]$alpha), 3))
plots_hist_AlamarCSFpT181 <- plot_rate_of_change(RofC_AlamarCSFpT181_Z, "RofC_AlamarCSFpT181_Z", rel_accum_AlamarCSFpT181_Z,
                                                 xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarCSFpT181",]$xi), 3), 
                                                 omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarCSFpT181",]$omega), 3), 
                                                 alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarCSFpT181",]$alpha), 3))
plots_hist_CSFLumipulsepT181 <- plot_rate_of_change(RofC_LumipulseCSFpT181_Z, "RofC_LumipulseCSFpT181_Z", rel_accum_LumipulseCSFpT181_Z,
                                                 xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "LumipulseCSFpT181",]$xi), 3), 
                                                 omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "LumipulseCSFpT181",]$omega), 3), 
                                                 alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "LumipulseCSFpT181",]$alpha), 3))
plots_hist_plasmapT217 <- plot_rate_of_change(RofC_plasmapT217_Nico_Z, "RofC_plasmapT217_Nico_Z", rel_accum_plasmapT217_Nico_Z,
                    xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT217_Nico",]$xi), 3), 
                    omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT217_Nico",]$omega), 3), 
                    alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT217_Nico",]$alpha), 3))
plots_hist_plasmapT181_Jucker <- plot_rate_of_change(RofC_plasmapT181_Jucker_Z, "RofC_plasmapT181_Jucker_Z", rel_accum_plasmapT181_Jucker_Z,
                    xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT181_Jucker",]$xi), 3), 
                    omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT181_Jucker",]$omega), 3), 
                    alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT181_Jucker",]$alpha), 3))
plots_hist_plasmapT181_Nico <- plot_rate_of_change(RofC_plasmapT181_Nico_Z, "RofC_plasmapT181_Nico_Z", rel_accum_plasmapT181_Nico_Z,
                                              xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT181_Nico",]$xi), 3), 
                                              omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT181_Nico",]$omega), 3), 
                                              alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "plasmapT181_Nico",]$alpha), 3))
plots_hist_AlamarpT181 <- plot_rate_of_change(RofC_AlamarpT181_Z, "RofC_AlamarpT181_Z", rel_accum_AlamarpT181_Z,
                    xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarpT181",]$xi), 3), 
                    omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarpT181",]$omega), 3), 
                    alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarpT181",]$alpha), 3))
plots_hist_AlamarpT217 <- plot_rate_of_change(RofC_AlamarpT217_Z, "RofC_AlamarpT217_Z", rel_accum_AlamarpT217_Z,
                    xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarpT217",]$xi), 3), 
                    omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarpT217",]$omega), 3), 
                    alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "AlamarpT217",]$alpha), 3))
plots_hist_pib <- plot_rate_of_change(RofC_pib_Z, "RofC_pib_Z", rel_accum_pib_Z,
                    xi = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "pib",]$xi), 3), 
                    omega = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "pib",]$omega), 3), 
                    alpha = round(as.numeric(RofC_fits[RofC_fits$biomarker_col == "pib",]$alpha), 3))
#panel d
plots_ARCvBaseline_CSFpT217 <- get_RofC_vs_baseline_scatter(
  df = RofC_CSFpT217_Nico_Z,
  RofC_col = "RofC_CSFpT217_Nico_Z",
  biomarker_col = "CSFpT217_Nico_Z",
  thresh_obj = CSFpT217_Nico_Z_PET_Thresh,
  XLAB = "Baseline CSF pTau217 (Z)",
  YLAB = "ARC CSF pTau217 (Z/year)",
  rel_accum_thresh = rel_accum_CSFpT217_Nico_Z
)

plots_ARCvBaseline_CSFpT181 <- get_RofC_vs_baseline_scatter(
  df = RofC_CSFpT181_Nico_Z,
  RofC_col = "RofC_CSFpT181_Nico_Z",
  biomarker_col = "CSFpT181_Nico_Z",
  thresh_obj = CSFpT181_Nico_Z_PET_Thresh,
  XLAB = "Baseline CSF pTau181 (Z)",
  YLAB = "ARC CSF pTau181 (Z/year)",
  rel_accum_thresh = rel_accum_CSFpT181_Nico_Z
)

plots_ARCvBaseline_AlamarpT181 <- get_RofC_vs_baseline_scatter(
  df = RofC_AlamarpT181_Z,
  RofC_col = "RofC_AlamarpT181_Z",
  biomarker_col = "Alamar_pTau181_Z",
  thresh_obj = Alamar_pTau181_Z_PET_Thresh,
  XLAB = "Baseline Alamar pTau181 (Z)",
  YLAB = "ARC Alamar pTau181 (Z/year)",
  rel_accum_thresh = rel_accum_AlamarpT181_Z
)

plots_ARCvBaseline_AlamarpT217 <- get_RofC_vs_baseline_scatter(
  df = RofC_AlamarpT217_Z,
  RofC_col = "RofC_AlamarpT217_Z",
  biomarker_col = "Alamar_pTau217_Z",
  thresh_obj = Alamar_pTau217_Z_PET_Thresh,
  XLAB = "Baseline Alamar pTau217 (Z)",
  YLAB = "ARC Alamar pTau217 (Z/year)",
  rel_accum_thresh = rel_accum_AlamarpT217_Z
)

plots_ARCvBaseline_AlamarCSFpT181 <- get_RofC_vs_baseline_scatter(
  df = RofC_AlamarCSFpT181_Z,
  RofC_col = "RofC_AlamarCSFpT181_Z",
  biomarker_col = "AlamarCSF_pTau181_Z",
  thresh_obj = AlamarCSF_pTau181_Z_PET_Thresh,
  XLAB = "Baseline CSF Alamar pTau181 (Z)",
  YLAB = "ARC CSF Alamar pTau181 (Z/year)",
  rel_accum_thresh = rel_accum_AlamarCSFpT181_Z
)

plots_ARCvBaseline_LumipulseCSFpT181 <- get_RofC_vs_baseline_scatter(
  df = RofC_LumipulseCSFpT181_Z,
  RofC_col = "RofC_LumipulseCSFpT181_Z",
  biomarker_col = "LUMIPULSE_CSF_pTau_Z",
  thresh_obj = CSFpT181_Lumipulse_Z_PET_Thresh,
  XLAB = "Baseline CSF Lumipulse pTau181 (Z)",
  YLAB = "ARC CSF Lumipulse pTau181 (Z/year)",
  rel_accum_thresh = rel_accum_LumipulseCSFpT181_Z
)

plots_ARCvBaseline_AlamarCSFpT217 <- get_RofC_vs_baseline_scatter(
  df = RofC_AlamarCSFpT217_Z,
  RofC_col = "RofC_AlamarCSFpT217_Z",
  biomarker_col = "AlamarCSF_pTau217_Z",
  thresh_obj = AlamarCSF_pTau217_Z_PET_Thresh,
  XLAB = "Baseline CSF Alamar pTau217 (Z)",
  YLAB = "ARC CSF Alamar pTau217 (Z/year)",
  rel_accum_thresh = rel_accum_AlamarCSFpT217_Z
)

plots_ARCvBaseline_pib <- get_RofC_vs_baseline_scatter(
  df = RofC_pib_Z,
  RofC_col = "RofC_pib_Z",
  biomarker_col = "CL_Z",
  thresh_obj = Apos_thresh_CL_Z,
  XLAB = "Baseline Amyloid PET (Z)",
  YLAB = "ARC Amyloid PET (Z/year)",
  rel_accum_thresh = rel_accum_pib_Z
)


plots_ARCvBaseline_plasmapT217 <- get_RofC_vs_baseline_scatter(
  df = RofC_plasmapT217_Nico_Z,
  RofC_col = "RofC_plasmapT217_Nico_Z",
  biomarker_col = "plasmapTau217_Nico_Z",
  thresh_obj = plasmapTau217_Nico_Z_PET_Thresh,
  XLAB = "Baseline Plasma pTau217 (Z)",
  YLAB = "ARC Plasma pTau217 (Z/year)",
  rel_accum_thresh = rel_accum_plasmapT217_Nico_Z
)


plots_ARCvBaseline_plasmapT181_Jucker <- get_RofC_vs_baseline_scatter(
  df = RofC_plasmapT181_Jucker_Z,
  RofC_col = "RofC_plasmapT181_Jucker_Z",
  biomarker_col = "plasmapTau181_jucker_Z",
  thresh_obj = plasmapTau181_jucker_Z_PET_Thresh,
  XLAB = "Baseline Plasma pTau181 - Jucker (Z)",
  YLAB = "ARC Plasma pTau181 - Jucker (Z/year)",
  rel_accum_thresh = rel_accum_plasmapT181_Jucker_Z
)

plots_ARCvBaseline_plasmapT181_Nico <- get_RofC_vs_baseline_scatter(
  df = RofC_plasmapT181_Nico_Z,
  RofC_col = "RofC_plasmapT181_Nico_Z",
  biomarker_col = "plasmapTau181_Nico_Z",
  thresh_obj = plasmapTau181_Nico_Z_PET_Thresh,
  XLAB = "Baseline Plasma pTau181 - Nico (Z)",
  YLAB = "ARC Plasma pTau181 - Nico (Z/year)",
  rel_accum_thresh = rel_accum_plasmapT181_Nico_Z
)

#panel e


CSFpT217_Nico_Z_RelAccum[[3]]

CSFpT181_Nico_Z_RelAccum[[3]]

AlamarpT181_Z_RelAccum[[3]]

AlamarpT217_Z_RelAccum[[3]]

AlamarCSFpT181_Z_RelAccum[[3]]

LumipulseCSFpT181_Z_RelAccum[[3]]

AlamarCSFpT217_Z_RelAccum[[3]]

plasmapT181_Jucker_Z_RelAccum[[3]]

plasmapT181_Nico_Z_RelAccum[[3]]

plasmapT217_Nico_Z_RelAccum[[3]]

pib_Z_RelAccum[[3]]


layout <- rbind(c(1, 1, 1, 2, 2, 2),
                c(1, 1, 1, 2, 2, 2),
                c(3, 3, 4, 4, 5, 6),
                c(3, 3, 4, 4, 5, 7))


grid.arrange(plots_CL$full, plots_hist_pib,
             plots_ARCvBaseline_pib$withoutThresh, pib_Z_RelAccum[[3]],
             plots_ARCvBaseline_pib$withThresh, plots_CL$low, plots_CL$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_PIB.pptx", width = 10, height = 8)



grid.arrange(plots_ReliableAccumThresh_CSFpT181_Nico_Z$full, plots_hist_CSFpT181,
             plots_ARCvBaseline_CSFpT181$withoutThresh,  CSFpT181_Nico_Z_RelAccum[[3]],
             plots_ARCvBaseline_CSFpT181$withThresh, plots_ReliableAccumThresh_CSFpT181_Nico_Z$low,
             plots_ReliableAccumThresh_CSFpT181_Nico_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_CSFpT181_Nico.pptx", width = 10, height = 8)


grid.arrange(plots_ReliableAccumThresh_CSFpT181_Lumipulse$full, plots_hist_CSFLumipulsepT181,
             plots_ARCvBaseline_LumipulseCSFpT181$withoutThresh,  CSFpT181_Lumipulse_Z_PET_Thresh[[3]],
             plots_ARCvBaseline_LumipulseCSFpT181$withThresh, plots_ReliableAccumThresh_CSFpT181_Lumipulse$low,
             plots_ReliableAccumThresh_CSFpT181_Lumipulse$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_CSFpT181_Lumipulse.pptx", width = 10, height = 8)

grid.arrange(plots_ReliableAccumThresh_CSFpT217_Nico_Z$full, plots_hist_CSFpT217,
             plots_ARCvBaseline_CSFpT217$withoutThresh,  CSFpT217_Nico_Z_RelAccum[[3]],
             plots_ARCvBaseline_CSFpT217$withThresh, plots_ReliableAccumThresh_CSFpT217_Nico_Z$low,
             plots_ReliableAccumThresh_CSFpT217_Nico_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_CSFpT217_Nico.pptx", width = 10, height = 8)

grid.arrange(plots_ReliableAccumThresh_plasmapTau181_Jucker_Z$full, plots_hist_plasmapT181_Jucker,
             plots_ARCvBaseline_plasmapT181_Jucker$withoutThresh,  plasmapT181_Jucker_Z_RelAccum[[3]],
             plots_ARCvBaseline_plasmapT181_Jucker$withThresh, plots_ReliableAccumThresh_plasmapTau181_Jucker_Z$low,
             plots_ReliableAccumThresh_plasmapTau181_Jucker_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_plasmapTau181_Jucker_Z.pptx", width = 10, height = 8)

grid.arrange(plots_ReliableAccumThresh_plasmapTau181_Nico_Z$full, plots_hist_plasmapT181_Nico,
             plots_ARCvBaseline_plasmapT181_Nico$withoutThresh,  plasmapT181_Nico_Z_RelAccum[[3]],
             plots_ARCvBaseline_plasmapT181_Nico$withThresh, plots_ReliableAccumThresh_plasmapTau181_Nico_Z$low,
             plots_ReliableAccumThresh_plasmapTau181_Nico_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_plasmapTau181_Nico_Z.pptx", width = 10, height = 8)

grid.arrange(plots_ReliableAccumThresh_plasmapTau217_Nico_Z$full, plots_hist_plasmapT217,
             plots_ARCvBaseline_plasmapT217$withoutThresh,  plasmapT217_Nico_Z_RelAccum[[3]],
             plots_ARCvBaseline_plasmapT217$withThresh, plots_ReliableAccumThresh_plasmapTau217_Nico_Z$low,
             plots_ReliableAccumThresh_plasmapTau217_Nico_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_plasmapTau217_Z.pptx", width = 10, height = 8)

grid.arrange(plots_ReliableAccumThresh_AlamarCSF_pTau181_Z$full, plots_hist_AlamarCSFpT181,
             plots_ARCvBaseline_AlamarCSFpT181$withoutThresh,  AlamarCSFpT181_Z_RelAccum[[3]],
             plots_ARCvBaseline_AlamarCSFpT181$withThresh, plots_ReliableAccumThresh_AlamarCSF_pTau181_Z$low,
             plots_ReliableAccumThresh_AlamarCSF_pTau181_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_AlamarCSFpT181.pptx", width = 10, height = 8)


grid.arrange(plots_ReliableAccumThresh_AlamarCSF_pTau217_Z$full, plots_hist_AlamarCSFpT217,
             plots_ARCvBaseline_AlamarCSFpT217$withoutThresh,  AlamarCSFpT217_Z_RelAccum[[3]],
             plots_ARCvBaseline_AlamarCSFpT217$withThresh, plots_ReliableAccumThresh_AlamarCSF_pTau217_Z$low,
             plots_ReliableAccumThresh_AlamarCSF_pTau217_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_AlamarCSFpT217.pptx", width = 10, height = 8)

grid.arrange(plots_ReliableAccumThresh_Alamar_pTau181_Z$full, plots_hist_AlamarpT181,
             plots_ARCvBaseline_AlamarpT181$withoutThresh,  AlamarpT181_Z_RelAccum[[3]],
             plots_ARCvBaseline_AlamarpT181$withThresh, plots_ReliableAccumThresh_Alamar_pTau181_Z$low,
             plots_ReliableAccumThresh_Alamar_pTau181_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_AlamarpT181.pptx", width = 10, height = 8)


grid.arrange(plots_ReliableAccumThresh_Alamar_pTau217_Z$full, plots_hist_AlamarpT217,
             plots_ARCvBaseline_AlamarpT217$withoutThresh,  AlamarpT217_Z_RelAccum[[3]],
             plots_ARCvBaseline_AlamarpT217$withThresh, plots_ReliableAccumThresh_Alamar_pTau217_Z$low,
             plots_ReliableAccumThresh_Alamar_pTau217_Z$high,
             layout_matrix = layout)
graph2ppt(file = "./Figures/FullCharacterization_AlamarpT217.pptx", width = 10, height = 8)

