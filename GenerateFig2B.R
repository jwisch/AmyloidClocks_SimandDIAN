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
source("./ADNIanalysis/ADNI_analysis_withCombat.R")


RofC_fits$modality <- c("CSF", "CSF", "CSF", "CSF", "Plasma", "Plasma",
                        "Plasma", "Plasma", "Plasma", 
                        "CSF", "CSF", "CSF", "PET")
RofC_fits$Analyte <- c("pT181", "pT217", "pT181", "pT217", "pT217", "pT181",
                       "pT181", "pT181", "pT217", "pT181",
                       "pT217", "pT181", "PIB")
RofC_fits$Assay <- c("IPMS", "IPMS", "Alamar", "Alamar","IPMS", "IPMS",
                     "Simoa", "Alamar", "Alamar", "Alamar",
                     "Alamar", "Lumiuplse", "PET")
RofC_fits$cohort <- "DIAN"
RofC_fits_ADNI$cohort <- "ADNI"

RofC_fits <- rbind(RofC_fits[, c("biomarker_col", "xi", "omega", "alpha", "modality", "Analyte", "Assay", "cohort")], 
                   RofC_fits_ADNI[, c("biomarker_col", "xi", "omega", "alpha", "modality", "Analyte", "Assay", "cohort")])
RofC_fits[RofC_fits$biomarker_col %in% c("pib", "PET"),]$Analyte <- "PET"
ggplot(RofC_fits, aes(x = omega, y = xi, label = biomarker_col, colour = Assay,
                            shape = Analyte)) +
  geom_point(aes(size = alpha, alpha = 0.6)) + #scale_colour_manual(values = c("#007B82", "#C77CFF",  "green", "black","#F8766D", "orange" )) +
  scale_shape_manual(values = c(15, 17,19)) +ylim(c(-4.5, 1)) + xlim(c(-0.05, 10)) +
  # geom_text(vjust = 1.3, hjust = -0.1) +
  facet_wrap(~modality) +
  theme_bw() +
  labs(y = "xi", x = "omega") +
  theme(legend.position = "bottom")


graph2pdf(file = "./Figures/3DPlot_RealData_Distribution_parameters.pdf", width = 5, height = 4)
library(ggbreak)
ggplot(RofC_fits, aes(x = omega, y = xi, label = biomarker_col, colour = Assay, shape = Analyte)) +
  geom_point(aes(size = alpha), alpha = 0.6) + 
  scale_shape_manual(values = c(15, 17, 19)) +
  facet_wrap(~modality) +
  theme_bw() +
  labs(y = "xi", x = "omega") +
  theme(legend.position = "bottom") + ylim(c(-4.5, 1)) + xlim(c(-0.05, 10)) +
  
  # ---- HERE IS THE AXIS BREAK ----
# This skips everything between -3 and -1 on the y-axis
scale_y_cut(breaks=c(-1, -2.9), which=c(1, 3), scales=c(3, 1), space=.5)
graph2pdf(file = "./Figures/3DPlot_RealData_Distribution_parameters_withBreak.pdf", width = 5, height = 4)

