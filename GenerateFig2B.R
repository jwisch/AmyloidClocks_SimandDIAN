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

RofC_fits$modality <- c("CSF", "CSF", "CSF", "CSF", "Plasma", "Plasma",
                        "Plasma", "Plasma", "Plasma", 
                        "CSF", "CSF", "CSF", "PET")
RofC_fits$Analyte <- c("pT181", "pT217", "pT181", "pT217", "pT217", "pT181",
                       "pT181", "pT181", "pT217", "pT181",
                       "pT217", "pT181", "PIB")
RofC_fits$Assay <- c("IPMS", "IPMS", "Alamar", "Alamar","IPMS", "IPMS",
                     "Simoa", "Alamar", "Alamar", "Alamar",
                     "Alamar", "Lumiuplse", "PET")

ggplot(RofC_fits, aes(x = omega, y = xi, label = biomarker_col, colour = Assay,
                            shape = Analyte)) +
  geom_point(aes(size = alpha, alpha = 0.6)) + scale_colour_manual(values = c("#007B82", "#C77CFF",  "green", "black","#F8766D" )) +
  scale_shape_manual(values = c(19, 17,19)) +
  # geom_text(vjust = 1.3, hjust = -0.1) +
  facet_wrap(~modality) +
  theme_bw() +
  labs(y = "xi", x = "omega") +
  theme(legend.position = "bottom")


graph2ppt(file = "./Figures/3DPlot_RealData_Distribution_parameters.pptx", width = 5, height = 4)

