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

df <- read.csv("./Data/cleaned_df_withBatchNorming_20260106.csv")
cols_to_zscore <- c("CSFpT181_Nico", "CSFpT217_Nico", "LUMIPULSE_CSF_pTau",
                    "plasmapTau217_Nico", "plasmapTau181_Nico", "plasmapTau181_jucker",
                    "Alamar_pTau181", "Alamar_pTau217", "AlamarCSF_pTau181",
                    "AlamarCSF_pTau217", "CL")

for (col in cols_to_zscore) {
  df <- compute_zscore(df, !!sym(col))
}

df <- df %>%
  group_by(newid18) %>%
  mutate(TimefromBaseline = VISITAGEc - min(VISITAGEc, na.rm = TRUE)) %>%
  ungroup()

#Alamar plasma values couldn't be normed to non-carriers because there weren't enough values.
df$Alamar_pTau217_Z <- (df$Alamar_pTau217 - mean(df[df$Mutation == 1 & !is.na(df$Alamar_pTau217) &
                                                      df$DIAN_EYO < -10,]$Alamar_pTau217))/
  sd(df[df$Mutation == 1 & !is.na(df$Alamar_pTau217) &
          df$DIAN_EYO < -10,]$Alamar_pTau217)

df$Alamar_pTau181_Z <- (df$Alamar_pTau217 - mean(df[df$Mutation == 1 & !is.na(df$Alamar_pTau181) &
                                                      df$DIAN_EYO < -10,]$Alamar_pTau181))/
  sd(df[df$Mutation == 1 & !is.na(df$Alamar_pTau181) &
          df$DIAN_EYO < -10,]$Alamar_pTau181)



###############################################################################
###############################################################################
##Fig 1: Biofluid vs CL and ROC curve
###############################################################################
###############################################################################
Apos_thresh_CL_Z <- approx(df$CL, df$CL_Z, 18)$y
df$Apos <- ifelse(df$CL_Z > Apos_thresh_CL_Z, 1, 0)

Alamar_pTau217_Z_PET_Thresh <- run_cutpointr_analysis(df, "Alamar_pTau217_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                  "Centiloids, Amyloid PET (Z)", "Alamar pTau217 (Z)")
Alamar_pTau181_Z_PET_Thresh <- run_cutpointr_analysis(df, "Alamar_pTau181_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                                      "Centiloids, Amyloid PET (Z)", "Alamar pTau181 (Z)")
AlamarCSF_pTau217_Z_PET_Thresh <- run_cutpointr_analysis(df, "AlamarCSF_pTau217_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                                      "Centiloids, Amyloid PET (Z)", "CSF Alamar pTau217 (Z)")
AlamarCSF_pTau181_Z_PET_Thresh <- run_cutpointr_analysis(df, "AlamarCSF_pTau181_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                                      "Centiloids, Amyloid PET (Z)", "CSF Alamar pTau181 (Z)")
CSFpT217_Nico_Z_PET_Thresh <- run_cutpointr_analysis(df, "CSFpT217_Nico_Z", "Apos", y1 = -3, y2 = 33, "CL_Z",
                                                      "Centiloids, Amyloid PET (Z)", "CSF pTau217 (Z)")
CSFpT181_Nico_Z_PET_Thresh <- run_cutpointr_analysis(df, "CSFpT181_Nico_Z", "Apos", y1 = -3, y2 = 23, "CL_Z",
                                                     "Centiloids, Amyloid PET (Z)", "CSF pTau181 (Z)")
plasmapTau181_Nico_Z_PET_Thresh <- run_cutpointr_analysis(df, "plasmapTau181_Nico_Z", "Apos", y1 = -3, y2 = 23, "CL_Z",
                                                     "Centiloids, Amyloid PET (Z)", "Plasma pTau181 - Nico (Z)")
plasmapTau181_jucker_Z_PET_Thresh <- run_cutpointr_analysis(df, "plasmapTau181_jucker_Z", "Apos", y1 = -3, y2 = 23, "CL_Z",
                                                     "Centiloids, Amyloid PET (Z)", "Plasma pTau181 - Jucker (Z)")
plasmapTau217_Nico_Z_PET_Thresh <- run_cutpointr_analysis(df, "plasmapTau217_Nico_Z", "Apos", y1 = -10, y2 =50, "CL_Z",
                                                     "Centiloids, Amyloid PET (Z)", "Plasma pTau217 - Nico (Z)")
CSFpT181_Lumipulse_Z_PET_Thresh <- run_cutpointr_analysis(df, "LUMIPULSE_CSF_pTau_Z", "Apos", y1 = -3, y2 = 23, "CL_Z",
                                                     "Centiloids, Amyloid PET (Z)", "CSF Lumipulse pTau181 (Z)")
# results$summary         # Precision, recall, TP/TN/FP/FN
# results$threshold_mean  # Mean optimal cutpoint
# results$roc_plot        # ROC curve
# results$classification_plot  # Classification visualization

layout <- rbind(c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4),
                c(5, 5, 5, 6, 6, 7, 7, 7, 8, 8),
                c(9, 9, 9, 10, 10, 11, 11, 11, 12, 12),
                c(13, 13, 13, 14, 14, 15, 15, 15, 16, 16),
                c(17, 17, 17, 18, 18, 19, 19, 19, 20, 20))

grid.arrange(CSFpT217_Nico_Z_PET_Thresh$classification_plot + ggtitle("A."), 
             CSFpT217_Nico_Z_PET_Thresh$roc_plot + ggtitle("B.") , 
             AlamarCSF_pTau217_Z_PET_Thresh$classification_plot + ggtitle("C."), 
             AlamarCSF_pTau217_Z_PET_Thresh$roc_plot + ggtitle("D.") ,              
             plasmapTau217_Nico_Z_PET_Thresh$classification_plot + ggtitle("E."), 
             plasmapTau217_Nico_Z_PET_Thresh$roc_plot  + ggtitle("F."),              
             Alamar_pTau217_Z_PET_Thresh$classification_plot + ggtitle("G."), 
             Alamar_pTau217_Z_PET_Thresh$roc_plot + ggtitle("H."), 
             CSFpT181_Nico_Z_PET_Thresh$classification_plot + ggtitle("I."), 
             CSFpT181_Nico_Z_PET_Thresh$roc_plot  + ggtitle("J."),             
             CSFpT181_Lumipulse_Z_PET_Thresh$classification_plot + ggtitle("K."),
             CSFpT181_Lumipulse_Z_PET_Thresh$roc_plot + ggtitle("L."),            
             AlamarCSF_pTau181_Z_PET_Thresh$classification_plot + ggtitle("M."), 
             AlamarCSF_pTau181_Z_PET_Thresh$roc_plot  + ggtitle("N."), 
             plasmapTau181_Nico_Z_PET_Thresh$classification_plot + ggtitle("O."), 
             plasmapTau181_Nico_Z_PET_Thresh$roc_plot  + ggtitle("P."), 
             plasmapTau181_jucker_Z_PET_Thresh$classification_plot + ggtitle("Q."), 
             plasmapTau181_jucker_Z_PET_Thresh$roc_plot  + ggtitle("R."),              
             Alamar_pTau181_Z_PET_Thresh$classification_plot + ggtitle("S."), 
             Alamar_pTau181_Z_PET_Thresh$roc_plot + ggtitle("T.") , 
             layout_matrix = layout)

graph2ppt(file = "./Figures/ThresholdsforApos.pptx", width = 18, height = 14)




###############################################################################
## Fig 2: Info on Reliable accumulation
###############################################################################
df <- df[!duplicated(df[, c("newid18", "visit"),]),]
RofC_CSFpT217_Nico_Z <- get_RofC_df(df, "newid18", "CSFpT217_Nico_Z", "TimefromBaseline")
RofC_plasmapT217_Nico_Z <- get_RofC_df(df, "newid18", "plasmapTau217_Nico_Z", "TimefromBaseline")
RofC_AlamarpT217_Z <- get_RofC_df(df, "newid18", "Alamar_pTau217_Z", "TimefromBaseline") 
RofC_AlamarCSFpT217_Z <- get_RofC_df(df, "newid18", "AlamarCSF_pTau217_Z", "TimefromBaseline")

RofC_pib_Z <- get_RofC_df(df, "newid18", "CL_Z", "TimefromBaseline")
RofC_CSFpT181_Nico_Z <- get_RofC_df(df, "newid18", "CSFpT181_Nico_Z", "TimefromBaseline")
RofC_plasmapT181_Nico_Z <- get_RofC_df(df, "newid18", "plasmapTau181_Nico_Z", "TimefromBaseline")
RofC_plasmapT181_Jucker_Z <- get_RofC_df(df, "newid18", "plasmapTau181_jucker_Z", "TimefromBaseline")
RofC_AlamarpT181_Z <- get_RofC_df(df, "newid18", "Alamar_pTau181_Z", "TimefromBaseline") 
RofC_AlamarCSFpT181_Z <- get_RofC_df(df, "newid18", "AlamarCSF_pTau181_Z", "TimefromBaseline")
RofC_LumipulseCSFpT181_Z <- get_RofC_df(df, "newid18", "LUMIPULSE_CSF_pTau", "TimefromBaseline")


rel_accum_CSFpT217_Nico_Z <- mean(max(RofC_CSFpT217_Nico_Z[RofC_CSFpT217_Nico_Z$classification == 1,]$rate_of_change),
                                  min(RofC_CSFpT217_Nico_Z[RofC_CSFpT217_Nico_Z$classification == 2,]$rate_of_change))
rel_accum_plasmapT217_Nico_Z <- mean(max(RofC_plasmapT217_Nico_Z[RofC_plasmapT217_Nico_Z$classification == 1,]$rate_of_change),
                                     min(RofC_plasmapT217_Nico_Z[RofC_plasmapT217_Nico_Z$classification == 2,]$rate_of_change))
rel_accum_AlamarpT217_Z <- mean(max(RofC_AlamarpT217_Z[RofC_AlamarpT217_Z$classification == 1,]$rate_of_change),
                                min(RofC_AlamarpT217_Z[RofC_AlamarpT217_Z$classification == 2,]$rate_of_change))
rel_accum_AlamarCSFpT217_Z <- mean(max(RofC_AlamarCSFpT217_Z[RofC_AlamarCSFpT217_Z$classification == 1,]$rate_of_change),
                                   min(RofC_AlamarCSFpT217_Z[RofC_AlamarCSFpT217_Z$classification == 2,]$rate_of_change))

rel_accum_pib_Z <- mean(max(RofC_pib_Z[RofC_pib_Z$classification == 1,]$rate_of_change),
                        min(RofC_pib_Z[RofC_pib_Z$classification == 2,]$rate_of_change))
rel_accum_CSFpT181_Nico_Z <- mean(max(RofC_CSFpT181_Nico_Z[RofC_CSFpT181_Nico_Z$classification == 1,]$rate_of_change),
                                  min(RofC_CSFpT181_Nico_Z[RofC_CSFpT181_Nico_Z$classification == 2,]$rate_of_change))
rel_accum_plasmapT181_Nico_Z <- mean(max(RofC_plasmapT181_Nico_Z[RofC_plasmapT181_Nico_Z$classification == 1,]$rate_of_change),
                                     min(RofC_plasmapT181_Nico_Z[RofC_plasmapT181_Nico_Z$classification == 2,]$rate_of_change))
rel_accum_plasmapT181_Jucker_Z <- mean(max(RofC_plasmapT181_Jucker_Z[RofC_plasmapT181_Jucker_Z$classification == 1,]$rate_of_change),
                                       min(RofC_plasmapT181_Jucker_Z[RofC_plasmapT181_Jucker_Z$classification == 2,]$rate_of_change))
rel_accum_AlamarpT181_Z <- mean(max(RofC_AlamarpT181_Z[RofC_AlamarpT181_Z$classification == 1,]$rate_of_change),
                                min(RofC_AlamarpT181_Z[RofC_AlamarpT181_Z$classification == 2,]$rate_of_change))
rel_accum_AlamarCSFpT181_Z <- mean(max(RofC_AlamarCSFpT181_Z[RofC_AlamarCSFpT181_Z$classification == 1,]$rate_of_change),
                                   min(RofC_AlamarCSFpT181_Z[RofC_AlamarCSFpT181_Z$classification == 2,]$rate_of_change))
rel_accum_LumipulseCSFpT181_Z <- mean(max(RofC_LumipulseCSFpT181_Z[RofC_LumipulseCSFpT181_Z$classification == 1,]$rate_of_change),
                                      min(RofC_LumipulseCSFpT181_Z[RofC_LumipulseCSFpT181_Z$classification == 2,]$rate_of_change))



colnames(RofC_CSFpT217_Nico_Z)[2] <- "RofC_CSFpT217_Nico_Z"
colnames(RofC_plasmapT217_Nico_Z)[2] <- "RofC_plasmapT217_Nico_Z"
colnames(RofC_AlamarpT217_Z)[2] <- "RofC_AlamarpT217_Z"
colnames(RofC_AlamarCSFpT217_Z)[2] <- "RofC_AlamarCSFpT217_Z"

colnames(RofC_pib_Z)[2] <- "RofC_pib_Z"
colnames(RofC_CSFpT181_Nico_Z)[2] <- "RofC_CSFpT181_Nico_Z"
colnames(RofC_plasmapT181_Nico_Z)[2] <- "RofC_plasmapT181_Nico_Z"
colnames(RofC_plasmapT181_Jucker_Z)[2] <- "RofC_plasmapT181_Jucker_Z"
colnames(RofC_AlamarpT181_Z)[2] <- "RofC_AlamarpT181_Z"
colnames(RofC_AlamarCSFpT181_Z)[2] <- "RofC_AlamarCSFpT181_Z"
colnames(RofC_LumipulseCSFpT181_Z)[2] <- "RofC_LumipulseCSFpT181_Z"


RofC_CSFpT217_Nico_Z <- combine_RofC_and_df(df, "CSFpT217_Nico_Z", RofC_CSFpT217_Nico_Z,"RofC_CSFpT217_Nico_Z", rel_accum_CSFpT217_Nico_Z )
RofC_plasmapT217_Nico_Z <- combine_RofC_and_df(df, "plasmapTau217_Nico_Z", RofC_plasmapT217_Nico_Z,"RofC_plasmapT217_Nico_Z", rel_accum_plasmapT217_Nico_Z )
RofC_AlamarpT217_Z <- combine_RofC_and_df(df, "Alamar_pTau217_Z", RofC_AlamarpT217_Z,"RofC_AlamarpT217_Z", rel_accum_AlamarpT217_Z )
RofC_AlamarCSFpT217_Z <- combine_RofC_and_df(df, "AlamarCSF_pTau217_Z", RofC_AlamarCSFpT217_Z,"RofC_AlamarCSFpT217_Z", rel_accum_AlamarCSFpT217_Z )
RofC_pib_Z <- combine_RofC_and_df(df, "CL_Z", RofC_pib_Z,"RofC_pib_Z", rel_accum_pib_Z )
RofC_CSFpT181_Nico_Z <- combine_RofC_and_df(df, "CSFpT181_Nico_Z", RofC_CSFpT181_Nico_Z,"RofC_CSFpT181_Nico_Z", rel_accum_CSFpT181_Nico_Z )
RofC_plasmapT181_Nico_Z <- combine_RofC_and_df(df, "plasmapTau181_Nico_Z", RofC_plasmapT181_Nico_Z,"RofC_plasmapT181_Nico_Z", rel_accum_plasmapT181_Nico_Z )
RofC_plasmapT181_Jucker_Z <- combine_RofC_and_df(df, "plasmapTau181_jucker_Z", RofC_plasmapT181_Jucker_Z,"RofC_plasmapT181_Jucker_Z", rel_accum_plasmapT181_Jucker_Z )
RofC_AlamarpT181_Z <- combine_RofC_and_df(df, "Alamar_pTau181_Z", RofC_AlamarpT181_Z,"RofC_AlamarpT181_Z", rel_accum_AlamarpT181_Z )
RofC_AlamarCSFpT181_Z <- combine_RofC_and_df(df, "AlamarCSF_pTau181_Z", RofC_AlamarCSFpT181_Z,"RofC_AlamarCSFpT181_Z", rel_accum_AlamarCSFpT181_Z )
RofC_LumipulseCSFpT181_Z <- combine_RofC_and_df(df, "LUMIPULSE_CSF_pTau_Z", RofC_LumipulseCSFpT181_Z,"RofC_LumipulseCSFpT181_Z", rel_accum_LumipulseCSFpT181_Z )


#Creating flag for reliable accumulation

cp_CSFpT217_Nico_Z <- cutpointr(RofC_CSFpT217_Nico_Z, CSFpT217_Nico_Z, ReliableAccumulator, 
                                method = maximize_metric, metric = youden, na.rm = TRUE)
cp_plasmapT217_Nico_Z <- cutpointr(RofC_plasmapT217_Nico_Z, plasmapTau217_Nico_Z, ReliableAccumulator, 
                                   method = maximize_metric, metric = youden, na.rm = TRUE)
cp_AlamarpT217_Z <- cutpointr(RofC_AlamarpT217_Z, Alamar_pTau217_Z, ReliableAccumulator, 
                              method = maximize_metric, metric = youden, na.rm = TRUE)
cp_AlamarCSFpT217_Z <- cutpointr(RofC_AlamarCSFpT217_Z, AlamarCSF_pTau217_Z, ReliableAccumulator, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)

cp_pib_Z <- cutpointr(RofC_pib_Z, CL_Z, ReliableAccumulator, 
                      method = maximize_metric, metric = youden, na.rm = TRUE)
cp_CSFpT181_Nico_Z <- cutpointr(RofC_CSFpT181_Nico_Z, CSFpT181_Nico_Z, ReliableAccumulator, 
                                method = maximize_metric, metric = youden, na.rm = TRUE)
cp_plasmapT181_Nico_Z <- cutpointr(RofC_plasmapT181_Nico_Z, plasmapTau181_Nico_Z, ReliableAccumulator, 
                                   method = maximize_metric, metric = youden, na.rm = TRUE)
cp_plasmapT181_Jucker_Z <- cutpointr(RofC_plasmapT181_Jucker_Z, plasmapTau181_jucker_Z, ReliableAccumulator, 
                                     method = maximize_metric, metric = youden, na.rm = TRUE)
cp_AlamarpT181_Z <- cutpointr(RofC_AlamarpT181_Z, Alamar_pTau181_Z, ReliableAccumulator, 
                              method = maximize_metric, metric = youden, na.rm = TRUE)
cp_AlamarCSFpT181_Z <- cutpointr(RofC_AlamarCSFpT181_Z, AlamarCSF_pTau181_Z, ReliableAccumulator, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
cp_LumipulseCSFpT181_Z <- cutpointr(RofC_LumipulseCSFpT181_Z, LUMIPULSE_CSF_pTau_Z, ReliableAccumulator, 
                                    method = maximize_metric, metric = youden, na.rm = TRUE)

CSFpT217_Nico_Z_RelAccum <- run_cutpointr_analysis(RofC_CSFpT217_Nico_Z, "CSFpT217_Nico_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                     "CSFpT217_Nico_Z",
                                                      "CSF pTau217 (Z)", "ARC CSF pTau217 (Z)")

CSFpT181_Nico_Z_RelAccum <- run_cutpointr_analysis(RofC_CSFpT181_Nico_Z, "CSFpT181_Nico_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                   "CSFpT181_Nico_Z",
                                                   "CSF pTau181 (Z)", "ARC CSF pTau181 (Z)")

AlamarpT181_Z_RelAccum <- run_cutpointr_analysis(RofC_AlamarpT181_Z, "Alamar_pTau181_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                   "Alamar_pTau181_Z",
                                                   "Alamar pTau181 (Z)", "ARC Alamar pTau181 (Z)")

AlamarpT217_Z_RelAccum <- run_cutpointr_analysis(RofC_AlamarpT217_Z, "Alamar_pTau217_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                 "Alamar_pTau217_Z",
                                                 "Alamar pTau217 (Z)", "ARC Alamar pTau217 (Z)")

AlamarCSFpT181_Z_RelAccum <- run_cutpointr_analysis(RofC_AlamarCSFpT181_Z, "AlamarCSF_pTau181_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                 "AlamarCSF_pTau181_Z",
                                                 "CSF Alamar pTau181 (Z)", "ARC CSF Alamar pTau181 (Z)")

AlamarCSFpT217_Z_RelAccum <- run_cutpointr_analysis(RofC_AlamarCSFpT217_Z, "AlamarCSF_pTau217_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                 "AlamarCSF_pTau217_Z",
                                                 "CSF Alamar pTau217 (Z)", "ARC CSF Alamar pTau217 (Z)")

LumipulseCSFpT181_Z_RelAccum <- run_cutpointr_analysis(RofC_LumipulseCSFpT181_Z, "LUMIPULSE_CSF_pTau_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                    "LUMIPULSE_CSF_pTau_Z",
                                                    "CSF Lumipulse pTau181 (Z)", "ARC CSF Lumipulse pTau217 (Z)")

plasmapT181_Jucker_Z_RelAccum <- run_cutpointr_analysis(RofC_plasmapT181_Jucker_Z, "plasmapTau181_jucker_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                   "plasmapTau181_jucker_Z",
                                                   "Plasma pTau181 - Jucker (Z)", "ARC Plasma pTau181 - Jucker (Z)")

plasmapT181_Nico_Z_RelAccum <- run_cutpointr_analysis(RofC_plasmapT181_Nico_Z, "plasmapTau181_Nico_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                        "plasmapTau181_Nico_Z",
                                                        "Plasma pTau181 - Nico (Z)", "ARC Plasma pTau181 - Nico (Z)")

plasmapT217_Nico_Z_RelAccum <- run_cutpointr_analysis(RofC_plasmapT217_Nico_Z, "plasmapTau217_Nico_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                        "plasmapTau217_Nico_Z",
                                                        "Plasma pTau217 - Nico (Z)", "ARC Plasma pTau217 - Nico (Z)")

pib_Z_RelAccum <- run_cutpointr_analysis(RofC_pib_Z, "CL_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                      "CL_Z",
                                                      "Amyloid PET (Z)", "ARC Amyloid PET (Z)")

df <- define_baseline_rel_to_RelAccum(df, "CL_Z", cp_pib_Z)
df <- define_baseline_rel_to_RelAccum(df, "Alamar_pTau217_Z", cp_AlamarpT217_Z, "ReliableAccumThresh_Alamar_pTau217_Z")
df <- define_baseline_rel_to_RelAccum(df, "AlamarCSF_pTau217_Z", cp_AlamarCSFpT217_Z, "ReliableAccumThresh_AlamarCSF_pTau217_Z")
df <- define_baseline_rel_to_RelAccum(df, "plasmapTau217_Nico_Z", cp_plasmapT217_Nico_Z, "ReliableAccumThresh_plasmapTau217_Nico_Z")
df <- define_baseline_rel_to_RelAccum(df, "plasmapTau181_Nico_Z", cp_plasmapT181_Nico_Z, "ReliableAccumThresh_plasmapTau181_Nico_Z")
df <- define_baseline_rel_to_RelAccum(df, "plasmapTau181_jucker_Z", cp_plasmapT181_Jucker_Z, "ReliableAccumThresh_plasmapTau181_Jucker_Z")
df <- define_baseline_rel_to_RelAccum(df, "Alamar_pTau181_Z", cp_AlamarpT181_Z, "ReliableAccumThresh_Alamar_pTau181_Z")
df <- define_baseline_rel_to_RelAccum(df, "AlamarCSF_pTau181_Z", cp_AlamarCSFpT181_Z, "ReliableAccumThresh_AlamarCSF_pTau181_Z")
df <- define_baseline_rel_to_RelAccum(df, "LUMIPULSE_CSF_pTau_Z", cp_LumipulseCSFpT181_Z, "ReliableAccumThresh_LumipulseCSF_pTau181_Z")
df <- define_baseline_rel_to_RelAccum(df, "CSFpT181_Nico_Z", cp_CSFpT181_Nico_Z, "ReliableAccumThresh_CSFpT181_Nico_Z")
df <- define_baseline_rel_to_RelAccum(df, "CSFpT217_Nico_Z", cp_CSFpT217_Nico_Z, "ReliableAccumThresh_CSFpT217_Nico_Z")


ggplot(df, aes(x = TimefromBaseline, y = CL_Z, group = newid18)) + geom_line(colour = "grey88") +
  geom_line(data = df[df$ReliableAccumThresh ==1,], 
            aes(x = TimefromBaseline, y = CL_Z, group = newid18), colour = "black") +
  theme_bw() + xlab("Time from Baseline") + ylab("Z")


# Map dataframes and corresponding columns by name
RofC_list <- list(
  RofC_CSFpT181_Nico_Z,
  RofC_CSFpT217_Nico_Z,
  RofC_plasmapT217_Nico_Z,
  RofC_plasmapT181_Nico_Z,
  RofC_plasmapT181_Jucker_Z,
  RofC_AlamarpT181_Z,
  RofC_AlamarpT217_Z,
  RofC_AlamarCSFpT181_Z,
  RofC_AlamarCSFpT217_Z,
  RofC_LumipulseCSFpT181_Z,
  RofC_pib_Z
)

RofC_cols <- c(
  "RofC_CSFpT181_Nico_Z",
  "RofC_CSFpT217_Nico_Z",
  "RofC_pTau217_Nico_Z",
  "RofC_pTau181_Nico_Z",
  "RofC_pTau181_Jucker_Z",
  "RofC_AlamarpT181_Z",
  "RofC_AlamarpT217_Z",
  "RofC_AlamarCSFpT181_Z",
  "RofC_AlamarCSFpT217_Z",
  "RofC_LumipulseCSFpT181_Z",
  "RofC_pib_Z"
)

# Fit skew-normal to each RofC column


RofC_fits <- purrr::map(RofC_list, function(df_i) {
  # Dynamically detect the "RofC_" column
  rofc_col <- grep("^RofC_", names(df_i), value = TRUE)
  
  # Handle cases with no matching or multiple matches
  if (length(rofc_col) != 1) {
    return(tibble(
      biomarker_col = paste(rofc_col, collapse = ", "),
      xi = NA_real_,
      omega = NA_real_,
      alpha = NA_real_,
      fit_status = "no or multiple RofC_ columns found",
      n_values = NA_integer_
    ))
  }
  
  vals <- df_i[[rofc_col]]
  vals <- vals[!is.na(vals)]
  
  biomarker_name <- gsub("^RofC_|_Z$", "", rofc_col)
  fit_status <- "ok"
  xi <- omega <- alpha <- NA_real_
  
  if (length(vals) <= 10) {
    fit_status <- "too few values"
  } else if (sd(vals) == 0) {
    fit_status <- "zero variance"
  } else {
    fit <- tryCatch(selm(vals ~ 1), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) {
      fit_status <- "fit failed"
    } else {
      params <- tryCatch(coef(fit, "DP"), error = function(e) NULL)
      if (is.null(params)) {
        fit_status <- "param extraction failed"
      } else {
        xi <- params["xi"]
        omega <- params["omega"]
        alpha <- params["alpha"]
      }
    }
  }
  
  tibble(
    biomarker_col = biomarker_name,
    xi = xi,
    omega = omega,
    alpha = alpha,
    fit_status = fit_status,
    n_values = length(vals)
  )
}) %>%
  bind_rows()



#################################################################################
##MAKING FIGURE 2 COMPONENTS
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

#TODO: Fix p2 so that i just have highlighted dots instead of different facets

p1 <- ggplot(RofC_pib_Z, aes(x = RofC_pib_Z)) +
  xlim(c(-2.5, 2.5)) +
  xlab("Rate of Change (Z/year)") +
  # Smoothed density scaled to histogram height
  geom_density(aes(y = after_stat(density * n * diff(range(x)) / 100)), colour = "#007B82", fill = "#007B82",
               alpha = 0.1) +
  geom_density(data = RofC_plasmapT217_Nico_Z, 
               aes(x = RofC_plasmapT217_Nico_Z, y = after_stat(density * n * diff(range(x)) / 100)),
               colour = "#F8766D", fill = "#F8766D",
               alpha = 0.1) +
  geom_density(data = RofC_AlamarCSFpT217_Z, 
               aes(x = RofC_AlamarCSFpT217_Z, y = after_stat(density * n * diff(range(x)) / 100)),
               colour = "#C77CFF", fill = "#C77CFF",
               alpha = 0.1) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) 

RofC_fits$highlighter <- 0
RofC_fits$highlighter[RofC_fits$biomarker_col == "plasmapT217_Nico"] <- 1
RofC_fits$highlighter[RofC_fits$biomarker_col == "AlamarCSFpT217"] <- 2
RofC_fits$highlighter[RofC_fits$biomarker_col == "pib"] <- 3
RofC_fits$highlighter <- as.factor(RofC_fits$highlighter)
RofC_fits$modality <- c("CSF", "CSF", "Plasma", "Plasma",
                        "Plasma", "Plasma", "Plasma", 
                        "CSF", "CSF", "CSF", "PET")
p2 <- ggplot(RofC_fits, aes(x = omega, y = xi, label = biomarker_col, colour = highlighter,
                            shape = modality)) +
  geom_point(size = 3) + scale_colour_manual(values = c("black", "#F8766D", "#C77CFF", "#007B82")) +
  scale_shape_manual(values = c(19, 3,17)) +
  geom_text(vjust = 1.3, hjust = -0.1) +
  facet_wrap(~cut(alpha, 4)) +
  theme_bw() +
  labs(y = "xi", x = "omega", title = "Relationships across alpha bins") +
  theme(legend.position = "none")

layout_matrix = rbind(c(1), c(2), c(2))

grid.arrange(p1, p2, layout_matrix = layout_matrix)

graph2ppt(file = "./Figures/3DPlot_RealData_Distribution_parameters.pptx", width = 5, height = 4)

tmp <- setDT(df[!is.na(df$plasmapTau181_Nico_Z),])[, .N, by = newid18]
nrow(tmp[tmp$N > 1,])
