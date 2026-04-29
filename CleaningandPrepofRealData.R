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
library(data.table)

library(emmeans) # for confidence intervals
library(export)

source("./functions.R")
source("./Simulation_functions.R")
source("./Fig1and2Funcs.R")

df <- read.csv("../Alamar_AmyloidClocks/Data/cleaned_df_withBatchNorming_20260420.csv") #not available to post on git because of data sharing restrictions
colnames(df)[16:17] <- c("CSFpT181_Nico", "CSFpT217_Nico")
df <- df[-15,]
cols_to_zscore <- c("CSFpT181_Nico", "CSFpT217_Nico", "LUMIPULSE_CSF_pTau",
                    "plasmapTau217_Nico", "plasmapTau181_Nico", "plasmapTau181_jucker",
                    "Alamar_pTau181", "Alamar_pTau217", "AlamarCSF_pTau181",
                    "AlamarCSF_pTau217", "BD.pTau.181", "BD.pTau.217", "CL")

for (col in cols_to_zscore) {
  df <- compute_zscore(df, !!sym(col))
}

df <- df[!is.na(df$VISITAGEc),]
df <- df %>%
  group_by(newid18) %>%
  mutate(TimefromBaseline = VISITAGEc - min(VISITAGEc, na.rm = TRUE)) %>%
  ungroup()



#Alamar values couldn't be normed to non-carriers because there weren't enough values.
df$Alamar_pTau217_Z <- (df$Alamar_pTau217 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau217, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau217, na.rm = TRUE)

df$Alamar_pTau181_Z <- (df$Alamar_pTau181 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau181, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau181, na.rm = TRUE)

df$AlamarCSF_pTau217_Z <- (df$Alamar_pTau217 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$AlamarCSF_pTau217, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau217, na.rm = TRUE)

df$AlamarCSF_pTau181_Z <- (df$Alamar_pTau181 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$AlamarCSF_pTau181, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau181, na.rm = TRUE)


df$BD.pTau.181_Z <- (df$BD.pTau.181 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.181, na.rm = TRUE)) /
  sd(df[!duplicated(df$newid18)  & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.181, na.rm = TRUE)

df$BD.pTau.217_Z <- (df$BD.pTau.217 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.217, na.rm = TRUE)) /
  sd(df[!duplicated(df$newid18)  & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.217, na.rm = TRUE)


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
BD.pTau.217_Z_PET_Thresh <- run_cutpointr_analysis(df, "BD.pTau.217_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                                   "Centiloids, Amyloid PET (Z)", "CSF Alamar BD pTau217 (Z)")
BD.pTau.181_Z_PET_Thresh <- run_cutpointr_analysis(df, "BD.pTau.181_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                                   "Centiloids, Amyloid PET (Z)", "CSF Alamar BD pTau181 (Z)")




# results$summary         # Precision, recall, TP/TN/FP/FN
# results$threshold_mean  # Mean optimal cutpoint
# results$roc_plot        # ROC curve
# results$classification_plot  # Classification visualization

layout <- rbind(c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4),
                c(5, 5, 5, 6, 6, 7, 7, 7, 8, 8),
                c(9, 9, 9, 10, 10, 11, 11, 11, 12, 12),
                c(13, 13, 13, 14, 14, 15, 15, 15, 16, 16),
                c(17, 17, 17, 18, 18, 19, 19, 19, 20, 20),
                c(21, 21, 21, 22, 22, 23, 23, 23, 24, 24))

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
             BD.pTau.217_Z_PET_Thresh$classification_plot + ggtitle("U."),
             BD.pTau.217_Z_PET_Thresh$roc_plot + ggtitle("V."),
             BD.pTau.181_Z_PET_Thresh$classification_plot + ggtitle("W."),
             BD.pTau.181_Z_PET_Thresh$roc_plot + ggtitle("X."),
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
RofC_BD.pTau.181_Z <- get_RofC_df(df[!is.na(df$BD.pTau.181_Z),], "newid18", "BD.pTau.181_Z", "TimefromBaseline", modelType = "E")
RofC_BD.pTau.217_Z <- get_RofC_df(df, "newid18", "BD.pTau.217_Z", "TimefromBaseline")


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
rel_accum_BD.pTau.181_Z <- mean(max(RofC_BD.pTau.181_Z[RofC_BD.pTau.181_Z$classification == 1,]$rate_of_change, na.rm = TRUE),
                                min(RofC_BD.pTau.181_Z[RofC_BD.pTau.181_Z$classification == 2,]$rate_of_change))
rel_accum_BD.pTau.217_Z <- mean(max(RofC_BD.pTau.217_Z[RofC_BD.pTau.217_Z$classification == 1,]$rate_of_change),
                                min(RofC_BD.pTau.217_Z[RofC_BD.pTau.217_Z$classification == 2,]$rate_of_change))

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
colnames(RofC_BD.pTau.181_Z)[2] <- "RofC_BD.pTau.181_Z"
colnames(RofC_BD.pTau.217_Z)[2] <- "RofC_BD.pTau.217_Z"


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
RofC_BD.pTau.181_Z <- combine_RofC_and_df(df, "BD.pTau.181_Z", RofC_BD.pTau.181_Z,"RofC_BD.pTau.181_Z", rel_accum_BD.pTau.181_Z )
RofC_BD.pTau.217_Z <- combine_RofC_and_df(df, "BD.pTau.217_Z", RofC_BD.pTau.217_Z,"RofC_BD.pTau.217_Z", rel_accum_BD.pTau.217_Z )


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
cp_BD.pTau.181_Z  <- cutpointr(RofC_BD.pTau.181_Z, BD.pTau.181_Z, ReliableAccumulator, 
                               method = maximize_metric, metric = youden, na.rm = TRUE)
cp_BD.pTau.217_Z  <- cutpointr(RofC_BD.pTau.217_Z, BD.pTau.217_Z, ReliableAccumulator, 
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
BD.pTau.181_Z_RelAccum <- run_cutpointr_analysis(RofC_BD.pTau.181_Z, "BD.pTau.181_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                 "BD.pTau.181_Z",
                                                 "BD CSF pTau 181 (Z)", "ARC Amyloid PET (Z)")
BD.pTau.217_Z_RelAccum <- run_cutpointr_analysis(RofC_BD.pTau.217_Z, "BD.pTau.217_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                                 "BD.pTau.217_Z",
                                                 "BD CSF pTau 217 (Z)", "ARC Amyloid PET (Z)")

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
df <- define_baseline_rel_to_RelAccum(df, "BD.pTau.181_Z", cp_BD.pTau.181_Z, "ReliableAccumThresh_BD.pTau.181_Z")
df <- define_baseline_rel_to_RelAccum(df, "BD.pTau.217_Z", cp_BD.pTau.217_Z, "ReliableAccumThresh_BD.pTau.217_Z")


ggplot(df, aes(x = TimefromBaseline, y = CL_Z, group = newid18)) + geom_line(colour = "grey88") +
  geom_line(data = df[df$ReliableAccumThresh ==1,], 
            aes(x = TimefromBaseline, y = CL_Z, group = newid18), colour = "black") +
  theme_bw() + xlab("Time from Baseline") + ylab("Z")


# Map dataframes and corresponding columns by name
RofC_list <- list(
  RofC_CSFpT181_Nico_Z,
  RofC_CSFpT217_Nico_Z,
  RofC_BD.pTau.181_Z,
  RofC_BD.pTau.217_Z,
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
  "RofC_BD.pTau.181_Z",
  "RofC_BD.pTau.217_Z",
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


RofC_CSFpT217_Nico_Z <- merge(RofC_CSFpT217_Nico_Z[, c("newid18", "RofC_CSFpT217_Nico_Z")], 
                              df[df$TimefromBaseline == 0 & !is.na(df$CSFpT217_Nico_Z),
                                 c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                   "CL_Z", "CSFpT217_Nico_Z")], by = "newid18", all = FALSE)
RofC_plasmapT217_Nico_Z <- merge(RofC_plasmapT217_Nico_Z[, c("newid18", "RofC_plasmapT217_Nico_Z")], 
                                 df[df$TimefromBaseline == 0 & !is.na(df$plasmapTau217_Nico_Z),
                                    c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                      "CL_Z", "plasmapTau217_Nico_Z")], by = "newid18", all = FALSE)
RofC_AlamarpT217_Z <- merge(RofC_AlamarpT217_Z[, c("newid18", "RofC_AlamarpT217_Z")], 
                            df[df$TimefromBaseline == 0 & !is.na(df$Alamar_pTau217_Z),
                               c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                 "CL_Z", "Alamar_pTau217_Z")], by = "newid18", all = FALSE)
RofC_AlamarCSFpT217_Z <- merge(RofC_AlamarCSFpT217_Z[, c("newid18", "RofC_AlamarCSFpT217_Z")], 
                               df[df$TimefromBaseline == 0 & !is.na(df$AlamarCSF_pTau217_Z),
                                  c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                    "CL_Z", "AlamarCSF_pTau217_Z")], by = "newid18", all = FALSE)
RofC_pib_Z <- merge(RofC_pib_Z[, c("newid18", "RofC_pib_Z")], 
                    df[df$TimefromBaseline == 0 & !is.na(df$CL_Z),
                       c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                         "CL_Z")], by = "newid18", all = FALSE)
RofC_CSFpT181_Nico_Z <- merge(RofC_CSFpT181_Nico_Z[, c("newid18", "RofC_CSFpT181_Nico_Z")], 
                              df[df$TimefromBaseline == 0 & !is.na(df$CSFpT181_Nico_Z),
                                 c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                   "CL_Z", "CSFpT181_Nico_Z")], by = "newid18", all = FALSE)
RofC_plasmapT181_Nico_Z <- merge(RofC_plasmapT181_Nico_Z[, c("newid18", "RofC_plasmapT181_Nico_Z")], 
                                 df[df$TimefromBaseline == 0 & !is.na(df$plasmapTau181_Nico_Z),
                                    c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                      "CL_Z", "plasmapTau181_Nico_Z")], by = "newid18", all = FALSE)
RofC_plasmapT181_Jucker_Z <- merge(RofC_plasmapT181_Jucker_Z[, c("newid18", "RofC_plasmapT181_Jucker_Z")], 
                                   df[df$TimefromBaseline == 0 & !is.na(df$plasmapTau181_jucker_Z),
                                      c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                        "CL_Z", "plasmapTau181_jucker_Z")], by = "newid18", all = FALSE)
RofC_AlamarpT181_Z <- merge(RofC_AlamarpT181_Z[, c("newid18", "RofC_AlamarpT181_Z")], 
                            df[df$TimefromBaseline == 0 & !is.na(df$Alamar_pTau181_Z),
                               c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                 "CL_Z", "Alamar_pTau181_Z")], by = "newid18", all = FALSE)
RofC_AlamarCSFpT181_Z <- merge(RofC_AlamarCSFpT181_Z[, c("newid18", "RofC_AlamarCSFpT181_Z")], 
                               df[df$TimefromBaseline == 0 & !is.na(df$AlamarCSF_pTau181_Z),
                                  c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                    "CL_Z", "AlamarCSF_pTau181_Z")], by = "newid18", all = FALSE)
RofC_LumipulseCSFpT181_Z <- merge(RofC_LumipulseCSFpT181_Z[, c("newid18", "RofC_LumipulseCSFpT181_Z")], 
                                  df[df$TimefromBaseline == 0 & !is.na(df$LUMIPULSE_CSF_pTau_Z),
                                     c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                       "CL_Z", "LUMIPULSE_CSF_pTau_Z")], by = "newid18", all = FALSE)
RofC_BD.pTau.181_Z <- merge(RofC_BD.pTau.181_Z[, c("newid18", "RofC_BD.pTau.181_Z")], 
                            df[df$TimefromBaseline == 0 & !is.na(df$BD.pTau.181_Z),
                               c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                 "CL_Z", "BD.pTau.181_Z")], by = "newid18", all = FALSE)
RofC_BD.pTau.217_Z <- merge(RofC_BD.pTau.217_Z[, c("newid18", "RofC_BD.pTau.217_Z")], 
                            df[df$TimefromBaseline == 0 & !is.na(df$BD.pTau.217_Z),
                               c("newid18", "VISITAGEc", "SEX", "RACE", "apoe", "Mutation", "fam_mutation",
                                 "CL_Z", "BD.pTau.217_Z")], by = "newid18", all = FALSE)


RofC_CSFpT217_Nico_Z <- get_relAccumFlag(RofC_CSFpT217_Nico_Z, "RofC_CSFpT217_Nico_Z",
                                         rel_accum_CSFpT217_Nico_Z)
RofC_plasmapT217_Nico_Z <- get_relAccumFlag(RofC_plasmapT217_Nico_Z, "RofC_plasmapT217_Nico_Z",
                                            rel_accum_plasmapT217_Nico_Z)
RofC_AlamarpT217_Z <- get_relAccumFlag(RofC_AlamarpT217_Z, "RofC_AlamarpT217_Z",
                                       rel_accum_AlamarpT217_Z)
RofC_AlamarCSFpT217_Z <- get_relAccumFlag(RofC_AlamarCSFpT217_Z, "RofC_AlamarCSFpT217_Z",
                                          rel_accum_AlamarCSFpT217_Z)
RofC_pib_Z <- get_relAccumFlag(RofC_pib_Z, "RofC_pib_Z",
                               rel_accum_pib_Z)
RofC_CSFpT181_Nico_Z <- get_relAccumFlag(RofC_CSFpT181_Nico_Z, "RofC_CSFpT181_Nico_Z",
                                         rel_accum_CSFpT181_Nico_Z)
RofC_plasmapT181_Nico_Z <- get_relAccumFlag(RofC_plasmapT181_Nico_Z, "RofC_plasmapT181_Nico_Z",
                                            rel_accum_plasmapT181_Nico_Z)
RofC_plasmapT181_Jucker_Z <- get_relAccumFlag(RofC_plasmapT181_Jucker_Z, "RofC_plasmapT181_Jucker_Z",
                                              rel_accum_plasmapT181_Jucker_Z)
RofC_AlamarpT181_Z <- get_relAccumFlag(RofC_AlamarpT181_Z, "RofC_AlamarpT181_Z",
                                       rel_accum_AlamarpT181_Z)
RofC_AlamarCSFpT181_Z <- get_relAccumFlag(RofC_AlamarCSFpT181_Z, "RofC_AlamarCSFpT181_Z",
                                          rel_accum_AlamarCSFpT181_Z)
RofC_LumipulseCSFpT181_Z <- get_relAccumFlag(RofC_LumipulseCSFpT181_Z, "RofC_LumipulseCSFpT181_Z",
                                             rel_accum_LumipulseCSFpT181_Z)
RofC_BD.pTau.181_Z <- get_relAccumFlag(RofC_BD.pTau.181_Z, "RofC_BD.pTau.181_Z",
                                       rel_accum_BD.pTau.181_Z)
RofC_BD.pTau.217_Z <- get_relAccumFlag(RofC_BD.pTau.217_Z, "RofC_BD.pTau.217_Z",
                                       rel_accum_BD.pTau.217_Z)



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
cp_BD.pTau.181_Z  <- cutpointr(RofC_BD.pTau.181_Z, BD.pTau.181_Z, ReliableAccumulator, 
                               method = maximize_metric, metric = youden, na.rm = TRUE)
cp_BD.pTau.217_Z  <- cutpointr(RofC_BD.pTau.217_Z, BD.pTau.217_Z, ReliableAccumulator, 
                               method = maximize_metric, metric = youden, na.rm = TRUE)


###############################################################################
## AIM 2: Derive estimates of time-from-amyloid positivity 
## based on each amyloid biomarker
###############################################################################
df$Apos_by_PET <- ifelse(df$CL > 18, 1, 0)

cp_PET_CSFpT217_Nico_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                     CSFpT217_Nico_Z, Apos_by_PET, 
                                                     "CSF pTau217 - Nico (Z)", "Amyloid Positive - 18 CL")

cp_PET_plasmapTau217_Nico_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                          plasmapTau217_Nico_Z, Apos_by_PET, 
                                                          "Plasma pTau217 - Nico (Z)", "Amyloid Positive - 18 CL")
cp_PET_Alamar_pTau217_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                      Alamar_pTau217_Z, Apos_by_PET, 
                                                      "Alamar pTau217 (Z)", "Amyloid Positive - 18 CL")
cp_PET_AlamarCSF_pTau217_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                         AlamarCSF_pTau217_Z, Apos_by_PET, 
                                                         "Alamar CSF pTau217 (Z)", "Amyloid Positive - 18 CL")
cp_PET_BD.pTau.181_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                   BD.pTau.181_Z, Apos_by_PET, 
                                                   "Alamar CSF BD pTau181 (Z)", "Amyloid Positive - 18 CL")
cp_PET_BD.pTau.217_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                   BD.pTau.217_Z, Apos_by_PET, 
                                                   "Alamar CSF BD pTau217 (Z)", "Amyloid Positive - 18 CL")
cp_PET_CSFpT181_Nico_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                     CSFpT181_Nico_Z, Apos_by_PET, 
                                                     "CSF pTau181 - Nico (Z)", "Amyloid Positive - 18 CL")
cp_PET_plasmapTau181_Nico_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                          plasmapTau181_Nico_Z, Apos_by_PET, 
                                                          "Plasma pTau181 - Nico (Z)", "Amyloid Positive - 18 CL")
cp_PET_plasmapTau181_Jucker_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                            plasmapTau181_jucker_Z, Apos_by_PET, 
                                                            "Plasma pTau181 - Jucker (Z)", "Amyloid Positive - 18 CL")
cp_PET_Alamar_pTau181_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                      Alamar_pTau181_Z, Apos_by_PET, 
                                                      "Alamar pTau181 (Z)", "Amyloid Positive - 18 CL")
cp_PET_AlamarCSF_pTau181_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                         AlamarCSF_pTau181_Z, Apos_by_PET, 
                                                         "Alamar CSF pTau181 (Z)", "Amyloid Positive - 18 CL")
cp_PET_LumipulseCSF_pTau181_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                            LUMIPULSE_CSF_pTau_Z, Apos_by_PET, 
                                                            "Lumipulse CSF pTau181 (Z)", "Amyloid Positive - 18 CL")

