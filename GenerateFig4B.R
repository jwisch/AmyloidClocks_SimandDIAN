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

df <- read.csv("./Data/cleaned_df_withBatchNorming_20260420.csv")#Not included for data sharing reasons
colnames(df)[16:17] <- c("CSFpT181_Nico", "CSFpT217_Nico")
df <- df[-15,]

cols <- which(names(df) == "CL") : which(names(df) == "BD.pTau.217")
cols <- cols[names(df)[cols] != "Study.Period"]


#Identifying cross-sectional and longituidnal counts
df_flags <- df %>%
  group_by(newid18) %>%
  summarize(
    across(
      all_of(cols),
      list(
        any = ~ sum(!is.na(.x)) >= 1,
        two_plus = ~ sum(!is.na(.x)) >= 2
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )


###############################################################################
## NORMING RELEVANT VALUES TO NON-CARRIERS TO GET Z SCORES
###############################################################################

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


#Alamar values couldn't be normed to non-carriers because there weren't enough values.
df$Alamar_pTau217_Z <- (df$Alamar_pTau217 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau217, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau217, na.rm = TRUE)

df$Alamar_pTau181_Z <- (df$Alamar_pTau181 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau181, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$Alamar_pTau181, na.rm = TRUE)

df$AlamarCSF_pTau217_Z <- (df$AlamarCSF_pTau217 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$AlamarCSF_pTau217, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$AlamarCSF_pTau217, na.rm = TRUE)

df$AlamarCSF_pTau181_Z <- (df$AlamarCSF_pTau181 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$AlamarCSF_pTau181, na.rm = TRUE))/
  sd(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$AlamarCSF_pTau181, na.rm = TRUE)


df$BD.pTau.181_Z <- (df$BD.pTau.181 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.181, na.rm = TRUE)) /
  sd(df[!duplicated(df$newid18)  & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.181, na.rm = TRUE)

df$BD.pTau.217_Z <- (df$BD.pTau.217 - mean(df[!duplicated(df$newid18) & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.217, na.rm = TRUE)) /
  sd(df[!duplicated(df$newid18)  & (df$Mutation == 0 | (df$Mutation == 1 & df$DIAN_EYO < -10)),]$BD.pTau.217, na.rm = TRUE)

###############################################################################
## AIM 1: 1.	Identify the initial thresholds for reliable positive accumulation 
## for each amyloid biomarker using imaging and fluid biomarker data.
###############################################################################

tmp <- df[!is.na(df$Alamar_pTau217_Z),]

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


GMM_Hist_plot <- function(df, title, rel_accum_thresh){
  p <- ggplot(df, aes(x = rate_of_change, fill = classification)) + 
    geom_histogram(alpha = 0.3, position = "identity") + theme_bw() +
    theme(legend.position = "bottom") + ggtitle(title) + xlab("Rate of Change") +
    geom_vline(xintercept = rel_accum_thresh, linetype = "dashed") +
    xlim(c(-2.5, 2.5))
  return(p)
}

GMM_Hist_plot(RofC_CSFpT217_Nico_Z, "CSF pTau217 - Nico Z", rel_accum_CSFpT217_Nico_Z)
GMM_Hist_plot(RofC_plasmapT217_Nico_Z, "Plasma pTau217 - Nico Z", rel_accum_plasmapT217_Nico_Z)
GMM_Hist_plot(RofC_AlamarpT217_Z, "Alamar pTau217 Z", rel_accum_AlamarpT217_Z)
GMM_Hist_plot(RofC_AlamarCSFpT217_Z, "CSF Alamar pTau217 Z", rel_accum_AlamarCSFpT217_Z)

GMM_Hist_plot(RofC_pib_Z, "PET - PiB Z", rel_accum_pib_Z)
GMM_Hist_plot(RofC_CSFpT181_Nico_Z, "CSF pTau181 - Nico Z", rel_accum_CSFpT181_Nico_Z)
GMM_Hist_plot(RofC_plasmapT181_Nico_Z, "Plasma pTau181 - Nico Z", rel_accum_plasmapT181_Nico_Z)
GMM_Hist_plot(RofC_plasmapT181_Jucker_Z, "Plasma pTau181 - Jucker Z", rel_accum_plasmapT181_Jucker_Z)
GMM_Hist_plot(RofC_AlamarpT181_Z, "Alamar pTau181 Z", rel_accum_AlamarpT181_Z)
GMM_Hist_plot(RofC_AlamarCSFpT181_Z, "CSF Alamar pTau181 Z", rel_accum_AlamarCSFpT181_Z)
GMM_Hist_plot(RofC_LumipulseCSFpT181_Z, "CSF Lumipulse pTau181 Z", rel_accum_LumipulseCSFpT181_Z)



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

cp_PET_BD.pTau.181_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                   BD.pTau.181_Z, Apos_by_PET, 
                                                   "Alamar CSF BD pTau181 (Z)", "Amyloid Positive - 18 CL")
cp_PET_BD.pTau.217_obj <- get_RofC_scatter_and_ROC(df[!duplicated(df$newid18),],
                                                   BD.pTau.217_Z, Apos_by_PET, 
                                                   "Alamar CSF BD pTau217 (Z)", "Amyloid Positive - 18 CL")




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
