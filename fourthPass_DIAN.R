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
library(cowplot)
source("./functions.R")

df <- read.csv("../Alamar_AmyloidClocks/Data/cleaned_df_withBatchNorming_20260420.csv")
colnames(df)[16:17] <- c("CSFpT181_Nico", "CSFpT217_Nico")
df <- df[-15,]

cols <- which(names(df) == "CL") : which(names(df) == "BD.pTau.217")
cols <- cols[names(df)[cols] != "Study.Period"]
# ggpairs(df[!duplicated(df$newid18), cols])
# 
# graph2ppt(file = "./Figures/ggpairPlotAllBaselineData.pptx", width = 9, height = 9)


#Identifying cross-sectional and longituidnal counts
df_flags <- df %>%
  group_by(newid18) %>%
  summarize(
    across(
      CL:DIAN_EYO,
      list(
        any = ~ sum(!is.na(.x)) >= 1,
        two_plus = ~ sum(!is.na(.x)) >= 2
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  select(-Study.Period_any, -Study.Period_two_plus)


sum(df_flags$plasmapTau217_Nico_two_plus)
df_flags <- merge(df[!duplicated(df$newid18), c("newid18", "VISITAGEc", "SEX", "RACE", "apoe",
                                                "Mutation", "fam_mutation", "DIAN_EYO")], df_flags, 
                  by = "newid18")

df_flags$SEX <- as.factor(df_flags$SEX)  
levels(df_flags$SEX) <- c("male", "female")
df_flags$fam_mutation <-as.factor(df_flags$fam_mutation)
levels(df_flags$fam_mutation) <- c("PSEN1", "PSEN2", "APP")

myVars <- c("VISITAGEc", "DIAN_EYO", "N", "SEX", "RACE", "apoe", "fam_mutation", 
            names(df_flags)[9:(length(df_flags) - 2)])
factorVars <- c("SEX", "RACE", "apoe", "fam_mutation")

tab <- CreateTableOne(vars = myVars, strata = "Mutation", factorVars = factorVars, data = df_flags)

## Convert TableOne → data frame
tab_df <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

## Write to Excel
write.csv(tab_df, file = "./Figures/TableOne.csv")
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

# tmp <- df[!is.na(df$Alamar_pTau217_Z),]

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

robust_bootstrap_get_Time_to_Positivity <- function(
    df,
    PET_pos_threshold,
    id_name,
    time_name,
    value_name,
    num_bootstraps = 1000,
    bootstrap_percent = 0.8,
    degree = 3,
    printIter = TRUE
) {
  
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  unique_ids <- unique(df[[id_name]])
  n_ids <- length(unique_ids)
  
  df_res <- vector("list", num_bootstraps)
  
  # --- SAFE ITERATION ---
  safe_iteration <- function(i) {
    tryCatch({
      
      # Sample IDs
      sampled_idx <- sample(
        seq_len(n_ids),
        size = floor(bootstrap_percent * n_ids),
        replace = TRUE
      )
      sampled_ids <- unique_ids[sampled_idx]
      sampled_data <- df[df[[id_name]] %in% sampled_ids, ]
      
      # Model
      result <- get_Time_to_Positivity(
        sampled_data,
        id_name,
        time_name,
        value_name,
        PET_pos_threshold,
        degree
      )
      
      # Validate
      if (is.null(result) ||
          length(result$Time_to_Positivity) < 2 ||
          length(result$actual_predicted_val) < 2 ||
          all(is.na(result$Time_to_Positivity)) ||
          all(is.na(result$actual_predicted_val))) {
        return(NULL)
      }
      
      # Interpolation
      interpolated_val <- suppressWarnings(
        sapply(Time_Window, function(t) {
          approx(
            result$Time_to_Positivity,
            result$actual_predicted_val,
            xout = t,
            rule = 1,
            ties = mean
          )$y
        })
      )
      
      if (all(is.na(interpolated_val))) return(NULL)
      
      data.frame(
        Time_Window = Time_Window,
        interpolated_val = interpolated_val
      )
      
    }, error = function(e) {
      if (printIter) message("Iteration ", i, " failed: ", e$message)
      return(NULL)
    })
  }
  
  # --- RUN LOOP (this WILL go to num_bootstraps) ---
  for (i in seq_len(num_bootstraps)) {
    if (printIter) print(i)
    
    # try() catches anything tryCatch misses (important)
    res_i <- try(safe_iteration(i), silent = TRUE)
    
    if (inherits(res_i, "try-error") || is.null(res_i)) {
      df_res[[i]] <- NULL
    } else {
      df_res[[i]] <- res_i
    }
  }
  
  # --- Diagnostics ---
  num_failed <- sum(sapply(df_res, is.null))
  message("Completed all iterations.")
  message("Failed iterations: ", num_failed, " / ", num_bootstraps)
  
  df_res_clean <- df_res[!sapply(df_res, is.null)]
  
  if (length(df_res_clean) == 0) {
    warning("All iterations failed.")
    return(data.frame(
      Time_to_Positivity = numeric(0),
      Estimate = numeric(0),
      CI_Lower = numeric(0),
      CI_Upper = numeric(0)
    ))
  }
  
  # --- Aggregate ---
  bootstrap_matrix <- do.call(rbind, df_res_clean)
  
  mean_result <- bootstrap_matrix %>%
    dplyr::group_by(Time_Window) %>%
    dplyr::summarize(
      interpolated_val_mean = median(interpolated_val, na.rm = TRUE),
      .groups = "drop"
    )
  
  sd_result <- bootstrap_matrix %>%
    dplyr::group_by(Time_Window) %>%
    dplyr::summarize(
      interpolated_val_sd = sd(interpolated_val, na.rm = TRUE),
      .groups = "drop"
    )
  
  ci_calc <- merge(mean_result, sd_result, by = "Time_Window")
  
  ci_df <- data.frame(
    Time_to_Positivity = ci_calc$Time_Window,
    Estimate = ci_calc$interpolated_val_mean,
    CI_Lower = ci_calc$interpolated_val_mean - 1.96 * ci_calc$interpolated_val_sd,
    CI_Upper = ci_calc$interpolated_val_mean + 1.96 * ci_calc$interpolated_val_sd
  )
  
  # --- Threshold alignment ---
  try({
    if (all(is.finite(ci_df$Estimate)) &&
        any(ci_df$Estimate < PET_pos_threshold) &&
        any(ci_df$Estimate > PET_pos_threshold)) {
      
      adjustment <- approx(
        ci_df$Estimate,
        ci_df$Time_to_Positivity,
        xout = PET_pos_threshold,
        ties = mean
      )$y
      
      if (!is.na(adjustment)) {
        ci_df$Time_to_Positivity <- ci_df$Time_to_Positivity - adjustment
      }
    }
  }, silent = TRUE)
  
  return(ci_df[, c("Time_to_Positivity", "Estimate", "CI_Lower", "CI_Upper")])
}
# 
# result_PET_BD.pTau.181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_BD.pTau.181_obj[[3]]$optimal_cutpoint,
#                                                                "newid18", "TimefromBaseline",
#                                                                "BD.pTau.181_Z", num_bootstraps = 1000)
# result_PET_BD.pTau.217_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_BD.pTau.217_obj[[3]]$optimal_cutpoint,
#                                                              "newid18", "TimefromBaseline",
#                                                              "BD.pTau.217_Z", num_bootstraps = 1000)
# saveRDS(result_PET_BD.pTau.181_Z, "./BootstrapResults/result_PET_BD.pTau.181_Z_20260420.RDS")
# saveRDS(result_PET_BD.pTau.217_Z, "./BootstrapResults/result_PET_BD.pTau.217_Z_20260420.RDS")
# 
#  result_PET_CSFpT217_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_CSFpT217_Nico_obj[[3]]$optimal_cutpoint,
#                                                               "newid18", "TimefromBaseline",
#                                                               "CSFpT217_Nico_Z", num_bootstraps = 1000)
#  result_PET_plasmapTau217_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_plasmapTau217_Nico_obj[[3]]$optimal_cutpoint,
#                                                                    "newid18", "TimefromBaseline",
#                                                                    "plasmapTau217_Nico_Z", num_bootstraps = 1000)
#  result_PET_Alamar_pTau217_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_Alamar_pTau217_obj[[3]]$optimal_cutpoint,
#                                                                "newid18", "TimefromBaseline",
#                                                                "Alamar_pTau217_Z", num_bootstraps = 1000)
 # result_PET_AlamarCSF_pTau217_Z <- robust_bootstrap_get_Time_to_Positivity(data.frame(df[!is.na(df$AlamarCSF_pTau217_Z),]), cp_PET_AlamarCSF_pTau217_obj[[3]]$optimal_cutpoint,
 #                                                                 "newid18", "TimefromBaseline",
 #                                                                 "AlamarCSF_pTau217_Z", num_bootstraps = 1000)
#  result_pib_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (18 - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE),
#                                                 "newid18", "TimefromBaseline",
#                                                 "CL_Z", num_bootstraps = 1000)
#  # 
#  saveRDS(result_PET_CSFpT217_Nico_Z, "./BootstrapResults/result_PET_CSFpT217_Nico_Z_20260420.RDS")
# saveRDS(result_PET_plasmapTau217_Nico_Z, "./BootstrapResults/result_PET_plasmapTau217_Nico_Z_20260420.RDS")
# saveRDS(result_PET_Alamar_pTau217_Z, "./BootstrapResults/result_PET_Alamar_pTau217_Z_batchCorrected_20260420.RDS")
# saveRDS(result_pib_Z, "./BootstrapResults/result_pib_Z_20260420.RDS")
 # saveRDS(result_PET_AlamarCSF_pTau217_Z, "./BootstrapResults/result_PET_AlamarCSF_pTau217_Z_batchCorrected_20260420.RDS")
#  # 
#  # 
#  # 
#  result_PET_CSFpT181_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_CSFpT181_Nico_obj[[3]]$optimal_cutpoint ,
#                                                               "newid18", "TimefromBaseline",
#                                                               "CSFpT181_Nico_Z", num_bootstraps = 1000)
#  result_PET_plasmapTau181_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_plasmapTau181_Nico_obj[[3]]$optimal_cutpoint,
#                                                                    "newid18", "TimefromBaseline",
#                                                                    "plasmapTau181_Nico_Z", num_bootstraps = 1000)
#  result_PET_plasmapTau181_Jucker_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_plasmapTau181_Jucker_obj[[3]]$optimal_cutpoint,
#                                                                      "newid18", "TimefromBaseline",
#                                                                      "plasmapTau181_jucker_Z", num_bootstraps = 1000)
#  result_PET_Alamar_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint,
#                                                              "newid18", "TimefromBaseline",
#                                                                "Alamar_pTau181_Z", num_bootstraps = 1000)
 # result_PET_AlamarCSF_pTau181_Z <- robust_bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint,
 #                                                             "newid18", "TimefromBaseline",
 #                                                               "AlamarCSF_pTau181_Z", num_bootstraps = 1000)
#  result_PET_LumipulseCSF_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), cp_PET_LumipulseCSF_pTau181_obj[[3]]$optimal_cutpoint,
#                                                              "newid18", "TimefromBaseline",
#                                                                "LUMIPULSE_CSF_pTau_Z", num_bootstraps = 1000)
#  # 
#  saveRDS(result_PET_CSFpT181_Nico_Z, "./BootstrapResults/result_PET_CSFpT181_Nico_Z_20260420.RDS")
#  saveRDS(result_PET_plasmapTau181_Nico_Z, "./BootstrapResults/result_PET_plasmapTau181_Nico_Z_20260420.RDS")
#  saveRDS(result_PET_plasmapTau181_Jucker_Z, "./BootstrapResults/result_PET_plasmapTau181_Jucker_Z_20260420.RDS")
#  saveRDS(result_PET_Alamar_pTau181_Z, "./BootstrapResults/result_PET_Alamar_pTau181_Z_batchCorrected_20260420.RDS")
# saveRDS(result_PET_AlamarCSF_pTau181_Z, "./BootstrapResults/result_PET_AlamarCSF_pTau181_Z_batchCorrected_20260420.RDS")
#  saveRDS( result_PET_LumipulseCSF_pTau181_Z, "./BootstrapResults/result_PET_LumipulseCSF_pTau181_Z_batchCorrected_20260420.RDS")
 
#######HERE. UPDATE WHAT WE READ IN
result_PET_CSFpT217_Nico_Z <- readRDS("./BootstrapResults/result_PET_CSFpT217_Nico_Z_20260420.RDS")
result_PET_Alamar_pTau217_Z <- readRDS("./BootstrapResults/result_PET_Alamar_pTau217_Z_batchCorrected_20260420.RDS")
result_pib_Z <- readRDS( "./BootstrapResults/result_pib_Z_20260420.RDS")
result_PET_AlamarCSF_pTau217_Z <- readRDS("./BootstrapResults/result_PET_AlamarCSF_pTau217_Z_batchCorrected_20260420.RDS")
result_PET_plasmapTau217_Nico_Z <- readRDS("./BootstrapResults/result_PET_plasmapTau217_Nico_Z_20260420.RDS")

result_PET_CSFpT181_Nico_Z <- readRDS("./BootstrapResults/result_PET_CSFpT181_Nico_Z_20260420.RDS")
result_PET_plasmapTau181_Nico_Z <- readRDS("./BootstrapResults/result_PET_plasmapTau181_Nico_Z_20260420.RDS")
result_PET_plasmapTau181_Jucker_Z <- readRDS("./BootstrapResults/result_PET_plasmapTau181_Jucker_Z_20260420.RDS")
result_PET_Alamar_pTau181_Z <- readRDS("./BootstrapResults/result_PET_Alamar_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_AlamarCSF_pTau181_Z <- readRDS("./BootstrapResults/result_PET_AlamarCSF_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_LumipulseCSF_pTau181_Z <- readRDS( "./BootstrapResults/result_PET_LumipulseCSF_pTau181_Z_batchCorrected_20260420.RDS")
result_PET_BD.pTau.181_Z <- readRDS("./BootstrapResults/result_PET_BD.pTau.181_Z_20260420.RDS")
result_PET_BD.pTau.217_Z <- readRDS("./BootstrapResults/result_PET_BD.pTau.217_Z_20260420.RDS")



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

get_TimefromAposPlot <- function(result_df, cp, YLAB){
 p <- ggplot(result_df[result_df$Estimate > cp,], 
         aes(x = Time_to_Positivity, y = Estimate, 
             ymin = CI_Lower, ymax = CI_Upper)) +
    geom_line() + geom_ribbon(alpha = 0.3) + theme_bw() + 
    xlab("Estimated time from A+") + ylab(YLAB) +
    xlim(c(-10, 20)) 
 return(p)
  
}

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

get_threshold_crossings <- function(df, id_col, age_col, value_col, threshold = 2.6) {
  # Ensure dplyr is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it with install.packages('dplyr').")
  }
  
  df %>%
    dplyr::arrange(.data[[id_col]], .data[[age_col]]) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      cross_age = {
        valid <- stats::complete.cases(.data[[value_col]], .data[[age_col]])
        x <- .data[[age_col]][valid]
        y <- .data[[value_col]][valid]
        
        if (any(y < threshold) & any(y > threshold)) {
          stats::approx(x = y, y = x, xout = threshold)$y
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )
}

plot_threshold_crossings <- function(df, id_col, age_col, value_col, threshold = 2.6) {
  
  # Get threshold crossings
  crossings <- get_threshold_crossings(
    df = df,
    id_col = id_col,
    age_col = age_col,
    value_col = value_col,
    threshold = threshold
  )
  
  # Merge crossings with original data
  df <- merge(df, crossings[!is.na(crossings$cross_age), ], 
              by.x = id_col, by.y = id_col, all = TRUE)
  
  # Compute observed crossing age
  df$observed_cross_age <- df[[age_col]] - df$TimefromApos_Z
  
  # Subset for unique IDs for plotting
  df_unique <- df[!duplicated(df[[id_col]]), ]
  
  # Define color limits
  eyo_vals <- df_unique$DIAN_EYO[!is.na(df_unique$cross_age)]
  color_limits <- c(floor(min(eyo_vals)), ceiling(max(eyo_vals)))
  # Calculate MAE and RMSE for annotation
  mae_val <- round(MLmetrics::MAE(df_unique$observed_cross_age[!is.na(df_unique$cross_age) & !is.na(df_unique$observed_cross_age)], 
                       df_unique$cross_age[!is.na(df_unique$cross_age) & !is.na(df_unique$observed_cross_age)]), 2)
  rmse_val <- round(MLmetrics::RMSE(df_unique$observed_cross_age[!is.na(df_unique$cross_age) & !is.na(df_unique$observed_cross_age)], 
                         df_unique$cross_age[!is.na(df_unique$cross_age) & !is.na(df_unique$observed_cross_age)]), 2)
  
  # Build plot
  p <- ggplot(df_unique, aes(x = observed_cross_age, y = cross_age, group = .data[[id_col]], colour = DIAN_EYO)) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +
    xlim(c(23, 57)) +
    ylim(c(23, 57)) +
    scale_color_gradient(low = "#00008B", high = "red", 
                         limits = color_limits,
                         oob = scales::squish,
                         name = "EYO") +
    ylab("Predicted Age at Conversion") +
    xlab("Observed Age at Conversion") +
    theme(legend.position = "bottom") +
    annotate(
      geom = "label",
      x = 38,
      y = 33,
      hjust = 0,
      vjust = 1,
      label = paste0("MAE = ", mae_val, " years\nRMSE = ", rmse_val, " years"),
      fill = "#FFD700",
      color = "black",
      label.size = 0.3,
      label.r = unit(0.15, "lines")
    ) + geom_point()
  
  return(list(p, mae_val))
}

id_col <- "newid18"
age_col <- "VISITAGEc"
value_col <- "Alamar_pTau181_Z"
threshold <- cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint

crossings <- get_threshold_crossings(
  df = df,
  id_col = id_col,
  age_col = age_col,
  value_col = value_col,
  threshold = threshold
)

crossings_CSF <- get_threshold_crossings(
  df = df,
  id_col = id_col,
  age_col = age_col,
  value_col = "AlamarCSF_pTau181_Z",
  threshold = cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint
)


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


lemon::grid_arrange_shared_legend(p_error_pib[[1]],
                                  p_error_CSFpT217_Nico[[1]],  p_error_AlamarCSF_pTau217[[1]], p_error_Alamar_pTau217[[1]], p_error_plasmapTau217_Nico[[1]],
                                  p_error_CSFpT181_Nico[[1]], p_error_AlamarCSF_pTau181[[1]], p_error_Alamar_pTau181[[1]],
                                  p_error_LumipulseCSF_pTau181[[1]],p_error_plasmapTau181_Nico[[1]], p_error_plasmapTau181_jucker[[1]],
                                  p_error_BD.pTau.181_Z[[1]], p_error_BD.pTau.217_Z[[1]],
                                  layout_matrix = layout_matrix)

graph2ppt(file = "./Figures/ObservedvsPredictedAgeatConversion.pptx", width = 10, height = 9)

results <- data.frame("Biomarker" = c("PiB", "CSFpT217_IPMS", "CSFpT217_Alamar",
                           "PlasmapT217_Alamar", "PlasmapT217_IPMS",
                           "CSFpT181_IPMS", "CSFpT181_Alamar", "PlasmapT181_Alamar",
                           "CSFpT181_Lumipulse", "PlasmapT181_IPMS", "PlasmapT181_Simoa",
                           "CSFBDpT181"),
           "MAE" = c(p_error_pib[[2]],
                     p_error_CSFpT217_Nico[[2]],  p_error_AlamarCSF_pTau217[[2]], p_error_Alamar_pTau217[[2]], p_error_plasmapTau217_Nico[[2]],
                     p_error_CSFpT181_Nico[[2]], p_error_AlamarCSF_pTau181[[2]], p_error_Alamar_pTau181[[2]],
                     p_error_LumipulseCSF_pTau181[[2]],p_error_plasmapTau181_Nico[[2]], p_error_plasmapTau181_jucker[[2]],
                     p_error_BD.pTau.217_Z[[2]]),
           "Modality" = c("PET", "CSF", "CSF", "Plasma", "Plasma", 
                          "CSF", "CSF", "Plasma", "CSF", "Plasma", "Plasma", "CSF"))

library(forcats)
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
  scale_shape_manual(values = c(19, 3,17)) + xlab("") + ylab("Mean Average Error (Years)")
graph2ppt(file = "./Figures/Evaluation_of_Real_Models_PointsvsMAE.pptx", width = 4.5, height = 5)


grid.arrange(p_pib, p_error_pib,
             p_CSF_217, p_error_CSFpT217_Nico, 
             p_AlamarCSF_217,p_error_AlamarCSF_pTau217,
             p_Alamar_217, p_error_Alamar_pTau217, 
             p_plasma_217, p_error_plasmapTau217_Nico, 
             
             p_CSF_181, p_error_CSFpT181_Nico, 
             p_AlamarCSF_181,p_error_AlamarCSF_pTau181, 
             p_Alamar_181, p_error_Alamar_pTau181,
             p_LumipulseCSF_181, p_error_LumipulseCSF_pTau181,
             p_plasma_181, p_error_plasmapTau181_Nico, 
             p_plasma_181_jucker,p_error_plasmapTau181_jucker, 
             p_BD.pTau.217_Z, p_error_BD.pTau.217_Z, ncol = 2)

# --- Step 1: Extract the legend from one error plot ---
shared_legend <- get_legend(
  p_error_pib + theme(legend.position = "bottom")
)

# --- Step 2: Remove legends from ALL individual plots ---
strip_legend <- function(p) p + theme(legend.position = "none")

plots_no_legend <- list(
  strip_legend(p_pib),                strip_legend(p_error_pib),
  strip_legend(p_CSF_217),            strip_legend(p_error_CSFpT217_Nico),
  strip_legend(p_AlamarCSF_217),      strip_legend(p_error_AlamarCSF_pTau217),
  strip_legend(p_Alamar_217),         strip_legend(p_error_Alamar_pTau217),
  strip_legend(p_plasma_217),         strip_legend(p_error_plasmapTau217_Nico),
  strip_legend(p_CSF_181),            strip_legend(p_error_CSFpT181_Nico),
  strip_legend(p_AlamarCSF_181),      strip_legend(p_error_AlamarCSF_pTau181),
  strip_legend(p_Alamar_181),         strip_legend(p_error_Alamar_pTau181),
  strip_legend(p_LumipulseCSF_181),   strip_legend(p_error_LumipulseCSF_pTau181),
  strip_legend(p_plasma_181),         strip_legend(p_error_plasmapTau181_Nico),
  strip_legend(p_plasma_181_jucker),  strip_legend(p_error_plasmapTau181_jucker),
  strip_legend(p_BD.pTau.217_Z), strip_legend(p_error_BD.pTau.217_Z)
)

# --- Step 3: Arrange all plots in a 2-column layout ---
plots_grid <- arrangeGrob(
  grobs = plots_no_legend,
  ncol = 2
)

# --- Step 4: Add the shared legend below the grid ---
final_plot <- grid.arrange(
  plots_grid,
  shared_legend,
  ncol = 1,
  heights = c(10, 0.7)   # adjust as needed
)

final_plot
graph2ppt(file = "./Figures/Evaluation_of_Real_Models_20260420.pptx", width = 4.5, height = 10)




plots_no_legend <- list(
  strip_legend(p_pib),                strip_legend(p_error_pib),
  strip_legend(p_AlamarCSF_217),      strip_legend(p_error_AlamarCSF_pTau217),
  strip_legend(p_plasma_217),         strip_legend(p_error_plasmapTau217_Nico)
)

# --- Step 3: Arrange all plots in a 2-column layout ---
plots_grid <- arrangeGrob(
  grobs = plots_no_legend,
  ncol = 2
)

# --- Step 4: Add the shared legend below the grid ---
final_plot <- grid.arrange(
  plots_grid,
  shared_legend,
  ncol = 1,
  heights = c(10, 0.7)   # adjust as needed
)

final_plot
graph2ppt(file = "./Figures/Evaluation_of_Real_Models_selectedFew.pptx", width = 4.5, height = 4)

################################################################################
##TRANSLATE Z SCORE TO TIME FROM A+
################################################################################
#WHY IS IT SO HARD FOR ME TO FIGURE OUT HOW TO STACK BOXES?
#TO DO: MAKE A NICE LOOKING PLOT HERE.....


result_PET_CSFpT217_Nico_Z$group <- "CSFpT217"
result_PET_Alamar_pTau217_Z$group <- "pT217_Alamar"
result_PET_AlamarCSF_pTau217_Z$group <- "pT217_AlamarCSF"
result_pib_Z$group <- "PET"
result_PET_CSFpT181_Nico_Z$group <- "CSFpT181"
result_PET_plasmapTau217_Nico_Z$group <- "PlasmapT217_Nico"
result_PET_plasmapTau181_Nico_Z$group <- "PlasmapT181_Nico"
result_PET_plasmapTau181_Jucker_Z$group <- "PlasmapT181_Jucker"
result_PET_Alamar_pTau181_Z$group <- "pT181_Alamar"
result_PET_AlamarCSF_pTau181_Z$group <- "pT181_AlamarCSF"
result_PET_LumipulseCSF_pTau181_Z$group <- "pT181_Lumipulse"
result_PET_BD.pTau.217_Z$group <- "BDpT217"
Time_Vec <- c(-4, -2, 0, 2, 4, 6, 8)

Time_to_Z <- rbind(Estimate_to_CI(Time_Vec, result_PET_CSFpT217_Nico_Z),
                    Estimate_to_CI(Time_Vec, result_PET_Alamar_pTau217_Z),
                    Estimate_to_CI(Time_Vec, result_PET_AlamarCSF_pTau217_Z),
                   Estimate_to_CI(Time_Vec, result_PET_plasmapTau217_Nico_Z),
                   
                    Estimate_to_CI(Time_Vec, result_pib_Z),
                    Estimate_to_CI(Time_Vec, result_PET_CSFpT181_Nico_Z),
                   Estimate_to_CI(Time_Vec, result_PET_plasmapTau181_Nico_Z),
                    Estimate_to_CI(Time_Vec, result_PET_plasmapTau181_Jucker_Z),
                    Estimate_to_CI(Time_Vec,  result_PET_Alamar_pTau181_Z),
                   Estimate_to_CI(Time_Vec,  result_PET_AlamarCSF_pTau181_Z),
                   Estimate_to_CI(Time_Vec,  result_PET_BD.pTau.217_Z))

colnames(Time_to_Z) <- c( "mean", "lower", "upper","range", "labeltext")

Time_to_Z$upper <- ifelse(Time_to_Z$upper < Time_to_Z$mean, 22, Time_to_Z$upper)

Time_to_Z |>
  forestplot(
    mean = mean,
    lower = lower,
    upper = upper,
    labeltext = labeltext,
    fn.ci_norm = fpDrawNormalCI,
    boxsize = 0.5,
    line.margin = 0.1,
    xlab = "Estimated Time from A+",
    ylab = "",
    txt_gp = fpTxtGp(
      xlab = gpar(fontsize = 18),
      ticks = gpar(fontsize = 16),
      title = gpar(fontsize = 14, fontface = "bold")
    )
  ) |>
  fp_set_style(
    box = c(
      rep("#133331", length(Time_Vec)), 
      rep("#aadedb", length(Time_Vec)), 
      rep("gold", length(Time_Vec)),
      rep("#440226", length(Time_Vec)),
      rep("#f62797", length(Time_Vec)), 
      rep("#FB9FD1", length(Time_Vec)),
      rep("#1F77B4", length(Time_Vec)),  # new color 7
      rep("#1F77B4", length(Time_Vec)),  # new color 7
      
      rep("#1F77B4", length(Time_Vec)),  # new color 7
      rep("#FF7F0E", length(Time_Vec)),  # new color 8
      rep("#2CA02C", length(Time_Vec))   # new color 9
    ) |> 
      lapply(function(x) gpar(fill = x, col = "#555555")),
    default = gpar(vertices = TRUE)
  ) |>
  fp_set_zebra_style("#F5F9F9") 

graph2ppt(file = "./Figures/ForestPlot_Time_CIs.pptx", width = 5, height = 9.4)


