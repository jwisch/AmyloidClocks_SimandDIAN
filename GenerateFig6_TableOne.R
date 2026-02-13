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

source("./functions.R")

df <- read.csv("./Data/cleaned_df_withBatchNorming_20260106.csv")


cols <- which(names(df) == "CL") : which(names(df) == "AlamarCSF_pTau217")

ggpairs(df[!duplicated(df$newid18), cols])

graph2ppt(file = "./Figures/ggpairPlotAllBaselineData.pptx", width = 9, height = 9)


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


 result_PET_CSFpT217_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_CSFpT217_Nico_obj[[3]]$optimal_cutpoint - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$CSFpT217_Nico, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$CSFpT217_Nico, na.rm = TRUE),
                                                              "newid18", "TimefromBaseline",
                                                              "CSFpT217_Nico_Z", num_bootstraps = 1000)
 result_PET_plasmapTau217_Nico_Z <- bootstrap_get_Time_to_Positivity_amended(data.frame(df), (cp_PET_plasmapTau217_Nico_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$plasmapTau217_Nico, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$plasmapTau217_Nico, na.rm = TRUE),
                                                                   "newid18", "TimefromBaseline",
                                                                   "plasmapTau217_Nico_Z", num_bootstraps = 1000)
 result_PET_Alamar_pTau217_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_Alamar_pTau217_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 1 & !duplicated(df$newid18) & df$DIAN_EYO < -10,]$Alamar_pTau217, na.rm = TRUE))/sd(df[df$Mutation == 1 & !duplicated(df$newid18) & df$DIAN_EYO < -10,]$Alamar_pTau217, na.rm = TRUE),
                                                               "newid18", "TimefromBaseline",
                                                               "Alamar_pTau217_Z", num_bootstraps = 1000)
 result_PET_AlamarCSF_pTau217_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_AlamarCSF_pTau217_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$Alamar_pTau217, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$AlamarCSF_pTau217, na.rm = TRUE),
                                                                 "newid18", "TimefromBaseline",
                                                                 "AlamarCSF_pTau217_Z", num_bootstraps = 1000)
 result_pib_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (18 - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$CL, na.rm = TRUE),
                                                "newid18", "TimefromBaseline",
                                                "CL_Z", num_bootstraps = 1000)
 result_PET_CSFpT181_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_CSFpT181_Nico_obj[[3]]$optimal_cutpoint - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$CSFpT181_Nico, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$CSFpT181_Nico, na.rm = TRUE),
                                                              "newid18", "TimefromBaseline",
                                                              "CSFpT181_Nico_Z", num_bootstraps = 1000)
 result_PET_plasmapTau181_Nico_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_plasmapTau181_Nico_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$plasmapTau181_Nico, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$plasmapTau181_Nico, na.rm = TRUE),
                                                                   "newid18", "TimefromBaseline",
                                                                   "plasmapTau181_Nico_Z", num_bootstraps = 1000)
 result_PET_plasmapTau181_Jucker_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_plasmapTau181_Jucker_obj[[3]]$optimal_cutpoint - mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$plasmapTau181_jucker, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$plasmapTau181_jucker, na.rm = TRUE),
                                                                     "newid18", "TimefromBaseline",
                                                                     "plasmapTau181_jucker_Z", num_bootstraps = 1000)
 result_PET_Alamar_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_Alamar_pTau181_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 1 & !duplicated(df$newid18) & df$DIAN_EYO < -10,]$Alamar_pTau181, na.rm = TRUE))/sd(df[df$Mutation == 1 & !duplicated(df$newid18) & df$DIAN_EYO < -10,]$Alamar_pTau181, na.rm = TRUE),
                                                             "newid18", "TimefromBaseline",
                                                               "Alamar_pTau181_Z", num_bootstraps = 1000)
 result_PET_AlamarCSF_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_AlamarCSF_pTau181_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$AlamarCSF_pTau181, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$AlamarCSF_pTau181, na.rm = TRUE),
                                                             "newid18", "TimefromBaseline",
                                                               "AlamarCSF_pTau181_Z", num_bootstraps = 1000)
 result_PET_LumipulseCSF_pTau181_Z <- bootstrap_get_Time_to_Positivity(data.frame(df), (cp_PET_LumipulseCSF_pTau181_obj[[3]]$optimal_cutpoint- mean(df[df$Mutation == 0 & !duplicated(df$newid18),]$LUMIPULSE_CSF_pTau, na.rm = TRUE))/sd(df[df$Mutation == 0 & !duplicated(df$newid18),]$LUMIPULSE_CSF_pTau, na.rm = TRUE),
                                                             "newid18", "TimefromBaseline",
                                                               "LUMIPULSE_CSF_pTau_Z", num_bootstraps = 1000)

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
df$TimefromApos_plasmapT181_Jucker_Z <- approx(y = result_PET_plasmapTau181_Jucker_Z$Time_to_Positivity, 
                                             x = result_PET_plasmapTau181_Jucker_Z$Estimate, xout = df$plasmapTau181_jucker_Z)$y
df$TimefromApos_AlamarpT181_Z <- approx(y = result_PET_Alamar_pTau181_Z$Time_to_Positivity, 
                                      x = result_PET_Alamar_pTau181_Z$Estimate, xout = df$Alamar_pTau181_Z)$y
df$TimefromApos_AlamarCSFpT181_Z <- approx(y = result_PET_AlamarCSF_pTau181_Z$Time_to_Positivity, 
                                        x = result_PET_AlamarCSF_pTau181_Z$Estimate, xout = df$AlamarCSF_pTau181_Z)$y
df$TimefromApos_LumipulseCSFpT181_Z <- approx(y = result_PET_LumipulseCSF_pTau181_Z$Time_to_Positivity, 
                                           x = result_PET_LumipulseCSF_pTau181_Z$Estimate, xout = df$LUMIPULSE_CSF_pTau_Z)$y
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



layout_matrix <- rbind(c(1, 1, 1, 1, 1, 1, 1, 1, 1), c(2, 2, 3, 3, 4,4, 5, 5, 5), c(6,6,6, 7,7,7, 8,8,8), c(9,9,9,10,10, 10,11,11, 11))
grid.arrange(p_pib,
              p_CSF_217, p_AlamarCSF_217, p_Alamar_217, p_plasma_217,
             p_CSF_181, p_AlamarCSF_181, p_Alamar_181,
             p_LumipulseCSF_181, p_plasma_181, p_plasma_181_jucker,
              
             layout_matrix = layout_matrix)
graph2ppt(file = "./Figures/BootstrappedAmyloidTimeResults_20260106.pptx", width = 10, height = 9)

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
  
  return(p)
}





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



lemon::grid_arrange_shared_legend(p_error_pib,
                                  p_error_CSFpT217_Nico,
                                  p_error_AlamarCSF_pTau217, p_error_Alamar_pTau217, p_error_plasmapTau217_Nico,
                                  p_error_CSFpT181_Nico, p_error_AlamarCSF_pTau181, p_error_Alamar_pTau181,
                                  p_error_LumipulseCSF_pTau181,p_error_plasmapTau181_Nico, p_error_plasmapTau181_jucker,
                                  
                                  layout_matrix = layout_matrix)

graph2ppt(file = "./Figures/ObservedvsPredictedAgeatConversion.pptx", width = 10, height = 9)







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
             p_plasma_181_jucker,p_error_plasmapTau181_jucker, ncol = 2)
library(cowplot)
library(gridExtra)
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
  strip_legend(p_plasma_181_jucker),  strip_legend(p_error_plasmapTau181_jucker)
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
graph2ppt(file = "./Figures/Evaluation_of_Real_Models.pptx", width = 4.5, height = 10)




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
