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
source("./Fig1and2Funcs.R")

df <- read.csv("../Alamar_AmyloidClocks/Data/cleaned_df_withBatchNorming_20260420.csv")


cols <- c(which(names(df) == "fam_mutation") : which(names(df) == "plasmapTau181_jucker" ),
          which(names(df) == "Study.Period") : which(names(df) == "BD.pTau.217"))

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

df_flags <- merge(df[!duplicated(df$newid18), c("newid18", "VISITAGEc", "SEX", "RACE", "apoe",
                                                "Mutation", "fam_mutation", "DIAN_EYO")], df_flags, 
                  by = "newid18")

df_flags$SEX <- as.factor(df_flags$SEX)  
levels(df_flags$SEX) <- c("male", "female")
df_flags$fam_mutation <-as.factor(df_flags$fam_mutation)
levels(df_flags$fam_mutation) <- c("PSEN1", "PSEN2", "APP")

df_flags <- df_flags[df_flags$CL_two_plus == TRUE |df_flags$CSFpT181_Nico_corrected_any == TRUE |
                       df_flags$LUMIPULSE_CSF_pTau_two_plus == TRUE | df_flags$plasmapTau217_Nico_two_plus == TRUE |
                       df_flags$plasmapTau181_Nico_two_plus == TRUE | df_flags$plasmapTau181_jucker_two_plus |
                       df_flags$CSFpT217_Nico_corrected_two_plus == TRUE | df_flags$Alamar_pTau181_two_plus |
                       df_flags$Alamar_pTau217_two_plus == TRUE | df_flags$AlamarCSF_pTau217_two_plus == TRUE |
                       df_flags$AlamarCSF_pTau181_two_plus == TRUE | df_flags$BD.pTau.217_two_plus == TRUE,]

myVars <- c("VISITAGEc", "DIAN_EYO", "N", "SEX", "RACE", "apoe", "fam_mutation", 
            "CL_two_plus", "CSFpT181_Nico_corrected_any", "LUMIPULSE_CSF_pTau_two_plus", 
            "plasmapTau217_Nico_two_plus", "plasmapTau181_jucker_two_plus",
            "CSFpT181_Nico_corrected_two_plus", "CSFpT217_Nico_corrected_two_plus",
            "Alamar_pTau181_two_plus","Alamar_pTau217_two_plus",
            "AlamarCSF_pTau217_two_plus", "AlamarCSF_pTau181_two_plus",
            "BD.pTau.217_any")
factorVars <- c("SEX", "RACE", "apoe", "fam_mutation")

tab <- CreateTableOne(vars = myVars, strata = "Mutation", factorVars = factorVars, data = df_flags)

## Convert TableOne → data frame
tab_df <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

## Write to Excel
write.csv(tab_df, file = "./Figures/TableOne.csv")
