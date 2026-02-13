#This script can only be run with DIAN datafreeze 18
#DIAN datafreeze 18 is only available with a data request

library(ggplot2)
library(GGally)
library(dplyr)
library(sva) #HAVE TO SWITCH TO R VERSION 4.5 TO MAKE THIS WORK
library(mclust)

source("./functions.R")

lumipulse <- read.csv(".././DIAN_DF18/df18_release/lumi_DCF.csv")
lumipulse <- lumipulse[, c("newid18", "visit", "LUMIPULSE_CSF_pTau")]
lumipulse <- lumipulse[!is.na(lumipulse$LUMIPULSE_CSF_pTau),]

ala_csf <- read.csv(".././DIAN_DF18/BMC_Alamar_CSF/combined_Alamar_after_qulificantion.csv")
ala_csf <- ala_csf[, c("newid18", "visit", "pTau.181", "pTau.217")]

ala_plasma <- read.csv(".././DIAN_DF18/BMC_Alamar_Plasma/Alamar_Plasma_after_QC.csv")
ala_plasma <- ala_plasma[, c("newid18", "visit", "pTau.181", "pTau.217")]

nico <- read.csv(".././DIAN_DF18/Project1&2/df18_pj2_mstau.csv")
nico <- nico[, c("newid18", "visit", "pT181.T181.", "pT217.T217.")]
colnames(nico)[3:4] <- c("CSFpT181_Nico", "CSFpT217_Nico")

plasma181 <- read.csv(".././DIAN_DF18/Project1&2/D2512_mj_newid18.csv")
plasma181 <- plasma181[, c("newid18", "visit", "Plasma_pTau__pg_ml_")]
colnames(plasma181)[3] <- c("plasmapTau181_jucker")

plasma_nico <- read.csv(".././DIAN_DF18/Project1&2/D2512_ptauMx_newid18_20260106.csv")
plasma_nico <- plasma_nico[, c("newid18", "visit", "X.p.tau217", "X.p.tau181")]
colnames(plasma_nico)[3:4] <- c("plasmapTau217_Nico", "plasmapTau181_Nico") 
plasma_nico <- plasma_nico[!plasma_nico$newid18 == "",]

#THIS IS ALL THE STANDARD DATA. I HAVE COMBINED IT INTO A SINGLE DATAFRAME
pib <- read.csv(".././DIAN_DF18/df18_release/pib_DCF.csv")
csf <-  read.csv(".././DIAN_DF18/df18_release/BIOMARKER_DCF.csv")

demogs <- read.csv(".././DIAN_DF18/df18_release/DEMOGRAPHICS_DCF.csv")
genetics <- read.csv(".././DIAN_DF18/df18_release/GENETIC_DCF.csv")

df <- demogs[, c("newid18", "imagid", "visit", "VISITAGEc", "SEX", "RACE")]
df <- merge(df, genetics[, c("newid18", "visit", "apoe", "Mutation", "fam_mutation")], by = c("newid18", "visit"), all = FALSE)

pib <- pib[, c("newid18", "visit", "PIB_fSUVR_rsf_TOT_CORTMEAN")]
pib <- pib[!is.na(pib$PIB_fSUVR_rsf_TOT_CORTMEAN),]
pib$CL <- pib$PIB_fSUVR_rsf_TOT_CORTMEAN * 45 - 47.5
pib <- pib[!is.na(pib$CL),]

df <- merge(df, pib[, c("newid18", "visit", "CL")], by = c("newid18", "visit"),
            all.x = TRUE, all.y = FALSE)



rm(demogs, genetics, pib, plasma, csf)

df <- merge(df, nico, by = c("newid18", "visit"), all.x = TRUE, all.y = FALSE)
rm(nico)

df <- merge(df, lumipulse, by = c("newid18", "visit"),
            all.x = TRUE, all.y = FALSE)
rm(lumipulse)
df <- merge(df, plasma_nico, by = c("newid18", "visit"),
            all.x = TRUE, all.y = FALSE)

rm(plasma_nico)
df <- merge(df, plasma181, by = c("newid18", "visit"),
            all.x = TRUE, all.y = FALSE)

rm(plasma181)

cols_to_convert <- c("CSFpT181_Nico", "CSFpT217_Nico", 
                     "plasmapTau181_Nico", "plasmapTau217_Nico")

df <- df %>%
  mutate(across(all_of(cols_to_convert), as.numeric))

###############################################################################
#Cleaning alamar data for batch effects from Wenjing
###############################################################################


run_batch_PCA_QC <- function(
    ala_plasma,
    col_start,
    col_end
) {
  library(dplyr)
  library(tibble)
  library(mclust)
  library(sva)
  
  # ------------------------------ #
  # 1. SELECT & CLEAN PROTEIN DATA
  # ------------------------------ #
  
  protein_data <- ala_plasma %>%
    dplyr::select(newid18, all_of(col_start):all_of(col_end)) %>%
    mutate(across(-newid18, as.numeric))
  
  protein_data <- protein_data[complete.cases(protein_data), ]
  
  # --------------------------------------------- #
  # 2. FIRST PCA → compute batch clusters (Mclust)
  # --------------------------------------------- #
  
  pca_obj <- prcomp(protein_data[, -1], scale. = TRUE)
  pca <- as.data.frame(pca_obj$x)
  pca$Participant_ID <- protein_data$newid18
  
  Vstdiv <- c(mean(pca$PC1) - 3 * sd(pca$PC1), mean(pca$PC1) + 3 * sd(pca$PC1))
  Hstdiv <- c(mean(pca$PC2) - 3 * sd(pca$PC2), mean(pca$PC2) + 3 * sd(pca$PC2))
  
  # Mclust classification
  fit <- Mclust(pca[, c("PC1", "PC2")], G = 2)
  pca$batch <- as.factor(fit$classification)
  
  # -------------------------- #
  # 3. PREPARE COMBAT MATRIX
  # -------------------------- #
  
  complete_rows <- ala_plasma %>%
    select(all_of(col_start):all_of(col_end)) %>%
    mutate(across(everything(), as.numeric)) %>%
    complete.cases()
  
  combat_data <- ala_plasma %>%
    select(all_of(col_start):all_of(col_end)) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(complete_rows) %>%
    as.matrix() %>%
    t()
  
  combat_data <- apply(combat_data, 2, as.numeric) %>%
    matrix(nrow = nrow(combat_data))
  
  batch <- factor(pca$batch)
  
  combat_corrected <- ComBat(
    dat = combat_data,
    batch = batch,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  # ---------------------------------------- #
  # 4. RECONSTRUCT CLEAN DATAFRAME AFTER COMBAT
  # ---------------------------------------- #
  
  protein_adjusted <- t(combat_corrected) %>% as.data.frame()
  colnames(protein_adjusted) <- names(ala_plasma)[match(col_start, names(ala_plasma)):
                                                    match(col_end, names(ala_plasma))]
  
  IDs <- ala_plasma$newid18[complete_rows]
  visits <- ala_plasma$visit[complete_rows]
  
  data_pca <- protein_adjusted %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Participant_ID")
  
  data_pca$Participant_ID <- IDs
  data_pca$visit <- visits
  
  # ---------------- #
  # 5. SECOND PCA QC
  # ---------------- #
  
  protein_clean <- data_pca %>%
    select(all_of(col_start):all_of(col_end)) %>%
    mutate(across(everything(), as.numeric))
  
  pca_obj2 <- prcomp(protein_clean, scale. = TRUE)
  pca2 <- as.data.frame(pca_obj2$x)
  
  pca2$Participant_ID <- data_pca$Participant_ID
  pca2$visit <- data_pca$visit
  
  Vstdiv2 <- c(mean(pca2$PC1) - 3 * sd(pca2$PC1), mean(pca2$PC1) + 3 * sd(pca2$PC1))
  Hstdiv2 <- c(mean(pca2$PC2) - 3 * sd(pca2$PC2), mean(pca2$PC2) + 3 * sd(pca2$PC2))
  
  # Remove outliers
  keep_ids2 <- pca2 %>%
    filter(
      PC1 >= Vstdiv2[1], PC1 <= Vstdiv2[2],
      PC2 >= Hstdiv2[1], PC2 <= Hstdiv2[2]
    ) %>%
    pull(Participant_ID)
  
  # ------------------------- #
  # 6. FINAL CLEANED DATAFRAME
  # ------------------------- #
  
  Ances_ds_after_QC <- data_pca %>%
    filter(Participant_ID %in% keep_ids2)
  
  colnames(Ances_ds_after_QC)[1] <- "newid18"
  
  Ances_ds_after_QC <- Ances_ds_after_QC[!duplicated(Ances_ds_after_QC),]
  
  return(Ances_ds_after_QC)
}
ala_plasma_combatted <- run_batch_PCA_QC(
  ala_plasma = ala_plasma,
  col_start  = "pTau.181",
  col_end    = "pTau.217"
)

ala_csf_combatted <- run_batch_PCA_QC(
  ala_plasma = ala_csf,
  col_start  = "pTau.181",
  col_end    = "pTau.217"
)
###############################################################################
###############################################################################
colnames(ala_plasma_combatted)[2:3] <- c("Alamar_pTau181", "Alamar_pTau217")
colnames(ala_csf_combatted)[2:3] <- c("AlamarCSF_pTau181", "AlamarCSF_pTau217")

df <- merge(df, ala_plasma_combatted,
            by = c("newid18", "visit"), all.x = TRUE, all.y = FALSE)

df <- merge(df, ala_csf_combatted,
            by = c("newid18", "visit"), all.x = TRUE, all.y = FALSE)

eyo <- read.csv(".././DIAN_DF18/df18_release/CLINICAL_DCF.csv")
eyo <- eyo[, c("newid18", "visit", "DIAN_EYO")]

df <- merge(df, eyo,
            by = c("newid18", "visit"), all.x = TRUE, all.y = FALSE)

#Dropping all values more than 5 sd from mean
df <- df %>%
  mutate(across(CL:AlamarCSF_pTau217, ~ {
    m  <- mean(.x, na.rm = TRUE)
    sd5 <- 5 * sd(.x, na.rm = TRUE)
    .x[.x < m - sd5 | .x > m + sd5] <- NA
    .x
  }))




write.csv(df, "./Data/cleaned_df_withBatchNorming_20260106.csv", row.names = FALSE)


