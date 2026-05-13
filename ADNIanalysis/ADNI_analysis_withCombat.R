library(data.table)
library(lubridate)
library(sva)
library(tableone)
source("./functions.R")
source("./Simulation_functions.R")
source("./Fig1and2Funcs.R")
source("./functions_DistributionSpecificity.R")
df <- read.csv(".././ADNI_Feb2026/PTDEMOG_24Feb2026.csv")

df <- df[, c("PTID", "VISCODE", "VISDATE", "PTDOBYY", "PTGENDER")]
df <- df[!(df$PTID) %like% "381_S_",]


# alamar <- read.csv(".././ADNI_Feb2026/BSHRI_PLA_CSF_NULISA_CNS_24Feb2026.csv")
# alamar <- alamar[alamar$Target %in% c("pTau-217", "pTau-181"),]
# alamar <- alamar[alamar$SampleMatrixType == "CSF",]
# longitudinal_ala <- setDT(alamar)[, .N, by = "PTID"]
# longitudinal_ala <- longitudinal_ala[longitudinal_ala$N > 1,]
# 
# 
# 
# alamar_181 <- alamar[alamar$Target == "pTau-181", c("PTID", "EXAMDATE", "VISCODE", "PlateID", "NPQ")]
# alamar_217 <- alamar[alamar$Target == "pTau-217", c("PTID", "VISCODE", "NPQ")]
# 
# ala <- merge(alamar_181, alamar_217, by = c("PTID", "VISCODE"))
# colnames(ala)[5:6] <- c("pT181", "pT217") 
# 
# ##COMBAT ON ALAMAR DATA
# # 1. Extract numeric biomarker matrix
# dat <- ala %>%
#   select(pT181, pT217) %>%
#   as.matrix()
# 
# # 2. Transpose: features x samples
# dat_t <- t(dat)
# 
# # 3. Define batch variable
# batch <- factor(ala$PlateID)
# 
# # 4. Run ComBat
# combat_corrected <- ComBat(
#   dat = dat_t,
#   batch = batch,
#   par.prior = TRUE,
#   prior.plots = FALSE
# )
# 
# # 5. Transpose back to original shape
# corrected_df <- t(combat_corrected) %>% as.data.frame()
# 
# # 6. Reattach metadata
# ala_corrected <- ala %>%
#   select(PTID, VISCODE, EXAMDATE, PlateID) %>%
#   bind_cols(corrected_df)
# 
# ala_corrected <- ala_corrected[!duplicated(ala_corrected),]
# 
# 
# rm(ala, alamar, alamar_181, alamar_217, combat_corrected, corrected_df, dat, dat_t)
# 
# colnames(ala_corrected)[5:6] <- c("Alamar_pT181", "Alamar_pT217")
apoe <- read.csv(".././ADNI_Feb2026/APOERES_24Feb2026.csv")
apoe <- apoe[, c("PTID", "GENOTYPE")]
apoe <- apoe[!duplicated(apoe),]

mmse <- read.csv(".././ADNI_Feb2026/MMSE_04Mar2026.csv")
mmse <- mmse[, c("PTID", "VISDATE", "MMSCORE")]
mmse <- mmse[!(mmse$PTID) %like% "381_S_",]

cdr <- read.csv(".././ADNI_Feb2026/CDR_04Mar2026.csv")
cdr <- cdr[, c("PTID", "VISDATE", "CDGLOBAL", "CDRSB")]
cdr <- cdr[!(cdr$PTID) %like% "381_S_",]

av45 <- read.csv(".././ADNI_Feb2026/HarmonizedADNI.csv")

filter_val <- mean(av45[!duplicated(av45$PTID),]$CENTILOIDS.combat, na.rm = TRUE) + 5 * sd(av45[!duplicated(av45$PTID),]$CENTILOIDS.combat, na.rm = TRUE)
av45 <- av45[av45$CENTILOIDS.combat < filter_val,]

c2n <- read.csv(".././ADNI_Feb2026/UPENN_PLASMA_FUJIREBIO_QUANTERIX_24Feb2026.csv")
c2n <- c2n[, c("PTID", "EXAMDATE", "pT217_AB42_F")]
c2n <- c2n[c2n$pT217_AB42_F > 0,]
c2n <- c2n[!(c2n$PTID) %like% "381_S_",]


filter_val <- mean(c2n[!duplicated(c2n$PTID),]$pT217_AB42_F, na.rm = TRUE) + 5 * sd(c2n[!duplicated(c2n$PTID),]$pT217_AB42_F, na.rm = TRUE)
c2n <- c2n[c2n$pT217_AB42_F < filter_val,]

df$VISDATE <- as.Date(df$VISDATE, format = "%Y-%m-%d")
df$PTDOBYY <- as.Date(df$PTDOBYY, format = "%Y-%m-%d")
df$age_at_visit <- round(as.numeric(df$VISDATE - df$PTDOBYY) / 365.25, 1)

av45 <- merge(av45, df[, c("PTID", "PTDOBYY")], by = "PTID",
              all.x = TRUE, all.y = FALSE)
c2n <- merge(c2n, df[, c("PTID", "PTDOBYY")], by = "PTID",
             all.x = TRUE, all.y = FALSE)

av45$SCANDATE <- as.Date(av45$SCANDATE, format = "%Y-%m-%d")

setDT(df)
setDT(av45)

# Keep only common PTIDs for efficiency
common_ids <- intersect(df$PTID, av45$PTID)
df_sub <- df[PTID %in% common_ids]
av45_sub <- av45[PTID %in% common_ids]

# Merge by PTID to get all possible matches for each VISDATE
merged_all <- av45_sub[df_sub, on = .(PTID), allow.cartesian = TRUE]

# Compute difference in days between VISDATE and SCANDATE
merged_all[, diff_days := abs(as.numeric(VISDATE - SCANDATE))]

# Keep only matches within 2 years (~730 days)
merged_all <- merged_all[diff_days <= 730]

# For each df row, keep only the nearest SCANDATE
merged_nearest <- merged_all[
  , .SD[which.min(diff_days)], by = .(PTID, VISDATE)
]

# Optional: remove helper column
merged_nearest[, diff_days := NULL]

# Convert to data.table
setDT(merged_nearest)
setDT(c2n)

# Convert EXAMDATE to Date if it's not already
c2n[, EXAMDATE := as.Date(EXAMDATE)]

# Keep only common PTIDs
common_ids <- intersect(merged_nearest$PTID, c2n$PTID)
merged_sub <- merged_nearest[PTID %in% common_ids]
c2n_sub <- c2n[PTID %in% common_ids]

# Cartesian join by PTID
merged_all2 <- c2n_sub[merged_sub, on = .(PTID), allow.cartesian = TRUE]

# Compute difference in days
merged_all2[, diff_days := abs(as.numeric(VISDATE - EXAMDATE))]

# Keep only matches within 2 years (~730 days)
merged_all2 <- merged_all2[diff_days <= 730]

# For each merged_nearest row, keep only the nearest EXAMDATE
final_merged <- merged_all2[, .SD[which.min(diff_days)], by = .(PTID, VISDATE)]

# Optional: remove helper column
final_merged[, diff_days := NULL]


# Calculate TimefromBaseline_PET in years
final_merged[, TimefromBaseline_plasma := as.numeric(EXAMDATE - min(EXAMDATE)) / 365.25, by = PTID]
final_merged[, TimefromBaseline_PET := as.numeric(SCANDATE - min(SCANDATE)) / 365.25, by = PTID]

mmse$VISDATE <- as.Date(mmse$VISDATE, format = "%Y-%m-%d")
cdr$VISDATE <- as.Date(cdr$VISDATE, format = "%Y-%m-%d")

matched <- merge(final_merged, mmse, by = c("PTID", "VISDATE"),
                 all = FALSE)
matched <- merge(matched, cdr, by = c("PTID", "VISDATE"),
                 all = FALSE)

matched$Apos <- ifelse(matched$CENTILOIDS.combat > 18, 1, 0)
counts <- setDT(matched)[, .N, by = PTID]
counts <- counts[counts$N > 1,]

tableone <- merge(matched[, c("PTID", "Age", "PTGENDER", "CENTILOIDS.combat", "pT217_AB42_F", "Apos", "CDGLOBAL")], counts, by = "PTID", all = FALSE)
tableone <- tableone[!duplicated(tableone$PTID),]

vars <- c("Age", "PTGENDER", "CENTILOIDS.combat", "pT217_AB42_F", "Apos", "CDGLOBAL")
factorVars <- c("PTGENDER", "Apos", "CDGLOBAL")
CreateTableOne(data = tableone, vars = vars, factorVars = factorVars)


# source("./ADNIanalysis/cleanADNIAlamardata.R")

mean_CL <- mean(matched[matched$CDGLOBAL == 0 & matched$TimefromBaseline == 0,]$CENTILOIDS.combat, na.rm = TRUE)
sd_CL <- sd(matched[matched$CDGLOBAL == 0 & matched$TimefromBaseline == 0,]$CENTILOIDS.combat, na.rm = TRUE)

mean_pT <- mean(matched[matched$CDGLOBAL == 0 & matched$TimefromBaseline == 0,]$pT217_AB42_F, na.rm = TRUE)
sd_pT <- sd(matched[matched$CDGLOBAL == 0 & matched$TimefromBaseline == 0,]$pT217_AB42_F, na.rm = TRUE)

# mean_Alamar_pT181 <- mean(ala[ala$CDGLOBAL == 0 & ala$TimefromBaseline == 0,]$Alamar_pT181, na.rm = TRUE)
# sd_Alamar_pT181 <- sd(ala[ala$CDGLOBAL == 0 & ala$TimefromBaseline == 0,]$Alamar_pT181, na.rm = TRUE)
# 
# mean_Alamar_pT217 <- mean(ala[ala$CDGLOBAL == 0 & ala$TimefromBaseline == 0,]$Alamar_pT217, na.rm = TRUE)
# sd_Alamar_pT217 <- sd(ala[ala$CDGLOBAL == 0 & ala$TimefromBaseline == 0,]$Alamar_pT217, na.rm = TRUE)

matched$CL_Z <- (matched$CENTILOIDS.combat - mean_CL) / sd_CL
matched$pTau217_Z <- (matched$pT217_AB42_F - mean_pT) / sd_pT

# ala$CL_Z <- (ala$CL.combat - mean_CL) / sd_CL
# ala$Alamar_pT181_Z <- (ala$Alamar_pT181 - mean_Alamar_pT181) / sd_Alamar_pT181
# ala$Alamar_pT217_Z <- (ala$Alamar_pT217 - mean_Alamar_pT217) / sd_Alamar_pT217


Apos_thresh_CL_Z <- approx(matched$CENTILOIDS.combat, matched$CL_Z, 18)$y
matched$Apos <- ifelse(matched$CL_Z > Apos_thresh_CL_Z, 1, 0)
# ala$Apos <- ifelse(ala$CL_Z > Apos_thresh_CL_Z, 1, 0)

###################
##PICK UP HERE AND RUN ALAMAR STUFF
###################
pT217_Z_PET_Thresh <- run_cutpointr_analysis(matched, "pTau217_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
                                             "Centiloids, Amyloid PET (Z)", "Plasma pTau217 / AB42 (Z)")
# Alamar_pT181_Z_PET_Thresh <- run_cutpointr_analysis(ala, "Alamar_pT181_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
#                                                     "Centiloids, Amyloid PET (Z)", "Alamar Plasma pTau181 (Z)")
# Alamar_pT217_Z_PET_Thresh <- run_cutpointr_analysis(ala, "Alamar_pT217_Z", "Apos", y1 = -3, y2 = 7, "CL_Z",
#                                              "Centiloids, Amyloid PET (Z)", "Alamar Plasma pTau217 (Z)")

RofC_PET_Z <- get_RofC_df(matched, "PTID", "CL_Z", "TimefromBaseline")
RofC_pT217_Z <- get_RofC_df(matched, "PTID", "pTau217_Z", "TimefromBaseline")

# longitudinal_ala <- setDT(ala)[, .N, by = "PTID"]
# longitudinal_ala <- longitudinal_ala[longitudinal_ala$N > 1,]
# 
# RofC_Alamar_pT181_Z <- get_RofC_df(data.frame(ala), "PTID", "Alamar_pT181_Z", "TimefromBaseline")
# RofC_Alamar_pT217_Z <- get_RofC_df(ala, "PTID", "Alamar_pT217_Z", "TimefromBaseline")



rel_accum_PET_Z <- mean(max(RofC_PET_Z[RofC_PET_Z$classification == 1,]$rate_of_change),
                        min(RofC_PET_Z[RofC_PET_Z$classification == 2,]$rate_of_change))

rel_accum_pT217_Z <- mean(max(RofC_pT217_Z[RofC_pT217_Z$classification == 1,]$rate_of_change),
                          min(RofC_pT217_Z[RofC_pT217_Z$classification == 2,]$rate_of_change))

colnames(RofC_PET_Z)[2] <- "RofC_PET_Z"
colnames(RofC_pT217_Z)[2] <- "RofC_pT217_Z"


RofC_PET_Z <- combine_RofC_and_df(matched, "CL_Z", RofC_PET_Z,"RofC_PET_Z", rel_accum_PET_Z, "PTID", "Age", "PTGENDER")
RofC_pT217_Z <- combine_RofC_and_df(matched, "pTau217_Z", RofC_pT217_Z,"RofC_pT217_Z", rel_accum_pT217_Z, "PTID", "Age", "PTGENDER")


cp_PET_Z <- cutpointr(RofC_PET_Z, CL_Z, ReliableAccumulator, 
                      method = maximize_metric, metric = youden, na.rm = TRUE)
cp_pT217_Z <- cutpointr(RofC_pT217_Z, pTau217_Z, ReliableAccumulator, 
                      method = maximize_metric, metric = youden, na.rm = TRUE)


PET_Z_RelAccum <- run_cutpointr_analysis(RofC_PET_Z, "CL_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                         "CL_Z",
                                         "Amyloid PET (Z)", "ARC Amyloid PET (Z)")
pT217_Z_RelAccum <- run_cutpointr_analysis(RofC_pT217_Z, "pTau217_Z", "ReliableAccumulator", y1 = -3, y2 = 7, 
                                         "pTau217_Z",
                                         "Plasma pTau217/AB42 (Z)", "ARC Plasma pTau217/AB42(Z)")

matched <- define_baseline_rel_to_RelAccum(matched, "CL_Z", cp_PET_Z, id_col = "PTID", age_col = "Age")
matched <- define_baseline_rel_to_RelAccum(matched, "pTau217_Z", cp_pT217_Z, "ReliableAccumThresh_pTau217_Z", id_col = "PTID", age_col = "Age")

RofC_list <- list(
  RofC_PET_Z,
  RofC_pT217_Z
)

RofC_cols <- c(
  "RofC_PET_Z",
  "RofC_pT217_Z"
)


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


RofC_PET_Z <- merge(RofC_PET_Z[, c("PTID", "RofC_PET_Z")], 
                    matched[matched$TimefromBaseline == 0 & !is.na(matched$CL_Z),
                       c("PTID", "Age", "PTGENDER", 
                         "CL_Z")], by = "PTID", all = FALSE)

RofC_pT217_Z <- merge(RofC_pT217_Z[, c("PTID", "RofC_pT217_Z")], 
                    matched[matched$TimefromBaseline == 0 & !is.na(matched$pTau217_Z),
                            c("PTID", "Age", "PTGENDER", 
                              "pTau217_Z")], by = "PTID", all = FALSE)


RofC_PET_Z <- get_relAccumFlag(RofC_PET_Z, "RofC_PET_Z",
                               rel_accum_PET_Z)
RofC_pT217_Z <- get_relAccumFlag(RofC_pT217_Z, "RofC_pT217_Z",
                               rel_accum_pT217_Z)


cp_PET_Z <- cutpointr(RofC_PET_Z, CL_Z, ReliableAccumulator, 
                      method = maximize_metric, metric = youden, na.rm = TRUE)
cp_pT217_Z <- cutpointr(RofC_pT217_Z, pTau217_Z, ReliableAccumulator, 
                      method = maximize_metric, metric = youden, na.rm = TRUE)


cp_PET_pT217_obj <- get_RofC_scatter_and_ROC(matched[!duplicated(matched$PTID),],
                                             pTau217_Z, Apos, 
                                                          "Plasma pTau217/AB42 (Z)", "Amyloid Positive - 18 CL")


RofC_fits$modality <- c("PET", "Plasma")
RofC_fits$Analyte <- c("PET", "pT217")
RofC_fits$Assay <- c( "PET", "Fujirubio")

p2b <- ggplot(RofC_fits, aes(x = omega, y = xi, label = biomarker_col, colour = Assay,
                      shape = Analyte)) +
  geom_point(aes(size = alpha, alpha = 0.6)) + scale_colour_manual(values = c("#007B82", "#C77CFF",  "green", "black","#F8766D" )) +
  scale_shape_manual(values = c(19, 17,19)) +
  # geom_text(vjust = 1.3, hjust = -0.1) +
  facet_wrap(~modality) +
  theme_bw() +
  labs(y = "xi", x = "omega") +
  theme(legend.position = "bottom") + 
  xlim(c(0, 7.5)) + ylim(c(-4, 0.6))



# result_PET_Z <- bootstrap_get_Time_to_Positivity(data.frame(matched), Apos_thresh_CL_Z,
#                                                  "PTID", "TimefromBaseline",
#                                                  "CL_Z", num_bootstraps = 1000)
# 
# result_pT217_Z <- bootstrap_get_Time_to_Positivity(data.frame(matched), cp_PET_pT217_obj[[3]]$optimal_cutpoint,
#                                                  "PTID", "TimefromBaseline",
#                                                  "pTau217_Z", num_bootstraps = 1000)
# 
# saveRDS(result_PET_Z, "./Data/result_PET_Z_ADNI.RDS")
# saveRDS(result_pT217_Z, "./Data/result_pT217_Z_ADNI.RDS")

result_PET_Z <- readRDS("./Data/result_PET_Z_ADNI.RDS")
result_pT217_Z <- readRDS("./Data/result_pT217_Z_ADNI.RDS")


matched$TimefromApos_Z <- approx(y = result_PET_Z$Time_to_Positivity, x = result_PET_Z$Estimate, xout = matched$CL_Z)$y
matched$TimefromApos_pT217_Z <- approx(y = result_pT217_Z$Time_to_Positivity, x = result_pT217_Z$Estimate, xout = matched$pTau217_Z)$y


p_PET  <- get_TimefromAposPlot(result_PET_Z, cp_PET_Z$optimal_cutpoint, 
                               "Amyloid PET -  Z")
p_pT217  <- get_TimefromAposPlot(result_pT217_Z, cp_pT217_Z$optimal_cutpoint, 
                               "Plasma pT217/AB42 Z")

grid.arrange(p_PET, p_pT217, nrow = 2, ncol = 1)



p_error_PET <- plot_threshold_crossings(data.frame(matched), id_col = "PTID", age_col = "Age", 
                                        value_col = "CL_Z", threshold = Apos_thresh_CL_Z, color_grad = NULL, 
                                        XLIM = c(50, 90), YLIM = c(50, 90),
                                        XPOS = 70, YPOS = 65)
p_error_pT217 <- plot_threshold_crossings(matched, id_col = "PTID", age_col = "Age", 
                                          value_col = "pTau217_Z", threshold = Apos_thresh_CL_Z, color_grad = NULL)#, 
                                          # XLIM = c(50, 90), YLIM = c(50, 90),
                                          # XPOS = 70, YPOS = 65)


results <- data.frame("Biomarker" = c("PET", "PlasmapT217/AB42"),
                      "MAE" = c(p_error_PET[[2]],
                                p_error_pT217[[2]]),
                      "Modality" = c("PET", "Plasma"))



p3b <- results %>%
  mutate(Biomarker = forcats::fct_reorder(Biomarker, MAE, .fun = mean)) %>%
  ggplot(aes(x = Biomarker, y = MAE, shape = Modality, colour = Modality)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + geom_hline(yintercept = 2, linetype = "dashed") +
  scale_colour_manual(values = c("black", "#007B82", "#F8766D")) +
  scale_shape_manual(values = c(19, 3,17)) + xlab("") + ylab("Mean Average Error (Years)")


# result_PET_dist_Z <- bs_DistSpecplot(data.frame(matched), Apos_thresh_CL_Z,
#                               "PTID", "TimefromBaseline",
#                                "CL_Z", num_bootstraps = 1000)
# 
# result_pT217_dist_Z <- bs_DistSpecplot(data.frame(matched), cp_PET_pT217_obj[[3]]$optimal_cutpoint,
#                                        "PTID", "TimefromBaseline",
#                                       "pTau217_Z", num_bootstraps = 1000)
# 
# saveRDS(result_PET_dist_Z, "./Data/result_PET_dist_Z.RDS")
# saveRDS(result_pT217_dist_Z, "./Data/result_pT217_dist_Z.RDS")

result_PET_dist_Z <- readRDS("./Data/result_PET_dist_Z.RDS")
result_pT217_dist_Z <- readRDS("./Data/result_pT217_dist_Z.RDS")


p_4_PET <- plot_hist_density_at_time(data = result_PET_dist_Z, rel_accum = -20, 
                                time_vals = c(0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c("#003B44",   # darker lead-in
                                                "#007B82",   # 2nd (your first required color)
                                                "#00BFC4",   # 3rd (your second required color)
                                                "#6FE3E6",   # lighter continuation
                                                "#CFF7F8") ,  # very light endpoint
                                "Amyloid PET")


p_4_plasma <- plot_hist_density_at_time(data = result_pT217_dist_Z, rel_accum = -20, 
                          time_vals = c(0, 2, 4, 6),
                          binwidth = 0.01,
                          color_vec =   c(
                            "#A93E36",  # darkest
                            "#F8766D",  # 2nd darkest (your anchor)
                            "#FB9A93",
                            "#FDC1BD",
                            "#FFE5E3"   # lightest
                          ) ,  # very light endpoint
                          "Plasma pT217/AB42 [Fujirubio]")
layout_matrix <- rbind(c(1, 1), c(2, 2), c(3, 4))
grid.arrange(p2b, p3b, p_4_PET, p_4_plasma, layout_matrix = layout_matrix)
grid.arrange(p_PET, p_error_PET[[1]],
             p_pT217, p_error_pT217[[1]], nrow = 2, ncol = 2)











