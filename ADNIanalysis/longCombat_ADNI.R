library(longCombat)
library(invgamma)
library(lme4)
library(readxl)
library(data.table)
library(dplyr)

df <- read.csv(".././ADNI_Feb2026/PTDEMOG_24Feb2026.csv")

df <- df[, c("PTID", "VISCODE", "VISDATE", "PTDOBYY", "PTGENDER")]
df$VISDATE <- as.Date(df$VISDATE, format = "%Y-%m-%d")
df$PTDOBYY <- as.Date(df$PTDOBYY, format = "%Y-%m-%d")
df$age_at_visit <- round(as.numeric(df$VISDATE - df$PTDOBYY) / 365.25, 1)

cdr <- read.csv(".././ADNI_Feb2026/CDR_04Mar2026.csv")
cdr <- cdr[, c("PTID", "VISDATE", "VISCODE", "CDGLOBAL", "CDRSB")]


av45 <- read.csv(".././ADNI_Feb2026/UCBERKELEY_AMY_6MM_24Feb2026.csv")
av45 <- av45[, c("LONIUID", "PTID", "VISCODE", "SITEID", 
                 "SCANDATE", "TRACER", "AMYLOID_STATUS", "CENTILOIDS", "SUMMARY_SUVR")]

filter_val <- mean(av45[!duplicated(av45$PTID),]$CENTILOIDS, na.rm = TRUE) + 5 * sd(av45[!duplicated(av45$PTID),]$CENTILOIDS, na.rm = TRUE)
av45 <- av45[av45$CENTILOIDS < filter_val,]
av45 <- merge(av45, df[, c("PTID", "PTDOBYY")], by = "PTID",
              all.x = TRUE, all.y = FALSE)

av45$SCANDATE <- as.Date(av45$SCANDATE, format = "%Y-%m-%d")
av45$Age <- as.numeric(av45$SCANDATE - av45$PTDOBYY) / 365.25

av45 <- unique(av45, by = c("PTID", "VISCODE"))
av45 <- av45[!is.na(av45$Age),]

setDT(av45)[, TimefromBaseline := as.numeric(Age- min(Age)), by = PTID]

av45$SITEID <- as.factor(av45$SITEID)
cdr$VISDATE <- as.Date(cdr$VISDATE, format = "%Y-%m-%d")

#Matching av45 to nearest possible cdr
setDT(cdr)
setDT(av45)
setkey(cdr, PTID, VISDATE)
av45_cdr <- cdr[av45, on = .(PTID, VISDATE = SCANDATE), roll = "nearest"]
av45_cdr <- av45_cdr[!is.na(av45_cdr$CDGLOBAL),]

av45_cdr <- data.frame(av45_cdr)[!(av45_cdr$PTID %like% "381_S_"),]

av45_combat <- longCombat(idvar='PTID', 
           timevar='TimefromBaseline',
           batchvar='SITEID', 
           features=c("CENTILOIDS", "SUMMARY_SUVR"), 
           formula='Age + CDGLOBAL*TimefromBaseline',
           ranef='(1|PTID)',
           data=data.frame(av45_cdr[, c("PTID", "TimefromBaseline",
                                        "SITEID", "CENTILOIDS", 
                                        "SUMMARY_SUVR", "Age", 
                                        "CDGLOBAL")]))

av45_harmonized <- av45_combat$data_combat


av45 <- merge(av45[, c("PTID", "SITEID", "TimefromBaseline",
                       "Age",
                       "SCANDATE", "CENTILOIDS")],
              av45_harmonized[, c("PTID", "TimefromBaseline", "CENTILOIDS.combat")], 
              by = c("PTID",  "TimefromBaseline"))




library(ggplot2)
p1 <- ggplot(av45_harmonized[av45_harmonized$TimefromBaseline < 20,], aes(x = TimefromBaseline, y = CENTILOIDS.combat, group = PTID, colour = SITEID)) +
  geom_point(alpha = 0.8) + geom_line(alpha = 0.2) + theme(legend.position = "none") + theme_bw() + guides(colour = "none")

p2 <- ggplot(av45[av45$TimefromBaseline < 20,], aes(x = TimefromBaseline, y = CENTILOIDS, group = PTID, colour = SITEID)) +
  geom_point(alpha = 0.8) + geom_line(alpha = 0.2) + theme_bw() + guides(colour = "none")

gridExtra::grid.arrange(p1, p2, nrow = 1)


# # get the harmonized data
# simdata_harmonized <- simdata_combat$data_combat
# # save combat feature names
# featurenames.combat <- names(simdata_harmonized)[4:23]
# # merge with original dataframe
# simdata <- merge(simdata, simdata_harmonized[,c(1,2,4:23)], by=c('subid', 'time'))

ARC_raw <- av45[
  , if (.N >= 2) {
    list(annual_rate = coef(lm(CENTILOIDS ~ TimefromBaseline))[2])
  } else {
    list(annual_rate = NA_real_)
  },
  by = PTID
]


ARC_combat <- setDT(av45_harmonized)[
  , if (.N >= 2) {
    list(annual_rate = coef(lm(CENTILOIDS.combat ~ TimefromBaseline))[2])
  } else {
    list(annual_rate = NA_real_)
  },
  by = PTID
]

hist(ARC_raw$annual_rate)
hist(ARC_combat$annual_rate)

ARC_raw$annual_rate_Z <- (ARC_raw$annual_rate - mean(ARC_raw$annual_rate, na.rm = TRUE))/ sd(ARC_raw$annual_rate, na.rm = TRUE)
ARC_combat$annual_rate_Z <- (ARC_combat$annual_rate - mean(ARC_combat$annual_rate, na.rm = TRUE))/ sd(ARC_combat$annual_rate, na.rm = TRUE)

hist(ARC_raw$annual_rate_Z)
hist(ARC_combat$annual_rate_Z)




ARC_raw$annual_rate_Z <- (ARC_raw$annual_rate - mean(ARC_raw$annual_rate, na.rm = TRUE))/ sd(ARC_raw$annual_rate, na.rm = TRUE)
ARC_combat$annual_rate_Z <- (ARC_combat$annual_rate - mean(ARC_combat$annual_rate, na.rm = TRUE))/ sd(ARC_combat$annual_rate, na.rm = TRUE)

hist(ARC_raw$annual_rate_Z)
hist(ARC_combat$annual_rate_Z)


#THESE ARE THE OUTLIERS TO DROP
ARC_raw[abs(ARC_raw$annual_rate_Z) > 5,]
ARC_combat[abs(ARC_combat$annual_rate_Z) > 5,]

bad_ids <- ARC_raw[abs(annual_rate_Z) > 5, PTID]
av45[PTID %in% bad_ids, CENTILOIDS := NA]
bad_ids <- ARC_combat[abs(annual_rate_Z) > 5, PTID]
av45[PTID %in% bad_ids, CENTILOIDS.combat := NA]



write.csv(av45, "./processedData/HarmonizedADNI.csv", row.names = FALSE)
