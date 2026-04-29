setDT(av45)
setDT(ala_corrected)
setDT(cdr)

# ensure Date type
ala_corrected[, EXAMDATE := as.Date(EXAMDATE)]
av45_u <- unique(av45, by = c("PTID", "SCANDATE"))
ala_u  <- unique(ala_corrected, by = c("PTID", "EXAMDATE"))
cdr_u  <- unique(cdr, by = c("PTID", "VISDATE"))

# candidate matches
pairs_ala_av45 <- av45_u[
  ala_u,
  on = .(PTID),
  allow.cartesian = TRUE
][
  , date_diff := abs(as.numeric(SCANDATE - EXAMDATE))
][
  date_diff <= 365
]

pairs_ala_av45 <- pairs_ala_av45[
  order(date_diff)
][
  !duplicated(paste(PTID, EXAMDATE))
]

pairs_ala_av45 <- pairs_ala_av45[
  !duplicated(paste(PTID, SCANDATE))
]

ala_av45 <- merge(
  pairs_ala_av45,
  ala_u,
  by = c("PTID", "EXAMDATE"),
  all.x = TRUE
)

pairs_all <- cdr_u[
  ala_av45,
  on = .(PTID),
  allow.cartesian = TRUE
][
  , date_diff_cdr := abs(as.numeric(VISDATE - EXAMDATE))
][
  date_diff_cdr <= 365
]

pairs_all <- pairs_all[
  order(date_diff_cdr)
][
  !duplicated(paste(PTID, EXAMDATE))
]

pairs_all <- pairs_all[
  !duplicated(paste(PTID, VISDATE))
]


final_df <- merge(
  pairs_all,
  cdr_u,
  by = c("PTID", "VISDATE"),
  all.x = TRUE
)

ala <- final_df[,c("PTID", "VISDATE", "EXAMDATE", "SCANDATE", "TimefromBaseline", 
                  "Age", "CDGLOBAL.x", "CENTILOIDS.combat", "Alamar_pT181.x", "Alamar_pT217.x")]

colnames(ala)[7:10] <- c("CDGLOBAL", "CL.combat", "Alamar_pT181", "Alamar_pT217")

rm(ala_av45, ala_corrected, ala_pet, ala_u, av45, av45_sub, av45_u,
   cdr_u, pairs_all, final_df)
