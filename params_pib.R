#These parameters were extracted from real world DIAN data when available
#And literature when it was not


minMeaningfulZ <- 0.732
AposThresh <- 1.82
annual_rate <- 0.295
var_annual_rate <- 0.934
mean_within_above <- 0.08 # PIB: 8% ± 7% from Tolboom 2009 paper
# # 10.3% ± 1.25% for pTau217 alzpath from brum 2023
# # 16.7% ± 1.8% for pTau181 gothenberg simoa from brum 2023
sd_within_above <- 0.03 * mean_within_above
mean_within_below <- mean_within_above / 2  # PIB: 4.4% ± 4.2% from Tolboom paper
sd_within_below <- 0.03 * mean_within_below

avg_interval <- 2.02
interval_noise <- 1.06

xi <- -0.41
omega <- 1.00
alpha <- 2.87