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
library(forcats)
library(emmeans) # for confidence intervals
library(export)
source("./functions.R")
source("./Simulation_functions.R")
source("./Fig1and2Funcs.R")
source("./ADNI_prep_for_analysis.R")

df_for_positivity <- readRDS("./Data/df_for_positivity_DIAN.RDS")

matched$Apos_by_plasma <- ifelse(matched$pTau217_Z > cp_PET_pT217_obj[[3]]$optimal_cutpoint, 1, 0)


get_nonConverters <- function(
    df,
    id_col,
    apos_col) {
  # Keep only IDs with >= 2 rows
  df_filtered <- df %>%
    group_by({{ id_col }}) %>%
    filter(n() >= 2) %>%
    ungroup()
# 1. Identify IDs that are strictly 0 (Never Positive)
zero_ids <- df %>%
  group_by({{ id_col }}) %>%
  filter(max({{ apos_col }}, na.rm = TRUE) == 0) %>% 
  pull({{ id_col }}) %>%
  unique()

# 2. Identify IDs that are strictly 1 (Always Positive in this window)
solo_ids <- df %>%
  group_by({{ id_col }}) %>%
  filter(min({{ apos_col }}, na.rm = TRUE) == 1) %>% 
  pull({{ id_col }}) %>%
  unique()
return(list(zero_ids, solo_ids))}



plot_nonConverter_histogram <- function(
    df,
    id_col,
    apos_col,
    time_col,
    timefrombaseline_col,
    plotTitle = NULL,
    seed = 123
) {
  
  #-----------------------------
  # Get non-converter IDs
  #-----------------------------
  nonConverter_obj <- get_nonConverters(
    df = df,
    id_col = {{ id_col }},
    apos_col = {{ apos_col }}
  )
  
  #-----------------------------
  # Compute TimefromApos_basedOnFirst
  #-----------------------------
  matched_firstApos <- df %>%
    group_by({{ id_col }}) %>%
    mutate(
      first_time_val = {{ time_col }}[
        which.min({{ timefrombaseline_col }})
      ],
      TimefromApos_basedOnFirst =
        first_time_val + {{ timefrombaseline_col }}
    ) %>%
    ungroup()
  #-----------------------------
  # Randomly select 1 row per PTID
  #-----------------------------
  set.seed(seed)
  
  zero_df <- matched_firstApos %>%
    filter(
      {{ id_col }} %in% nonConverter_obj[[1]],
      {{ timefrombaseline_col }} > 1
    ) %>%
    group_by({{ id_col }}) %>%
    slice_sample(n = 1) %>%
    ungroup()

  solo_df <- matched_firstApos %>%
    filter(
      {{ id_col }} %in% nonConverter_obj[[2]],
      {{ timefrombaseline_col }} > 1
    ) %>%
    group_by({{ id_col }}) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  #-----------------------------
  # Plot
  #-----------------------------
  p <- ggplot() +
    geom_histogram(
      data = zero_df,
      aes(x = TimefromApos_basedOnFirst),
      fill = "#E07529",
      alpha = 0.4,
      position = "identity"
    ) +
    geom_histogram(
      data = solo_df,
      aes(x = TimefromApos_basedOnFirst),
      fill = "#7F7991",
      alpha = 0.4,
      position = "identity"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    theme_bw() +
    labs(
      x = "Forecasted Time from Positivity",
      y = "Count",
      title = plotTitle
    )
  
  return(p)
}

plot_nonConverter_histogram(
  df = matched,
  id_col = PTID,
  apos_col = Apos,
  time_col = TimefromApos_Z,
  timefrombaseline_col = TimefromBaseline,
  plotTitle = "A. Amyloid PET, ADNI"
)

plot_nonConverter_histogram(
  df = matched,
  id_col = PTID,
  apos_col = Apos_by_plasma,
  time_col = TimefromApos_pT217_Z,
  timefrombaseline_col = TimefromBaseline,
  plotTitle = "B. Plasma pTau217/AB42, ADNI"
)


plot_nonConverter_histogram(
  df = df_for_positivity,
  id_col = newid18,
  apos_col = Apos_by_CL_Z,
  time_col = TimefromApos_Z,
  timefrombaseline_col = TimefromBaseline,
  plotTitle = "C. Amyloid PET, DIAN"
)

plot_nonConverter_histogram(
  df = df_for_positivity,
  id_col = newid18,
  apos_col = Apos_by_CSFpT181_Nico_Z,
  time_col = TimefromApos_CSFpT181_Z,
  timefrombaseline_col = TimefromBaseline,
  plotTitle = "C. CSF pT181 - IPMS, DIAN"
)





