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
library(gganimate)
library(patchwork)
library(magick)
library(purrr)
source("./Simulation_functions.R")
source("./functions.R")
seed <- 10
#TODO: Fix the mean_within_above at a singl evalue and figure out how to alter histogram properties
mean_within_above <- 0.03 #fixed within individual variance

n_ids <- 300
num_of_scans_for_model <- 2
minMeaningfulZ <- 0.9022
AposThresh <- 2.600
annual_rate <- 0.601
var_annual_rate <- 0.934
# mean_within_above <- 0.08 # PIB: 8% ± 7% from Tolboom 2009 paper
# # 10.3% ± 1.25% for pTau217 alzpath from brum 2023
# # 16.7% ± 1.8% for pTau181 gothenberg simoa from brum 2023
sd_within_above <- 0.03 * mean_within_above
mean_within_below <- mean_within_above / 2  # PIB: 4.4% ± 4.2% from Tolboom paper
sd_within_below <- 0.03 * mean_within_below

avg_interval <- 2.02
interval_noise <- 1.06

xi <- 0
omega <- 0.5
alpha <- 1

num_frames <- 100
#Generating synthetic data
tmp <- generate_synthetic_data_with_noise(n_ids = n_ids, 
                                          n_visits = 8,
                                          avg_interval = avg_interval,
                                          interval_noise = interval_noise,
                                          minMeaningfulZ = minMeaningfulZ,
                                          AposThresh = AposThresh,
                                          maxMeaningfulZ = 25,
                                          annual_rate = annual_rate,
                                          var_annual_rate = var_annual_rate,
                                          mean_within_above = mean_within_above,  
                                          sd_within_above   = sd_within_above,
                                          mean_within_below = mean_within_below, 
                                          sd_within_below   = sd_within_below,
                                          plateau_decay = 0.15,
                                          prob_neg = 0.5,
                                          xi = xi,
                                          omega = omega,
                                          alpha = alpha,
                                          seed = seed)

df <- tmp[[1]]
converters <- tmp[[2]] #identify converters

df <- merge(df, converters, by = "ID", all.x = TRUE)
colnames(df)[3] <- "Z"
df <- df[, c("ID", "TimefromBaseline", "Z", "Z_true", "ever_cross", "first_cross")]
df[is.na(df$ever_cross),]$ever_cross <- FALSE


#Applying a manual shift to make it so not everyone converts within the first 5 years
df <- adjust_crossing_times(df,
                            mean_within_below = mean_within_below,
                            sd_within_below = sd_within_below)

##randomly sampling consecutive visits from the originally simulated data
#I will use the excluded data to get observed conversion dates outside the bounds
#of the data that I was training

df <- df[df$TimefromBaseline < 16,]

df$Observed_Conversion_Time <- approx(df$Z, df$TimefromBaseline, 2.6)$y
df$Observed_Conversion_Time <- ifelse(df$ever_cross == TRUE, df$Observed_Conversion_Time, NA)
df$Observed_Time_from_Positivity <- df$TimefromBaseline - df$Observed_Conversion_Time


 df_sample <- df %>%
    group_by(ID) %>%
    arrange(TimefromBaseline, .by_group = TRUE) %>%  # ensures chronological order
    group_modify(~ {
      n <- nrow(.x)
      if (n < num_of_scans_for_model) {
        return(.x)  # not enough rows, return all
      } else {
        start_idx <- sample(1:(n - (num_of_scans_for_model - 1)), 1)
        return(.x[start_idx:(start_idx + (num_of_scans_for_model - 1)), ])
      }
    }) %>%
    ungroup()


 #Linear regression per id per iteration
 roc_df_all <- df_sample %>%
   filter(!is.na(Z), !is.na(TimefromBaseline)) %>%
   group_by(ID) %>%
   filter(n() >= 2) %>%   # must have at least 2 points per ID
   summarise(
     model = list(lm(Z ~ TimefromBaseline)),
     .groups = "drop"
   )
 
 roc_df_all <- roc_df_all %>%
   mutate(
     slope = map_dbl(model, ~ coef(.x)[["TimefromBaseline"]]),
     intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]])
   )
 
roc_df_all <- merge(roc_df_all[, c("ID", "slope", "intercept")], df_sample[!duplicated(df_sample$ID), c("ID", "Z")],
                    by = "ID", all = FALSE)
 
 
 df_res <- get_Time_to_Positivity(df , PET_pos_threshold = AposThresh,  id_name = "ID", 
                                  time_name = "TimefromBaseline", value_name = "Z", degree = 3)


df_sample$Predicted_Time_to_Positivity <- approx(  df_res$actual_predicted_val, df_res$Time_to_Positivity,as.numeric(df_sample$Z))$y
ID1 <- "22"
ID2 <- "16"
ID3 <- "32"
ID4 <- "12"
df_sample$Observed_Time_to_Positivity <- df_sample$TimefromBaseline - df_sample$first_cross

ggplot(df_sample, aes(x = Predicted_Time_to_Positivity, y = Observed_Time_to_Positivity)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() + xlim(c(-10, 18)) + ylim(c(-10, 18))

ggplot(df_sample, aes(x = TimefromBaseline, y = Z, group = ID)) +
  geom_point() +
  geom_line(alpha = 0.4) +
  theme_classic() +
  labs(y = "Biomarker (Z)", x = "TimefromBaseline") +
  geom_point(data = df_sample[df_sample$ID == ID1,], aes(x = TimefromBaseline, y = Z, group = 1), colour = "#fc8d59", size = 3) +
  geom_point(data = df_sample[df_sample$ID == ID2,], aes(x = TimefromBaseline, y = Z, group = 1), colour = "#007B82", size = 3) +
  geom_point(data = df_sample[df_sample$ID == ID3,], aes(x = TimefromBaseline, y = Z,  group = 1), colour = "#C77CFF", size = 3) 
  
  
 
ggplot(df_sample, aes(x = Predicted_Time_to_Positivity, y = Z, group = ID)) +
  geom_point() +
  geom_line(alpha = 0.4) +
  theme_classic() +
  labs(y = "Biomarker (Z)", x = "Time from Biomarker Positivity") +
  geom_point(data = df_sample[df_sample$ID == ID1,], aes(x = Predicted_Time_to_Positivity, y = Z, group = 1), colour = "#fc8d59", size = 3) +
  geom_point(data = df_sample[df_sample$ID == ID2,], aes(x = Predicted_Time_to_Positivity, y = Z, group = 1), colour = "#007B82", size = 3) +
  geom_point(data = df_sample[df_sample$ID == ID3,], aes(x = Predicted_Time_to_Positivity, y = Z, group = 1), colour = "#C77CFF", size = 3) 



get_bs <- function(df, PET_pos_threshold, id_name, time_name, value_name,
                   num_bootstraps = 1000, bootstrap_percent = 0.8, degree = 3, printIter = TRUE) {
  df_res <- list()
  df_input_data <- list()
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  # List to store bootstrapped results
  df_bs <- vector("list", num_bootstraps)
  
  # Perform bootstrap resampling
  for (i in 1:num_bootstraps) {
    if(printIter == TRUE){
      print(i)
    }
    # Sample IDs with replacement
    sampled_ids <- sample(seq(from = 1, to = length(unique(df[[ id_name]]))), 
                          size = floor(bootstrap_percent * length(unique(df[[id_name]]))), replace = TRUE)
    ID_vec <- unique(df[[id_name]])
    ID_vec <- ID_vec[sampled_ids]
    # Subset the data for the sampled IDs
    df <- data.frame(df)
    sampled_data <- df[df[[id_name]] %in% (ID_vec),] 
    
    # Apply get_Time_to_Positivity on the sampled data
    result <- get_Time_to_Positivity(sampled_data,id_name, time_name, value_name, PET_pos_threshold, degree)
    
    # Store the actual_predicted_val from this iteration
    df_bs[[i]] <- data.frame("val_bs" = result$actual_predicted_val, 
                             "time_bs" = result$Time_to_Positivity)
    interpolated_val <- vector("list", length(Time_Window))
    
    for(j in 1:length(Time_Window)){
      
      interpolated_val[[j]] <- approx( result$Time_to_Positivity, result$actual_predicted_val, as.numeric(Time_Window[j]))$y
    }
    df_res[[i]] <- data.frame("Time_Window" = Time_Window)
    
    df_res[[i]]$interpolated_val <- unlist(interpolated_val)
    df_res[[i]]$iteration <- i
    df_input_data[[i]] <- sampled_data
    df_input_data[[i]]$iteration <- i
    
  }
  
  return(list(df_res, df_input_data))}
res <- get_bs(df_sample, PET_pos_threshold = AposThresh, id_name = "ID", time_name = "TimefromBaseline", value_name = "Z",
              num_bootstraps = num_frames, bootstrap_percent = 0.8, degree = 3, printIter = TRUE)

input_data <- res[[2]]
res <- res[[1]]



# bind all data + track which list element it came from
df_all <- bind_rows(res, .id = "frame") %>%
  mutate(step = as.integer(frame))

df_all_inputs <- bind_rows(input_data, .id = "frame") %>%
  mutate(step = as.integer(frame))

#Linear regression per id per iteration
roc_df <- df_all_inputs %>%
  filter(!is.na(Z), !is.na(TimefromBaseline)) %>%
  group_by(iteration, ID) %>%
  filter(n() >= 2) %>%   # must have at least 2 points per ID
  summarise(
    model = list(lm(Z ~ TimefromBaseline)),
    .groups = "drop"
  )

roc_df <- roc_df %>%
  mutate(
    slope = map_dbl(model, ~ coef(.x)[["TimefromBaseline"]]),
    intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]])
  )







ids <- unique(df_sample$ID)

interp_df <- map_dfr(1:num_frames, function(iter) {

  df_iter <- df_all %>% filter(iteration == iter)

  map_dfr(ids, function(id) {

    target_x <- df_sample %>%
      filter(ID == id) %>%
      slice(1) %>%
      pull(Z) %>%
      as.numeric()

    approx(
      x = df_iter$interpolated_val,
      y = df_iter$Time_Window,
      xout = target_x,
      rule = 2
    ) %>%
      as.data.frame() %>%
      mutate(
        iteration = iter,
        ID = id
      )
  })
})

interp_df <- interp_df %>%
  rename(Z = x, Time_Window_interp = y)

df_sample_baseline <- df_sample[!duplicated(df_sample$ID),c("ID", "Observed_Time_to_Positivity")]
interp_df <- merge(interp_df, df_sample_baseline, by = "ID")


id_summary <- interp_df %>%
  group_by(ID) %>%
  summarise(
    obs_ttp = first(Observed_Time_to_Positivity)
  )

id_summary <- id_summary %>%
  filter(
    !is.na(obs_ttp),
    is.finite(obs_ttp)
  )

id_summary <- id_summary %>%
  mutate(rank = ntile(obs_ttp, 4))

selected_ids <- id_summary %>%
  group_by(rank) %>%
  slice_sample(n = 1) %>%
  ungroup()


interp_cum <- map_dfr(1:num_frames, function(i) {
  interp_df %>%
    filter(iteration <= i) %>%
    mutate(frame = i)
})


ids <- as.factor(selected_ids$ID)

tmp <- interp_cum[interp_cum$ID %in% ids, c("ID", "Observed_Time_to_Positivity")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[order(tmp$Observed_Time_to_Positivity),]
ids <- factor(tmp$ID,
              levels = tmp$ID)



make_originalData_plot <- function(f) {
  
  df_f <- df_all_inputs[df_all_inputs$frame == f,]

  ggplot(df_sample, aes(x = TimefromBaseline, y = Z, group = ID), alpha = 0.3) +
    geom_point(data = df_sample, aes(x = TimefromBaseline, y = Z, group = ID), alpha = 0.3) +
    geom_line(data = df_sample, aes(x = TimefromBaseline, y = Z, group = ID), alpha = 0.2) +
    geom_point(data = df_f, aes(x = TimefromBaseline, y = Z, group = ID)) + 
      geom_line(data = df_f, aes(x = TimefromBaseline, y = Z, group = ID), alpha = 0.8) + theme_classic() + xlab("Time from Baseline") +
    geom_point(data = df_sample[df_sample$ID == ids[1],], 
               aes(x = TimefromBaseline, y = Z), colour = "#fc8d59", size = 3) +
    geom_line(data = df_sample[df_sample$ID == ids[1], ], 
              aes(x = TimefromBaseline, y = Z), colour = "#fc8d59", linewidth = 1) +
    geom_point(data = df_sample[df_sample$ID == ids[2],],  
               aes(x = TimefromBaseline, y = Z), colour = "#C77CFF", size = 3) +
    geom_line(data = df_sample[df_sample$ID == ids[2],],  
              aes(x = TimefromBaseline, y = Z), colour = "#C77CFF", linewidth = 1) +
    geom_point(data = df_sample[df_sample$ID == ids[3],],  
               aes(x = TimefromBaseline, y = Z), colour = "#007B82", size = 3) +
    geom_line(data = df_sample[df_sample$ID == ids[3],],  
              aes(x = TimefromBaseline, y = Z), colour = "#007B82", linewidth = 1) +
    geom_point(data = df_sample[df_sample$ID == ids[4],],  
               aes(x = TimefromBaseline, y = Z), colour = "#F2C14E", size = 3) +
    geom_line(data = df_sample[df_sample$ID == ids[4],],  
              aes(x = TimefromBaseline, y = Z), colour = "#F2C14E", linewidth = 1) +
    ylab("Biomarker Level (Z)")
}


make_ROC_plot <- function(f) {
  
  roc_df_tmp <- merge(df_all_inputs[df_all_inputs$iteration == f, c("ID", "iteration", "Z", "TimefromBaseline")],
                      roc_df[roc_df$iteration == f, c("ID", "iteration", "slope")], by = "ID", all = FALSE)
  roc_df_tmp <- roc_df_tmp[!duplicated(roc_df_tmp$ID),]
  
  
  ggplot(roc_df_tmp, aes(x = Z, y = slope)) + geom_point() + theme_classic() +
    geom_point(data = roc_df_all, aes(x = Z, y = slope), alpha = 0.3) + 
    xlab("Baseline Biomarker Level (Z)") + ylab("Annualized Rate of Change (Z/yr)") +
    geom_point(data = roc_df_all[roc_df_all$ID == ids[1],], 
               aes(x = Z, y = slope), colour = "#fc8d59", size = 3) +
    geom_point(data = roc_df_all[roc_df_all$ID == ids[2],], 
               aes(x = Z, y = slope), colour = "#C77CFF", size = 3) +
    geom_point(data = roc_df_all[roc_df_all$ID == ids[3],], 
               aes(x = Z, y = slope), colour = "#007B82", size = 3) +
    geom_point(data = roc_df_all[roc_df_all$ID == ids[4],], 
               aes(x = Z, y = slope), colour = "#F2C14E", size = 3) +
    ylim(c(mean(roc_df_all$slope, na.rm = TRUE) - 3 *sd(roc_df_all$slope, na.rm = TRUE),  
           mean(roc_df_all$slope, na.rm = TRUE) + 3 *sd(roc_df_all$slope, na.rm = TRUE)))

}







make_top_plot <- function(f) {
  
  df_f <- df_all %>% filter(frame <= f)
  
  df_line <- interp_cum %>% filter(frame <= f)
  
  vline_df <- df_line %>%
    group_by(ID) %>%
    summarise(
      Z = first(Z, na.rm = TRUE),
      .groups = "drop"
    ) 
  vline_df$ID <- as.factor(vline_df$ID)
  vline_df <- vline_df[vline_df$ID %in% tmp$ID,]
  vline_df$ID <- factor(tmp$ID,
                    levels = tmp$ID)
  ggplot(df_f, aes(x = Time_Window, y = interpolated_val)) +
    geom_line(aes(group = frame), colour = "black", alpha = 0.2) +
    # geom_hline(yintercept = vline_df[vline_df$ID == ids[1],]$Z, linetype = "dashed", colour = "#C77CFF") +
    # geom_hline(yintercept = vline_df[vline_df$ID == ids[2],]$Z, linetype = "dashed", colour = "#fc8d59") +
    # geom_hline(yintercept = vline_df[vline_df$ID == ids[3],]$Z, linetype = "dashed", colour = "#007B82") +
    # geom_hline(yintercept = vline_df[vline_df$ID == ids[4],]$Z, linetype = "dashed", colour = "#F2C14E") +
    theme_classic() +
    labs(
      x = "Time from Biomarker Positivity",
      y = "Biomarker (Z)"
    ) + xlim(c(-10, 20)) + ylim(c(-9, 50))
}

make_bottom_plot <- function(f) {
  interp_cum$ID <- as.character(interp_cum$ID)
  
  df_f <- interp_cum %>% filter(frame <= f)
 
  vline_df <- df_f %>%
    group_by(ID) %>%
    summarise(
      Observed_Time_to_Positivity = median(Observed_Time_to_Positivity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(Observed_Time_to_Positivity))
  
  vline_df <- vline_df[vline_df$ID %in% ids,]
  
  label_df <- interp_cum %>%
    distinct(ID, Observed_Time_to_Positivity) %>%
    mutate(
      facet_label = as.numeric(
        round(Observed_Time_to_Positivity, 2)
      )
    )
  
  label_df <- label_df[label_df$ID %in% ids,]
  
  label_lookup <- setNames(
    paste0(
      #label_df$ID,
      #"\n
      "Obs TFP = ",
      round(label_df$Observed_Time_to_Positivity, 2)
    ),
    label_df$ID
  )
  
  df_f <- interp_cum %>%
    filter(frame <= f) %>%
    left_join(label_df, by = "ID")
  
  df_f <- df_f[df_f$ID %in% ids,]
  df_f$ID <- factor(df_f$ID, levels = rev(ids))
  label_df$ID <- factor(label_df$ID, levels = rev(ids))
  
  ggplot(df_f, aes(x = Time_Window_interp, fill = ID)) +
    geom_histogram(  binwidth = binwidth,
                     boundary = x_min, alpha = 0.7) + 
    scale_fill_manual(values = c("#F2C14E", "#C77CFF", "#007B82", "#fc8d59")) + #, "#4C6FE7"
    geom_vline(data = label_df, aes(xintercept = Observed_Time_to_Positivity),
      linetype = "dashed",
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ ID, labeller = labeller(ID = label_lookup), nrow = 1) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.line.y = element_blank()) +
    labs(
      x = "Time from Biomarker Positivity (TFP)",
      y = ""
    ) + scale_x_continuous(limits = c(-10, 20)) + scale_y_continuous(limits = c(0, max_count))
}


n_frames <- max(as.numeric(df_all$frame))
x_min <- min(interp_cum$Time_Window_interp, na.rm = TRUE)
x_max <- max(interp_cum$Time_Window_interp, na.rm = TRUE)
binwidth <- (x_max - x_min) / 30

max_count <- interp_cum %>%
  count(ID, 
        bin = cut(Time_Window_interp, breaks = seq(x_min, x_max, by = binwidth))) %>%
  summarise(max_n = max(n, na.rm = TRUE)) %>%
  pull(max_n)
n_frames <- num_frames


for (f in 1:n_frames) {
  
  originalData_plot <- make_originalData_plot(f)
  ROC_plot <- make_ROC_plot(f)
  top_plot <- make_top_plot(f)
  bottom_plot <- make_bottom_plot(f)
  
  layout_matrix <- rbind(c(1, 2, 2),
                         c(3, 2, 2),
                         c(4, 4, 4))
  
  # ---- draw grobs into a captured raster ----
  img <- grid.grabExpr(
    grid.arrange(
      originalData_plot,
      top_plot,
      ROC_plot,
      bottom_plot,
      layout_matrix = layout_matrix
    )
  )
  
  # convert to magick image
  combo <- image_graph(width = 1500, height = 900, res = 150)
  grid.draw(img)
  dev.off()
  
  # IMPORTANT: image_graph already produced a magick object
  image_write(combo,
              sprintf("./animations_simulation/4panelimage_03percwithin_xi-omegapoint5alpha1/frame_%04d.png", f))
}
##Wait to do this till after manually manipulating the frames that have been written to frames_combo

frames <- list.files("./animations_simulation/4panelimage_03percwithin_xi-omegapoint5alpha1", full.names = TRUE)
frames <- sort(frames)
animation <- image_read(frames)
animation <- image_animate(animation, fps = 10)

image_write(animation, "./animations_simulation/4panelimage_03percwithin_xi-omegapoint5alpha1.gif")
