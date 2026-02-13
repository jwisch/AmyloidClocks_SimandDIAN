#This figure generates the top and bottom of Figure 3 separately
#These figures were combined in powerpoint post hoc
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(export)
library(gridExtra)
library(sn)

# Folder containing your RDS files
folder <- "./RDS_Simulation_v02_fullBS/"

# List all RDS files in the folder
all_files <- list.files(folder, pattern = "\\.RDS$", full.names = TRUE)


###Generates the cartoon distributions for the top of the plot
simulate_dist_single <- function(
    xi,
    omega,
    alpha,
    n = 5000
) {
  
  # --- Simulate data ---
  sim_data <- data.frame(
    values = rsn(n, xi = xi, omega = omega, alpha = alpha)
  )
  
  # --- Plot ---
  p <- ggplot(sim_data, aes(x = values)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 50,
      fill = "#00BFC4",
      color = "#00BFC4",
      alpha = 0.7
    ) +
    geom_density(color = "#007B82", linewidth = 1) +
    theme_bw(base_size = 13) +
    theme(
      # remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # remove outer box / panel border
      panel.border = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    scale_x_continuous(limits = c(-5, 5)) +
    ylim(c(0, 1.2)) +
    geom_vline(xintercept = 0, linetype = "dashed")
  
  return(p)
}

p1 <- simulate_dist_single(xi = 0, omega = 0.5, alpha = 0)
p2 <- simulate_dist_single(xi = 1, omega = 0.5, alpha = 0)
p3 <- simulate_dist_single(xi = 0, omega = 0.5, alpha = 2)
p4 <- simulate_dist_single(xi = 1, omega = 0.5, alpha = 2)
grid.arrange(p1, p2, p3, p4, nrow = 1, ncol = 4)

graph2ppt(file = "./figures/MAEDistribution_withinIndividVariation_FixedXiAlphaOmega_Distributions.pptx", width = 9, height = 2)



#Reading in generated data and making corresponding plots

# Define a regex pattern to extract the parts
# Pattern breakdown:
# ^(.*?)_X(.*?)Omega(.*?)Alpha(.*?)_mean_within_above_(.*?)\\.RDS$
# 1: METRIC
# 2: VALUE1
# 3: VALUE2
# 4: VALUE3
# 5: VALUE4
pattern <- "^(.*?)_X(.*?)Omega(.*?)Alpha(.*?)_mean_within_above_(.*?)\\.RDS$"

# Filter files that match the pattern
matching_files <- all_files[str_detect(basename(all_files), pattern)]

# Function to read a file and extract info
read_file_with_metadata <- function(file) {
  fname <- basename(file)
  
  # Extract info from filename
  matches <- str_match(fname, pattern)
  
  METRIC <- matches[2]
  VALUE1 <- as.character(matches[3])
  VALUE2 <- as.character(matches[4])
  VALUE3 <- as.character(matches[5])
  VALUE4 <- as.numeric(matches[6])
  
  # Read the RDS file (assuming it contains a vector)
  values <- readRDS(file)
  
  # Return a data frame
  data.frame(
    METRIC = METRIC,
    VALUE1 = VALUE1,
    VALUE2 = VALUE2,
    VALUE3 = VALUE3,
    VALUE4 = VALUE4,
    Value = values
  )
}

# Read and combine all files
combined_df <- map_dfr(matching_files, read_file_with_metadata)

# Check the result
head(combined_df)


combined_df <- combined_df[combined_df$METRIC %in% c("MAE", "R2", "RMSE"),]
colnames(combined_df)[2:5] <- c("Xi", "Omega", "Alpha", "Mean_within_above")

combined_df <- combined_df[, c("METRIC", "Xi", "Omega", "Alpha", "Mean_within_above", "Value")]


results <- setDT(combined_df)[, .(
  median_Value = median(Value),
  mean_Value   = mean(Value),
  sd_Value     = sd(Value)
), by = .(METRIC, Xi, Omega, Alpha, Mean_within_above)]



p1 <- ggplot(results[results$METRIC == "MAE" & results$Alpha == 0,], aes(x = Mean_within_above, y = median_Value,
                                              ymin = mean_Value - 1.96 * sd_Value,
                                              ymax = mean_Value + 1.96 * sd_Value)) +
  geom_point() + geom_errorbar() + theme_bw() + facet_wrap(~Xi) + ylim(c(-6, 9)) +
  xlab("Mean within-individual Variation") + ylab("Mean Average Error") + geom_hline(yintercept = 2, linetype = "dashed")

p2<- ggplot(results[results$METRIC == "MAE" & results$Alpha == 2,], aes(x = Mean_within_above, y = median_Value,
                                                                   ymin = mean_Value - 1.96 * sd_Value,
                                                                   ymax = mean_Value + 1.96 * sd_Value)) +
  geom_point() + geom_errorbar() + theme_bw() + facet_wrap(~Xi) + ylim(c(-6, 9))+
  xlab("Mean within-individual Variation") + ylab("Mean Average Error")+ geom_hline(yintercept = 2, linetype = "dashed")

grid.arrange(p1, p2, nrow = 1)
graph2ppt(file = "./figures/MAEDistribution_withinIndividVariation_FixedXiAlphaOmega.pptx", width = 10, height = 5.78)

