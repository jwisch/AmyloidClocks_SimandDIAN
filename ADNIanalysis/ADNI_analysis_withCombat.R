library(data.table)
library(lubridate)
library(sva)
library(tableone)
source("./functions.R")
source("./Simulation_functions.R")
source("./Fig1and2Funcs.R")
source("./functions_DistributionSpecificity.R")

source("./ADNI_prep_for_analysis.R")

p2b <- ggplot(RofC_fits_ADNI, aes(x = omega, y = xi, label = biomarker_col, colour = Assay,
                      shape = Analyte)) +
  geom_point(aes(size = alpha, alpha = 0.6)) + scale_colour_manual(values = c("#007B82", "#C77CFF",  "green", "black","#F8766D" )) +
  scale_shape_manual(values = c(19, 17,19)) +
  # geom_text(vjust = 1.3, hjust = -0.1) +
  facet_wrap(~modality) +
  theme_bw() +
  labs(y = "xi", x = "omega") +
  theme(legend.position = "bottom") + 
  xlim(c(0, 7.5)) + ylim(c(-4, 0.6))




p_PET  <- get_TimefromAposPlot(result_PET_Z, cp_PET_Z$optimal_cutpoint, 
                               "Amyloid PET -  Z")
p_pT217  <- get_TimefromAposPlot(result_pT217_Z, cp_pT217_Z$optimal_cutpoint, 
                               "Plasma pT217/AB42 Z")

grid.arrange(p_PET, p_pT217, nrow = 2, ncol = 1)






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


p_4_PET <- plot_hist_density_at_time(data = result_PET_dist_Z, 
                                time_vals = c(0, 2, 4, 6),
                                binwidth = 0.01,
                                color_vec =   c("#003B44",   # darker lead-in
                                                "#007B82",   # 2nd (your first required color)
                                                "#00BFC4",   # 3rd (your second required color)
                                                "#6FE3E6",   # lighter continuation
                                                "#CFF7F8") ,  # very light endpoint
                                "Amyloid PET")
p_4_PET <- p_4_PET + xlim(c(-0.5, 2))
colnames(result_pT217_dist_Z) <- c("Time_to_Positivity", "Estimate")
p_4_plasma <- plot_hist_density_at_time(data = result_pT217_dist_Z, 
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
p_4_plasma <- p_4_plasma + xlim(c(-0.5, 2))

layout_matrix <- rbind(c(1, 1), c(2, 2), c(3, 4))
grid.arrange(p2b, p3b, p_4_PET, p_4_plasma, layout_matrix = layout_matrix)
grid.arrange(p_PET, p_error_PET[[1]],
             p_pT217, p_error_pT217[[1]], nrow = 2, ncol = 2)

grid.arrange(p_4_PET, p_4_plasma, nrow = 1, ncol = 2)
graph2ppt(file = "./Figures/Fig4c.pptx", width = 6.5, height = 5)





p1 <- ggplot(matched[!duplicated(matched$PTID) & matched$PTID %in% counts$PTID,], aes(x = CENTILOIDS.combat, y = pT217_AB42_F)) +
  geom_point() + theme_bw() + geom_vline(xintercept = cp_PET_Z$optimal_cutpoint * sd_CL + mean_CL, linetype = "dashed") +
  geom_hline(yintercept = cp_pT217_Z$optimal_cutpoint * sd_pT + mean_pT, linetype = "dashed") + 
  ggtitle("A.") + xlab("Cortical Amyloid (CL)") + ylab("Plasma pTau217/AB42 - Fujirebio")

p2 <- cp_PET_pT217_obj[[1]] + ggtitle("B.")
grid.arrange(p1, p2, nrow = 1, ncol = 2)
graph2ppt(file = "./Figures/ADNIBasicDataandCutpoints.pptx", width = 6.5, height = 4)

matched$pTpos <- ifelse(matched$pT217_AB42_F > (cp_pT217_Z$optimal_cutpoint * sd_pT + mean_pT), 1, 0)
matched$CLpos <- ifelse(matched$CENTILOIDS.combat > (cp_PET_Z$optimal_cutpoint * sd_CL + mean_CL),1 , 0)
table(matched[!duplicated(matched$PTID) & matched$PTID %in% counts$PTID,]$CLpos, matched[!duplicated(matched$PTID) & matched$PTID %in% counts$PTID,]$pTpos)


plot_df <- matched[matched$PTID %in% counts$PTID,]

