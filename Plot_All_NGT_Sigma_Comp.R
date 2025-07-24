#####
#'Creates a Plot comparing the calculated and estimated diffusions lengths of d18O with errors based on the 5th and 95th percentile 
#'
#'
#--------------------------
#import packages
library(tidyr)
library(ggplot2)
library(dplyr)
library(viridis)

# load in calculated and estimated diffusion lengths
source("your_filepath/Calc_Firn_Diffusion_Length.R")
source("your_filepath/NGT_rechnungen/ALL_ExNGT_GRIP_Diffusion_Est.R")
source("your_filepath/Cleaned_Estimation_B22_Diffusionlength.R")
source("your_filepath/B22_Bag24_PSD.R")
source("your_filepath/Estimation_B22_Diffusionlength.R")
#----------------------------------------------------
# define calculated sigma and rename columns to use later
sim_sig <- mean.calc.sig_df %>%
  rename(site = site_i, sim_diffusion_length = mean.sigma)

# define estimated sigma for different Sites and convert values from m to cm 
est_sigma_summary_cm_ALL_NGT <- sigma_summary_ALL_NGT %>%
  mutate(sigma_est = sigma_est_ALL_NGT * 100, q5 = q5 * 100, q95 = q95 * 100) 

est_sigma_summary_cm_B22 <- B22_sigma_summary %>%
  mutate(sigma_est = B22_sigma_est * 100, q5 = q5 * 100, q95 = q95 * 100) 

est_sigma_summary_cm_cleaned_B22 <- cleaned_B22_sigma_summary %>%
  mutate(sigma_est = cleaned_B22_sigma_est * 100, q5 = q5 * 100, q95 = q95 * 100) 

# add all estimated sigmas into a data frame
est_sigma_summary_cm <- bind_rows(
  est_sigma_summary_cm_ALL_NGT,
  est_sigma_summary_cm_B22,
  est_sigma_summary_cm_cleaned_B22
)


# make combined dataframe with simulated and estimated diffusion lengths
comparison_sig <- sim_sig %>%
  inner_join(est_sigma_summary_cm, by = "site") 

print(comparison_sig)


Site <- comparison_sig$site # define sites


# plot the comparison with q5 and q95 as errorbars
comp_plot <- ggplot(comparison_sig, aes(x = sim_diffusion_length, y = sigma_est, label = Site, color = Site)) +
  geom_point(size = 3) +
  geom_errorbar(aes(
    ymin = q5,
    ymax = q95
  ), width = 0.2) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  #geom_text(nudge_y = 0.4, nudge_x = -0.1, check_overlap = TRUE) +
  labs(
    x = "Simulated Diffusion Length in cm",
    y = "Observed Diffusion Length in cm",
  ) +
  xlim(1, 16) +
  ylim(1, 16)+
  scale_color_viridis_d(option = "viridis") +
  theme_minimal()

print(comp_plot)

# save plot as pdf and png
ggsave(
  filename = paste0("plots/", depth_range_str, "final_comparison_est_calc_sigma.png"),
  plot = comp_plot,
  width = 6,
  height = 4)

ggsave(
  filename = paste0("plots/", depth_range_str, "final_comparison_est_calc_sigma.pdf"),
  plot = comp_plot,
  width = 6,
  height = 4)


elevation_df_site<- elevation_df_renamed %>% rename(site = site_i)
comparison_sig_el <- comparison_sig %>%
  left_join(elevation_df_renamed, by = c("site" = "site_i"))
print(comparison_sig_el)

comparison_sig %>% count(site) %>% filter(n > 1)
elevation_df_renamed %>% count(site_i) %>% filter(n > 1)



# Now plot with elevation color coding
comp_plot2 <- ggplot(comparison_sig_el, aes(x = sim_diffusion_length, y = sigma_est, label = site, color = elevation)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_text(hjust = -0.13, vjust = 0.7, size = 3, check_overlap = FALSE, color = "black") +
  labs(
    x = "Simulated Diffusion Length in cm",
    y = "Observed Diffusion Length in cm",
    color = "Elevation [m]"
  ) +
  xlim(1, 16) +
  ylim(1, 16) +
  scale_color_viridis_c(option = "cviridis") +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
  )+
  theme_minimal()

print(comp_plot)


ggsave(
  filename = paste0("plots/", depth_range_str, "elv_final_comparison_est_calc_sigma.pdf"),
  plot = comp_plot2,
  width = 6,
  height = 4)
