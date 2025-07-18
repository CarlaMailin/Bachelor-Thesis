## calculating the correlation between melt layer frequency as well as melt layer thickness and difference between theoretical and estimated diffusion lengths using data from Samira Zander

##----------------------------------------------
# import library
library(dplyr)
library(ggplot2)
#-----------------------------------

# define dataframe with melt data for each core
#melt_layer <- data.frame(site = c( "exB22", "exB16", "cleaned_B22", "exB18","exB21", "exB23", "exB26", "not_cleaned_B22" ),

melt_layer <- data.frame(
  site = c(
    "B16-2019",
    "cleaned_B22-2019",
    "B18-2012",
    "B21-2012",
    "B23-2012",
    "B26-2012",
    "n_c_B22-2019"
  ),
  
  m_freq = c(
    11.666667,
    3.441860,
    16.444444,
    5.727273,
    16.375000,
    14.000000,
    3.441860
  ),
  m_mean_thickness = c(
    0.570667,
    0.610500,
    0.411111,
    0.429455,
    0.587500,
    0.241500,
    0.610500
  ) # in cm
)


# calculate difference between simulated and estimated diffusion length
comparison_sig$difference <- comparison_sig$sim_diffusion_length - comparison_sig$sigma_est
#comparison_sig$q5 <- comparison_sig$q5 + comparison_sig$difference
#comparison_sig$q95 <- comparison_sig$q95 + comparison_sig$difference

print(comparison_sig$difference) #debug
print(comparison_sig) #debug


# add difference to data frame
difference_sigma <- comparison_sig[, c("site", "sim_diffusion_length", "sigma_est", "difference")]
print(difference_sigma)

# join all data - estimated and calculated as well as melt data in a single data frame,each with a corresponding site
comparison_data <- left_join(comparison_sig, melt_layer, by = "site")
print(comparison_data)

# plot melt layer frequency vs difference between estimated and calculated diffusion length
freq_plot <- ggplot(comparison_data, aes(y = m_freq, x = difference, color = site)) +
  geom_point(size = 3, shape = 16) +  # shape = 4 is "+"
  geom_errorbar(aes(xmin = sim_diffusion_length - q5, xmax = sim_diffusion_length - q95),
                width = 0.2) +
  labs(y = "Melt Frequency", x = "Deviation of diffusion length in cm", color = "Site") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
  )+
  theme(plot.background = element_rect(fill = "white", color = NA))
theme(legend.position = "right")


thickness_plot <-  ggplot(comparison_data,
                          aes(x = difference, y = m_mean_thickness, color = site)) +
  geom_point(size = 3, shape = 16) +  # shape = 4 is "+"
  geom_errorbar(aes(xmin = sim_diffusion_length - q5, xmax = sim_diffusion_length - q95),
                width = 0.01) +
  labs(y = "Mean Thickness of Melt Layer in cm", x = "Deviation of diffusion length in cm", color = "Site") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
  )+
  theme(plot.background = element_rect(fill = "white", color = NA))
theme(legend.position = "right")


print(freq_plot)
print(thickness_plot)

ggsave(
  filename = paste0("plots/", depth_range_str, "melt_frequency_plot.pdf"),
  plot = freq_plot,
  bg = "white",
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  filename = paste0("plots/", depth_range_str, "thickness_plot.pdf"),
  plot = thickness_plot,
  bg = "white",
  width = 6,
  height = 4,
  dpi = 300
)


ggsave(
  filename = paste0("plots/", depth_range_str, "melt_frequency_plot.png"),
  plot = freq_plot,
  width = 6,
  height = 4
)

ggsave(
  filename = paste0("plots/", depth_range_str, "thickness_plot.png"),
  plot = thickness_plot,
  width = 6,
  height = 4
)

elevation_df_site<- elevation_df_renamed %>% rename(site = site_i)
comparison_data_el <- comparison_data %>%
  left_join(elevation_df_renamed, by = c("site" = "site_i"))
print(comparison_data_el)

comparison_data %>% count(site) %>% filter(n > 1)
elevation_df_renamed %>% count(site_i) %>% filter(n > 1)



# Now plot with elevation color coding
comp_plot2 <- ggplot(comparison_data_el, aes(x = sim_diffusion_length, y = sigma_est, label = site, color = elevation)) +
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
  theme_minimal()

print(comp_plot)


ggsave(
  filename = paste0("plots/", depth_range_str, "elv_final_comparison_est_calc_sigma.pdf"),
  plot = comp_plot2,
  width = 6,
  height = 4)
