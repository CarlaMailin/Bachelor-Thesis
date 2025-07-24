##################################################################
#' Calculations based on DensityHL function from the FirnR package, where the Johnsen Correction factor was changed to TRUE
#'Date compiled using the ClimPar function from the exNGT R package (https://github.com/EarthSystemDiagnostics/exNGT) using data from (doi.org/10.5194/cp-12-171-2016) Weichbach et al.
#' calculates the density profiles for firn in NGT cores

#-------------------------------------
# import packages
library(dplyr)
library(viridis)  # for scale_color_viridis_d()
library(ggplot2)
library(patchwork)
source("your_filepath/FirnR/R/DensityHL.R")
#-----------------------------------------------
# read in data with temperature and accumulation rate for all sites
climdata <- read.csv(
  "your_filepath/climatePar.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)
#-----------------------------------------------------------------
# define constants
depth <- 0:100  # in m - can be adjusted - value is around deepest depth value of NGT core
rho.s <- 320 #329.6133  # surface density value in kg/m³ rounded from 329.613, which is the mean of surface density in Schaller (2016).
ngt_sites <- c("B16",
               "B18",
               "B21",
               "B22",
               "B23",
               "B26", 
               "GRIP")
elevation <- climdata$Elevation   
print(elevation)
# create lists for objects created in loop
density_results <- list() # density profile result of individual core
mean.den <- list() # mean firn density for depth range in [kg/m^3]
max.den <- list() # maximum firn density in depth interval in [kg/m^3]
min.den <- list() # minimum firn density in depth interval in [kg/m^3]
f_m_mean <- list() # mean frequency of annual cycle [1/m] with density profile in the depth range used in calculation
f_m_max <- list() # max. frequency of annual cycle [1/m] with density profile in the depth range used in calculation
f_m_min <- list() # min. frequency of annual cycle [1/m] with density profile in the depth range used in calculation
h_mean <- list() # mean core depth [m] of annual cylce with density profile in depth range used in calculation
h_max <- list() # max. core depth [m] of annual cylce with density profile in depth range used in calculation
h_min <- list() # min. core depth [m] of annual cylce with density profile in depth range used in calculation
df_density_all  <- data.frame()

# calculate density for each site with loop
for (i in seq_along(ngt_sites)) {
  site_j <- ngt_sites[i] # select site from site list
  site_data <- climdata %>% filter(Site == site_j) # filter through sites for selected site
  
  t.mean  <- site_data$meanTemperature + 273.15 # mean temperature in Kelvin
  bdot    <- site_data$accRate # mean accumulation Rate in kg/m^2/yr
  elevation_site <- site_data$Elevation
  #densification model from FirnR package DensityHL
  density_results <- DensityHL(
    depth = depth,
    rho.surface = rho.s,
    T = t.mean,
    bdot = bdot,
   JohnsenCorr = TRUE
  )
  
  # Save density results to dataframe with site and depth data
  df_density <- data.frame(site_j = site_j,
                           rho = density_results$rho,
                           # calc.
                           depth = density_results$depth.we,
                           elevation = elevation_site
                           )
  
  depth_range <- df_density$depth >= 23.2 &
    df_density$depth <= 30 # define a depth range to calculate for specific depth intevals
  filtered_df_density <- df_density[depth_range, ] # filtered density profile for specific depth interval
  mean.den[[site_j]] <- mean(filtered_df_density$rho, na.rm = TRUE) # mean density in depth interval
  max.den[[site_j]] <- max(filtered_df_density$rho, na.rm = TRUE) # max. density in depth interval
  min.den[[site_j]] <- min(filtered_df_density$rho, na.rm = TRUE) # min. density in depth interval
  #------------------------------------
  # define function for annual cycle
  h_mean[[site_j]] <- bdot / mean.den[[site_j]] # depth/height of firn mean accumulation per year
  h_min[[site_j]] <- bdot / min.den[[site_j]] # depth/height min firn accumulation in used depth interval
  h_max[[site_j]] <- bdot / max.den[[site_j]] # depth/height max firn accumulation in used depth interval
  f_m_mean[[site_j]] <- 1 / h_mean[[site_j]] # mean annual cycle freq in depth domain
  f_m_max[[site_j]] <- 1 / h_max[[site_j]] # max. annual cycle freq in depth domain
  f_m_min[[site_j]] <- 1 / h_min[[site_j]] # min. annual cycle freq in depth domain
  
  # save each density profiles of core after calculation in data frame
  df_density_all <- rbind(df_density_all, df_density)
}

#define function to make individul lists with calculated parameters into dataframe
list_to_df <- function(lst, colname) {
  data.frame(site = names(lst),
             value = unlist(lst),
             row.names = NULL) %>%
    rename(!!colname := value)
}
print(df_density_all)
# define data frame with all results included
summary_density_df <- bind_cols(
  list_to_df(mean.den, "mean_density"),
  list_to_df(min.den, "min_density")["min_density"],
  list_to_df(max.den, "max_density")["max_density"],
  list_to_df(f_m_mean, "f_m_mean")["f_m_mean"],
  list_to_df(f_m_max, "f_m_max")["f_m_max"],
  list_to_df(f_m_min, "f_m_min")["f_m_min"]
)

print(summary_density_df)


# plot firn density profiles over entire calculated depth
densityplot <- ggplot(df_density_all, aes(x = rho, y = depth, color = site_j)) +
  geom_line(linewidth = 1) +
  scale_y_reverse() +  # so depth increases downward
  labs(x = "Density [kg/m³]", y = "Depth [m]", #title = "Density Profiles by Site",
       color = "Site") +
  theme_minimal()

print(densityplot)

ggsave(
  filename = paste0("plots/", "NGT_Density_Profile.pdf"),
  plot = densityplot,
  width = 6,
  height = 4
    )


ggsave(
  filename = paste0("plots/", "NGT_Density_Profile.png"),
  plot = densityplot,
  width = 6,
  height = 4
)



# plot using an elevation color scale 
site_name_map <- c("B16"= "B16-2019", "B18"= "B18-2012", "B21"= "B21-2012", "B22"= "B22-2019", "B23" = "B23-2012", "B26" = "B26-2012", "GRIP" =  "GRIP")

# Rename sites in df_density_all using site_name_map
df_density_all <- df_density_all %>%
  mutate(site_j = if_else(site_j %in% names(site_name_map), site_name_map[site_j], site_j))

# Summarize to get one elevation and mean depth per site (to position labels properly)
site_labels <- df_density_all %>%
  group_by(site_j, elevation) %>%
  summarise(depth = mean(depth), .groups = "drop") %>%
  arrange(desc(elevation)) %>%
  mutate(label_y = seq(from = max(depth), to = min(depth), length.out = n()))  # Add label_y

# Main plot
p1 <- ggplot(df_density_all, aes(x = rho, y = depth, group = site_j)) +
  geom_line(aes(color = elevation), linewidth = 1) +
  scale_color_viridis_c(option = "viridis", name = "Elevation [m]") +
  scale_y_reverse(limits = c(max(df_density_all$depth), min(df_density_all$depth))) +
  labs(x = "Density [kg/m³]", y = "Depth [m]") +
  theme_minimal() +
  theme(legend.position = "right")

# Label plot next to main plot
p2 <- ggplot(site_labels, aes(y = label_y)) +
  geom_point(aes(x = 1, color = elevation), size = 4) +       # colored dots
  geom_text(aes(x = 1.1, label = site_j), color = "black", hjust = 0, size = 4) +  # black text next to dots
  scale_color_viridis_c(option = "viridis", guide = "none") +  # no legend for points here
  scale_y_continuous(breaks = site_labels$label_y, labels = NULL, expand = c(0.1, 0.1)) +
  scale_x_continuous(limits = c(0.9, 1.5)) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine side by side
elevation_density_plot <- p1 + p2 + plot_layout(widths = c(4, 2))
print(elevation_density_plot)
ggsave(
  filename = paste0("plots/", "elevation_density.pdf"),
  plot = elevation_density_plot,
  width = 6,
  height = 4)
