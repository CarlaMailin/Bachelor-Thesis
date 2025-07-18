########################################################
#' calculations are made using the FirnR package written by Thomas Münch and Thomas Laepple
#' calculate diffusion lengths from the NGT cores
#'climatological data compiled using the ClimPar function from the exNGT R package (https://github.com/EarthSystemDiagnostics/exNGT) using data from (doi.org/10.5194/cp-12-171-2016) Weichbach et al.
#'
#' * B16-B30 elevation, temperature and accumulation rate data from
#' Weissbach et al., Clim. Past, https://doi.org/10.5194/cp-12-171-2016, 2016.
#'
#' * GRIP and NGRIP elevation, temperature and accumulation rate data from
#' Vinther et al., J. Geophys. Res., https://doi.org/10.1029/2005JD006921, 2006.
#' 
#' #' * Average surface snow density (0-1 m) from
#' Schaller et al., PANGAEA, https://doi.org/10.1594/PANGAEA.867874, 2016.
#--------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
source("C:/Users/maili/Documents/Bachelorarbeit/FirnR/R/DensityHL.R")
source("C:/Users/maili/Documents/Bachelorarbeit/FirnR/R/Diffusivity.R")
source("C:/Users/maili/Documents/Bachelorarbeit/FirnR/R/DiffusionLength.R")
#-----------------------------------------------------------------------
# read in Data, compiled using exNGT package
climdata <- read.csv(
  "C:/Users/maili/Documents/Bachelorarbeit/NGT_rechnungen/climatePar.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# a list of all sites, for which the Diffusionlength will be calculated - not all of these are used later - potentially remove the unused ones?
all_sites <- c(
  "B16",
  "B17",
  "B18",
  "B20",
  "B21",
  "B22",
  "B23",
  "B26",
  "B27.28",
  "B29",
  "B30",
  "GRIP",
  "GISP2",
  "NEGIS",
  "NEEM",
  "NGRIP"
)

elevation <- climdata$Elevation  
print(elevation)
# define lists to store values in, that are calc. in loop
mean.sig.site <- list()
median.sig.site <- list()
filtered_sigma_site <- list()

# calculate Diffusionlength using loop
for (site_i in all_sites) {
  site_data <- climdata %>% filter(Site == site_i) # load in site specific data
  t.mean <- site_data$meanTemperature + 273.15 # mean temperature in Kelvin
  bdot   <- site_data$accRate  # mean accumulation Rate in kg/m^2/yr
  press.mean  <- site_data$surfacePressure # mean surface pressure in mbar
  elevation <- site_data$Elevation
  kRho.s <- 320 # surface density value in kg/m³ rounded from 329.613, which is the mean of surface density in Schaller et al.
  
  depth <- 0:100 # depth of core in m
  elevation_site <- site_data$Elevation
  # calculate the density profiles of the cores using DenstiyHL function from FirnR package
  rho <- DensityHL(
    depth = depth,
    rho.surface = kRho.s,
    T = t.mean,
    bdot = bdot,
    JohnsenCorr = TRUE
  )$rho
  # calculate diffusion length of d180
  sigma.d18O <- DiffusionLength(
    depth,
    rho,
    T = t.mean,
    P = press.mean,
    bdot = bdot,
    dD = FALSE,
  )
  # calculate diffusion length of dD
  sigma.dD <- DiffusionLength(
    depth,
    rho,
    T = t.mean,
    P = press.mean,
    bdot = bdot,
    dD = TRUE
  )
  
  # save diffusionlength for d18O with corresponding site and depth in dataframe
  sigma_df <- data.frame(site_i, sigma.d18O, depth, elevation)
  
  # calculated maxiumum diffusion length in core for chosen depth
  print(max(sigma_df$sigma.d18O, na.rm = TRUE))
  max_index <- which.max(sigma_df$sigma.d18O)
  depth_at_max_sigma <- sigma_df$depth[max_index]
  print(depth_at_max_sigma)
  
  # apply a depth filter, to only calculate diffusion length in set depth range
  depth_range <- sigma_df$depth >= 23.2 & sigma_df$depth <=30
  filtered_sigma_df <- sigma_df[depth_range, ]
  
  
  # calculate mean and median sigma values in depth range
  mean.sigma <- mean(filtered_sigma_df$sigma.d18O, na.rm = TRUE)
  median.sigma <- median(filtered_sigma_df$sigma.d18O, na.rm = TRUE)
  # create dataframes for the mean and medians values for sigma with corresponding site
  mean.sig.site[[site_i]] <- data.frame(site_i, mean.sigma)
  median.sig.site[[site_i]] <- data.frame(site_i, median.sigma)
  filtered_sigma_site[[site_i]] <- data.frame(site_i, filtered_sigma_df)
  
  # plot diffusion length profiles
  pdf(
    paste0("plots/", site_i, "_diffusion_profile.png"),
    width = 6,
    height = 4
  )
  # plot diffusion profile
  plot(
    sigma.d18O,
    depth,
    ylim = c(100, 0),
    xlim = c(0, 20),
    type = "l",
    xlab = "diffusion length (cm)",
    ylab = "depth (m)",
    main = "",
    lwd = 2,
    las = 1,
    col = "cadetblue"
  )
  lines(sigma.dD, depth, lwd = 2, col = "deeppink4")
  legend(
    "topleft",
    c("d18O", "dD"),
    col = c("cadetblue", "deeppink4"),
    lwd = 2,
    bty = "n"
  )
  dev.off()
  
}

# add all mean and median diffusionslengths to a dataframe
mean.calc.sig_df <- bind_rows(mean.sig.site)
median.calc.sig_df <- bind_rows(median.sig.site)
str(filtered_sigma_site)
elevation_df <- bind_rows(
  lapply(names(filtered_sigma_site), function(site_name) {
    df <- filtered_sigma_site[[site_name]]
    tibble(site_i = site_name, elevation = unique(df$elevation))
  })
)
print(elevation_df)
# add duplicate values for exNGT, so site column can be used later for different data of the same cores
dup_map <- tibble::tibble(
  original = c("B27.28", "B22", "B16", "B22", "B18", "B21", "B23", "B26", "B22"),
  new = c("exB2728", "B22-2019", "B16-2019", "cleaned_B22-2019", "B18-2012", "B21-2012", "B23-2012", "B26-2012", "n_c_B22-2019")
)

duplicates <- purrr::map2_dfr(dup_map$original, dup_map$new, function(orig, new_name) {
  mean.calc.sig_df %>%
    filter(site_i == orig) %>%
    mutate(site_i = new_name)
})

mean.calc.sig_df <- bind_rows(mean.calc.sig_df, duplicates)
print(mean.calc.sig_df) #debug



filtered_sigma_site_df <- bind_rows(filtered_sigma_site)
filtered_sigma_site_df <- filtered_sigma_site_df %>%
  filter(site_i %in% c("B16", "B18", "B21", "B23", "B26", "B22", "GRIP"))

diff_length_plot <- ggplot(filtered_sigma_site_df, aes(y = sigma.d18O, x = depth, group = site_i,  color = site_i)) +
  geom_line(size =  1) +
  labs(
    x = "Depth [m]",
    y = expression(sigma~"[cm]"),
    color = "Site"
  ) +
  #scale_x_reverse() +  # so depth increases downward
  theme_minimal() +
  scale_color_viridis_d(option = "C") +  # oder eine andere Palette
  theme(
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )

print(diff_length_plot)


ggsave(
  filename = paste0("plots/", "diffusionprofiles.pdf"),
  plot = diff_length_plot,
  width = 6,
  height = 4)

ggsave(
  filename = paste0("plots/", "diffusionprofiles.png"),
  plot = diff_length_plot,
  width = 6,
  height = 4)


# plot with elevation data

site_name_map <- c(
  "B16" = "B16-2019",
  "B18" = "B18-2012",
  "B21" = "B21-2012",
  "B22" = "B22-2019",
  "B23" = "B23-2012",
  "B26" = "B26-2012",
  "GRIP" = "GRIP"
)

# Apply renaming
filtered_sigma_site_df <- filtered_sigma_site_df %>%
  mutate(site_i = recode(as.character(site_i), !!!site_name_map))

site_labels <- filtered_sigma_site_df %>%
  group_by(site_i) %>%
  summarise(
    elevation = first(elevation),
    depth = mean(depth),
    .groups = "drop"
  ) %>%
  arrange(desc(elevation)) %>%
  mutate(label_y = seq(from = n(), to = 1, by = -1))  # Even spacing, top to bottom

# Main plot
p1 <- ggplot(filtered_sigma_site_df, aes(x = depth, y = sigma.d18O, group = site_i)) +
  geom_line(aes(color = elevation), linewidth = 1) +
  #scale_color_viridis_c(option = "viridis", name = "Elevation [m]") +
  scale_color_viridis_c(option = "viridis", name = "Elevation [m]") +
  labs(y = expression(sigma~"[cm]"), x = "Depth [m]") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
  )+
  theme(legend.position = "right")

# Label plot
p2 <- ggplot(site_labels, aes(y = label_y)) +
  geom_point(aes(x = 1, color = elevation), size = 4) +
  geom_text(aes(x = 1.1, label = site_i), color = "black", hjust = 0, size = 4) +
  scale_color_viridis_c(option = "viridis", guide = "none") +
  scale_y_continuous(breaks = site_labels$label_y, labels = NULL, expand = c(0.1, 0.1)) +
  scale_x_continuous(limits = c(0.9, 1.5)) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine
elevation_diffusivity_plot <- p1 + p2 + plot_layout(widths = c(4, 2))
print(elevation_diffusivity_plot)

ggsave(
  filename = paste0("plots/", "elevation_diffusivity.pdf"),
  plot = elevation_diffusivity_plot,
  width = 6,
  height = 4)





