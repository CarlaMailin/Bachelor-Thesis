# Evaluation of B22 CFA Data in comparison to B22 discrete data (1 cm) resolution measured in April 2025 for Bag 16, 17 and 24, which all contain melt layers, which were identified, measured and documented by Samira Zander
#'
#-------------------------------------
#import packages
library(dplyr)
library(ggplot2)
#-------------------------------------------------------
# define depths of melt layers
melt_layers_16_17 <- data.frame(pos = c(16.0955, 16.0765, 16.009, 15.988))
melt_layers_24 <- data.frame(pos = c(23.314, 23.174))
#--------------------------------------------------------------------
# read in CFA Data from B22 core
B22_cfa_data <- read.csv(
  "your_filepath//ExNGT_B22_Full-core_processed_cleaned_vs_NOT-cleaned_070425(2)(ExNGT_B22_039_002_online data).csv",
  skip = 1,
  header = TRUE,
  fileEncoding = "latin1"
  ,
  check.names = TRUE,
  stringsAsFactors = FALSE
)
#names(B22_data) <- make.names(names(B22_data), unique = TRUE) # assign names
site_name <- "B22_cfa" # assign site name - not to be confused with discretly measured B22

# create data frame with d018 isotope data and corresponding depth
f_B22_cleaned_data <- data.frame(iso = B22_cfa_data$`d18O....`, depth = B22_cfa_data$`Depth..m.`) # cleaned isotope Data
f_B22_data <- data.frame(iso = B22_cfa_data$`d18O.....1`, depth = B22_cfa_data$`Depth..m..1`) # not cleaned isotope Data

# reverse data so it starts at the lowest depth
f_B22_data_rev <- f_B22_data[nrow(f_B22_data):1, ]
# create depth interval to calculate psd over a specific depth for Bag 16/17
f_B22_data_rev1617 <- f_B22_data_rev[f_B22_data_rev$depth >= 15 &
                                       f_B22_data_rev$depth <= 17, ]
# create depth interval to calculate psd over a specific depth for Bag 24
f_B22_data_rev24 <- f_B22_data_rev[f_B22_data_rev$depth >= 23 &
                                     f_B22_data_rev$depth <= 24, ]
# reverse data so it starts at the lowest depth
f_B22_cleaned_data_rev <- f_B22_cleaned_data[nrow(f_B22_cleaned_data):1, ]
# create depth interval to calculate psd over a specific depth for Bag 16/17
f_B22_cleaned_data_rev1617 <- f_B22_cleaned_data_rev[f_B22_cleaned_data_rev$depth >= 15 &
                                                       f_B22_data_rev$depth <= 17, ]
# create depth interval to calculate psd over a specific depth for Bag 24
f_B22_cleaned_data_rev24 <- f_B22_cleaned_data_rev[f_B22_cleaned_data_rev$depth >= 23 &
                                                     f_B22_data_rev$depth <= 24, ]

# extract d018 isotope and depth data for Bag 16/17
B22_iso_1617 <- f_B22_data_rev1617$iso
B22_depth_1617 <- f_B22_data_rev1617$depth
# extract d018 isotope and depth data for Bag 24
B22_iso_24 <- f_B22_data_rev24$iso
B22_depth_24 <- f_B22_data_rev24$depth
# extract d018 isotope and data Bag 16/17
B22_cleaned_iso_1617 <- f_B22_cleaned_data_rev1617$iso
B22_cleaned_depth_1617 <- f_B22_cleaned_data_rev1617$depth
# extract d018 isotope and depth data Bag 24
B22_cleaned_iso_24 <- f_B22_cleaned_data_rev24$iso
B22_cleaned_depth_24 <- f_B22_cleaned_data_rev24$depth

# create dataframe for Bags 16/17 and Bag 24 with both cleaned and uncleaned data
cfa_24 <- data.frame(
  depth_c = B22_cleaned_depth_24,
  d18O_c = B22_cleaned_iso_24,
  depth = B22_depth_24,
  d18O = B22_iso_24,
  bag = "Bag24"
)


cfa_1617 <- data.frame(
  depth_c = B22_cleaned_depth_1617,
  d18O_c = B22_cleaned_iso_1617,
  depth = B22_depth_1617,
  d18O = B22_iso_1617,
  bag = "Bag16/17"
)

#------------------------------------------------------------------------------
# read in discrete B22 data

# read in Bag 16
B22_data_16 <- read.csv2(
  "C:/Users/maili/Documents/Bachelorarbeit/Isotopendaten/B22_1cm/Probenliste_ExNGTB22_diskrete_melttest_mbeh(ExNGTB22 Bag 16).csv",
  header = TRUE,
  row.names = NULL
)

# read in Bag 17
B22_data_17 <- read.csv2(
  "your_filepath/Probenliste_ExNGTB22_diskrete_melttest_mbeh(ExNGTB22 Bag 17).csv",
  header = TRUE,
  row.names = NULL
)

# read in Bag 24
B22_data_24 <- read.csv2(
  "your_filepath/Probenliste_ExNGTB22_diskrete_melttest_mbeh(ExNGTB22 Bag 24).csv",
  header = TRUE,
  row.names = NULL,
  fileEncoding = "latin1"
) # needed because of Ü comment

# add needed paramters from Bags into dataframes
df_16 <- data.frame(
  depth = B22_data_16[["Tiefe..mittlere."]],
  d18O  = B22_data_16[["d18O"]],
  dD    = B22_data_16[["dD"]],
  bag      = rep("Bag16", nrow(B22_data_16))
)

df_17 <- data.frame(
  depth = B22_data_17[["Tiefe..mittlere."]],
  d18O  = B22_data_17[["d18O"]],
  dD   = B22_data_17[["dD"]],
  bag      = rep("Bag17", nrow(B22_data_17))
)

df_24 <- data.frame(
  depth = B22_data_24[["Tiefe..mittlere."]],
  d18O = B22_data_24[["d18O"]],
  dD   = B22_data_24[["dD"]],
  bag      = rep("Bag24", nrow(B22_data_24))
)

# min depth as offset for depth assignment
off_16 <- 15
off_17 <- 16
off_24 <- 23
df_16 <- df_16 %>% mutate(depth_off = depth + off_16)
df_17 <- df_17 %>% mutate(depth_off = depth + off_17)
df_24 <- df_24 %>% mutate(depth_off = depth + off_24)

# calc max and min values of depth based time series                           
max(df_24$d18O)
min(df_24$d18O)
max(df_16$d18O)
min(df_16$d18O)
max(df_17$d18O)
min(df_17$d18O)

df_16_17 <- bind_rows(df_16, df_17)
#-----------------------------------------------------------
# create combined data frames for each Bag with discrete and CFA Data
bag24_combined <- bind_rows(
  df_24 %>% 
    select(depth = depth_off, d18O) %>% 
    mutate(source = "Bag24", type = "Discrete"),
  
  cfa_24 %>% 
    select(depth = depth_c, d18O = d18O_c) %>% 
    mutate(source = "Bag24", type = "CFA (cleaned)"),
  
  cfa_24 %>% 
    select(depth = depth, d18O = d18O) %>% 
    mutate(source = "Bag24", type = "CFA"))


bag1617_combined <- bind_rows(
  df_16_17 %>% 
    select(depth = depth_off, d18O) %>% 
    mutate(source = "Bag16/17", type = "Discrete"),
  
  cfa_1617 %>% 
    select(depth = depth_c, d18O = d18O_c) %>% 
    mutate(source = "Bag16/17", type = "CFA (cleaned)"),
  
  cfa_1617 %>% 
    select(depth = depth, d18O = d18O) %>% 
    mutate(source = "Bag16/17", type = "CFA")
)
#---------------------------------------------------------------
# plot results
bag24_plot <- ggplot(bag24_combined, aes(x = depth, y = d18O, color = type)) +
  geom_point(data = subset(bag24_combined, type == "Discrete"), size = 3, shape = 3) +
  geom_line(data = subset(bag24_combined, type == "CFA (cleaned)"), linewidth = 0.5) +
  geom_line(data = subset(bag24_combined, type == "CFA"), linewidth = 0.5, linetype = "dashed") +
  geom_vline(data = melt_layers_24, aes(xintercept = pos), 
             color = "black", linetype = "dotted", linewidth = 0.5) +
  scale_color_viridis_d(option = "C", begin = 0, end = 0.85) +
  labs(
    x = "Depth [m]",
    y = expression(delta^18*O~"[\u2030]"),
   color = "Data Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.7),  # x, y coordinates between 0 and 1
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),  
    legend.box.background = element_rect(color = "black"), 
    legend.text = element_text(size = 14),       
    legend.title = element_text(size = 14),
    text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.key.size = unit(0.8, "lines")        
  
  )


bag1617_plot <- ggplot(bag1617_combined, aes(x = depth, y = d18O, color = type)) +
  geom_point(data = subset(bag1617_combined, type == "Discrete"), size = 3, shape = 3) +
  geom_line(data = subset(bag1617_combined, type == "CFA (cleaned)"), linewidth = 0.5) +
  geom_line(data = subset(bag1617_combined, type == "CFA"), linewidth = 0.5, linetype = "dashed") +
  geom_vline(data = melt_layers_16_17, aes(xintercept = pos), 
             color = "black", linetype = "dotted", linewidth = 0.5) +
  scale_color_viridis_d(option = "C", begin = 0, end = 0.85) +
  labs(
    x = "Depth [m]",
    y = expression(delta^18*O~"[\u2030]"),
    color = "Data Type"
  ) +
  theme_minimal()+
  theme(
    legend.position = c(0.85, 0.7), 
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),  
    legend.box.background = element_rect(color = "black"),
    legend.text = element_text(size = 14),       
    legend.title = element_text(size = 14),   
    text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.key.size = unit(0.8, "lines")         
  )
print(bag1617_plot)
print(bag24_plot)
#----------------------------------------------------------------------
# save plots as pdf and png
ggsave(
  filename = paste0("plots/", "B22_Bag1617timeseries.pdf"),
  plot = bag1617_plot,
  width = 6,
  height = 4
)
ggsave(
  filename = paste0("plots/", "B22_Bag24timeseries.pdf"),
  plot = bag24_plot,
  width = 6,
  height = 4
)


ggsave(
  filename = paste0("plots/", "B22_Bag1617timeseries.png"),
  plot = bag1617_plot,
  width = 6,
  height = 4
)
ggsave(
  filename = paste0("plots/", "B22_Bag24timeseries.png"),
  plot = bag24_plot,
  width = 6,
  height = 4
)







