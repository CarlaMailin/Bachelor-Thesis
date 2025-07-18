#####
#' This estimates the diffusionlength for the B22 core with isotope data obtained with CFA measurements, which include "cleaned" and "not cleaned" data - where all data that results in a loss of signal was removed

#-------------------------------------
#import packages
library(PaleoSpec)
library(dplyr)
library(cmdstanr)
library(posterior)
library(ggplot2)
source("C:/Users/maili/Documents/Bachelorarbeit/IceDiffusionTools-main/all_functions.R")
source("C:/Users/maili/Documents/Bachelorarbeit/paleospec-master/paleospec-master/R/gg_spec_depth.R") # ggspec with renamed axis
source("C:/Users/maili/Documents/Bachelorarbeit/paleospec-master/paleospec-master/R/MakeEquidistant_ts_depth.R") # added a return of the depth in MakeEquidistant
#-------------------------------------------------------
# define path for fitting model from Ice Diffusion Tools
diffusion_length_path <- "C:/Users/maili/Documents/Bachelorarbeit/IceDiffusionTools-main/diffusion_length_fit.stan"
#------------------------------------------------------
# read in isotopedata
B22_data <- read.csv("C:\\Users\\maili\\Documents\\Bachelorarbeit\\Isotopendaten\\B22_CFA\\ExNGT_B22_Full-core_processed_cleaned_vs_NOT-cleaned_070425(2)(ExNGT_B22_039_002_online data).csv",
                     skip = 1, header = TRUE,  fileEncoding = "latin1"
                     , check.names = TRUE, stringsAsFactors = FALSE    )
#names(B22_data) <- make.names(names(B22_data), unique = TRUE) # assign names 
site_name <- "not_cleaned_B22" # assign site name - not to be confused with discretly measured B22

# create data frame with d018 isotope data and corresponding depth
#f_B22_data <- data.frame(iso = B22_data$`d18O....`,depth = B22_data$`Depth..m.`) # cleaned isotope Data
f_B22_data <- data.frame(iso = B22_data$`d18O.....1`,depth = B22_data$`Depth..m..1`) # not cleaned isotope Data


f_B22_data_rev <- f_B22_data[nrow(f_B22_data):1, ] # reverse data so it starts at the lowest depth
# calculate psd over a specific depth
depth_range <-f_B22_data_rev$depth >=23.2 &f_B22_data_rev$depth  <= 30 # define a depth range to be used to calc. spectra

# Extract subset from isotope data
B22_iso <- subset(f_B22_data_rev$iso, depth_range)
B22_depth <- subset(f_B22_data_rev$depth, depth_range)

# create filter that removes NA values from data - only necessary for cleaned data
#valid_idx <- !is.na(B22_iso) & !is.na(B22_depth)
#B22_iso <- B22_iso[valid_idx] # create subset of data without NA values
#B22_depth <- B22_depth[valid_idx] # create subset of data without NA values


# evenly space data 
c.iso <- MakeEquidistant(B22_depth, B22_iso, dt = 0.026)
str(c.iso)
#print(c.iso$depth)
#print(c.iso$ts)

# compare interpolated and original data - debug
pdf(paste0("plots/", site_name, "B22_cfainterpolationcomp.pdf"), width = 6, height = 4)

plot(B22_depth, B22_iso, type = "l", col = "#2a9d8f", 
     xlab = "Depth [m]", ylab = "δO18 [‰]", main = "")
lines(c.iso$depth, c.iso$ts, col = "#6A4C93")
legend("topright", legend = c("original data", "interpolated data"), col = c("#2a9d8f", "#6A4C93"), 
       lty = c(1, 1), lwd = 1, cex = 0.5)

dev.off()

B22_iso <- c.iso$ts
B22_depth <- c.iso$depth

sigma = 0.1235610# define sigma prior
#-----------------------------------------------------
B22_freq <- freq_axis(B22_depth) # create frequency axis
dz <- mean(diff(c.iso$depth)) # define spacing
print(dz) # debug
B22_raw.spec <- raw_spec(c.iso$depth, c.iso$ts) # use raw_spec function to create psd (with detrend = TRUE)
#----------------------------------------------------------
# define variables and constants
alpha <- mean(B22_raw.spec$spec[3:5], na.rm = TRUE) # average the first 5 points without the slowest frequency for the inital spectrum P0 before diffusion
beta <- 0
noise <- 0.1 # noise in [‰]
#print(sigma) # debug
#---------------------------------------------------
B22_theo.iso <- diffusion_spec(B22_freq, alpha, beta, sigma, dz = dz, noise = 0.1) # calculate theoretical PSD after diffusion
B22_theo.spec <- list(freq = B22_freq, spec = B22_theo.iso)# add theoretical frequency to list with theoretical PSD
#----------------------------------------------------------------
#fit data with diffusion_length_fit.stan for the diffusion length fit estimate
B22_results <- new_diffusion_length_estimate(
  model_path = diffusion_length_path,
  depth = B22_depth,
  rec = B22_iso,
  alpha_mean = alpha ,
  alpha_sd <- sd(B22_raw.spec$spec[3:5], na.rm = TRUE)*100, #standard deviation is the sd of the alpha value to the power of 10
  beta_mean = 0 ,
  beta_sd = 0.001 ,
  sigma_mean = sigma*2 ,
  sigma_sd = 1,
  noise_mean = noise * 1.1,
  noise_sd = noise / 5
)

saveRDS(B22_results, file = paste0( "results/", "B22_cfa", "_fit_result.rds")) # save fit results

# save fit spectral data in dataframe
B22_fit_spec <- data.frame(
  freq = B22_results$specs$freq,
  spec = B22_results$specs$fit,
  type = "Model Fit")


B22_fit.stats <- # # summarise all variables from fit
  B22_results$fit_stats$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975 ))
  )


# plot power spectra from the data, the fit and the theoretical expectation
B22_spec_plot <- gg_spec(list(raw = B22_raw.spec, fit = B22_fit_spec, theory = B22_theo.spec), time_unit = "m")
print(B22_spec_plot)
min_depth <- floor(min(B22_depth))  
max_depth <- ceiling(max(B22_depth))
depth_range_str <- paste0(min_depth, "to", max_depth)
print(depth_range_str)
site_name_ <- site_name
ggsave(
  filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.pdf"),
  plot = B22_spec_plot,
  width = 6,
  height = 4)
ggsave(
  filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.png"),
  plot = B22_spec_plot,
  width = 6,
  height = 4)  
# make fit results into dataframe and print
B22_diffusion_est_params <- bind_rows(B22_fit.stats,  .id = "not_cleaned_B22")
print(B22_diffusion_est_params)

# define noise from fit
B22_noise <- B22_diffusion_est_params %>% 
  filter(variable == "noise") %>% 
  pull(mean)

print(B22_noise)

# add all diffusion lengths with site and q5 and q95 into dataframe
B22_sigma_summary <- B22_diffusion_est_params %>%
  filter(variable == "sigma") %>%
  select( mean, q5, q95) %>%
  mutate(site = "n_c_B22-2019")

B22_sigma_est <- B22_sigma_summary$mean # estimated diffusion lengths
B22_q5 <- B22_sigma_summary$q5 # q5 from estimated diffusion length fit
B22_q25 <- B22_sigma_summary$q95 # q95 from estimated diffusion length fit

print(B22_q5)
print(B22_q25)
print(B22_sigma_est)


# add cleaned and uncleaned data in one plot - requires to run Cleaned_Estimation_B22_Diffusionlength.R as well beforehand

# souce("C:/Users/maili/Documents/Bachelorarbeit/NGT_rechnungen/Cleaned_Estimation_B22_Diffusionlength.R")
# source("C:/Users/maili/Documents/Bachelorarbeit/NGT_rechnungen/B22_Bag24_PSD.R")

B22_colors <- c(
  "Experimental n.c." = "#1b9e77",
  "Experimental clean" = "#d95f02",
  "Fit n.c." = "#7570b3",
  "Fit clean" = "#e7298a",
  "Theory" = "#66a61e",
  "Discrete" = "#e6ab02"
)


B22_comp_spec_plot <- gg_spec(list(
  "Experimental n.c." = B22_raw.spec,
  "Experimental clean" = cleaned_B22_raw.spec,
  "Fit n.c." = B22_fit_spec,
  "Theory" = B22_theo.spec,
  "Fit clean" = cleaned_B22_fit_spec
), time_unit = "m") +
  scale_color_manual(values = B22_colors) +
  labs(color = NULL) +  
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
  )


print(B22_comp_spec_plot)

# save plot
ggsave(
  filename = paste0("plots/", "Z_spectrum_plot.pdf"),
  plot = B22_comp_spec_plot,
  width = 6,
  height = 4)

B22_comp_spec_plot_Bag24 <- gg_spec(list(
  "Experimental n.c." = B22_raw.spec,
  "Experimental clean" = cleaned_B22_raw.spec,
  "Discrete" = Bag24_raw.spec
), time_unit = "m") +
  scale_color_manual(values = B22_colors) +
  labs(color = NULL) +  
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
  )


print(B22_comp_spec_plot_Bag24)

# save plot
ggsave(
  filename = paste0("plots/", "Z_spectrum_plot.pdf"),
  plot = B22_comp_spec_plot,
  width = 6,
  height = 4)

ggsave(
  filename = paste0("plots/", "Z24_spectrum_plot.pdf"),
  plot = B22_comp_spec_plot_Bag24,
  width = 6,
  height = 4)


