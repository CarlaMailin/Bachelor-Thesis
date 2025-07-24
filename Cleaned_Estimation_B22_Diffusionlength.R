#####
#' This estimates the diffusion length for the B22 core with isotope data obtained with CFA measurements, which include "cleaned" and "not cleaned" data - where all data that results in a loss of signal was removed

#-------------------------------------
#import packages
library(PaleoSpec)
library(dplyr)
library(cmdstanr)
library(posterior)
library(ggplot2)
source("your_filepath/IceDiffusionTools-main/all_functions.R")
source("your_filepath/paleospec-master/paleospec-master/R/gg_spec_depth.R") # ggspec with renamed axis
source("your_filepath/paleospec-master/paleospec-master/R/MakeEquidistant_ts_depth.R") # added a return of the depth in MakeEquidistant
#-------------------------------------------------------
# define path for fitting model from Ice Diffusion Tools
diffusion_length_path <- "Cyour_filepath/IceDiffusionTools-main/diffusion_length_fit.stan"
#------------------------------------------------------
# read in isotopedata
cleaned_B22_data <- read.csv("your_filepath/ExNGT_B22_Full-core_processed_cleaned_vs_NOT-cleaned_070425(2)(ExNGT_B22_039_002_online data).csv",
                     skip = 1, header = TRUE,  fileEncoding = "latin1"
                     , check.names = TRUE, stringsAsFactors = FALSE    )
#names(B22_data) <- make.names(names(B22_data), unique = TRUE) # assign names 
site_name <- "cleaned_B22_cfa" # assign site name - not to be confused with discretly measured cleaned_B22

# create data frame with d018 isotope data and corresponding depth
f_cleaned_B22_data <- data.frame(iso = cleaned_B22_data$`d18O....`,depth = cleaned_B22_data$`Depth..m.`) # cleaned isotope Data

f_cleaned_B22_data_rev <- f_cleaned_B22_data[nrow(f_cleaned_B22_data):1, ] # reverse data so it starts at the lowest depth
depth_range <-f_cleaned_B22_data_rev$depth >=23.2 &f_cleaned_B22_data_rev$depth  <= 30 # define a depth range to be used to calc. spectra

# Extract subset from isotope data
cleaned_B22_iso <- subset(f_cleaned_B22_data_rev$iso, depth_range)
cleaned_B22_depth <- subset(f_cleaned_B22_data_rev$depth, depth_range)

# create filter that removes NA values from data - only necessary for cleaned data#
valid_idx <- !is.na(cleaned_B22_iso) & !is.na(cleaned_B22_depth)
cleaned_B22_iso <- cleaned_B22_iso[valid_idx] # create subset of data without NA values
cleaned_B22_depth <- cleaned_B22_depth[valid_idx] # create subset of data without NA values


# evenly space data 
c.iso <- MakeEquidistant(cleaned_B22_depth, cleaned_B22_iso, dt = 0.026)
str(c.iso)
#print(c.iso$depth)
#print(c.iso$ts)

# compare interpolated and original data - debug
pdf(paste0("plots/", site_name, "cleaned_B22_cfainterpolationcomp.pdf"), width = 6, height = 4)

plot(cleaned_B22_depth, cleaned_B22_iso, type = "l", col = "blue", 
     xlab = "Depth [m]", ylab = "δO18 [‰]", main = "")
lines(c.iso$depth, c.iso$ts, col = "red")
legend("topright", legend = c("original data", "interpolated data"), col = c("blue", "red"), 
       lty = c(1, 1), lwd = 1, cex = 0.5)

dev.off()

cleaned_B22_iso <- c.iso$ts
cleaned_B22_depth <- c.iso$depth

sigma = 0.1235610# define sigma prior
#-----------------------------------------------------
cleaned_B22_freq <- freq_axis(cleaned_B22_depth) # create frequency axis
dz <- mean(diff(c.iso$depth)) # define spacing
print(dz) # debug
cleaned_B22_raw.spec <- raw_spec(c.iso$depth, c.iso$ts) # use raw_spec function to create psd (with detrend = TRUE)
#----------------------------------------------------------
# define variables and constants
alpha <- mean(cleaned_B22_raw.spec$spec[3:5], na.rm = TRUE) # average the first 5 points without the slowest frequency for the inital spectrum P0 before diffusion
beta <- 0
noise <- 0.1 # noise in [‰]
#print(sigma) # debug
#---------------------------------------------------
cleaned_B22_theo.iso <- diffusion_spec(cleaned_B22_freq, alpha, beta, sigma, dz = dz, noise = 0.1) # calculate theoretical PSD after diffusion
cleaned_B22_theo.spec <- list(freq = cleaned_B22_freq, spec = cleaned_B22_theo.iso)# add theoretical frequency to list with theoretical PSD
#----------------------------------------------------------------
#fit data with diffusion_length_fit.stan for the diffusion length fit estimate
cleaned_B22_results <- new_diffusion_length_estimate(
  model_path = diffusion_length_path,
  depth = cleaned_B22_depth,
  rec = cleaned_B22_iso,
  alpha_mean = alpha ,
  alpha_sd <- sd(cleaned_B22_raw.spec$spec[3:5], na.rm = TRUE)*100, #standard deviation is the sd of the alpha value to the power of 10
  beta_mean = 0 ,
  beta_sd = 0.001 ,
  sigma_mean = sigma*2 ,
  sigma_sd = 1,
  noise_mean = noise * 1.1,
  noise_sd = noise / 5
)

saveRDS(cleaned_B22_results, file = paste0( "results/", "cleaned_B22_cfa", "_fit_result.rds")) # save fit results

# save fit spectral data in dataframe
cleaned_B22_fit_spec <- data.frame(
  freq = cleaned_B22_results$specs$freq,
  spec = cleaned_B22_results$specs$fit,
  type = "Model Fit")


cleaned_B22_fit.stats <- # # summarise all variables from fit
  cleaned_B22_results$fit_stats$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975 ))
  )


# plot power spectra from the data, the fit and the theoretical expectation
cleaned_B22_spec_plot <- gg_spec(list(raw = cleaned_B22_raw.spec, fit = cleaned_B22_fit_spec, theory = cleaned_B22_theo.spec), time_unit = "m")
print(cleaned_B22_spec_plot)
# save plot
min_depth <- floor(min(cleaned_B22_depth))  
max_depth <- ceiling(max(cleaned_B22_depth))
depth_range_str <- paste0(min_depth, "to", max_depth)
print(depth_range_str)
site_name_ <- site_name
ggsave(
  filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.pdf"),
  plot = cleaned_B22_spec_plot,
  width = 6,
  height = 4)
ggsave(
  filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.png"),
  plot = cleaned_B22_spec_plot,
  width = 6,
  height = 4)  
# make fit results into dataframe and print
cleaned_B22_diffusion_est_params <- bind_rows(cleaned_B22_fit.stats,  .id = "cleaned_B22_cfa")
print(cleaned_B22_diffusion_est_params)

# define noise from fit
cleaned_B22_noise <- cleaned_B22_diffusion_est_params %>% 
  filter(variable == "noise") %>% 
  pull(mean)

print(cleaned_B22_noise)

# add all diffusion lengths with site and q5 and q95 into dataframe
cleaned_B22_sigma_summary <- cleaned_B22_diffusion_est_params %>%
  filter(variable == "sigma") %>%
  select( mean, q5, q95) %>%
  mutate(site = "cleaned_B22-2019")

cleaned_B22_sigma_est <- cleaned_B22_sigma_summary$mean # estimated diffusion lengths
cleaned_B22_q5 <- cleaned_B22_sigma_summary$q5 # q5 from estimated diffusion length fit
cleaned_B22_q25 <- cleaned_B22_sigma_summary$q95 # q95 from estimated diffusion length fit

print(cleaned_B22_q5)
print(cleaned_B22_q25)
print(cleaned_B22_sigma_est)



