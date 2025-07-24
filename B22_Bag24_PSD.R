
#####
#' This estimates the diffusion length for the Bag24 core with isotope data obtained with CFA measurements, which include "cleaned" and "not cleaned" data - where all data that results in a loss of signal was removed

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
diffusion_length_path <- "your_filepath/IceDiffusionTools-main/diffusion_length_fit.stan"
#------------------------------------------------------
# read in isotopedata
Bag24_data_24 <- read.csv2(
  "your_filepath/Probenliste_ExNGTB22_diskrete_melttest_mbeh(ExNGTB22 Bag 24).csv",
  header = TRUE,
  row.names = NULL,
  fileEncoding = "latin1"
) # needed because of Ü comment

#names(Bag24_data) <- make.names(names(Bag24_data), unique = TRUE) # assign names 
site_name <- "Bag24_Bag24" # assign site name 

# create data frame with d018 isotope data and corresponding depth
#f_Bag24_data <- data.frame(iso = Bag24_data$`d18O....`,depth = Bag24_data$`Depth..m.`) # cleaned isotope Data
f_Bag24_data <- data.frame( depth = Bag24_data_24[["Tiefe..mittlere."]],
                           iso = Bag24_data_24[["d18O"]]) 
print(f_Bag24_data)



# Extract subset from isotope data
Bag24_iso <- f_Bag24_data$iso
Bag24_depth <- f_Bag24_data$depth

# create filter that removes NA values from data - only necessary for cleaned data
#valid_idx <- !is.na(Bag24_iso) & !is.na(Bag24_depth)
#Bag24_iso <- Bag24_iso[valid_idx] # create subset of data without NA values
#Bag24_depth <- Bag24_depth[valid_idx] # create subset of data without NA values


# evenly space data 
c.iso <- MakeEquidistant(Bag24_depth, Bag24_iso, dt = 0.026)
str(c.iso)
#print(c.iso$depth)
#print(c.iso$ts)

# compare interpolated and original data - debug
pdf(paste0("plots/", site_name, "Bag24_cfainterpolationcomp.pdf"), width = 6, height = 4)

plot(Bag24_depth, Bag24_iso, type = "l", col = "#2a9d8f", 
     xlab = "Depth [m]", ylab = "δO18 [‰]", main = "")
lines(c.iso$depth, c.iso$ts, col = "#6A4C93")
legend("topright", legend = c("original data", "interpolated data"), col = c("#2a9d8f", "#6A4C93"), 
       lty = c(1, 1), lwd = 1, cex = 0.5)

dev.off()

Bag24_iso <- c.iso$ts
Bag24_depth <- c.iso$depth

sigma = 0.1235610# define sigma prior
#-----------------------------------------------------
Bag24_freq <- freq_axis(Bag24_depth) # create frequency axis
dz <- mean(diff(c.iso$depth)) # define spacing
print(dz) # debug
Bag24_raw.spec <- raw_spec(c.iso$depth, c.iso$ts) # use raw_spec function to create psd (with detrend = TRUE)
#----------------------------------------------------------
# define variables and constants
alpha <- mean(Bag24_raw.spec$spec[3:5], na.rm = TRUE) # average the first 5 points without the slowest frequency for the inital spectrum P0 before diffusion
beta <- 0
noise <- 0.1 # noise in [‰]
#print(sigma) # debug
#---------------------------------------------------
#----------------------------------------------------------------
#fit data with diffusion_length_fit.stan for the diffusion length fit estimate
Bag24_results <- new_diffusion_length_estimate(
  model_path = diffusion_length_path,
  depth = Bag24_depth,
  rec = Bag24_iso,
  alpha_mean = alpha ,
  alpha_sd <- sd(Bag24_raw.spec$spec[3:5], na.rm = TRUE)*100, #standard deviation is the sd of the alpha value to the power of 10
  beta_mean = 0 ,
  beta_sd = 0.001 ,
  sigma_mean = sigma*2 ,
  sigma_sd = 1,
  noise_mean = noise * 1.1,
  noise_sd = noise / 5
)

saveRDS(Bag24_results, file = paste0( "results/", "Bag24_cfa", "_fit_result.rds")) # save fit results

# save fit spectral data in dataframe
Bag24_fit_spec <- data.frame(
  freq = Bag24_results$specs$freq,
  spec = Bag24_results$specs$fit,
  type = "Model Fit")


Bag24_fit.stats <- # # summarise all variables from fit
  Bag24_results$fit_stats$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975 ))
  )


# plot power spectra from the data, the fit and the theoretical expectation
Bag24_spec_plot <- gg_spec(list(raw = Bag24_raw.spec, fit = Bag24_fit_spec), time_unit = "m")
print(Bag24_spec_plot)
min_depth <- floor(min(Bag24_depth))  
max_depth <- ceiling(max(Bag24_depth))
depth_range_str <- paste0(min_depth, "to", max_depth)
print(depth_range_str)
site_name_ <- site_name
ggsave(
  filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.pdf"),
  plot = Bag24_spec_plot,
  width = 6,
  height = 4)
ggsave(
  filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.png"),
  plot = Bag24_spec_plot,
  width = 6,
  height = 4)  
# make fit results into dataframe and print
Bag24_diffusion_est_params <- bind_rows(Bag24_fit.stats,  .id = "not_cleaned_Bag24")
print(Bag24_diffusion_est_params)

# define noise from fit
Bag24_noise <- Bag24_diffusion_est_params %>% 
  filter(variable == "noise") %>% 
  pull(mean)

print(Bag24_noise)

# add all diffusion lengths with site and q5 and q95 into dataframe
Bag24_sigma_summary <- Bag24_diffusion_est_params %>%
  filter(variable == "sigma") %>%
  select( mean, q5, q95) %>%
  mutate(site = "B22 Bag24")

Bag24_sigma_est <- Bag24_sigma_summary$mean # estimated diffusion lengths
Bag24_q5 <- Bag24_sigma_summary$q5 # q5 from estimated diffusion length fit
Bag24_q25 <- Bag24_sigma_summary$q95 # q95 from estimated diffusion length fit

print(Bag24_q5)
print(Bag24_q25)
print(Bag24_sigma_est)
