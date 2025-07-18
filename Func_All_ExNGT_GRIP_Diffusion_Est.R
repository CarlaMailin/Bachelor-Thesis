###########
#' Function estimates diffusion length of d18O for exNGT-2019 cores in a given depth interval
#' Estimation and spectrum calculation methods and functions used from Shaw et al. "Novel approach to estimate the water isotope diffusion length in  deep ice cores with an application to Marine Isotope Stage 19  in the Dome C ice core" (doi.org/10.5194/tc-18-3685-2024)
#' 
#' ---------------------------------------------------------------
# define function to estimate diffusion lengths at all sites
EstSigma_exNGT <- function(data, site_name, alpha, beta, sigma) {
  cat("Processing", site_name, "\n")
  cat("theo diffusion length", sigma, "\n")
  depth_range <-  data$`depth` >=23.2 & data$`depth` <= 30# define a depth range to be used to calc. spectra
  
  # Extract subset from isotope data
  iso <- subset(data$`d18O`, depth_range)
  depth <- subset(data$`depth`, depth_range)
  cat("number of valid values:", length(iso), "\n")
  print(summary(iso))
  # create filter that removes NA values from data - only necessary for cleaned data
  valid_idx <- !is.na(iso) & !is.na(depth)
  iso <- iso[valid_idx] # create subset of data without NA values
  depth <- depth[valid_idx] 
  cat("number of valid values:", length(iso), "\n")
  print(summary(iso))
  #-----------------------------------------------------
  #-----------------------------------------------------
  # evenly space data 
  c.iso <- MakeEquidistant(depth, iso, dt = 0.026)
  iso <- c.iso$ts
  print(is.ts(iso))
  print(summary(iso))
  depth <- c.iso$depth
  freq <- freq_axis(depth) # make frequency axis
  dz <- mean(diff(c.iso$depth)) # define spacing of depth data
  #----------------------------------------------------
  # compare interpolated and original data as debug in plot
  interpcomp <- plot(depth, iso, col = "blue", pch = 20, xlab = "Depth [m]", ylab = "δO18 [‰]", main = "")
  lines(c.iso$depth, c.iso$ts, col = "red")
  legend("topright", legend = c("original data", "interpolated data"), col = c("blue", "red"), 
         lty = c(1, 1), lwd = 1, cex = 0.7)
  
  print(interpcomp)
  #--------------------------------
  raw.spec <- raw_spec(c.iso$depth, c.iso$ts)# use raw_spec function to calc. power spectrum with detrend = TRUE
  #-------------------------------------------------------------
  # define variables and constants
  alpha <- mean(raw.spec$spec[3:5], na.rm = TRUE) # average the first 5 points without the slowest frequency for the inital spectrum P0 before diffusion
  beta <- 0
  noise <- 0.1 # noise = 0.1 ‰ 
  #print(sigma) # debug
  #---------------------------------------------------------------
  theo.iso <- diffusion_spec(freq, alpha, beta, sigma, dz = dz, noise = 0.1)   # calculate theoretical PSD after diffusion
  theo.spec <- list(freq = freq, spec = theo.iso) # add theoretical frequency to list with theoretical PSD
  cat("Processing", site_name, "\n")
  #-----------------------------------------------------------------
  #fit data with diffusion_length_fit.stan for the diffusion length fit estimate
  message("used diffusion length prior: ", sigma)
  results <- new_diffusion_length_estimate(
    model_path = diffusion_length_path,
    depth = depth, # evenly spaced depth
    rec = iso, # evenly spaced isotope data
    alpha_mean = alpha ,
    alpha_sd <- sd(raw.spec$spec[3:5], na.rm = TRUE)*100 , # standard deviation is the sd of the alpha value to the power of 10
    beta_mean = 0 ,
    beta_sd = 0.001 ,
    sigma_mean = sigma*2 ,
    sigma_sd = 1, 
    noise_mean = noise * 1.1,
    noise_sd = noise / 5
  )
  
  saveRDS(results, file = paste0( "results/", site_name, "_fit_result.rds")) # save fit results
  
  # save fit spectral data in dataframe
  fit_spec <- data.frame(
    freq = results$specs$freq,
    spec = results$specs$fit,
    type = "Model Fit")
  
  
  fit.stats <- # summarise all variables from fit
    results$fit_stats$summary(
      variables = NULL,
      posterior::default_summary_measures(),
      extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975 ))
    )
  
  # plot power spectra from the data, the fit and the theoretical expecation
  spec_plot <- gg_spec(list(Experimental = raw.spec, Fit = fit_spec, Theory = theo.spec), time_unit = "m") +
    theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
  )
  print(spec_plot)
  
  # save plot incl. site name
  min_depth <- floor(min(depth))  
  max_depth <- ceiling(max(depth))
  depth_range_str <- paste0(min_depth, "to", max_depth)
  print(depth_range_str)
  site_name_ <- site_name
  ggsave(
    filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.pdf"),
    plot = spec_plot,
    width = 6,
    height = 4)
  ggsave(
    filename = paste0("plots/", site_name_, depth_range_str, "_spectrum_plot.png"),
    plot = spec_plot,
    width = 6,
    height = 4)
  
  
  cat("everything worked") # debug
  
  return(list(
    #sigma_est = sigma_est,
    fit = results,
    fit_spec = fit_spec,
    fit.stats = fit.stats
  ))
  
}







