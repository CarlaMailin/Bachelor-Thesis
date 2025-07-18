### 
#'calculate array of estimated diffusion lengths from d18O with the EstSigma Function for the GRIP core and exNGT cores drilled in 2019 and 2012, without B19 and B27.28 with the calculated diffusionlengths as priors and using for the estimation of the diffusion length Fyntan Shaw's code and methods from "Novel approach to estimate the water isotope diffusion length in  deep ice cores with an application to Marine Isotope Stage 19  in the Dome C ice core" (doi.org/10.5194/tc-18-3685-2024)
#' The prior for the diffusion length sigma is the calculated diffusion length using functions from the exNGT package
#-------------------------------------
#import packages
library(PaleoSpec)
library(ggplot2)
library(PaleoSpec)
library(dplyr)
library(cmdstanr)
library(posterior)
source("C:/Users/maili/Documents/Bachelorarbeit/IceDiffusionTools-main/all_functions.R")
source("C:/Users/maili/Documents/Bachelorarbeit/NGT_rechnungen/Func_All_ExNGT_GRIP_Diffusion_Est.R")
source("C:/Users/maili/Documents/Bachelorarbeit/NGT_rechnungen/Calc_Firn_Diffusion_Length.R") # have to manually set depth range still for calculated one in Estimations_Firn_Diffusion_Length
source("C:/Users/maili/Documents/Bachelorarbeit/paleospec-master/paleospec-master/R/gg_spec_depth.R") # ggspec with renamed axis
source("C:/Users/maili/Documents/Bachelorarbeit/paleospec-master/paleospec-master/R/MakeEquidistant_ts_depth.R") # added a return of the depth in MakeEquidistant
#-------------------------------------------------------
# define path for fitting model from Ice Diffusion Tools
diffusion_length_path <- "C:/Users/maili/Documents/Bachelorarbeit/IceDiffusionTools-main/diffusion_length_fit.stan"
#------------------------------------------------------
# define folder path where Isotope data is
folder_path <- "C:/Users/maili/Documents/Bachelorarbeit/Isotopendaten/all_exNGT"

# List all .tab files in the folder
ALL_NGT_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
ALL_NGT_file_titles <- tools::file_path_sans_ext(basename(ALL_NGT_files)) # define names of files that belong to the cores

# convert datasets to dataframes and add all dataframes to one list
data_list_ALL_NGT <- Map(function(file, title) {
  isodata_df <- read.csv(file, skip = 0, header = TRUE, row.names = NULL)
  return(isodata_df)
}, ALL_NGT_files, ALL_NGT_file_titles)

names(data_list_ALL_NGT) <- ALL_NGT_file_titles # add corresponding file name to core in data list
#ALL_NGT_site_names <- c("exB16", "exB18", "exB21", "exB22", "exB23", "exB26", "GRIP") # define site names 
print(names(data_list_ALL_NGT))
ALL_NGT_site_names <- c("B16-2019", "B18-2012", "B21-2012", "B23-2012","B26-2012", "GRIP") # double checked with data heads, that that is correct
names(data_list_ALL_NGT) <- ALL_NGT_site_names # Assign site names
str(data_list_ALL_NGT)
names(data_list_ALL_NGT)


#---------------------------------------------------------------------
# define calculated diffusionlengths to use later as a prior
sigma <- mean.calc.sig_df$mean/ 100  # convert from cm to m
names(sigma) <- mean.calc.sig_df$site_i # add site names to prior diffusionlength
print(sigma)

# loop for estimation of sites using the calc. sigma as a prior for the estimated diffusion length value
results_all <- mapply(
  EstSigma_exNGT,
  data = data_list_ALL_NGT,
  site_name = names(data_list_ALL_NGT), 
  sigma = sigma[names(data_list_ALL_NGT)], # reorders sigma (calc. value as prior) to match sitenames
  SIMPLIFY = FALSE
)
#str(results_all)

# sort results with site name and estimated sigma
results_named_ALL_NGT <- lapply(names(results_all), function(site) {
  fit.stats <- results_all[[site]]$fit.stats
  fit.stats$site <- site
  return(fit.stats)
})

names(results_named_ALL_NGT) <- names(results_all) # sync names 

# make into dataframe and print
diffusion_est_params_ALL_NGT <- bind_rows(results_named_ALL_NGT, .id = "site")
print(diffusion_est_params_ALL_NGT)

noise_ALL_NGT <- diffusion_est_params_ALL_NGT %>% 
  filter(variable == "noise") %>% 
  pull(mean)
print(noise_ALL_NGT)

# add all diffusion lengths with site and q5 and q95 into dataframe
sigma_summary_ALL_NGT <- diffusion_est_params_ALL_NGT %>%
  filter(variable == "sigma") %>%
  select(site, mean, q5, q95) %>%
  arrange(site)

sigma_est_ALL_NGT <- sigma_summary_ALL_NGT$mean # estimated diffusion lengths
q5_ALL_NGT <- sigma_summary_ALL_NGT$q5 # q5 from estimated diffusion length fit
q25_ALL_NGT <- sigma_summary_ALL_NGT$q95 # q95 from estimated diffusion length fit

print(sigma_summary_ALL_NGT)
#debug
print(sigma[names(data_list_ALL_NGT)])
#










