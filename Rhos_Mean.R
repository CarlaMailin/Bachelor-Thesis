### calculate rho.s

density_data <- read.table("C:/Users/maili/Documents/Bachelorarbeit/Isotopendaten/Schaller_2016/datasets/N2E_density.tab", 
                   skip = 38, header = TRUE, sep = "\t")

# calculate mean density (7. column  = density)
mean_density <- mean(density_data[[7]], na.rm = TRUE) * 1000 # convert from cm^-3 to m^-3

print(mean_density)



