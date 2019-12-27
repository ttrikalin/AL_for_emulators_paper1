#########################################
#### Benchmark simulator (2-D)
#########################################

# get functions to populate the seeding designs 
source('./00_algorithm1_functions.R')

#Gramacy, R. B., & Lee, H. K. (2008). 
#Gaussian processes and limiting linear models. 
#Computational Statistics & Data Analysis, 53(1), 123-136.

sim_GL2008 <- function(x_in_0_1) {    
  #map [0,1]^2 -> [-2,6]^2
  
  x1 <- 8*x_in_0_1[,1] - 2
  x2 <- 8*x_in_0_1[,2] - 2
  y <- x1*exp(-x1^2 - x2^2)
  return(y)
}

D_GL2008 <- get_seeding_design(number_of_interior_points = 8, input_dimensions = 2)
D_GL2008$t <- 0
D_GL2008 <- D_GL2008[,c('point', 't', 'x1', 'x2')]
D_GL2008$y <- sim_GL2008(D_GL2008[,c('x1', 'x2')])

# The application dependent threshold for resampling (~5% of range)
T_resample_GL2008 <- 0.05
T_SE_GL2008 <- T_resample_GL2008
#########################################

simulator <- sim_GL2008
simulator_name <- "2D_GL2008"
simulator_name_plots <- "Gramacy & Lee (2008)"
input_names <- c('x1', 'x2')
D <- D_GL2008
T_resample <- T_resample_GL2008
T_SE <- T_SE_GL2008

# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5

# package in an object:
simulator_info_list_GL2008_2D <- list(
  'simulator' = simulator,
  'input_dimensions' = 2,
  'input_names' = input_names,
  'simulator_name' = simulator_name,
  'simulator_name_plots' = simulator_name_plots, 
  'D' = D, 
  'T_resample' = T_resample, 
  'T_SE' = T_SE, 
  'T_SE_ALM' = T_SE_ALM,
  'T_SE_ALC' = T_SE_ALC
)




