#########################################
#### Benchmark simulator (2-D)
#########################################

# get functions to populate the seeding designs 
source('./00_algorithm1_functions.R')

# The Branin, or Branin-Hoo, function has three global minima; rescaled for input at [0,1]^2
# Picheny, V., Wagner, T., & Ginsbourger, D. (2012). 
#         A benchmark of kriging-based infill criteria for noisy optimization


sim_BHP2012 <- function(x_in_0_1) {                 
  # map x [0, 1]^2 as per Picheny
  x1 <- 15*x_in_0_1[,1] - 5 
  x2 <- 15*x_in_0_1[,2] 
  y <- (1/51.95) * (
          (x2 - 5.1*(x1^2)/(4*pi^2) + 5*x1/pi -6)^2 +
          (10 - 10/(8*pi))*cos(x1) - 
          44.81 )
  return(y)
}

D_BHP2012 <- get_seeding_design(number_of_interior_points = 8, input_dimensions = 2)
D_BHP2012$t <- 0
D_BHP2012 <- D_BHP2012[,c('point', 't', 'x1', 'x2')]
D_BHP2012$y <- sim_BHP2012(D_BHP2012[,c('x1', 'x2')])

# The application dependent threshold for resampling (~3% of range)
T_resample_BHP2012 <- 0.15
T_SE_BHP2012 <- T_resample_BHP2012
#########################################

simulator <- sim_BHP2012
simulator_name <- "2D_BHP2012"
simulator_name_plots <- "Branin-Hoo-Picheny (2012)"
input_names <- c('x1', 'x2')
D <- D_BHP2012
T_resample <- T_resample_BHP2012
T_SE <- T_SE_BHP2012

# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5

input_dimensions <- 2


# package in an object:
simulator_info_list_BHP2012_2D <- list(
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




