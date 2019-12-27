#########################################
#### Benchmark simulator (3-D)
#########################################

# get functions to populate the seeding designs 
source('./00_algorithm1_functions.R')

# Dette, H., & Pepelyshev, A. (2010). Generalized Latin 
# hypercube design for computer experiments. Technometrics, 52(4).


sim_DP2010 <- function(x_in_0_1) {                 
  x1 <- x_in_0_1[,1] 
  x2 <- x_in_0_1[,2] 
  x3 <- x_in_0_1[,3] 
  y <- 2.463019  * (exp(-2/(x1^1.75)) + exp(-2/(x2^1.5)) + exp(-2/(x3^1.25)))  # normalized between 0 and 1
  return(y)
}

D_DP2010 <- get_seeding_design(number_of_interior_points = 8, input_dimensions = 3)
D_DP2010$t <- 0
D_DP2010 <- D_DP2010[,c('point', 't', 'x1', 'x2', 'x3')]
D_DP2010$y <- sim_DP2010(D_DP2010[,c('x1', 'x2', 'x3')])

# The application dependent threshold for resampling (~5% of range)
T_resample_DP2010 <- 0.03
T_SE_DP2010 <- T_resample_DP2010
#########################################

simulator <- sim_DP2010
simulator_name <- "3D_DP2010"
simulator_name_plots <- "Dette & Pepelyshev (2010)"
input_names <- c('x1', 'x2', 'x3')
D <- D_DP2010
T_resample <- T_resample_DP2010
T_SE <- T_SE_DP2010

# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5

# package in an object:
simulator_info_list_DP2010_3D <- list(
  'simulator' = simulator,
  'input_dimensions' = 3,
  'input_names' = input_names,
  'simulator_name' = simulator_name,
  'simulator_name_plots' = simulator_name_plots, 
  'D' = D, 
  'T_resample' = T_resample, 
  'T_SE' = T_SE, 
  'T_SE_ALM' = T_SE_ALM,
  'T_SE_ALC' = T_SE_ALC
)




