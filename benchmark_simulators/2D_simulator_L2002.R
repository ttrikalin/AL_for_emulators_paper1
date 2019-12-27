#########################################
#### Benchmark simulator (2-D)
#########################################

# get functions to populate the seeding designs 
source('./00_algorithm1_functions.R')

# Lim, Y. B., Sacks, J., Studden, W. J., & Welch, W. J. (2002). 
# Design and analysis of computer experiments when the output is 
# highly correlated over the input space. Canadian Journal of Statistics, 30(1), 109-126.

sim_L2002 <- function(x_in_0_1) {                 
  x1 <- x_in_0_1[,1] 
  x2 <- x_in_0_1[,2] 
  y <- 9 + 2.5*x1 - 17.5*x2 + 
       2.5*x1*x2 + 19*x2^2 - 
       7.5*x1^3 -
       2.5*x1*x2^2 - 5.5*x2^4 +
       x1^3 * x2^2
  return(y)
}

D_L2002 <- get_seeding_design(number_of_interior_points = 8, input_dimensions = 2)
D_L2002$t <- 0
D_L2002 <- D_L2002[,c('point', 't', 'x1', 'x2')]
D_L2002$y <- sim_L2002(D_L2002[,c('x1', 'x2')])

# The application dependent threshold for resampling (~3% of range)
T_resample_L2002 <- 0.15
T_SE_L2002 <- T_resample_L2002
#########################################

simulator <- sim_L2002
simulator_name <- "2D_L2002"
simulator_name_plots <- "Lim (2002)"
input_names <- c('x1', 'x2')
D <- D_L2002
T_resample <- T_resample_L2002
T_SE <- T_SE_L2002

# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5

# package in an object:
simulator_info_list_L2002_2D <- list(
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




