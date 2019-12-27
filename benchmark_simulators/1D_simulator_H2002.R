#########################################
#### Benchmark simulator (1-D)
#########################################

# Higdon 2002. Space and space-time modeling using process convolutions. 
#              In Quantitative methods for current environmental issues (pp. 37-56). Springer London.
sim_H2002 <- function(x_in_0_1) {                 
  # map x from [0,1] to [0, 10]
  x <- 10 * x_in_0_1
  y <- sin(2*pi*x/10) + 0.2*sin(2*pi*x/2.5)
  return(y)
}

# The hull of the input polytope
X_H2002<- list()
X_H2002[['h1']] <- list('h1', 0, 0)
X_H2002[['h2']] <- list('h2', 0, 1)

# add some interior points -- a seeding design 
X_H2002[['i1']] <- list('i1', 0, 0.25)
X_H2002[['i2']] <- list('i2', 0, 0.50)
X_H2002[['i3']] <- list('i3', 0, 0.75)

# The starting set of the design points, as a dataframe 
D_H2002 <- setNames(do.call(rbind.data.frame, X_H2002), c('point', 't', 'x'))
D_H2002$y <- sim_H2002(D_H2002$x)

# The application dependent threshold for resampling (~7% of range)
T_resample_H2002 <- 0.15
T_SE_H2002 <- T_resample_H2002
#########################################


simulator <- sim_H2002
simulator_name <- "1D_H2002"
simulator_name_plots <- "Higdon (2002)"
D <- D_H2002
T_resample <- T_resample_H2002
T_SE <- T_SE_H2002

input_names <- 'x'
input_dimensions <- 1


# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5

# package in an object:
simulator_info_list_H2002_1D <- list(
  'simulator' = simulator,
  'input_dimensions' = 1,
  'input_names' = 'x',
  'simulator_name' = simulator_name,
  'simulator_name_plots' = simulator_name_plots, 
  'D' = D, 
  'T_resample' = T_resample, 
  'T_SE' = T_SE, 
  'T_SE_ALM' = T_SE_ALM,
  'T_SE_ALC' = T_SE_ALC
)
