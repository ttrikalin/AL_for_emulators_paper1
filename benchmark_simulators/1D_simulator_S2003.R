#########################################
#### Benchmark simulator (1-D)
#########################################

# Santner 2003. Santner, T. J., Williams, B. J., & Notz, W. I. (2003). 
#               The design and analysis of computer experiments. Springer Verlag.

sim_S2003 <- function(x_in_0_1) {                 
  # map x from [0,1] to [0, 10]
  x <- x_in_0_1
  y <- exp(-1.4*x) * cos(3.5*pi*x)
  return(y)
}

# The hull of the input polytope
X_S2003<- list()
X_S2003[['h1']] <- list('h1', 0, 0)
X_S2003[['h2']] <- list('h2', 0, 1)

# add some interior points -- a seeding design 
X_S2003[['i1']] <- list('i1', 0, 0.25)
X_S2003[['i2']] <- list('i2', 0, 0.50)
X_S2003[['i3']] <- list('i3', 0, 0.75)

# The starting set of the design points, as a dataframe 
D_S2003 <- setNames(do.call(rbind.data.frame, X_S2003), c('point', 't', 'x'))
D_S2003$y <- sim_S2003(D_S2003$x)

# The application dependent threshold for resampling (~7% of range)
T_resample_S2003 <- 0.15
T_SE_S2003 <- T_resample_S2003
#########################################

simulator <- sim_S2003
simulator_name <- "1D_S2003"
simulator_name_plots <- "Santner (2003)"
D <- D_S2003
T_resample <- T_resample_S2003
T_SE <- T_SE_S2003

# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5

input_names <- 'x'
input_dimensions <- 1 

# package in an object:
simulator_info_list_S2003_1D <- list(
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




