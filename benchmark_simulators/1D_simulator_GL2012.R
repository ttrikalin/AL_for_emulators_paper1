#########################################
#### Benchmark simulator (1-D)
#########################################

# Gramacy & Lee 2012. Cases for the nugget in modeling computer experiments. Statistics and Computing, 22(3), 713-22 
sim_GL2012 <- function(x_in_0_1) {                 
  # map x from [0,1] to [0.5, 2.5]
  x <- 2* x_in_0_1 + 0.5 
  y <- sin(10*pi*x)/(2*x) + (x-1)^4 
  return(y)
}

# The hull of the input polytope
X_GL2012 <- list()
X_GL2012[['h1']] <- list('h1', 0, 0)
X_GL2012[['h2']] <- list('h2', 0, 1)

# add some interior points -- a seeding design
X_GL2012[['i1']] <- list('i1', 0, 0.25)
X_GL2012[['i2']] <- list('i2', 0, 0.50)
X_GL2012[['i3']] <- list('i3', 0, 0.75)

# The starting set of the design points, as a dataframe
D_GL2012 <- setNames(do.call(rbind.data.frame, X_GL2012), c('point', 't', 'x'))
D_GL2012$y <- sim_GL2012(D_GL2012$x)

# The application dependent threshold for resampling (~7.5% of range)
T_resample_GL2012 <- 0.45
T_SE_GL2012 <- T_resample_GL2012
#########################################


simulator <- sim_GL2012
simulator_name <- "1D_GL2012"
simulator_name_plots <- "Gramacy & Lee (2012)"
input_names <- c('x')
D <- D_GL2012
T_resample <- T_resample_GL2012
T_SE <- T_SE_GL2012

# For the Gramacy BTGP algorithm 
T_SE_ALM <- T_SE
T_SE_ALC <- 10^-5



# package in an object:
simulator_info_list_GL2012_1D <- list(
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
