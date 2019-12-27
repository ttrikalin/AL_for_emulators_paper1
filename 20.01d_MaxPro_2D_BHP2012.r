### code to implement the benchmark tests in the manuscript 
library(GPfit)
library(geometry)
library(matrixStats)
library(stringr)
options(stringsAsFactors = FALSE)


### Read useful functions 
source('./00_algorithm1_functions.r')  

### Read simulator(), D, T_resample
# The Branin, or Branin-Hoo, function has three global minima; rescaled for input at [0,1]^2
# Picheny, V., Wagner, T., & Ginsbourger, D. (2012).
#         A benchmark of kriging-based infill criteria for noisy optimization
source('./benchmark_simulators/2D_simulator_BHP2012.R')


D_0 <- D
number_of_interior_points_0 <- length(grep('i', D_0$point))
T_max <- 100 


results_trace <- list()
results_trace[['D_list']] <- list()
results_trace[['f']] <- list()
results_trace[['simulator']] <- simulator
results_trace[['simulator_name']] <- simulator_name
results_trace[['simulator_name_plots']] <- simulator_name_plots

results_trace[['input_names']] <- input_names
results_trace[['input_dimensions']] <- length(input_names)
results_trace[['beta']] <- input_names

cat(paste('\nNow running MaxPro designs over ',T_max, ' steps:\n' , sep ='')) 
for (t in 1:T_max) {
  cat(paste(t, '.', sep=''))
  X <- get_seeding_design_MaxPro(number_of_interior_points = number_of_interior_points_0+t, input_dimensions = length(input_names))
  Y <- simulator(X[,input_names])
  
  results_trace[['D_list']][[t]] <- cbind(X, Y)
  if (t>1) {
    last_beta <- results_trace[['f']][[t-1]][['beta']]
    results_trace[['f']][[t]] <- GP_fit(X=X[,input_names], Y=Y, optim_start = last_beta)
  } else if(t==1) {
    results_trace[['f']][[t]] <- GP_fit(X=X[,input_names], Y=Y)
  }
}
cat('\n')


saveRDS(results_trace, file='./benchmark_results/2D_BHP2012/2D_BHP2012_SpaceFill.rds')

