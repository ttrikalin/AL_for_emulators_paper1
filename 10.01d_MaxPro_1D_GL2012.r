### code to implement the benchmark tests in the manuscript 
library(GPfit)
library(geometry)
library(matrixStats)
library(stringr)
options(stringsAsFactors = FALSE)


### Read useful functions 
source('./00_algorithm1_functions.r')  

### Read simulator(), D, T_resample
# Gramacy & Lee 2012. Cases for the nugget in modeling computer experiments. Statistics and Computing, 22(3), 713-22 
source('./benchmark_simulators/1D_simulator_GL2012.r')  

D_0 <- D
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
    X <- seq(0, 1, length.out = dim(D_0)[1] + t)
    Y <- simulator(X)
    
    results_trace[['D_list']][[t]] <- cbind(X, Y)
    if (t>1) {
      last_beta <- results_trace[['f']][[t-1]][['beta']]
      results_trace[['f']][[t]] <- GP_fit(X=X, Y=Y, optim_start = last_beta)
    } else if(t==1) {
      results_trace[['f']][[t]] <- GP_fit(X=X, Y=Y)
    }
}
cat('\n')


saveRDS(results_trace, file='./benchmark_results/1D_GL2012/1D_GL2012_SpaceFill.rds')

