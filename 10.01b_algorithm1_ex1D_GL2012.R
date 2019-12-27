### code to implement the benchmark tests in the manuscript 
library(GPfit)
library(geometry)
library(matrixStats)
library(stringr)
options(stringsAsFactors = FALSE)


### Read useful functions 
#append_jk_emulators <- function(Design, f_list, current_beta=NULL)
#get_predictions_at_C <- function(f_list, candidate_points)
#get_prediction_range_at_C <- function(predictions) 
#add_design_point <- function(Design, new_x, iteration, simulator)
source('./00_algorithm1_functions.r')  

#Read algorithm 1 function 
#get_sequential_design_1D <- function(simulator_info_list, 
#                                     n_consecutive_T_resample=5,
#                                     n_consecutive_T_SE=5, 
#                                     save_as_you_go_path=NULL, 
#                                     max_iter = 100 ) 
source('./01_algorithm1_1D.r')  

#Read algorithm 2 functions 
#get_candidate_points_1D <- function(X)
source('./02_algorithm2.r')  

### Read simulator(), D, T_resample
# Gramacy & Lee 2012. Cases for the nugget in modeling computer experiments. Statistics and Computing, 22(3), 713-22 
source('./benchmark_simulators/1D_simulator_GL2012.r')  

algorithm_trace <- get_sequential_design_1D(simulator_info_list = simulator_info_list_GL2012_1D, 
                                            save_as_you_go_path = 'interim_results')

# save the algorithm_trace 
saveRDS(algorithm_trace, file='./benchmark_results/1D_GL2012/1D_GL2012_AL.rds')

