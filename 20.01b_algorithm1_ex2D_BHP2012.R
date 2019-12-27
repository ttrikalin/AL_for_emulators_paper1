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
#get_sequential_design_pD <- function(simulator_info_list, 
#                                     n_consecutive_T_resample=5,
#                                     n_consecutive_T_SE=5, 
#                                     save_as_you_go_path=NULL, 
#                                     continue_from_list =NULL, # list(algorithm_trace, T_resample, T_SE),
#                                     max_iter = 100 ) 

source('./01_algorithm1_pD.r')  

#Read algorithm 2 functions 
#get_candidate_points_1D <- function(X)
source('./02_algorithm2.r')  

### Read simulator(), D, T_resample
# The Branin, or Branin-Hoo, function has three global minima; rescaled for input at [0,1]^2
# Picheny, V., Wagner, T., & Ginsbourger, D. (2012).
#         A benchmark of kriging-based infill criteria for noisy optimization
source('./benchmark_simulators/2D_simulator_BHP2012.R')


# Start from scratch 
algorithm_trace <- get_sequential_design_pD(simulator_info_list = simulator_info_list_BHP2012_2D, 
                                            save_as_you_go_path = 'interim_results_2D')


# #Example where we start from a previous algorithm trace 
# old_algorithm_trace <- readRDS(file='./interim_results/__interim_results_iteration_018.rds')
# continue_from_list <- list('algorithm_trace' = old_algorithm_trace, 
#                            'T_SE'=0.15, 
#                            'T_resample'=0.15)
# algorithm_trace <- get_sequential_design_pD(simulator_info_list = simulator_info_list_BHP2012_2D, 
#                                             continue_from_list = continue_from_list,
#                                             save_as_you_go_path = 'interim_results')




# save the algorithm_trace 
saveRDS(algorithm_trace, file='./benchmark_results/2D_BHP2012/2D_BHP2012_AL.rds')

