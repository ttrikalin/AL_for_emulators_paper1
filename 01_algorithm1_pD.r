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

#Read algorithm 2 functions 
#get_candidate_points_1D <- function(X)
source('./02_algorithm2.r')  



### Example Usage 
# Santner 2003. Santner, T. J., Williams, B. J., & Notz, W. I. (2003). 
#               The design and analysis of computer experiments. Springer Verlag.
#source('./benchmark_simulators/1D_simulator_S2003.r')  
#my_trace <- get_sequential_design_1D(simulator_info_list = simulator_info_list_S2003_1D, 
#                                     save_as_you_go_path = 'interim_results')

get_sequential_design_pD <- function(simulator_info_list, 
                                     n_consecutive_T_resample=5,
                                     n_consecutive_T_SE=5, 
                                     save_as_you_go_path=NULL, 
                                     continue_from_list =NULL, # list(algorithm_trace, T_resample, T_SE)
                                     max_iter = 100 ) {

  if(is.null(continue_from_list)) {
        ### Set up the lists we need to trace the whole history of Algorithm 1 
        f <- list()                # list of lists of emulators developed at the end of each iteration  
        C <- list()                # list of sets of candidate points in ieach iteration
        Predicted_at_C <- list()   # list of predictions with each emulator for candidate points 
        SE_predicted_at_C <- list()   # list of SE of predictions the the current emulator [all] for candidate points 
        Prediction_range_at_C <- list() # list of ranges of prediction with each emulator for candidate points 
        largest_prediction_range <- list()  # list of the largest range of predictions at candidate points 
        largest_prediction_SE <- list()  # list of the largest prediction SE at candidate points with current emulator
        C_star1_index <-list()      # list of the index of the candidate point with the maximum prediction range
        C_star1 <- list()           # list of candidate points with the maximum prediction range
        C_star2_index <-list()      # list of the index of the candidate point with the maximum predicted SE
        C_star2 <- list()           # list of candidate points with the maximum predicted SE
        C_star_index <-list()      # list of the index of the selected candidate points (either C_star1_index or C_star2_index)
        C_star <- list()           # list of the selected candidate points (either C_star1 or C_star2)
        meets_T_resample <- list()  # was the T_resample met for the required consecutive iterations?
        meets_T_SE <- list()        #  was the T_SE met for the required consecutive iterations?
        
        t <- 0                      # initialize the iteration counter (will be overwritten if we continue from a previous algorithm_trace )
        go_on <- TRUE
        
        ### unpack the information on the simulator. E.g., for S2003_1D, and store it in the 
        ### return object algorithm_trace 
        # simulator_info_list_S2003_1D <- list(
        #   'simulator' = simulator,
        #   'input_dimensions' = 1,
        #   'simulator_name' = simulator_name,
        #   'simulator_name_plots' = simulator_name_plots, 
        #   'D' = D, 
        #   'T_resample' = T_resample, 
        #   'T_SE' = T_SE
        # )
        
        algorithm_trace <- list()
        for (var_name in names(simulator_info_list)) {
          assign(var_name, simulator_info_list[[var_name]])
          algorithm_trace[[var_name]] <- simulator_info_list[[var_name]]
        }
  } else {
    # record newly-specified options (the var-names will be overwritten by the contents of algorithm_trace)
    new_n_consecutive_T_resample <- n_consecutive_T_resample
    new_n_consecutive_T_SE <-  n_consecutive_T_SE
    new_max_iter <- max_iter 
    ### Over-write the algorithm_trace and various lists 
    cat('"continue_from_list" option was specified; examining whether stoping criteria are more stringent\n')
    algorithm_trace <- continue_from_list[['algorithm_trace']] 
    for (var_name in names(algorithm_trace)) {
      assign(var_name, algorithm_trace[[var_name]])
    }
    # restore option values 
    n_consecutive_T_resample <- new_n_consecutive_T_resample
    n_consecutive_T_SE <-  new_n_consecutive_T_SE
    max_iter <- new_max_iter 

    #use any newly provided T_SE or T_resample 
    tmp_T_resample <- continue_from_list[['T_resample']]
    if(!is.null(tmp_T_resample)) { T_resample <- tmp_T_resample }
    tmp_T_SE <- continue_from_list[['T_SE']]
    if(!is.null(tmp_T_SE)) { T_SE <- tmp_T_SE }
    
    # add list definitions that are not packaged in algorithm_trace
    SE_predicted_at_C <- list()   # list of SE of predictions the the current emulator [all] for candidate points 
    meets_T_resample <- list()  # was the T_resample met for the required consecutive iterations?
    meets_T_SE <- list()        #  was the T_SE met for the required consecutive iterations?
    
    # check that we have tightened some criterion: T_resample, T_SE, n_consecutive_T_resample, n_consecutive_T_SE
    is_T_resample_more_stringent <-  !is_T_resample_met(cutoff=T_resample, 
                                                       largest_prediction_range=largest_prediction_range, 
                                                       n_consecutive=n_consecutive_T_resample)
    cat(paste('is_T_resample_more_stringent = ', is_T_resample_more_stringent, '\n',sep=''))
    
    is_T_SE_more_stringent <-  !is_T_SE_met(cutoff=T_SE, 
                                           largest_prediction_SE=largest_prediction_SE, 
                                           n_consecutive=n_consecutive_T_SE)
    cat(paste('is_T_SE_more_stringent = ', is_T_SE_more_stringent, '\n',sep=''))
    
    if(is_T_resample_more_stringent | is_T_SE_more_stringent) { 
        cat(paste('Convergence criteria are not met; continuing from iteration ', t, '\n', sep=''))  
        go_on <- TRUE 
    } else {
      cat('You need stricter convergence criteria\n')
      cat('Returning the algorithm trace from which I was asked to continue\n')
      return(algorithm_trace)
    }
  }
  
  
  algorithm_trace[['n_consecutive_T_resample']] <- n_consecutive_T_resample
  algorithm_trace[['n_consecutive_T_SE']] <- n_consecutive_T_SE
  algorithm_trace[['max_iter']] <- max_iter
  
  ### start the list that will be returned:  
  

  while(go_on==TRUE & t <= max_iter ) {
    t<- t+1
    cat(paste("Now running iteration ", t, '\n', sep=''))
    
    #### fit an emulator to set of design points D
    f[[t]] <- list()
    f[[t]][['all']] <- GP_fit(X=D[,input_names], Y=D$y)
    current_beta <- f[[t]][['all']]$beta 
    
    #### Append jack-knifed emulators by dropping interior points 
    #### f[[t]][['all']]: f is fit on all points 
    #### f[[t]][['i1']]: f_(-i1) is fit dropping point i1 and so on. 
    f[[t]] <- append_jk_emulators(Design=D, f_list=f[[t]], current_beta = current_beta) 
    
    #### Get candidate points 
    C[[t]] <- get_candidate_points_pD(emulator=f[[t]][['all']], method='maxSE') 
    
    #### Get predictions Y_hat, SE at each candidate point
    Predicted_at_C[[t]] <- get_predictions_at_C(f_list=f[[t]], candidate_points=C[[t]])
    
    
    #### Get the ranges of the predictions Y_hat at each candidate point
    Prediction_range_at_C[[t]] <- get_prediction_range_at_C(Predicted_at_C[[t]])
    
    #### Obtain the candidate input with the maximum prediction range 
    C_star1_index[[t]] <- order(Prediction_range_at_C[[t]], decreasing = TRUE)[1]
    C_star1[[t]] <- C[[t]][C_star1_index[[t]], ]
    largest_prediction_range[[t]] <- Prediction_range_at_C[[t]][C_star1_index[[t]]]
    meets_T_resample[[t]] <- is_T_resample_met(cutoff=T_resample, 
                                               largest_prediction_range=largest_prediction_range, 
                                               n_consecutive=n_consecutive_T_resample)
    
    #### Obtain the candidate input with the maximum prediction SE 
    SE_predicted_at_C[[t]] <- sqrt(predict.GP(f[[t]][['all']], xnew=C[[t]])$MSE)
    C_star2_index[[t]] <- order(SE_predicted_at_C[[t]], decreasing = TRUE)[1]
    C_star2[[t]] <- C[[t]][C_star2_index[[t]], ]
    largest_prediction_SE[[t]] <- SE_predicted_at_C[[t]][C_star2_index[[t]]]
    meets_T_SE[[t]] <- is_T_SE_met(cutoff=T_SE, 
                                   largest_prediction_SE=largest_prediction_SE, 
                                   n_consecutive=n_consecutive_T_SE)
    
    #evaluate T_resample, T_SE stopping criteria for n_consecutive iterations
    if(!meets_T_resample[[t]]) {
      D <- add_design_point(Design=D, new_x=C_star1[[t]], iteration=t, simulator=simulator) 
      C_star[[t]] <- C_star1[[t]]
      C_star_index[[t]] <- C_star1_index[[t]]
    } else if(!meets_T_SE[[t]]) { 
      D <- add_design_point(Design=D, new_x=C_star2[[t]], iteration=t, simulator=simulator)
      C_star[[t]] <- C_star2[[t]]
      C_star_index[[t]] <- C_star2_index[[t]]
    } else{
      go_on <- FALSE
      C_star[[t]] <- C_star1[[t]]                # just save a value on exit, say from resample
      C_star_index[[t]] <- C_star1_index[[t]]    # just save a value on exit, say from resample
    }
    
    
    # overwrite current algorithm_trace lists
    algorithm_trace[['t']] <- t
    algorithm_trace[['D']] <- D
    algorithm_trace[['f']] <- f
    algorithm_trace[['C']] <- C
    algorithm_trace[['Predicted_at_C']] <- Predicted_at_C
    algorithm_trace[['Prediction_range_at_C']] <- Prediction_range_at_C
    algorithm_trace[['largest_prediction_range']] <- largest_prediction_range
    algorithm_trace[['largest_prediction_SE']] <- largest_prediction_SE
    algorithm_trace[['C_star1_index']] <- C_star1_index
    algorithm_trace[['C_star1']] <- C_star1 
    algorithm_trace[['C_star2_index']] <- C_star2_index
    algorithm_trace[['C_star2']] <- C_star2 
    algorithm_trace[['C_star_index']] <-C_star_index
    algorithm_trace[['C_star']] <- C_star
    
    
    if(!is.null(save_as_you_go_path)){
      if(!dir.exists(save_as_you_go_path)) { dir.create(save_as_you_go_path) }
      save_path <- file.path(save_as_you_go_path, sprintf('__interim_results_iteration_%03d.rds',t))
      saveRDS(algorithm_trace, file=save_path)
      cat(paste("\tInterim results saved in ", save_path, '\n', sep=''))
    }
  }
  
  return(algorithm_trace)
}
  
  

