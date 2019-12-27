##### Gramacy & Lee treed-GP stochastic emulator
library(tgp)


#### read needed functions 
# is_AL_stopping_criterion_met<- function(cutoff, 
#                                         largest_AL_statistic, 
#                                         n_consecutive=1, 
#                                         algorithm='Gramacy_Lee_BTGPLLM')
source('./00_algorithm1_functions.R')


# ### read the simulator 
# source('benchmark_simulators/2D_simulator_GL2008.R')
# 
# ### unpack the simulator input list 
# for(var_name in names(simulator_info_list_GL2008_2D)) {
#   assign(var_name, simulator_info_list_GL2008_2D[[var_name]])
# }
# T_SE_ALM <- T_SE
# T_SE_ALC <- 10^-5
# simulator_info_list_GL2008_2D[['T_SE_ALM']] <- T_SE_ALM
# simulator_info_list_GL2008_2D[['T_SE_ALC']] <- T_SE_ALC
# 
# 
# ## FOR DEBUGGING:
# simulator_info_list = simulator_info_list_GL2012_1D
# n_consecutive_AL=5
# save_as_you_go_path='./interim_results_2D_gl'
# continue_from_list =NULL # list(algorithm_trace, TSE_ALM, TSE_ALC)
# n_candidate_points = 200
# max_iter = 100
# AL_criterion = 'ALM'
# #


get_sequential_design_GL_BTGPLL <- function(simulator_info_list, 
                                     n_consecutive_AL=5, 
                                     save_as_you_go_path=NULL, 
                                     continue_from_list =NULL, # list(algorithm_trace, TSE_ALM, TSE_ALC)
                                     n_candidate_points = 200, 
                                     max_iter = 100, 
                                     AL_criterion = 'ALM',      # ALM or ALC
                                     seed=0) {
  
  
  
  if(is.null(continue_from_list)) {
      
    #setup lists 
      algorithm_trace <- list()     # algorithm_object
      for (var_name in names(simulator_info_list)) {
        assign(var_name, simulator_info_list[[var_name]])
        algorithm_trace[[var_name]] <- simulator_info_list[[var_name]]
      }
      algorithm_trace[['algorithm_name']] <- 'GL_BTGPLMM'  # Gramacy & Lee Bayesian Treed GP LLM models 
      algorithm_trace[['AL_criterion']] <- AL_criterion
      
      
      algorithm_trace[['max_iter']] <- max_iter
      algorithm_trace[['n_candidate_points']] <- n_candidate_points
      
      f <- list()          # The list of sequential Bayesian Tree GP LLM models per iteration 
      ALC <- list()        # The list of ALC values over C per iteration 
      ALM <- list()        # The list of ALC values over C per iteration 
      largest_ALC <- list()  # The largest statistic ALC over the input points (list over iteration)
      C_star_ALC <- list()        # The candidate point with the max ALC per iteration (list over iterations)
      C_star_ALC_index <- list()  # The index of the candidate point with the max ALC per iteration (list over iterations)
      largest_ALM <- list()       # The largest statistic ALM over the input points (list over iteration)
      C_star_ALM <- list()        # The candidate point with the max ALM per iteration (list over iterations)
      C_star_ALM_index <- list()  # The index of the candidate point with the max ALM per iteration (list over iterations)
      largest_AL_statistic <- list()  # The largest statistic (ALM or ALC, depending on which is acive) over the input points (list over iteration)
      C_star <- list()            # The selected next candidate point (list over iterations)
      C_star_index <- list()      # The index of the next selected candidate point (list over iterations)
  
      t<- 0
      go_on <- TRUE      
      
      
      
  } else {
    
    ###### START EDITS 
    
    # record newly-specified options (the var-names will be overwritten by the contents of algorithm_trace)
    new_n_consecutive_AL <- n_consecutive_AL
    new_AL_criterion <-  AL_criterion 
    new_max_iter <- max_iter 
    new_n_candidate_points <- n_candidate_points
    ### Over-write the algorithm_trace and various lists 
    cat('"continue_from_list" option was specified; examining whether stoping criteria are more stringent\n')
    algorithm_trace <- continue_from_list[['algorithm_trace']] 
    for (var_name in names(algorithm_trace)) {
      assign(var_name, algorithm_trace[[var_name]])
    }
    # restore option values 
    n_consecutive_AL <- new_n_consecutive_AL
    AL_criterion <-  new_AL_criterion
    max_iter <- new_max_iter 
    n_candidate_points <- new_n_candidate_points
    
    #use any newly provided T_SE_ALC or T_SE_ALM 
    tmp_T_SE_ALM <- continue_from_list[['T_SE_ALM']]
    if(!is.null(tmp_T_SE_ALM)) { T_SE_ALM <- tmp_T_SE_ALM }
    tmp_T_SE_ALC <- continue_from_list[['T_SE_ALC']]
    if(!is.null(tmp_T_SE_ALC)) { T_SE_ALC <- tmp_T_SE_ALC }
    
    # add list definitions that are not packaged in algorithm_trace
    
    # check that we have tightened some criterion: T_SE_[ALM | ALC],n_consecutive_AL
    tmp_is_ALC_met <- is_AL_stopping_criterion_met(cutoff = T_SE_ALC, 
                                                              largest_AL_statistic = largest_ALC, 
                                                              n_consecutive=n_consecutive_AL, 
                                                              algorithm='Gramacy_Lee_BTGPLLM')
    tmp_is_ALM_met <- is_AL_stopping_criterion_met(cutoff = T_SE_ALM, 
                                                   largest_AL_statistic = largest_ALM, 
                                                   n_consecutive=n_consecutive_AL, 
                                                   algorithm='Gramacy_Lee_BTGPLLM')
    cat(paste('is_T_ALC_more_stringent = ', !tmp_is_ALC_met, '\n',sep=''))
    cat(paste('is_T_ALM_more_stringent = ', !tmp_is_ALM_met, '\n',sep=''))

    cat(paste('Active criterion = ', AL_criterion,'\n',sep=''))
    
    if((AL_criterion=='ALC' & ! tmp_is_ALC_met ) | 
       (AL_criterion=='ALM' & ! tmp_is_ALM_met ) ) { 
      cat(paste('Convergence criteria for ', AL_criterion, ' are not met; continuing from iteration ', t, '\n', sep=''))  
      go_on <- TRUE 
    } else {
      cat('You need stricter convergence criteria\n')
      cat('Returning the algorithm trace from which I was asked to continue\n')
      return(algorithm_trace)
    }
    
    
    ### END EDITS 
  }
      
      
      if(AL_criterion=='ALM') {   
        T_SE_criterion <- T_SE_ALM 
      } else if (AL_criterion=='ALC'){
        T_SE_criterion <- T_SE_ALC 
      } else {
        cat('AL_criterion should be either ALC or ALM')
        return (NULL)
      }

  
      while (go_on == TRUE & t <= max_iter) {
        
        if(!is.null(seed)) { set.seed(seed)}
      
        t <- t +1 
        cat(paste('Now running iteration ', t, '\n', sep=''))  
      
        # Every 10 new design points, see if we need to update the D optimal design
        # trees split regions when you have at least 10 points in a region
        how_many_design_points <- dim(D)[1]
        if(t==1) {
          C<- get_candidate_points_GL_BTGPLLM(seeding_design = D, n_candidate_points=200, seed=seed) 
        } else if (how_many_design_points%%10 == 0) {
          C<- get_candidate_points_GL_BTGPLLM(seeding_design = D, n_candidate_points=200, seed=seed) 
        }
        
        
        
        ## Now get predictions at the candidate points C; 
        ##The easiest way is to predict at the same time we fit the GP LLM model 
        f[[t]] <- btgpllm(X= D[, input_names], Z=D$y, XX=C, pred.n =TRUE, corr='expsep', verb=0, 
                           Ds2x=TRUE) # Get the ALC criterion (points that most reduce the mean estimated variance)
        
        ALC[[t]] <- f[[t]][['Ds2x']]  # The Active-Learning Cohn criterion: reduction in squared error averaged over the input space 
        ALM[[t]] <- f[[t]][['ZZ.q']]  # The Active Learning McKay criterion: Vector of quantile norms (ranges)
        
        # find and record the point that is selected with the ALC criterion 
        largest_ALC[[t]] <- max(ALC[[t]])
        C_star_ALC_index[[t]] <- which(ALC[[t]] == largest_ALC[[t]])[1]
        C_star_ALC[[t]] <- matrix(C[C_star_ALC_index[[t]],], ncol=input_dimensions)
        
        # find and record the point that is selected with the ALM criterion
        largest_ALM[[t]] <- max(ALM[[t]])
        C_star_ALM_index[[t]] <- which(ALM[[t]] == largest_ALM[[t]])[1]
        C_star_ALM[[t]] <- matrix(C[C_star_ALM_index[[t]], ], ncol=input_dimensions)
        
        if(AL_criterion =='ALM') {                     # ALM 
          largest_AL_statistic[[t]] <- largest_ALM[[t]]
          C_star[[t]] = C_star_ALM[[t]] 
          C_star_index[[t]] = C_star_ALM_index[[t]] 
        } else {                                       # ALC
          largest_AL_statistic[[t]] <- largest_ALC[[t]]
          C_star[[t]] = C_star_ALC[[t]] 
          C_star_index[[t]] = C_star_ALC_index[[t]] 
        }
        
      
      
        is_stopping_criterion_met <- is_AL_stopping_criterion_met(cutoff = T_SE_criterion, 
                                                                  largest_AL_statistic = largest_AL_statistic, 
                                                                  n_consecutive=n_consecutive_AL, 
                                                                  algorithm='Gramacy_Lee_BTGPLLM')
        if(!is_stopping_criterion_met){
          D <- add_design_point(Design =D, new_x = C_star[[t]], iteration=t, simulator=simulator)
        } else {
          go_on <- FALSE
        }
        
      
  
  
  
      
      algorithm_trace[['t']] <- t
      algorithm_trace[['D']] <- D
      algorithm_trace[['f']] <- f
      algorithm_trace[['C']] <- C
      algorithm_trace[['ALC']] <- ALC
      algorithm_trace[['ALM']] <- ALM
      algorithm_trace[['largest_ALC']] <- largest_ALC
      algorithm_trace[['C_star_ALC']] <- C_star_ALC
      algorithm_trace[['C_star_ALC_index']] <- C_star_ALC_index
      algorithm_trace[['largest_ALM']] <- largest_ALM
      algorithm_trace[['C_star_ALM']] <-  C_star_ALM
      algorithm_trace[['C_star_ALM_index']] <- C_star_ALM_index
      algorithm_trace[['largest_AL_statistic']] <- largest_AL_statistic
      algorithm_trace[['C_star']] <- C_star
      algorithm_trace[['C_star_index']] <- C_star_index
      
      
      
      if(!is.null(save_as_you_go_path)){
        if(!dir.exists(save_as_you_go_path)) { dir.create(save_as_you_go_path) }
        save_path <- file.path(save_as_you_go_path, sprintf('__interim_results_iteration_%03d.rds',t))
        saveRDS(algorithm_trace, file=save_path)
        cat(paste("\tInterim results saved in ", save_path, '\n', sep=''))
      }
      
  }
      
  return(algorithm_trace)
      

}