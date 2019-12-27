##### Gramacy & Lee treed-GP stochastic emulator
library(tgp)


#### read needed functions 
# is_AL_stopping_criterion_met<- function(cutoff, 
#                                         largest_AL_statistic, 
#                                         n_consecutive=1, 
#                                         algorithm='Gramacy_Lee_BTGPLLM')
source('./00_algorithm1_functions.R')

#Read the algorithm
source('./03_Gramacy_Lee_tgp_algorithm.R')

### read the simulator 
source('benchmark_simulators/2D_simulator_BHP2012.R')



tmp_dir_stem <- 'tmp_2D_B'

for (seed_value in 5:9) {
  
  for (stop_criterion in c('ALC')) {
    if(!dir.exists(tmp_dir_stem)) {dir.create(tmp_dir_stem)}
    setwd(tmp_dir_stem)
    if(!dir.exists(as.character(seed_value))) {dir.create(as.character(seed_value))}
    setwd(as.character(seed_value))
    
    
    cat(paste('*********** RUNNING FOR SEED = ', seed_value, ' and ', stop_criterion, '\n',sep=''))  
    algorithm_trace <- get_sequential_design_GL_BTGPLL(simulator_info_list = simulator_info_list_BHP2012_2D,
                                                     n_consecutive_AL=5,
                                                     save_as_you_go_path=paste('./interim_results_2D_B', seed_value, sep=''),
                                                     continue_from_list =NULL, # list(algorithm_trace, TSE_ALM, TSE_ALC)
                                                     n_candidate_points = 200,
                                                     max_iter = 100,
                                                     AL_criterion = stop_criterion,      # ALM or ALC, 
                                                     seed = seed_value  # NULL to ignore
    )
    setwd('../..')
    saveRDS(algorithm_trace, file = paste('./benchmark_results/2D_BHP2012/2D_BHP2012_AL_BTGP_', stop_criterion, '_', seed_value, '.rds', sep=''))
  }
}

