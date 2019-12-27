library(GPfit)
library(tgp)


### function to get the metrics in 1D
get_metrics_1D <- function(simulator, emulator, step=0.005) {
  x_list <- seq(0, 1, step)
  y_hat <- predict.GP(emulator, xnew=x_list)$Y_hat
  y_simulator <- simulator(x_list)
  SIMULATOR_RANGE <- max(y_simulator) - min(y_simulator)
  ABS_DIFF <- abs(y_simulator - y_hat)
  RMSE <- sqrt(mean(as.matrix(ABS_DIFF)^2))
  MAX <- max(ABS_DIFF)  
  
  x_MAX <- x_list[order(ABS_DIFF, decreasing = TRUE)[1]]
  return(list('MAX'=MAX, 'x_MAX'=x_MAX, 'RMSE'=RMSE,'sMAX'=MAX/SIMULATOR_RANGE, 'sRMSE'=RMSE/SIMULATOR_RANGE))
}


get_input_grid <- function(input_names, step = NULL) {
  input_dimensions <- length(input_names)
  
  if (is.null(step)) {
    number_of_points <- min(10000, 100^input_dimensions)
    grid_point_list <- seq(from=0, to=1, length.out=number_of_points^(1/input_dimensions))
  } else {
    grid_point_list <- seq(from=0, to=1, by=step)
  }
  expand_list <- list()
  for(p in 1:input_dimensions) {
    expand_list[[paste('x', p, sep='')]] <- grid_point_list
  }
  x_list <- expand.grid(expand_list)
  return(x_list)
}



### function to get the metrics in p>=1 dimensions
get_metrics_over_emulator_list_pD <- function(simulator, emulator_list, step=NULL) {
  
  number_of_emulators <- length(emulator_list) 
  input_dimensions <- dim(emulator_list[[1]]$X)[2]
  
  #input_names <- grep('x',colnames(emulator_list[[1]]$X), value = TRUE)
  if (input_dimensions==1) {
    input_names <- 'x'
  } else {
      input_names <- paste('x', 1:input_dimensions, sep = '')
  }
  
  x_list <- get_input_grid(input_names = input_names, step = step)
  
  y_simulator <- simulator(x_list)
  SIMULATOR_RANGE <- max(y_simulator) - min(y_simulator)
  
  results <- data.frame()
  cat(paste('\nExamining ', number_of_emulators, ' emulators:\n', sep=''))
  for(f_i in 1:number_of_emulators) {
    cat(paste(f_i, '.', sep=''))
    emulator <- emulator_list[[f_i]]
    if(class(emulator)=='GP') {
      y_hat <- predict.GP(emulator, xnew=x_list)$Y_hat
      results[f_i, 'n_points'] <- length(emulator$Y)
    } else if(class(emulator)=='tgp') {
      y_hat <- tryCatch( { predict(object=emulator, XX = x_list, pred.n = FALSE, MAP=TRUE)$ZZ.km },
                        error=function(cond) { return(rep(NA, dim(x_list)[1]))}
                       )
      results[f_i, 'partitions'] <- length(emulator$trees)
      results[f_i, 'n_points'] <- emulator$n
    }
      
    emulator_name <- names(emulator_list)[[f_i]]
    results[f_i, 'iteration'] <- f_i
    if(is.null(emulator_name)) {
    
      results[f_i, 'emulator'] <- f_i
    } else { 
      results[f_i, 'emulator'] <- emulator_name
    }
    results[f_i, 'simulator_range'] <- SIMULATOR_RANGE
    
    ABS_DIFF <- as.matrix(abs(y_simulator - y_hat))
    x_MAX <- x_list[order(ABS_DIFF, decreasing = TRUE)[1],]
    results[f_i, 'MAX'] <- max(ABS_DIFF)
    results[f_i, 'sMAX'] <- results[f_i, 'MAX'] / SIMULATOR_RANGE
    results[f_i, paste(input_names, '_MAX', sep='')] <- x_list[order(ABS_DIFF, decreasing = TRUE)[1],]
    results[f_i, 'RMSE'] <- sqrt(mean(ABS_DIFF^2))
    results[f_i, 'sRMSE'] <- results[f_i, 'RMSE'] / SIMULATOR_RANGE
    
  }
  return(results)
}


### 
### The above function to accommodate the PSAPC emulator
get_metrics_over_emulator_list_pD_for_PSAPC <- function(simulator, PSAPC_grid, emulator_list, step=NULL) {
  
  
  number_of_emulators <- length(emulator_list) 
  input_dimensions <- dim(emulator_list[[1]]$X)[2]
  
  #input_names <- grep('x',colnames(emulator_list[[1]]$X), value = TRUE)
  if (input_dimensions==1) {
    input_names <- 'x'
  } else {
    input_names <- paste('x', 1:input_dimensions, sep = '')
  }
  
  #x_list <- get_input_grid(input_names = input_names, step = step)
  x_list <- PSAPC_grid[, input_names]
  y_simulator <- simulator(x_list, PSAPC_grid=PSAPC_grid)
  SIMULATOR_RANGE <- max(y_simulator) - min(y_simulator)
  
  
  results <- data.frame()
  cat(paste('Running over ', number_of_emulators, ' iterations:\n', sep=''))
  for(f_i in 1:number_of_emulators) {
    cat(paste(f_i, '.', sep=''))
    emulator <- emulator_list[[f_i]]
    if(class(emulator)=='GP') {
      y_hat <- predict.GP(emulator, xnew=x_list)$Y_hat
      results[f_i, 'n_points'] <- length(emulator$Y)
    } else if(class(emulator)=='tgp') {
      y_hat <- predict(object=emulator, XX = x_list, pred.n = FALSE, MAP=TRUE)$ZZ.km
      results[f_i, 'partitions'] <- length(emulator$trees)
      results[f_i, 'n_points'] <- emulator$n
    }
    
    emulator_name <- names(emulator_list)[[f_i]]
    results[f_i, 'iteration'] <- f_i
    if(is.null(emulator_name)) {
      
      results[f_i, 'emulator'] <- f_i
    } else { 
      results[f_i, 'emulator'] <- emulator_name
    }
    results[f_i, 'simulator_range'] <- SIMULATOR_RANGE
    
    ABS_DIFF <- as.matrix(abs(y_simulator - y_hat))
    x_MAX <- x_list[order(ABS_DIFF, decreasing = TRUE)[1],]
    results[f_i, 'MAX'] <- max(ABS_DIFF)
    results[f_i, 'sMAX'] <- results[f_i, 'MAX'] / SIMULATOR_RANGE
    results[f_i, paste(input_names, '_MAX', sep='')] <- x_MAX
    results[f_i, 'RMSE'] <- sqrt(mean(ABS_DIFF^2))
    results[f_i, 'sRMSE'] <- results[f_i, 'RMSE'] / SIMULATOR_RANGE
    
  }
  return(results)
}








### function to get predicted values and CIs in 1D
get_meanY_and_95CI_1D <- function(emulator, step=0.005) {
    x_list <- seq(0, 1, step)
    prediction_list <- predict.GP(emulator, xnew=x_list)
    
    prediction <- cbind(prediction_list$Y_hat, 
                        prediction_list$Y_hat - 1.96*sqrt(prediction_list$MSE), 
                        prediction_list$Y_hat + 1.96*sqrt(prediction_list$MSE))
    colnames(prediction) <- c('mean', 'lower', 'upper')
    return(prediction)
}



### function to get predicted values and CIs in 2D
get_meanY_and_95CI_2D <- function(emulator, step=0.01) {
  x_list <- expand.grid(list('x1'=seq(0, 1, step), 'x2'=seq(0, 1, step)))
  prediction_list <- predict.GP(emulator, xnew=x_list)
  
  prediction <- cbind(prediction_list$Y_hat, 
                      prediction_list$Y_hat - 1.96*sqrt(prediction_list$MSE), 
                      prediction_list$Y_hat + 1.96*sqrt(prediction_list$MSE))
  colnames(prediction) <- c('mean', 'lower', 'upper')

  prediction <- cbind(x_list, prediction)
  return(prediction)
}





str_eval <- function(x) {
  return(eval(parse(text=x)))
  }



