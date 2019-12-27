library(GPfit)
library(geometry)
library(matrixStats)
library(stringr)
library(MaxPro)  
options(stringsAsFactors = FALSE)   #very important to avoid promoting strings to factors anywhere in this set of functions. 



#########################################
### Seeding design helper functions   
#########################################


####  Get a seeding design for a simulator using Max Projections algorithm with Euclidian norm
####  This gives us the internal points (but also some external points)
get_seeding_design_MaxPro <- function(number_of_interior_points, input_dimensions) {
  if( input_dimensions==1 ) {
    D <- matrix(data=seq(from=0, to=1, length.out = number_of_interior_points+2), ncol=1)
  } else {
    D <- MaxProLHD(n=number_of_interior_points, p=input_dimensions)$Design
    D <- MaxPro(InitialDesign = D)$Design  
  }
  D<- as.data.frame(D)
  colnames(D) <- paste('x', seq(1:dim(D)[2]), sep='')
  point  <- paste('i', seq(1:dim(D)[1]), sep='')
  D<- cbind(point , D)
  return(D)
}

### Get the design hull 
get_design_hull <- function(input_dimensions) {
  expand_list <- list()
  for(p in 1:input_dimensions) {
    expand_list[[paste('x', p, sep='')]] <- 0:1
  }
  H <- expand.grid(expand_list)
  point <- paste('h', seq(1, dim(H)[1]), sep='')
  return(cbind(point, H))
}

# Combine the design hull with the interior points, and drop any redundancies from the interior points 
get_seeding_design <-  function(number_of_interior_points, input_dimensions) {
    D_hull <- get_design_hull(input_dimensions = input_dimensions)
    D_interior <- get_seeding_design_MaxPro(number_of_interior_points = number_of_interior_points, input_dimensions = input_dimensions)
    
    x_indices <- 2:dim(D_interior)[2]
    # define the radius of a small ball: say, 1/10 of the Euclidean distance between the MaxPro points 1, 2
    # If you use a different distance in the MaxPro, fix here as well 
    s <- 2 # Euclidean distance 
    small_distance <- (sum((D_interior[1,x_indices] - D_interior[2,x_indices])^s))^(1/s) / 10
    
    
    n_hull <- dim(D_hull)[1]
    n_interior <- dim(D_interior)[1]
    to_keep <- rep(TRUE, n_interior)
    for (i in 1:n_interior) {
      examine_next_hull_point <- TRUE 
      h <- 0 
      while(examine_next_hull_point & h<n_hull) {
        h<- h+1
        tmp_distance <- (sum((D_interior[i,x_indices] - D_hull[h,x_indices])^s))^(1/s)
        if (tmp_distance<=small_distance) {
          examine_next_hull_point <- FALSE 
          to_keep[i] <- FALSE
        }
      }
    }

    D <- rbind(D_hull, D_interior[to_keep,])
    return(D)
}



#########################################
### Algorithm 1 helper functions   
#########################################



#### Append jack-knifed emulators by dropping interior points 
#### f[[t]][['all']]: f will all points used. 
#### f[[t]][['i1']]: f_(-i1) and so on. 
append_jk_emulators <- function(Design, f_list, current_beta=NULL) {
  interior_points <- Design$point[which(substr(Design$point, 1, 1)=='i')]
  x_indices <- grep('x', names(Design))
  for (point in interior_points) {
    Design_minus_i <- Design[which(Design$point!=point),]
    
    f_list[[point]] <- GP_fit(X=Design_minus_i[,x_indices], Y=Design_minus_i$y, optim_start = current_beta)
  }
  return(f_list)
}



#### Get predictions Y_hat at each candidate point
get_predictions_at_C <- function(f_list, candidate_points){
  # if the 1D case breaks check this recasting of vector into a matrix:
  input_dimensions <- dim(f_list[[1]]$X)[2]
  if(input_dimensions>1 & is.vector(candidate_points)) { 
    candidate_points <- matrix(candidate_points, nrow=1)  
  }
  emulator_names_list <- names(f_list)
  i <- 0
  for(name in emulator_names_list){
    i<- i+1
    tmp <- NULL
    tmp <- predict.GP(f_list[[name]], xnew = candidate_points)$Y_hat
    if (i==1) { Y_hat <- tmp }
    else      { Y_hat <- cbind(Y_hat, tmp) }
  }
  colnames(Y_hat) <- emulator_names_list
  rownames(Y_hat) <- paste('c', seq(1:length(tmp)), sep='')
  return(Y_hat)
}



#### Get predictions Y_hat at each candidate point
get_prediction_SE_at_C <- function(emulator, candidate_points) {
  input_dimensions <- dim(emulator$X)[2]
  if(input_dimensions>1 & is.vector(candidate_points)) { 
    candidate_points <- matrix(candidate_points, nrow=1)  
  }
  SE <- sqrt(predict.GP(emulator, xnew = candidate_points)$MSE)
  return(SE)
}




#### Get the ranges of the predictions Y_hat at each candidate point
get_prediction_range_at_C <- function(predictions) {
  ranges <- rowMaxs(predictions) - rowMins(predictions)
  names(ranges) <- rownames(predictions)
  return(ranges)
}





#### Add the design point 
add_design_point <- function(Design, new_x, iteration, simulator) {
  cat('\t* Adding a new design point\n')
  current_point_names <- grep('i', Design$point, value=TRUE)
  
  input_names <- grep('x',colnames(Design), value = TRUE)
  if(length(input_names)>1) { new_x <- as.matrix(new_x, ncol=length(input_names)) }
  
  new_point_number <- 1 + max(as.numeric(str_replace(current_point_names, 'i', '')))  
  new_point_name <- paste('i', new_point_number, sep='')
  
  r <- nrow(Design)+1
  Design[r,'point'] <- new_point_name
  Design[r,'t'] <- iteration
  Design[r,input_names] <- new_x
  Design[r,'y'] <- simulator(new_x)
  
  rownames(Design) <- Design$point
  return(Design)
}




#### Examine if the resampling threshold has been met 
is_T_resample_met <- function(cutoff, largest_prediction_range, n_consecutive=1, K=NULL) {
     if(is.null(K)) {K <- length(largest_prediction_range)}
     if (K<n_consecutive) { 
         return(FALSE) 
     } else {
          cutoff_met <- TRUE 
          for (t in (K-n_consecutive+1):K) {
              cutoff_met <- ((largest_prediction_range[[t]] <= cutoff) & cutoff_met)
          }
          names(cutoff_met) <- 'is_T_resample_met'
          return(cutoff_met)  
     }
}

#### Examine if the T_SE threshold has been met 
is_T_SE_met <- function(cutoff, largest_prediction_SE, n_consecutive=1, K=NULL) {
  if(is.null(K)) {K <- length(largest_prediction_SE)}
  if (K<n_consecutive) { 
    return(FALSE) 
  } else {
    cutoff_met <- TRUE 
    for (t in (K-n_consecutive+1):K) {
      cutoff_met <- ((largest_prediction_SE[[t]] <= cutoff) & cutoff_met)
    }
    names(cutoff_met) <- 'is_T_SE_met'
    return(cutoff_met)  
  }
}




#########################################
### Algorithm 2 helper functions   
#########################################

get_simplex_centroid <- function(simplex) {
  return(colMeans(simplex))
}

get_centroids_from_tesselation <- function(triangulation_full_object){
    X <- triangulation_full_object$p
    tri <- triangulation_full_object$tri
    
    number_simplexes <- dim(tri)[1]
    input_dimensions <- dim(X)[2]
    
    X_centroids <- matrix(rep(NA, number_simplexes*input_dimensions), ncol=input_dimensions)
    for (s in 1:number_simplexes) {
      X_centroids[s,] <- get_simplex_centroid(X[tri[s,], ] )
    }
    colnames(X_centroids) <- colnames(X)
    return(X_centroids)
}





get_simplex_centroids_1D <-function(X) {
  indexes_hull <- delaunayn(X)    
  return ((X[indexes_hull[,1]] + X[indexes_hull[,2]]) / 2)
}


get_simplex_centroids_pD <- function(X) {
  T_obj <- delaunayn(X, output.options=c('Fn','Fa'))    
  X_centroids <- get_centroids_from_tesselation(triangulation_full_object = T_obj)
  
  # order them in decreasing measure for each simplex 
  decreasing_area_order <- order(T_obj$areas, decreasing = TRUE)
  X_centroids <- X_centroids[decreasing_area_order, ]
  
  return (X_centroids)
}


## get points with max prediction error
optim_prediction_SE_at_par <- function(par, emulator){
    SE <- get_prediction_SE_at_C(emulator = emulator, candidate_points = par)
    return(SE)
}

find_local_maxSE <- function(emulator, starting_point) {
  
    input_dimensions <- dim(emulator$X)[2]
    optimization_result <- optim(par=starting_point, 
                                 fn=optim_prediction_SE_at_par, 
                                 method = 'L-BFGS-B', 
                                 lower=rep(0, input_dimensions), 
                                 upper=rep(1, input_dimensions),
                                 emulator = emulator, 
                                 control=list('fnscale'=-1,'factr'=10^10))  

  return(optimization_result$par)
}


find_local_maxSE_with_lin_constr <- function(emulator, starting_point, constraints_list) {
  
  input_dimensions <- dim(emulator$X)[2]
  optimization_result <- constrOptim(theta=starting_point, 
                               f=optim_prediction_SE_at_par, 
                               grad=NULL,
                               ui=constraints_list[['constraints_matrix']],
                               ci=constraints_list[['constraints_rhs']], 
                               emulator = emulator, 
                               control=list('fnscale'=-1)  )
  
  return(optimization_result$par)
}


#########################################
### Helper functions for Benchmark algorithms
#########################################




#### Examine if the stopping threshold has been met 
is_AL_stopping_criterion_met<- function(cutoff, 
                                        largest_AL_statistic, 
                                        n_consecutive=1, 
                                        algorithm='Gramacy_Lee_BTGPLLM', 
                                        K=NULL) {
  if (algorithm == 'Gramacy_Lee_BTGPLLM') {
    is_met <- is_AL_stopping_criterion_met_GLBTGPLMM(cutoff=cutoff, largest_AL_statistic=largest_AL_statistic, n_consecutive=n_consecutive, K=K)
  } else {
    cat(paste('algorithm ', algorithm, ' is not implemented. Specify a valid option\n', sep=''))
  }
  return(is_met)

}


is_AL_stopping_criterion_met_GLBTGPLMM <- function(cutoff, largest_AL_statistic, n_consecutive=1, K=NULL) {
  if(is.null(K)) {K <- length(largest_AL_statistic)}
  if (K<n_consecutive) { 
    return(FALSE) 
  } else {
    cutoff_met <- TRUE 
    for (t in (K-n_consecutive+1):K) {
      cutoff_met <- ((largest_AL_statistic[[t]] <= cutoff) & cutoff_met)
    }
    names(cutoff_met) <- 'is_AL_stopping_criterion_met'
    return(cutoff_met)  
  }
}

### Get D-optimal candidate points For the Gramacy & Lee BTGPLLM algorithm 
get_candidate_points_GL_BTGPLLM <- function(seeding_design, n_candidate_points=200, seed=0) {
  ## Gramacy and Lee 
  ## create a large number of random design points 
  ## and then select the n_candidate_points that are D-optimal 

  input_names <- grep('x', colnames(seeding_design), value = TRUE)
  input_dimensions <- length(input_names)
  
  x_mat <- c(0,1)
  if(input_dimensions>1) {
    for (j in 2:input_dimensions) {
      x_mat <- rbind(x_mat, c(0, 1))
    }
  }
  
  if(!is.null(seed)) {set.seed(seed)}
  Xcand <- lhs(5000, x_mat)
 
  
  ## learn an initial tree partitioning by using a GP LLM
  tmp_mdl <- btgpllm(X= seeding_design[, input_names], Z=seeding_design$y, pred.n =FALSE, corr='expsep', verb=0)
  
  # Use tgp.design() to get a set of D-optimal candidate design points 
  # given the current partitioning; the sequential algorithm will samle from these
  # These are the candidate points
  C <- tgp.design(howmany=n_candidate_points, Xcand = Xcand, tmp_mdl)  # XX in Gramacy's notation
  colnames(C) <- input_names
  rownames(C) <- paste('c', 1:(dim(C)[1]), sep='')
  return(C)
}



### Get D-optimal candidate points For the Gramacy & Lee BTGPLLM algorithm 
get_candidate_points_GL_BTGPLLM_for_PSAPC <- function(seeding_design, PSAPC_grid, n_candidate_points=200, seed=0) {
  ## Gramacy and Lee 
  
  input_names <- grep('x', colnames(seeding_design), value = TRUE)
  input_dimensions <- length(input_names)
  
  x_mat <- c(0,1)
  if(input_dimensions>1) {
    for (j in 2:input_dimensions) {
      x_mat <- rbind(x_mat, c(0, 1))
    }
  }
  
  if(!is.null(seed)) {set.seed(seed)}
  Xcand <-  PSAPC_grid[, input_names]  
  
  
  ## learn an initial tree partitioning by using a GP LLM
  tmp_mdl <- btgpllm(X= seeding_design[, input_names], Z=seeding_design$y, pred.n =FALSE, corr='expsep', verb=0)
  
  # Use tgp.design() to get a set of D-optimal candidate design points 
  # given the current partitioning; the sequential algorithm will samle from these
  # These are the candidate points
  C <- tgp.design(howmany=n_candidate_points, Xcand = Xcand, tmp_mdl)  # XX in Gramacy's notation
  colnames(C) <- input_names
  rownames(C) <- paste('c', 1:(dim(C)[1]), sep='')
  return(C)
}



#### Add the design point 
add_design_point_for_PSAPC <- function(Design, new_x, iteration, simulator, PSAPC_grid) {
  cat('\t* Adding a new design point\n')
  current_point_names <- grep('i', Design$point, value=TRUE)
  
  input_names <- grep('x',colnames(Design), value = TRUE)
  if(length(input_names)>1) { new_x <- as.matrix(new_x, ncol=length(input_names)) }
  
  new_point_number <- 1 + max(as.numeric(str_replace(current_point_names, 'i', '')))  
  new_point_name <- paste('i', new_point_number, sep='')
  
  r <- nrow(Design)+1
  Design[r,'point'] <- new_point_name
  Design[r,'t'] <- iteration
  Design[r,input_names] <- new_x
  Design[r,'y'] <- simulator(new_x, PSAPC_grid = PSAPC_grid)
  
  rownames(Design) <- Design$point
  return(Design)
}

