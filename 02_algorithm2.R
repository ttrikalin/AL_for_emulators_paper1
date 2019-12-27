


## Temp: In 1 dimension -- just the midpoint is a pretty good approximation of the 
##       point that maximizes the SE of the prediction 
##       We use this instead of optimizing the SE of the predictions 

get_candidate_points_1D <- function(X) {
  return(get_centroids_1D(X))
}


## In p>=1 dimensions -- the Centroid of each simplex 
get_candidate_points_pD <- function(emulator, method='centroids', PSAPC_grid=NULL, constraints_list=NULL) {
  
    X <- emulator$X
    input_names <- colnames(X)
    candidate_points_df <- data.frame()
    
    if (method=='centroids') {
      candidate_points_df <- get_simplex_centroids_pD(X)
    } 
    
    if (method %in% c('maxSE', 'PSAPC_grid')) {
      optim_start_points <- get_simplex_centroids_pD(X)
      number_simplexes <- dim(optim_start_points)[1]
      for (s in 1:number_simplexes) {
          if(is.null(constraints_list)) {
              candidate_points_df[s, input_names] <- find_local_maxSE(emulator = emulator, 
                                                                    starting_point = optim_start_points[s,])  
          } else {
              candidate_points_df[s, input_names] <- find_local_maxSE_with_lin_constr(emulator = emulator, 
                                                                    starting_point = optim_start_points[s,], 
                                                                    constraints = constraints_list)  
          }
          
      }
    } 
    
    candidate_points_df <- unique(candidate_points_df)
    
    # Snap the candidate points to the closest point in the PSAPC_grid
    if (method=='PSAPC_grid') {
      if(is.null(PSAPC_grid)) {cat('\n, specify option PSAPC_grid')}
      for (point in 1:dim(candidate_points_df)[1]) {
          distance <- rep(0, dim(PSAPC_grid)[1])
          for(x in input_names) {
              distance <- distance + (candidate_points_df[point, x] - PSAPC_grid[, x])^2
          }
          candidate_points_df[point, input_names] <- PSAPC_grid[(distance == min(distance)), input_names ]
      
          ### remove the candidate points that are already in the emulator's design 
          X_included <- emulator$X
          for (existing_point in 1:dim(X_included)[1]) {
            if (identical(candidate_points_df[point, input_names], X_included[existing_point, input_names])) {
                candidate_points_df[point, input_names] <- rep(NA, length(input_names))
            }
          }
          candidate_points_df <- candidate_points_df[complete.cases(candidate_points_df),]
      }
      
      
      
    }
    
    return(candidate_points_df)
}









