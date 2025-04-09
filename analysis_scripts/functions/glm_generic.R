get_diagonal_from_event_start_and_end <- function(event_start, event_end, mat_length, max_beta){
  if(sum(event_start) > 0){
    
    # Get lengths of consecutive ones
    ones_lengths <- event_end - event_start
    ones_starts <-  event_start
    
    # ~~~ left aligned ~~~
    # create an empty matrix with number of rows equal to mat_length and columns equal to max_beta
    m_left <- matrix(0, nrow = mat_length, ncol = max_beta)
    
    # fill the first column of the matrix with 1s at the starting positions of sequences of ones
    m_left[event_start, 1] <- 1
    
    # loop through 2 to the maximum length of consecutive ones
    for(n_beta in 2:max(ones_lengths)) {
      # create a mask to filter sequences of ones that are at least n_beta long
      beta_mask <- ones_lengths >= n_beta
      # apply the mask to get the truncated start positions of the sequences of ones
      ones_starts_truc <- event_start[beta_mask]
      
      # fill the matrix at the appropriate positions to form diagonals
      m_left[ones_starts_truc + n_beta - 1, n_beta] <- 1
    }
    
    # ~~~ right aligned ~~~
    # create an empty matrix with number of rows equal to mat_length and columns equal to max_beta
    m_right <- matrix(0, nrow = mat_length, ncol = max_beta)
    
    # cill the last column of the matrix with 1s at the starting positions of sequences of ones
    m_right[ones_starts + ones_lengths - 1, max_beta] <- 1
    
    # loop through 2 to the maximum length of consecutive ones
    for(n_beta in 2:max(ones_lengths)) {
      # create a mask to filter sequences of ones that are at least n_beta long
      beta_mask <- ones_lengths >= n_beta
      # apply the mask to get the truncated start positions of the sequences of ones
      ones_stops_truc <- (ones_starts + ones_lengths - 1)[beta_mask]
      
      # fill the matrix at the appropriate positions to form diagonals
      m_right[ones_stops_truc - n_beta + 1, max_beta - n_beta + 1] <- 1
    }
    
    
  } else { # if there are no events, return a matrix of 0s
    m_left <- matrix(0, nrow = mat_length, ncol = max_beta)
    m_right <- matrix(0, nrow = mat_length, ncol = max_beta)
  }
  
  return(list(diag_matrix_left = m_left, diag_matrix_right = m_right))
}

get_diagonal_from_trial_ids <- function(df_trial_summary, var_id, n_samples_per_trial) {
  # Function to convert trial data into a diagonal matrix
  # This function takes a dataframe with trial summary (1 observation / trial), a variable ID, and the number of samples per trial,
  # and returns a matrix where each trial's values are represented as diagonal matrices, stacked row-wise.
  
  
  # Initialize an empty matrix to store the results
  matrix_session <- matrix(nrow = 0, ncol = n_samples_per_trial)
  
  # Loop through each unique trial ID
  for(value in df_trial_summary %>% pull(!!as.name(var_id))) {
    
    # Create a diagonal matrix with the extracted value and specified number of samples
    matrix_trial <- diag(1, n_samples_per_trial, n_samples_per_trial) * value
    
    # Append the diagonal matrix to the session matrix
    matrix_session <- rbind(matrix_session, matrix_trial)
  }
  
  # Return the resulting matrix
  return(matrix_session)
}


matrix_shift_horizontal <- function(matrix, previous_coord_start, previous_coord_end, new_coord_start, new_coord_end){
  
  matrix[new_coord_start[1]:new_coord_end[1], new_coord_start[2]:new_coord_end[2]]  <- 
    matrix[previous_coord_start[1]:previous_coord_end[1], previous_coord_start[2]:previous_coord_end[2]] 
  
  matrix[previous_coord_start[1]:previous_coord_end[1], previous_coord_start[2]:previous_coord_end[2]] <- 0
  
  return(matrix)
}
