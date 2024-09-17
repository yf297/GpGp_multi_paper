# parameter_initialization.R

#------------------------------------------
# Parameter Initialization Functions
#------------------------------------------

# Function to get parameter indices for the flexible-A model
get_parameter_indices <- function(ncomp) {
  neach <- ncomp * (ncomp + 1) / 2
  
  # Indices for variance parameters (sig)
  var_indices <- 1:neach
  
  # Indices for range parameters (ran)
  ran_indices <- (neach + 1):(2 * neach)
  
  # Index for Delta_B
  delta_B_index <- 2 * neach + 1
  
  # Indices for smoothness parameters (smo)
  smo_indices <- (2 * neach + 2):(3 * neach + 1)
  
  # Index for Delta_A
  delta_A_index <- 3 * neach + 2
  
  # Indices for nugget parameters (nug)
  nug_indices <- (3 * neach + 3):(4 * neach + 2)
  
  # Create index mapping for parameters
  index_mat <- matrix(NA, ncomp, ncomp)
  index_vec <- 1:neach
  index_mat[upper.tri(index_mat, diag = TRUE)] <- index_vec
  index_mat <- t(index_mat)
  
  # Diagonal indices (marginal parameters)
  diag_indices <- which(row(index_mat) == col(index_mat), arr.ind = TRUE)
  marginal_var_indices <- var_indices[index_mat[diag_indices]]
  marginal_ran_indices <- ran_indices[index_mat[diag_indices]]
  marginal_smo_indices <- smo_indices[index_mat[diag_indices]]
  marginal_nug_indices <- nug_indices[index_mat[diag_indices]]
  
  # Off-diagonal indices (cross parameters)
  cross_indices <- which(row(index_mat) != col(index_mat) & !is.na(index_mat), arr.ind = TRUE)
  cross_var_indices <- var_indices[index_mat[cross_indices]]
  cross_ran_indices <- ran_indices[index_mat[cross_indices]]
  cross_smo_indices <- smo_indices[index_mat[cross_indices]]
  cross_nug_indices <- nug_indices[index_mat[cross_indices]]
  
  return(list(
    var_indices = var_indices,
    ran_indices = ran_indices,
    delta_B_index = delta_B_index,
    smo_indices = smo_indices,
    delta_A_index = delta_A_index,
    nug_indices = nug_indices,
    cross_var_indices = cross_var_indices,
    cross_ran_indices = cross_ran_indices,
    cross_smo_indices = cross_smo_indices,
    cross_nug_indices = cross_nug_indices,
    marginal_var_indices = marginal_var_indices,
    marginal_ran_indices = marginal_ran_indices,
    marginal_smo_indices = marginal_smo_indices,
    marginal_nug_indices = marginal_nug_indices
  ))
}


# Function to initialize starting parameters for the flexible-A model
initialize_parameters <- function(fit1, param_indices, ncomp) {
  # Extract marginal parameter estimates from fit1
  marginal_ran_estimates <- fit1$logparms[param_indices$marginal_ran_indices]
  marginal_smo_estimates <- fit1$logparms[param_indices$marginal_smo_indices]
  marginal_var_estimates <- fit1$logparms[param_indices$marginal_var_indices]
  marginal_nug_estimates <- fit1$logparms[param_indices$marginal_nug_indices]
  
  # Start with marginal parameters from the independent fit
  start_logparms <- fit1$logparms
  
  # Initialize Deltas
  start_logparms[param_indices$delta_A_index] <- log(0.001)  # Small positive value
  start_logparms[param_indices$delta_B_index] <- log(0.01)   # Small positive value
  
  M <- matrix(0.5, ncomp, ncomp)
  diag(M) <- 1
  start_logparms[param_indices$cross_ran_inds]  <- log(finv(M)) 
  start_logparms[param_indices$cross_smo_inds]  <-  log(finv(M))  

  
  return(start_logparms)
}
