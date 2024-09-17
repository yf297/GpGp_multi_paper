# Load required packages and source necessary scripts
library("GpGpm")

# Source helper functions
source("R/helper.R")
source("R/fisher_multi.R")
source("R/fit_multi.R")
source("R/link_multi.R")
source("R/check.R")

#------------------------------------------
# Load and Format Data
#------------------------------------------

# Load the weather data
data <- get(load("data/formatted/Weather.RData"))

# Extract response vector y and locations locs
y <- as.vector(data[, 1])
locs <- as.matrix(data[, -1])

# Convert the component labels in locs to numeric factor levels
component_labels <- as.numeric(as.factor(locs[, ncol(locs)]))
locs[, ncol(locs)] <- component_labels

# Create design matrix X (excluding intercept)
X <- model.matrix(~ -1 + factor(component_labels))

#------------------------------------------
# Define Order and Neighbor Functions
#------------------------------------------

# Define order and neighbor functions
order_fun <- get("order_completely_random")
neighbor_fun <- get("nearest_multi_any")

# Set number of neighbors for Vecchia approximation
m <- 40

# Indicator for whether to include correlated nugget terms
corrnug <- FALSE

#------------------------------------------
# Define Number of Components
#------------------------------------------

# Number of components
ncomp <- length(unique(locs[, ncol(locs)]))

#------------------------------------------
# First Fit: Independent Components
#------------------------------------------

# Obtain starting parameters
start_parms <- get_start_parms(y, X, locs, "matern_multi")$start_parms

# Take logarithm of starting parameters
start_logparms <- log(start_parms)

# Append zeros for Delta_B and Delta_A at appropriate positions
neach <- ncomp * (ncomp + 1) / 2
start_logparms <- append(start_logparms, 0, after = 2 * neach)
start_logparms <- append(start_logparms, 0, after = 3 * neach + 1)

# Get parameter indices
param_indices <- get_parameter_indices(ncomp)

# Define link functions
linkfuns <- list(
  link = link_flex_a,
  dlink = dlink_flex_a
)

# Set cross variance and cross nugget parameters to zero in start_logparms
start_logparms[param_indices$cross_var_indices] <- 0
start_logparms[param_indices$cross_nug_indices] <- 0

# Set active parameters (those to be estimated)
active_logparms <- rep(TRUE, length(start_logparms))

# Deactivate cross parameters (set to FALSE)
active_logparms[c(
  param_indices$cross_var_indices,
  param_indices$cross_ran_indices,
  param_indices$delta_B_index,
  param_indices$cross_smo_indices,
  param_indices$delta_A_index,
  param_indices$cross_nug_indices
)] <- FALSE

cat("Starting independent fit...\n")

# Fit the model with independent components
fit1 <- fit_multi(
  y = y,
  locs = locs,
  X = X,
  neighbor_fun = neighbor_fun,
  order_fun = order_fun,
  m = m,
  start_logparms = start_logparms,
  linkfuns = linkfuns,
  active_logparms = active_logparms
)

#------------------------------------------
# Second Fit: Flexible-A Model
#------------------------------------------

# Initialize starting parameters for flexible-A model
start_logparms <- initialize_parameters(fit1, param_indices, ncomp)

# Deactivate cross nugget parameters if correlated nugget is not used
active_logparms <- rep(TRUE, length(start_logparms))
if (!corrnug) {
  active_logparms[param_indices$cross_nug_indices] <- FALSE
}

cat("Starting flexible-A model fit...\n")

# Fit the model with flexible-A cross-covariance structure
fit2 <- fit_multi(
  y = y,
  locs = locs,
  X = X,
  neighbor_fun = neighbor_fun,
  order_fun = order_fun,
  m = m,
  start_logparms = start_logparms,
  linkfuns = linkfuns,
  active_logparms = active_logparms
)



