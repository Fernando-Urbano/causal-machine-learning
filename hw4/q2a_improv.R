###########################################################
# Simulation Code for High-Dimensional Parametric Case
#
# The code below outlines a modular simulation study in R.
# The study examines two estimators (plugin and doubly robust)
# for the ATT under various data-generating processes (DGPs),
# using penalized regression methods (LASSO and Ridge).
#
# The code:
# 1. Defines functions for:
#    - Data generation
#    - Fitting models (mu0_hat and p_hat)
#    - Computing plug-in and DR estimators
#    - Running the simulation across combinations of parameters
#    - Saving results to a LaTeX table and CSV
#
# 2. Loops over sample sizes, dimensions, sparsity levels, and penalty types.
# 3. Stores individual replication results in a data.frame, converts it to a LaTeX table, and saves it.
#
# Note: 
# - You do not need to run this code. Just place it in an R script and run.
# - Make sure necessary packages are installed: glmnet, xtable, dplyr
# - The code as given is a template and may need minor adjustments.
#
###########################################################

# Load necessary packages
library(glmnet)
library(xtable)
library(dplyr)  # Added for write_csv function

###########################################################
# Function to generate true parameter vectors with given sparsity
###########################################################
generate_parameters <- function(d, s, type = c("beta","theta"), seed = NULL) {
  # type: "beta" or "theta" just for clarity, they are handled similarly.
  # s = number of non-zero entries.
  if (!is.null(seed)) set.seed(seed)
  
  # Nonzero elements are drawn from say N(0,1), zero otherwise
  param <- numeric(d)
  nonzero_indices <- sample(1:d, size = s)
  param[nonzero_indices] <- rnorm(s, mean = 1, sd = 1)  # can adjust if needed
  return(param)
}

###########################################################
# Function to generate data
###########################################################
generate_data <- function(n, d, beta0, theta0, tau = 1, seed = NULL) {
  # Data generating process:
  # X ~ N(0, I_d)
  # p(x) = 1 / (1 + exp(-theta0' x))
  # mu_0(x) = beta0' x
  # T ~ Bernoulli(p(x))
  # Y = T * (mu_0(X) + tau) + (1 - T)*mu0(X) + noise
  # Here we add some noise to Y. Let's say noise ~ N(0,1)
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- matrix(rnorm(n*d), nrow = n, ncol = d)
  p_x <- 1 / (1 + exp(- X %*% theta0))
  T_ind <- rbinom(n, 1, p_x)
  
  mu0_x <- X %*% beta0
  Y <- T_ind * (mu0_x + tau) + (1 - T_ind)*mu0_x + rnorm(n,0,1)
  
  return(list(X = X, Y = Y, T = T_ind, p_true = p_x, mu0_true = mu0_x))
}

###########################################################
# Function to estimate p(x) using penalized logistic regression
# penalty_type = "lasso" or "ridge"
###########################################################
estimate_p <- function(X, T, penalty_type = "lasso") {
  # For logistic regression: family = "binomial"
  # LASSO: alpha=1, Ridge: alpha=0
  alpha_val <- ifelse(penalty_type == "lasso", 1, 0)
  fit_p <- glmnet(X, T, family = "binomial", alpha = alpha_val)
  
  # Use cross-validation to choose lambda
  cvfit_p <- cv.glmnet(X, T, family = "binomial", alpha = alpha_val)
  lambda_best <- cvfit_p$lambda.min
  
  # final model
  p_hat_func <- function(newX) {
    preds <- predict(cvfit_p, newX, s = "lambda.min", type = "response")
    return(preds)
  }
  
  return(p_hat_func)
}

###########################################################
# Function to estimate mu0(x) using penalized linear regression
# penalty_type = "lasso" or "ridge"
###########################################################
estimate_mu0 <- function(X, Y, T, penalty_type = "lasso") {
  # We only use the control units (T=0) to estimate mu0
  X0 <- X[T==0, , drop=FALSE]
  Y0 <- Y[T==0]
  
  alpha_val <- ifelse(penalty_type == "lasso", 1, 0)
  fit_mu0 <- glmnet(X0, Y0, family = "gaussian", alpha = alpha_val)
  
  # Use cross-validation to choose lambda
  cvfit_mu0 <- cv.glmnet(X0, Y0, family = "gaussian", alpha = alpha_val)
  
  mu0_hat_func <- function(newX) {
    preds <- predict(cvfit_mu0, newX, s = "lambda.min", type = "response")
    return(preds)
  }
  
  return(mu0_hat_func)
}

###########################################################
# Function to compute the plug-in estimator:
#   hat{tau}_PI = hat{mu}_1 - hat{mu}_0
# where:
#   hat{mu}_1 = (1/n) * sum_{i=1}^n [T_i Y_i / hat{p}]
#   hat{mu}_0 = (1/n) * sum_{i=1}^n [(1 - T_i) * hat{p}(X_i)*Y_i / (1 - hat{p}(X_i))]
###########################################################
compute_plugin_estimator <- function(Y, T, p_hat_func, X) {
  p_hat_vals <- p_hat_func(X)
  p_bar      <- mean(p_hat_vals)
  
  mu1_hat <- mean( (T * Y) / p_bar )
  mu0_hat <- mean( ((1 - T) * p_hat_vals * Y) / (1 - p_hat_vals) )
  
  tau_pi <- mu1_hat - mu0_hat
  return(tau_pi)
}

###########################################################
# Function to compute doubly robust estimator:
#   hat{tau}_DR = hat{mu}_1 - hat{mu}_0
#   where:
#   hat{mu}_1 = (1/n) * sum_{i} (T_i Y_i / hat{p})
#   hat{mu}_0 = (1/(n * hat{p})) * sum_{i} [ T_i hat{mu}_0(X_i) + ((1 - T_i)*hat{p}(X_i)*Y_i)/(1 - hat{p}(X_i)) ]
###########################################################
compute_DR_estimator <- function(Y, T, p_hat_func, mu0_hat_func, X) {
  p_hat_vals <- p_hat_func(X)
  p_bar      <- mean(p_hat_vals)
  
  mu0_hat_vals <- mu0_hat_func(X)
  
  mu1_hat <- mean((T * Y) / p_bar)
  term    <- T * mu0_hat_vals + ((1 - T) * p_hat_vals * Y)/(1 - p_hat_vals)
  mu0_hat <- (1/(p_bar)) * mean(term)
  
  tau_dr <- mu1_hat - mu0_hat
  return(tau_dr)
}

###########################################################
# Main simulation function:
# Runs simulations over combinations of (n, d, sBeta, sTheta, penalty_type),
# and repeats for multiple replicates to get individual estimates.
###########################################################
run_simulation <- function(
    ns = c(1000, 5000),
    ds = c(10, 50, 500, 5000),
    sparsity_levels = c("d/10", "d/2", "d"),
    penalty_types = c("lasso","ridge"),
    reps = 100,   # number of Monte Carlo reps per combination
    tau = 1,
    seed = 123
) {
  set.seed(seed)
  
  # Initialize an empty data.frame to store individual replication results
  results <- data.frame(
    n = integer(),
    d = integer(),
    sBeta = character(),
    sTheta = character(),
    penalty = character(),
    estimator = character(),
    estimate = numeric(),
    replicate = integer(),
    stringsAsFactors = FALSE
  )
  
  # Helper to parse sparsity strings like "d/10", "d/2", "d"
  sparsity_value <- function(sp, d) {
    if (sp == "d") {
      return(d)
    } else if (sp == "d/10") {
      return(max(1, floor(d/10)))  # Ensure at least 1
    } else if (sp == "d/2") {
      return(floor(d/2))
    } else {
      stop("Unknown sparsity level")
    }
  }
  
  # Total number of combinations for progress tracking (optional)
  total_combinations <- length(ns) * length(ds) * length(sparsity_levels) * length(sparsity_levels) * length(penalty_types)
  comb_counter <- 1
  
  # Loop over all combinations
  for (n_val in ns) {
    for (d_val in ds) {
      for (sBeta_str in sparsity_levels) {
        sBeta_val <- sparsity_value(sBeta_str, d_val)
        for (sTheta_str in sparsity_levels) {
          sTheta_val <- sparsity_value(sTheta_str, d_val)
          
          for (penalty_type in penalty_types) {
            
            # Optionally, skip combinations where sample size < dimension
            if (n_val >= d_val) {
              message(paste("Skipping combination: n =", n_val, ", d =", d_val, 
                            "(n < d)"))
              next
            }
            
            message(paste("Running combination", comb_counter, "of", total_combinations, 
                          ": n =", n_val, ", d =", d_val, 
                          ", sBeta =", sBeta_str, ", sTheta =", sTheta_str, 
                          ", penalty =", penalty_type))
            comb_counter <- comb_counter + 1
            
            # Run replicates
            for (r in 1:reps) {
              # Generate true parameters
              beta0 <- generate_parameters(d_val, sBeta_val, type = "beta")
              theta0 <- generate_parameters(d_val, sTheta_val, type = "theta")
              
              # Generate data
              data_list <- generate_data(n_val, d_val, beta0, theta0, tau = tau)
              X <- data_list$X
              Y <- data_list$Y
              T_vec <- data_list$T
              
              # Estimate models
              p_hat_func <- estimate_p(X, T_vec, penalty_type = penalty_type)
              mu0_hat_func <- estimate_mu0(X, Y, T_vec, penalty_type = penalty_type)
              
              # Compute estimators
              tau_pi <- compute_plugin_estimator(Y, T_vec, p_hat_func, X)
              tau_dr <- compute_DR_estimator(Y, T_vec, p_hat_func, mu0_hat_func, X)
              
              # Append plugin estimator result
              results <- rbind(
                results,
                data.frame(
                  n = n_val,
                  d = d_val,
                  sBeta = sBeta_str,
                  sTheta = sTheta_str,
                  penalty = penalty_type,
                  estimator = "plugin",
                  estimate = tau_pi,
                  replicate = r,
                  stringsAsFactors = FALSE
                )
              )
              
              # Append DR estimator result
              results <- rbind(
                results,
                data.frame(
                  n = n_val,
                  d = d_val,
                  sBeta = sBeta_str,
                  sTheta = sTheta_str,
                  penalty = penalty_type,
                  estimator = "DR",
                  estimate = tau_dr,
                  replicate = r,
                  stringsAsFactors = FALSE
                )
              )
              
              # Optional: Progress tracking within replicates
              if (r %% 1 == 0) {
                timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
                message(paste("  ", timestamp, "Completed replicate", r, "of", reps))
              }
            } # end replicates
            
          } # end penalty_types
        } # end sTheta
      } # end sBeta
    } # end d
  } # end n
  
  return(results)
}

###########################################################
# Run the simulation
###########################################################

# Example run with fewer reps for demonstration.
# **Important**: Running with large 'ds' and 'reps = 100' can be computationally intensive.
# Adjust 'reps' as needed based on your computational resources.

# generate random number between 0 and 1000
seed_opt <- sample(0:10000, 1)

sim_results <- run_simulation(
  ns = c(1000, 5000),
  ds = c(10, 50, 500, 5000),
  sparsity_levels = c("d/10", "d/2", "d"),
  penalty_types = c("lasso", "ridge"),
  reps = 1,  # Set to desired number of replications
  seed = seed_opt
)

###########################################################
# Save results to a LaTeX table and CSV
###########################################################

# Save as LaTeX table
# Note: Due to potentially large size, consider filtering or summarizing before converting to LaTeX.
# Here, we save the entire results, but it might not be practical for very large datasets.

# Uncomment the following lines to generate the LaTeX table
# out_table <- xtable(sim_results, 
#                     caption = "Simulation Results for ATT Estimators", 
#                     label = "tab:sim_results")
# print(out_table, file = "q2a_simulation_results.tex", include.rownames = FALSE)

# Save as CSV for easier data manipulation and analysis
# Write with the timestamp of the time and date before the csv extension
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

write.csv(sim_results, paste0("q2a-data/q2a_all_simulation_results_n_less_d_", timestamp, "_", seed_opt, ".csv"))

# Optional: View the first few rows of the results
head(sim_results)