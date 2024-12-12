###########################################################
# Simulation Code for Nonparametric Low-Dimensional Case
#
# The code below outlines a modular simulation study in R.
#
# Changes from part (a):
#  - Nonparametric data-generating processes for mu0(x) and p(x).
#  - We now consider functions mu0 and p that are nonlinear.
#  - We use methods like random forests and deep neural networks to estimate mu0 and p.
#  - Consider sparsity in the sense that only a subset of the D covariates actually matter.
#
# Like before:
#  - We'll try sample sizes n = {1000, 5000, 15000}
#  - Dimensions d = {1, 3, 5, 10}
#  - For sparsity, we specify that only a subset s of the variables enter into the nonlinear functions.
#
# The code:
# 1. Defines functions for:
#    - Data generation (nonlinear mu0 and p)
#    - Fitting models (mu0_hat and p_hat) using random forest and deep nets
#    - Computing plugin and DR estimators
#    - Running the simulation across combinations of parameters and methods
#    - Saving results to a LaTeX table
#
# Notes:
# - Make sure necessary packages are installed: randomForest, keras, tensorflow, xtable
# - For deep nets, ensure keras is configured properly.
# - The code is a template and may require adjustments.
# - You do not need to run this code, just provide it.
#
###########################################################

library(randomForest)
library(keras)
library(xtable)

###########################################################
# Nonlinear functions for mu0 and p
#
# We'll define mu0 and p as follows:
# mu0(x): A nonlinear function of a subset of x. For example:
#   If only the first s variables matter:
#     mu0(x) = sin(x_1) + 0.5 * x_2^2 + ... (add complexity as desired)
#
# p(x): A nonlinear logistic function:
#   p(x) = 1 / (1 + exp(-g(x)))
# where g(x) is a nonlinear function, e.g.:
#   g(x) = 0.5 * sin(x_1) + (x_2^2)/4 - cos(x_3) + ...
#
# We'll allow s < d for sparsity. If s < d, then variables beyond s are ignored in mu0 and g.
###########################################################

mu0_function <- function(X, s) {
  # X: n-by-d matrix
  # Only first s columns matter.
  # Example nonlinear form:
  # mu0(x) = sum_{j=1 to s} [ sin(x_j) ] if s>0,
  # plus maybe some polynomial terms for variety.
  
  if (s == 0) {
    # If no active covariates, just return 0
    return(rep(0, nrow(X)))
  }
  
  # Let's construct a simple nonlinear combination:
  # mu0(x) = sin(x_1) + 0.5 * x_2^2 + ... (we'll vary depending on s)
  # For simplicity, let's just sum sin(x_j) for j=1,...,s
  val <- rep(0, nrow(X))
  for (j in 1:s) {
    val <- val + sin(X[,j])  # simple nonlinear term
  }
  return(val)
}

g_function <- function(X, s) {
  # Nonlinear function for p(x) inside logistic.
  # We'll do something similar:
  # g(x) = 0.5*sin(x_1) + (x_2^2)/4 - cos(x_3) + ...
  
  if (s == 0) {
    return(rep(0, nrow(X)))
  }
  
  val <- rep(0, nrow(X))
  if (s >= 1) val <- val + 0.5 * sin(X[,1])
  if (s >= 2) val <- val + (X[,2]^2)/4
  if (s >= 3) val <- val - cos(X[,3])
  # If s > 3, add more arbitrary nonlinearities:
  for (j in 4:s) {
    val <- val + 0.3 * (X[,j]^3)  # just an example of a different nonlinearity
  }
  return(val)
}

###########################################################
# Data generation function
# 
# Data generating process:
#   X ~ N(0, I_d)
#   p(x) = 1 / (1 + exp(-g(x)))
#   T ~ Bernoulli(p(x))
#   mu0(x) as defined
#   Y = T*(mu0(x) + tau) + (1-T)*mu0(x) + noise (noise ~ N(0,1))
#
# Sparsity: only first s variables matter. The rest do not enter mu0 or p.
###########################################################
generate_data <- function(n, d, s, tau = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- matrix(rnorm(n*d), nrow = n, ncol = d)
  mu0_vals <- mu0_function(X, s)
  
  g_vals <- g_function(X, s)
  p_x <- 1/(1+exp(-g_vals))
  
  T_ind <- rbinom(n, 1, p_x)
  Y <- T_ind*(mu0_vals + tau) + (1 - T_ind)*mu0_vals + rnorm(n,0,1)
  
  return(list(X = X, Y = Y, T = T_ind, p_true = p_x, mu0_true = mu0_vals))
}

###########################################################
# Estimation functions
#
# We provide two methods: randomForest and a deep neural net.
# Feel free to add more.
#
# We'll define generic estimation functions that take method arguments:
# For p(x): we fit a model to predict T from X.
# For mu0(x): we fit a model to predict Y from X where T=0.
###########################################################

fit_p_randomForest <- function(X, T) {
  rf_p <- randomForest(x = X, y = as.factor(T))
  p_hat_func <- function(newX) {
    pred_probs <- predict(rf_p, newX, type="prob")[,2]
    return(pred_probs)
  }
  return(p_hat_func)
}

fit_mu0_randomForest <- function(X, Y, T) {
  X0 <- X[T==0, , drop=FALSE]
  Y0 <- Y[T==0]
  rf_mu0 <- randomForest(x = X0, y = Y0)
  mu0_hat_func <- function(newX) {
    preds <- predict(rf_mu0, newX)
    return(preds)
  }
  return(mu0_hat_func)
}

# For neural nets using Keras:
# We'll create a simple feed-forward network.
# Note: This requires a working Keras/TensorFlow backend.
# For simplicity, we define a helper that builds a small MLP.

build_mlp <- function(input_dim) {
  model <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu", input_shape = input_dim) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_dense(units = 1, activation = "linear")
  model %>% compile(
    loss = "mse",
    optimizer = "adam"
  )
  return(model)
}

build_mlp_classification <- function(input_dim) {
  model <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu", input_shape = input_dim) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_dense(units = 1, activation = "sigmoid")
  model %>% compile(
    loss = "binary_crossentropy",
    optimizer = "adam",
    metrics = c("accuracy")
  )
  return(model)
}

fit_p_deepnet <- function(X, T, epochs=10, batch_size=32) {
  model <- build_mlp_classification(dim(X)[2])
  model %>% fit(X, T, epochs = epochs, batch_size = batch_size, verbose = 0)
  
  p_hat_func <- function(newX) {
    preds <- model %>% predict(newX)
    return(as.numeric(preds))
  }
  return(p_hat_func)
}

fit_mu0_deepnet <- function(X, Y, T, epochs=10, batch_size=32) {
  X0 <- X[T==0, , drop=FALSE]
  Y0 <- Y[T==0]
  
  model <- build_mlp(dim(X0)[2])
  model %>% fit(X0, Y0, epochs = epochs, batch_size = batch_size, verbose = 0)
  
  mu0_hat_func <- function(newX) {
    preds <- model %>% predict(newX)
    return(as.numeric(preds))
  }
  return(mu0_hat_func)
}

###########################################################
# Estimators
###########################################################
compute_plugin_estimator <- function(Y, T, p_hat_func, X) {
  p_hat_vals <- p_hat_func(X)
  p_bar      <- mean(p_hat_vals)
  
  mu1_hat <- mean( (T * Y) / p_bar )
  mu0_hat <- mean( ((1 - T) * p_hat_vals * Y) / (1 - p_hat_vals) )
  
  tau_pi <- mu1_hat - mu0_hat
  return(tau_pi)
}

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
# Main simulation function
#
# Runs simulations over combinations of (n, d, s), methods (rf, deepnet),
# and repeats for multiple replicates to get sampling distributions.
###########################################################
run_simulation <- function(
  ns = c(1000, 5000, 15000),
  ds = c(1, 3, 5, 10),
  sparsities = c(1,3,5,10),  # s must be <= d, represents how many covariates matter
  methods = c("rf","deepnet"),
  reps = 50,
  tau = 1,
  seed = 123
) {
  set.seed(seed)
  
  results <- data.frame(
    n = integer(0),
    d = integer(0),
    s = integer(0),
    method = character(0),
    estimator = character(0),
    mean_est = numeric(0),
    sd_est = numeric(0)
  )
  
  for (n_val in ns) {
    for (d_val in ds) {
      # Ensure s <= d
      valid_s <- sparsities[sparsities <= d_val]
      for (s_val in valid_s) {
        for (method_val in methods) {
          tau_pi_vals <- numeric(reps)
          tau_dr_vals <- numeric(reps)
          
          for (r in 1:reps) {
            data_list <- generate_data(n_val, d_val, s_val, tau = tau)
            X <- data_list$X
            Y <- data_list$Y
            T_vec <- data_list$T
            
            # Fit models for p and mu0
            if (method_val == "rf") {
              p_hat_func <- fit_p_randomForest(X, T_vec)
              mu0_hat_func <- fit_mu0_randomForest(X, Y, T_vec)
            } else if (method_val == "deepnet") {
              p_hat_func <- fit_p_deepnet(X, T_vec)
              mu0_hat_func <- fit_mu0_deepnet(X, Y, T_vec)
            } else {
              stop("Unknown method")
            }
            
            tau_pi_vals[r] <- compute_plugin_estimator(Y, T_vec, p_hat_func, X)
            tau_dr_vals[r] <- compute_DR_estimator(Y, T_vec, p_hat_func, mu0_hat_func, X)
          }
          
          # Summarize
          pi_mean <- mean(tau_pi_vals)
          pi_sd   <- sd(tau_pi_vals)
          dr_mean <- mean(tau_dr_vals)
          dr_sd   <- sd(tau_dr_vals)
          
          results <- rbind(
            results,
            data.frame(
              n = n_val,
              d = d_val,
              s = s_val,
              method = method_val,
              estimator = "plugin",
              mean_est = pi_mean,
              sd_est = pi_sd
            ),
            data.frame(
              n = n_val,
              d = d_val,
              s = s_val,
              method = method_val,
              estimator = "DR",
              mean_est = dr_mean,
              sd_est = dr_sd
            )
          )
        } # end method
      } # end s
    } # end d
  } # end n
  
  return(results)
}

###########################################################
# Run the simulation
###########################################################

# This might be very computationally heavy (especially with deepnets).
# Adjust reps or skip deepnet if this is too large.

sim_results <- run_simulation(
  ns = c(1000, 5000, 15000),
  ds = c(1, 3, 5, 10),
  sparsities = c(1, 3, 5, 10),
  methods = c("rf","deepnet"),
  reps = 50,
  seed = 123
)

###########################################################
# Save results to a LaTeX table
###########################################################
out_table <- xtable(sim_results, 
                    caption = "Nonparametric Nonlinear Simulation Results for ATT Estimators", 
                    label = "tab:sim_results_nonlinear")
print(out_table, file = "q2b_simulation_results_nonlinear.tex", include.rownames = FALSE)