# Load necessary packages

library(tidyverse)
library(glmnet)

# Function to generate data
generate_data <- function(n, d, beta_0, theta_0, tau, noise_sd) {
  # n: sample size
  # d: dimension of X
  # beta_0: true coefficients for Y(0)
  # theta_0: true coefficients for propensity score
  # tau: true treatment effect
  # noise_sd: standard deviation of the noise in Y(0)
  
  # Simulate covariates X ~ N(0, I_d)
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  # Compute true propensity scores
  p_true <- 1 / (1 + exp(-X %*% theta_0))
  
  # Ensure probabilities are within a reasonable range to avoid extreme values
  p_true <- pmin(pmax(p_true, 0.01), 0.99)
  
  # Generate treatment assignments T ~ Bernoulli(p_true)
  T <- rbinom(n, 1, prob = p_true)
  
  # Generate potential outcomes
  epsilon_0 <- rnorm(n, mean = 0, sd = noise_sd)
  Y0 <- X %*% beta_0 + epsilon_0
  Y1 <- Y0 + tau
  
  # Observed outcome
  Y <- T * Y1 + (1 - T) * Y0
  
  # Return a data frame
  data <- tibble(Y = as.vector(Y),
                 T = T,
                 X = I(X),
                 p_true = as.vector(p_true))
  
  return(data)
}

# Function to estimate propensity score using MLE with regularization (ridge regression)
estimate_ps_mle <- function(data) {
  library(glmnet)
  
  # Prepare data
  X <- as.matrix(data$X)
  T <- data$T
  
  # Fit logistic regression with ridge penalty
  cv_fit <- cv.glmnet(X, T, family = "binomial", alpha = 0)
  
  # Use the lambda that minimizes cross-validated error
  best_lambda <- cv_fit$lambda.min
  
  # Predict propensity scores
  p_hat <- predict(cv_fit, newx = X, s = best_lambda, type = "response")
  
  # Ensure estimated probabilities are within a reasonable range
  p_hat <- pmin(pmax(p_hat, 0.01), 0.99)
  
  return(p_hat)
}

# Function to estimate propensity score using NLS with optim()
estimate_ps_nls <- function(data) {
  # Prepare data for NLS
  X <- as.matrix(data$X)
  T <- data$T
  
  # Define logistic function
  logistic_function <- function(theta, X) {
    1 / (1 + exp(-X %*% theta))
  }
  
  # Define objective function (sum of squared errors)
  sse_function <- function(theta) {
    p_hat <- logistic_function(theta, X)
    sum((T - p_hat)^2)
  }
  
  # Initial guess for theta
  theta_start <- rep(0, ncol(X))
  
  # Use optim to minimize sse_function
  optim_result <- optim(par = theta_start,
                        fn = sse_function,
                        method = "BFGS",
                        control = list(maxit = 1000))
  
  # Estimated coefficients
  theta_hat <- optim_result$par
  
  # Estimated propensity scores
  p_hat <- logistic_function(theta_hat, X)
  
  # Ensure estimated probabilities are within a reasonable range
  p_hat <- pmin(pmax(p_hat, 0.01), 0.99)
  
  return(p_hat)
}

# Function to compute ATT estimator
compute_att <- function(data, p_hat) {
  # Number of treated and control units
  n_T <- sum(data$T)
  n_C <- sum(1 - data$T)
  
  # Avoid division by zero
  if (n_T == 0 || n_C == 0) {
    return(NA)
  }
  
  # Compute mu1_hat (average outcome among treated units)
  mu1_hat <- sum(data$T * data$Y) / n_T
  
  # Compute weights for control units
  weight <- p_hat / (1 - p_hat)
  
  # Compute mu0_hat (weighted average outcome among control units)
  mu0_hat <- sum((1 - data$T) * weight * data$Y) / sum((1 - data$T) * weight)
  
  # Compute ATT estimator
  tau_hat <- mu1_hat - mu0_hat
  
  return(tau_hat)
}

# Function to run a single simulation iteration
run_simulation <- function(n, d, beta_0, theta_0, tau, noise_sd) {
  # Generate data
  data <- generate_data(n, d, beta_0, theta_0, tau, noise_sd)
  
  # Estimate propensity scores using MLE
  p_hat_mle <- estimate_ps_mle(data)
  
  # Estimate propensity scores using NLS
  p_hat_nls <- estimate_ps_nls(data)
  
  # Compute ATT estimators
  tau_hat_mle <- compute_att(data, p_hat_mle)
  tau_hat_nls <- compute_att(data, p_hat_nls)
  
  # True ATT
  tau_true <- tau
  
  # Return results as a tibble
  result <- tibble(
    n = n,
    d = d,
    tau_true = tau_true,
    tau_hat_mle = tau_hat_mle,
    tau_hat_nls = tau_hat_nls
  )
  
  return(result)
}

# Parameters for the simulation
set.seed(42)  # For reproducibility
num_simulations <- 100
sample_sizes <- c(50, 100, 250, 500, 1000)
dimensions <- c(5, 10, 50, 100)
noise_levels <- c(.1, .3, .5)
beta_0 <- c(1, -1, 0.5, rep(0, max(dimensions) - 3))  # Adjust length as per d
theta_0 <- c(0.5, -0.5, 0.3, rep(0, max(dimensions) - 3))
tau <- 2  # True treatment effect

summary_stats <- tibble()

# Run simulations
for (n in sample_sizes) {
  for (d in dimensions) {
    simulations <- tibble()
    for (noise_sd in noise_levels) {
      # Adjust beta_0 and theta_0 to match current dimension d
      beta_0_current <- beta_0[1:d]
      theta_0_current <- theta_0[1:d]
      
      for (sim in 1:num_simulations) {
        result <- run_simulation(n, d, beta_0_current, theta_0_current, tau, noise_sd)
        result <- result %>%
          mutate(noise_sd = noise_sd,
                 sim = sim)
        simulations <- bind_rows(simulations, result)
        print(paste0("Simulations ", sim, " for noise: ", noise_sd, "; d:", d, "; n: ", n))
      }
    }
    new_summary_stats <- simulations %>%
      group_by(n, d, noise_sd) %>%
      summarise(
        tau_true = first(tau_true),
        bias_mle = mean(tau_hat_mle - tau_true),
        bias_nls = mean(tau_hat_nls - tau_true),
        sd_mle = sd(tau_hat_mle),
        sd_nls = sd(tau_hat_nls),
        mse_mle = mean((tau_hat_mle - tau_true)^2),
        mse_nls = mean((tau_hat_nls - tau_true)^2),
        n_simulations = n(),
        .groups = 'keep'
      )
    simulations %>% write_csv(paste0('q1g-data/simulations_d', d, '_n', n, '.csv'))
    summary_stats <- summary_stats %>% rbind(new_summary_stats)
    summary_stats %>% write_csv('q1g-data/simulations_summary_stats.csv')
  }
}