########################################
# Load Libraries
########################################
library(dplyr)
library(glmnet)        # For LASSO
library(randomForest)  # For Random Forests
library(nnet)          # For neural networks (simple example)
# library(keras)       # For deep neural networks if desired
library(MASS)          # For parametric models, etc.

########################################
# Data Preparation
########################################

# Load data (assuming CSV file in working directory)
census <- read.csv("census2000.csv")

# Construct outcome and treatment
# Define treatment: T = 1 if female, 0 if male
census <- census %>% 
  mutate(T = ifelse(sex == "F", 1, 0),
         wage_rate = income / hours,
         Y = log(wage_rate))

# Define covariates X (excluding treatment and outcome)
# Here we include: age, marital, race, education as controls
X_vars <- c("age", "marital", "race", "education")
census <- census %>% 
  filter(!is.na(Y)) %>%  # Ensure no missing Y
  filter(!is.na(T))      # Ensure no missing T

X <- census[, X_vars]
T_ind <- census$T
Y <- census$Y
n <- nrow(census)

########################################
# Functions to Fit Models
########################################

# 1) Propensity Score Estimation
fit_propensity <- function(X, T_ind, method = "logistic") {
  if (method == "logistic") {
    # Parametric logistic regression
    model <- glm(T_ind ~ ., data = data.frame(T_ind, X), family = binomial)
    p_hat <- predict(model, newdata = data.frame(X), type = "response")
  } else if (method == "lasso") {
    # LASSO logistic regression
    x_mat <- model.matrix(~ . - 1, data = X)
    cvfit <- cv.glmnet(x_mat, T_ind, family = "binomial", alpha = 1)
    p_hat <- predict(cvfit, newx = x_mat, s = "lambda.min", type = "response")
  } else if (method == "forest") {
    # Random forest for classification
    rf_mod <- randomForest(as.factor(T_ind) ~ ., data = data.frame(T_ind=as.factor(T_ind), X))
    p_hat <- predict(rf_mod, newdata = X, type = "prob")[,2]
  } else if (method == "nn") {
    # Simple neural network for classification using nnet
    nn_mod <- nnet(as.factor(T_ind) ~ ., data = data.frame(T_ind=as.factor(T_ind), X), 
                   size = 5, maxit = 500, trace = FALSE)
    p_hat <- predict(nn_mod, newdata = X, type = "raw")
  } else {
    stop("Unknown propensity method")
  }
  return(p_hat)
}

# 2) Outcome Model for mu0(X)
#    We estimate E[Y | X, T=0]. This typically requires subsetting T=0 data for training.
fit_mu0 <- function(X, Y, T_ind, method = "linear") {
  # Train only on control units
  X0 <- X[T_ind == 0, , drop = FALSE]
  Y0 <- Y[T_ind == 0]
  
  if (method == "linear") {
    model <- lm(Y0 ~ ., data = data.frame(Y0, X0))
    mu0_hat <- predict(model, newdata = X)
  } else if (method == "lasso") {
    x_mat0 <- model.matrix(~ . - 1, data = X0)
    cvfit <- cv.glmnet(x_mat0, Y0, alpha = 1)
    x_mat_full <- model.matrix(~ . - 1, data = X)
    mu0_hat <- as.numeric(predict(cvfit, newx = x_mat_full, s = "lambda.min"))
  } else if (method == "forest") {
    rf_mod <- randomForest(Y0 ~ ., data = data.frame(Y0, X0))
    mu0_hat <- predict(rf_mod, newdata = X)
  } else if (method == "nn") {
    # Simple neural network for regression
    # Here we might standardize Y and X, or use a single hidden layer
    nn_mod <- nnet(Y0 ~ ., data = data.frame(Y0, X0), size = 5, linout = TRUE, trace = FALSE, maxit=500)
    mu0_hat <- predict(nn_mod, newdata = X)
  } else {
    stop("Unknown mu0 method")
  }
  
  return(mu0_hat)
}

########################################
# Estimators (PI and DR)
########################################

# Plug-in estimator of the ATT
# tau_PI = mu1_hat - mu0_hat
# where mu1_hat = (1/n)*sum(T_i*Y_i/p_hat)
# and mu0_hat = (1/n)*sum[(1 - T_i)*p_hat(X_i)*Y_i/(1 - p_hat(X_i))]
tau_PI <- function(Y, T_ind, p_hat) {
  mu1_hat <- mean((T_ind * Y) / mean(p_hat))
  mu0_hat <- mean(((1 - T_ind) * p_hat * Y) / (1 - p_hat))
  tau <- mu1_hat - mu0_hat
  return(tau)
}

# Doubly Robust estimator of the ATT
# tau_DR = [ (1/n)*sum(T_i*Y_i/p_hat) ] 
#          - (1/p_hat)* (1/n)*sum[ T_i*mu0_hat(X_i) + ((1 - T_i)*p_hat(X_i)*Y_i/(1 - p_hat(X_i))) ]
tau_DR <- function(Y, T_ind, p_hat, mu0_hat) {
  p_bar <- mean(T_ind)
  part1 <- mean((T_ind * Y) / p_bar)
  part2 <- mean( T_ind * mu0_hat + ((1 - T_ind) * p_hat * Y)/(1 - p_hat) )
  tau <- part1 - (1/p_bar)*part2
  return(tau)
}

########################################
# Cross-Fitting Utility
########################################

# Cross-fitting involves splitting the data into K folds, 
# fitting models on K-1 folds, and predicting on the held-out fold.
# This is just a conceptual example with K=2 for brevity.

cross_fit_estimates <- function(X, Y, T_ind, p_method, mu0_method, K=2) {
  set.seed(123) # For reproducibility
  fold_ids <- sample(rep(1:K, length.out = nrow(X)))
  tau_PI_vals <- numeric(K)
  tau_DR_vals <- numeric(K)
  
  for (k in 1:K) {
    train_idx <- which(fold_ids != k)
    test_idx  <- which(fold_ids == k)
    
    # Fit p and mu0 on training set
    p_hat_train <- fit_propensity(X[train_idx,], T_ind[train_idx], method = p_method)
    mu0_hat_train <- fit_mu0(X[train_idx,], Y[train_idx], T_ind[train_idx], method = mu0_method)
    
    # Now we need predictions on the test set.
    # For simplicity, we re-use the same fitting functions on the training set to get models, 
    # but we should store the models. Here, since our fit_ functions return predictions, 
    # we can just re-fit. A better approach would be to modify fit_ functions to return models too.
    
    # Refit the propensity model to train data and predict on test
    # We'll use a logistic model as an example. If using other methods, adapt accordingly.
    df_train <- data.frame(T_ind=T_ind[train_idx], X[train_idx,])
    # Ensure columns are properly named
    prop_mod <- glm(T_ind ~ ., data = df_train, family = binomial)
    p_hat_test <- predict(prop_mod, newdata = X[test_idx,], type = "response")
    
    # Refit the mu0 model on training data and predict on test data
    # Note: The mu0_fit functions returns predictions for full sample, let's do a direct lm here:
    df_train_mu0 <- data.frame(Y0=Y[train_idx], X[train_idx,])
    mu0_test_model <- lm(Y0 ~ ., data = df_train_mu0)
    mu0_test_predictions <- predict(mu0_test_model, newdata=X[test_idx,])
    
    # Compute PI and DR using test set predictions
    tau_PI_vals[k] <- tau_PI(Y[test_idx], T_ind[test_idx], p_hat_test)
    tau_DR_vals[k] <- tau_DR(Y[test_idx], T_ind[test_idx], p_hat_test, mu0_test_predictions)
  }
  
  return(list(PI = mean(tau_PI_vals), DR = mean(tau_DR_vals)))
}

########################################
# Main Analysis
########################################

# Example without cross-fitting using parametric models
p_hat_logistic <- fit_propensity(X, T_ind, method = "logistic")
mu0_hat_linear <- fit_mu0(X, Y, T_ind, method = "linear")

tau_pi_param <- tau_PI(Y, T_ind, p_hat_logistic)
tau_dr_param <- tau_DR(Y, T_ind, p_hat_logistic, mu0_hat_linear)

cat("Results (No Cross-Fit, Parametric):\n")
cat("PI Estimator:", tau_pi_param, "\n")
cat("DR Estimator:", tau_dr_param, "\n\n")

# Example with cross-fitting using the same methods
results_cf_param <- cross_fit_estimates(X, Y, T_ind, p_method = "logistic", mu0_method = "linear")
cat("Results (Cross-Fit, Parametric):\n")
cat("PI Estimator:", results_cf_param$PI, "\n")
cat("DR Estimator:", results_cf_param$DR, "\n\n")

# Example using LASSO for p and mu0 without cross-fitting
p_hat_lasso <- fit_propensity(X, T_ind, method = "lasso")
mu0_hat_lasso <- fit_mu0(X, Y, T_ind, method = "lasso")

tau_pi_lasso <- tau_PI(Y, T_ind, p_hat_lasso)
tau_dr_lasso <- tau_DR(Y, T_ind, p_hat_lasso, mu0_hat_lasso)

cat("Results (No Cross-Fit, LASSO):\n")
cat("PI Estimator:", tau_pi_lasso, "\n")
cat("DR Estimator:", tau_dr_lasso, "\n\n")

# Similarly, you can run for forests, neural nets, etc.
# Example: Random Forests
p_hat_forest <- fit_propensity(X, T_ind, method = "forest")
mu0_hat_forest <- fit_mu0(X, Y, T_ind, method = "forest")

tau_pi_forest <- tau_PI(Y, T_ind, p_hat_forest)
tau_dr_forest <- tau_DR(Y, T_ind, p_hat_forest, mu0_hat_forest)

cat("Results (No Cross-Fit, Random Forest):\n")
cat("PI Estimator:", tau_pi_forest, "\n")
cat("DR Estimator:", tau_dr_forest, "\n\n")

# For neural networks, you'd do similarly:
# p_hat_nn <- fit_propensity(X, T_ind, method = "nn")
# mu0_hat_nn <- fit_mu0(X, Y, T_ind, method = "nn")
# tau_pi_nn <- tau_PI(Y, T_ind, p_hat_nn)
# tau_dr_nn <- tau_DR(Y, T_ind, p_hat_nn, mu0_hat_nn)

# cat("Results (No Cross-Fit, Neural Network):\n")
# cat("PI Estimator:", tau_pi_nn, "\n")
# cat("DR Estimator:", tau_dr_nn, "\n\n")

########################################
# After obtaining all the results:
# - Compare PI vs DR
# - Compare with and without cross-fitting
# - Compare across different methods (parametric, LASSO, forest, NN)
# - Discuss the findings as requested.
########################################