# Load required packages
library(tidyverse)
library(randomForest)
library(gbm)

# Read the data
q3 <- read_csv("data_for_HW4.csv")

# Split the data into the two sources
q3_e1 <- q3 %>% filter(e == 1)
q3_e2 <- q3 %>% filter(e == 2)

# Define a function to perform Doubly Robust ATE estimation using flexible models
doubly_robust_ate_ml <- function(data, 
                                 ps_model_method = "randomForest", 
                                 outcome_model_method = "randomForest") {
  # Extract variables
  Y <- data$y
  T_ind <- data$t
  X <- data %>% select(x.1, x.2, x.3, x.4, x.5)
  
  # Step 1: Estimate propensity score p(x) using specified method
  if (ps_model_method == "randomForest") {
    ps_model <- randomForest(as.factor(T_ind) ~ ., data = X, ntree = 100)
    p_hat <- predict(ps_model, newdata = X, type = "prob")[,2]
  } else if (ps_model_method == "gbm") {
    ps_model <- gbm(as.factor(T_ind) ~ ., data = cbind(X, T_ind = as.factor(T_ind)), 
                    distribution = "bernoulli", n.trees = 100, interaction.depth = 3, 
                    shrinkage = 0.1, verbose = FALSE)
    p_hat <- predict(ps_model, newdata = X, n.trees = 100, type = "response")
  } else {
    stop("Unsupported propensity score model method.")
  }
  
  # Ensure propensity scores are bounded away from 0 and 1 to avoid division by zero
  p_hat <- pmin(pmax(p_hat, 0.01), 0.99)
  
  # Step 2: Estimate outcome models mu_1(x) and mu_0(x) using specified method
  if (outcome_model_method == "randomForest") {
    # Model for treated (T=1)
    outcome_model_t1 <- randomForest(y ~ ., data = X %>% mutate(y = Y, t = T_ind), 
                                     subset = T_ind == 1, ntree = 100)
    mu1_hat <- predict(outcome_model_t1, newdata = X)
    
    # Model for control (T=0)
    outcome_model_t0 <- randomForest(y ~ ., data = X %>% mutate(y = Y, t = T_ind), 
                                     subset = T_ind == 0, ntree = 100)
    mu0_hat <- predict(outcome_model_t0, newdata = X)
    
  } else if (outcome_model_method == "gbm") {
    # Model for treated (T=1)
    outcome_model_t1 <- gbm(y ~ ., data = X %>% mutate(y = Y, t = T_ind), 
                            distribution = "gaussian", n.trees = 100, interaction.depth = 3, 
                            shrinkage = 0.1, verbose = FALSE, subset = T_ind == 1)
    mu1_hat <- predict(outcome_model_t1, newdata = X, n.trees = 100)
    
    # Model for control (T=0)
    outcome_model_t0 <- gbm(y ~ ., data = X %>% mutate(y = Y, t = T_ind), 
                            distribution = "gaussian", n.trees = 100, interaction.depth = 3, 
                            shrinkage = 0.1, verbose = FALSE, subset = T_ind == 0)
    mu0_hat <- predict(outcome_model_t0, newdata = X, n.trees = 100)
    
  } else {
    stop("Unsupported outcome model method.")
  }
  
  # Step 3: Compute the Doubly Robust ATE
  tau_dr <- mean((mu1_hat - mu0_hat) + 
                   (T_ind * (Y - mu1_hat) / p_hat) - 
                   ((1 - T_ind) * (Y - mu0_hat) / (1 - p_hat)))
  
  # Step 4: Compute the influence function for variance estimation
  IF <- (mu1_hat - mu0_hat) + 
    (T_ind * (Y - mu1_hat) / p_hat) - 
    ((1 - T_ind) * (Y - mu0_hat) / (1 - p_hat)) - tau_dr
  
  # Compute standard error and 95% confidence interval
  se_tau <- sd(IF) / sqrt(nrow(data))
  ci_lower <- tau_dr - 1.96 * se_tau
  ci_upper <- tau_dr + 1.96 * se_tau
  
  list(ATE = tau_dr, CI = c(ci_lower, ci_upper))
}

# Apply the Doubly Robust Estimator using different model combinations
# Example 1: Propensity Score - Random Forest; Outcome Models - Random Forest
res_e1_rf <- doubly_robust_ate_ml(q3_e1, 
                                  ps_model_method = "randomForest", 
                                  outcome_model_method = "randomForest")

res_e2_rf <- doubly_robust_ate_ml(q3_e2, 
                                  ps_model_method = "randomForest", 
                                  outcome_model_method = "randomForest")

# Example 2: Propensity Score - GBM; Outcome Models - GBM
res_e1_gbm <- doubly_robust_ate_ml(q3_e1, 
                                   ps_model_method = "gbm", 
                                   outcome_model_method = "gbm")

res_e2_gbm <- doubly_robust_ate_ml(q3_e2, 
                                   ps_model_method = "gbm", 
                                   outcome_model_method = "gbm")

# Example 3: Propensity Score - Random Forest; Outcome Models - GBM
res_e1_rf_gbm <- doubly_robust_ate_ml(q3_e1, 
                                      ps_model_method = "randomForest", 
                                      outcome_model_method = "gbm")

res_e2_rf_gbm <- doubly_robust_ate_ml(q3_e2, 
                                      ps_model_method = "randomForest", 
                                      outcome_model_method = "gbm")

# Example 4: Propensity Score - GBM; Outcome Models - Random Forest
res_e1_gbm_rf <- doubly_robust_ate_ml(q3_e1, 
                                      ps_model_method = "gbm", 
                                      outcome_model_method = "randomForest")

res_e2_gbm_rf <- doubly_robust_ate_ml(q3_e2, 
                                      ps_model_method = "gbm", 
                                      outcome_model_method = "randomForest")

# Compile all results into a tibble for comparison
results <- tibble(
  Data_Source = rep(c("e=1", "e=2"), each = 4),
  PS_Model = rep(c("Random Forest", "GBM", "Random Forest", "GBM"), times = 2),
  Outcome_Model = rep(c("Random Forest", "GBM", "GBM", "Random Forest"), times = 2),
  ATE = c(res_e1_rf$ATE, res_e1_gbm$ATE, res_e1_rf_gbm$ATE, res_e1_gbm_rf$ATE,
          res_e2_rf$ATE, res_e2_gbm$ATE, res_e2_rf_gbm$ATE, res_e2_gbm_rf$ATE),
  CI_Lower = c(res_e1_rf$CI[1], res_e1_gbm$CI[1], res_e1_rf_gbm$CI[1], res_e1_gbm_rf$CI[1],
               res_e2_rf$CI[1], res_e2_gbm$CI[1], res_e2_rf_gbm$CI[1], res_e2_gbm_rf$CI[1]),
  CI_Upper = c(res_e1_rf$CI[2], res_e1_gbm$CI[2], res_e1_rf_gbm$CI[2], res_e1_gbm_rf$CI[2],
               res_e2_rf$CI[2], res_e2_gbm$CI[2], res_e2_rf_gbm$CI[2], res_e2_gbm_rf$CI[2])
)

# Print the compiled results
print(results)

# Optional: Visualize the ATE estimates with confidence intervals
results %>%
  mutate(Model = paste("PS:", PS_Model, "| Outcome:", Outcome_Model)) %>%
  ggplot(aes(x = Model, y = ATE, ymin = CI_Lower, ymax = CI_Upper, color = Data_Source)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(position = position_dodge(width = 0.7), width = 0.2) +
  labs(title = "Doubly Robust ATE Estimates with Flexible Models",
       x = "Model Combination",
       y = "ATE Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
