# Load required packages
library(tidyverse)
library(broom)
library(moments)

# Read data
q3 <- read_csv("data_for_HW4.csv")

# Split by data source
q3_e1 <- q3 %>% filter(e == 1)
q3_e2 <- q3 %>% filter(e == 2)

######################################
# Fit linear models with interactions #
######################################

# Model with treatment and interactions for e=1
model_e1 <- lm(y ~ t*(x.1 + x.2 + x.3 + x.4 + x.5), data = q3_e1)
# Model with treatment and interactions for e=2
model_e2 <- lm(y ~ t*(x.1 + x.2 + x.3 + x.4 + x.5), data = q3_e2)

#########################################
# Compute the CATE for each observation #
#########################################

# For each model, CATE_i = beta_t + x_i * (beta_{t:x})
# Extract coefficients
coefs_e1 <- coef(model_e1)
coefs_e2 <- coef(model_e2)

# Create a function to compute CATE given a dataset and coefficients
compute_cate <- function(data, coefs) {
  # Extract covariates
  X_mat <- data %>% select(x.1, x.2, x.3, x.4, x.5)
  
  # Construct CATE = coefs["t"] + sum_j x_j[i]*coefs[paste0("t:x.j")]
  # Identify the relevant coefficient names for interaction terms
  int_terms <- names(coefs)[grepl("t:", names(coefs))]
  
  # Compute linear predictor for T=1 minus T=0 on each observation
  # T=0 part: y0 = (intercept + sum(b_x_j * x_j))
  # T=1 part: y1 = y0 + b_t + sum(b_{t:x_j} * x_j)
  # So CATE = (b_t + sum(b_{t:x_j} * x_j))
  
  cate <- coefs["t"]
  # Add interaction terms
  for (term in int_terms) {
    # term looks like "t:x.j"
    varname <- sub("t:", "", term)
    cate <- cate + X_mat[[varname]] * coefs[term]
  }
  return(cate)
}

cate_e1 <- compute_cate(q3_e1, coefs_e1)
cate_e2 <- compute_cate(q3_e2, coefs_e2)

######################################
# Plot the distribution of the CATEs #
######################################

# Plot distributions
par(mfrow = c(1,2))
hist(cate_e1, main = "CATE Distribution (e=1)", xlab = "CATE", breaks = 30)
hist(cate_e2, main = "CATE Distribution (e=2)", xlab = "CATE", breaks = 30)

tibble(
  stats = c(
    "Mean",
    "Std",
    "Skewness",
    "Kurtosis",
    "Q10",
    "Q25",
    "Median",
    "Q75",
    "Q90"
  ),
  e1 = c(
    cate_e1 %>% mean(),
    cate_e1 %>% sd(),
    cate_e1 %>% skewness(),
    cate_e1 %>% kurtosis(),
    cate_e1 %>% quantile(.1),
    cate_e1 %>% quantile(.25),
    cate_e1 %>% quantile(.5),
    cate_e1 %>% quantile(.75),
    cate_e1 %>% quantile(.9)
  ),
  e2 = c(
    cate_e2 %>% mean(),
    cate_e2 %>% sd(),
    cate_e2 %>% skewness(),
    cate_e2 %>% kurtosis(),
    cate_e2 %>% quantile(.1),
    cate_e2 %>% quantile(.25),
    cate_e2 %>% quantile(.5),
    cate_e2 %>% quantile(.75),
    cate_e2 %>% quantile(.9)
  )
)

########################################
# Compute the ATE and its confidence CI #
########################################

# The ATE is the expectation of the CATE. Since we have a linear model:
# ATE = average of CATE_i = mean over i of [b_t + sum_j x_j[i]*b_{t:x_j}]
# This is the same as: ATE = b_t + sum_j [E(x_j)*b_{t:x_j}]

# Let's compute ATE directly as mean of predicted CATE (empirical average)
ate_e1 <- mean(cate_e1)
ate_e2 <- mean(cate_e2)

# To get a confidence interval for the ATE, we need the variance estimate.
# Since ATE is a linear combination of coefficients, we can use the delta method:
# Var(ATE) = Var(b_t + sum_j mean(x_j)*b_{t:x_j})
# Let's compute for each dataset.

get_ate_ci <- function(model, data) {
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  
  # Mean of x variables
  x_means <- colMeans(data[,c("x.1","x.2","x.3","x.4","x.5")])
  
  # Identify t and t:x coefficients
  b_t <- coefs["t"]
  t_int_terms <- names(coefs)[grepl("t:", names(coefs))]
  
  # Compute ATE = b_t + sum_j mean(x_j)*b_{t:x_j}
  ATE <- b_t
  # Create a vector of partial derivatives for ATE w.r.t. the parameters
  grad <- rep(0, length(coefs))
  names(grad) <- names(coefs)
  grad["t"] <- 1
  
  for (term in t_int_terms) {
    varname <- sub("t:", "", term)
    ATE <- ATE + x_means[varname]*coefs[term]
    grad[term] <- x_means[varname]
  }
  
  # Variance of ATE
  var_ATE <- t(grad) %*% vcov_mat %*% grad
  se_ATE <- sqrt(var_ATE)
  
  # 95% CI
  CI <- c(ATE - 1.96*se_ATE, ATE + 1.96*se_ATE)
  
  return(list(ATE = ATE, CI = CI))
}

# Compute ATE and CI for e=1
res_e1 <- get_ate_ci(model_e1, q3_e1)

# Compute ATE and CI for e=2
res_e2 <- get_ate_ci(model_e2, q3_e2)

# Print results
cat("Results for data source e=1:\n")
cat("ATE:", res_e1$ATE, "\n")
cat("95% CI:", res_e1$CI[1], "to", res_e1$CI[2], "\n\n")

cat("Results for data source e=2:\n")
cat("ATE:", res_e2$ATE, "\n")
cat("95% CI:", res_e2$CI[1], "to", res_e2$CI[2], "\n\n")

# Compare the ATEs here to those obtained when we did not control for the x variables.
# The previous part had simple difference-in-means ATEs, and now we have adjusted ATEs
# that incorporate the interaction terms, providing a conditional treatment effect model.