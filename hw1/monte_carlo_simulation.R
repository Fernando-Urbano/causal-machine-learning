# Set seed for reproducibility
set.seed(91)

# Simulation parameters
n_sim <- 1000         # Number of simulations
n <- 1000             # Sample size per simulation
p_S <- 0.5            # Probability of S=1
p_T <- 0.5            # Probability of T=1

# True parameters
beta_0 <- 0
beta_S <- 1
beta_T <- 2
beta_ST <- 0.5
gamma <- 1
sigma <- 1

# Initialize vectors to store estimates
ATE_T_reg <- numeric(n_sim)
ATE_T_dim <- numeric(n_sim)

for (i in 1:n_sim) {
  # Generate treatments
  S <- rbinom(n, 1, p_S)
  T <- rbinom(n, 1, p_T)
  
  # Generate covariate
  X <- rnorm(n, 0, 1)
  
  # Generate outcome
  epsilon <- rnorm(n, 0, sigma)
  Y <- beta_0 + beta_S * S + beta_T * T + beta_ST * (S * T) + gamma * X + epsilon
  
  # Regression-based estimator
  model <- lm(Y ~ S + T + S:T)
  beta_hat <- coef(model)
  S_bar <- mean(S)
  ATE_T_reg[i] <- beta_hat["T"] + beta_hat["S:T"] * S_bar
  
  # Difference-in-means estimator
  Y_T1 <- Y[T == 1]
  Y_T0 <- Y[T == 0]
  ATE_T_dim[i] <- mean(Y_T1) - mean(Y_T0)
}

# Calculate true ATE_T
# E[Y(T=1)] - E[Y(T=0)] = (beta_T + beta_ST * E[S]) - (0) = beta_T + beta_ST * p_S
true_ATE_T <- beta_T + beta_ST * p_S

# Calculate Bias
bias_reg <- mean(ATE_T_reg) - true_ATE_T
bias_dim <- mean(ATE_T_dim) - true_ATE_T

# Calculate Variance
var_reg <- var(ATE_T_reg)
var_dim <- var(ATE_T_dim)

# Calculate Mean Squared Error
mse_reg <- mean((ATE_T_reg - true_ATE_T)^2)
mse_dim <- mean((ATE_T_dim - true_ATE_T)^2)

# Summary of results
results <- data.frame(
  Estimator = c("Regression-Based", "Difference-in-Means"),
  Bias = c(bias_reg, bias_dim),
  Variance = c(var_reg, var_dim),
  MSE = c(mse_reg, mse_dim)
)

print(results)

# Consistency
n_max <- 2450

ATE_T_dim <- numeric(n_max)
ATE_T_reg <- numeric(n_max)

for (n in 1:n_max){
    # Generate treatments
    S <- rbinom(n*10, 1, p_S)
    T <- rbinom(n*10, 1, p_T)
    
    # Generate covariate
    X <- rnorm(n*10, 0, 1)
    
    # Generate outcome
    epsilon <- rnorm(n*10, 0, sigma)
    Y <- beta_0 + beta_S * S + beta_T * T + beta_ST * (S * T) + gamma * X + epsilon
    
    # Regression-based estimator
    model <- lm(Y ~ S + T + S:T)
    beta_hat <- coef(model)
    S_bar <- mean(S)
    ATE_T_reg[n] <- beta_hat["T"] + beta_hat["S:T"] * S_bar
    
    # Difference-in-means estimator
    Y_T1 <- Y[T == 1]
    Y_T0 <- Y[T == 0]
    ATE_T_dim[n] <- mean(Y_T1) - mean(Y_T0)
}

tibble(
  "ATE Dim" = ATE_T_dim,
  "ATE Reg" = ATE_T_reg
) %>% 
  mutate(n=(seq(1, n())) * 10) %>% 
  gather(id, value, -n) %>% 
  ggplot(aes(n, value)) +
  geom_line(aes(color=id)) +
  geom_hline(aes(yintercept=2.25, color="True ATE"), size=1) +
  theme_minimal() +
  scale_color_manual(name = "", values = c(
    "ATE Dim" = "blue", "ATE Reg" = "red", "True ATE" = "black"
  )) +
  theme(legend.position = "bottom") +
  labs(
    title = "Consistence of ATE",
    y="ATE Estimate in Simulation",
    x="Sample Size (n)"
  )
ggsave("monte_carlo_ate_consistence.png", width = 8, height = 5, dpi = 300)
