library(tidyverse)

q3 <- read_csv("data_for_HW4.csv")

q3_analysis <- q3 %>% 
  mutate(x.1 = case_when(
    x.1 < quantile(x.1, .25) ~ "P0-P24",
    x.1 < quantile(x.1, .50) ~ "P25-P49",
    x.1 < quantile(x.1, .75) ~ "P50-P74",
    TRUE ~ "P75-P100",
  )) %>% 
  mutate(x.2 = case_when(
    x.2 < quantile(x.2, .25) ~ "P0-P24",
    x.2 < quantile(x.2, .50) ~ "P25-P49",
    x.2 < quantile(x.2, .75) ~ "P50-P74",
    TRUE ~ "P75-P100",
  )) %>% 
  mutate(x.3 = case_when(
    x.3 < quantile(x.3, .25) ~ "P0-P24",
    x.3 < quantile(x.3, .50) ~ "P25-P49",
    x.3 < quantile(x.3, .75) ~ "P50-P74",
    TRUE ~ "P75-P100",
  )) %>% 
  mutate(x.4 = case_when(
    x.4 < quantile(x.4, .25) ~ "P0-P24",
    x.4 < quantile(x.4, .50) ~ "P25-P49",
    x.4 < quantile(x.4, .75) ~ "P50-P74",
    TRUE ~ "P75-P100",
  )) %>% 
  mutate(x.5 = case_when(
    x.5 < quantile(x.5, .25) ~ "P0-P24",
    x.5 < quantile(x.5, .50) ~ "P25-P49",
    x.5 < quantile(x.5, .75) ~ "P50-P74",
    TRUE ~ "P75-P100",
  ))

q3_analysis %>% 
  group_by(e, x.1, x.2, x.3, x.4, x.5) %>% 
  summarise(t = mean(t))

# rbind(
#   q3_analysis %>% 
#     group_by(e, x.1) %>% 
#     summarise(t = mean(t)) %>% 
#     spread(x.1, t) %>% 
#     mutate(x = "x1"),
#   q3_analysis %>% 
#     group_by(e, x.2) %>% 
#     summarise(t = mean(t)) %>% 
#     spread(x.2, t) %>% 
#     mutate(x = "x2"),
#   q3_analysis %>% 
#     group_by(e, x.3) %>% 
#     summarise(t = mean(t)) %>% 
#     spread(x.3, t) %>% 
#     mutate(x = "x3"),
#   q3_analysis %>% 
#     group_by(e, x.4) %>% 
#     summarise(t = mean(t)) %>% 
#     spread(x.4, t) %>% 
#     mutate(x = "x4"),
#   q3_analysis %>% 
#     group_by(e, x.5) %>% 
#     summarise(t = mean(t)) %>% 
#     spread(x.5, t) %>% 
#     mutate(x = "x5")
# ) %>% 
#   arrange(e)

doubly_robust_ate_corrected <- function(data) {
  Y <- data$y
  T_ind <- data$t
  X <- data %>% select(x.1, x.2, x.3, x.4, x.5)
  
  ps_model <- glm(T_ind ~ x.1 + x.2 + x.3 + x.4 + x.5,data = data, family = binomial)
  p_hat <- predict(ps_model, type = "response")
  
  outcome_model <- lm(y ~ t * (x.1 + x.2 + x.3 + x.4 + x.5), data = data)
  
  data_T1 <- data
  data_T1$t <- 1
  mu1_hat <- predict(outcome_model, newdata = data_T1)
  
  data_T0 <- data
  data_T0$t <- 0
  mu0_hat <- predict(outcome_model, newdata = data_T0)
  
  tau_dr <- mean(mu1_hat - mu0_hat + (T_ind * (Y - mu1_hat) / p_hat) - 
                   ((1 - T_ind) * (Y - mu0_hat) / (1 - p_hat)))
  
  IF <- (mu1_hat - mu0_hat) + (T_ind * (Y - mu1_hat) / p_hat) - 
    ((1 - T_ind) * (Y - mu0_hat) / (1 - p_hat)) - tau_dr
  
  se_tau <- sd(IF) / sqrt(nrow(data))
  ci_lower <- tau_dr - 1.96 * se_tau
  ci_upper <- tau_dr + 1.96 * se_tau
  
  list(ATE = tau_dr, CI = c(ci_lower, ci_upper))
}

q3_e1 <- q3 %>% filter(e == 1)
q3_e2 <- q3 %>% filter(e == 2)

res_e1_corrected <- doubly_robust_ate_corrected(q3_e1)
res_e2_corrected <- doubly_robust_ate_corrected(q3_e2)

cat("Results for data source e=1 (doubly robust, corrected):\n")
cat("ATE:", res_e1_corrected$ATE, "\n")
cat("95% CI:", res_e1_corrected$CI[1], "to", res_e1_corrected$CI[2], "\n\n")

cat("Results for data source e=2 (doubly robust, corrected):\n")
cat("ATE:", res_e2_corrected$ATE, "\n")
cat("95% CI:", res_e2_corrected$CI[1], "to", res_e2_corrected$CI[2], "\n\n")

