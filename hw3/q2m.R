library(MASS)
library(stats)
library(tidyverse)
library(kableExtra)
library(xtable)
library(scales)

calc_p_x <- function(data, x, h) {
  return(length(data[(data > x) & (data <= x + h)]) / length(data))
}

calc_f_x <- function(data, x, h) {
  return(calc_p_x(data, x, h) / h)
}

calc_var_f_x <- function(data, x, h) {
  p <- calc_p_x(data, x, h)
  return((p * (1 - p)) / (length(data) * (h**2)))
}

build_ci <- function(ci_level, variance, mean, n) {
  alpha <- 1 - ci_level
  critical_value <- qt(1 - alpha / 2, df = n - 1)
  std_error <- sqrt(variance) * critical_value
  return(c(mean - std_error, mean + std_error))
}

prod_simulation <- function(distribution_name, true_pdf, data, x, h, ci_level = 0.95) {
  sample_pdf <- calc_f_x(data, x, h)
  sample_var_pdf <- calc_var_f_x(data, x, h)
  ci <- build_ci(.95, sample_var_pdf, sample_pdf, length(data))
  return(tibble(
    'distribution_name' = distribution_name,
    'x' = x, 'h' = h, 'n' = length(data),
    'true_pdf' = true_pdf,
    'sample_pdf' = sample_pdf, 'sample_var_pdf' = sample_var_pdf,
    'ci_lower_bound' = ci[1], 'ci_upper_bound' = ci[2],
    'ci_level' = ci_level
  ))
}

N_SIM <- 1000
X <- 0.5

set.seed(42)

N_OPS = c(10, seq(25, 75, 25), seq(100, 1000, 100))

VAR_EPSILON = 0.0001

H_VALUES = (1 / (N_OPS**(1/3)) + VAR_EPSILON)[c(1:8)]

for (n in N_OPS) {
  simulations <- tibble()
  for (h in H_VALUES) {
    for (i in 1:N_SIM) {
      dist_information <- list(
        list('T-Student(30)', dt(X, df = 30), rt(n, df = 30)),
        list('T-Student(50)', dt(X, df = 50), rt(n, df = 50)),
        list('T-Student(100)', dt(X, df = 100), rt(n, df = 100)),
        list('Normal(0, 1)', dnorm(X, mean = 0, sd = 1), rnorm(n, mean = 0, sd = 1)),
        list('Normal(1, 1)', dnorm(X, mean = 1, sd = 1), rnorm(n, mean = 1, sd = 1)),
        list('Gamma(1, 1)', dgamma(X, shape = 1, rate = 1), rgamma(n, shape = 1, rate = 1)),
        list('Gamma(2, 3)', dgamma(X, shape = 2, rate = 3), rgamma(n, shape = 2, rate = 3)),
        list('Log-Normal(1, 1)', dlnorm(X, meanlog = 1, sdlog = 1), rlnorm(n, meanlog = 1, sdlog = 1))
      )
      for (dist_info in dist_information) {
        simulations <- simulations %>%
          rbind(prod_simulation(
            distribution_name = dist_info[[1]], dist_info[[2]], dist_info[[3]], x = X, h = h
          ))
      }
      print(paste0("Simulation: ", i, "; n: ", n, "; h: ", h))
    }
  }
  simulations %>% write_csv(paste0('q2m-data/x', X*10, '/simulations_n', n, '_q2m.csv'))
}

for (x in c(.5, 1)) {
  simulation_files <- list.files(path = paste0("q2m-data/x", X * 10, "/"), pattern = "\\.csv$", full.names = TRUE)
  
  simulations <- tibble()
  for (simulation_file in simulation_files) {
    simulations <- simulations %>% 
      rbind(read.csv(simulation_file))
  }
  
  simulations_agg <- simulations %>% 
    mutate(inside_ci = true_pdf >= ci_lower_bound & true_pdf <= ci_upper_bound) %>% 
    group_by(x, h, n, distribution_name) %>% 
    summarise(
      pct_inside_ci = mean(inside_ci),
      n_sim = n(),
      .groups = 'keep'
    ) %>% 
    ungroup()
  
  distributions <- list(
    list('T-Student(30)', 't_student_30'),
    list('T-Student(50)', 't_student_50'),
    list('T-Student(100)', 't_student_100'),
    list('Normal(0, 1)', 'normal_0_1'),
    list('Normal(1, 1)', 'normal_1_1'),
    list('Gamma(1, 1)', 'gamma_1_1'),
    list('Gamma(2, 3)', 'gamma_2_3'),
    list('Log-Normal(1, 1)', 'log_normal_1_1')
  )
  
  color_by_value <- function(value) {
    if (is.na(value)) return("")
    scaled_value <- abs(value - 0.95) ** (1/2)
    color <- gradient_n_pal(c("red", "white", "green"))(1 - scaled_value / 0.95)
    sprintf("\\cellcolor[HTML]{%s} %.0f", substring(color, 2), value * 100)
  }
  
  for (distribution in distributions) {
    
    simulation_dist_agg <- simulations_agg %>% 
      filter(distribution_name == distribution[[1]]) %>% 
      dplyr::select(n, h, pct_inside_ci) %>% 
      spread(n, pct_inside_ci)
    
    formatted_table <- simulation_dist_agg %>%
      mutate(h = sprintf("\\textbf{%.2f}", h)) %>%
      mutate(across(-h, ~sapply(., color_by_value)))
    
    latex_table <- xtable(formatted_table, align = "c|l|" %>% paste0(paste(rep("c", 5), collapse = ""), "|"))
    
    latex_file <- paste0("q2m-data/simulations_f_x", X * 10, "_", distribution[[2]], "_comp_h_n", ".tex")
    sink(latex_file)
    print(latex_table, include.rownames = FALSE, sanitize.text.function = identity)
    sink()
  }
}