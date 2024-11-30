library(tidyverse)
library(knitr)


simulations <- read_csv("q1g-data/simulations_summary_stats.csv")

n_opts <- simulations$n %>% unique()
noise_sd_opts <- simulations$noise_sd %>% unique()

for (noise_sd_opt in noise_sd_opts) {
  for (n_opt in n_opts){
    simulations %>% 
      filter(noise_sd == noise_sd_opt, n == n_opt) %>% 
      select(n, d, noise_sd, bias_mle, bias_nls) %>% 
      rename("Bias MLE" = bias_mle, "Bias NLS" = bias_nls) %>% 
      gather(id, value, -n, -d, -noise_sd) %>% 
      ggplot(aes(d, value)) +
      geom_line(aes(color = id)) +
      geom_point(aes(color = id)) +
      theme_minimal() +
      labs(
        title = 'Bias MLE and NLS by X Dimension with Noise SD: ' %>%
          paste0(noise_sd, " and N: ", n),
        x = 'd (X Dimension)',
        y = 'Bias'
      ) +
      theme(legend.position = 'bottom') +
      scale_color_manual(
        name =  '',
        values = c('red', 'blue')
      )
    ggsave('q1g-plots/bias_per_d' %>% paste0('_n', n_opt, '_noise', noise_sd_opt * 10, '.png'), width = 8, height = 5)
    
    simulations %>% 
      filter(noise_sd == noise_sd_opt, n == n_opt) %>% 
      select(n, d, noise_sd, mse_mle, mse_nls) %>% 
      rename("MSE MLE" = mse_mle, "MSE NLS" = mse_nls) %>% 
      gather(id, value, -n, -d, -noise_sd) %>% 
      ggplot(aes(d, value)) +
      geom_line(aes(color = id)) +
      geom_point(aes(color = id)) +
      theme_minimal() +
      labs(
        title = 'MSE MLE and NLS by X Dimension',
        x = 'd (X Dimension)',
        y = 'Bias'
      ) +
      theme(legend.position = 'bottom') +
      scale_color_manual(
        name =  '',
        values = c('red', 'blue')
      )
    ggsave('q1g-plots/mse_per_d' %>% paste0('_n', n_opt, '_noise', noise_sd_opt * 10, '.png'), width = 8, height = 5)
    
    simulations %>% 
      filter(noise_sd == noise_sd_opt, n == n_opt) %>% 
      select(n, d, noise_sd, sd_mle, sd_nls) %>% 
      rename("Std MLE" = sd_mle, "Std NLS" = sd_nls) %>% 
      gather(id, value, -n, -d, -noise_sd) %>% 
      ggplot(aes(d, value)) +
      geom_line(aes(color = id)) +
      geom_point(aes(color = id)) +
      theme_minimal() +
      labs(
        title = 'Std MLE and NLS by X Dimension',
        x = 'd (X Dimension)',
        y = 'Bias'
      ) +
      theme(legend.position = 'bottom') +
      scale_color_manual(
        name =  '',
        values = c('red', 'blue')
      )
    ggsave('q1g-plots/std_per_d' %>% paste0('_n', n_opt, '_noise', noise_sd_opt * 10, '.png'), width = 8, height = 5)
  }
}

simulations %>% 
  group_by(d) %>% 
  select(-n, -noise_sd, -tau_true, -n_simulations) %>% 
  summarise_all(., ~ mean(.)) %>% 
  purrr::set_names("Dimensions", "Bias MLE", "Bias NLS", "SD MLE", "SD NLS", "MSE MLE", "MSE NLS") %>% 
  kable(format = "latex", booktabs = TRUE) %>% 
  writeLines("q1g-plots/stats_by_d.tex")

simulations %>% 
  group_by(noise_sd) %>% 
  select(-n, -d, -tau_true, -n_simulations) %>% 
  summarise_all(., ~ mean(.)) %>% 
  purrr::set_names("Noise", "Bias MLE", "Bias NLS", "SD MLE", "SD NLS", "MSE MLE", "MSE NLS") %>% 
  kable(format = "latex", booktabs = TRUE) %>% 
  writeLines("q1g-plots/stats_by_noise_sd.tex")

simulations %>% 
  group_by(n) %>% 
  select(-d, -noise_sd, -tau_true, -n_simulations) %>% 
  summarise_all(., ~ mean(.)) %>% 
  purrr::set_names("N Samples", "Bias MLE", "Bias NLS", "SD MLE", "SD NLS", "MSE MLE", "MSE NLS") %>% 
  kable(format = "latex", booktabs = TRUE) %>% 
  writeLines("q1g-plots/stats_by_n.tex")
  

