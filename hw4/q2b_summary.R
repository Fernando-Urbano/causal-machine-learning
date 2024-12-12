# check all the files inside the folder q2a-data
library(tidyverse)
library(knitr)
library(kableExtra)
# Read all the files in the folder
files <- list.files("q2b-data", full.names = TRUE)

simulations <- tibble()

for (file in files) {
  simulations <- simulations %>% 
    rbind(read_csv(file))
}

simulations <- simulations %>%
  mutate(
    plugin_error = abs(plugin_mean - 1),
    dr_error = abs(dr_mean - 1),
    n_div_d = n / d
  ) %>% 
  dplyr::select(-dr_mean, -plugin_mean)

simulations_summary <- simulations %>% 
  group_by(n, d, s, method, n_div_d) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), std_plugin_error = sd(plugin_error, na.rm = TRUE),
  ) %>% 
  purrr::set_names("N", "D", "S", "Estimator", "N / D", "Avg DR Error", "Std DR Error", "Avg Plug-In Error", "Std Plug-In Error")

n <- 2

for (i in 1:n) {
  
  # divide the simulations_summary in 5 chunks
  n_rows <- nrow(simulations_summary)
  chunk_size <- n_rows %/% n
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, n_rows)
  
  simulations_summary_chunk <- simulations_summary[start:end, ]
  
  # Create LaTeX table with kable and style it with kableExtra
  latex_table <- simulations_summary_chunk %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      caption = "Simulation Results Part " %>% paste0(i),
      label = "tab:simulation_results_part" %>% paste0(i),
      longtable = TRUE,             # Use longtable for tables that span multiple pages
      linesep = "",                  # Remove extra lines between rows
      escape = FALSE                 # Allow LaTeX commands in the table
    ) %>%
    # Optionally, you can add additional styling or formatting here
    # For example, aligning numeric columns to the right
    column_spec(1, bold = TRUE) %>%  # Bold the first column (n)
    column_spec(2, width = "3cm")    # Adjust the width of the second column (d)
  
  # Define the file path where the LaTeX table will be saved
  output_file <- "report/q2b_simulations_summary_part" %>% paste0(i) %>% paste0(".tex")
  
  # Write the LaTeX table to the file
  writeLines(latex_table, con = output_file)
  
  # Optional: Print a message upon successful saving
  cat("LaTeX table saved to", output_file, "\n")
  
}

simulations_summary <- simulations %>% 
  group_by(n, d, s, method, n_div_d) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  )

# Ensure the output directory exists
if(!dir.exists("q2b-plots")) {
  dir.create("q2b-plots")
}

# 1. Average Error of the Estimators per N
simulations %>% 
  group_by(n) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = n, y = avg_dr_error, color = "Avg DR Error")) +
  geom_line(aes(x = n, y = avg_plugin_error, color = "Avg Plug-In Error")) +
  geom_point(aes(x = n, y = avg_dr_error, color = "Avg DR Error")) +
  geom_point(aes(x = n, y = avg_plugin_error, color = "Avg Plug-In Error")) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators per N',
    x = 'Sample Size (N)',
    y = 'Average Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_estimator_and_n.png', width = 8, height = 5)

# 2. Standard Error of the Estimators per N
simulations %>% 
  group_by(n) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = n, y = std_dr_error, color = "Std DR Error")) +
  geom_line(aes(x = n, y = std_plugin_error, color = "Std Plug-In Error")) +
  geom_point(aes(x = n, y = std_dr_error, color = "Std DR Error")) +
  geom_point(aes(x = n, y = std_plugin_error, color = "Std Plug-In Error")) +
  theme_minimal() +
  labs(
    title = 'Standard Error of Estimators per N',
    x = 'Sample Size (N)',
    y = 'Standard Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/std_error_by_estimator_and_n.png', width = 8, height = 5)

# 3. Average Error of the Estimators by D
simulations %>% 
  group_by(d) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = d, y = avg_dr_error, color = "Avg DR Error")) +
  geom_line(aes(x = d, y = avg_plugin_error, color = "Avg Plug-In Error")) +
  geom_point(aes(x = d, y = avg_dr_error, color = "Avg DR Error")) +
  geom_point(aes(x = d, y = avg_plugin_error, color = "Avg Plug-In Error")) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators by Dimensions (D)',
    x = 'Dimensions (D)',
    y = 'Average Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_estimator_and_d.png', width = 8, height = 5)

# 4. Standard Error of the Estimators by D
simulations %>% 
  group_by(d) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = d, y = std_dr_error, color = "Std DR Error")) +
  geom_line(aes(x = d, y = std_plugin_error, color = "Std Plug-In Error")) +
  geom_point(aes(x = d, y = std_dr_error, color = "Std DR Error")) +
  geom_point(aes(x = d, y = std_plugin_error, color = "Std Plug-In Error")) +
  theme_minimal() +
  labs(
    title = 'Standard Error of Estimators by Dimensions (D)',
    x = 'Dimensions (D)',
    y = 'Standard Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/std_error_by_estimator_and_d.png', width = 8, height = 5)

# 5. Average Error of the Method by N in DR Estimator
simulations %>% 
  group_by(n, method) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = n, y = avg_dr_error, color = method)) +
  geom_point(aes(x = n, y = avg_dr_error, color = method)) +
  theme_minimal() +
  labs(
    title = 'Average Error of Methods by N in DR Estimator',
    x = 'Sample Size (N)',
    y = 'Average DR Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_method_and_n_dr_estimator.png', width = 8, height = 5)

# 6. Average Error of the Estimators by D in DR Estimator
simulations %>% 
  group_by(d, method) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = d, y = avg_dr_error, color = method)) +
  geom_point(aes(x = d, y = avg_dr_error, color = method)) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators by D in DR Estimator',
    x = 'Dimensions (D)',
    y = 'Average DR Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_method_and_d_dr_estimator.png', width = 8, height = 5)

# 7. Average Error of the Estimators by D in Plug-In Estimator
simulations %>% 
  group_by(d, method) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = d, y = avg_plugin_error, color = method)) +
  geom_point(aes(x = d, y = avg_plugin_error, color = method)) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators by D in Plug-In Estimator',
    x = 'Dimensions (D)',
    y = 'Average Plug-In Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_method_and_d_plug_in_estimator.png', width = 8, height = 5)

# 8. Average Error of the Estimators by N in Plug-In Estimator
simulations %>% 
  group_by(n, method) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    std_dr_error = sd(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = n, y = avg_plugin_error, color = method)) +
  geom_point(aes(x = n, y = avg_plugin_error, color = method)) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators by N in Plug-In Estimator',
    x = 'Sample Size (N)',
    y = 'Average Plug-In Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_method_and_n_plug_in_estimator.png', width = 8, height = 5)

# 9. Average Error of Estimators by Scenario (s)
simulations %>% 
  group_by(s, method) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE), 
    avg_plugin_error = mean(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = s, y = avg_dr_error, color = method)) +
  geom_point(aes(x = s, y = avg_dr_error, color = method)) +
  geom_line(aes(x = s, y = avg_plugin_error, color = method), linetype = "dashed") +
  geom_point(aes(x = s, y = avg_plugin_error, color = method), shape = 17) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators by Scenario (s)',
    x = 'Scenario (s)',
    y = 'Average Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_by_scenario_and_method.png', width = 8, height = 5)

# 10. Standard Deviation of Estimators by Scenario (s)
simulations %>% 
  group_by(s, method) %>% 
  summarise(
    std_dr_error = sd(dr_error, na.rm = TRUE), 
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = s, y = std_dr_error, color = method)) +
  geom_point(aes(x = s, y = std_dr_error, color = method)) +
  geom_line(aes(x = s, y = std_plugin_error, color = method), linetype = "dashed") +
  geom_point(aes(x = s, y = std_plugin_error, color = method), shape = 17) +
  theme_minimal() +
  labs(
    title = 'Standard Deviation of Estimators by Scenario (s)',
    x = 'Scenario (s)',
    y = 'Standard Deviation',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/std_dev_by_scenario_and_method.png', width = 8, height = 5)

# 11. Average Error vs. n_div_d
simulations %>% 
  group_by(n_div_d, method) %>% 
  summarise(
    avg_dr_error = mean(dr_error, na.rm = TRUE),
    avg_plugin_error = mean(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = n_div_d, y = avg_dr_error, color = method)) +
  geom_point(aes(x = n_div_d, y = avg_dr_error, color = method)) +
  geom_line(aes(x = n_div_d, y = avg_plugin_error, color = method), linetype = "dashed") +
  geom_point(aes(x = n_div_d, y = avg_plugin_error, color = method), shape = 17) +
  theme_minimal() +
  labs(
    title = 'Average Error of Estimators vs. n_div_d',
    x = 'n / d',
    y = 'Average Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/avg_error_vs_n_div_d.png', width = 8, height = 5)

# 12. Standard Deviation vs. n_div_d
simulations %>% 
  group_by(n_div_d, method) %>% 
  summarise(
    std_dr_error = sd(dr_error, na.rm = TRUE),
    std_plugin_error = sd(plugin_error, na.rm = TRUE)
  ) %>% 
  ggplot() +
  geom_line(aes(x = n_div_d, y = std_dr_error, color = method)) +
  geom_point(aes(x = n_div_d, y = std_dr_error, color = method)) +
  geom_line(aes(x = n_div_d, y = std_plugin_error, color = method), linetype = "dashed") +
  geom_point(aes(x = n_div_d, y = std_plugin_error, color = method), shape = 17) +
  theme_minimal() +
  labs(
    title = 'Standard Deviation of Estimators vs. n_div_d',
    x = 'n / d',
    y = 'Standard Deviation',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/std_dev_vs_n_div_d.png', width = 8, height = 5)

# 13. Relationship Between plugin_sd and plugin_error
simulations %>% 
  ggplot(aes(x = plugin_sd, y = plugin_error, color = method)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = 'Relationship Between Plugin SD and Plugin Error',
    x = 'Plugin Standard Deviation',
    y = 'Plugin Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/plugin_sd_vs_plugin_error.png', width = 8, height = 5)

# 14. Relationship Between dr_sd and dr_error
simulations %>% 
  ggplot(aes(x = dr_sd, y = dr_error, color = method)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = 'Relationship Between DR SD and DR Error',
    x = 'DR Standard Deviation',
    y = 'DR Error',
    color = 'Method'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2b-plots/dr_sd_vs_dr_error.png', width = 8, height = 5)