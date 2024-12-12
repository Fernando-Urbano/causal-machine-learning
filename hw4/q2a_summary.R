# check all the files inside the folder q2a-data
library(tidyverse)
library(knitr)
library(kableExtra)
# Read all the files in the folder
files <- list.files("q2a-data", full.names = TRUE)

simulations <- tibble()

for (file in files) {
  simulations <- simulations %>% 
    rbind(read_csv(file))
}

simulations <- simulations %>%
  mutate(
    error = abs(estimate - 1),
    n_div_d = n / d
  )

simulations_summary <- simulations %>% 
  group_by(n, d, sBeta, sTheta, estimator, n_div_d) %>% 
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  purrr::set_names("N", "D", "sBeta", "sTheta", "Estimator", "N / D", "Avg Error", "Std Error")

n <- 4

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
  output_file <- "report/q2a_simulations_summary_part" %>% paste0(i) %>% paste0(".tex")
  
  # Write the LaTeX table to the file
  writeLines(latex_table, con = output_file)
  
  # Optional: Print a message upon successful saving
  cat("LaTeX table saved to", output_file, "\n")

}

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  group_by(estimator, n_div_d) %>% 
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n_div_d, avg_error)) +
  geom_line(aes(color = estimator)) +
  geom_point(aes(color = estimator)) +
  theme_minimal() +
  labs(
    title = 'Avg Error of the Estimators per N / D',
    x = 'N Sample / Dimensions',
    y = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/avg_error_of_the_estimator_per_n_div_d.png', width = 8, height = 5)


simulations %>% 
  
  

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  group_by(estimator, n_div_d) %>% 
  filter(n_div_d >= 2) %>% 
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n_div_d, avg_error)) +
  geom_line(aes(color = estimator)) +
  geom_point(aes(color = estimator)) +
  theme_minimal() +
  labs(
    title = 'Avg Error of the Estimators per N / D (filtering Ratio >= 2)',
    x = 'N Sample / Dimensions',
    y = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/avg_error_of_the_estimator_per_n_div_d_filtering_2.png', width = 8, height = 5)

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  group_by(estimator, n_div_d) %>% 
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n_div_d, std_error)) +
  geom_line(aes(color = estimator)) +
  geom_point(aes(color = estimator)) +
  theme_minimal() +
  labs(
    title = 'Std Error of the Estimators per N / D',
    x = 'N Sample / Dimensions',
    y = 'Std Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/std_error_of_the_estimator_per_n_div_d.png', width = 8, height = 5)


simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  group_by(estimator, n_div_d) %>% 
  filter(n_div_d >= 2) %>% 
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n_div_d, std_error)) +
  geom_line(aes(color = estimator)) +
  geom_point(aes(color = estimator)) +
  theme_minimal() +
  labs(
    title = 'Std Error of the Estimators per N / D (filtering Ratio >= 2)',
    x = 'N Sample / Dimensions',
    y = 'Std Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/std_error_of_the_estimator_per_n_div_d_filtering_2.png', width = 8, height = 5)

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  group_by(estimator, n) %>% 
  mutate(n = as.character(n)) %>%
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n, avg_error)) +
  geom_bar(aes(fill = estimator), stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Avg Error of the Estimators per N',
    x = 'N',
    y = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/avg_error_of_the_estimator_per_n.png', width = 8, height = 5)

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  filter(n_div_d >= 2) %>% 
  group_by(estimator, n) %>% 
  mutate(n = as.character(n)) %>% 
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n, avg_error)) +
  geom_bar(aes(fill = estimator), stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Avg Error of the Estimators per N (filtering N / D >= 2)',
    x = 'N',
    y = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/avg_error_of_the_estimator_per_n_filtering_2.png', width = 8, height = 5)

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  group_by(estimator, n) %>% 
  mutate(n = as.character(n)) %>%
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n, std_error)) +
  geom_bar(aes(fill = estimator), stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Std Error of the Estimators per N',
    x = 'N',
    y = 'Std Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/avg_error_of_the_estimator_per_n.png', width = 8, height = 5)

simulations %>% 
  mutate(error = abs(estimate - 1)) %>% 
  mutate(n_div_d = n / d) %>% 
  filter(n_div_d >= 2) %>% 
  group_by(estimator, n) %>% 
  mutate(n = as.character(n)) %>%
  summarise(avg_error = mean(error, na.rm = TRUE), std_error = sd(error, na.rm = TRUE)) %>% 
  ggplot(aes(n, std_error)) +
  geom_bar(aes(fill = estimator), stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Std Error of the Estimators per N (filtering N / D >= 2)',
    x = 'N',
    y = 'Std Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/std_error_of_the_estimator_per_n.png', width = 8, height = 5)


# -----------------------------------------------------------------------------
# Plot 1a: Avg Error of Estimators by Penalty Type
# -----------------------------------------------------------------------------
avg_error_penalty <- simulations %>% 
  group_by(estimator, penalty) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE),
    std_error = sd(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(avg_error_penalty, aes(penalty, avg_error, fill = estimator)) +
  geom_bar(stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Avg Error of the Estimators by Penalty Type',
    x = 'Penalty Type',
    y = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')
ggsave('q2a-plots/avg_error_of_estimators_by_penalty.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 1b: Avg Error of Estimators by Penalty Type (n_div_d >= 2)
# -----------------------------------------------------------------------------
avg_error_penalty_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, penalty) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE),
    std_error = sd(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(avg_error_penalty_filtered, aes(penalty, avg_error, fill = estimator)) +
  geom_bar(stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Avg Error of the Estimators by Penalty Type (n/d >= 2)',
    x = 'Penalty Type',
    y = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/avg_error_of_estimators_by_penalty_filter2.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 2a: Std Error of Estimators by Penalty Type
# -----------------------------------------------------------------------------
std_error_penalty <- simulations %>% 
  group_by(estimator, penalty) %>% 
  summarise(
    std_error = sd(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(std_error_penalty, aes(penalty, std_error, fill = estimator)) +
  geom_bar(stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Std Error of the Estimators by Penalty Type',
    x = 'Penalty Type',
    y = 'Std Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/std_error_of_estimators_by_penalty.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 2b: Std Error of Estimators by Penalty Type (n_div_d >= 2)
# -----------------------------------------------------------------------------
std_error_penalty_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, penalty) %>% 
  summarise(
    std_error = sd(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(std_error_penalty_filtered, aes(penalty, std_error, fill = estimator)) +
  geom_bar(stat = 'identity', position = "dodge") +
  theme_minimal() +
  labs(
    title = 'Std Error of the Estimators by Penalty Type (n/d >= 2)',
    x = 'Penalty Type',
    y = 'Std Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/std_error_of_estimators_by_penalty_filter2.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 3a: Avg Error by sBeta and sTheta
# -----------------------------------------------------------------------------
avg_error_sBeta_sTheta <- simulations %>% 
  group_by(estimator, sBeta, sTheta) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(avg_error_sBeta_sTheta, aes(sBeta, sTheta, fill = avg_error)) +
  geom_tile() +
  facet_wrap(~ estimator) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = 'Avg Error of Estimators by sBeta and sTheta',
    x = 'sBeta',
    y = 'sTheta',
    fill = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/avg_error_by_sBeta_sTheta.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 3b: Avg Error by sBeta and sTheta (n_div_d >= 2)
# -----------------------------------------------------------------------------
avg_error_sBeta_sTheta_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, sBeta, sTheta) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(avg_error_sBeta_sTheta_filtered, aes(sBeta, sTheta, fill = avg_error)) +
  geom_tile() +
  facet_wrap(~ estimator) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = 'Avg Error of Estimators by sBeta and sTheta (n/d >= 2)',
    x = 'sBeta',
    y = 'sTheta',
    fill = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/avg_error_by_sBeta_sTheta_filter2.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 4a: Performance of Estimators Across Different n/d Ratios and Penalty Types
# -----------------------------------------------------------------------------
avg_error_n_div_d_penalty <- simulations %>% 
  group_by(estimator, penalty, n_div_d) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(avg_error_n_div_d_penalty, aes(n_div_d, avg_error, color = estimator)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ penalty) +
  theme_minimal() +
  labs(
    title = 'Avg Error Across n/d Ratios by Penalty Type',
    x = 'N / D Ratio',
    y = 'Avg Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/avg_error_across_n_div_d_by_penalty.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 4b: Performance of Estimators Across Different n/d Ratios and Penalty Types (n_div_d >= 2)
# -----------------------------------------------------------------------------
avg_error_n_div_d_penalty_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, penalty, n_div_d) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(avg_error_n_div_d_penalty_filtered, aes(n_div_d, avg_error, color = estimator)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ penalty) +
  theme_minimal() +
  labs(
    title = 'Avg Error Across n/d Ratios by Penalty Type (n/d >= 2)',
    x = 'N / D Ratio',
    y = 'Avg Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/avg_error_across_n_div_d_by_penalty_filter2.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 5a: Boxplot of Errors by Estimator
# -----------------------------------------------------------------------------
ggplot(simulations, aes(estimator, error, fill = estimator)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = 'Distribution of Errors by Estimator',
    x = 'Estimator',
    y = 'Error'
  ) +
  theme(legend.position = 'none')

ggsave('q2a-plots/error_distribution_by_estimator.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 5b: Boxplot of Errors by Estimator (n_div_d >= 2)
# -----------------------------------------------------------------------------
ggplot(simulations %>% filter(n_div_d >= 2), aes(estimator, error, fill = estimator)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = 'Distribution of Errors by Estimator (n/d >= 2)',
    x = 'Estimator',
    y = 'Error'
  ) +
  theme(legend.position = 'none')

ggsave('q2a-plots/error_distribution_by_estimator_filter2.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 6a: Interaction Effect Between Estimator and Penalty on Avg Error
# -----------------------------------------------------------------------------
interaction_avg_error <- simulations %>% 
  group_by(estimator, penalty) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(interaction_avg_error, aes(estimator, avg_error, fill = penalty)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme_minimal() +
  labs(
    title = 'Interaction Effect of Estimator and Penalty on Avg Error',
    x = 'Estimator',
    y = 'Avg Error',
    fill = 'Penalty Type'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/interaction_estimator_penalty_avg_error.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 6b: Interaction Effect Between Estimator and Penalty on Avg Error (n_div_d >= 2)
# -----------------------------------------------------------------------------
interaction_avg_error_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, penalty) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(interaction_avg_error_filtered, aes(estimator, avg_error, fill = penalty)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme_minimal() +
  labs(
    title = 'Interaction Effect of Estimator and Penalty on Avg Error (n/d >= 2)',
    x = 'Estimator',
    y = 'Avg Error',
    fill = 'Penalty Type'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/interaction_estimator_penalty_avg_error_filter2.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 7a: Trend of Avg Error Over Replicates
# -----------------------------------------------------------------------------
trend_avg_error_replicates <- simulations %>% 
  group_by(estimator, replicate) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(trend_avg_error_replicates, aes(replicate, avg_error, color = estimator)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = 'Trend of Avg Error Over Replicates',
    x = 'Replicate',
    y = 'Avg Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/trend_avg_error_over_replicates.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 7b: Trend of Avg Error Over Replicates (n_div_d >= 2)
# -----------------------------------------------------------------------------
trend_avg_error_replicates_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, replicate) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(trend_avg_error_replicates_filtered, aes(replicate, avg_error, color = estimator)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = 'Trend of Avg Error Over Replicates (n/d >= 2)',
    x = 'Replicate',
    y = 'Avg Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/trend_avg_error_over_replicates_filter2.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 8a: Heatmap of Avg Error by Estimator and n/d Ratio
# -----------------------------------------------------------------------------
heatmap_avg_error <- simulations %>% 
  group_by(estimator, n_div_d) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(heatmap_avg_error, aes(factor(n_div_d), estimator, fill = avg_error)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = 'Heatmap of Avg Error by Estimator and N/D Ratio',
    x = 'N / D Ratio',
    y = 'Estimator',
    fill = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/heatmap_avg_error_estimator_n_div_d.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 8b: Heatmap of Avg Error by Estimator and n/d Ratio (n_div_d >= 2)
# -----------------------------------------------------------------------------
heatmap_avg_error_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(estimator, n_div_d) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(heatmap_avg_error_filtered, aes(factor(n_div_d), estimator, fill = avg_error)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = 'Heatmap of Avg Error by Estimator and N/D Ratio (n/d >= 2)',
    x = 'N / D Ratio',
    y = 'Estimator',
    fill = 'Avg Error'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/heatmap_avg_error_estimator_n_div_d_filter2.png', width = 10, height = 6)

# -----------------------------------------------------------------------------
# Plot 9a: Faceted Line Plots for Each Penalty Type Showing Avg Error Across n/d Ratios
# -----------------------------------------------------------------------------
faceted_avg_error <- simulations %>% 
  group_by(penalty, estimator, n_div_d) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(faceted_avg_error, aes(n_div_d, avg_error, color = estimator)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ penalty) +
  theme_minimal() +
  labs(
    title = 'Avg Error Across N/D Ratios by Penalty Type and Estimator',
    x = 'N / D Ratio',
    y = 'Avg Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/faceted_avg_error_n_div_d_penalty.png', width = 12, height = 8)

# -----------------------------------------------------------------------------
# Plot 9b: Faceted Line Plots for Each Penalty Type Showing Avg Error Across n/d Ratios (n_div_d >= 2)
# -----------------------------------------------------------------------------
faceted_avg_error_filtered <- simulations %>% 
  filter(n_div_d >= 2) %>%
  group_by(penalty, estimator, n_div_d) %>% 
  summarise(
    avg_error = mean(error, na.rm = TRUE)
  ) %>% 
  ungroup()

ggplot(faceted_avg_error_filtered, aes(n_div_d, avg_error, color = estimator)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ penalty) +
  theme_minimal() +
  labs(
    title = 'Avg Error Across N/D Ratios by Penalty Type and Estimator (n/d >= 2)',
    x = 'N / D Ratio',
    y = 'Avg Error',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/faceted_avg_error_n_div_d_penalty_filter2.png', width = 12, height = 8)

# -----------------------------------------------------------------------------
# Plot 10a: Scatter Plot of N vs. D Colored by Estimator
# -----------------------------------------------------------------------------
ggplot(simulations, aes(n, d, color = estimator)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = 'Scatter Plot of N vs. D Colored by Estimator',
    x = 'Sample Size (N)',
    y = 'Dimensions (D)',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/scatter_n_vs_d_by_estimator.png', width = 8, height = 5)

# -----------------------------------------------------------------------------
# Plot 10b: Scatter Plot of N vs. D Colored by Estimator (n_div_d >= 2)
# -----------------------------------------------------------------------------
ggplot(simulations %>% filter(n_div_d >= 2), aes(n, d, color = estimator)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = 'Scatter Plot of N vs. D Colored by Estimator (n/d >= 2)',
    x = 'Sample Size (N)',
    y = 'Dimensions (D)',
    color = 'Estimator'
  ) +
  theme(legend.position = 'bottom')

ggsave('q2a-plots/scatter_n_vs_d_by_estimator_filter2.png', width = 8, height = 5)

