# Load required packages
library(tidyverse)

# Read in the data
q3 <- read_csv("data_for_HW4.csv")

# Separate the data into the two sources
q3_e1 <- q3 %>% filter(e == 1)
q3_e2 <- q3 %>% filter(e == 2)

# Compute the Average Treatment Effect (ATE) in each data source by simple difference in means
# For e=1
mean_y_t1_e1 <- mean(q3_e1$y[q3_e1$t == 1])
mean_y_t0_e1 <- mean(q3_e1$y[q3_e1$t == 0])
ate_e1 <- mean_y_t1_e1 - mean_y_t0_e1

# For e=2
mean_y_t1_e2 <- mean(q3_e2$y[q3_e2$t == 1])
mean_y_t0_e2 <- mean(q3_e2$y[q3_e2$t == 0])
ate_e2 <- mean_y_t1_e2 - mean_y_t0_e2

# Get confidence intervals using a two-sample t-test for each source
ttest_e1 <- t.test(q3_e1$y[q3_e1$t == 1], q3_e1$y[q3_e1$t == 0])
ttest_e2 <- t.test(q3_e2$y[q3_e2$t == 1], q3_e2$y[q3_e2$t == 0])

# Print out the results
cat("Results for data source e=1 (possibly targeted):\n")
cat("ATE estimate:", ate_e1, "\n")
cat("95% CI:", ttest_e1$conf.int[1], "to", ttest_e1$conf.int[2], "\n\n")

cat("Results for data source e=2 (fully randomized):\n")
cat("ATE estimate:", ate_e2, "\n")
cat("95% CI:", ttest_e2$conf.int[1], "to", ttest_e2$conf.int[2], "\n\n")

# Comment on findings:
# The first data source (e=1) may show a different ATE than the fully randomized second data source (e=2),
# potentially due to the non-random treatment assignment. If the treatment in e=1 was targeted based on certain 
# characteristics, that may cause bias in the simple difference-in-means estimate. 
# In contrast, e=2 is fully randomized, so the difference in means should yield an unbiased estimate of the ATE.