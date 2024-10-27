library(tidyverse)
library(readxl)

SUBTRACT_CITY = 0

SMALL_CITY = 1 - SUBTRACT_CITY
BIG_CITY = 2 - SUBTRACT_CITY

marketing_data <- read.csv("marketing_data.csv", sep=",") %>% 
  as_tibble() %>% 
  select(-X) %>% 
  mutate(city = city - SUBTRACT_CITY)

# Statistics of Data
# Mean of the Variables
marketing_data %>% 
  summarise(across(c(visitors, spend, city), mean, na.rm = TRUE))

# Mean Spending and Visitors per Group
marketing_data %>% 
  group_by(city) %>% 
  summarise(across(c(visitors, spend), mean, na.rm = TRUE))

# Covariance Matrix
marketing_data %>% 
  cov()

# Cov city * spending, spending
cov(marketing_data$city * marketing_data$spend, marketing_data$spend)

# Cov city, spending
cov(marketing_data$city, marketing_data$spend)

# Var spending
var(marketing_data$spend)

# Var city
var(marketing_data$city)

# Var visitors
var(marketing_data$visitors)

# Cov city * demeaned spending, spending
cov((marketing_data$spend - mean(marketing_data$spend)) * marketing_data$city, marketing_data$spend)

# General Regression and Plot
ols_visitors_spending <- lm(visitors ~ spend, data = marketing_data)
summary(ols_visitors_spending)

# Call:
# lm(formula = visitors ~ spend, data = marketing_data)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -5.1353 -2.1523  0.3405  2.1699  5.2667 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  12.0359     0.3813   31.57   <2e-16 ***
# spend        -1.3877     0.6768   -2.05   0.0416 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.477 on 198 degrees of freedom
# Multiple R-squared:  0.02079,	Adjusted R-squared:  0.01585 
# F-statistic: 4.205 on 1 and 198 DF,  p-value: 0.04163

# By City Regression and Plot
marketing_data %>% 
  ggplot(aes(x=spend, y=visitors)) +
  geom_point() +
  labs(
    title="Spend vs Visitors",
    x="Spend",
    y="Visitors"
  ) +
  geom_smooth(aes(color = "Linear Regression"), method = "lm", se = FALSE) +
  theme_minimal() +
  scale_color_manual(name = "", values = c("Linear Regression" = "blue")) +
  theme(legend.position = "bottom")
ggsave("marketing_spending_vs_visitors.png", width = 8, height = 5, dpi = 300)

marketing_data_small_cities <- marketing_data %>% 
  filter(city == SMALL_CITY)

marketing_data_big_cities <- marketing_data %>% 
  filter(city == BIG_CITY)

marketing_data %>% 
  mutate(city=ifelse(city == SMALL_CITY, "Small", "Big")) %>% 
  ggplot(aes(x=spend, y=visitors, color=city)) +
  geom_point() +
  labs(
    title="Spend vs Visitors",
    x="Spend",
    y="Visitors"
  ) +
  geom_smooth(aes(color = city), method = "lm", se = FALSE) +
  theme_minimal() +
  scale_color_manual(name = "", values = c("Big" = "blue", "Small" = "red")) +
  theme(legend.position = "bottom")
ggsave("marketing_spending_vs_visitors_by_city.png", width = 8, height = 5, dpi = 300)

ols_visitors_spending_small_cities <- lm(
  visitors ~ spend, data = marketing_data_small_cities
)
summary(ols_visitors_spending_small_cities)

# Call:
#   lm(formula = visitors ~ spend, data = marketing_data_small_cities)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.2205 -0.6139 -0.0231  0.8664  2.3934 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  12.2342     0.2196  55.704  < 2e-16 ***
#   spend         2.6232     0.4873   5.383 4.37e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.151 on 107 degrees of freedom
# Multiple R-squared:  0.2131,	Adjusted R-squared:  0.2057 
# F-statistic: 28.97 on 1 and 107 DF,  p-value: 4.37e-07


ols_visitors_spending_big_cities <- lm(
  visitors ~ spend, data = marketing_data_big_cities
)
summary(ols_visitors_spending_big_cities)

# Call:
#   lm(formula = visitors ~ spend, data = marketing_data_big_cities)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.84717 -0.85156 -0.05242  0.85707  2.60109 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.7741     0.3517  19.259  < 2e-16 ***
#   spend         3.5916     0.5219   6.882 7.95e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.153 on 89 degrees of freedom
# Multiple R-squared:  0.3473,	Adjusted R-squared:   0.34 
# F-statistic: 47.37 on 1 and 89 DF,  p-value: 7.952e-10

# Multivariate Regression
marketing_data_with_interaction <- marketing_data %>%
  mutate(interaction_spend_city = (spend) * city)

ols_visitors_with_interaction <- lm(
  visitors ~ ., data = marketing_data_with_interaction
)
summary(ols_visitors_with_interaction)

# Call:
#   lm(formula = visitors ~ ., data = marketing_data_with_interaction)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.2205 -0.7589 -0.0294  0.8690  2.6011 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             17.6943     0.5628  31.439   <2e-16 ***
#   spend                    1.6549     1.1062   1.496    0.136    
# city                    -5.4601     0.4144 -13.176   <2e-16 ***
#   interaction_spend_city   0.9683     0.7139   1.356    0.177    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.152 on 196 degrees of freedom
# Multiple R-squared:  0.7903,	Adjusted R-squared:  0.7871 
# F-statistic: 246.2 on 3 and 196 DF,  p-value: < 2.2e-16