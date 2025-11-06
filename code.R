library(tidyverse)
library(readxl)
library(fitdistrplus)
library(viridis)
library(dplyr)
library(quantreg)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")
data$sex <- as.factor(data$Sex)
data$beta60 <- data$`Beta60 (g/kg/h)`
data$weight <- data$`Weight (kg)`
data$height <- data$`Height (cm)`
data$age <- data$`Age (years)`
data$beta <- -data$beta60

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BETA EDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# beta vs weight
weight_plot <- ggplot(data, aes(x=weight, y=beta, colour = sex)) +
  geom_point() +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(title="Figure 1: β vs Weight", x="Weight (kg)", y="β Elimination Rate (g/kg/h)")

# beta vs height
height_plot <- ggplot(data, aes(x=height, y=beta, colour = sex)) +
  geom_point() +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(title="Figure 2: β vs Height", x="Weight (cm)", y="β Elimination Rate (g/kg/h)")

# beta vs age
age_plot <- ggplot(data, aes(x=age, y=beta, colour = sex)) +
  geom_point() +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(title="Figure 3: β vs Age", x="Age (years)", y="β Elimination Rate (g/kg/h)")

# beta vs gender
sex_plot <- ggplot(data, aes(x = sex, y = beta60,
                             fill = sex)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("lightblue", "steelblue")) +
  labs(
    title = "Figure 4: β vs Gender",
    x = "Gender",
    y = "β Elimination Rate (g/kg/h)"
  ) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# pairwise correlations between characteristics and beta
corr <- data %>%
  dplyr::select(beta, weight, age, height) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(3)

#%%%%%%%%%%%%%%%%%%%% BETA DISTRIBUTION TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Density plot of beta.
density_plot <- ggplot(data, aes(x = beta)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "navy",
                 alpha = 0.55) +
  geom_density(color = "navy", 
               size = 1) +
  geom_vline(aes(xintercept=mean(beta)),
             linetype = "dashed",
             color = "navy",
             size = 1) +
  geom_vline(aes(xintercept=quantile(beta, 0.025)),
             linetype = "dotted",
             color = "navy",
             size = 1) +
  geom_vline(aes(xintercept=quantile(beta, 0.975)),
             linetype = "dotted",
             color = "navy",
             size = 1) +
  labs(x = expression(""*beta), y = "Density", title = expression("Distribution of "*beta))

# compute group statistics by sex
stats <- data %>%
  group_by(sex) %>%
  summarize(
    mean = mean(beta),
    q025 = quantile(beta, 0.025),
    q975 = quantile(beta, 0.975)
  )

# beta density plot (by sex)
sex_density_plot <- ggplot(data, aes(x = beta, color = sex, fill = sex)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005, alpha = 0.4, position = "identity") +
  geom_density(size = 1, alpha = 0) +
  geom_vline(data = stats, aes(xintercept = mean, color = sex),
             linetype = "dashed", size = 1) +
  geom_vline(data = stats, aes(xintercept = q025, color = sex),
             linetype = "dotted", size = 1) +
  geom_vline(data = stats, aes(xintercept = q975, color = sex),
             linetype = "dotted", size = 1) +
  scale_color_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  scale_fill_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(x = expression(""*beta), y = "Density", title = expression("Distribution of "*beta*" by Sex"))

# Residual plots of fitted distributions
normal_fit <- fitdist(data$beta, distr = "norm", method = "mle")
beta_fit <- fitdist(data$beta, distr = "beta", method = "mle")
gamma_fit <- fitdist(data$beta, distr = "gamma", method = "mle")

# plot each fitted distribution on residual plot
plot.legend <- c("Normal", "Beta", "Gamma")
denscomp(list(normal_fit, beta_fit, gamma_fit), legendtext = plot.legend)
qqcomp(list(normal_fit, beta_fit, gamma_fit), legendtext = plot.legend)
cdfcomp(list(normal_fit, beta_fit, gamma_fit), legendtext = plot.legend)
ppcomp(list(normal_fit, beta_fit, gamma_fit), legendtext = plot.legend)

# modelling as gamma and finding 2.5% quantile
gamma_fit <- fitdist(-data$beta60, distr = "gamma", method = "mle")
summary(gamma_fit)
par(mar=c(1, 1, 1, 1))
plot(gamma_fit)
gamma_shape <- gamma_fit$estimate[[1]]
gamma_rate <- gamma_fit$estimate[[2]]
q025 <- -qgamma(0.975, gamma_shape, gamma_rate)
q025

# modelling as gamma by sex
m_data <- data %>% filter(Sex == "male")
m_gamma_fit <- fitdist(-m_data$beta60, distr = "gamma", method = "mle")
m_gamma_shape <- m_gamma_fit$estimate[[1]]
m_gamma_rate <- m_gamma_fit$estimate[[2]]
m_q025 <- -qgamma(0.975, m_gamma_shape, m_gamma_rate)

f_data <- data %>% filter(Sex == "female")
f_gamma_fit <- fitdist(-f_data$beta60, distr = "gamma", method = "mle")
f_gamma_shape <- f_gamma_fit$estimate[[1]]
f_gamma_rate <- f_gamma_fit$estimate[[2]]
f_q025 <- -qgamma(0.975, f_gamma_shape, f_gamma_rate)

# prepare data frame by sex for plotting.
gamma_sex_df <- data.frame(
  x = data$beta,
  male_density   = dgamma(data$beta, shape = m_gamma_shape, rate = m_gamma_rate),
  female_density = dgamma(data$beta, shape = f_gamma_shape, rate = f_gamma_rate)
)

# convert dataframe to long format.
gamma_long <- gamma_sex_df %>%
  pivot_longer(cols = c(male_density, female_density),
               names_to = "Sex",
               values_to = "density") %>%
  mutate(Sex = ifelse(Sex == "male_density", "male", "female"))

# extract gamma quantile statistics by sex.
stats_gamma <- data.frame(
  Sex = rep(c("male", "female"), each = 3),
  Quantile = rep(c("q025", "q50", "q975"), times = 2),
  value = c(m_quantiles, f_quantiles)
)

# plot the gamma distribution by sex.
sex_gamma_plot <- ggplot(data, aes(x = beta, color = Sex, fill = Sex)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005, alpha = 0.4, position = "identity") +
  geom_line(data = gamma_long,
            aes(x = x, y = density, color = Sex),
            size = 1.2) +
  geom_vline(
    data = stats_gamma,
    aes(xintercept = value, color = Sex, linetype = Quantile),
    size = 1
  ) +
  scale_color_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  scale_fill_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  scale_linetype_manual(
    values = c("q025" = "dotted", "q50" = "dashed", "q975" = "dotted"),
    labels = c("2.5%", "50%", "97.5%")
  ) +
  labs(
    x = expression(beta),
    y = "Density",
    title = expression("Gamma Distributions of "*beta*" by Sex"),
    linetype = "Quantile"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold")
  )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODELLING BETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- linear modelling beta 60 attempt ---

# model incorporating all characteristics
lm_model <- lm(beta ~ weight + age + height + sex, data = data)

# Residual plots
par(mfrow = c(2,2))
plot(lm_model)

# Try modelling square root of beta (residual plots indicate
# a slight quadratic mean-variance relationship)

sqrt_lm_model <- lm(sqrt(beta) ~ weight + age + height + sex, data = data)

# Residual plots
par(mfrow = c(2,2))
plot(sqrt_lm_model)

# Create predicted values
pred_data <- data

# quantiles for lm
resid_q025 <- quantile(residuals(lm_model), probs = 0.025)
pred_data$lm_pred <- predict(lm_model)
pred_data$lm_q025 <- pred_data$lm_pred + resid_q025

# using quantile regression instead
rq_model <- rq(beta ~ weight + age + height + sex, data = data, tau = 0.025)
pred_data$rq_pred <- predict(rq_model)

# Combine into a long format for ggplot
plot_data <- pred_data %>%
  tidyr::pivot_longer(cols = c(lm_q025, rq_pred),
                      names_to = "model",
                      values_to = "predicted_beta")

# Plot
quantile_v_lm_plot <- ggplot(plot_data, aes(x = beta, y = predicted_beta, color = model)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = expression("Comparing Models for the 2.5% Quantile for "*beta),
    x = expression("Actual "*beta),
    y = expression("Predicted "*beta),
    color = "Model Type"
  ) +
  scale_color_manual(values = c("lm_q025" = "firebrick", "rq_pred" = "forestgreen"),
                     labels = c("Linear Regression", "Quantile Regression"))



# compute group statistics
stats <- data %>%
  group_by(Sex) %>%
  summarize(
    mean = mean(beta60),
    q025 = quantile(beta60, 0.025),
    q975 = quantile(beta60, 0.975)
  )

stats

# plot
sex_density_plot <- ggplot(data, aes(x = beta60, color = Sex, fill = Sex)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005, alpha = 0.55, position = "identity") +
  geom_density(size = 1, alpha= 0) +
  geom_vline(data = stats, aes(xintercept = mean, color = Sex),
             linetype = "dashed", size = 1) +
  geom_vline(data = stats, aes(xintercept = q025, color = Sex),
             linetype = "dotted", size = 1) +
  geom_vline(data = stats, aes(xintercept = q975, color = Sex),
             linetype = "dotted", size = 1) +
  theme_minimal() +
  labs(x = "beta_60", y = "Density", title = "Distribution of beta_60 by Sex")

sex_density_plot

# quantile regression to increase accuracy of quantile prediction
library(quantreg)
rqfit <- rq(sqrt((-1)*beta60) ~ weight + age + height + sex, data = data, tau = 0.025)
summary(rqfit)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# setup from question
test_person <- data.frame(weight = 70, height = 160, age = 70, sex = "female")
Ct <- 0.15
t  <- 2

# empirical quantiles for Beta60 (population approach)
beta_pop <- quantile(data$beta, probs = c(0.025, 0.5, 0.975))

# C0 values using population quantiles
C0_pop <- Ct + beta_pop * t

# gamma fitted whole population distribution
gamma_quantiles <- qgamma(c(0.025, 0.5, 0.975), gamma_shape, gamma_rate)
C0_gamma <- Ct + gamma_quantiles * t

# gamma fitted to males
m_gamma_quantiles <- qgamma(c(0.025, 0.5, 0.975), m_gamma_shape, m_gamma_rate)
C0_gamma_m <- Ct + m_gamma_quantiles * t

# gamma fitted to females
f_gamma_quantiles <- qgamma(c(0.025, 0.5, 0.975), f_gamma_shape, f_gamma_rate)
C0_gamma_f <- Ct + f_gamma_quantiles * t

# predicted Beta using regression model
lm_pred_vals <- predict(lm_model, newdata = test_person, interval = "prediction", level = 0.95)

# C0 values from regression model
beta_pred_lm <- as.numeric(lm_pred_vals)
names(beta_pred_lm) <- c("Fitted", "Lower (PI 2.5%)", "Upper (PI 97.5%)")
C0_lm <- Ct + beta_pred_lm * t

# predicted beta using quantile regression model
rq_pred_vals <- predict(rq_model, newdata = test_person, interval = "confidence", level = 0.95)
C0_rq <- Ct + rq_pred_vals * t

# combine results into one comparison table 
results_table <- tibble(
  Approach = c("Population (empirical quantiles)", "Gamma Fitted", "Male Gamma Fitted", "Female Gamma Fitted", "Linear Regression", "Quantile Regression"),
  Lower_2.5 = c(C0_pop[1], C0_gamma[1], C0_gamma_m[1], C0_gamma_f[1], C0_lm[2], C0_rq[2]),  # note ordering: lwr=2nd col in pred
  Central   = c(C0_pop[2], C0_gamma[2], C0_gamma_m[2], C0_gamma_f[2], C0_lm[1], C0_rq[1]),
  Upper_97.5 = c(C0_pop[3], C0_gamma[3], C0_gamma_m[3], C0_gamma_f[3], C0_lm[3], C0_rq[3])
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# define A, Co, Vd
data$A <- data$`Amount of Alcohol Consumed (g)`
data$Co <- data$`Co (g/Kg)`
data$Vd <- data$A/(data$Co *data$weight)

# view summary and quantiles of Vd
summary(data$Vd)
quantile(data$Vd, probs = c(0.025, 0.5, 0.975))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vd EDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# histogram of Vd
ggplot(data, aes(x = Vd)) +
  geom_histogram(bins = 20, colour = "navy", fill = "steelblue") +
  labs(title = "Distribution of Volume of Distribution (Vd)",
       x = "Vd (L/kg)", y = "Count")

# Vd vs weight
weight_plot2 <- ggplot(data, aes(x = weight, y = Vd, colour = sex)) +
  geom_point() +
  geom_smooth(method = "lm", colour = "navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(title = "Figure 5: Vd vs Weight", x = "Weight (kg)", y = "Vd (L/kg)")

# Vd vs height
height_plot2 <- ggplot(data, aes(x = height, y = Vd, colour = sex)) +
  geom_point() +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  geom_smooth(method = "lm", colour = "navy") +
  labs(title = "Figure 6: Vd vs Height", x = "Height (cm)", y = "Vd (L/kg)")

# Vd vs age
age_plot2 <- ggplot(data, aes(x = age, y = Vd, colour = sex)) +
  geom_point() +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  geom_smooth(method = "lm", colour = "navy") +
  labs(title = "Figure 7: Vd vs Age", x = "Age (years)", y = "Vd (L/kg)")

# Vd vs gender
sex_plot2 <- ggplot(data, aes(x = sex, y = Vd,
                              fill = sex)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("lightblue", "steelblue")) +
  labs(
    title = "Figure 8: Vd vs Gender",
    x = "Gender",
    y = "Vd (L/kg)"
  ) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

(weight_plot2 + height_plot2) /
  (age_plot2 + sex_plot2) /
  plot_layout(ncol = 1) +
  plot_annotation(
    title = "Vd Characteristics Plots",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Test correlation to investigate independent assumption
cor.test(data$beta, data$Vd, use = "complete.obs")

# View quantiles of each coefficient for comparison
beta_range <- quantile(data$beta, probs = c(0.025, 0.975))
Vd_range   <- quantile(data$Vd, probs = c(0.025, 0.975))

# Scatter plot between beta and V_d
ggplot(data, aes(x = beta, y = Vd, colour = sex)) +
  geom_point() +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  geom_smooth(method = "lm", color = "navy") +
  labs(
    title = "Figure 10: Relationship between β and Vd",
    subtitle = paste("Correlation =", round(cor(data$beta, data$Vd, use = "complete.obs"), 3)),
    x = "β elimination rate (g/kg/h)",
    y = "Vd (L/kg)"
  )

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALT APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Visualize the joint distribution

ggplot(data, aes(x = beta, y = Vd)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_density_2d(color = "navy") +
  # Add the marginal 97.5th percentiles
  geom_vline(xintercept = quantile(data$beta, 0.975, na.rm = TRUE), 
             linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = quantile(data$Vd, 0.975, na.rm = TRUE), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Joint Distribution of β and Vd",
    subtitle = "Red lines show marginal 97.5th percentiles",
    x = "β (g/kg/h)",
    y = "Vd (L/kg)"
  ) +
  annotate("text", 
           x = quantile(data$beta, 0.975, na.rm = TRUE), 
           y = min(data$Vd, na.rm = TRUE),
           label = "97.5th percentile β", 
           hjust = -0.1, color = "red")

