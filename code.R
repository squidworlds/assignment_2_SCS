suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(afex))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(cv))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(broom))

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
  geom_point(size = 3) +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  labs(title="Figure 1: β vs Weight", x="Weight (kg)", y="β Elimination Rate (g/kg/h)")

# beta vs height
height_plot <- ggplot(data, aes(x=height, y=beta, colour = sex)) +
  geom_point(size = 3) +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  labs(title="Figure 2: β vs Height", x="Height (cm)", y="β Elimination Rate (g/kg/h)")

# beta vs age
age_plot <- ggplot(data, aes(x=age, y=beta, colour = sex)) +
  geom_point(size = 3) +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  labs(title="Figure 3: β vs Age", x="Age (years)", y="β Elimination Rate (g/kg/h)")

# beta vs gender
sex_plot <- ggplot(data, aes(x = sex, y = beta60,
                             fill = sex)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("seagreen", "skyblue")) +
  labs(
    title = "Figure 4: β vs Gender",
    x = "Gender",
    y = "β Elimination Rate (g/kg/h)"
  ) +
  guides(fill = "none")

# Compute correlation and p-value between beta and each characteristic.
corr_table <- data %>%
  dplyr::select(beta, weight, age, height) %>%
  pivot_longer(cols = c(weight, age, height), names_to = "Characteristic", values_to = "Value") %>%
  group_by(Characteristic) %>%
  summarise(
    Correlation = cor(beta, Value, use = "pairwise.complete.obs"),
    P_value = cor.test(beta, Value)$p.value
  ) %>%
  mutate(
    Correlation = round(Correlation, 3),
    P_value = signif(P_value, 3)
  )

#%%%%%%%%%%%%%%%%%%%% BETA DISTRIBUTION TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Density plot of beta.
density_plot <- ggplot(data, aes(x = beta)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "steelblue",
                 alpha = 0.55) +
  geom_density(color = "navy", 
               size = 1) +
  geom_vline(aes(xintercept=mean(beta)),
             linetype = "dashed",
             color = "red",
             size = 2) +
  geom_vline(aes(xintercept=quantile(beta, 0.025)),
             linetype = "dotted",
             color = "red",
             size = 2) +
  geom_vline(aes(xintercept=quantile(beta, 0.975)),
             linetype = "dotted",
             color = "red",
             size = 2) +
  labs(x = expression(""*beta), y = "Density",
       title = expression("Distribution of "*beta))

# compute group statistics
stats <- data %>%
  group_by(sex) %>%
  summarize(
    mean = mean(beta),
    q025 = quantile(beta, 0.025),
    q975 = quantile(beta, 0.975)
  )

# helper: color mapping for male/female lines
stats <- stats %>%
  mutate(line_colour = ifelse(sex == "female", "red", "blue"))

# custom legend labels for quantile lines
legend_lines <- data.frame(
  sex = c("male", "female"),
  color = c("blue", "red"),
  linetype = c("dotted", "dotted"),
  label = c("97.5th quantile (male)", "97.5th quantile (female)")
)

# main plot
sex_density_plot <- ggplot(data, aes(x = beta, fill = sex)) +
  geom_histogram(aes(y = after_stat(density), color = sex),
                 binwidth = 0.005, alpha = 0.5, position = "identity") +
  geom_density(aes(color = sex), size = 1, alpha = 0) +
  
  # mean lines (dashed, red/blue)
  geom_vline(data = stats,
             aes(xintercept = mean, color = line_colour, linetype = "Mean"),
             size = 2) +
  
  # lower quantile (dotted)
  geom_vline(data = stats,
             aes(xintercept = q025, color = line_colour, linetype = "2.5th quantile"),
             size = 2) +
  
  # upper quantile (dotted)
  geom_vline(data = stats,
             aes(xintercept = q975, color = line_colour, linetype = "97.5th quantile"),
             size = 2) +
  
  # consistent manual color/fill mapping
  scale_fill_manual(values = c("male" = "skyblue", "female" = "seagreen"),
                    name = "Sex") +
  scale_color_manual(
    name = "Quantile Colour",
    values = c("male" = "skyblue", "female" = "seagreen",
               "red" = "red", "blue" = "blue"),
    breaks = c("blue", "red"),
    labels = c("male quantiles (blue)", "female quantiles (red)")
  ) +
  scale_linetype_manual(
    name = "Statistic",
    values = c("Mean" = "dashed", "2.5th quantile" = "dotted", "97.5th quantile" = "dotted")
  ) +
  labs(
    x = expression(beta),
    y = "Density",
    title = expression("Distribution of "*beta*" by Sex")
  )

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

gamma_df <- data.frame(
  x = data$beta,
  density = dgamma(data$beta, shape = gamma_shape, rate = gamma_rate)
)

quantiles <- qgamma(c(0.025, 0.5, 0.975), gamma_shape, gamma_rate)

stats_gamma <- data.frame(
  Quantile = rep(c("q025", "q50", "q975")),
  value = quantiles
)

# plot gamma fit for beta.
gamma_plot <- ggplot(data, aes(x = beta)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.005,
    fill = "steelblue",
    alpha = 0.5,
    position = "identity") +
  geom_line(
    data = gamma_df,
    aes(x = x,
        y = density
    ),
    color = "navy",
    size = 1.2) +
  geom_vline(
    data = stats_gamma,
    aes(xintercept = value,
        linetype = Quantile),
    color = "red",
    size = 2
  ) +
  scale_linetype_manual(
    values = c("q025" = "dotted", "q50" = "dashed", "q975" = "dotted"),
    labels = c("2.5%", "50%", "97.5%")
  ) +
  labs(
    x = expression(beta),
    y = "Density",
    title = expression("Gamma Distribution of "*beta),
    linetype = "Quantile"
  ) 

m_data <- data %>% filter(Sex == "male")
m_gamma_fit <- fitdist(m_data$beta, distr = "gamma", method = "mle")
m_gamma_shape <- m_gamma_fit$estimate[[1]]
m_gamma_rate <- m_gamma_fit$estimate[[2]]
m_quantiles <- qgamma(c(0.025, 0.5, 0.975), m_gamma_shape, m_gamma_rate)

f_data <- data %>% filter(Sex == "female")
f_gamma_fit <- fitdist(f_data$beta, distr = "gamma", method = "mle")
f_gamma_shape <- f_gamma_fit$estimate[[1]]
f_gamma_rate <- f_gamma_fit$estimate[[2]]
f_quantiles <- qgamma(c(0.025, 0.5, 0.975), f_gamma_shape, f_gamma_rate)

gamma_sex_df <- data.frame(
  x = data$beta,
  male_density   = dgamma(data$beta, shape = m_gamma_shape, rate = m_gamma_rate),
  female_density = dgamma(data$beta, shape = f_gamma_shape, rate = f_gamma_rate)
)

gamma_sex_long <- gamma_sex_df %>%
  pivot_longer(cols = c(male_density, female_density),
               names_to = "Sex",
               values_to = "density") %>%
  mutate(Sex = ifelse(Sex == "male_density", "male", "female"))

stats_gamma_sex <- data.frame(
  Sex = rep(c("male", "female"), each = 3),
  Quantile = rep(c("q025", "q50", "q975"), times = 2),
  value = c(m_quantiles, f_quantiles)
)

# gamma fit of beta by sex
sex_gamma_plot <- ggplot(data, aes(x = beta)) +
  # histogram with correct fill colors
  geom_histogram(
    aes(y = after_stat(density), fill = Sex),
    binwidth = 0.005, alpha = 0.5, position = "identity"
  ) +
  
  # gamma density curves (match fill colours)
  geom_line(
    data = gamma_sex_long,
    aes(x = x, y = density, color = Sex, linetype = "Gamma density"),
    size = 1.2
  ) +
  
  # quantile lines (male = blue, female = red)
  geom_vline(
    data = stats_gamma_sex %>%
      mutate(line_colour = ifelse(Sex == "female", "red", "blue")),
    aes(xintercept = value, color = line_colour, linetype = Quantile),
    size = 2
  ) +
  
  # color scales (fill = distribution colour, line = quantile/gamma)
  scale_fill_manual(
    values = c("male" = "skyblue", "female" = "seagreen"),
    name = "Sex"
  ) +
  scale_color_manual(
    name = "Quantile Colour",
    values = c("male" = "skyblue", "female" = "seagreen",
               "red" = "red", "blue" = "blue"),
    breaks = c("skyblue", "seagreen", "blue", "red"),
    labels = c("Male gamma curve", "Female gamma curve",
               "Male quantiles (blue)", "Female quantiles (red)")
  ) +
  
  # line type labels for quantiles
  scale_linetype_manual(
    name = "Statistic",
    values = c("Gamma density" = "solid",
               "q025" = "dotted",
               "q50"  = "dashed",
               "q975" = "dotted"),
    labels = c("Gamma density" = "Gamma density curve",
               "q025" = "2.5th quantile",
               "q50"  = "50th quantile (median)",
               "q975" = "97.5th quantile")
  ) +
  
  labs(
    x = expression(beta),
    y = "Density",
    title = expression("Gamma Distributions of "*beta*" by Sex")
  )

#%%%%%%%%%%%%%%%%%%%%%%% LINEAR MODELLING BETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  geom_point( size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = expression("Comparing Models for the 2.5% Quantile for "*beta),
    x = expression("Actual "*beta),
    y = expression("Predicted "*beta),
    color = "Model Type"
  ) +
  scale_color_manual(values = c("lm_q025" = "orange", "rq_pred" = "darkred"),
                     labels = c("Linear Regression", "Quantile Regression"))

# compute group statistics
stats <- data %>%
  group_by(Sex) %>%
  summarize(
    mean = mean(beta60),
    q025 = quantile(beta60, 0.025),
    q975 = quantile(beta60, 0.975)
  )

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

# %%%%%%%%%%%%%%%%%%%%%%%%%% ALTERNATE RESULT APPROACHES %%%%%%%%%%%%%%%%%%%%

# Given parameters
Ct <- 0.15         # measured concentration (g/kg)
t  <- 2            # time since arrest (hours)
x  <- 0.47         # legal limit (g/kg)

# Fitted Gamma parameters for beta
gamma_shape <- 31.9538
gamma_rate  <- 173.7079

# 1. Compute 2.5th and 97.5th percentiles of beta
beta_q025 <- qgamma(0.025, shape = gamma_shape, rate = gamma_rate)
beta_q975 <- qgamma(0.975, shape = gamma_shape, rate = gamma_rate)

# 2. Compute 95% interval for C0
C0_lower <- Ct + beta_q025 * t
C0_upper <- Ct + beta_q975 * t
C0_interval <- c(C0_lower, C0_upper)

# 3. Compute probability that true C0 exceeded legal limit
beta_threshold <- (x - Ct) / t
p_over_limit <- 1 - pgamma(beta_threshold, shape = gamma_shape, rate = gamma_rate)

# 4. Monte Carlo simulation
set.seed(123)
n_sim <- 1000
beta_samples <- rgamma(n_sim, shape = gamma_shape, rate = gamma_rate)
C0_samples <- Ct + beta_samples * t
p_over_limit_sim <- mean(C0_samples > x)

# Summaries
C0_mean <- mean(C0_samples)
C0_sd   <- sd(C0_samples)

# Create summary table
results <- data.frame(
  Statistic = c(
    "Lower 95% bound for C₀",
    "Upper 95% bound for C₀",
    "Analytical P(C₀ > 0.47 g/kg)",
    "Simulated P(C₀ > 0.47 g/kg)",
    "Mean of simulated C₀",
    "SD of simulated C₀"
  ),
  Value = c(
    round(C0_lower, 3),
    round(C0_upper, 3),
    round(p_over_limit, 3),
    round(p_over_limit_sim, 3),
    round(C0_mean, 3),
    round(C0_sd, 3)
  )
)

# Parameters
gamma_shape <- 31.9538
gamma_rate <- 173.7079

# Compute 2.5th and 97.5th percentiles of beta
beta_q025 <- qgamma(0.025, shape = gamma_shape, rate = gamma_rate)
beta_q975 <- qgamma(0.975, shape = gamma_shape, rate = gamma_rate)

# Current concentration
C_t <- 0.15
t <- 2

# Legal limit for C0
C_legal <- 0.47        

# Derived beta that would reach legal limit
beta_limit <- (C_legal - C_t)/t

# Define range of beta
# beta <- seq(0, 0.5, length.out = 1000)
beta <- seq(beta_q025, beta_q975, length.out = 1000)

# Calculate C0
C0 <- C_t + beta * t

# Gamma PDF for beta
pdf_beta <- dgamma(beta, shape = gamma_shape, rate = gamma_rate)

# Plot C0 vs beta
plot(beta, C0, type = "l", lwd = 2, col = "blue",
     xlab = expression(beta), ylab = expression(C[0]),
     main = expression(paste(C[0], " vs ", beta)))

# Add vertical line at derived beta_limit
abline(v = beta_limit, col = "darkgreen", lwd = 2, lty = 2)

# Shading for beta values that make C0 > C_legal
beta_shade <- beta[beta > beta_limit]
C0_shade <- C0[beta > beta_limit]

polygon(c(beta_shade, rev(beta_shade)),
        c(rep(min(C0), length(beta_shade)), rev(C0_shade)),
        col = rgb(1, 0, 0, 0.3), border = NA)

# Add Gamma PDF on secondary y-axis
par(new = TRUE)
plot(beta, pdf_beta, type = "l", lwd = 2, col = "red",
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, col = "red", col.axis = "red")
mtext("Density", side = 4, line = 3, col = "red")

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
density_plot <- ggplot(data, aes(x = Vd)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "steelblue",
                 alpha = 0.55) +
  geom_density(color = "navy", 
               size = 1) +
  geom_vline(aes(xintercept=mean(Vd)),
             linetype = "dashed",
             color = "red",
             size = 2) +
  geom_vline(aes(xintercept=quantile(Vd, 0.025)),
             linetype = "dotted",
             color = "red",
             size = 2) +
  geom_vline(aes(xintercept=quantile(Vd, 0.975)),
             linetype = "dotted",
             color = "red",
             size = 2) +
  labs(x = "Vd (L/kg)", y = "Density", title = expression("Distribution of Vd"))

# Vd vs weight
weight_plot2 <- ggplot(data, aes(x = weight, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", colour = "navy") +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  labs(title = "Figure 6: Vd vs Weight", x = "Weight (kg)", y = "Vd (L/kg)")

# Vd vs height
height_plot2 <- ggplot(data, aes(x = height, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  geom_smooth(method = "lm", colour = "navy") +
  labs(title = "Figure 7: Vd vs Height", x = "Height (cm)", y = "Vd (L/kg)")

# Vd vs age
age_plot2 <- ggplot(data, aes(x = age, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  geom_smooth(method = "lm", colour = "navy") +
  labs(title = "Figure 8: Vd vs Age", x = "Age (years)", y = "Vd (L/kg)")

# Vd vs gender
sex_plot2 <- ggplot(data, aes(x = sex, y = Vd,
                              fill = sex)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("seagreen", "skyblue")) +
  labs(
    title = "Figure 9: Vd vs Gender",
    x = "Gender",
    y = "Vd (L/kg)"
  ) +
  guides(fill = "none")

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
  geom_point(alpha = 0.6, color = "orange", size = 3) +
  geom_density_2d(color = "red") +
  # Add the marginal 97.5th percentiles
  geom_vline(xintercept = quantile(data$beta, 0.975, na.rm = TRUE), 
             linetype = "dashed", color = "blue", linewidth = 1) +
  geom_hline(yintercept = quantile(data$Vd, 0.975, na.rm = TRUE), 
             linetype = "dashed", color = "blue", linewidth = 1) +
  labs(
    title = "Joint Distribution of β and Vd",
    subtitle = "Blue lines show marginal 97.5th percentiles",
    x = "β (g/kg/h)",
    y = "Vd (L/kg)"
  ) +
  annotate("text", 
           x = quantile(data$beta, 0.975, na.rm = TRUE), 
           y = min(data$Vd, na.rm = TRUE),
           label = "97.5th percentile β", size = 7, 
           hjust = -0.1, color = "blue")


# Test person (from task 2)
A <- mean(data$A) # Chose this kinda randomly but made the most sense to me 
weight <- 70
t <- 2

# Calculate Ct using empirical joint distribution
data$Ct_joint <- (A / (weight * data$Vd)) - data$beta * t

# Compute the quantiles of C_t
quantile(data$Ct_joint, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# Current (independent) approach
beta_ind <- quantile(data$beta, 0.975, na.rm = TRUE)
Vd_ind   <- quantile(data$Vd, 0.975, na.rm = TRUE)
Ct_independent <- (A / (weight * Vd_ind)) - beta_ind * t

# Table comparison
results_compare <- tibble(
  Method = c("Empirical joint (β,Vd)", "Independent 97.5th percentiles"),
  Lower_2.5 = c(round(quantile(data$Ct_joint, 0.025, na.rm = TRUE), 3), ""),
  Median     = c(round(quantile(data$Ct_joint, 0.5, na.rm = TRUE), 3), ""),
  Upper_97.5 = c(round(quantile(data$Ct_joint, 0.975, na.rm = TRUE), 3), 
                 round(Ct_independent, 3))
)
