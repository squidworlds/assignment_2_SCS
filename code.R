knitr::opts_chunk$set(echo = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(afex)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(viridis)
  library(gridExtra)
  library(patchwork)
  library(kableExtra)
  library(cv)
  library(stringr)
  library(fitdistrplus)
  library(readxl)
  library(quantreg)
  library(purrr)
  library(broom)
  library(triangle)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%% DATA WRANGLING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- data wrangling ---
data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")
data$sex <- as.factor(data$Sex)
data$beta60 <- data$`Beta60 (g/kg/h)`
data$weight <- data$`Weight (kg)`
data$height <- data$`Height (cm)`
data$age <- data$`Age (years)`
data$beta <- -data$beta60

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BETA EDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- beta EDA ---

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
    title = "Figure 4: β vs Sex",
    x = "Sex",
    y = "β Elimination Rate (g/kg/h)"
  ) +
  guides(fill = "none")


# Compute correlation and p-value between beta and each characteristic.
corr_table <- data %>%
  dplyr::select(beta, weight, age, height) %>%
  pivot_longer(cols = c(weight, age, height), 
               names_to = "Characteristic", values_to = "Value") %>%
  group_by(Characteristic) %>%
  summarise(
    Correlation = cor(beta, Value, use = "pairwise.complete.obs"),
    P_value = cor.test(beta, Value)$p.value
  ) %>%
  mutate(
    Characteristic = str_to_title(Characteristic),
    Correlation = round(Correlation, 3),
    P_value = signif(P_value)
  )

#%%%%%%%%%%%%%%%%%%%% BETA DISTRIBUTION TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

density_plot <- ggplot(data, aes(x = beta)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "steelblue",
                 alpha = 0.55) +
  geom_density(color = "navy", 
               size = 3) +
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
  labs(x = expression(""*beta), y = "Density", title = expression("Figure 5: Distribution of "*beta))

# define colours 
sex_cols <- c("male" = "skyblue", "female" = "seagreen")

# compute group statistics
stats <- data %>%
  group_by(sex) %>%
  summarise(
    mean  = mean(beta),
    q025  = quantile(beta, 0.025),
    q975  = quantile(beta, 0.975)
  )

# main plot
sex_density_plot <- ggplot(data, aes(x = beta, colour = sex)) +
  
  # density curves
  geom_density(size = 3) +
  
  # mean lines
  geom_vline(data = stats,
             aes(xintercept = mean, colour = sex, linetype = "Mean"),
             size = 2.2) +
  
  # lower quantile
  geom_vline(data = stats,
             aes(xintercept = q025, colour = sex, linetype = "2.5th quantile"),
             size = 2) +
  
  # upper quantile
  geom_vline(data = stats,
             aes(xintercept = q975, colour = sex, linetype = "97.5th quantile"),
             size = 2) +
  
  # unified colour control
  scale_colour_manual(values = sex_cols, name = "Sex") +
  scale_fill_manual(values = sex_cols, name = "Sex") +
  
  # consistent linetypes
  scale_linetype_manual(
    name   = "Statistic",
    values = c("Mean" = "dashed",
               "2.5th quantile" = "dotted",
               "97.5th quantile" = "dotted")
  ) +
  
  labs(
    x = expression(beta),
    y = "Density",
    title = expression("Figure 6: Distribution of " * beta * " by Sex")
  )

# Fit distributions
normal_fit <- fitdist(data$beta, distr = "norm", method = "mle")
beta_fit   <- fitdist(data$beta, distr = "beta", method = "mle")
gamma_fit  <- fitdist(data$beta, distr = "gamma", method = "mle")

plot.legend <- c("Normal", "Beta", "Gamma")

# Layout
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2), oma = c(0, 0, 5, 0), cex = 2)

# --- DENSITY (custom lines) ---
cols <- c("red", "green", "blue")
ltys <- c(1, 2, 3)
lwds <- rep(5, 5)

denscomp(
  list(normal_fit, beta_fit, gamma_fit),
  legendtext = plot.legend,
  fitcol = cols,
  fitlty = ltys,
  fitlwd = lwds,
  xlab = "beta",
  ylab = "Density",
  main = "Density Comparison",
  xlegend = "topright"
)

# --- QQ PLOT (use defaults) ---
qqcomp(
  list(normal_fit, beta_fit, gamma_fit),
  legendtext = plot.legend,
  fitcol = cols,
  main = "QQ Plot Comparison",
  xlegend = "bottomright"
)

# --- CDF PLOT (defaults) ---
cdfcomp(
  list(normal_fit, beta_fit, gamma_fit),
  legendtext = plot.legend,
  fitcol = cols,
  fitlty = ltys,
  fitlwd = lwds,
  main = "CDF Comparison",
  xlegend = "bottomright"
)

# --- PP PLOT (defaults) ---
ppcomp(
  list(normal_fit, beta_fit, gamma_fit),
  legendtext = plot.legend,
  fitcol = cols,
  main = "PP Plot Comparison",
  xlegend = "bottomright"
)

# --- Add main title across all plots ---
mtext("Figure 7: Comparison of Fitted Distributions", outer = TRUE, cex = 3, font = 2)

# Gamma approach

# Extract fitted parameters
gamma_shape <- gamma_fit$estimate[[1]]
gamma_rate  <- gamma_fit$estimate[[2]]

# Build smooth gamma curve over the observed range
gamma_curve <- data.frame(
  x = seq(min(data$beta), max(data$beta), length.out = 1000)
) %>%
  mutate(density = dgamma(x, shape = gamma_shape, rate = gamma_rate))

# Gamma quantiles
stats_gamma <- data.frame(
  Quantile = c("2.5%", "50%", "97.5%"),
  value = qgamma(c(0.025, 0.50, 0.975), gamma_shape, gamma_rate)
)

# Plot
gamma_plot <- ggplot(data, aes(x = beta)) +
  
  # Empirical histogram
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.005,
    fill = "steelblue", alpha = 0.5
  ) +
  
  # Fitted gamma density
  geom_line(
    data = gamma_curve,
    aes(x = x, y = density),
    colour = "navy",
    linewidth = 3
  ) +
  
  # Gamma quantile lines
  geom_vline(
    data = stats_gamma,
    aes(xintercept = value, linetype = Quantile),
    colour = "red",
    linewidth = 1.2
  ) +
  
  scale_linetype_manual(
    values = c("2.5%" = "dotted", "50%" = "dashed", "97.5%" = "dotted")
  ) +
  
  labs(
    x = expression(beta),
    y = "Density",
    title = expression("Figure 8: Empirical Distribution of " * beta * 
                         " with Fitted Gamma Model"),
    linetype = "Gamma Quantile"
  )

# Gamma by Sex

# Fit gamma models for each sex
m_data <- data %>% filter(Sex == "male")
m_gamma_fit  <- fitdist(m_data$beta, "gamma", method = "mle")
m_gamma_shape <- m_gamma_fit$estimate[[1]]
m_gamma_rate  <- m_gamma_fit$estimate[[2]]
m_quantiles <- qgamma(c(0.025, 0.5, 0.975), m_gamma_shape, m_gamma_rate)

f_data <- data %>% filter(Sex == "female")
f_gamma_fit  <- fitdist(f_data$beta, "gamma", method = "mle")
f_gamma_shape <- f_gamma_fit$estimate[[1]]
f_gamma_rate  <- f_gamma_fit$estimate[[2]]
f_quantiles <- qgamma(c(0.025, 0.5, 0.975), f_gamma_shape, f_gamma_rate)


# build density curves
m_curve <- tibble(
  x = seq(min(m_data$beta), max(m_data$beta), length.out = 1000),
  density = dgamma(x, m_gamma_shape, m_gamma_rate)
)

f_curve <- tibble(
  x = seq(min(f_data$beta), max(f_data$beta), length.out = 1000),
  density = dgamma(x, f_gamma_shape, f_gamma_rate)
)

# Quantile df for geom_vline 
m_stats <- tibble(
  Quantile = c("2.5%", "50%", "97.5%"),
  value = m_quantiles
)

f_stats <- tibble(
  Quantile = c("2.5%", "50%", "97.5%"),
  value = f_quantiles
)


# Male plot
male_plot <- ggplot(m_data, aes(beta)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "skyblue",
                 alpha = 0.4) +
  
  geom_line(data = m_curve,
            aes(x = x, y = density),
            colour = "skyblue4",
            linewidth = 3) +
  
  geom_vline(data = m_stats,
             aes(xintercept = value, linetype = Quantile),
             colour = "skyblue4",
             linewidth = 2) +
  
  scale_linetype_manual(
    values = c("2.5%" = "dotted",
               "50%"  = "dashed",
               "97.5%" = "dotted")
  ) +
  
  labs(
    title = "Male",
    x = expression(beta),
    y = "Density"
  ) 


# Female plot
female_plot <- ggplot(f_data, aes(beta)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "seagreen",
                 alpha = 0.4) +
  
  geom_line(data = f_curve,
            aes(x = x, y = density),
            colour = "seagreen4",
            linewidth = 3) +
  
  geom_vline(data = f_stats,
             aes(xintercept = value, linetype = Quantile),
             colour = "seagreen4",
             linewidth = 2) +
  
  scale_linetype_manual(
    values = c("2.5%" = "dotted",
               "50%"  = "dashed",
               "97.5%" = "dotted")
  ) +
  
  labs(
    title = "Female",
    x = expression(beta),
    y = "Density"
  ) 


#%%%%%%%%%%%%%%%%%%%%%%% LINEAR MODELLING BETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- linear modelling beta attempt ---

# model incorporating all characteristics
lm_model <- lm(beta ~ weight + height + sex, data = data)

# Residual plots
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2), oma = c(0, 0, 5, 0), cex = 2)

plot(lm_model, which = 1, caption = "", main = "Residuals vs Fitted", sub.caption = "", cex.main = 1.5)

plot(lm_model, which = 2, caption = "", main = "Q-Q Residuals", sub.caption = "", cex.main = 1.5)

plot(lm_model, which = 3, caption = "", main = "Scale Location", sub.caption = "", cex.main = 1.5)

plot(lm_model, which = 5, caption = "", main = "Residuals vs Leverage", sub.caption = "", cex.main = 1.5)

# quantile regression approach 

# fit models
lm_model <- lm(beta ~ weight + height + sex, data = data)
rq_model <- rq(beta ~ weight + height + sex, data = data, tau = 0.025)

# predictions for linear model
pred_int <- predict(lm_model, newdata = data, interval = "prediction", level = 0.95)

# create df for comparing predicted vs actual
pred_data <- data %>%
  mutate(
    lm_pred_mean = pred_int[, "fit"],   # mean OLS prediction
    lm_pred_lwr  = pred_int[, "lwr"],   # 2.5% prediction interval (lower)
    lm_pred_upr  = pred_int[, "upr"]    # 97.5% prediction interval (upper)
  )

# add empirical 2.5% quantile of residuals to find 2.5% quantile of linear model
resid_q025 <- quantile(residuals(lm_model), probs = 0.025)
pred_data <- pred_data %>%
  mutate(lm_q025 = lm_pred_mean + resid_q025)

# quantile regression predictions
pred_data$rq_pred <- predict(rq_model, newdata = data)

# reshape
plot_data <- pred_data %>%
  pivot_longer(
    cols = c(lm_q025, rq_pred),
    names_to = "model",
    values_to = "predicted_beta"
  )

# plot
quantile_v_lm_plot <- ggplot() +
  # Points for models (quantile LM, quantile regression)
  geom_point(
    data = plot_data,
    aes(x = beta, y = predicted_beta, color = model),
    alpha = 0.6, size = 2
  ) +
  # Linear trend lines for clarity
  geom_smooth(
    data = plot_data,
    aes(x = beta, y = predicted_beta, color = model),
    method = lm, se = TRUE, linewidth = 1
  ) +
  # Manual colors and labels
  scale_color_manual(
    values = c(
      "lm_q025"      = "gold",
      "rq_pred"      = "darkred"
    ),
    labels = c(
      "lm_q025"      = "Linear Model 2.5%",
      "rq_pred"      = "Quantile Linear Model 2.5%"
    )
  ) +
  # Labels and theme
  labs(
    title = expression("Figure 11: Actual vs Predicted 2.5% Quantile of " * beta),
    x = expression("Actual " * beta),
    y = expression("Predicted " * beta),
    color = "Model Type"
  )

# %%%%%%%%%%%%%%%%%%%%%% APPLYING VARIOUS APPROACHES %%%%%%%%%%%%%%%%%%%%%%%%%

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
# fit quantiles
rq_fit_quantiles <- rq(beta ~ weight + height + sex, data = data, tau = c(0.025, 0.5, 0.975))
# predict on new data
rq_pred_vals <- predict(rq_fit_quantiles, newdata = test_person, interval = "none")
C0_rq <- Ct + rq_pred_vals * t

# combine results into one comparison table 
results_table <- tibble(
  Approach = c("Empirical Quantiles", "Gamma Fitted", "Female Gamma Fitted", "Quantile Regression"),
  Lower_2.5 = c(C0_pop[1], C0_gamma[1], C0_gamma_f[1], C0_rq[1]),
  Central   = c(C0_pop[2], C0_gamma[2], C0_gamma_f[2], C0_rq[2]),
  Upper_97.5 = c(C0_pop[3], C0_gamma[3], C0_gamma_f[3], C0_rq[3])
)

# setup from question
test_person <- data.frame(weight = 70, height = 160, age = 70, sex = "female")
Ct <- 0.15
t  <- 2
x_limit <- 0.47

# 1. Current approach (fixed beta)
beta_fixed <- quantile(data$beta, probs = 0.025)
C0_fixed   <- Ct + beta_fixed * t
p_fixed    <- as.numeric(C0_fixed > x_limit) 

# 2. Gamma fitted approach to whole population
beta_q025_gamma <- qgamma(0.025, shape = gamma_shape, rate = gamma_rate)
beta_q50_gamma <- qgamma(0.5, shape = gamma_shape, rate = gamma_rate)
beta_q975_gamma <- qgamma(0.975, shape = gamma_shape, rate = gamma_rate)

C0_gamma_lower <- Ct + beta_q025_gamma * t
C0_gamma_med <- Ct + beta_q50_gamma * t
C0_gamma_upper <- Ct + beta_q975_gamma * t

# Probability analytically from Gamma CDF
beta_threshold <- (x_limit - Ct) / t
p_gamma <- 1 - pgamma(beta_threshold, shape = gamma_shape, rate = gamma_rate)

# 2. Gamma fitted approach to female population
beta_q025_gamma_f <- qgamma(0.025, shape = f_gamma_shape, rate = f_gamma_rate)
beta_q50_gamma_f <- qgamma(0.5, shape = f_gamma_shape, rate = f_gamma_rate)
beta_q975_gamma_f <- qgamma(0.975, shape = f_gamma_shape, rate = f_gamma_rate)

C0_gamma_lower_f <- Ct + beta_q025_gamma_f * t
C0_gamma_med_f <- Ct + beta_q50_gamma_f * t
C0_gamma_upper_f <- Ct + beta_q975_gamma_f * t

# Probability analytically from Gamma CDF
beta_threshold <- (x_limit - Ct) / t
p_gamma_f <- 1 - pgamma(beta_threshold, shape = f_gamma_shape, rate = f_gamma_rate)

# 4. Quantile regression approach
# Fit lower/upper quantile regression for beta
rq_lower <- rq(beta ~ weight + height + sex, data = data, tau = 0.025)
rq_med <- rq(beta ~ weight + height + sex, data = data, tau = 0.5)
rq_upper <- rq(beta ~ weight + height + sex, data = data, tau = 0.975)

# Predict beta quantiles for each observation in the dataset
beta_q025_qr <- as.numeric(predict(rq_lower, newdata = test_person))
beta_q50_qr <- as.numeric(predict(rq_med, newdata = test_person))
beta_q975_qr <- as.numeric(predict(rq_upper, newdata = test_person))

# Summarize C_0 interval
C0_qr_lower <- Ct + beta_q025_qr * t
C0_qr_med  <- Ct + beta_q50_qr * t
C0_qr_upper <- Ct + beta_q975_qr * t

# Approximate P(C_0 > limit) under QR for the given suspect information
set.seed(123)
n_sim <- 10000
# draw from triangular distribution
beta_samples_qr <- rtriangle(
  n = n_sim,
  a = beta_q025_qr,
  b = beta_q975_qr,
  c = beta_q50_qr
)
# apply formula for C_0
C0_samples_qr <- Ct + beta_samples_qr * t
p_qr <- mean(C0_samples_qr > x_limit, na.rm = TRUE)

# plot
results_df <- tibble(
  Method = c("Current (Fixed β = 0.126)", "Gamma Model", "Female Gamma Model", "Quantile Regression"),
  Lower  = c(C0_fixed, C0_gamma_lower, C0_gamma_lower_f, C0_qr_lower),
  Upper  = c(C0_fixed, C0_gamma_upper, C0_gamma_upper_f, C0_qr_upper),
  Mean   = c(C0_fixed, C0_gamma_med, C0_gamma_med_f, C0_qr_med),
  Prob   = c(p_fixed, p_gamma, p_gamma_f, p_qr)
)

# Plot
plot <- ggplot(results_df, aes(y = Method)) +
  # Interval bands (thin lines)
  geom_segment(aes(x = Lower, xend = Upper, y = Method, yend = Method,
                   colour = Method), linewidth = 2) +
  # Endpoints for clarity
  geom_point(aes(x = Lower, colour = Method), size = 3) +
  geom_point(aes(x = Upper, colour = Method), size = 3) +
  # Mean point
  geom_point(aes(x = Mean, colour = Method), size = 10, shape = 18) +
  # Legal limit line
  geom_vline(xintercept = x_limit, linetype = "dashed",
             colour = "darkgreen", linewidth = 1) +
  annotate("text", x = x_limit + 0.01, y = 3.3,
           label = "Legal limit = 0.47 g/kg",
           colour = "darkgreen", hjust = 0, size = 9) +
  annotate("rect", xmin = x_limit - 0.005, xmax = x_limit + 0.005,
           ymin = 0.5, ymax = 4.5, alpha = 0.1, fill = "darkgreen") +
  geom_text(
    data = results_df %>% filter(Method != first(Method)),
    aes(x = Upper, y = Method,
        label = sprintf("P(C[0] > 0.47) == %.3f", Prob)),
    parse = TRUE,
    hjust = 1,
    vjust = 1.5,
    size = 9,
    colour = "black"
  ) +
  scale_x_continuous(name = expression(C[0]~"(g/kg)"),
                     limits = c(0.35, 0.72),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("Current (Fixed β = 0.126)" = "grey60",
                                 "Gamma Model"             = "firebrick",
                                 "Female Gamma Model"      = "seagreen",
                                 "Quantile Regression"     = "hotpink"),
                      name = "Method") +
  ggtitle(expression("Figure 12: Comparison of Uncertainty Intervals for Reporting Frameworks for "*C[0]*" in Case Study"))


# Recommended Framework

# Create summary table
summary_table <- data.frame(
  Quantity = c(
    "Measured BAC (Ct)",
    "Time since driving (t)",
    "Legal limit",
    "95% uncertainty interval",
    "Probability C₀ exceeded legal limit"
  ),
  Value = c(
    "0.15 g/kg",
    "2 hours",
    "0.47 g/kg",
    "(0.43, 0.69) g/kg",
    "94%"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vd EDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# define A, Co, Vd
data$A <- data$`Amount of Alcohol Consumed (g)`
data$Co <- data$`Co (g/Kg)`
data$Vd <- data$A/(data$Co *data$weight)

# view summary and quantiles of Vd
#summary(data$Vd)
#quantile(data$Vd, probs = c(0.025, 0.5, 0.975))

# histogram of Vd
density_plot <- ggplot(data, aes(x = Vd)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.02,
                 fill = "steelblue",
                 alpha = 0.55) +
  geom_density(color = "navy", 
               size = 3) +
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
  labs(x = "Vd (L/kg)", y = "Density", title = expression("Figure 13: Distribution of Vd"))

# Vd vs weight
weight_plot2 <- ggplot(data, aes(x = weight, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", colour = "navy") +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  labs(title = "Figure 14: Vd vs Weight", x = "Weight (kg)", y = "Vd (L/kg)")

# Vd vs height
height_plot2 <- ggplot(data, aes(x = height, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  geom_smooth(method = "lm", colour = "navy") +
  labs(title = "Figure 15: Vd vs Height", x = "Height (cm)", y = "Vd (L/kg)")

# Vd vs age
age_plot2 <- ggplot(data, aes(x = age, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  geom_smooth(method = "lm", colour = "navy") +
  labs(title = "Figure 16: Vd vs Age", x = "Age (years)", y = "Vd (L/kg)")

# Vd vs gender
sex_plot2 <- ggplot(data, aes(x = sex, y = Vd,
                              fill = sex)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("seagreen", "skyblue")) +
  labs(
    title = "Figure 17: Vd vs Gender",
    x = "Gender",
    y = "Vd (L/kg)"
  ) +
  guides(fill = "none")

# Test correlation to investigate independence assumption
corr_test <- cor.test(data$beta, data$Vd, use = "complete.obs")

# Extract correlation coefficient and p-value
corr_coef <- corr_test$estimate
p_value   <- corr_test$p.value

# View quantiles of each coefficient for comparison
beta_range <- quantile(data$beta, probs = c(0.025, 0.975))
Vd_range   <- quantile(data$Vd, probs = c(0.025, 0.975))

# Scatter plot between beta and Vd
scatter <- ggplot(data, aes(x = beta, y = Vd, colour = sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("male" = "skyblue", "female" = "seagreen")) +
  geom_smooth(method = "lm", color = "navy") +
  labs(
    title = "Figure 18: Relationship between β and Vd",
    subtitle = paste0(
      "Correlation = ", round(corr_coef, 3),
      ",  p-value = ", format.pval(p_value, digits = 5, eps = .00001)
    ),
    x = "β elimination rate (g/kg/h)",
    y = "Vd (L/kg)"
  )

# Visualize the joint distribution

# Compute thresholds
beta_thr <- quantile(data$beta, 0.975, na.rm = TRUE)
Vd_thr   <- quantile(data$Vd, 0.975, na.rm = TRUE)

# Add shaded region
annotate("rect",
         xmin = beta_thr,
         xmax = Inf,
         ymin = Vd_thr,
         ymax = Inf,
         alpha = 0.1,
         fill = "blue")


joint <- ggplot(data, aes(x = beta, y = Vd)) +
  geom_point(alpha = 0.6, color = "orange", size = 3) +
  geom_density_2d(color = "red") +
  
  # Highlight top-right corner
  annotate("rect",
           xmin = beta_thr,
           xmax = Inf,
           ymin = Vd_thr,
           ymax = Inf,
           alpha = 0.1,
           fill = "blue") +
  
  # Add marginal 97.5th percentiles
  geom_vline(xintercept = beta_thr, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_hline(yintercept = Vd_thr,   linetype = "dashed", color = "blue", linewidth = 1) +
  
  labs(
    title = "Figure 19: Joint Distribution of β and Vd",
    subtitle = "Blue region shows values above joint 97.5th percentiles",
    x = "β (g/kg/h)",
    y = "Vd (L/kg)"
  ) +
  annotate("text", 
           x = beta_thr, 
           y = min(data$Vd, na.rm = TRUE),
           label = "97.5th percentile β", size = 7, 
           hjust = -0.1, color = "blue")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALT APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test person (from task 2)
A <- mean(data$A) # Chose this kinda randomly but made the most sense to me 
weight <- 70
t <- 2

# Calculate Ct using empirical joint distribution
data$Ct_joint <- (A / (weight * data$Vd)) - data$beta * t

# Compute the quantiles of C_t
Ct_quant <- quantile(data$Ct_joint, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# Current (independent) approach
beta_ind <- quantile(data$beta, 0.975, na.rm = TRUE)
Vd_ind   <- quantile(data$Vd, 0.975, na.rm = TRUE)
Ct_independent <- (A / (weight * Vd_ind)) - beta_ind * t

# Table comparison
results_compare <- tibble(
  Method = c("Empirical joint (β,Vd)", "Independent 97.5th percentiles"),
  Lower_2.5 = c(round(Ct_quant[1], 3), ""),
  Median     = c(round(Ct_quant[2], 3), ""),
  Upper_97.5 = c(round(Ct_quant[3], 3), 
                 round(Ct_independent, 3))
)
