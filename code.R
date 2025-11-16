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
lm_model <- lm(beta ~ weight + height + sex, data = data)

# Residual plots
par(mfrow = c(2,2))
plot(lm_model)

# Try modelling square root of beta (residual plots indicate
# a slight quadratic mean-variance relationship)

sqrt_lm_model <- lm(sqrt(beta) ~ weight + height + sex, data = data)

# Residual plots
par(mfrow = c(2,2))
plot(sqrt_lm_model)

# compute group statistics
stats <- data %>%
  group_by(Sex) %>%
  summarize(
    mean = mean(beta60),
    q025 = quantile(beta60, 0.025),
    q975 = quantile(beta60, 0.975)
  )

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
    cols = c(lm_pred_mean, lm_q025, rq_pred),
    names_to = "model",
    values_to = "predicted_beta"
  )

# plot
quantile_v_lm_plot <- ggplot() +
  # 95% prediction interval shaded area (between lm_pred_lwr and lm_pred_upr)
  geom_ribbon(
    data = pred_data,
    aes(x = beta, ymin = lm_pred_lwr, ymax = lm_pred_upr),
    fill = "lightblue", alpha = 0.25
  ) +
  # Points for models (mean LM, quantile LM, quantile regression)
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
  # 1:1 line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Manual colors and labels
  scale_color_manual(
    values = c(
      "lm_pred_mean" = "blue",
      "lm_q025"      = "gold",
      "rq_pred"      = "darkred"
    ),
    labels = c(
      "lm_pred_mean" = "Linear Model (Mean Prediction)",
      "lm_q025"      = "2.5% Quantile of Linear Model Residuals",
      "rq_pred"      = "Quantile Regression (τ = 0.025)"
    )
  ) +
  # Labels and theme
  labs(
    title = expression("Figure 12: Actual vs Predicted " * beta * " with 95% Prediction Interval"),
    x = expression("Actual " * beta),
    y = expression("Predicted " * beta),
    color = "Model Type"
  )


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    "Interval length",
    "Lower 95% bound for C₀",
    "Upper 95% bound for C₀",
    "Analytical P(C₀ > 0.47 g/kg)",
    "Simulated P(C₀ > 0.47 g/kg)",
    "Mean of simulated C₀",
    "SD of simulated C₀"
  ),
  Value = c(
    round(C0_upper - C0_lower, 3),
    round(C0_lower, 3),
    round(C0_upper, 3),
    round(p_over_limit, 3),
    round(p_over_limit_sim, 3),
    round(C0_mean, 3),
    round(C0_sd, 3)
  )
)

kable(
  results,
  caption = "Table 2: Summary of Back-Calculated Blood-Alcohol Concentration (C₀) Using Gamma Approach",
  col.names = c("Statistic", "Value"),
  align = c("l", "c"),
  format = "markdown"
)

# Fit quantile regression models for beta
rq_lower  <- rq(beta ~ weight + age + sex, data = data, tau = 0.025)
rq_upper  <- rq(beta ~ weight + age + sex, data = data, tau = 0.975)

# Predict quantile values of beta for each observation
beta_q025 = predict(rq_lower,  newdata = data)
beta_q975 = predict(rq_upper,  newdata = data)

# Compute 95% interval for C0 for each observation
C0_lower = Ct + mean(beta_q025) * t
C0_upper = Ct + mean(beta_q975) * t
C0_interval <- c(C0_lower, C0_upper)

# Monte Carlo-style resampling from the quantile predictions
set.seed(123)
n_sim <- 1000
# approximate sampling of beta between 2.5% and 97.5% quantiles
beta_samples <- runif(n_sim, min = mean(beta_q025), max = mean(beta_q975))
C0_samples <- Ct + beta_samples * t

# Probability that C0 exceeded the legal limit x
p_over_limit_sim <- mean(C0_samples > x)

# Summaries
C0_mean <- mean(C0_samples)
C0_sd   <- sd(C0_samples)

# Create summary table
results <- data.frame(
  Statistic = c(
    "Interval Length",
    "Lower 95% bound for C₀",
    "Upper 95% bound for C₀",
    "Simulated P(C₀ > 0.47 g/kg)",
    "Mean of simulated C₀",
    "SD of simulated C₀"
  ),
  Value = c(
    round(C0_upper - C0_lower, 3),
    round(C0_lower, 3),
    round(C0_upper, 3),
    round(p_over_limit_sim, 3),
    round(C0_mean, 3),
    round(C0_sd, 3)
  )
)

kable(
  results,
  caption = "Table 3: Summary of Back-Calculated Blood-Alcohol Concentration using Quantile Regression (C₀)",
  col.names = c("Statistic", "Value"),
  align = c("l", "c"),
  format = "markdown"
)

# Parameters
gamma_shape <- 31.9538
gamma_rate  <- 173.7079

# Quantiles for beta
beta_q025 <- qgamma(0.025, shape = gamma_shape, rate = gamma_rate)
beta_q975 <- qgamma(0.975, shape = gamma_shape, rate = gamma_rate)

# Case-specific parameters
C_t <- 0.15
t <- 2
C_legal <- 0.47

# Derived beta that yields legal limit for C0
beta_limit <- (C_legal - C_t) / t

# Sequence of beta values
beta <- seq(beta_q025, beta_q975, length.out = 1000)

# Corresponding C0 values
C0 <- C_t + beta * t

# Gamma PDF for beta
pdf_beta <- dgamma(beta, shape = gamma_shape, rate = gamma_rate)

# Plot C0 vs beta
plot(beta, C0, type = "l", lwd = 2, col = "blue",
     xlab = expression(beta), ylab = expression(C[0]),
     main = expression(paste("Figure 13: BAC at C_0 vs Elimination Rate (", beta, ")")))

# Add horizontal line for legal limit on C0 scale
abline(h = C_legal, col = "darkgreen", lwd = 2, lty = 2)

# Add vertical line at beta corresponding to legal limit
abline(v = beta_limit, col = "darkgreen", lwd = 2, lty = 2)

# Add marker for (beta_limit, C0 = 0.47)
points(beta_limit, C_legal, pch = 19, col = "darkgreen", cex = 1.5)
text(beta_limit, C_legal - 0.015, labels = expression(paste("(", beta[limit], ", ", C[legal], ")")), 
     pos = 4, col = "darkgreen")

# --- Overlay PDF curve (secondary y-axis) ---
par(new = TRUE)
plot(beta, pdf_beta, type = "n", axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(pdf_beta)))

# Shade area under the PDF curve for beta > beta_limit
shade_idx <- beta >= beta_limit
polygon(c(beta[shade_idx], rev(beta[shade_idx])),
        c(pdf_beta[shade_idx], rep(0, sum(shade_idx))),
        col = rgb(1, 0, 0, 0.3), border = NA)

# Draw the full PDF line
lines(beta, pdf_beta, col = "red", lwd = 2)

# # Add secondary axis for density
# axis(side = 4, col = "red", col.axis = "red")
# mtext("Density of β", side = 4, line = 3, col = "red")

# Add legend
legend("topleft",
       legend = c(expression(C[0] == C[t] + beta*t),
                  expression("Legal limit"),
                  expression("P(C"[0]*"> limit)")),
       col = c("blue", "darkgreen", "red"),
       lwd = c(2, 2, NA),
       pch = c(NA, NA, 15),
       pt.cex = 1.5,
       bty = "n",
       cex = 0.9,
       fill = c(NA, NA, rgb(1,0,0,0.3)))


# 1. Current approach (fixed beta)
beta_fixed <- 0.126
C0_fixed   <- Ct + beta_fixed * t
p_fixed    <- as.numeric(C0_fixed > x_limit) 

# 2. Gamma model (use report MLEs)
gamma_shape <- 31.9538
gamma_rate  <- 173.7079

beta_q025_gamma <- qgamma(0.025, shape = gamma_shape, rate = gamma_rate)
beta_q975_gamma <- qgamma(0.975, shape = gamma_shape, rate = gamma_rate)

C0_gamma_lower <- Ct + beta_q025_gamma * t
C0_gamma_upper <- Ct + beta_q975_gamma * t
C0_gamma_mean  <- mean(c(C0_gamma_lower, C0_gamma_upper))

# Probability analytically from Gamma CDF
beta_threshold <- (x_limit - Ct) / t
p_gamma <- 1 - pgamma(beta_threshold, shape = gamma_shape, rate = gamma_rate)

# 3. Quantile regression approach (population-average bands)
# Fit lower/upper quantile regression for beta
rq_lower <- rq(beta ~ weight + age + sex, data = data, tau = 0.025)
rq_upper <- rq(beta ~ weight + age + sex, data = data, tau = 0.975)

# Predict β quantiles for each observation in the dataset
beta_q025_qr <- as.numeric(predict(rq_lower, newdata = data))
beta_q975_qr <- as.numeric(predict(rq_upper, newdata = data))

# Summarize Co interval by averaging β quantiles across the population 
C0_qr_lower <- Ct + mean(beta_q025_qr, na.rm = TRUE) * t
C0_qr_upper <- Ct + mean(beta_q975_qr, na.rm = TRUE) * t
C0_qr_mean  <- mean(c(C0_qr_lower, C0_qr_upper))

p_qr <- p_over_limit_sim

# plot
results_df <- tibble(
  Method = c("Current (Fixed β=0.126)", "Gamma model", "Quantile regression"),
  Lower  = c(C0_fixed, C0_gamma_lower, C0_qr_lower),
  Upper  = c(C0_fixed, C0_gamma_upper, C0_qr_upper),
  Mean   = c(C0_fixed, C0_gamma_mean, C0_qr_mean),
  Prob   = c(p_fixed, p_gamma, p_qr)
)

# Plot
p <- ggplot(results_df, aes(y = Method)) +
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
           ymin = 0.5, ymax = 3.5, alpha = 0.1, fill = "darkgreen") +
  # Probability labels (skip Current approach)
  # Probability labels (skip Current approach)
  geom_text(
    data = subset(results_df, Method != "Current (Fixed β=0.126)"),
    aes(x = Upper, 
        y = Method, 
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
  scale_colour_manual(values = c("Current (Fixed β=0.126)" = "grey60",
                                 "Gamma model"             = "firebrick",
                                 "Quantile regression"     = "hotpink"),
                      name = "Method") +
  ggtitle(expression("Comparison of reporting frameworks for " * C[0]))

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
