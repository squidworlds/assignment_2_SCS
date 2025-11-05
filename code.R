library(tidyverse)
library(readxl)
library(fitdistrplus)
library(viridis)
library(dplyr)
library(quantreg)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")
data$sex <- as.factor(data$Sex)
data$beta60 <- data$`Beta60 (g/kg/h)`
data$weight <- data$`Weight (kg)`
data$height <- data$`Height (cm)`
data$age <- data$`Age (years)`
data$beta <- -data$beta60

# --- beta EDA ---

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

quantile(data$beta60, 0.025)

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

gamma_quantiles <- -qgamma(c(0.975, 0.5, 0.025), gamma_shape, gamma_rate)#
gamma_quantiles[2]


# beta?
beta_fit <- fitdist(-data$beta60, distr = "beta", method = "mle")
summary(beta_fit)

# sex density plot
data <- read_excel("C:/Users/Saioa/OneDrive - University of Edinburgh/y4s1/SCS/assignment_2_SCS/SCS_BAC_and_BrAC_split_TOP.xlsx")
data$beta60 <- df$`Beta60 (g/kg/h)`
data$Sex <- as.factor(data$Sex)
data

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


# quantile regression
library(quantreg)
rqfit <- rq(sqrt((-1)*beta60) ~ weight + age + height + sex, data = data, tau = 0.025)
summary(rqfit)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# define A, Co, Vd
data$A <- data$`Amount of Alcohol Consumed (g)`
data$Co <- data$`Co (g/Kg)`
data$Vd <- data$A/(data$Co *data$weight)

# view summary and quantiles of Vd
summary(data$Vd)
quantile(data$Vd, probs = c(0.025, 0.5, 0.975))

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

# test correlation to investigate independent assumption
cor.test(data$beta60, data$Vd, use = "complete.obs")

# view quantiles of each coefficient for comparison
beta_range <- quantile(data$beta, probs = c(0.025, 0.975))
Vd_range   <- quantile(data$Vd,     probs = c(0.025, 0.975))
