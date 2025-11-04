library(tidyverse)
library(readxl)
library(fitdistrplus)
library(viridis)
library(dplyr)
library(quantreg)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# data <- read_excel("C:/Users/Saioa/OneDrive - University of Edinburgh/y4s1/SCS/assignment_2_SCS/SCS_BAC_and_BrAC_split_TOP.xlsx")
data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")

data$Sex <- as.factor(data$Sex)

# beta60 vs weight
weight_plot <- ggplot(data, aes(x=`Weight (kg)`, y=`Beta60 (g/kg/h)`)) +
  geom_point(alpha=0.5, color=viridis(5)[3]) +
  geom_smooth(method="lm", color="navy") +
  labs(title="Figure 1: Beta60 vs Weight", x="Weight (kg)", y="Beta60 (g/kg/h)")

# beta60 vs height
height_plot <- ggplot(data, aes(x=`Height (cm)`, y=`Beta60 (g/kg/h)`)) +
  geom_point(alpha=0.5, color=viridis(5)[3]) +
  geom_smooth(method="lm", color="navy") +
  labs(title="Figure 2: Beta60 vs Height", x="Weight (cm)", y="Beta60 (g/kg/h)")

# beta60 vs age
age_plot <- ggplot(data, aes(x=`Age (years)`, y=`Beta60 (g/kg/h)`)) +
  geom_point(alpha=0.5, color=viridis(5)[3]) +
  geom_smooth(method="lm", color="navy") +
  labs(title="Figure 3: Beta60 vs Age", x="Age (years)", y="Beta60 (g/kg/h)")

# beta60 vs gender
sex_plot <- ggplot(data, aes(x = Sex, y = `Beta60 (g/kg/h)`,
                             fill = Sex)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("lightblue", "steelblue")) +
  labs(
    title = "Figure 4: Beta60 vs Gender",
    x = "Gender",
    y = "Beta60 (g/kg/h)"
  ) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# (weight_plot + height_plot) /
#   (age_plot + sex_plot) /
#   plot_layout(ncol = 1) +
#   plot_annotation(
#     title = "Beta60 vs Characteristics Plots",
#     theme = theme(plot.title = element_text(size = 16, face = "bold"))
#   )

#%%%%%%%%%%%%%%%%%%%% BETA DISTRIBUTION TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data$beta60 <- df$`Beta60 (g/kg/h)`

density_plot <- ggplot(data, aes(x = beta60)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "lightgreen",
                 alpha = 0.55) +
  geom_density(color = "seagreen", 
               size = 1) +
  geom_vline(aes(xintercept=mean(beta60)),
             linetype = "dashed",
             color = "seagreen",
             size = 1) +
  geom_vline(aes(xintercept=quantile(beta60, 0.025)),
             linetype = "dashed",
             color = "seagreen",
             size = 1) +
  geom_vline(aes(xintercept=quantile(beta60, 0.975)),
             linetype = "dashed",
             color = "seagreen",
             size = 1)

density_plot

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
