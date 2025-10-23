library(tidyverse)
library(readxl)
library(fitdistrplus)

df <- read_excel("C:/Users/Saioa/OneDrive - University of Edinburgh/y4s1/SCS/assignment_2_SCS/SCS_BAC_and_BrAC_split_TOP.xlsx")

df$beta_60 <- df$`Beta60 (g/kg/h)`

density_plot <- ggplot(df, aes(x = beta_60)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.005,
                 fill = "lightgreen",
                 alpha = 0.55) +
  geom_density(color = "seagreen", 
               size = 1) +
  geom_vline(aes(xintercept=mean(beta_60)),
             linetype = "dashed",
             color = "seagreen",
             size = 1) +
  geom_vline(aes(xintercept=quantile(beta_60, 0.025)),
             linetype = "dashed",
             color = "seagreen",
             size = 1) +
  geom_vline(aes(xintercept=quantile(beta_60, 0.975)),
             linetype = "dashed",
             color = "seagreen",
             size = 1)

density_plot

quantile(df$beta_60, 0.025)

gamma_fit <- fitdist(-df$beta_60, distr = "gamma", method = "mle")
summary(gamma_fit)
par(mar=c(1, 1, 1, 1))
plot(gamma_fit)

beta_fit <- fitdist(df$beta_60, distr = "beta", method = "mle")
summary(beta_fit)
