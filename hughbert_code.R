# load packages
library(tidyverse)
library(readxl)
library(ggplot2)
library(viridis)
library(patchwork)
library(kableExtra)

# --- data wrangling ---
data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")
data$sex <- as.factor(data$Sex)
data$beta60 <- data$`Beta60 (g/kg/h)`
data$weight <- data$`Weight (kg)`
data$height <- data$`Height (cm)`
data$age <- data$`Age (years)`
data$beta <- -data$beta60

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- beta60 EDA ---

# beta60 vs weight
weight_plot <- ggplot(data, aes(x=weight, y=beta60, colour = sex)) +
  geom_point() +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(title="Figure 1: Beta60 vs Weight", x="Weight (kg)", y="Beta60 (g/kg/h)")

# beta60 vs height
height_plot <- ggplot(data, aes(x=height, y=beta60, colour = sex)) +
  geom_point() +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) +
  labs(title="Figure 2: Beta60 vs Height", x="Weight (cm)", y="Beta60 (g/kg/h)")

# beta60 vs age
age_plot <- ggplot(data, aes(x=age, y=beta60, colour = sex)) +
  geom_point() +
  geom_smooth(method="lm", color="navy") +
  scale_colour_manual(values = c("male" = "steelblue", "female" = "lightblue")) 
  labs(title="Figure 3: Beta60 vs Age", x="Age (years)", y="Beta60 (g/kg/h)")

# beta60 vs gender
sex_plot <- ggplot(data, aes(x = sex, y = beta60,
                             fill = sex)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, color = "navy", alpha = 0.7) +
  scale_fill_manual(values = c("lightblue", "steelblue")) +
  labs(
    title = "Figure 4: Beta60 vs Gender",
    x = "Gender",
    y = "Beta60 (g/kg/h)"
  ) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

(weight_plot + height_plot) /
  (age_plot + sex_plot) /
  plot_layout(ncol = 1) +
  plot_annotation(
    title = "Beta60 vs Characteristics Plots",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# pairwise correlations between characteristics and beta60
data %>% 
  select(beta60, weight, age, height) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(3)

# quantiles of beta60
quantile(data$beta60, probs = c(0.025, 0.5, 0.975))

# --- linear modelling beta 60 attempt ---

# model incorporating all characteristics
model <- lm(beta60 ~ weight + age + height + sex, data = data)
summary(model)
http://127.0.0.1:39811/graphics/ad1037ac-b58b-43d8-8acf-eb01ec5f6dea.png
# Residual plots
par(mfrow = c(2,2))
plot(model)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# test person info (given in question)
test_person <- data.frame(weight = 70, height = 160, age = 70, sex = "female")

# predicted Beta60 using from regression model
pred_vals <- predict(model, newdata = test_person, interval = "prediction", level = 0.95)
pred_vals

# Ct and t from question
Ct <- 0.15
t  <- 2

# empirical quantiles for Beta60 (population approach)
beta_pop <- quantile(data$beta60, probs = c(0.025, 0.5, 0.975))
beta_pop

# C0 values using population quantiles
C0_pop <- Ct + beta_pop * t
C0_pop_abs <- Ct + abs(beta_pop) * t   # absolute values since beta60 is negative

# C0 values using model prediction
beta_pred <- as.numeric(pred_vals)
names(beta_pred) <- c("Fitted", "Lower (PI 2.5%)", "Upper (PI 97.5%)")

C0_model <- Ct + beta_pred * t
C0_model_abs <- Ct + abs(beta_pred) * t

# combine results into one comparison table 
results_table <- tibble(
  Approach = c("Population (empirical quantiles)", "Model-based (predicted)"),
  Lower_2.5 = c(C0_pop_abs[1], C0_model_abs[2]),  # note ordering: lwr=2nd col in pred
  Central   = c(C0_pop_abs[2], C0_model_abs[1]),
  Upper_97.5 = c(C0_pop_abs[3], C0_model_abs[3])
)

kable(results_table, digits = 3,
      col.names = c("Approach", "Lower (2.5%)", "Central", "Upper (97.5%)"),
      caption = "Comparison of estimated C₀ values for the population and model-based approaches")

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
beta_range <- quantile(data$beta60, probs = c(0.025, 0.975))
Vd_range   <- quantile(data$Vd,     probs = c(0.025, 0.975))


