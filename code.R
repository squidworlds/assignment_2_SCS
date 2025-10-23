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


(weight_plot + height_plot) /
  (age_plot + sex_plot) /
  plot_layout(ncol = 1) +
  plot_annotation(
    title = "Beta60 vs Characteristics Plots",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )