# Load necessary libraries
library(ggplot2)

# Define data
samples <- c("ZC1_S4", "ZC2_S5", "ZC3_S6", "ZC4_S7", "ZP1_S20", "ZP2_S21", "ZP3_S22")
unclassified_reads <- c(29.59, 40.58, 37.58, 41.50, 33.24, 38.00, 39.30)
missing_reads <- c(25.74, 36.72, 33.73, 37.64, 29.37, 34.14, 35.42)

# Create a dataframe for correlation plot
df_correlation <- data.frame(Sample = samples, 
                             Unclassified_Reads = unclassified_reads, 
                             Missing_Reads = missing_reads)

# Generate scatter plot with regression line
correlation_plot <- ggplot(df_correlation, aes(x = Unclassified_Reads, y = Missing_Reads, label = Sample)) +
  geom_point(size = 4, color = "blue") +  # Scatter plot points
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
  geom_text(vjust = -1, size = 5) +  # Add sample labels
  labs(title = "Correlation between Unclassified and Missing Reads",
       x = "Unclassified Reads",
       y = "Missing Reads") +
  theme_minimal()

# Save the correlation plot
ggsave("correlation_plot.png", correlation_plot, width = 8, height = 6, dpi = 300)
print(correlation_plot)  # Display plot

# ---- Dual Line Plot ----
# Transform data into long format for ggplot2
df_long <- data.frame(
  Sample = rep(samples, 2),
  Reads = c(unclassified_reads, missing_reads),
  Type = rep(c("Unclassified Reads", "Missing Reads"), each = length(samples))
)

# Ensure Sample is a factor to maintain correct order
df_long$Sample <- factor(df_long$Sample, levels = samples)

# Plot both variables on the same graph
dual_line_plot <- ggplot(df_long, aes(x = Sample, y = Reads, group = Type, color = Type)) +
  geom_line(linewidth = 1) +  # Use linewidth instead of size
  geom_point(size = 3) +  # Points for each data value
  labs(title = "Unclassified vs Missing Reads",
       x = "Sample",
       y = "Read Count",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_y_continuous(limits = c(20, 45))  # Fixed y-axis range

# Save the dual-line plot
ggsave("dual_line_plot.png", dual_line_plot, width = 8, height = 6, dpi = 300)
print(dual_line_plot)  # Display plot