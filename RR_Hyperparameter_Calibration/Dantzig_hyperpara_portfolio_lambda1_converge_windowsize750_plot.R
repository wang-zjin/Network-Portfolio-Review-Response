# Clear the workspace 
rm(list = ls())

# Set the working directory 
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

# Load Functions and packages
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

# Define the rolling window size
window_size = 750

# Load the precomputed Dantzig lambda values for eigenvector centrality
load(paste0("dantzig_lambda1_portfolio_parameter_WS",window_size,"_rollingwindow_20250126.RData"))

# Unlist the Dantzig lambda values and store them in a vector
lmd1.EC.Dantzig = unlist(dantzig_lambda1_portfolio_parameter)

# Create a sequence for the x-axis
x_750 = 1:length(lmd1.EC.Dantzig)

# Compute the rolling median of the Dantzig lambda values
median_lambda1_750 = sapply(1:length(lmd1.EC.Dantzig), function(i) median(lmd1.EC.Dantzig[1:i]))

# Convert your data into a data frame
plot_data <- data.frame(
  x = x_750,
  lambda1_Dantzig = lmd1.EC.Dantzig,  # Dantzig lambda values
  median_lambda1_750 = median_lambda1_750   # Rolling median values
)

# Create the plot using ggplot2
plot_dantzig <- ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = lambda1_Dantzig, color = "lambda_1 (WS750)"), size = 1) +  # First line
  geom_line(aes(y = median_lambda1_750, color = "median(lambda_1) (WS750)"), size = 1) +  # Second line
  scale_color_manual(values = c("lambda_1 (WS750)" = "blue", "median(lambda_1) (WS750)" = "red")) +  # Explicit colors
  labs(
    x = "",
    y = "Hyperparameter lambda_1 Value",
    color = "Legend"
  ) +
  theme(
    panel.grid = element_blank(), # Remove all grid lines
    panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
    plot.background = element_rect(fill = "transparent", color = NA),   # Make plot background transparent
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.title.y = element_text(size = 16),
    legend.position = "none"  # Remove the legend
  )

# Save the plot
ggsave(filename = "Dantzig_hyper_portfolio_lambda1_converge_rollingwindow_750.png",
       plot = plot_dantzig,
       width = 10,
       height = 6,
       dpi = 300)
