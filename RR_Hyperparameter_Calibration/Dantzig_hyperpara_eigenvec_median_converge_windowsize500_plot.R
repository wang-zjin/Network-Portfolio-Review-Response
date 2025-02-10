# Clear the workspace 
rm(list = ls())

# Set the working directory 
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

# Load Functions and packages
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

# Define the rolling window size
window_size = 500

# Load the precomputed Dantzig lambda values for eigenvector centrality
load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))

# Unlist the Dantzig lambda values and store them in a vector
lmd.EC.Dantzig = unlist(lmd.EC.Dantzig.list)

# Create a sequence for the x-axis
x_500 = 1:length(lmd.EC.Dantzig)

# Compute the rolling median of the Dantzig lambda values
median_rho_500 = sapply(1:length(lmd.EC.Dantzig), function(i) median(lmd.EC.Dantzig[1:i]))

# Convert your data into a data frame
plot_data <- data.frame(
  x = x_500,
  lmd_EC_Dantzig = lmd.EC.Dantzig,  # Dantzig lambda values
  median_rho_500 = median_rho_500   # Rolling median values
)

# Create the plot using ggplot2
plot_dantzig <- ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = median_rho_500, color = "median(rho^E) (WS500)"), size = 1) +  
  scale_color_manual(values = c("median(rho^E) (WS500)" = "red")) +  # Explicit colors
  labs(
    x = "",
    y = "",
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
ggsave(filename = "Dantzig_hyper_eigenvector_median_converge_rollingwindow_500.png",
       plot = plot_dantzig,
       width = 10,
       height = 6,
       dpi = 300)
