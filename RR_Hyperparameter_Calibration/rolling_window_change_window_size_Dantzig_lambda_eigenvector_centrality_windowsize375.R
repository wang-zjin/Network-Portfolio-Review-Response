# We want to adjust the window size, make it 375, 500, 750, 1000
setwd("~/Documents/GitHub/Network-Portfolio/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

estimate_rollingwindow_dantzig_lambda_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=375)
# Tuning parameter of Dantzig estimation: 232623.287 sec elapsed
estimate_rollingwindow_dantzig_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=375)
# Estmate Rolling Window Eigenvector Centrality by Dantzig-type selector: 346.918 sec elapsed

estimate_dantzig_lambda1_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=375)
estimate_dantzig_lambda2_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=375)
# CV hyperparameter lambda2 for Portfolio Parameter by Dantzig-type selector: 2057.196 sec elapsed
estimate_dantzig_lambda3_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=375)
estimate_glasso_rho_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=375)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector: 103.486 sec elapsed

