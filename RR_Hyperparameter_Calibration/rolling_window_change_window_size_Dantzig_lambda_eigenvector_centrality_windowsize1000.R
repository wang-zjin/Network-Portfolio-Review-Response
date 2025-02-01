# We want to adjust the window size, make it 375, 500, 750, 1000
setwd("~/Documents/GitHub/Network-Portfolio/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

estimate_rollingwindow_dantzig_lambda_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# Tuning parameter of Dantzig estimation: 128616.124 sec elapsed (updated)
estimate_rollingwindow_dantzig_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# Estmate Rolling Window Eigenvector Centrality by Dantzig-type selector: 160.066 sec elapsed (updated)

estimate_dantzig_lambda1_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# CV hyperparameter lambda1 for Portfolio Parameter by Dantzig-type selector: 2300.064 sec elapsed (updated)
estimate_dantzig_lambda2_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# CV hyperparameter lambda2 for Portfolio Parameter by Dantzig-type selector: 1934.276 sec elapsed (updated)
estimate_dantzig_lambda3_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector: 1793.048 sec elapsed (updated)
estimate_glasso_rho_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector: 94.922 sec elapsed (updated)

estimate_portfolio_parameter(EC_file_name = "eigenvector_centrality_WS1000_20250117.RData", 
                             DS_lambda1_file_name = "dantzig_lambda1_portfolio_parameter_WS1000_20250117.RData", 
                             DS_lambda2_file_name = "dantzig_lambda2_portfolio_parameter_WS750_20250117.RData", 
                             DS_lambda3_file_name = "dantzig_lambda3_portfolio_parameter_WS750_20250117.RData", 
                             window_size=1000)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector: 94.922 sec elapsed (updated)

