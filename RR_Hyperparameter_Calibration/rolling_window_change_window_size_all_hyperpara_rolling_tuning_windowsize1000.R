rm(list = ls())

# We want to adjust the window size, make it 375, 500, 750, 1000
setwd("~/Documents/GitHub/Network-Portfolio/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

estimate_rollingwindow_dantzig_lambda_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# Tuning parameter of Dantzig estimation: 128616.124 sec elapsed (updated)
estimate_rollingwindow_dantzig_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=1000)
# Estmate Rolling Window Eigenvector Centrality by Dantzig-type selector: 160.066 sec elapsed (updated)

estimate_rollingwindow_dantzig_lambda1_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=1000)
# Tuning parameter of Dantzig estimation: 50247.999 sec elapsed (updated)
estimate_rollingwindow_dantzig_lambda2_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=1000)
# Tuning parameter of Dantzig estimation: 58047.227 sec elapsed (updated)
estimate_rollingwindow_dantzig_lambda3_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=1000)
# Tuning parameter of Dantzig estimation: 38946.124 sec elapsed (updated)
estimate_rollingwindow_glasso_rho_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                      window_size=1000)
# Tuning parameter of glasso estimation: 3755.784 sec elapsed (updated)

estimate_rollingwindow_portfolio_parameter(EC_file_name = "eigenvector_centrality_WS1000_20250117.RData", 
                             DS_lambda1_file_name = "dantzig_lambda1_portfolio_parameter_WS1000_rollingwindow_20250126.RData", 
                             DS_lambda2_file_name = "dantzig_lambda2_portfolio_parameter_WS1000_rollingwindow_20250126.RData", 
                             DS_lambda3_file_name = "dantzig_lambda3_portfolio_parameter_WS1000_rollingwindow_20250126.RData", 
                             window_size=1000)
# Estimate Portfolio Parameters rho1, rho2, rho3 by Dantzig-type selector: 502.721 sec elapsed (updated)
