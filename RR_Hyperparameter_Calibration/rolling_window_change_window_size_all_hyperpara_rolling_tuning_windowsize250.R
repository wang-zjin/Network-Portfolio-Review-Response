rm(list = ls())

# We want to adjust the window size, make it 375, 500, 750, 1000
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

estimate_rollingwindow_dantzig_lambda_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=250)
# Tuning parameter of Dantzig estimation: 15877.521 sec elapsed (updated)
estimate_rollingwindow_dantzig_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=250)
# Estmate Rolling Window Eigenvector Centrality by Dantzig-type selector: 161.398 sec elapsed (updated)

estimate_rollingwindow_dantzig_lambda1_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv",
                                                           window_size=250)
# CV hyperparameter lambda1 for Portfolio Parameter by Dantzig-type selector:
estimate_rollingwindow_dantzig_lambda2_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=250)
# CV hyperparameter lambda2 for Portfolio Parameter by Dantzig-type selector: 56329.52 sec elapsed (updated)
estimate_rollingwindow_dantzig_lambda3_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=250)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector:
estimate_rollingwindow_glasso_rho_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                      window_size=250)
# CV hyperparameter rho for Portfolio Parameter by Dantzig-type selector: 12297.907 sec elapsed (updated)

estimate_rollingwindow_portfolio_parameter(EC_file_name = "eigenvector_centrality_W250_20250117.RData", 
                             DS_lambda1_file_name = "dantzig_lambda1_portfolio_parameter_WS250_rollingwindow_20250126.RData", 
                             DS_lambda2_file_name = "dantzig_lambda2_portfolio_parameter_WS250_rollingwindow_20250126.RData", 
                             DS_lambda3_file_name = "dantzig_lambda3_portfolio_parameter_WS250_rollingwindow_20250126.RData", 
                             window_size=250)
# Estimate Portfolio Parameters rho1, rho2, rho3 by Dantzig-type selector: 551.884 sec elapsed
