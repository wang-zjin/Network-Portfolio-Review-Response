rm(list = ls())

# We want to adjust the window size, make it 375, 500, 750, 1000
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

estimate_rollingwindow_dantzig_lambda_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=750)
# Tuning parameter of Dantzig estimation: 119563.322 sec elapsed (updated)
estimate_rollingwindow_dantzig_eigenvectorcentrality(file_name = "SP500 securities_up_20230306.csv", window_size=750)
# Estmate Rolling Window Eigenvector Centrality by Dantzig-type selector: 275.584 sec elapsed (updated)

estimate_rollingwindow_dantzig_lambda1_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=750)
# CV hyperparameter lambda1 for Portfolio Parameter by Dantzig-type selector: 2506.124 sec elapsed
estimate_rollingwindow_dantzig_lambda2_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=750)
# CV hyperparameter lambda2 for Portfolio Parameter by Dantzig-type selector: 1733.059 sec elapsed
estimate_rollingwindow_dantzig_lambda3_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                           window_size=750)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector: 1808.611 sec elapsed 
estimate_rollingwindow_glasso_rho_portfolio_parameter(file_name = "SP500 securities_up_20230306.csv", 
                                                      window_size=750)
# CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector: 104.225 sec elapsed 

estimate_rollingwindow_portfolio_parameter(EC_file_name = "eigenvector_centrality_WS750_20250117.RData", 
                             DS_lambda1_file_name = "dantzig_lambda1_portfolio_parameter_WS750_20250117.RData", 
                             DS_lambda2_file_name = "dantzig_lambda2_portfolio_parameter_WS750_20250117.RData", 
                             DS_lambda3_file_name = "dantzig_lambda3_portfolio_parameter_WS750_20250117.RData", 
                             window_size=750)
# Estimate Portfolio Parameters rho1, rho2, rho3 by Dantzig-type selector: 551.884 sec elapsed
