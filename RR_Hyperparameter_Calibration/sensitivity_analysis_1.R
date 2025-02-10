# Clear the workspace 
rm(list = ls())

# Set the working directory 
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

# Function to compute the L2 (Euclidean) norm difference between two vectors
l2_norm <- function(x) {
  sqrt(sum((x)^2))
}
relative_weight <- function(x,y){
  z=(y-x)/x
  return(z)
}
relative_error <- function(x,y){
  z=l2_norm(y-x)/l2_norm(x)
  return(z)
}

weight_2constraint_3constraint_network_port = function(theta1, theta2, theta3, phi_star, mu_star){
  
  ###### Dantzig 2 constraints network portfolio #####
  if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
    alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
    w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
  } else{
    w = theta1/sum(theta1)
  }
  w[is.na(w)] <- 0
  w_2constraint_dantzig<-w
  
  ###### Dantzig 3 constraints network portfolio #####
  if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
    if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
      # 3 constraint case
      M1 <- cbind(rbind(1,phi_star,mean(estimated_mu)),
                  rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
                  rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
      M2 <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
                  rbind(1,phi_star,mean(estimated_mu)),
                  rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
      M <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
                 rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
                 rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
      gamma1 <- det(M1)/det(M)
      gamma2 <- det(M2)/det(M)
      alpha1 <- c(gamma1*sum(theta3))
      alpha2 <- c(gamma2*sum(theta2))
      w1 <- (1-alpha1-alpha2)*theta1/sum(theta1)
      w2 <- alpha1*theta3/sum(theta3)
      w3 <- alpha2*theta2/sum(theta2)
      w = w1 +  w2 + w3
    }
    else{
      # 1 centrality constraint case
      alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
      w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
    }} else{
      if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
        # 1 expected return constraint case 
        alpha=c((sum(theta2)*sum(theta1)*mean(estimated_mu)-(sum(theta2))^2)/(estimated_mu%*%theta2*sum(theta1)-(sum(theta2))^2))
        w = alpha*theta2/sum(theta2)+(1-alpha)*theta1/sum(theta1)
      }
      else{
        # global minimum variance case 
        w = theta1/sum(theta1)
      }
    }
  w[is.na(w)] <- 0
  w_3constraint_dantzig<-w
  
  return(list("w_2constraint_dantzig"=w_2constraint_dantzig,
              "w_3constraint_dantzig"=w_3constraint_dantzig))
}

# load data
prices<-read.csv("SP500 securities_up_20230306.csv")
ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))

#return
return<- Return.calculate(ZOO, method="log")
set.seed(42) # for reproducibility
p = 3
stock_select =sample(1:454, p)
return<- return[-1, stock_select]
returnstd<-xts(return)

CORR = cor(returnstd[,1:p])
# Check positive-definiteness of CORR
eigenvalues = eig(CORR)
isPositiveDefinite = all(eigenvalues > 0)
# So CORR is positive-definite

# Given data
A <- CORR - diag(1,p,p)
mu <- apply(returnstd, 2, mean)
sigma <- apply(returnstd, 2, sd)
Sigma <- diag(sigma) %*% CORR %*% diag(sigma)

network_assets = network.correlation(returnstd)
phi =linfun3_1(A-diag(max(eigen(A)$value),p),
                 rep(0,p),
                 lambda=0.1,
                 abs(eigen(A-diag(max(eigen(A)$value),p))$vector[1,1])
)
phi=phi/max(phi)
# phi <- eigen_centrality(network_assets, directed = FALSE, scale = TRUE)$vector

phi_star = 0.96
mu_star = mean(mu)

###### PLug-In 2 constraints network portfolio #####
net.gmin.port = network.efficient.portfolio(phi, Sigma, phi_star,TRUE)
w =net.gmin.port$weights
w_network_2constraint_plugin<-w

###### PLug-In 3 constraints network portfolio #####
net.gmin.port = network.3constraint.portfolio(phi, mu, Sigma, phi_star, mu_star, TRUE)
w =net.gmin.port$weights
w_network_3constraint_plugin<-w

# Number of observations to generate
n_simu <- 10000

# Generate multivariate normal data
set.seed(42) # for reproducibility
data <- mvrnorm(n = n_simu, mu = mu, 
                Sigma = Sigma)

# Estimate mean vector
estimated_mu <- colMeans(data)

# Estimate correlation matrix
estimated_A <- cor(data)

# Estimate covariance matrix
estimated_Sigma <- cov(data)

# Error between True covariance Sigma and estimated covariance estimated_A
# Using Frobenius norm
error_norm <- norm(estimated_Sigma - Sigma, type = "F")
cat("Frobenius norm of the error:", error_norm, "\n")
# Relative error
relative_error_Sigma <- error_norm / norm(Sigma, type = "F")
cat("Relative error:", relative_error_Sigma, "\n")

error_norm <- norm(inv(estimated_Sigma) - inv(Sigma), type = "F")
cat("Frobenius norm of the error:", error_norm, "\n")
relative_error_Sigma_inv <- error_norm / norm(inv(estimated_Sigma), type = "F")
cat("Relative error:", relative_error_Sigma_inv, "\n")

# Estimate Eigenvector centrality
# EC_DS =linfun3_1(estimated_A-diag(max(eigen(estimated_A)$value),p),
#                  rep(0,p),
#                  lambda=0.1,
#                  abs(eigen(estimated_A--diag(max(eigen(estimated_A)$value),p))$vector[1,1])
# )
# EC_DS=EC_DS/max(EC_DS)
EC_DS = phi

##### Dantzig #######
# # lmd1 = estimate_dantzig_lambda1_portfolio_parameter_simplified(data, window_size=500)
# # lmd2 = estimate_dantzig_lambda2_portfolio_parameter_simplified(data, window_size=500)
# # lmd3 = estimate_dantzig_lambda3_portfolio_parameter_simplified(data, window_size=500)
# lmd1 = 0.3358
# lmd2 = 0.0001
# lmd3 = 0.1
# theta1 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
# theta2 =linfun1(estimated_Sigma,estimated_mu,lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
# theta3 =linfun1(estimated_Sigma,EC_DS,lambda=lmd3) # lambda <= 0.1 will lead to be infeasible


lmd1 = 0.001
theta1 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd1)
relative_error_SigmaInv_theta1 = relative_error(rep(1,p), estimated_Sigma %*% theta1)
relative_error_theta1 = relative_error(inv(Sigma) %*% rep(1,p), theta1)
cat("Relative error of Sigma * theta1:", relative_error_SigmaInv_theta1, "\n")
cat("Relative error of theta1:", relative_error_theta1, "\n")

lmd2 = 0.00001
theta2 =linfun1(estimated_Sigma,estimated_mu,lambda=lmd2)
relative_error_SigmaInv_theta2 = relative_error(estimated_mu, estimated_Sigma %*% theta2)
relative_error_theta2 = relative_error(inv(Sigma) %*% estimated_mu, theta2)
cat("Relative error of Sigma * theta2:", relative_error_SigmaInv_theta2, "\n")
cat("Relative error of theta2:", relative_error_theta2, "\n")

lmd3 = 0.001
theta3 =linfun1(estimated_Sigma,EC_DS,lambda=lmd3)
relative_error_SigmaInv_theta3 = relative_error(EC_DS, estimated_Sigma %*% theta3)
relative_error_theta3 = relative_error(inv(Sigma) %*% EC_DS, theta3)
cat("Relative error of Sigma * theta3:", relative_error_SigmaInv_theta3, "\n")
cat("Relative error of theta3:", relative_error_theta3, "\n")

# EC_DS
# lmd1=0.1
# 
# relative_error(estimated_mu, estimated_Sigma %*% theta2)
# 
# relative_error(EC_DS, estimated_Sigma %*% theta3)
# 
# estimated_Sigma %*% theta2
# 
# estimated_Sigma %*% theta3

###### Dantzig 2 constraints network portfolio #####
if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
  # gamma = (sum(theta1) + 1)/sum(theta3)
  # w = gamma * theta3 - theta1
  # sum(w)
  
  # gamma1 = (sum((inv(Sigma)) %*% rep(1,p)) + 1)/sum((inv(Sigma)) %*% EC_DS) 
  # w1 = gamma1 * (inv(Sigma)) %*% EC_DS - (inv(Sigma)) %*% rep(1,p)
  # sum(w1)
  
  alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
  w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
} else {
  w = theta1/sum(theta1)
  
  w1 = inv(Sigma) %*% rep(1,p) / sum((inv(Sigma)) %*% rep(1,p))
}
w[is.na(w)] <- 0
w_2constraint_dantzig<-w

relative_error(w_network_2constraint_plugin, w_2constraint_dantzig)
relative_error(w_network_2constraint_plugin, w1)
relative_error(w_2constraint_dantzig, w1)


###### Dantzig 3 constraints network portfolio #####
if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
  if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
    # 3 constraint case
    M1 <- cbind(rbind(1,phi_star,mean(estimated_mu)),
                rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
                rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
    M2 <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
                rbind(1,phi_star,mean(estimated_mu)),
                rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
    M <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
               rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
               rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
    gamma1 <- det(M1)/det(M)
    gamma2 <- det(M2)/det(M)
    alpha1 <- c(gamma1*sum(theta3))
    alpha2 <- c(gamma2*sum(theta2))
    w1 <- (1-alpha1-alpha2)*theta1/sum(theta1)
    w2 <- alpha1*theta3/sum(theta3)
    w3 <- alpha2*theta2/sum(theta2)
    w = w1 +  w2 + w3
  }
  else{
    # 1 centrality constraint case
    alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
    w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
  }} else{
    if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
      # 1 expected return constraint case
      alpha=c((sum(theta2)*sum(theta1)*mean(estimated_mu)-(sum(theta2))^2)/(estimated_mu%*%theta2*sum(theta1)-(sum(theta2))^2))
      w = alpha*theta2/sum(theta2)+(1-alpha)*theta1/sum(theta1)
    }
    else{
      # global minimum variance case
      w = theta1/sum(theta1)
    }
  }
w[is.na(w)] <- 0
w_3constraint_dantzig<-w

l2_norm_reletive_error_2constraint = relative_error(w_network_2constraint_plugin, w_2constraint_dantzig)
l2_norm_reletive_error_3constraint = relative_error(w_network_3constraint_plugin, w_3constraint_dantzig)
cat("L2 norm reletive error for 2-constraint portfolios:", l2_norm_reletive_error_2constraint, "\n")
cat("L2 norm reletive error for 3-constraint portfolios:", l2_norm_reletive_error_3constraint, "\n")

# ----- Simulation over different theta values -----
# Create an empty data frame to store the errors for each scenario
results <- data.frame(
  lmd1 = numeric(),
  lmd2 = numeric(),
  lmd3 = numeric(),
  error_2constraint = numeric(),
  error_3constraint = numeric()
)
num_lambda = 20
for (i1 in 1:num_lambda) {
  for (i2 in 1:num_lambda) {
    for (i3 in 1:num_lambda) {
      lmd1 = lmd1.list[i1]
      lmd2 = lmd2.list[i2]
      lmd3 = lmd3.list[i3]
      
      theta1 = theta1.list[[i1]]
      theta2 = theta2.list[[i2]]
      theta3 = theta3.list[[i3]]
      w = weight_2constraint_3constraint_network_port(theta1, theta2, theta3, phi_star, mu_star)
      
      l2_norm_reletive_error_2constraint = relative_error(w_network_2constraint_plugin, w$w_2constraint_dantzig)
      l2_norm_reletive_error_3constraint = relative_error(w_network_3constraint_plugin, w$w_3constraint_dantzig)
      
      # Store the results for this scenario
      results <- rbind(results, data.frame(
        lmd1 = lmd1,
        lmd2 = lmd2,
        lmd3 = lmd3,
        error_2constraint = l2_norm_reletive_error_2constraint,
        error_3constraint = l2_norm_reletive_error_3constraint
      ))
    }
  }
}

#### First case #####

# Install and load plotly if you haven't already
# install.packages("plotly")
library(plotly)
library(htmlwidgets)

# Create a 3D scatter plot for error_2constraint
fig <- plot_ly(
  data = results,
  x = ~lmd1,
  y = ~lmd2,
  z = ~lmd3,
  color = ~error_2constraint,   
  colors = colorRamp(c("blue", "red")),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
)
fig  # This will open an interactive window in your viewer/browser.
# Save the Plotly figure as an HTML file
saveWidget(as_widget(fig), "plotly_3d_scatter_2constraint.html")

# Create a 3D scatter plot for error_2constraint
fig <- plot_ly(
  data = results,
  x = ~lmd1,
  y = ~lmd2,
  z = ~lmd3,
  color = ~error_3constraint,   
  colors = colorRamp(c("blue", "red")),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
)
fig  # This will open an interactive window in your viewer/browser.
# Save the Plotly figure as an HTML file
saveWidget(as_widget(fig), "plotly_3d_scatter_3constraint.html")


# # Compare the 2-constraint portfolios
# l2_norm_2constraint <- l2_norm(w_2constraint_dantzig-w_network_2constraint_plugin)
# cat("L2 norm difference for 2-constraint portfolios:", l2_norm_2constraint, "\n")
# cat("L2 norm relative difference for 2-constraint portfolios:", l2_norm_2constraint/l2_norm(w_network_2constraint_plugin), "\n")
# 
# # Compare the 3-constraint portfolios
# l2_norm_3constraint <- l2_norm(w_3constraint_dantzig-w_network_3constraint_plugin)
# cat("L2 norm difference for 3-constraint portfolios:", l2_norm_3constraint, "\n")
# cat("L2 norm relative difference for 3-constraint portfolios:", l2_norm_3constraint/l2_norm(w_network_3constraint_plugin), "\n")
# 
# weight_lambdas = function(lmd1 = 0.3358, lmd2 = 0.0001, lmd3 = 0.1){
#   
#   theta1 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
#   theta2 =linfun1(estimated_Sigma,estimated_mu,lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
#   theta3 =linfun1(estimated_Sigma,EC_DS,lambda=lmd3) # lambda <= 0.1 will lead to be infeasible
#   
#   ###### Dantzig 2 constraints network portfolio #####
#   if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
#     alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
#     w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
#   } else{
#     w = theta1/sum(theta1)
#   }
#   w[is.na(w)] <- 0
#   w_2constraint_dantzig<-w
#   
#   ###### Dantzig 3 constraints network portfolio #####
#   if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
#     if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
#       # 3 constraint case
#       M1 <- cbind(rbind(1,phi_star,mean(estimated_mu)),
#                   rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
#                   rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
#       M2 <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
#                   rbind(1,phi_star,mean(estimated_mu)),
#                   rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
#       M <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
#                  rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
#                  rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
#       gamma1 <- det(M1)/det(M)
#       gamma2 <- det(M2)/det(M)
#       alpha1 <- c(gamma1*sum(theta3))
#       alpha2 <- c(gamma2*sum(theta2))
#       w1 <- (1-alpha1-alpha2)*theta1/sum(theta1)
#       w2 <- alpha1*theta3/sum(theta3)
#       w3 <- alpha2*theta2/sum(theta2)
#       w = w1 +  w2 + w3
#     }
#     else{
#       # 1 centrality constraint case
#       alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
#       w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
#     }} else{
#       if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
#         # 1 expected return constraint case 
#         alpha=c((sum(theta2)*sum(theta1)*mean(estimated_mu)-(sum(theta2))^2)/(estimated_mu%*%theta2*sum(theta1)-(sum(theta2))^2))
#         w = alpha*theta2/sum(theta2)+(1-alpha)*theta1/sum(theta1)
#       }
#       else{
#         # global minimum variance case 
#         w = theta1/sum(theta1)
#       }
#     }
#   w[is.na(w)] <- 0
#   w_3constraint_dantzig<-w
#   
#   return(list(e))
# }
# 
# 
# 
# sensitivity_Monte_Carlo = function(){
#   data <- mvrnorm(n = n_simu, mu = mu, 
#                   Sigma = Sigma)
#   
#   # Estimate mean vector
#   estimated_mu <- colMeans(data)
#   
#   # Estimate correlation matrix
#   estimated_A <- cor(data)
#   
#   # Estimate covariance matrix
#   estimated_Sigma <- cov(data)
#   
#   # Estimate Eigenvector centrality
#   EC_DS =linfun3_1(estimated_A-diag(1,p)-diag(max(eigen(estimated_A)$value),p),
#                    rep(0,p),
#                    lambda=0.45,
#                    abs(eigen(estimated_A-diag(1,p)-diag(max(eigen(estimated_A)$value),p))$vector[1,1])
#   )
#   EC_DS=EC_DS/max(EC_DS)
#   
#   ##### Dantzig #######
#   # lmd1 = estimate_dantzig_lambda1_portfolio_parameter_simplified(data, window_size=500)
#   lmd1 = 0.3358
#   lmd2 = 0.0001
#   lmd3 = 0.1
#   theta1 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
#   theta2 =linfun1(estimated_Sigma,estimated_mu,lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
#   theta3 =linfun1(estimated_Sigma,EC_DS,lambda=lmd3) # lambda <= 0.1 will lead to be infeasible
#   
#   ###### Dantzig 2 constraints network portfolio #####
#   if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
#     alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
#     w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
#   } else{
#     w = theta1/sum(theta1)
#   }
#   w[is.na(w)] <- 0
#   w_2constraint_dantzig<-w
#   
#   ###### Dantzig 3 constraints network portfolio #####
#   if((c(theta1/sum(theta1))%*%c(EC_DS))>phi_star){
#     if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
#       # 3 constraint case
#       M1 <- cbind(rbind(1,phi_star,mean(estimated_mu)),
#                   rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
#                   rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
#       M2 <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
#                   rbind(1,phi_star,mean(estimated_mu)),
#                   rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
#       M <- cbind(rbind(sum(theta1),EC_DS%*%theta1,estimated_mu%*%theta1),
#                  rbind(sum(theta3),EC_DS%*%theta3,estimated_mu%*%theta3),
#                  rbind(sum(theta2),EC_DS%*%theta2,estimated_mu%*%theta2))
#       gamma1 <- det(M1)/det(M)
#       gamma2 <- det(M2)/det(M)
#       alpha1 <- c(gamma1*sum(theta3))
#       alpha2 <- c(gamma2*sum(theta2))
#       w1 <- (1-alpha1-alpha2)*theta1/sum(theta1)
#       w2 <- alpha1*theta3/sum(theta3)
#       w3 <- alpha2*theta2/sum(theta2)
#       w = w1 +  w2 + w3
#     }
#     else{
#       # 1 centrality constraint case
#       alpha=c((sum(theta3)*sum(theta1)*phi_star-(sum(theta3))^2)/(EC_DS%*%theta3*sum(theta1)-(sum(theta3))^2))
#       w = alpha*theta3/sum(theta3)+(1-alpha)*theta1/sum(theta1)
#     }} else{
#       if((c(theta1/sum(theta1))%*%c(estimated_mu))<mean(estimated_mu)){
#         # 1 expected return constraint case 
#         alpha=c((sum(theta2)*sum(theta1)*mean(estimated_mu)-(sum(theta2))^2)/(estimated_mu%*%theta2*sum(theta1)-(sum(theta2))^2))
#         w = alpha*theta2/sum(theta2)+(1-alpha)*theta1/sum(theta1)
#       }
#       else{
#         # global minimum variance case 
#         w = theta1/sum(theta1)
#       }
#     }
#   w[is.na(w)] <- 0
#   w_3constraint_dantzig<-w
#   
#   # Function to compute the L2 (Euclidean) norm difference between two vectors
#   l2_norm <- function(x) {
#     sqrt(sum((x)^2))
#   }
#   relative_error <- function(x,y){
#     z=(y-x)/x
#     return(z)
#   }
#   
#   l2_norm_reletive_error_2constraint = l2_norm(relative_error(w_network_2constraint_plugin, w_2constraint_dantzig))
#   l2_norm_reletive_error_3constraint = l2_norm(relative_error(w_network_3constraint_plugin, w_3constraint_dantzig))
#   cat("L2 norm reletive error for 2-constraint portfolios:", l2_norm_reletive_error_2constraint, "\n")
#   cat("L2 norm reletive error for 3-constraint portfolios:", l2_norm_reletive_error_3constraint, "\n")
#   
#   # # Compare the 2-constraint portfolios
#   # l2_norm_2constraint <- l2_norm(w_2constraint_dantzig, w_network_2constraint_plugin)
#   # cat("L2 norm difference for 2-constraint portfolios:", l2_norm_2constraint, "\n")
#   # cat("L2 norm relative difference for 2-constraint portfolios:", l2_norm_2constraint/sqrt(sum(w_network_2constraint_plugin^2)), "\n")
#   # 
#   # # Compare the 3-constraint portfolios
#   # l2_norm_3constraint <- l2_norm(w_3constraint_dantzig, w_network_3constraint_plugin)
#   # cat("L2 norm difference for 3-constraint portfolios:", l2_norm_3constraint, "\n")
#   # cat("L2 norm relative difference for 3-constraint portfolios:", l2_norm_3constraint/sqrt(sum(w_network_3constraint_plugin^2)), "\n")
#   
#   return(list(
#     error_2constraint = l2_norm_reletive_error_2constraint,
#     error_3constraint = l2_norm_reletive_error_3constraint
#   ))
# }
# 
# 
# # Number of bootstrap replications
# B <- 10  
# l2_errors_2constraint <- numeric(B)
# l2_errors_3constraint <- numeric(B)
# 
# set.seed(42)
# for (b in 1:B) {
#   errors = sensitivity_Monte_Carlo()
#   
#   # Compute L2 error for the 2-constraint portfolios
#   l2_errors_2constraint[b] <- errors$error_2constraint
#   
#   # Compute L2 error for the 3-constraint portfolios
#   l2_errors_3constraint[b] <- errors$error_3constraint
# }
# 
# # Now compute summary statistics:
# error_mean_2constraint <- mean(l2_errors_2constraint)
# error_sd_2constraint <- sd(l2_errors_2constraint)
# error_ci_2constraint <- quantile(l2_errors_2constraint, c(0.025, 0.975))
# 
# cat("Bootstrap mean L2 error (2-constraint):", error_mean_2constraint, "\n")
# cat("Bootstrap error 95% CI:", error_ci_2constraint, "\n")
# 
# # Now compute summary statistics:
# error_mean_3constraint <- mean(l2_errors_3constraint)
# error_sd_3constraint <- sd(l2_errors_3constraint)
# error_ci_3constraint <- quantile(l2_errors_3constraint, c(0.025, 0.975))
# 
# cat("Bootstrap mean L2 error (3-constraint):", error_mean_3constraint, "\n")
# cat("Bootstrap error 95% CI:", error_ci_3constraint, "\n")