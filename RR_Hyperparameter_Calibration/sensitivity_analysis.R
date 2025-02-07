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

# load data
prices<-read.csv("SP500 securities_up_20230306.csv")
ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))

#return
return<- Return.calculate(ZOO, method="log")
return<- return[-1, 1:10]
returnstd<-xts(return)
# p=dim(return)[2]
p=dim(return)[2]

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
phi <- eigen_centrality(network_assets,directed = FALSE, scale = TRUE)$vector

phi_star = quantile(phi, 0.1)
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
n_simu <- 1000000

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
relative_error <- error_norm / norm(Sigma, type = "F")
cat("Relative error:", relative_error, "\n")

error_norm <- norm(inv(estimated_Sigma) - inv(Sigma), type = "F")
cat("Frobenius norm of the error:", error_norm, "\n")
relative_error <- error_norm / norm(Sigma, type = "F")
cat("Relative error:", relative_error, "\n")

# Estimate Eigenvector centrality
EC_DS =linfun3_1(estimated_A-diag(1,p)-diag(max(eigen(estimated_A)$value),p),
                 rep(0,p),
                 lambda=0.45,
                 abs(eigen(estimated_A-diag(1,p)-diag(max(eigen(estimated_A)$value),p))$vector[1,1])
)
EC_DS=EC_DS/max(EC_DS)

##### Dantzig #######
# lmd1 = estimate_dantzig_lambda1_portfolio_parameter_simplified(data, window_size=500)
# lmd2 = estimate_dantzig_lambda2_portfolio_parameter_simplified(data, window_size=500)
# lmd3 = estimate_dantzig_lambda3_portfolio_parameter_simplified(data, window_size=500)
lmd1 = 0.3358
lmd2 = 0.0001
lmd3 = 0.1
theta1 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
theta2 =linfun1(estimated_Sigma,estimated_mu,lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
theta3 =linfun1(estimated_Sigma,EC_DS,lambda=lmd3) # lambda <= 0.1 will lead to be infeasible


relative_error_list1 = c()
lmd1.list = seq(0.01,0.9,length=100)
l.lmd1 = length(lmd1.list)
for(i in 1:l.lmd1){
  lmd1 = lmd1.list[i]
  theta1 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd1)
  relative_error_list1[i] = relative_error(rep(1,p), estimated_Sigma %*% theta1)
}

relative_error_list2 = c()
lmd2.list = seq(0.000001,0.0001,length=100)
l.lmd2 = length(lmd2.list)
for(i in 1:l.lmd2){
  lmd2 = lmd2.list[i]
  theta2 =linfun1(estimated_Sigma,rep(1,p),lambda=lmd2)
  relative_error_list2[i] = relative_error(rep(1,p), estimated_Sigma %*% theta2)
}

lmd1=0.1

relative_error(estimated_mu, estimated_Sigma %*% theta2)

relative_error(EC_DS, estimated_Sigma %*% theta3)

estimated_Sigma %*% theta2

estimated_Sigma %*% theta3

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



l2_norm_reletive_error_2constraint = relative_error(w_network_2constraint_plugin, w_2constraint_dantzig)
l2_norm_reletive_error_3constraint = relative_error(w_network_3constraint_plugin, w_2constraint_dantzig)
cat("L2 norm reletive error for 2-constraint portfolios:", l2_norm_reletive_error_2constraint, "\n")
cat("L2 norm reletive error for 3-constraint portfolios:", l2_norm_reletive_error_3constraint, "\n")

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