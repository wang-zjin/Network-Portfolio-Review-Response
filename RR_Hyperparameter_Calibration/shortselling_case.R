# Clear the workspace 
rm(list = ls())

# Set the working directory 
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

# load data
prices<-read.csv("SP500 securities_up_20230306.csv")
ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))

#return
return<- Return.calculate(ZOO, method="log")
return<- return[-1, ]
returnstd<-xts(return)
p=dim(return)[2]

mean(returnstd[,1])
mean(returnstd[,3])
mean(returnstd[,6])
std(returnstd[,1])
std(returnstd[,3])
std(returnstd[,6])
cor(returnstd[,c(1,3,6)])

# Given data
A <- matrix(c(1, 0.3075, 0.3465, 0.3075, 1, 0.2993, 0.3465, 0.2993, 1), nrow = 3, ncol = 3)
mu <- c(0.00068, -0.00005, 0.000230)
sigma <- c(0.0168, 0.0212, 0.0178)
Sigma <- diag(sigma) %*% A %*% diag(sigma)
phi <- eigen(A)$vectors[,1] / max(eigen(A)$vectors[,1])
phi_star <- 0.98
mu_star <- 0

# Number of Monte Carlo samples
n_samples <- 10000

# Function to generate a feasible weight vector
generate_feasible_w <- function() {
  repeat {
    w <- runif(3)
    w <- w / sum(w)
    if (sum(w) == 1 && 
        phi %*% w <= phi_star &&
        mu %*% w >= mu_star &&
        all(w >= 0)) {
      return(w)
    }
  }
}

# Monte Carlo simulation
best_w <- NULL
best_obj <- Inf

for (i in 1:n_samples) {
  w <- generate_feasible_w()
  obj <- 0.5 * t(w) %*% Sigma %*% w
  if (obj < best_obj) {
    best_obj <- obj
    best_w <- w
  }
}

# Output the best weights and their corresponding objective value
cat("Best weights:\n", best_w, "\n")
cat("Best objective value:\n", best_obj, "\n")

# Generate observations based on the given parameters

# Number of observations to generate
n <- 500
n_outsample <- 250

# Generate multivariate normal data
set.seed(42) # for reproducibility
data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
data_out <- mvrnorm(n = n_outsample, mu = mu, Sigma = Sigma)

# Estimate mean vector
estimated_mu <- colMeans(data)

# Estimate correlation matrix
estimated_A <- cor(data)

# Estimate covariance matrix
estimated_Sigma <- cov(data)

# Print results
cat("Estimated mean vector:\n", estimated_mu, "\n")
cat("Estimated correlation matrix:\n")
print(estimated_A)

# Estimate Eigenvector centrality
EC_DS =linfun3_1(estimated_A-diag(1,3,3)-diag(max(eigen(estimated_A)$value),3,3),
                      rep(0,3),
                      lambda=0.618,
                      abs(eigen(estimated_A-diag(1,3,3)-diag(max(eigen(estimated_A)$value),3,3))$vector[1,1])
)
EC_DS=EC_DS/max(EC_DS)

# PLug-In #####
net.gmin.port = network.3constraint.portfolio(EC_DS, estimated_mu, estimated_Sigma, phi_star, mu_star, TRUE)
w =net.gmin.port$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_3constraint_plugin<-rowSums(aus)
cumureturn_3constraint_plugin<-cumprod(return_3constraint_plugin+1)

# Dantzig #####
lmd1 = 0.4636
lmd2 = 0.0001
lmd3 = 0.1
theta1 =linfun1(estimated_Sigma,rep(1,3),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
theta2 =linfun1(estimated_Sigma,estimated_mu,lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
theta3 =linfun1(estimated_Sigma,EC_DS,lambda=lmd3) # lambda <= 0.1 will lead to be infeasible

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
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_3constraint_dantzig<-rowSums(aus)
cumureturn_3constraint_dantzig<-cumprod(return_3constraint_dantzig+1)

# Glasso #####
glasso.icov=glasso(estimated_Sigma,rho=0.0001)$wi
if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS))>phi_star){
  if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(estimated_mu))<mean(estimated_mu)){
    # 3 constraint case
    M1 <- cbind(rbind(1,phi_star,mean(estimated_mu)),
                rbind(sum(glasso.icov%*%EC_DS),EC_DS%*%glasso.icov%*%EC_DS,estimated_mu%*%glasso.icov%*%EC_DS),
                rbind(sum(glasso.icov%*%estimated_mu),EC_DS%*%glasso.icov%*%estimated_mu,estimated_mu%*%glasso.icov%*%estimated_mu))
    M2 <- cbind(rbind(sum(glasso.icov),EC_DS%*%row_sums(glasso.icov),estimated_mu%*%row_sums(glasso.icov)),
                rbind(1,phi_star,mean(estimated_mu)),
                rbind(sum(glasso.icov%*%estimated_mu),EC_DS%*%glasso.icov%*%estimated_mu,estimated_mu%*%glasso.icov%*%estimated_mu))
    M <- cbind(rbind(sum(glasso.icov),EC_DS%*%theta1,estimated_mu%*%theta1),
               rbind(sum(glasso.icov%*%EC_DS),EC_DS%*%glasso.icov%*%EC_DS,estimated_mu%*%glasso.icov%*%EC_DS),
               rbind(sum(glasso.icov%*%estimated_mu),EC_DS%*%glasso.icov%*%estimated_mu,estimated_mu%*%glasso.icov%*%estimated_mu))
    gamma1 <- det(M1)/det(M)
    gamma2 <- det(M2)/det(M)
    alpha1 <- c(gamma1*sum(glasso.icov%*%EC_DS))
    alpha2 <- c(gamma2*sum(glasso.icov%*%estimated_mu))
    w1 <- (1-alpha1-alpha2)*row_sums(glasso.icov)/sum(glasso.icov)
    w2 <- alpha1*glasso.icov%*%EC_DS/sum(glasso.icov%*%EC_DS)
    w3 <- alpha2*glasso.icov%*%estimated_mu/sum(glasso.icov%*%estimated_mu)
    w = c(w1 +  w2 + w3)
  }
  else{
    # 1 centrality constraint case
    alpha=c((sum(glasso.icov%*%EC_DS)*sum(glasso.icov)*phi_star-(sum(glasso.icov%*%EC_DS))^2)/(EC_DS%*%glasso.icov%*%EC_DS*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS))^2))
    w = c(alpha[1]*glasso.icov%*%EC_DS/sum(glasso.icov%*%EC_DS)+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
  }} else{
  if(((row_sums(glasso.icov)/sum(glasso.icov))%*%estimated_mu)<mean(estimated_mu)){
    # 1 expected return constraint case 
    alpha=c((sum(glasso.icov%*%estimated_mu)*sum(glasso.icov)*mean(estimated_mu)-(sum(glasso.icov%*%estimated_mu))^2)/(estimated_mu%*%glasso.icov%*%estimated_mu*sum(glasso.icov)-(sum(glasso.icov%*%estimated_mu))^2))
    w = c(alpha[1]*glasso.icov%*%estimated_mu/sum(glasso.icov%*%estimated_mu)+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
  }
  else{
    # global minimum variance case 
    w=row_sums(glasso.icov)/sum(glasso.icov)
  }
}
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_3constraint_glasso<-rowSums(aus)
cumureturn_3constraint_glasso<-cumprod(return_3constraint_glasso+1)
