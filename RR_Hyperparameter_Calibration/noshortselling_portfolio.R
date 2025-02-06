# 任意选三只股票，求GMV,MV,EW,Plug-In,MC

# Clear the workspace 
rm(list = ls())

# Set the working directory 
setwd("~/Documents/GitHub/Network-Portfolio-Review-Response/RR_Hyperparameter_Calibration")

n_samples = 10000
stock_select = c(1,3,6)
# Number of observations to generate
n <- 500
n_outsample <- 250

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

p=3

MonteCarlo_network_2constraint_noshort_portfolio = function(n_samples = 10000, 
                                                    Sigma,
                                                    phi, phi_star){
  # Function to generate a feasible weight vector
  generate_feasible_w <- function() {
    repeat {
      w <- runif(3)
      w <- w / sum(w)
      if (sum(w) == 1 && 
          phi %*% w <= phi_star &&
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
  return(best_w)
}

MonteCarlo_network_3constraint_noshort_portfolio = function(n_samples = 10000, 
                                                    Sigma,
                                                    phi, phi_star,
                                                    mu, mu_star){
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
  return(best_w)
}

# Generate data

data = returnstd[1:n, stock_select]
data_out = returnstd[(n+1):(n+n_outsample), stock_select]

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
                 lambda=0.45,
                 abs(eigen(estimated_A-diag(1,3,3)-diag(max(eigen(estimated_A)$value),3,3))$vector[1,1])
)
EC_DS=EC_DS/max(EC_DS)

phi_star = mean(EC_DS)

###### minimum variance portfolio  #####
portf_minVar =globalMin.portfolio(estimated_mu,estimated_Sigma)
w =portf_minVar$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_minVar<-rowSums(aus)+1
cumureturn_minVar<-cumprod(return_minVar)
w_minVar<-w

###### mean variance portfolio  ######
portf_meanVar =efficient.portfolio(estimated_mu,estimated_Sigma,mean(estimated_mu))
w =portf_meanVar$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_meanVar<-rowSums(aus)+1
cumureturn_meanVar<-cumprod(return_meanVar)
w_meanVar<-w

###### equally weighted portfolio #####
w =matrix(1/p,1,p)
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_equal<-rowSums(aus)+1
cumureturn_equal<-cumprod(return_equal)
w_equal<-w

###### PLug-In 2 constraints network portfolio #####
net.gmin.port = network.efficient.portfolio(EC_DS, estimated_Sigma, phi_star,TRUE)
w =net.gmin.port$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_network_vary_with_phi<-rowSums(aus)+1
cumureturn_network_vary_with_phi<-cumprod(return_network_vary_with_phi)
w_network_vary_with_phi<-w

###### PLug-In 3 constraints network portfolio #####
net.gmin.port = network.3constraint.portfolio(EC_DS, estimated_mu, estimated_Sigma, phi_star, mu_star, TRUE)
w =net.gmin.port$weights
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_network_3constraint_plugin<-rowSums(aus)+1
cumureturn_network_3constraint_plugin<-cumprod(return_network_3constraint_plugin)
w_network_3constraint_plugin<-w

###### Monte Carlo 2 constraints portfolio #####
w = MonteCarlo_network_2constraint_noshort_portfolio(n_samples, Sigma, EC_DS, phi_star)
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_network_MC_2constraint<-rowSums(aus)+1
cumureturn_network_MC_2constraint<-cumprod(return_network_MC_2constraint)
w_network_MC_2constraint<-w

###### Monte Carlo 3 constraints portfolio #####
w = MonteCarlo_network_3constraint_noshort_portfolio(n_samples, Sigma, EC_DS, phi_star, estimated_mu, mu_star)
aus<-as.matrix(repmat(w,n_outsample,1)*data_out)
return_network_MC_3constraint<-rowSums(aus)+1
cumureturn_network_MC_3constraint<-cumprod(return_network_MC_3constraint)
w_network_MC_3constraint<-w

##### Compare results ######

tb = rbind(cbind(cumureturn_minVar[n_outsample]/n_outsample*252*100,
                 cumureturn_meanVar[n_outsample]/n_outsample*252*100,
                 cumureturn_equal[n_outsample]/n_outsample*252*100,
                 cumureturn_network_vary_with_phi[n_outsample]/n_outsample*252*100,
                 cumureturn_network_3constraint_plugin[n_outsample]/n_outsample*252*100,
                 cumureturn_network_MC_2constraint[n_outsample]/n_outsample*252*100,
                 cumureturn_network_MC_3constraint[n_outsample]/n_outsample*252*100),
           cbind(std(return_minVar)*sqrt(252)*100,
                 std(return_meanVar)*sqrt(252)*100,
                 std(return_equal)*sqrt(252)*100,
                 std(return_network_vary_with_phi)*sqrt(252)*100,
                 std(return_network_3constraint_plugin)*sqrt(252)*100,
                 std(return_network_MC_2constraint)*sqrt(252)*100,
                 std(return_network_MC_3constraint)*sqrt(252)*100),
           cbind(cumureturn_minVar[n_outsample]/n_outsample*252/(std(return_minVar)*sqrt(252))*100,
                 cumureturn_meanVar[n_outsample]/n_outsample*252/(std(return_meanVar)*sqrt(252))*100,
                 cumureturn_equal[n_outsample]/n_outsample*252/(std(return_equal)*sqrt(252))*100,
                 cumureturn_network_vary_with_phi[n_outsample]/n_outsample*252/(std(return_network_vary_with_phi)*sqrt(252))*100,
                 cumureturn_network_3constraint_plugin[n_outsample]/n_outsample*252/(std(return_network_3constraint_plugin)*sqrt(252))*100,
                 cumureturn_network_MC_2constraint[n_outsample]/n_outsample*252/(std(return_network_MC_2constraint)*sqrt(252))*100,
                 cumureturn_network_MC_3constraint[n_outsample]/n_outsample*252/(std(return_network_MC_3constraint)*sqrt(252))*100),
           cbind(maxDrawdown(return_minVar-1),
                 maxDrawdown(return_meanVar-1),
                 maxDrawdown(return_equal-1),
                 maxDrawdown(return_network_vary_with_phi-1),
                 maxDrawdown(return_network_3constraint_plugin-1),
                 maxDrawdown(return_network_MC_2constraint-1),
                 maxDrawdown(return_network_MC_3constraint-1))
)
rownames(tb) = c()
xtable(tb,digits = 2)

