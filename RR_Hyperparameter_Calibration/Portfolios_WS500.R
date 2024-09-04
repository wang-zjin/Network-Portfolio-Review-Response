rm(list = ls())

# setwd("~/Documents/GitHub/Network-Portfolio/Empirical_results/Monthly_rebalance")
setwd("~/Documents/GitHub/Network-Portfolio/RR_Hyperparameter_Calibration")

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

# set label
node.label=colnames(returnstd)
names(returnstd) = node.label

# rolling window
W<-list()
window_size = 500
for(t in 0: (floor((1857-window_size)/22)-1)){
  W[[(t+1)]]=returnstd[(1+t*22):(window_size+22+t*22),]
}
W_in<-list()
W_out<-list()
for(t in 0: (floor((1857-window_size)/22)-1)){
  W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
  W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
}
T.windows<-length(W)

#### Load Dantzig estimated parameters ####
load("eigenvector_centrality_WS500_20240723.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

load("theta_Dantzig_WS500_20240723.RData")
theta1=theta$theta1
theta2=theta$theta2
theta3=theta$theta3

load("Glasso_rho.RData")
rho=rho.WS500

#### Weights of different portfolios in each rolling window ####
##### minimum variance portfolio  #####

###### minimum variance portfolio with Dantzig estimation  ######
w<-list()
cumureturn_minVar_Dantzig<-list()
for(t in 1: length(W_in)){
  w[[t]]=theta1[[t]]/sum(theta1[[t]])
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_minVar_Dantzig[[t]]<-rowSums(aus)
}
return_minVar_Dantzig<-as.matrix(cbind(unlist(cumureturn_minVar_Dantzig))+1)
cumureturn_minVar_Dantzig<-cumprod(return_minVar_Dantzig)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_minVar_Dantzig<-w

###### minimum variance portfolio with plug in  ######
w<-list()
cumureturn_minVar<-list()
for(t in 1: length(W_in)){
  portf_minVar =globalMin.portfolio(ER_in[[(t)]],COV_in[[(t)]])
  w[[(t)]] =portf_minVar$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_minVar[[t]]<-rowSums(aus)
}
return_minVar<-as.matrix(cbind(unlist(cumureturn_minVar))+1)
cumureturn_minVar<-cumprod(return_minVar)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_minVar<-w

###### minimum variance portfolio with glasso #####
w<-list()
cumureturn<-list()
for(t in 1: length(W_in)){
  glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
  w[[t]]=row_sums(glasso.icov)/sum(glasso.icov)
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn[[t]]<-rowSums(aus)
}
return_minVar_glasso<-as.matrix(cbind(unlist(cumureturn))+1)
cumureturn_minVar_glasso<-cumprod(return_minVar_glasso)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_minVar_glasso<-w

###### minimum variance portfolio with plug in no short ######
w<-list()
cumureturn_minVar_noshort<-list()
for(t in 1: length(W_in)){
  portf_minVar =globalMin.portfolio(ER_in[[(t)]],COV_in[[(t)]],FALSE)
  w[[(t)]] =portf_minVar$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_minVar_noshort[[t]]<-rowSums(aus)
}
return_minVar_noshort<-as.matrix(cbind(unlist(cumureturn_minVar_noshort))+1)
cumureturn_minVar_noshort<-cumprod(return_minVar_noshort)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_minVar_noshort<-w

##### mean variance portfolio  #####
###### mean variance portfolio with plug in  ######
w<-list()
cumureturn_meanVar<-list()
for(t in 1: length(W_in)){
  portf_meanVar =efficient.portfolio(ER_in[[(t)]],COV_in[[(t)]],mean(ER_in[[t]]))
  w[[(t)]] =portf_meanVar$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_meanVar[[t]]<-rowSums(aus)
}
return_meanVar<-as.matrix(cbind(unlist(cumureturn_meanVar))+1)
cumureturn_meanVar<-cumprod(return_meanVar)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_meanVar<-w

###### mean variance portfolio with Dantzig estimation  ######
w<-list()
cumureturn_meanVar_Dantzig<-list()
for(t in 1: length(W_in)){
  alpha=(sum(theta2[[t]])*sum(theta1[[t]])*mean(ER_in[[t]])-(sum(theta2[[t]]))^2)/(ER_in[[(t)]]%*%theta2[[t]]*sum(theta1[[t]])-(sum(theta2[[t]]))^2)
  w[[(t)]] = alpha*theta2[[t]]/sum(theta2[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_meanVar_Dantzig[[t]]<-rowSums(aus)
}
return_meanVar_Dantzig<-as.matrix(cbind(unlist(cumureturn_meanVar_Dantzig))+1)
cumureturn_meanVar_Dantzig<-cumprod(return_meanVar_Dantzig)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_meanVar_Dantzig<-w

###### mean variance portfolio with glasso ######
w<-list()
cumureturn<-list()
for(t in 1: length(W_in)){
  glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
  alpha=(sum(glasso.icov%*%ER_in[[t]])*sum(glasso.icov)*mean(ER_in[[t]])-(sum(glasso.icov%*%ER_in[[t]]))^2)/(ER_in[[(t)]]%*%glasso.icov%*%ER_in[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%ER_in[[t]]))^2)
  w[[(t)]] = c(alpha[1]*glasso.icov%*%ER_in[[t]]/sum(glasso.icov%*%ER_in[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn[[t]]<-rowSums(aus)
}
return_meanVar_glasso<-as.matrix(cbind(unlist(cumureturn))+1)
cumureturn_meanVar_glasso<-cumprod(return_meanVar_glasso)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_meanVar_glasso<-w

###### mean variance portfolio with plug in no short ######
w<-list()
cumureturn_meanVar_noshort<-list()
for(t in 1: length(W_in)){
  portf_meanVar =efficient.portfolio(ER_in[[(t)]],COV_in[[(t)]],mean(ER_in[[t]]),shorts = FALSE)
  w[[(t)]] =portf_meanVar$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_meanVar_noshort[[t]]<-rowSums(aus)
}
return_meanVar_noshort<-as.matrix(cbind(unlist(cumureturn_meanVar_noshort))+1)
cumureturn_meanVar_noshort<-cumprod(return_meanVar_noshort)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_meanVar_noshort<-w
##### equally weighted portfolio #####
w<-list()
cumureturn_temporal<-list()
centrality_equal_portfolio<-list()
for(t in 1: length(W_in)){
  w[[(t)]] =matrix(1/p,1,p)
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
  centrality_equal_portfolio[[t]]<-as.double(w[[(t)]]%*%EC_in[[(t)]])
}
return_equal<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_equal<-cumprod(return_equal)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_equal<-w

##### network portfolio  constraint #####
### default setting centrality constraint as mean centrality ###
###### network portfolio fixed constraint as mean centrality with plug in ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio without short ##
  net.gmin.port = network.efficient.portfolio(EC_DS[[(t)]], COV_in[[(t)]],mean(EC_DS[[(t)]]))
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_2constraint<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_2constraint<-cumprod(return_network_2constraint)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_2constraint<-w

###### network portfolio fixed constraint as mean centrality with Dantzig ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio without short ##
  if(c(theta1[[t]]/sum(theta1[[t]]))%*%c(EC_DS[[(t)]])>mean(EC_DS[[t]])){
    alpha=(sum(theta3[[t]])*sum(theta1[[t]])*mean(EC_DS[[t]])-(sum(theta3[[t]]))^2)/(EC_DS[[(t)]]%*%theta3[[t]]*sum(theta1[[t]])-(sum(theta3[[t]]))^2)
    w[[(t)]] = c(alpha)*theta3[[t]]/sum(theta3[[t]])+(1-c(alpha))*theta1[[t]]/sum(theta1[[t]])
  }
  else{
    w[[(t)]] = theta1[[t]]/sum(theta1[[t]])
  }
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_2constraint_Dantzig<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_2constraint_Dantzig<-cumprod(return_network_2constraint_Dantzig)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_2constraint_Dantzig<-w

###### network portfolio fixed constraint as mean centrality with glasso ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio without short ##
  glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
  if(c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS[[t]])>mean(EC_DS[[t]])){
    alpha=c((sum(glasso.icov%*%EC_DS[[t]])*sum(glasso.icov)*mean(EC_DS[[t]])-(sum(glasso.icov%*%EC_DS[[t]]))^2)/(EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS[[t]]))^2))
    w[[(t)]] = c(alpha[1]*glasso.icov%*%EC_DS[[t]]/sum(glasso.icov%*%EC_DS[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
  }
  else{
    w[[(t)]] = row_sums(glasso.icov)/sum(glasso.icov)
  }
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_2constraint_glasso<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_2constraint_glasso<-cumprod(return_network_2constraint_glasso)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_2constraint_glasso<-w

###### network portfolio fixed constraint as mean centrality with plug in no short ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio without short ##
  net.gmin.port = network.efficient.portfolio(EC_DS[[(t)]], COV_in[[(t)]],mean(EC_DS[[(t)]]),FALSE)
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
  # print(t)
}
return_network_2constraint_noshort<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_2constraint_noshort<-cumprod(return_network_2constraint_noshort)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_2constraint_noshort<-w

###### network portfolio varying constraint with plug in ######
cumureturn_network_vary_with_phi<-list()
return_network_vary_with_phi<-list()
quantl<-seq(0.1,0.9,0.1)
w_network_vary_with_phi<-list()
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## compute 1-constraint network portfolio ##
    net.gmin.port = network.efficient.portfolio(EC_DS[[(t)]], COV_in[[(t)]],quantile(EC_DS[[(t)]],quantl[i]),TRUE)
    w[[(t)]] =net.gmin.port$weights
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi[[i]]<-cumprod(return_network_vary_with_phi[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi[[i]]<-w
}

###### network portfolio varying constraint with Dantzig ######
cumureturn_network_vary_with_phi_Dantzig<-list()
return_network_vary_with_phi_Dantzig<-list()
quantl<-seq(0.1,0.9,0.1)
w_network_vary_with_phi_Dantzig<-list()
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    if((c(theta1[[t]]/sum(theta1[[t]]))%*%c(EC_DS[[(t)]]))>quantile(EC_DS[[(t)]],quantl[i])){
      alpha=c((sum(theta3[[t]])*sum(theta1[[t]])*quantile(EC_DS[[(t)]],quantl[i])-(sum(theta3[[t]]))^2)/(EC_DS[[(t)]]%*%theta3[[t]]*sum(theta1[[t]])-(sum(theta3[[t]]))^2))
      w[[(t)]] = alpha*theta3[[t]]/sum(theta3[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
    }
    else{
      w[[(t)]] = theta1[[t]]/sum(theta1[[t]])
    }
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_Dantzig[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_Dantzig[[i]]<-cumprod(return_network_vary_with_phi_Dantzig[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_Dantzig[[i]]<-w
}

###### network portfolio varying constraint with glasso ######
cumureturn_network_vary_with_phi_glasso<-list()
return_network_vary_with_phi_glasso<-list()
quantl<-seq(0.1,0.9,0.1)
w_network_vary_with_phi_glasso<-list()
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## compute global minimum variance portfolio ##glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
    glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
    if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS[[t]]))>quantile(EC_DS[[(t)]],quantl[i])){
      alpha=c((sum(glasso.icov%*%EC_DS[[t]])*sum(glasso.icov)*quantile(EC_DS[[(t)]],quantl[i])-(sum(glasso.icov%*%EC_DS[[t]]))^2)/(EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS[[t]]))^2))
      w[[(t)]] = c(alpha[1]*glasso.icov%*%EC_DS[[t]]/sum(glasso.icov%*%EC_DS[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
    }
    else{
      w[[(t)]] = row_sums(glasso.icov)/sum(glasso.icov)
    }
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_glasso[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_glasso[[i]]<-cumprod(return_network_vary_with_phi_glasso[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_glasso[[i]]<-w
}

###### network portfolio varying constraint with plug in no short ######
cumureturn_network_vary_with_phi_noshort<-list()
return_network_vary_with_phi_noshort<-list()
w_network_vary_with_phi_noshort<-list()
quantl<-seq(0.1,0.9,0.1)
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## compute global minimum variance portfolio ##
    net.gmin.port = network.efficient.portfolio(EC_DS[[(t)]], COV_in[[(t)]],quantile(EC_DS[[(t)]],quantl[i]),FALSE)
    w[[(t)]] =net.gmin.port$weights
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_noshort[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_noshort[[i]]<-cumprod(return_network_vary_with_phi_noshort[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_noshort[[i]]<-w
}
###### network portfolio data-driven constraint with plug in ######
cumureturn_network_datadriven_phistar<-list()
return_network_datadriven_phistar<-list()
phi_star<-centrality_equal_portfolio
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  net.gmin.port = network.efficient.portfolio(EC_DS[[(t)]], COV_in[[(t)]],phi_star[[t]],TRUE)
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_datadriven_phistar<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar<-cumprod(return_network_datadriven_phistar)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar<-w

###### network portfolio data-driven constraint with Dantzig ######
cumureturn_network_datadriven_phistar_Dantzig<-list()
return_network_datadriven_phistar_Dantzig<-list()
phi_star<-centrality_equal_portfolio
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  if(c(theta1[[t]]/sum(theta1[[t]]))%*%c(EC_DS[[(t)]])>phi_star[[t]]){
    alpha=c((sum(theta3[[t]])*sum(theta1[[t]])*phi_star[[t]]-(sum(theta3[[t]]))^2)/(EC_DS[[(t)]]%*%theta3[[t]]*sum(theta1[[t]])-(sum(theta3[[t]]))^2))
    w[[(t)]] = alpha*theta3[[t]]/sum(theta3[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
  }
  else{
    w[[(t)]] = theta1[[t]]/sum(theta1[[t]])
  }
  aus <- as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]] <- rowSums(aus)
}
return_network_datadriven_phistar_Dantzig<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_Dantzig<-cumprod(return_network_datadriven_phistar_Dantzig)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_Dantzig<-w

###### network portfolio data-driven constraint with glasso ######
cumureturn_network_datadriven_phistar_glasso<-list()
return_network_datadriven_phistar_glasso<-list()
phi_star<-centrality_equal_portfolio
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
  if(c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS[[t]])>phi_star[[t]]){
    alpha=c((sum(glasso.icov%*%EC_DS[[t]])*sum(glasso.icov)*phi_star[[t]]-(sum(glasso.icov%*%EC_DS[[t]]))^2)/(EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS[[t]]))^2))
    w[[(t)]] = c(alpha[1]*glasso.icov%*%EC_DS[[t]]/sum(glasso.icov%*%EC_DS[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
  }
  else{
    w[[(t)]] = row_sums(glasso.icov)/sum(glasso.icov)
  }
  aus <- as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]] <- rowSums(aus)
}
return_network_datadriven_phistar_glasso<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_glasso<-cumprod(return_network_datadriven_phistar_glasso)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_glasso<-w

###### network portfolio data-driven constraint with plug in no short ######
cumureturn_network_datadriven_phistar_noshort<-list()
return_network_datadriven_phistar_noshort<-list()
cumureturn_temporal<-list()
phi_star<-centrality_equal_portfolio
w<-list()
for(t in 1: length(W_in)){
  ## compute network portfolio data-driven constraint no short ##
  net.gmin.port = network.efficient.portfolio(EC_DS[[(t)]], COV_in[[(t)]],phi_star[[t]],FALSE)
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_datadriven_phistar_noshort<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_noshort<-cumprod(return_network_datadriven_phistar_noshort)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_noshort<-w

##### network portfolio 3 constraint #####
### default setting centrality constraint as mean centrality ###

###### network portfolio fixed constraint with plug in ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute 3 constraint network portfolio ##
  net.gmin.port = network.3constraint.portfolio(EC_DS[[(t)]],ER_in[[(t)]], COV_in[[(t)]],mean(EC_DS[[(t)]]),mean(ER_in[[(t)]]))
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_3constraint<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_3constraint<-cumprod(return_network_3constraint)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_3constraint<-w

###### network portfolio fixed constraint with Dantzig ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute 3 constraint network portfolio ##
  if(c(theta1[[t]]/sum(theta1[[t]]))%*%c(EC_DS[[(t)]])>mean(EC_DS[[t]])){
    if(c(theta1[[t]]/sum(theta1[[t]]))%*%c(ER_in[[(t)]])<mean(ER_in[[t]])){
      # 3 constraint case
      M1 <- cbind(rbind(1,mean(EC_DS[[(t)]]),mean(ER_in[[(t)]])),
                  rbind(sum(theta3[[(t)]]),EC_DS[[(t)]]%*%theta3[[(t)]],ER_in[[(t)]]%*%theta3[[(t)]]),
                  rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
      M2 <- cbind(rbind(sum(theta1[[(t)]]),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                  rbind(1,mean(EC_DS[[(t)]]),mean(ER_in[[(t)]])),
                  rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
      M <- cbind(rbind(sum(theta1[[(t)]]),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                 rbind(sum(theta3[[(t)]]),EC_DS[[(t)]]%*%theta3[[(t)]],ER_in[[(t)]]%*%theta3[[(t)]]),
                 rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
      gamma1 <- det(M1)/det(M)
      gamma2 <- det(M2)/det(M)
      alpha1 <- c(gamma1*sum(theta3[[t]]))
      alpha2 <- c(gamma2*sum(theta2[[t]]))
      w1 <- (1-alpha1-alpha2)*theta1[[t]]/sum(theta1[[t]])
      w2 <- alpha1*theta3[[t]]/sum(theta3[[t]])
      w3 <- alpha2*theta2[[t]]/sum(theta2[[t]])
      w[[(t)]] = w1 +  w2 + w3
    }
    else{
      # 1 centrality constraint case
      alpha=(sum(theta3[[t]])*sum(theta1[[t]])*mean(EC_DS[[t]])-(sum(theta3[[t]]))^2)/(EC_DS[[(t)]]%*%theta3[[t]]*sum(theta1[[t]])-(sum(theta3[[t]]))^2)
      w[[(t)]] = alpha*theta3[[t]]/sum(theta3[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
    }}
  else{
    if(c(theta1[[t]]/sum(theta1[[t]]))%*%c(ER_in[[(t)]])<mean(ER_in[[t]])){
      # 1 expected return constraint case 
      alpha=c((sum(theta2[[t]])*sum(theta1[[t]])*mean(ER_in[[t]])-(sum(theta2[[t]]))^2)/(ER_in[[(t)]]%*%theta2[[t]]*sum(theta1[[t]])-(sum(theta2[[t]]))^2))
      w[[(t)]] = alpha*theta2[[t]]/sum(theta2[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
    }
    else{
      # global minimum variance case 
      w[[(t)]] = theta1[[t]]/sum(theta1[[t]])
    }
  }
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_3constraint_Dantzig<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_3constraint_Dantzig<-cumprod(return_network_3constraint_Dantzig)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_3constraint_Dantzig<-w

###### network portfolio fixed constraint with glasso ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute 3 constraint network portfolio ##
  glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
  if(c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS[[(t)]])>mean(EC_DS[[t]])){
    if(c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(ER_in[[(t)]])<mean(ER_in[[t]])){
      # 3 constraint case
      M1 <- cbind(rbind(1,mean(EC_DS[[(t)]]),mean(ER_in[[(t)]])),
                  rbind(sum(glasso.icov%*%EC_DS[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]]),
                  rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
      M2 <- cbind(rbind(sum(glasso.icov),EC_DS[[(t)]]%*%row_sums(glasso.icov),ER_in[[(t)]]%*%row_sums(glasso.icov)),
                  rbind(1,mean(EC_DS[[(t)]]),mean(ER_in[[(t)]])),
                  rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
      M <- cbind(rbind(sum(glasso.icov),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                 rbind(sum(glasso.icov%*%EC_DS[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]]),
                 rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
      gamma1 <- det(M1)/det(M)
      gamma2 <- det(M2)/det(M)
      alpha1 <- c(gamma1*sum(glasso.icov%*%EC_DS[[(t)]]))
      alpha2 <- c(gamma2*sum(glasso.icov%*%ER_in[[(t)]]))
      w1 <- (1-alpha1-alpha2)*row_sums(glasso.icov)/sum(glasso.icov)
      w2 <- alpha1*glasso.icov%*%EC_DS[[(t)]]/sum(glasso.icov%*%EC_DS[[(t)]])
      w3 <- alpha2*glasso.icov%*%ER_in[[(t)]]/sum(glasso.icov%*%ER_in[[(t)]])
      w[[(t)]] = c(w1 +  w2 + w3)
    }
    else{
      # 1 centrality constraint case
      alpha=c((sum(glasso.icov%*%EC_DS[[t]])*sum(glasso.icov)*mean(EC_DS[[t]])-(sum(glasso.icov%*%EC_DS[[t]]))^2)/(EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS[[t]]))^2))
      w[[(t)]] = c(alpha[1]*glasso.icov%*%EC_DS[[t]]/sum(glasso.icov%*%EC_DS[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
    }}
  else{
    if(c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(ER_in[[(t)]])<mean(ER_in[[t]])){
      # 1 expected return constraint case 
      alpha=c((sum(glasso.icov%*%ER_in[[t]])*sum(glasso.icov)*mean(ER_in[[t]])-(sum(glasso.icov%*%ER_in[[t]]))^2)/(ER_in[[(t)]]%*%glasso.icov%*%ER_in[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%ER_in[[t]]))^2))
      w[[(t)]] = c(alpha[1]*glasso.icov%*%ER_in[[t]]/sum(glasso.icov%*%ER_in[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
    }
    else{
      # global minimum variance case 
      w[[t]]=row_sums(glasso.icov)/sum(glasso.icov)
    }
  }
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_3constraint_glasso<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_3constraint_glasso<-cumprod(return_network_3constraint_glasso)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_3constraint_glasso<-w

###### network portfolio fixed constraint with plug in no short ######
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio without short ##
  net.gmin.port = network.3constraint.portfolio(EC_DS[[(t)]],ER_in[[t]], COV_in[[(t)]],mean(EC_DS[[(t)]]),mean(ER_in[[t]]),FALSE)
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_3constraint_noshort<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_3constraint_noshort<-cumprod(return_network_3constraint_noshort)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_3constraint_noshort<-w
###### network portfolio varying constraint with plug in ######
cumureturn_network_vary_with_phi_3constraint<-list()
return_network_vary_with_phi_3constraint<-list()
w_network_vary_with_phi_3constraint<-list()
quantl<-seq(0.1,0.9,0.1)
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## compute global minimum variance portfolio ##
    net.gmin.port = network.3constraint.portfolio(EC_DS[[(t)]],ER_in[[t]], COV_in[[(t)]],quantile(EC_DS[[(t)]],quantl[i]),mean(ER_in[[t]]),TRUE)
    w[[(t)]] =net.gmin.port$weights
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_3constraint[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_3constraint[[i]]<-cumprod(return_network_vary_with_phi_3constraint[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_3constraint[[i]]<-w
}

###### network portfolio varying constraint with Dantzig ######
cumureturn_network_vary_with_phi_3constraint_Dantzig<-list()
return_network_vary_with_phi_3constraint_Dantzig<-list()
w_network_vary_with_phi_3constraint_Dantzig<-list()
quantl<-seq(0.1,0.9,0.1)
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## network portfolio varying constraint with Dantzig estimation ##
    if((c(theta1[[t]]/sum(theta1[[t]]))%*%c(EC_DS[[(t)]]))>quantile(EC_DS[[(t)]],quantl[i])){
      if((c(theta1[[t]]/sum(theta1[[t]]))%*%c(ER_in[[(t)]]))<mean(ER_in[[t]])){
        # 3 constraint case
        M1 <- cbind(rbind(1,quantile(EC_DS[[(t)]],quantl[i]),mean(ER_in[[(t)]])),
                    rbind(sum(theta3[[(t)]]),EC_DS[[(t)]]%*%theta3[[(t)]],ER_in[[(t)]]%*%theta3[[(t)]]),
                    rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
        M2 <- cbind(rbind(sum(theta1[[(t)]]),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                    rbind(1,quantile(EC_DS[[(t)]],quantl[i]),mean(ER_in[[(t)]])),
                    rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
        M <- cbind(rbind(sum(theta1[[(t)]]),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                   rbind(sum(theta3[[(t)]]),EC_DS[[(t)]]%*%theta3[[(t)]],ER_in[[(t)]]%*%theta3[[(t)]]),
                   rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
        gamma1 <- det(M1)/det(M)
        gamma2 <- det(M2)/det(M)
        alpha1 <- c(gamma1*sum(theta3[[t]]))
        alpha2 <- c(gamma2*sum(theta2[[t]]))
        w1 <- (1-alpha1-alpha2)*theta1[[t]]/sum(theta1[[t]])
        w2 <- alpha1*theta3[[t]]/sum(theta3[[t]])
        w3 <- alpha2*theta2[[t]]/sum(theta2[[t]])
        w[[(t)]] = w1 +  w2 + w3
      }
      else{
        # 1 centrality constraint case
        alpha=c((sum(theta3[[t]])*sum(theta1[[t]])*quantile(EC_DS[[(t)]],quantl[i])-(sum(theta3[[t]]))^2)/(EC_DS[[(t)]]%*%theta3[[t]]*sum(theta1[[t]])-(sum(theta3[[t]]))^2))
        w[[(t)]] = alpha*theta3[[t]]/sum(theta3[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
      }}
    else{
      if((c(theta1[[t]]/sum(theta1[[t]]))%*%c(ER_in[[(t)]]))<mean(ER_in[[t]])){
        # 1 expected return constraint case 
        alpha=c((sum(theta2[[t]])*sum(theta1[[t]])*mean(ER_in[[t]])-(sum(theta2[[t]]))^2)/(ER_in[[(t)]]%*%theta2[[t]]*sum(theta1[[t]])-(sum(theta2[[t]]))^2))
        w[[(t)]] = alpha*theta2[[t]]/sum(theta2[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
      }
      else{
        # global minimum variance case 
        w[[(t)]] = theta1[[t]]/sum(theta1[[t]])
      }
    }
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_3constraint_Dantzig[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_3constraint_Dantzig[[i]]<-cumprod(return_network_vary_with_phi_3constraint_Dantzig[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_3constraint_Dantzig[[i]]<-w
}

###### network portfolio varying constraint with glasso ######
cumureturn_network_vary_with_phi_3constraint_glasso<-list()
return_network_vary_with_phi_3constraint_glasso<-list()
w_network_vary_with_phi_3constraint_glasso<-list()
quantl<-seq(0.1,0.9,0.1)
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## network portfolio varying constraint with glasso estimation ##
    glasso.icov=glasso(COV_in[[(t)]],rho=rho)$wi
    if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS[[(t)]]))>quantile(EC_DS[[(t)]],quantl[i])){
      if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(ER_in[[(t)]]))<mean(ER_in[[t]])){
        # 3 constraint case
        M1 <- cbind(rbind(1,quantile(EC_DS[[(t)]],quantl[i]),mean(ER_in[[(t)]])),
                    rbind(sum(glasso.icov%*%EC_DS[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]]),
                    rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
        M2 <- cbind(rbind(sum(glasso.icov),EC_DS[[(t)]]%*%row_sums(glasso.icov),ER_in[[(t)]]%*%row_sums(glasso.icov)),
                    rbind(1,quantile(EC_DS[[(t)]],quantl[i]),mean(ER_in[[(t)]])),
                    rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
        M <- cbind(rbind(sum(glasso.icov),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                   rbind(sum(glasso.icov%*%EC_DS[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]]),
                   rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
        gamma1 <- det(M1)/det(M)
        gamma2 <- det(M2)/det(M)
        alpha1 <- c(gamma1*sum(glasso.icov%*%EC_DS[[(t)]]))
        alpha2 <- c(gamma2*sum(glasso.icov%*%ER_in[[(t)]]))
        w1 <- (1-alpha1-alpha2)*row_sums(glasso.icov)/sum(glasso.icov)
        w2 <- alpha1*glasso.icov%*%EC_DS[[(t)]]/sum(glasso.icov%*%EC_DS[[(t)]])
        w3 <- alpha2*glasso.icov%*%ER_in[[(t)]]/sum(glasso.icov%*%ER_in[[(t)]])
        w[[(t)]] = c(w1 +  w2 + w3)
      }
      else{
        # 1 centrality constraint case
        alpha=c((sum(glasso.icov%*%EC_DS[[t]])*sum(glasso.icov)*quantile(EC_DS[[(t)]],quantl[i])-(sum(glasso.icov%*%EC_DS[[t]]))^2)/(EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS[[t]]))^2))
        w[[(t)]] = c(alpha[1]*glasso.icov%*%EC_DS[[t]]/sum(glasso.icov%*%EC_DS[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
      }}
    else{
      if(((row_sums(glasso.icov)/sum(glasso.icov))%*%ER_in[[(t)]])<mean(ER_in[[t]])){
        # 1 expected return constraint case 
        alpha=c((sum(glasso.icov%*%ER_in[[t]])*sum(glasso.icov)*mean(ER_in[[t]])-(sum(glasso.icov%*%ER_in[[t]]))^2)/(ER_in[[(t)]]%*%glasso.icov%*%ER_in[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%ER_in[[t]]))^2))
        w[[(t)]] = c(alpha[1]*glasso.icov%*%ER_in[[t]]/sum(glasso.icov%*%ER_in[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
      }
      else{
        # global minimum variance case 
        w[[t]]=row_sums(glasso.icov)/sum(glasso.icov)
      }
    }
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_3constraint_glasso[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_3constraint_glasso[[i]]<-cumprod(return_network_vary_with_phi_3constraint_glasso[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_3constraint_glasso[[i]]<-w
}

###### network portfolio varying constraint with plug in no short ######
cumureturn_network_vary_with_phi_3constraint_noshort<-list()
return_network_vary_with_phi_3constraint_noshort<-list()
w_network_vary_with_phi_3constraint_noshort<-list()
quantl<-seq(0.1,0.9,0.1)
for (i in 1:length(quantl)) {
  w<-list()
  cumureturn_temporal<-list()
  for(t in 1: length(W_in)){
    ## compute global minimum variance portfolio ##
    net.gmin.port = network.3constraint.portfolio(EC_DS[[(t)]],ER_in[[t]], COV_in[[(t)]],quantile(EC_DS[[(t)]],quantl[i]),mean(ER_in[[t]]),FALSE)
    w[[(t)]] =net.gmin.port$weights
    aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
    cumureturn_temporal[[t]]<-rowSums(aus)
  }
  return_network_vary_with_phi_3constraint_noshort[[i]]<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
  cumureturn_network_vary_with_phi_3constraint_noshort[[i]]<-cumprod(return_network_vary_with_phi_3constraint_noshort[[i]])
  w<-t(matrix(unlist(w),p,T.windows))
  colnames(w) = node.label
  w_network_vary_with_phi_3constraint_noshort[[i]]<-w
}

###### network portfolio varying constraint with glasso no short ######

###### network portfolio data-driven constraint with plug in ######
cumureturn_network_datadriven_phistar_3constraint<-list()
return_network_datadriven_phistar_3constraint<-list()
phi_star<-centrality_equal_portfolio
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  net.gmin.port = network.3constraint.portfolio(EC_DS[[(t)]],ER_in[[t]], COV_in[[(t)]],phi_star[[t]],mean(ER_in[[t]]),TRUE)
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_datadriven_phistar_3constraint<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_3constraint<-cumprod(return_network_datadriven_phistar_3constraint)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_3constraint<-w

###### network portfolio data-driven constraint with Dantzig ######
cumureturn_network_datadriven_phistar_3constraint_Dantzig<-list()
return_network_datadriven_phistar_3constraint_Dantzig<-list()
phi_star<-centrality_equal_portfolio
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  if((c(theta1[[t]]/sum(theta1[[t]]))%*%c(EC_DS[[(t)]]))>phi_star[[t]]){
    if((c(theta1[[t]]/sum(theta1[[t]]))%*%c(ER_in[[(t)]]))<phi_star[[t]]){
      # 3 constraint case
      M1 <- cbind(rbind(1,phi_star[[t]],mean(ER_in[[(t)]])),
                  rbind(sum(theta3[[(t)]]),EC_DS[[(t)]]%*%theta3[[(t)]],ER_in[[(t)]]%*%theta3[[(t)]]),
                  rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
      M2 <- cbind(rbind(sum(theta1[[(t)]]),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                  rbind(1,phi_star[[t]],mean(ER_in[[(t)]])),
                  rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
      M <- cbind(rbind(sum(theta1[[(t)]]),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                 rbind(sum(theta3[[(t)]]),EC_DS[[(t)]]%*%theta3[[(t)]],ER_in[[(t)]]%*%theta3[[(t)]]),
                 rbind(sum(theta2[[(t)]]),EC_DS[[(t)]]%*%theta2[[(t)]],ER_in[[(t)]]%*%theta2[[(t)]]))
      gamma1 <- det(M1)/det(M)
      gamma2 <- det(M2)/det(M)
      alpha1 <- c(gamma1*sum(theta3[[t]]))
      alpha2 <- c(gamma2*sum(theta2[[t]]))
      w1 <- (1-alpha1-alpha2)*theta1[[t]]/sum(theta1[[t]])
      w2 <- alpha1*theta3[[t]]/sum(theta3[[t]])
      w3 <- alpha2*theta2[[t]]/sum(theta2[[t]])
      w[[(t)]] = w1 +  w2 + w3
    }
    else{
      # 1 centrality constraint case
      alpha=c((sum(theta3[[t]])*sum(theta1[[t]])*phi_star[[t]]-(sum(theta3[[t]]))^2)/(EC_DS[[(t)]]%*%theta3[[t]]*sum(theta1[[t]])-(sum(theta3[[t]]))^2))
      w[[(t)]] = alpha*theta3[[t]]/sum(theta3[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
    }}
  else{
    if(((theta1[[t]]/sum(theta1[[t]]))%*%ER_in[[(t)]])<mean(ER_in[[t]])){
      # 1 expected return constraint case 
      alpha=c((sum(theta2[[t]])*sum(theta1[[t]])*mean(ER_in[[t]])-(sum(theta2[[t]]))^2)/(ER_in[[(t)]]%*%theta2[[t]]*sum(theta1[[t]])-(sum(theta2[[t]]))^2))
      w[[(t)]] = alpha*theta2[[t]]/sum(theta2[[t]])+(1-alpha)*theta1[[t]]/sum(theta1[[t]])
    }
    else{
      # global minimum variance case 
      w[[(t)]] = theta1[[t]]/sum(theta1[[t]])
    }
  }
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_datadriven_phistar_3constraint_Dantzig<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_3constraint_Dantzig<-cumprod(return_network_datadriven_phistar_3constraint_Dantzig)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_3constraint_Dantzig<-w

###### network portfolio data-driven constraint with glasso ######
cumureturn_network_datadriven_phistar_3constraint_glasso<-list()
return_network_datadriven_phistar_3constraint_glasso<-list()
phi_star<-centrality_equal_portfolio
w<-list()
cumureturn_temporal<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(EC_DS[[(t)]]))>phi_star[[t]]){
    if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%c(ER_in[[(t)]]))<phi_star[[t]]){
      # 3 constraint case
      M1 <- cbind(rbind(1,phi_star[[t]],mean(ER_in[[(t)]])),
                  rbind(sum(glasso.icov%*%EC_DS[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]]),
                  rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
      M2 <- cbind(rbind(sum(glasso.icov),EC_DS[[(t)]]%*%row_sums(glasso.icov),ER_in[[(t)]]%*%row_sums(glasso.icov)),
                  rbind(1,phi_star[[t]],mean(ER_in[[(t)]])),
                  rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
      M <- cbind(rbind(sum(glasso.icov),EC_DS[[(t)]]%*%theta1[[(t)]],ER_in[[(t)]]%*%theta1[[(t)]]),
                 rbind(sum(glasso.icov%*%EC_DS[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%EC_DS[[(t)]]),
                 rbind(sum(glasso.icov%*%ER_in[[(t)]]),EC_DS[[(t)]]%*%glasso.icov%*%ER_in[[(t)]],ER_in[[(t)]]%*%glasso.icov%*%ER_in[[(t)]]))
      gamma1 <- det(M1)/det(M)
      gamma2 <- det(M2)/det(M)
      alpha1 <- c(gamma1*sum(glasso.icov%*%EC_DS[[(t)]]))
      alpha2 <- c(gamma2*sum(glasso.icov%*%ER_in[[(t)]]))
      w1 <- (1-alpha1-alpha2)*row_sums(glasso.icov)/sum(glasso.icov)
      w2 <- alpha1*glasso.icov%*%EC_DS[[(t)]]/sum(glasso.icov%*%EC_DS[[(t)]])
      w3 <- alpha2*glasso.icov%*%ER_in[[(t)]]/sum(glasso.icov%*%ER_in[[(t)]])
      w[[(t)]] = c(w1 +  w2 + w3)
    }
    else{
      # 1 centrality constraint case
      alpha=(sum(glasso.icov%*%EC_DS[[t]])*sum(glasso.icov)*phi_star[[t]]-(sum(glasso.icov%*%EC_DS[[t]]))^2)/(EC_DS[[(t)]]%*%glasso.icov%*%EC_DS[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%EC_DS[[t]]))^2)
      w[[(t)]] = c(alpha[1]*glasso.icov%*%EC_DS[[t]]/sum(glasso.icov%*%EC_DS[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
    }}
  else{
    if((c(row_sums(glasso.icov)/sum(glasso.icov))%*%ER_in[[(t)]])<mean(ER_in[[t]])){
      # 1 expected return constraint case 
      alpha=(sum(glasso.icov%*%ER_in[[t]])*sum(glasso.icov)*mean(ER_in[[t]])-(sum(glasso.icov%*%ER_in[[t]]))^2)/(ER_in[[(t)]]%*%glasso.icov%*%ER_in[[t]]*sum(glasso.icov)-(sum(glasso.icov%*%ER_in[[t]]))^2)
      w[[(t)]] = c(alpha[1]*glasso.icov%*%ER_in[[t]]/sum(glasso.icov%*%ER_in[[t]])+(1-alpha[1])*row_sums(glasso.icov)/sum(glasso.icov))
    }
    else{
      # global minimum variance case 
      w[[t]]=row_sums(glasso.icov)/sum(glasso.icov)
    }
  }
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_datadriven_phistar_3constraint_glasso<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_3constraint_glasso<-cumprod(return_network_datadriven_phistar_3constraint_glasso)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_3constraint_glasso<-w

###### network portfolio data-driven constraint with plug in no short ######
cumureturn_network_datadriven_phistar_3constraint_noshort<-list()
return_network_datadriven_phistar_3constraint_noshort<-list()
cumureturn_temporal<-list()
phi_star<-centrality_equal_portfolio
w<-list()
for(t in 1: length(W_in)){
  ## compute global minimum variance portfolio ##
  net.gmin.port = network.3constraint.portfolio(EC_DS[[(t)]],ER_in[[t]], COV_in[[(t)]],phi_star[[t]],mean(ER_in[[t]]),FALSE)
  w[[(t)]] =net.gmin.port$weights
  aus<-as.matrix(repmat(w[[(t)]],22,1)*W_out[[t]])
  cumureturn_temporal[[t]]<-rowSums(aus)
}
return_network_datadriven_phistar_3constraint_noshort<-as.matrix(cbind(unlist(cumureturn_temporal))+1)
cumureturn_network_datadriven_phistar_3constraint_noshort<-cumprod(return_network_datadriven_phistar_3constraint_noshort)
w<-t(matrix(unlist(w),p,T.windows))
colnames(w) = node.label
w_network_datadriven_phistar_3constraint_noshort<-w

###### network portfolio data-driven constraint with Dantzig no short ######

#### Save portfolios ####
Portfolio.Scenario <- list("return_minVar"=return_minVar,"cumureturn_minVar"=cumureturn_minVar,"w_minVar"=w_minVar,
                           "return_minVar_Dantzig"=return_minVar_Dantzig,"cumureturn_minVar_Dantzig"=cumureturn_minVar_Dantzig,"w_minVar_Dantzig"=w_minVar_Dantzig,
                           "return_minVar_glasso"=return_minVar_glasso,"cumureturn_minVar_glasso"=cumureturn_minVar_glasso,"w_minVar_glasso"=w_minVar_glasso,
                           "return_minVar_noshort"=return_minVar_noshort,"cumureturn_minVar_noshort"=cumureturn_minVar_noshort,"w_minVar_noshort"=w_minVar_noshort,
                           "return_meanVar"=return_meanVar,"cumureturn_meanVar"=cumureturn_meanVar,"w_meanVar"=w_meanVar,
                           "return_meanVar_Dantzig"=return_meanVar_Dantzig,"cumureturn_meanVar_Dantzig"=cumureturn_meanVar_Dantzig,"w_meanVar_Dantzig"=w_meanVar_Dantzig,
                           "return_meanVar_glasso"=return_meanVar_glasso,"cumureturn_meanVar_glasso"=cumureturn_meanVar_glasso,"w_meanVar_glasso"=w_meanVar_glasso,
                           "return_meanVar_noshort"=return_meanVar_noshort,"cumureturn_meanVar_noshort"=cumureturn_meanVar_noshort,"w_meanVar_noshort"=w_meanVar_noshort,
                           "return_equal"=return_equal,"cumureturn_equal"=cumureturn_equal,"w_equal"=w_equal,
                           "return_network_2constraint"=return_network_2constraint,"cumureturn_network_2constraint"=cumureturn_network_2constraint,"w_network_2constraint"=w_network_2constraint,
                           "return_network_2constraint_Dantzig"=return_network_2constraint_Dantzig,"cumureturn_network_2constraint_Dantzig"=cumureturn_network_2constraint_Dantzig,"w_network_2constraint_Dantzig"=w_network_2constraint_Dantzig,
                           "return_network_2constraint_glasso"=return_network_2constraint_glasso,"cumureturn_network_2constraint_glasso"=cumureturn_network_2constraint_glasso,"w_network_2constraint_glasso"=w_network_2constraint_glasso,
                           "return_network_2constraint_noshort"=return_network_2constraint_noshort,"cumureturn_network_2constraint_noshort"=cumureturn_network_2constraint_noshort,"w_network_2constraint_noshort"=w_network_2constraint_noshort,
                           "return_network_vary_with_phi"=return_network_vary_with_phi,"cumureturn_network_vary_with_phi"=cumureturn_network_vary_with_phi,"w_network_vary_with_phi"=w_network_vary_with_phi,
                           "return_network_vary_with_phi_Dantzig"=return_network_vary_with_phi_Dantzig,"cumureturn_network_vary_with_phi_Dantzig"=cumureturn_network_vary_with_phi_Dantzig,"w_network_vary_with_phi_Dantzig"=w_network_vary_with_phi_Dantzig,
                           "return_network_vary_with_phi_glasso"=return_network_vary_with_phi_glasso,"cumureturn_network_vary_with_phi_glasso"=cumureturn_network_vary_with_phi_glasso,"w_network_vary_with_phi_glasso"=w_network_vary_with_phi_glasso,
                           "return_network_vary_with_phi_noshort"=return_network_vary_with_phi_noshort,"cumureturn_network_vary_with_phi_noshort"=cumureturn_network_vary_with_phi_noshort,"w_network_vary_with_phi_noshort"=w_network_vary_with_phi_noshort,
                           "return_network_datadriven_phistar"=return_network_datadriven_phistar,"cumureturn_network_datadriven_phistar"=cumureturn_network_datadriven_phistar,"w_network_datadriven_phistar"=w_network_datadriven_phistar,
                           "return_network_datadriven_phistar_Dantzig"=return_network_datadriven_phistar_Dantzig,"cumureturn_network_datadriven_phistar_Dantzig"=cumureturn_network_datadriven_phistar_Dantzig,"w_network_datadriven_phistar_Dantzig"=w_network_datadriven_phistar_Dantzig,
                           "return_network_datadriven_phistar_glasso"=return_network_datadriven_phistar_glasso,"cumureturn_network_datadriven_phistar_glasso"=cumureturn_network_datadriven_phistar_glasso,"w_network_datadriven_phistar_glasso"=w_network_datadriven_phistar_glasso,
                           "return_network_datadriven_phistar_noshort"=return_network_datadriven_phistar_noshort,"cumureturn_network_datadriven_phistar_noshort"=cumureturn_network_datadriven_phistar_noshort,"w_network_datadriven_phistar_noshort"=w_network_datadriven_phistar_noshort,
                           "return_network_3constraint"=return_network_3constraint,"cumureturn_network_3constraint"=cumureturn_network_3constraint,"w_network_3constraint"=w_network_3constraint,
                           "return_network_3constraint_Dantzig"=return_network_3constraint_Dantzig,"cumureturn_network_3constraint_Dantzig"=cumureturn_network_3constraint_Dantzig,"w_network_3constraint_Dantzig"=w_network_3constraint_Dantzig,
                           "return_network_3constraint_glasso"=return_network_3constraint_glasso,"cumureturn_network_3constraint_glasso"=cumureturn_network_3constraint_glasso,"w_network_3constraint_glasso"=w_network_3constraint_glasso,
                           "return_network_3constraint_noshort"=return_network_3constraint_noshort,"cumureturn_network_3constraint_noshort"=cumureturn_network_3constraint_noshort,"w_network_3constraint_noshort"=w_network_3constraint_noshort,
                           "return_network_vary_with_phi_3constraint"=return_network_vary_with_phi_3constraint,"cumureturn_network_vary_with_phi_3constraint"=cumureturn_network_vary_with_phi_3constraint,"w_network_vary_with_phi_3constraint"=w_network_vary_with_phi_3constraint,
                           "return_network_vary_with_phi_3constraint_Dantzig"=return_network_vary_with_phi_3constraint_Dantzig,"cumureturn_network_vary_with_phi_3constraint_Dantzig"=cumureturn_network_vary_with_phi_3constraint_Dantzig,"w_network_vary_with_phi_3constraint_Dantzig"=w_network_vary_with_phi_3constraint_Dantzig,
                           "return_network_vary_with_phi_3constraint_glasso"=return_network_vary_with_phi_3constraint_glasso,"cumureturn_network_vary_with_phi_3constraint_glasso"=cumureturn_network_vary_with_phi_3constraint_glasso,"w_network_vary_with_phi_3constraint_glasso"=w_network_vary_with_phi_3constraint_glasso,
                           "return_network_vary_with_phi_3constraint_noshort"=return_network_vary_with_phi_3constraint_noshort,"cumureturn_network_vary_with_phi_3constraint_noshort"=cumureturn_network_vary_with_phi_3constraint_noshort,"w_network_vary_with_phi_3constraint_noshort"=w_network_vary_with_phi_3constraint_noshort,
                           "return_network_datadriven_phistar_3constraint"=return_network_datadriven_phistar_3constraint,"cumureturn_network_datadriven_phistar_3constraint"=cumureturn_network_datadriven_phistar_3constraint,"w_network_datadriven_phistar_3constraint"=w_network_datadriven_phistar_3constraint,
                           "return_network_datadriven_phistar_3constraint_Dantzig"=return_network_datadriven_phistar_3constraint_Dantzig,"cumureturn_network_datadriven_phistar_3constraint_Dantzig"=cumureturn_network_datadriven_phistar_3constraint_Dantzig,"w_network_datadriven_phistar_3constraint_Dantzig"=w_network_datadriven_phistar_3constraint_Dantzig,
                           "return_network_datadriven_phistar_3constraint_glasso"=return_network_datadriven_phistar_3constraint_glasso,"cumureturn_network_datadriven_phistar_3constraint_glasso"=cumureturn_network_datadriven_phistar_3constraint_glasso,"w_network_datadriven_phistar_3constraint_glasso"=w_network_datadriven_phistar_3constraint_glasso,
                           "return_network_datadriven_phistar_3constraint_noshort"=return_network_datadriven_phistar_3constraint_noshort,"cumureturn_network_datadriven_phistar_3constraint_noshort"=cumureturn_network_datadriven_phistar_3constraint_noshort,"w_network_datadriven_phistar_3constraint_noshort"=w_network_datadriven_phistar_3constraint_noshort)
save(Portfolio.Scenario,file="Portfolios_WS500_20240723.RData")
