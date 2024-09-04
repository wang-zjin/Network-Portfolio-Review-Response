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
for(t in 0: (floor((1857-500)/22)-1)){
  W[[(t+1)]]=returnstd[(1+t*22):(522+t*22),]
}
W_in<-list()
W_out<-list()
for(t in 0: (floor((1857-500)/22)-1)){
  W_in[[(t+1)]]=W[[t+1]][c(1:500),]
  W_out[[(t+1)]]=W[[t+1]][c(501:522),]
}
T.windows<-length(W)
# correlation matrix, Expected return, covariance matrix
C_in <- list()
ER_in <- list()
COV_in <- list()
EC_in <- list()
for(t in 1: length(W_in)){
  C_in[[(t)]] =cor(W_in[[(t)]])
  ER_in[[(t)]] = colMeans(W_in[[(t)]])
  COV_in[[(t)]] = cov(W_in[[(t)]])
  network_port = network.correlation(W_in[[(t)]])
  EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
  max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
}

##### Dantzig selector estimation for eigenvector centrality #####


###### CV for tuning lambda in estimation eigenvector centrality using the first 500 data 
###### Task 1: tune lambda in estimation eigenvector centrality using rolling window of size 500 ######

tic("Tuning parameter of Dantzig estimation")

lmd.EC.Dantzig.list = list()
for (t in c(0:60)) {
  n<-dim(W_in[[(t+1)]])[1]
  B=100
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  lmd.i=c()
  for (valid.block in 1:n.block) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=W_in[[(t+1)]][train.ind,]
    returnstd.valid=W_in[[(t+1)]][valid.ind,]
    mu.train=rep(0,p)
    mu.valid=rep(0,p)
    
    corr.train=cor(returnstd.train)
    corr.train[is.na(corr.train)]=0
    A.train=corr.train-diag(1,p,p)
    corr.valid=cor(returnstd.valid)
    corr.valid[is.na(corr.valid)]=0
    A.valid=corr.valid-diag(1,p,p)
    
    cov.train=A.train-diag(max(eigen(A.train)$value),p,p)
    cov.valid=A.valid-diag(max(eigen(A.valid)$value),p,p)
    lambda.grid=seq(0,1,length=101)[2:101]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l.lmd=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      # lmd=0.2
      lin.train=linfun3(cov.train,mu.train,lmd,abs(eigen(cov.valid)$vector[1,1]))
      # max(lin.train)
      # max(cov.train%*%lin.train)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l.lmd=c(cv.l.lmd, lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l.lmd[which(cv.l.error==min(cv.l.error))])
  }
  lmd=mean(lmd.i[lmd.i<Inf])
  # lmd.EC.Dantzig <- 0.682
  lmd.EC.Dantzig <- lmd
  lmd.EC.Dantzig.list[[t+1]] = lmd.EC.Dantzig
}
toc()

save(lmd.EC.Dantzig.list,file="Dantzig_lambda_rolling_window.RData")
load("Dantzig_lambda_rolling_window.RData")

mean(unlist(lmd.EC.Dantzig.list))
std(unlist(lmd.EC.Dantzig.list))

#   EC_DS <- eigenvector centrality estimated by Dantzig selector
lmd.EC.Dantzig = 0.682
EC_DS<-list()
a=c()
for (t in 1: length(W_in)) {
  print(t)
  EC_DS[[t]] =linfun3(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p))$vector[1,1]))
  EC_DS[[t]]=EC_DS[[t]]/max(EC_DS[[t]])
  # boxplot(EC_DS[[t]])
  a[t]=sum(EC_DS[[t]]==0)
}
a
eigenvector_centrality = list("eigenvector_absolute_value"=EC_in,"eigenvector_centrality_Dantzig"=EC_DS,
                              "zeors_in_EC_DS"=a,"covariance_matrix"=COV_in,"correlation_matrix"=C_in,
                              "expected_return"=ER_in)
save(eigenvector_centrality,file = "eigenvector_centrality_20230513.RData")

load("eigenvector_centrality_20230513.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

############ Different window size: 750, 1000 
###### Task 2: window size 750 ######
tic("Tuning parameter of Dantzig estimation")
# rolling window
W<-list()
window_size = 750
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
# correlation matrix, Expected return, covariance matrix
C_in <- list()
ER_in <- list()
COV_in <- list()
EC_in <- list()
for(t in 1: length(W_in)){
  C_in[[(t)]] =cor(W_in[[(t)]])
  ER_in[[(t)]] = colMeans(W_in[[(t)]])
  COV_in[[(t)]] = cov(W_in[[(t)]])
  network_port = network.correlation(W_in[[(t)]])
  EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
  max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
}

lmd.EC.Dantzig.list = list()
for (t in c(0:49)) {
  n<-dim(W_in[[(t+1)]])[1]
  B=100
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  lmd.i=c()
  for (valid.block in 1:n.block) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=W_in[[(t+1)]][train.ind,]
    returnstd.valid=W_in[[(t+1)]][valid.ind,]
    mu.train=rep(0,p)
    mu.valid=rep(0,p)
    
    corr.train=cor(returnstd.train)
    corr.train[is.na(corr.train)]=0
    A.train=corr.train-diag(1,p,p)
    corr.valid=cor(returnstd.valid)
    corr.valid[is.na(corr.valid)]=0
    A.valid=corr.valid-diag(1,p,p)
    
    cov.train=A.train-diag(max(eigen(A.train)$value),p,p)
    cov.valid=A.valid-diag(max(eigen(A.valid)$value),p,p)
    lambda.grid=seq(0,1,length=101)[2:101]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l.lmd=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      # lmd=0.2
      lin.train=linfun3(cov.train,mu.train,lmd,abs(eigen(cov.valid)$vector[1,1]))
      # max(lin.train)
      # max(cov.train%*%lin.train)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l.lmd=c(cv.l.lmd, lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l.lmd[which(cv.l.error==min(cv.l.error))])
  }
  lmd=mean(lmd.i[lmd.i<Inf])
  # lmd.EC.Dantzig <- 0.682
  lmd.EC.Dantzig <- lmd
  lmd.EC.Dantzig.list[[t+1]] = lmd.EC.Dantzig
}
toc()
save(lmd.EC.Dantzig.list,file="Dantzig_lambda_rolling_window_750.RData")
# Tuning parameter of Dantzig estimation: 133072.148 sec elapsed
mean(unlist(lmd.EC.Dantzig.list)) # 0.5449429
std(unlist(lmd.EC.Dantzig.list)) # 0.1774072
load("Dantzig_lambda_rolling_window_750.RData")
lmd.EC.Dantzig.WS750=lmd.EC.Dantzig.list[[1]] # 0.5085714
lmd.EC.Dantzig = lmd.EC.Dantzig.WS750

#   EC_DS <- eigenvector centrality estimated by Dantzig selector
EC_DS<-list()
a=c()
for (t in 1: length(W_in)) {
  print(t)
  EC_DS[[t]] =linfun3(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p))$vector[1,1]))
  EC_DS[[t]]=EC_DS[[t]]/max(EC_DS[[t]])
  # boxplot(EC_DS[[t]])
  a[t]=sum(EC_DS[[t]]==0)
}
a
eigenvector_centrality = list("eigenvector_absolute_value"=EC_in,"eigenvector_centrality_Dantzig"=EC_DS,
                              "zeors_in_EC_DS"=a,"covariance_matrix"=COV_in,"correlation_matrix"=C_in,
                              "expected_return"=ER_in)
save(eigenvector_centrality,file = "eigenvector_centrality_WS750_20240723.RData")

load("eigenvector_centrality_WS750_20240723.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

###### Task 3: window size 1000 ######
tic("Tuning parameter of Dantzig estimation")
# rolling window
W<-list()
window_size = 1000
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
# correlation matrix, Expected return, covariance matrix
C_in <- list()
ER_in <- list()
COV_in <- list()
EC_in <- list()
for(t in 1: length(W_in)){
  C_in[[(t)]] =cor(W_in[[(t)]])
  ER_in[[(t)]] = colMeans(W_in[[(t)]])
  COV_in[[(t)]] = cov(W_in[[(t)]])
  network_port = network.correlation(W_in[[(t)]])
  EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
  max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
}

lmd.EC.Dantzig.list = list()
for (t in c(0:37)) {
  n<-dim(W_in[[(t+1)]])[1]
  B=100
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  lmd.i=c()
  for (valid.block in 1:n.block) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=W_in[[(t+1)]][train.ind,]
    returnstd.valid=W_in[[(t+1)]][valid.ind,]
    mu.train=rep(0,p)
    mu.valid=rep(0,p)
    
    corr.train=cor(returnstd.train)
    corr.train[is.na(corr.train)]=0
    A.train=corr.train-diag(1,p,p)
    corr.valid=cor(returnstd.valid)
    corr.valid[is.na(corr.valid)]=0
    A.valid=corr.valid-diag(1,p,p)
    
    cov.train=A.train-diag(max(eigen(A.train)$value),p,p)
    cov.valid=A.valid-diag(max(eigen(A.valid)$value),p,p)
    lambda.grid=seq(0,1,length=101)[2:101]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l.lmd=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      # lmd=0.2
      lin.train=linfun3(cov.train,mu.train,lmd,abs(eigen(cov.valid)$vector[1,1]))
      # max(lin.train)
      # max(cov.train%*%lin.train)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l.lmd=c(cv.l.lmd, lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l.lmd[which(cv.l.error==min(cv.l.error))])
  }
  lmd=mean(lmd.i[lmd.i<Inf])
  # lmd.EC.Dantzig <- 0.682
  lmd.EC.Dantzig <- lmd
  lmd.EC.Dantzig.list[[t+1]] = lmd.EC.Dantzig
}
toc()
save(lmd.EC.Dantzig.list,file="Dantzig_lambda_rolling_window_1000.RData")
# Tuning parameter of Dantzig estimation: 133072.148 sec elapsed
load("Dantzig_lambda_rolling_window_windowsize1000.RData")
mean(unlist(lmd.EC.Dantzig.list)) # 0.5449429
std(unlist(lmd.EC.Dantzig.list)) # 0.1774072
lmd.EC.Dantzig.WS1000=lmd.EC.Dantzig.list[[1]]
lmd.EC.Dantzig = lmd.EC.Dantzig.WS1000

#   EC_DS <- eigenvector centrality estimated by Dantzig selector
EC_DS<-list()
a=c()
for (t in 1: length(W_in)) {
  print(t)
  EC_DS[[t]] =linfun3(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p))$vector[1,1]))
  EC_DS[[t]]=EC_DS[[t]]/max(EC_DS[[t]])
  # boxplot(EC_DS[[t]])
  a[t]=sum(EC_DS[[t]]==0)
}
a
eigenvector_centrality = list("eigenvector_absolute_value"=EC_in,"eigenvector_centrality_Dantzig"=EC_DS,
                              "zeors_in_EC_DS"=a,"covariance_matrix"=COV_in,"correlation_matrix"=C_in,
                              "expected_return"=ER_in)
save(eigenvector_centrality,file = "eigenvector_centrality_WS1000_20240723.RData")

load("eigenvector_centrality_WS1000_20240723.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

#### Dantzig selector estimation for portfolio parameters ####

##### CV for tuning lambdas in portfolio solutions

##### CV for tuning lambda in estimation Sigma^-1 1 using the first 500/750/1000 data ######
###### Task 1: window size 500  ######

# Use for recurrence
window_size = 500
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=rep(1,p)
  mu.valid=rep(1,p)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    lin.train=linfun1(cov.train,mu.train,lmd)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd1.WS500=mean(lmd.i)
lmd1=0.4636

###### Task 2: window size 750  ######

# Use for recurrence
window_size = 750
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=rep(1,p)
  mu.valid=rep(1,p)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    lin.train=linfun1(cov.train,mu.train,lmd)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd1.WS750=mean(lmd.i) # = 0.3592
lmd1=0.3592

###### Task 3: window size 1000  ######

# Use for recurrence
window_size = 1000
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=rep(1,p)
  mu.valid=rep(1,p)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    lin.train=linfun1(cov.train,mu.train,lmd)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd1.WS1000=mean(lmd.i) # =0.2746
lmd1=0.2746

##### CV for tuning lambda in estimation Sigma^-1 mu, using the first 500/750/1000 data ######
###### Task 1: window size 500  ######
window_size = 500
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=colMeans(returnstd.train)
  mu.valid=colMeans(returnstd.valid)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  lambda.grid=seq(min(max(abs(mu.train))/100,min(abs(mu.train))), 0.01,length=101)[2:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    lin.train=linfun1(cov.train,mu.train,lmd)
    sum(lin.train==0)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd2.WS500=mean(lmd.i)
lmd2=0.002142681

###### Task 2: window size 750  ######
window_size = 750
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=colMeans(returnstd.train)
  mu.valid=colMeans(returnstd.valid)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  lambda.grid=seq(min(max(abs(mu.train))/100,min(abs(mu.train))), 0.01,length=101)[2:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    lin.train=linfun1(cov.train,mu.train,lmd)
    sum(lin.train==0)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd2.WS750=mean(lmd.i) # = 0.001682136
lmd2=0.001682136

###### Task 3: window size 1000  ######
window_size = 1000
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=colMeans(returnstd.train)
  mu.valid=colMeans(returnstd.valid)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  lambda.grid=seq(min(max(abs(mu.train))/100,min(abs(mu.train))), 0.01,length=101)[2:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    lin.train=linfun1(cov.train,mu.train,lmd)
    sum(lin.train==0)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd2.WS1000=mean(lmd.i) # =0.001122977
lmd2=0.001122977

##### CV for tuning lambda in estimation Sigma^-1 phi, using the first 500/750/1000 data ######
###### Task 1: window size 500  ######
window_size = 500
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  mu.train=EC_DS[[1]]
  mu.valid=EC_DS[[1]]
  # mu.train=linfun3(cor(returnstd.train)-diag(1,p,p)-diag(max(eigen(cor(returnstd.train))$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(cor(returnstd.train)-diag(1,p,p)-diag(max(eigen(cor(returnstd.train))$value),p,p))$vector[1,1]))
  # mu.train=mu.train/max(mu.train)
  # sum(mu.train==0)
  # mu.valid=linfun3(cor(returnstd.valid)-diag(1,p,p)-diag(max(eigen(cor(returnstd.valid))$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(cor(returnstd.valid)-diag(1,p,p)-diag(max(eigen(cor(returnstd.valid))$value),p,p))$vector[1,1]))
  # mu.valid=mu.valid/max(mu.valid)
  # sum(mu.valid==0)
  lambda.grid=seq(0.1, max(mu.train),length=101)[1:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    # lmd=0.1
    lin.train=linfun1(cov.train,mu.train,lmd)
    # sum(lin.train==0)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd3.WS500=mean(lmd.i)
lmd3=0.4132

###### Task 2: window size 750  ######
window_size = 750
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  mu.train=EC_DS[[1]]
  mu.valid=EC_DS[[1]]
  # mu.train=linfun3(cor(returnstd.train)-diag(1,p,p)-diag(max(eigen(cor(returnstd.train))$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(cor(returnstd.train)-diag(1,p,p)-diag(max(eigen(cor(returnstd.train))$value),p,p))$vector[1,1]))
  # mu.train=mu.train/max(mu.train)
  # sum(mu.train==0)
  # mu.valid=linfun3(cor(returnstd.valid)-diag(1,p,p)-diag(max(eigen(cor(returnstd.valid))$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(cor(returnstd.valid)-diag(1,p,p)-diag(max(eigen(cor(returnstd.valid))$value),p,p))$vector[1,1]))
  # mu.valid=mu.valid/max(mu.valid)
  # sum(mu.valid==0)
  lambda.grid=seq(0.1, max(mu.train),length=101)[1:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    # lmd=0.1
    lin.train=linfun1(cov.train,mu.train,lmd)
    # sum(lin.train==0)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd3.WS750=mean(lmd.i) # =0.2782
lmd3=0.2782

###### Task 3: window size 1000  ######
window_size = 1000
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
lmd.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  mu.train=EC_DS[[1]]
  mu.valid=EC_DS[[1]]
  # mu.train=linfun3(cor(returnstd.train)-diag(1,p,p)-diag(max(eigen(cor(returnstd.train))$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(cor(returnstd.train)-diag(1,p,p)-diag(max(eigen(cor(returnstd.train))$value),p,p))$vector[1,1]))
  # mu.train=mu.train/max(mu.train)
  # sum(mu.train==0)
  # mu.valid=linfun3(cor(returnstd.valid)-diag(1,p,p)-diag(max(eigen(cor(returnstd.valid))$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(cor(returnstd.valid)-diag(1,p,p)-diag(max(eigen(cor(returnstd.valid))$value),p,p))$vector[1,1]))
  # mu.valid=mu.valid/max(mu.valid)
  # sum(mu.valid==0)
  lambda.grid=seq(0.1, max(mu.train),length=101)[1:100]
  l.lambda=length(lambda.grid)
  cv.l.error=NULL
  cv.l=NULL
  for(i in 1:l.lambda){
    lmd=lambda.grid[i]
    print(i)
    # lmd=0.1
    lin.train=linfun1(cov.train,mu.train,lmd)
    # sum(lin.train==0)
    if(!(all(lin.train==0))){
      error=sum((cov.valid%*%lin.train-mu.valid)^2)
      cv.l.error=c(cv.l.error,error)
      cv.l=c(cv.l,lmd)
    }
  }
  lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
}
lmd.i
lmd3.WS1000=mean(lmd.i) # =0.2656
lmd3=0.2656

##### CV for tuning parameter in glasso, using the first 500/750/1000 data #####
###### Task 1: window size 500  ######
window_size = 500
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
rho.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=colMeans(returnstd.train)
  mu.valid=colMeans(returnstd.valid)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  rho.grid=seq(0,0.8,length=101)[2:101]
  l.rho=length(rho.grid)
  cv.rho.error=NULL
  cv.rho=NULL
  for (i in 1:l.rho){
    rho=rho.grid[i]
    prec.glasso=glasso(cov.train,rho=rho)$wi
    error=sum((prec.glasso%*%cov.valid-diag(rep(1,p)))^2)
    cv.rho.error=c(cv.rho.error,error)
    cv.rho=c(cv.rho,rho)
    print(i)
  }
  rho.i[valid.block]=cv.rho[which(cv.rho.error==min(cv.rho.error))]
}
rho.WS500=mean(rho.i)
rho<-0.016

###### Task 2: window size 750  ######
window_size = 750
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
rho.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=colMeans(returnstd.train)
  mu.valid=colMeans(returnstd.valid)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  rho.grid=seq(0,0.8,length=101)[2:101]
  l.rho=length(rho.grid)
  cv.rho.error=NULL
  cv.rho=NULL
  for (i in 1:l.rho){
    rho=rho.grid[i]
    prec.glasso=glasso(cov.train,rho=rho)$wi
    error=sum((prec.glasso%*%cov.valid-diag(rep(1,p)))^2)
    cv.rho.error=c(cv.rho.error,error)
    cv.rho=c(cv.rho,rho)
    print(i)
  }
  rho.i[valid.block]=cv.rho[which(cv.rho.error==min(cv.rho.error))]
}
rho.WS750=mean(rho.i) # =0.016
rho<-0.016

###### Task 3: window size 1000  ######
window_size = 1000
n<-dim(returnstd[1:window_size,])[1]
B=100
n.block=floor(n/B)
block.start=1+(0:(n.block-1))*B
# valid.block=sort(sample(1:n.block,floor(n.block/4)))
rho.i=c()
for (valid.block in 1:5) {
  valid.ind=NULL
  for(k in valid.block){
    valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
  }
  n.valid=length(valid.ind)
  train.ind=setdiff(1:n,valid.ind)
  n.train=length(train.ind)
  returnstd.train=returnstd[1:window_size,][train.ind,]
  returnstd.valid=returnstd[1:window_size,][valid.ind,]
  mu.train=colMeans(returnstd.train)
  mu.valid=colMeans(returnstd.valid)
  cov.train=cov(returnstd.train)
  cov.valid=cov(returnstd.valid)
  rho.grid=seq(0,0.8,length=101)[2:101]
  l.rho=length(rho.grid)
  cv.rho.error=NULL
  cv.rho=NULL
  for (i in 1:l.rho){
    rho=rho.grid[i]
    prec.glasso=glasso(cov.train,rho=rho)$wi
    error=sum((prec.glasso%*%cov.valid-diag(rep(1,p)))^2)
    cv.rho.error=c(cv.rho.error,error)
    cv.rho=c(cv.rho,rho)
    print(i)
  }
  rho.i[valid.block]=cv.rho[which(cv.rho.error==min(cv.rho.error))]
}
rho.WS1000=mean(rho.i) # =0.016
rho<-0.016

###### Save tuned parameters in glasso  ######

save(rho.WS500,rho.WS750,rho.WS1000,file="Glasso_rho.RData")

## Dantzig selector estimation for theta1, theta2, theta3, for window size 500/750/1000 #####
# theta1 <- Sigma^-1 1
# theta2 <- Sigma^-1 mu
# theta3 <- Sigma^-1 phi

###### Task 1: window size 500  ######
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
load("eigenvector_centrality_20230513.RData")

theta1<-list()
theta2<-list()
theta3<-list()
for(t in 1: length(W_in)){
  print(t)
  ptm<-proc.time()
  ## compute global minimum variance portfolio ##
  theta1[[t]] =linfun1(COV_in[[t]],rep(1,p),lambda=lmd1.WS500) # lambda <= 0.1 will lead to be infeasible
  # print('theta1')
  theta2[[t]] =linfun1(COV_in[[t]],ER_in[[t]],lambda=lmd2.WS500) # lambda <= 0.1 will lead to be infeasible
  # print('theta2')
  theta3[[t]] =linfun1(COV_in[[t]],EC_DS[[t]],lambda=lmd3.WS500) # lambda <= 0.1 will lead to be infeasible
  ptm<-proc.time()-ptm
  print(ptm)
}
theta<-list("theta1"=theta1,
            "theta2"=theta2,
            "theta3"=theta3)
save(theta,file="theta_Dantzig_20230513.RData")
load("theta_Dantzig_20230513.RData")
theta1=theta$theta1
theta2=theta$theta2
theta3=theta$theta3

###### Task 2: window size 750  ######
W<-list()
window_size = 750
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
load("eigenvector_centrality_WS750_20240723.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

theta1<-list()
theta2<-list()
theta3<-list()
for(t in 1: length(W_in)){
  print(t)
  ptm<-proc.time()
  ## compute global minimum variance portfolio ##
  theta1[[t]] =linfun1(COV_in[[t]],rep(1,p),lambda=lmd1.WS750) # lambda <= 0.1 will lead to be infeasible
  # print('theta1')
  theta2[[t]] =linfun1(COV_in[[t]],ER_in[[t]],lambda=lmd2.WS750) # lambda <= 0.1 will lead to be infeasible
  # print('theta2')
  theta3[[t]] =linfun1(COV_in[[t]],EC_DS[[t]],lambda=lmd3.WS750) # lambda <= 0.1 will lead to be infeasible
  ptm<-proc.time()-ptm
  print(ptm)
}
theta<-list("theta1"=theta1,
            "theta2"=theta2,
            "theta3"=theta3)
save(theta,file="theta_Dantzig_WS750_20240723.RData")
load("theta_Dantzig_WS750_20240723.RData")
theta1=theta$theta1
theta2=theta$theta2
theta3=theta$theta3

###### Task 3: window size 1000  ######
W<-list()
window_size = 1000
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
load("eigenvector_centrality_WS1000_20240723.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

theta1<-list()
theta2<-list()
theta3<-list()
for(t in 1: length(W_in)){
  print(t)
  ptm<-proc.time()
  ## compute global minimum variance portfolio ##
  theta1[[t]] =linfun1(COV_in[[t]],rep(1,p),lambda=lmd1.WS1000) # lambda <= 0.1 will lead to be infeasible
  # print('theta1')
  theta2[[t]] =linfun1(COV_in[[t]],ER_in[[t]],lambda=lmd2.WS1000) # lambda <= 0.1 will lead to be infeasible
  # print('theta2')
  theta3[[t]] =linfun1(COV_in[[t]],EC_DS[[t]],lambda=lmd3.WS1000) # lambda <= 0.1 will lead to be infeasible
  ptm<-proc.time()-ptm
  print(ptm)
}
theta<-list("theta1"=theta1,
            "theta2"=theta2,
            "theta3"=theta3)
save(theta,file="theta_Dantzig_WS1000_20240723.RData")
load("theta_Dantzig_WS1000_20240723.RData")
theta1=theta$theta1
theta2=theta$theta2
theta3=theta$theta3


