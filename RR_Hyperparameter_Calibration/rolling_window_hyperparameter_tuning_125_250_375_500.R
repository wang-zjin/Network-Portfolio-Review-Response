rm(list = ls())

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

#### Calibrate the hyperparameter of Dantzig Selector for Eigenvector Centrality, using rolling window, with different window size: 125/250/375/500 #####
##### Task 1: window size 125 #####

window_size = 125
Dantzig_eigencent_rollwind_hypertune(returnstd, window_size = window_size, lambda_max=4)
# Tuning parameter of Dantzig estimation in Eigenvector centrality: 45187.047 sec elapsed

# Report error, because the range of lambda is not wide enough for window_size = 125, as window_size decreases,
# the lambda value (hyperparameter) should also increases.

##### Task 2: window size: 250 ######

window_size = 250
Dantzig_eigencent_rollwind_hypertune(returnstd, window_size = window_size, lambda_max=2)

##### Task 3: window size: 375 ######

window_size = 375
Dantzig_eigencent_rollwind_hypertune(returnstd, window_size = window_size)

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
lmd1.WS500=0.4636
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
lmd1.WS750 = 0.3592
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
lmd2.WS500=0.002142681
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
lmd2.WS750= 0.001682136
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
lmd3.WS500=0.4132
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
lmd3.WS750 =0.2782
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
rho.WS500=0.016
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
rho.WS750 =0.016
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
load("eigenvector_centrality_WS500_20240723.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

lmd1.WS500=0.4636
lmd2.WS500=0.002142681
lmd3.WS500=0.4132

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
save(theta,file="theta_Dantzig_WS500_20240723.RData")
load("theta_Dantzig_WS500_20240723.RData")
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

lmd1.WS750 = 0.3592
lmd2.WS750 = 0.001682136
lmd3.WS750 = 0.2782

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


