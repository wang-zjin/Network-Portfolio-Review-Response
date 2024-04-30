# The codes visualize two networks with the same absolute values of weights.  
# In Network 1, all weights are positive
# Network 2 also includes negative weights which are represented by blue lines. 
# Each vertex size represents its magnitude of eigenvector centrality.


rm(list = ls())
setwd("/Users/LvB/Documents/GitHub/Network-Portfolio-Review-Response/Figure_for_review_response")
# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')
library(xtable)

###### Network 1 with all positive weights ######
prices<-read.csv("SP500 securities_up_20230306.csv")
selected_columns <- c("Dates","MSFT","AKAM","NLOK","INTU","AAPL","EBAY","ORCL","CSCO")
prices1<-prices[,selected_columns]
ZOO <- zoo(prices1[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
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

C = C_in[[(1)]]
C[1,4] = C[4,1] = -C[1,4]
C[abs(C)<0.23] <- -2*C[abs(C)<0.23]
C = abs(C)
I = diag(rep(1,8))
A = C - I
print(round(A,2))

g <- graph_from_adjacency_matrix(C, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
EC = eigen_centrality(g,directed = FALSE, scale = TRUE)$vector
print(round(EC,2))
print(round(sum(EC),2))
Gamma = diag(eigen_centrality(g)$vector)
print(round(rowSums(A%*%Gamma),2))
print(round(A%*%EC,2))
print(round(sum(A%*%EC),2))

#figure
gd <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
eigen_centr = exp(evcent(g)$vector)
colors = rep(adjustcolor("Gray", alpha.f = .8), length(E(gd)))
colors = ifelse(E(gd)$weight < 0, 'blue', colors) 
png(paste0("Network","_1", ".png"), width = 750, height = 600, bg = "transparent")
plot(gd, layout = layout_in_circle, vertex.label = c(1:8), vertex.size = eigen_centr*25,
     edge.color = colors, edge.width = abs(E(g)$weight)^3*80,  vertex.label.cex = 2
)
dev.off()





###### Network 2 with positive and negative weights ######
prices<-read.csv("SP500 securities_up_20230306.csv")
selected_columns <- c("Dates","MSFT","AKAM","NLOK","INTU","AAPL","EBAY","ORCL","CSCO")
prices1<-prices[,selected_columns]
ZOO <- zoo(prices1[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
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

C = C_in[[(1)]]
C[1,4] = C[4,1] = -C[1,4]
C[abs(C)<0.23] <- -2*C[abs(C)<0.23]
I = diag(rep(1,8))
A = C - I
print(round(A,2))
# xtable= xtable(A)
# print(xtable, include.rownames = FALSE,inlcude.colnames = FALSE) 

g <- graph_from_adjacency_matrix(C, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
EC = eigen_centrality(g,directed = FALSE, scale = TRUE)$vector
print(round(EC,2))
print(round(sum(EC),2))
Gamma = diag(eigen_centrality(g)$vector)
print(round(rowSums(A%*%Gamma),2))
print(round(A%*%EC,2))
print(round(sum(A%*%EC),2))

gd <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
eigen_centr = exp(evcent(g)$vector)
colors = rep(adjustcolor("Gray", alpha.f = .8), length(E(gd)))
colors = ifelse(E(gd)$weight < 0, 'blue', colors) 
png(paste0("Network","_2", ".png"), width = 750, height = 600, bg = "transparent")
plot(gd, layout = layout_in_circle, vertex.label = c(1:8), vertex.size = eigen_centr*25,
     edge.color = colors, edge.width = abs(E(g)$weight)^3*80, vertex.label.cex = 2
)
dev.off()

