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
#write.csv(A,file="Mat_A.csv",row.names = FALSE)
A = as.matrix(read.csv("Mat_A.csv"))
g <- graph_from_adjacency_matrix(C, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
EC = eigen_centrality(g,directed = FALSE, scale = TRUE)$vector
print(round(EC,2))
print(round(sum(EC),2))
Gamma = diag(eigen_centrality(g)$vector)
w1 = rep(0.125,8)
w1%*%Gamma%*%A
print(round(rowSums(w1%*%A%*%Gamma),3))
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
#write.csv(A,file="Mat_A2.csv",row.names = FALSE)
A = as.matrix(read.csv("Mat_A2.csv"))
print(round(A,2))

g <- graph_from_adjacency_matrix(C, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
EC = eigen_centrality(g,directed = FALSE, scale = TRUE)$vector
print(round(EC,2))
print(round(sum(EC),2))
Gamma = diag(eigen_centrality(g)$vector)
w1 = rep(0.125,8)
w1%*%A%*%Gamma
print(round(rowSums(w1%*%A%*%Gamma),3))
print(round(rowSums(A%*%Gamma),2))
print(round(A%*%EC,2))
print(round(sum(A%*%EC),2))
print(round(sum(abs(A%*%EC)),2))

gd <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
eigen_centr = exp(evcent(g)$vector)
colors = rep(adjustcolor("Gray", alpha.f = .8), length(E(gd)))
colors = ifelse(E(gd)$weight < 0, 'blue', colors) 
png(paste0("Network","_2", ".png"), width = 750, height = 600, bg = "transparent")
plot(gd, layout = layout_in_circle, vertex.label = c(1:8), vertex.size = eigen_centr*25,
     edge.color = colors, edge.width = abs(E(g)$weight)^3*80, vertex.label.cex = 2
)
dev.off()

