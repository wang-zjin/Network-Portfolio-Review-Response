# The codes visualize two networks with the same absolute values of weights.  
# In Network 1, all weights are positive
# Network 2 also includes negative weights which are represented by blue lines. 
# Each vertex size represents its magnitude of eigenvector centrality.


rm(list = ls())
# setwd("/Users/LvB/Documents/GitHub/Network-Portfolio-Review-Response/Figure_for_review_response")
setwd("~/Documents/GitHub/Network-Portfolio/Figure_for_review_response")
# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')
library(xtable)

###### Network 1 with all positive weights ######
#write.csv(A,file="Mat_A.csv",row.names = FALSE)
A = as.matrix(read.csv("Mat_A.csv"))
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
EC = eigen_centrality(g,directed = FALSE, scale = TRUE)$vector

print(round(EC,2))
print(round(sum(EC),2))
Gamma = diag(eigen_centrality(g)$vector)

w = matrix(repmat(1/8, 1, 8), nrow = 8)
print(round(t(w),3))
print(round(t(w) * rowSums(A%*%Gamma),3))
print(round(sum(t(w) * rowSums(A%*%Gamma)),3))

# w = matrix(c(-0.015, 0.145, 0.145, 0.145, 0.145, 0.145, 0.145, 0.145), nrow = 8)
w = matrix(c(0.145, 0.145, 0.145, 0.145, 0.145, 0.145, 0.145, -0.015), nrow = 8)
print(round(t(w),3))
print(round(t(w) * rowSums(A%*%Gamma),3))
print(round(sum(t(w) * rowSums(A%*%Gamma)),3))

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
A = as.matrix(read.csv("Mat_A.csv"))
A[1,4] = A[4,1] = -A[1,4]
A[2,3] = A[3,2] = -A[2,3]
A[2,5] = A[5,2] = -A[2,5]
A[2,7] = A[7,2] = -A[2,7]
A[3,5] = A[5,3] = -A[3,5]
A[3,6] = A[6,3] = -A[3,6]
A
# g <- graph_from_adjacency_matrix(C, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
EC = eigen_centrality(g,directed = FALSE, scale = TRUE)$vector

print(round(EC,2))
print(round(sum(EC),2))
Gamma = diag(eigen_centrality(g)$vector)

w = matrix(repmat(1/8, 1, 8), nrow = 8)
print(round(t(w),3))
print(round(t(w) * rowSums(A%*%Gamma),3))
print(round(sum(t(w) * rowSums(A%*%Gamma)),3))

w = matrix(c(0.145, 0.145, 0.145, 0.145, 0.145, 0.145, 0.145, -0.015), nrow = 8)
print(round(t(w),3))
print(round(t(w) * rowSums(A%*%Gamma),3))
print(round(sum(t(w) * rowSums(A%*%Gamma)),3))

gd <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
eigen_centr = exp(evcent(g)$vector)
colors = rep(adjustcolor("Gray", alpha.f = .8), length(E(gd)))
colors = ifelse(E(gd)$weight < 0, 'blue', colors) 
png(paste0("Network","_2", ".png"), width = 750, height = 600, bg = "transparent")
plot(gd, layout = layout_in_circle, vertex.label = c(1:8), vertex.size = eigen_centr*25,
     edge.color = colors, edge.width = abs(E(g)$weight)^3*80, vertex.label.cex = 2
)
dev.off()

