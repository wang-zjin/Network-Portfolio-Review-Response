# ===================================================
# Numeric study for network portfolio
# ===================================================

rm(list = ls())
library(quadprog)
library(clusterGeneration)
library(igraph)
library(stats)

# Set the random seed for reproducibility
#set.seed(2)
set.seed(15)
# Generate a random positive definite matrix
random_pd_matrix <- genPositiveDefMat(dim = 8)$Sigma
random_pd_matrix = abs(random_pd_matrix)
# Convert the positive definite matrix to a correlation matrix
A <- cov2cor(random_pd_matrix)
diag(A) = 0

# Define diagonal matrix Gamma 
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
Gamma <- diag(eigen_centrality(g)$vector)  
nc <- eigen_centrality(g)$vector
## Now we consider optimization problem (1)  
# min t(w)%*%Sigma_xi%*%w + t(w)%*%Q%*%w
# s.t. sum(w)=1
# is equivalent to
# min t(w)%*%(Sigma_xi+Q)%*%w
# s.t. sum(w)=1
Sigma_xi = random_pd_matrix
Q = A %*% Gamma %*% Sigma_xi %*% t(Gamma) %*% t(A)

# Define matrices
Dmat <- Sigma_xi + Q  # Dmat is the sum of Sigma_xi and Q
dvec <- rep(0,8)  # Coefficients of the linear term in the objective function (zero in this case)
Amat <- matrix(1, 8, 1)  # Constraint matrix for the sum of weights equal to 1
bvec <- 1  # Right-hand side of the constraint (1 in this case)
# Solve quadratic programming problem
sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)

# Extract solution vector
w_1 <- sol$solution
SS1 = t(w_1)%*%(Sigma_xi+Q)%*%w_1


## Now we consider network portfolio
phi = t(w_1)%*%nc
qtl=quantile(nc,seq(0.05,0.95,0.05))
N = 8

# lw = quantile(nc,0.35)
# hg = quantile(nc,0.65)
# sq = seq(lw,hg,by=0.001)

SS2 = rep(NA,length(qtl))
i=1
for ( phi in qtl ) {
  cov.mat = random_pd_matrix
  target.nc = phi
  Dmat = 2*cov.mat
  dvec = rep.int(0, N)
  Amat = cbind(rep(1,N), -nc) #, er
  bvec = c(1, -target.nc) #,target.er
  result = solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
  w = result$solution
  SS2[i] = t(w)%*%(Sigma_xi+Q)%*%w
  i=i+1
}
rets = cbind(qtl,SS2)

