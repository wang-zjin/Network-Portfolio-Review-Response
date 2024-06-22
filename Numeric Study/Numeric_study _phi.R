# ===================================================
# Numeric study for network portfolio
# ===================================================

rm(list = ls())
library(quadprog)
library(clusterGeneration)
library(igraph)



  # Set the random seed for reproducibility
  set.seed(2*i+6)
  # Generate a random positive definite matrix
  random_pd_matrix <- genPositiveDefMat(dim = 8)$Sigma
  # Convert the positive definite matrix to a correlation matrix
  A <- cov2cor(random_pd_matrix)
  A = abs(A)
  diag(A) = 0
  
  # Define diagonal matrix Gamma 
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  eigen_centrality(g)$vector
  Gamma <- diag(eigen_centrality(g)$vector)  
  
  ## Now we consider network portfolio
  cov.mat = random_pd_matrix
  Dmat = 2*cov.mat
  N=8
  dvec = rep.int(0, N)
  Amat = cbind(rep(1,N), -nc, diag(1,N), er)
  bvec = c(1, -target.nc, rep(0,N), target.er)
  result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
  w = result$solution
  
  
  
  






