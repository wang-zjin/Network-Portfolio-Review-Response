# ===================================================
# Numeric study for optimization (11) and (12)
# M_setting
# ===================================================

rm(list = ls())
library(quadprog)
library(clusterGeneration)
library(igraph)
library(xtable)

# ###### M ######
# Set the random seed for reproducibility
set.seed(200)
# Generate a random positive definite matrix
random_pd_matrix <- genPositiveDefMat(dim = 8)$Sigma
# Convert the positive definite matrix to a correlation matrix
A <- cov2cor(random_pd_matrix)
random_pd_matrix = abs(random_pd_matrix)
diag(A) = 0

# Define diagonal matrix Gamma
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
Gamma <- diag(eigen_centrality(g)$vector)

## Now we consider optimization problem (1)
# min t(w)%*%Sigma_xi%*%w + t(w)%*%Q%*%w
# s.t. sum(w)=1
# is equivalent to
# min t(w)%*%(Sigma_xi+Q)%*%w
# s.t. sum(w)=1
Sigma_xi = random_pd_matrix
Q = A %*% Gamma %*% Sigma_xi %*% t(Gamma) %*% t(A)+ A %*% Gamma%*% Sigma_xi+Sigma_xi%*%t(Gamma)%*% t(A)


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
M =  t(w_1) %*%A%*% Gamma %*%matrix(1, 8, 1) 

Dmat <- Sigma_xi
dvec <- rep(0,8)  # Coefficients of the linear term in the objective function (zero in this case)
Amat <- cbind(matrix(1, 8, 1), -A %*% Gamma %*%matrix(1, 8, 1))
bvec <- c(1,-M)  # Right-hand side of the constraint (1 and M in this case), lets set M = 100 for example
# Solve quadratic programming problem
sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
w_2 <- sol$solution
problem2_min_value = sol$value
SS2 = t(w_2)%*%(Sigma_xi+Q)%*%w_2
SS2 -SS1

## Optimization Problem (2)
# min t(w)%*%Sigma_xi%*%w
# s.t. sum(w)=1
#      t(w_1) %*%A%*% Gamma %*%matrix(1, 8, 1) <= M

SS2 = rep(NA,length(seq(0.05,0.12,by=0.002)))
i=1
for ( M in seq(0.05,0.12,by=0.002) ){
  Dmat <- Sigma_xi
  dvec <- rep(0,8)  # Coefficients of the linear term in the objective function (zero in this case)
  Amat <- cbind(matrix(1, 8, 1), -A %*% Gamma %*%matrix(1, 8, 1))
  bvec <- c(1,-M)  # Right-hand side of the constraint (1 and M in this case), lets set M = 100 for example
  # Solve quadratic programming problem
  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  w_2 <- sol$solution
  problem2_min_value = sol$value
  SS2[i] = t(w_2)%*%(Sigma_xi+Q)%*%w_2
  print(M)
  print(t(w_2)%*%(Sigma_xi+Q)%*%w_2)
  i=i+1
}

M = seq(0.05,0.12,by=0.002)
rets = cbind(M,SS2)
M_1 = t(w_1) %*%A%*% Gamma %*%matrix(1, 8, 1)

#relative error
re = (SS2-rep(SS1,length(M)))/c(SS1)
ret = cbind(M,re)

mean((SS2-rep(SS1,length(M)))/c(SS1))

colnames(ret) = c("$M$", "RE")
xtable= xtable(ret,digits=c(0,3,6))
print(xtable, include.rownames = FALSE,sanitize.text.function = function(x) {x})

