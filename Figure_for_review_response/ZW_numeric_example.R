rm(list = ls())
setwd("~/Documents/GitHub/Network-Portfolio/Figure_for_review_response")
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

library(quadprog)

# Define the covariance matrix Sigma_mu (replace with actual values)
# Sigma_xi <- matrix(c(
#   1.0 , 0 , 0 , 0 , 0 , 0 , 0 , 0,
#   0 , 0.9 , 0 , 0 , 0 , 0 , 0 , 0,
#   0 , 0 , 0.9 , 0 , 0 , 0 , 0 , 0,
#   0 , 0 , 0 , 0.9 , 0 , 0 , 0 , 0,
#   0 , 0 , 0 , 0 , 1.2 , 0 , 0 , 0,
#   0 , 0 , 0 , 0 , 0 , 1.1 , 0 , 0,
#   0 , 0 , 0 , 0 , 0 , 0 , 0.8 , 0,
#   0 , 0 , 0 , 0 , 0 , 0 , 0 , 1.3
# ), nrow = 8, byrow = TRUE)
Sigma_xi <- matrix(c(
  1.0 , 0.5, 0.3, 0.2, 0.1, 0, 0, 0,
  0.5 , 0.9, 0.5, 0.3, 0.2, 0.1, 0, 0,
  0.3 , 0.5, 0.9, 0.5, 0.3, 0.2, 0.1, 0,
  0.2 , 0.3, 0.5, 0.9, 0.5, 0.3, 0.2, 0.1,
  0.1 , 0.2, 0.3, 0.5, 1.2, 0.5, 0.3, 0.2,
  0 ,   0.1, 0.2, 0.3, 0.5, 1.1, 0.5, 0.3,
  0 ,   0 ,  0.1, 0.2, 0.3, 0.5, 0.8, 0.5,
  0 ,   0 ,  0 ,  0.1, 0.2, 0.3, 0.5, 1.3
), nrow = 8, byrow = TRUE)

# Generate random samples from the normal distribution
set.seed(123)  # for reproducibility
n <- 1000  # length of time series
xi <- MASS::mvrnorm(n, mu = rep(0, 8), Sigma = Sigma_xi)

# xi is now a matrix with n rows and 8 columns, each row representing a sample from the 8-dimensional normal distribution

# # Define correlation matrix A 
# A <- matrix(c(
#   0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
#   0.2, 0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
#   0.3, 0.3, 0, 0.4, 0.5, 0.6, 0.7, 0.8,
#   0.4, 0.4, 0.4, 0, 0.5, 0.6, 0.7, 0.8,
#   0.5, 0.5, 0.5, 0.5, 0, 0.6, 0.7, 0.8,
#   0.6, 0.6, 0.6, 0.6, 0.6, 0, 0.7, 0.8,
#   0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0, 0.8,
#   0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0
# ), nrow = 8, byrow = TRUE)

# We can also consider the following A example
A = as.matrix(read.csv("Mat_A.csv"))
A[1,4] = A[4,1] = -A[1,4]
A[2,3] = A[3,2] = -A[2,3]
A[2,5] = A[5,2] = -A[2,5]
A[2,7] = A[7,2] = -A[2,7]
A[3,5] = A[5,3] = -A[3,5]
A[3,6] = A[6,3] = -A[3,6]
print(A)

# Define diagonal matrix Gamma 
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
eigen_centrality(g)$vector
Gamma <- diag(eigen_centrality(g)$vector)  

# # Calculate r: returns
# r <- matrix(0, n, 8)  # Initialize r matrix
# for (i in 1:n) {
#   r[i,] <- solve(diag(8) - A %*% Gamma) %*% xi[i,]
# }


# r1 = xi[i,]
# r2 = A %*% Gamma %*% xi[i,]
# r3 = (A %*% Gamma) %*% (A %*% Gamma) %*% xi[i,]
# r4 = (A %*% Gamma) %*% (A %*% Gamma) %*% (A %*% Gamma) %*% xi[i,]

# # Expression of sum( (A %*% Gamma)^i %*% xi[1,] )
# # Initialize a list to store the expressions
# expressions <- list()
# xi_expr <- "xi[1,]"
# # Define A * Gamma expression
# A_Gamma <- "A %*% Gamma "
# # Loop through the number of desired expressions
# for (k in 1:4) {
#   # Generate the expression by repeating the base expression k times and insert A * Gamma expression
#   expr <- paste(rep( A_Gamma, each = k), collapse = " %*% ")
#   expr <- paste(c(expr, xi_expr), collapse = " %*% ")
#   # Store the expression in the list
#   expressions[[paste0("r", k)]] <- expr
# }
# # Print the elegant expressions
# for (k in 1:4) {
#   cat("r", k, ":", expressions[[paste0("r", k)]], "\n")
# }


# r is now a matrix with n rows and 8 columns, each row representing a simulated vector r

##### Now we consider optimization problem (1)  #####
# min t(w)%*%Sigma_xi%*%w + t(w)%*%Q%*%w
# s.t. sum(w)=1
# is equivalent to
# min t(w)%*%(Sigma_xi+Q)%*%w
# s.t. sum(w)=1
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
print(w_1)
global_min_value = sol$value
print(2*global_min_value)
print(t(w_1)%*%Sigma_xi%*%w_1)
print(t(w_1)%*%Q%*%w_1)
print(t(w_1)%*%(Sigma_xi+Q)%*%w_1)
B_1 = t(w_1) %*% A %*% Gamma 
print(B_1)
print(B_1 %*% rep(1,8))

######## Optimization Problem (2)
# min t(w)%*%Sigma_xi%*%w
# s.t. sum(w)=1
#      t( A %*% Gamma %*%matrix(1, 8, 1) )%*%w <= M
M = 0.8
Dmat <- Sigma_xi  
dvec <- rep(0,8)  # Coefficients of the linear term in the objective function (zero in this case)
Amat <- cbind(matrix(1, 8, 1), -A %*% Gamma %*%matrix(1, 8, 1))  
bvec <- c(1,-M)  # Right-hand side of the constraint (1 and M in this case), lets set M = 100 for example
# Solve quadratic programming problem
sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
w_2 <- sol$solution
# print(w_2)
problem2_min_value = sol$value
# print(2*problem2_min_value)
print(t(w_2)%*%Sigma_xi%*%w_2)
print(t(w_2)%*%Q%*%w_2)
print(t(w_2)%*%(Sigma_xi+Q)%*%w_2)
B_2 = t(w_2) %*% A %*% Gamma 
# print(B_2)
print(B_2 %*% rep(1,8))
print(sum(abs(B_2)))


#### Problem 3
# min t(w)%*%Sigma_xi%*%w 
# s.t. sum(w)=1
#      t(w)%*%M%*%w <= M

