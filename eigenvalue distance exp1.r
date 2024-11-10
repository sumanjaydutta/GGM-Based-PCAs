# Load necessary libraries
library(MASS)           # For mvrnorm (simulating returns)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(CVglasso)       # For Glasso covariance estimation
library(ggplot2)        # For visualization

# Set seed for reproducibility
set.seed(123)

# Experiment parameters
n <- 100          # Number of observations
p <- 50           # Number of variables
num_simulations <- 100  # Number of simulations

# Generate the true covariance matrix and simulate data
true_cov_matrix <- diag(runif(p, 0.5, 2))  # True covariance matrix with random variances
true_data <- mvrnorm(n, mu = rep(0, p), Sigma = true_cov_matrix)

# Perform PCA on true covariance to get true eigenvalues
true_eigen <- eigen(true_cov_matrix)
true_eigenvalues <- sort(true_eigen$values)

# Initialize lists to store distances for each method across simulations
glasso_distances <- numeric(num_simulations)
nlshrink_distances <- numeric(num_simulations)

# Simulation loop
for (i in 1:num_simulations) {
  # Step 1: Generate training data based on true covariance matrix
  train_data <- mvrnorm(n, mu = rep(0, p), Sigma = true_cov_matrix)
  
  # Step 2: Covariance estimation using Glasso
  glasso_fit <- CVglasso(X = train_data)
  precision_matrix_glasso <- glasso_fit$Omega
  eigen_glasso <- eigen(precision_matrix_glasso)
  estimated_eigenvalues_glasso <- sort(1 / eigen_glasso$values)
  
  # Step 3: Covariance estimation using Non-linear Shrinkage
  cov_matrix_nlshrink <- nlshrink_cov(train_data)
  eigen_nlshrink <- eigen(cov_matrix_nlshrink)
  estimated_eigenvalues_nlshrink <- sort(eigen_nlshrink$values)
  
  # Step 4: Calculate the distance between true and estimated eigenvalues for each method
  glasso_distances[i] <- sqrt(sum((true_eigenvalues - estimated_eigenvalues_glasso)^2))
  nlshrink_distances[i] <- sqrt(sum((true_eigenvalues - estimated_eigenvalues_nlshrink)^2))
}

# Create a data frame for visualization
results_df <- data.frame(
  Method = rep(c("Glasso", "Non-linear Shrinkage"), each = num_simulations),
  Distance = c(glasso_distances, nlshrink_distances)
)

# Visualization of the distance distributions for each method
ggplot(results_df, aes(x = Method, y = Distance)) +
  geom_boxplot() +
  labs(title = "Eigenvalue Distance Comparison",
       y = "Distance between True and Estimated Eigenvalues") +
  theme_minimal()
  
  
  #incorporating linear shrinkage
  
  # Load necessary libraries
library(MASS)           # For mvrnorm (simulating returns)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(CVglasso)       # For Glasso covariance estimation
library(ggplot2)        # For visualization

# Set seed for reproducibility
set.seed(123)

# Experiment parameters
n <- 100          # Number of observations
p <- 50           # Number of variables
num_simulations <- 100  # Number of simulations

# Generate the true covariance matrix and simulate data
true_cov_matrix <- diag(runif(p, 0.5, 2))  # True covariance matrix with random variances
true_data <- mvrnorm(n, mu = rep(0, p), Sigma = true_cov_matrix)

# Perform PCA on true covariance to get true eigenvalues
true_eigen <- eigen(true_cov_matrix)
true_eigenvalues <- sort(true_eigen$values)

# Initialize lists to store distances for each method across simulations
glasso_distances <- numeric(num_simulations)
linshrink_distances <- numeric(num_simulations)
nlshrink_distances <- numeric(num_simulations)

# Simulation loop
for (i in 1:num_simulations) {
  # Step 1: Generate training data based on true covariance matrix
  train_data <- mvrnorm(n, mu = rep(0, p), Sigma = true_cov_matrix)
  
  # Step 2: Covariance estimation using Glasso (Inverse Covariance Estimation)
  glasso_fit <- CVglasso(X = train_data)
  precision_matrix_glasso <- glasso_fit$Omega
  eigen_glasso <- eigen(precision_matrix_glasso)
  estimated_eigenvalues_glasso <- sort(1 / eigen_glasso$values)
  glasso_distances[i] <- sqrt(sum((true_eigenvalues - estimated_eigenvalues_glasso)^2))
  
  # Step 3: Covariance estimation using Linear Shrinkage
  cov_matrix_linshrink <- linshrink_cov(train_data)
  eigen_linshrink <- eigen(cov_matrix_linshrink)
  estimated_eigenvalues_linshrink <- sort(eigen_linshrink$values)
  linshrink_distances[i] <- sqrt(sum((true_eigenvalues - estimated_eigenvalues_linshrink)^2))
  
  # Step 4: Covariance estimation using Non-linear Shrinkage
  cov_matrix_nlshrink <- nlshrink_cov(train_data)
  eigen_nlshrink <- eigen(cov_matrix_nlshrink)
  estimated_eigenvalues_nlshrink <- sort(eigen_nlshrink$values)
  nlshrink_distances[i] <- sqrt(sum((true_eigenvalues - estimated_eigenvalues_nlshrink)^2))
}

# Combine results into a data frame for visualization
results_df <- data.frame(
  Method = rep(c("Glasso", "Linear Shrinkage", "Non-linear Shrinkage"), each = num_simulations),
  Distance = c(glasso_distances, linshrink_distances, nlshrink_distances)
)

# Plot the results
ggplot(results_df, aes(x = Method, y = Distance)) +
  geom_boxplot() +
  labs(title = "Eigenvalue Distance Comparison", y = "Distance from True Eigenvalues") +
  theme_minimal()

