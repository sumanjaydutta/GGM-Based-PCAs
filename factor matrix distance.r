# Load necessary libraries
library(MASS)           # For mvrnorm (simulating returns)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(CVglasso)       # For Glasso covariance estimation
library(ggplot2)        # For visualization

# Set seed for reproducibility
set.seed(123)

# Experiment parameters
n <- 100                # Number of observations
p <- 50                 # Number of variables
num_simulations <- 100  # Number of simulations
num_top_components <- 5 # Number of top components (statistical factors)

# Generate the true covariance matrix and simulate data
true_cov_matrix <- diag(runif(p, 0.5, 2))  # True covariance matrix with random variances
true_data <- mvrnorm(n, mu = rep(0, p), Sigma = true_cov_matrix)

# Perform PCA on true covariance to get true factors (principal components)
true_eigen <- eigen(true_cov_matrix)
true_factors <- true_data %*% true_eigen$vectors[, order(true_eigen$values, decreasing = TRUE)[1:num_top_components]]

# Initialize lists to store Frobenius norms for each method across simulations
glasso_frob_distances <- numeric(num_simulations)
linshrink_frob_distances <- numeric(num_simulations)
nlshrink_frob_distances <- numeric(num_simulations)

# Simulation loop
for (i in 1:num_simulations) {
  # Step 1: Generate training data based on true covariance matrix
  train_data <- mvrnorm(n, mu = rep(0, p), Sigma = true_cov_matrix)
  
  # Step 2: Statistical factors estimation using Glasso
  glasso_fit <- CVglasso(X = train_data)
  precision_matrix_glasso <- glasso_fit$Omega
  eigen_glasso <- eigen(precision_matrix_glasso)
  
  # Invert eigenvalues and sort for covariance context
  estimated_eigenvectors_glasso <- eigen_glasso$vectors[, order(1 / eigen_glasso$values, decreasing = TRUE)[1:num_top_components]]
  estimated_factors_glasso <- train_data %*% estimated_eigenvectors_glasso
  glasso_frob_distances[i] <- norm(true_factors - estimated_factors_glasso, type = "F")
  
  # Step 3: Statistical factors estimation using Linear Shrinkage
  cov_matrix_linshrink <- linshrink_cov(train_data)
  eigen_linshrink <- eigen(cov_matrix_linshrink)
  estimated_eigenvectors_linshrink <- eigen_linshrink$vectors[, order(eigen_linshrink$values, decreasing = TRUE)[1:num_top_components]]
  estimated_factors_linshrink <- train_data %*% estimated_eigenvectors_linshrink
  linshrink_frob_distances[i] <- norm(true_factors - estimated_factors_linshrink, type = "F")
  
  # Step 4: Statistical factors estimation using Non-linear Shrinkage
  cov_matrix_nlshrink <- nlshrink_cov(train_data)
  eigen_nlshrink <- eigen(cov_matrix_nlshrink)
  estimated_eigenvectors_nlshrink <- eigen_nlshrink$vectors[, order(eigen_nlshrink$values, decreasing = TRUE)[1:num_top_components]]
  estimated_factors_nlshrink <- train_data %*% estimated_eigenvectors_nlshrink
  nlshrink_frob_distances[i] <- norm(true_factors - estimated_factors_nlshrink, type = "F")
}

# Combine results into a data frame for visualization
results_df <- data.frame(
  Method = rep(c("Glasso", "Linear Shrinkage", "Non-linear Shrinkage"), each = num_simulations),
  Frobenius_Distance = c(glasso_frob_distances, linshrink_frob_distances, nlshrink_frob_distances)
)

# Plot the results
ggplot(results_df, aes(x = Method, y = Frobenius_Distance)) +
  geom_boxplot() +
  labs(title = "Statistical Factors Frobenius Norm Comparison", y = "Frobenius Norm of Difference") +
  theme_minimal()
