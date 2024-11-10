# Load necessary libraries
library(MASS)           # For mvrnorm (simulating returns)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(ggplot2)        # For visualization
library(CVglasso)       # For Glasso covariance estimation

# Set parameters
set.seed(123)
n_assets <- 50         # Number of assets
n_periods <- 20        # Number of periods (low sample case)

# Simulate returns data (50 assets, 20 periods)
mean_returns <- rep(0.001, n_assets)
cov_matrix <- matrix(0.02, n_assets, n_assets) + diag(0.05, n_assets)
returns <- mvrnorm(n = n_periods, mu = mean_returns, Sigma = cov_matrix)

# ---- Glasso Covariance Estimation ----
# Estimate the inverse covariance matrix using CVglasso
glasso_fit <- CVglasso(X = returns)
precision_matrix_glasso <- glasso_fit$Omega

# Eigen decomposition for Glasso (Inverse Covariance)
eigen_glasso <- eigen(precision_matrix_glasso)

# Sort eigenvalues in ascending order (smallest eigenvalue first)
sorted_indices_glasso <- order(eigen_glasso$values)
eigenvalues_glasso <- eigen_glasso$values[sorted_indices_glasso]
eigenvectors_glasso <- eigen_glasso$vectors[, sorted_indices_glasso]

# Invert the eigenvalues
inverted_eigenvalues_glasso <- 1 / eigenvalues_glasso

# Select eigenvectors corresponding to the smallest eigenvalues (last ones)
top_components_glasso <- eigenvectors_glasso[, 1:5]  # Select top 5 components

# ---- Linear Shrinkage Covariance ----
# Estimate the covariance matrix using linear shrinkage
cov_matrix_linshrink <- linshrink_cov(returns)

# Eigen decomposition for linear shrinkage covariance
eigen_linshrink <- eigen(cov_matrix_linshrink)
top_components_linshrink <- eigen_linshrink$vectors[, 1:5]  # Select top 5 components

# ---- Non-Linear Shrinkage Covariance ----
# Estimate the covariance matrix using non-linear shrinkage
cov_matrix_nlshrink <- nlshrink_cov(returns)

# Eigen decomposition for non-linear shrinkage covariance
eigen_nlshrink <- eigen(cov_matrix_nlshrink)
top_components_nlshrink <- eigen_nlshrink$vectors[, 1:5]  # Select top 5 components

# ---- Portfolio Construction Based on PCA Components ----
# Define function to calculate portfolio weights from top components
calculate_portfolio_weights <- function(top_components, inverted_eigenvalues = NULL) {
  if (!is.null(inverted_eigenvalues)) {
    # Adjust top components by the inverted eigenvalues (weighting the components)
    weights <- rowSums(top_components * inverted_eigenvalues[1:nrow(top_components)]) / sum(inverted_eigenvalues[1:nrow(top_components)])
  } else {
    # Simple equal weighting if no eigenvalue inversion is required
    weights <- rowSums(top_components) / sum(rowSums(top_components))
  }
  return(weights)
}

# Calculate portfolio weights for each method
portfolio_weights_glasso <- calculate_portfolio_weights(top_components_glasso, inverted_eigenvalues_glasso)
portfolio_weights_linshrink <- calculate_portfolio_weights(top_components_linshrink)
portfolio_weights_nlshrink <- calculate_portfolio_weights(top_components_nlshrink)

# ---- Portfolio Return Calculation ----
# Portfolio returns = weighted sum of asset returns
portfolio_return_glasso <- returns %*% portfolio_weights_glasso
portfolio_return_linshrink <- returns %*% portfolio_weights_linshrink
portfolio_return_nlshrink <- returns %*% portfolio_weights_nlshrink

# ---- Plot Portfolio Returns ----
portfolio_returns_data <- data.frame(
  Period = 1:n_periods,
  Glasso = portfolio_return_glasso,
  Linear_Shrinkage = portfolio_return_linshrink,
  Non_linear_Shrinkage = portfolio_return_nlshrink
)

# Plot portfolio returns
ggplot(portfolio_returns_data, aes(x = Period)) +
  geom_line(aes(y = Glasso, color = "Glasso")) +
  geom_line(aes(y = Linear_Shrinkage, color = "Linear Shrinkage")) +
  geom_line(aes(y = Non_linear_Shrinkage, color = "Non-linear Shrinkage")) +
  labs(title = "Portfolio Returns Comparison",
       x = "Period", y = "Portfolio Return") +
  theme_minimal() +
  scale_color_manual(values = c("Glasso" = "blue", 
                                "Linear Shrinkage" = "red", 
                                "Non-linear Shrinkage" = "green"))

#out sample variance
# Load necessary libraries
library(MASS)           # For mvrnorm (simulating returns)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(ggplot2)        # For visualization
library(CVglasso)       # For Glasso covariance estimation

# Set parameters
set.seed(123)
n_assets <- 50         # Number of assets
n_periods <- 20        # Number of periods (low sample case)

# Simulate returns data (50 assets, 20 periods)
mean_returns <- rep(0.001, n_assets)
cov_matrix <- matrix(0.02, n_assets, n_assets) + diag(0.05, n_assets)
returns <- mvrnorm(n = n_periods, mu = mean_returns, Sigma = cov_matrix)

# Split the data into in-sample (training) and out-of-sample (testing)
train_size <- floor(0.7 * n_periods)  # 70% for training
train_data <- returns[1:train_size, ]
test_data <- returns[(train_size + 1):n_periods, ]

# ---- Glasso Covariance Estimation ----
# Estimate the inverse covariance matrix using CVglasso
glasso_fit <- CVglasso(X = train_data)
precision_matrix_glasso <- glasso_fit$Omega

# Eigen decomposition for Glasso (Inverse Covariance)
eigen_glasso <- eigen(precision_matrix_glasso)

# Sort eigenvalues in ascending order (smallest eigenvalue first)
sorted_indices_glasso <- order(eigen_glasso$values)
eigenvalues_glasso <- eigen_glasso$values[sorted_indices_glasso]
eigenvectors_glasso <- eigen_glasso$vectors[, sorted_indices_glasso]

# Invert the eigenvalues
inverted_eigenvalues_glasso <- 1 / eigenvalues_glasso

# Select eigenvectors corresponding to the smallest eigenvalues (last ones)
top_components_glasso <- eigenvectors_glasso[, 1:5]  # Select top 5 components

# ---- Linear Shrinkage Covariance ----
# Estimate the covariance matrix using linear shrinkage
cov_matrix_linshrink <- linshrink_cov(train_data)

# Eigen decomposition for linear shrinkage covariance
eigen_linshrink <- eigen(cov_matrix_linshrink)
top_components_linshrink <- eigen_linshrink$vectors[, 1:5]  # Select top 5 components

# ---- Non-Linear Shrinkage Covariance ----
# Estimate the covariance matrix using non-linear shrinkage
cov_matrix_nlshrink <- nlshrink_cov(train_data)

# Eigen decomposition for non-linear shrinkage covariance
eigen_nlshrink <- eigen(cov_matrix_nlshrink)
top_components_nlshrink <- eigen_nlshrink$vectors[, 1:5]  # Select top 5 components

# ---- Portfolio Construction Based on PCA Components ----
# Define function to calculate portfolio weights from top components
calculate_portfolio_weights <- function(top_components, inverted_eigenvalues = NULL) {
  if (!is.null(inverted_eigenvalues)) {
    # Adjust top components by the inverted eigenvalues (weighting the components)
    weights <- rowSums(top_components * inverted_eigenvalues[1:nrow(top_components)]) / sum(inverted_eigenvalues[1:nrow(top_components)])
  } else {
    # Simple equal weighting if no eigenvalue inversion is required
    weights <- rowSums(top_components) / sum(rowSums(top_components))
  }
  return(weights)
}

# Calculate portfolio weights for each method
portfolio_weights_glasso <- calculate_portfolio_weights(top_components_glasso, inverted_eigenvalues_glasso)
portfolio_weights_linshrink <- calculate_portfolio_weights(top_components_linshrink)
portfolio_weights_nlshrink <- calculate_portfolio_weights(top_components_nlshrink)

# ---- Portfolio Return Calculation for Out-of-Sample Data ----
# Portfolio returns = weighted sum of asset returns (using out-of-sample data)
portfolio_return_glasso_out_of_sample <- test_data %*% portfolio_weights_glasso
portfolio_return_linshrink_out_of_sample <- test_data %*% portfolio_weights_linshrink
portfolio_return_nlshrink_out_of_sample <- test_data %*% portfolio_weights_nlshrink

# ---- Compute Out-of-Sample Portfolio Variance ----
# Function to compute out-of-sample portfolio variance
out_of_sample_variance <- function(returns, portfolio_weights) {
  portfolio_returns <- returns %*% portfolio_weights
  return(var(portfolio_returns))
}

# Compute out-of-sample variances
variance_glasso <- out_of_sample_variance(test_data, portfolio_weights_glasso)
variance_linshrink <- out_of_sample_variance(test_data, portfolio_weights_linshrink)
variance_nlshrink <- out_of_sample_variance(test_data, portfolio_weights_nlshrink)

# ---- Compute Reconstruction Error ----
# Function to compute reconstruction error (Mean Squared Error)
reconstruction_error <- function(original, reconstructed) {
  return(mean((original - reconstructed)^2))
}

# Reconstruct data using top components for each method
reconstructed_glasso <- top_components_glasso %*% t(top_components_glasso) %*% t(train_data)
reconstructed_linshrink <- top_components_linshrink %*% t(top_components_linshrink) %*% t(train_data)
reconstructed_nlshrink <- top_components_nlshrink %*% t(top_components_nlshrink) %*% t(train_data)

# Compute reconstruction errors
error_glasso <- reconstruction_error(train_data, t(reconstructed_glasso))
error_linshrink <- reconstruction_error(train_data, t(reconstructed_linshrink))
error_nlshrink <- reconstruction_error(train_data, t(reconstructed_nlshrink))

# ---- Plot Comparison of Reconstruction Errors and Out-of-Sample Variance ----
# Plot Reconstruction Error Comparison
ggplot(data.frame(Method = c("Glasso", "Linear Shrinkage", "Non-linear Shrinkage"),
                  Error = c(error_glasso, error_linshrink, error_nlshrink)),
       aes(x = Method, y = Error, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(title = "Reconstruction Error Comparison",
       x = "Method", y = "Reconstruction Error (MSE)") +
  theme_minimal()

# Plot Out-of-Sample Variance Comparison
ggplot(data.frame(Method = c("Glasso", "Linear Shrinkage", "Non-linear Shrinkage"),
                  Variance = c(variance_glasso, variance_linshrink, variance_nlshrink)),
       aes(x = Method, y = Variance, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(title = "Out-of-Sample Portfolio Variance Comparison",
       x = "Method", y = "Out-of-Sample Variance") +
  theme_minimal()

