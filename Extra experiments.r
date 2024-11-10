# Load necessary libraries
library(glasso)   # For graphical lasso estimation of precision matrix
library(MASS)     # For generating sample data

# Step 1: Generate or load a dataset
set.seed(123)
n <- 100  # Number of observations
p <- 10   # Number of variables

# Generate synthetic data with a multivariate normal distribution
# Create a random covariance matrix
Sigma <- matrix(runif(p * p, min = -1, max = 1), p, p)
Sigma <- Sigma %*% t(Sigma)  # Make it symmetric positive-definite
data <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

# Step 2: Estimate the precision matrix using Graphical Lasso
lambda <- 0.1  # Regularization parameter for Graphical Lasso
glasso_result <- glasso(Sigma, rho = lambda)
precision_matrix <- glasso_result$wi  # Precision matrix (inverse covariance)

# Step 3: Compute eigenvalues and eigenvectors of the precision matrix
precision_eigen <- eigen(precision_matrix)

# Invert eigenvalues to get "principal component" ordering similar to standard PCA
eigenvalues_inv <- 1 / precision_eigen$values
eigenvectors <- precision_eigen$vectors

# Sort by descending order of inverse eigenvalues to align with PCA's convention
sorted_indices <- order(eigenvalues_inv, decreasing = TRUE)
precision_pcs <- eigenvectors[, sorted_indices]
precision_explained_variance <- eigenvalues_inv[sorted_indices] / sum(eigenvalues_inv)

# Step 4: Project the data onto the principal components
# Transform the original data to the precision-based PCA space
precision_pca_scores <- scale(data) %*% precision_pcs

# Step 5: Compute standard PCA for comparison
standard_pca <- prcomp(data, scale. = TRUE)
standard_pcs <- standard_pca$rotation  # Eigenvectors of the covariance matrix
standard_explained_variance <- standard_pca$sdev^2 / sum(standard_pca$sdev^2)  # Proportion of variance explained
standard_pca_scores <- standard_pca$x  # Standard PCA scores

# Step 6: Print and compare the results
cat("Precision-based PCA Principal Components:\n")
print(precision_pcs)
cat("\nPrecision-based PCA Explained Variance:\n")
print(precision_explained_variance)

cat("\nStandard PCA Principal Components:\n")
print(standard_pcs)
cat("\nStandard PCA Explained Variance:\n")
print(standard_explained_variance)

# Compare Explained Variance visually
barplot(rbind(precision_explained_variance, standard_explained_variance),
        beside = TRUE, col = c("blue", "red"),
        names.arg = paste("PC", 1:p),
        legend = c("Precision-based PCA", "Standard PCA"),
        main = "Comparison of Explained Variance: Precision-based vs Standard PCA",
        ylab = "Proportion of Variance Explained")

# Step 7: Summary of difference in scores (correlation between scores)
correlation_scores <- cor(precision_pca_scores[,1:min(p, ncol(precision_pca_scores))], 
                          standard_pca_scores[,1:min(p, ncol(standard_pca_scores))])
cat("\nCorrelation between Precision-based PCA and Standard PCA Scores:\n")
print(correlation_scores)



#low-sample application 

# Load necessary libraries
library(glasso)   # For graphical lasso estimation of precision matrix
library(MASS)     # For generating sample data

# Step 1: Generate or load a dataset with fewer samples than variables
set.seed(123)
n <- 50   # Number of observations (samples)
p <- 100  # Number of variables (high-dimensional setting)

# Generate synthetic data with a multivariate normal distribution
# Create a random covariance matrix for high-dimensional data
Sigma <- matrix(runif(p * p, min = -1, max = 1), p, p)
Sigma <- Sigma %*% t(Sigma)  # Make it symmetric positive-definite
data <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

# Step 2: Estimate the precision matrix using Graphical Lasso
lambda <- 0.1  # Regularization parameter for Graphical Lasso
glasso_result <- glasso(Sigma, rho = lambda)
precision_matrix <- glasso_result$wi  # Precision matrix (inverse covariance)

# Step 3: Compute eigenvalues and eigenvectors of the precision matrix
precision_eigen <- eigen(precision_matrix)

# Invert eigenvalues to get "principal component" ordering similar to standard PCA
eigenvalues_inv <- 1 / precision_eigen$values
eigenvectors <- precision_eigen$vectors

# Sort by descending order of inverse eigenvalues to align with PCA's convention
sorted_indices <- order(eigenvalues_inv, decreasing = TRUE)
precision_pcs <- eigenvectors[, sorted_indices]
precision_explained_variance <- eigenvalues_inv[sorted_indices] / sum(eigenvalues_inv)

# Step 4: Project the data onto the precision-based PCA components
precision_pca_scores <- scale(data) %*% precision_pcs

# Step 5: Compute standard PCA for comparison
standard_pca <- prcomp(data, scale. = TRUE)
standard_pcs <- standard_pca$rotation  # Eigenvectors of the covariance matrix
standard_explained_variance <- standard_pca$sdev^2 / sum(standard_pca$sdev^2)  # Proportion of variance explained
standard_pca_scores <- standard_pca$x  # Standard PCA scores

# Step 6: Print and compare the results
cat("Precision-based PCA Principal Components (Top 10):\n")
print(precision_pcs[, 1:10])
cat("\nPrecision-based PCA Explained Variance (Top 10):\n")
print(precision_explained_variance[1:10])

cat("\nStandard PCA Principal Components (Top 10):\n")
print(standard_pcs[, 1:10])
cat("\nStandard PCA Explained Variance (Top 10):\n")
print(standard_explained_variance[1:10])

# Compare Explained Variance visually (Top 10 components)
barplot(rbind(precision_explained_variance[1:10], standard_explained_variance[1:10]),
        beside = TRUE, col = c("blue", "red"),
        names.arg = paste("PC", 1:10),
        legend = c("Precision-based PCA", "Standard PCA"),
        main = "Comparison of Explained Variance (Top 10): Precision-based vs Standard PCA",
        ylab = "Proportion of Variance Explained")

# Step 7: Summary of difference in scores (correlation between scores)
correlation_scores <- cor(precision_pca_scores[, 1:min(10, ncol(precision_pca_scores))], 
                          standard_pca_scores[, 1:min(10, ncol(standard_pca_scores))])
cat("\nCorrelation between Precision-based PCA and Standard PCA Scores (Top 10):\n")
print(correlation_scores)


#portfolio based comparison (one synthetic sample set)
# Load necessary libraries
library(MASS)     # For generating multivariate normal data
library(glasso)   # For graphical lasso estimation of precision matrix

# Define number of assets and observations (high-dimensional case)
n_assets <- 100
n_obs <- 50

# Generate synthetic data for returns
set.seed(123)
Sigma <- matrix(runif(n_assets * n_assets, min = -1, max = 1), n_assets, n_assets)
Sigma <- Sigma %*% t(Sigma)  # Ensure positive-definiteness
returns_data <- mvrnorm(n_obs, mu = rep(0, n_assets), Sigma = Sigma)

# Apply Standard PCA
standard_pca <- prcomp(returns_data, scale. = TRUE)
standard_pca_components <- standard_pca$rotation
standard_explained_variance <- standard_pca$sdev^2 / sum(standard_pca$sdev^2)

# Apply Precision-based PCA using Graphical Lasso for precision matrix
glasso_result <- glasso(Sigma, rho = 0.1)
precision_matrix <- glasso_result$wi
precision_eigen <- eigen(precision_matrix)
precision_pca_components <- precision_eigen$vectors[, order(1 / precision_eigen$values, decreasing = TRUE)]
precision_explained_variance <- 1 / precision_eigen$values / sum(1 / precision_eigen$values)

# Construct minimum variance portfolios based on both PCA methods
standard_weights <- apply(standard_pca_components[, 1:5], 1, function(pc) sum(pc * returns_data) / sum(returns_data))
precision_weights <- apply(precision_pca_components[, 1:5], 1, function(pc) sum(pc * returns_data) / sum(returns_data))

# Calculate portfolio returns and risks
portfolio_return_standard <- colMeans(returns_data %*% standard_weights)
portfolio_risk_standard <- sqrt(var(returns_data %*% standard_weights))

portfolio_return_precision <- colMeans(returns_data %*% precision_weights)
portfolio_risk_precision <- sqrt(var(returns_data %*% precision_weights))

# Print results for comparison
cat("Standard PCA Portfolio Return:", portfolio_return_standard, "\n")
cat("Standard PCA Portfolio Risk:", portfolio_risk_standard, "\n")
cat("Precision-based PCA Portfolio Return:", portfolio_return_precision, "\n")
cat("Precision-based PCA Portfolio Risk:", portfolio_risk_precision, "\n")

#Risk Analysis (Volatility Estimation)
#In this experiment, we'll estimate portfolio volatility using both Standard PCA and Precision-based PCA. 
#We'll simulate a set of asset returns, apply both methods to extract principal components, and then compute portfolio volatility.

# Load necessary libraries
library(MASS)        # For mvrnorm
library(glasso)      # For graphical lasso (graphical Lasso package)
library(quantmod)    # For financial data

# Simulate returns (Example: 10 assets, 500 periods)
set.seed(123)
n_assets <- 10
n_periods <- 500
mean_returns <- rep(0.001, n_assets)    # Mean returns for each asset
cov_matrix <- matrix(0.02, n_assets, n_assets) + diag(0.05, n_assets) # Covariance matrix
returns <- mvrnorm(n = n_periods, mu = mean_returns, Sigma = cov_matrix)

# Calculate Standard PCA (using covariance matrix)
pca_standard <- prcomp(returns, scale = TRUE)

# Precision Matrix using Graphical Lasso (to estimate the inverse covariance)
cov_matrix_est <- cov(returns)
glasso_result <- glasso(cov_matrix_est, rho = 0.1)  # Estimate the precision matrix
precision_matrix <- glasso_result$wi  # Extract the precision matrix

# PCA based on the precision matrix
eigen_precision <- eigen(precision_matrix)
principal_components_precision <- eigen_precision$vectors

# Portfolio Risk Estimation using Principal Components
# Using first 3 PCs for both methods
w_standard <- pca_standard$rotation[, 1:3]  # Weights for Standard PCA
w_precision <- principal_components_precision[, 1:3]  # Weights for Precision-based PCA

# Portfolio variance (risk) for both methods
portfolio_variance_standard <- t(w_standard) %*% cov_matrix %*% w_standard
portfolio_variance_precision <- t(w_precision) %*% cov_matrix %*% w_precision

# Portfolio Volatility
portfolio_volatility_standard <- sqrt(portfolio_variance_standard)
portfolio_volatility_precision <- sqrt(portfolio_variance_precision)

# Print results
cat("Portfolio Volatility (Standard PCA): ", portfolio_volatility_standard, "\n")
cat("Portfolio Volatility (Precision-based PCA): ", portfolio_volatility_precision, "\n")


#Factor Exposure Analysis (Factor Loadings Comparison)
#In this experiment, we will compare the factor loadings extracted by Standard PCA and Precision-based PCA. 
#We'll look at how each method captures the underlying systematic risk factors.

# Factor loadings from Standard PCA
factor_loadings_standard <- pca_standard$rotation[, 1:3]  # First 3 factors

# Factor loadings from Precision-based PCA
factor_loadings_precision <- principal_components_precision[, 1:3]  # First 3 factors

# Compare factor loadings between methods
factor_comparison <- data.frame(
  Factor = 1:3,
  Standard_PCA = apply(factor_loadings_standard, 2, mean),
  Precision_PCA = apply(factor_loadings_precision, 2, mean)
)

# Print factor loadings comparison
print(factor_comparison)

#Portfolio Performance Comparison (Tracking Error)
#In this experiment, we will compute the tracking error for portfolios constructed using both Standard PCA and Precision-based PCA. 
#The tracking error measures how well a portfolio tracks a benchmark.

# Define a benchmark (e.g., a market index return, using synthetic data)
benchmark_returns <- rnorm(n_periods, mean = 0.001, sd = 0.02)

# Portfolio returns for Standard PCA and Precision-based PCA
portfolio_returns_standard <- returns %*% w_standard
portfolio_returns_precision <- returns %*% w_precision

# Compute Tracking Error (standard deviation of the difference between portfolio and benchmark returns)
tracking_error_standard <- sd(portfolio_returns_standard - benchmark_returns)
tracking_error_precision <- sd(portfolio_returns_precision - benchmark_returns)

# Print tracking errors
cat("Tracking Error (Standard PCA): ", tracking_error_standard, "\n")
cat("Tracking Error (Precision-based PCA): ", tracking_error_precision, "\n")




#In this experiment, we will generate data where the number of observations 
#n is smaller than the number of assets (variables)  p, and apply both Standard PCA and Precision-based PCA.
#We will compare how stable or reliable the principal components are as the sample size decreases. This is typically referred to as the "small sample problem" in high-dimensional statistics.

library(glasso)               # For graphical lasso
library(MASS)                 # For mvrnorm
library(stats)                # For PCA

# Function to simulate data with n < p
simulate_data <- function(n, p) {
  # Simulate n observations of p assets (variables)
  mean_returns <- rep(0.001, p)
  cov_matrix <- matrix(0.02, p, p) + diag(0.05, p)
  returns <- mvrnorm(n = n, mu = mean_returns, Sigma = cov_matrix)
  return(returns)
}

# Function to run both PCA methods and compare stability
run_pca_comparison <- function(n, p) {
  # Simulate returns with n samples and p assets
  returns <- simulate_data(n, p)
  
  # Standard PCA
  pca_standard <- prcomp(returns, scale = TRUE)
  standard_pc1 <- pca_standard$x[, 1]  # First principal component

  # Graphical Lasso to estimate precision matrix
  cov_matrix <- cov(returns)
  glasso_result <- glasso(cov_matrix, rho = 0.1)  # Regularization parameter
  precision_matrix <- glasso_result$wi
  
  # PCA based on precision matrix (eigen decomposition of precision matrix)
  eigen_precision <- eigen(precision_matrix)
  precision_pc1 <- eigen_precision$vectors[, 1]  # First principal component from precision matrix
  
  # Return the results
  return(list(standard_pc1 = standard_pc1, precision_pc1 = precision_pc1))
}

# Simulating for different sample sizes and analyzing stability
n_values <- c(10, 20, 50, 100)  # Different sample sizes (n < p)
p <- 200  # Number of assets (p)

# Store results for comparison
results_standard <- list()
results_precision <- list()

for (n in n_values) {
  cat("Running PCA for sample size n =", n, "\n")
  results <- run_pca_comparison(n, p)
  
  results_standard[[as.character(n)]] <- results$standard_pc1
  results_precision[[as.character(n)]] <- results$precision_pc1
}

# Analyzing the stability of principal components by calculating pairwise correlation
# between components at different sample sizes
stability_standard <- sapply(n_values, function(n) {
  if (n > 10) {
    # Compare the first PC across different sample sizes for Standard PCA
    cor(results_standard[[as.character(n)]], results_standard[[as.character(n-10)]], use = "pairwise.complete.obs")
  } else {
    return(NA)
  }
})

stability_precision <- sapply(n_values, function(n) {
  if (n > 10) {
    # Compare the first PC across different sample sizes for Precision-based PCA
    cor(results_precision[[as.character(n)]], results_precision[[as.character(n-10)]], use = "pairwise.complete.obs")
  } else {
    return(NA)
  }
})

# Print out the stability analysis
cat("Stability of Standard PCA Components: \n")
print(stability_standard)

cat("Stability of Precision-based PCA Components: \n")
print(stability_precision)


#experimenting with shrinkage

# Load necessary libraries
library(MASS)        # For mvrnorm (simulating returns)
library(glasso)      # For graphical lasso (precision matrix estimation)
library(nlshrink)    # For non-linear shrinkage
library(quantmod)    # For financial data

# Simulate returns data (50 assets, 30 periods)
set.seed(123)
n_assets <- 50     # Number of assets
n_periods <- 30    # Number of periods (samples)
mean_returns <- rep(0.001, n_assets)
cov_matrix <- matrix(0.02, n_assets, n_assets) + diag(0.05, n_assets) 
returns <- mvrnorm(n = n_periods, mu = mean_returns, Sigma = cov_matrix)

# ---- Precision-based PCA (using Graphical Lasso for precision matrix estimation) ----
# Estimate precision matrix using glasso (graphical lasso)
cov_matrix_sample <- cov(returns)
glasso_fit <- glasso(cov_matrix_sample, rho = 0.1)  # rho controls the regularization strength
precision_matrix <- glasso_fit$wi  # The estimated precision (inverse covariance) matrix

# Eigen decomposition of the precision matrix
eigen_precision <- eigen(precision_matrix)  
principal_components_precision <- eigen_precision$vectors[, 1]  # First principal component (PC1)

# ---- PCA using Non-Linear Shrinkage for Covariance Matrix ----
cov_shrinked <- nlshrink_cov(returns)  # Non-linear shrinkage estimator for covariance matrix
cov_matrix_shrinked <- cov_shrinked  # Shrinked covariance matrix
eigen_shrinked <- eigen(cov_matrix_shrinked)  # Eigen decomposition of shrinked covariance matrix
principal_components_shrinked <- eigen_shrinked$vectors[, 1]  # First principal component (PC1)

# ---- Portfolio Weights and Performance ----

# Portfolio weights (using first principal component)
weights_precision <- principal_components_precision / sum(abs(principal_components_precision))  # Normalize to sum to 1
weights_shrinked <- principal_components_shrinked / sum(abs(principal_components_shrinked))  # Normalize to sum to 1

# Portfolio variance (risk)
portfolio_variance_precision <- t(weights_precision) %*% cov_matrix %*% weights_precision
portfolio_variance_shrinked <- t(weights_shrinked) %*% cov_matrix %*% weights_shrinked

# Portfolio volatility (standard deviation)
portfolio_volatility_precision <- sqrt(portfolio_variance_precision)
portfolio_volatility_shrinked <- sqrt(portfolio_variance_shrinked)

# Portfolio return (mean of weighted returns)
portfolio_return_precision <- mean(returns %*% weights_precision)
portfolio_return_shrinked <- mean(returns %*% weights_shrinked)

# ---- Results ----
cat("Portfolio Return (Precision-based PCA): ", portfolio_return_precision, "\n")
cat("Portfolio Volatility (Precision-based PCA): ", portfolio_volatility_precision, "\n")
cat("Portfolio Return (Non-linear Shrinkage PCA): ", portfolio_return_shrinked, "\n")
cat("Portfolio Volatility (Non-linear Shrinkage PCA): ", portfolio_volatility_shrinked, "\n")

# Compare the results
# Higher portfolio return and lower volatility would indicate better performance


