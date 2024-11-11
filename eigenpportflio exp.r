# Load necessary libraries
library(MASS)           # For mvrnorm (simulating returns)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(ggplot2)        # For visualization
library(CVglasso)       # For Glasso covariance estimation

# Set parameters
set.seed(123)
n_assets <- 100         # Number of assets
sample_sizes <- c(20, 50, 80, 110, 130, 160)  # Varying sample sizes (periods)
n_simulations <- 100   # Number of simulations per sample size

# Simulate returns data
mean_returns <- rep(0.001, n_assets)
cov_matrix <- matrix(0.02, n_assets, n_assets) + diag(0.05, n_assets)

# ---- Function to Perform the Experiment ----
run_experiment <- function(n_periods) {
  # Simulate returns data for a given sample size
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
  
  # Invert the eigenvalues
  inverted_eigenvalues_glasso <- 1 / eigen_glasso$values
  
  # Select eigenvectors corresponding to the smallest eigenvalues (last ones)
  top_components_glasso <- eigen_glasso$vectors[, 1:5]  # Select top 5 components
  
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
  
  return(c(variance_glasso, variance_linshrink, variance_nlshrink))
}

# ---- Run the Experiment for Different Sample Sizes ----
results <- data.frame()
for (n_periods in sample_sizes) {
  # Store variances across simulations for each sample size
  variances_glasso <- numeric(n_simulations)
  variances_linshrink <- numeric(n_simulations)
  variances_nlshrink <- numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    # Run the experiment for each simulation
    variances <- run_experiment(n_periods)
    
    # Store the results for each method
    variances_glasso[i] <- variances[1]
    variances_linshrink[i] <- variances[2]
    variances_nlshrink[i] <- variances[3]
  }
  
  # Calculate average variance for each method
  avg_variance_glasso <- mean(variances_glasso)
  avg_variance_linshrink <- mean(variances_linshrink)
  avg_variance_nlshrink <- mean(variances_nlshrink)
  
  # Append to results
  results <- rbind(results, data.frame(SampleSize = n_periods,
                                       Glasso = avg_variance_glasso,
                                       LinearShrinkage = avg_variance_linshrink,
                                       NonLinearShrinkage = avg_variance_nlshrink))
}

# ---- Plot Average Variance Comparison ----
results_long <- reshape(results, 
                       varying = c("Glasso", "LinearShrinkage", "NonLinearShrinkage"), 
                       v.names = "Variance", 
                       timevar = "Method", 
                       times = c("Glasso", "LinearShrinkage", "NonLinearShrinkage"), 
                       direction = "long")

ggplot(results_long, aes(x = SampleSize, y = Variance, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Out-of-Sample Portfolio Variance Comparison",
       x = "Sample Size (Periods)", y = "Average Out-of-Sample Variance") +
  theme_minimal() +
  scale_color_manual(values = c("Glasso" = "blue", 
                                "LinearShrinkage" = "red", 
                                "NonLinearShrinkage" = "green"))


#with student-t distribution

# Load necessary libraries
library(MASS)           # For mvrnorm
library(mvtnorm)        # For rmvt (multivariate t-distribution)
library(nlshrink)       # For Non-linear Shrinkage Covariance
library(ggplot2)        # For visualization
library(CVglasso)       # For Glasso covariance estimation

# Set parameters
set.seed(123)
n_assets <- 50         # Number of assets
sample_sizes <- c(20, 50, 100)  # Varying sample sizes (periods)
n_simulations <- 100   # Number of simulations per sample size
df <- 5                # Degrees of freedom for t-distribution

# Simulate returns data using multivariate Student's t-distribution
mean_returns <- rep(0.001, n_assets)
cov_matrix <- matrix(0.02, n_assets, n_assets) + diag(0.05, n_assets)

# Function to simulate returns using multivariate Student's t-distribution
simulate_returns_t <- function(n, df, mu, Sigma) {
  rmvt(n = n, delta = mu, sigma = Sigma, df = df)
}

# ---- Function to Perform the Experiment ----
run_experiment <- function(n_periods) {
  # Simulate returns data for a given sample size using multivariate t-distribution
  returns <- simulate_returns_t(n = n_periods, df = df, mu = mean_returns, Sigma = cov_matrix)
  
  # Split the data into in-sample (training) and out-of-sample (testing)
  train_size <- floor(0.7 * n_periods)  # 70% for training
  train_data <- returns[1:train_size, ]
  test_data <- returns[(train_size + 1):n_periods, ]
  
  # ---- Glasso Covariance Estimation ----
  glasso_fit <- CVglasso(X = train_data)
  precision_matrix_glasso <- glasso_fit$Omega
  
  # Eigen decomposition for Glasso (Inverse Covariance)
  eigen_glasso <- eigen(precision_matrix_glasso)
  inverted_eigenvalues_glasso <- 1 / eigen_glasso$values
  top_components_glasso <- eigen_glasso$vectors[, 1:5]  # Select top 5 components
  
  # ---- Linear Shrinkage Covariance ----
  cov_matrix_linshrink <- linshrink_cov(train_data)
  eigen_linshrink <- eigen(cov_matrix_linshrink)
  top_components_linshrink <- eigen_linshrink$vectors[, 1:5]  # Select top 5 components
  
  # ---- Non-Linear Shrinkage Covariance ----
  cov_matrix_nlshrink <- nlshrink_cov(train_data)
  eigen_nlshrink <- eigen(cov_matrix_nlshrink)
  top_components_nlshrink <- eigen_nlshrink$vectors[, 1:5]  # Select top 5 components
  
  # ---- Portfolio Construction Based on PCA Components ----
  calculate_portfolio_weights <- function(top_components, inverted_eigenvalues = NULL) {
    if (!is.null(inverted_eigenvalues)) {
      weights <- rowSums(top_components * inverted_eigenvalues[1:nrow(top_components)]) / sum(inverted_eigenvalues[1:nrow(top_components)])
    } else {
      weights <- rowSums(top_components) / sum(rowSums(top_components))
    }
    return(weights)
  }
  
  # Calculate portfolio weights for each method
  portfolio_weights_glasso <- calculate_portfolio_weights(top_components_glasso, inverted_eigenvalues_glasso)
  portfolio_weights_linshrink <- calculate_portfolio_weights(top_components_linshrink)
  portfolio_weights_nlshrink <- calculate_portfolio_weights(top_components_nlshrink)
  
  # ---- Portfolio Return Calculation for Out-of-Sample Data ----
  portfolio_return_glasso_out_of_sample <- test_data %*% portfolio_weights_glasso
  portfolio_return_linshrink_out_of_sample <- test_data %*% portfolio_weights_linshrink
  portfolio_return_nlshrink_out_of_sample <- test_data %*% portfolio_weights_nlshrink
  
  # ---- Compute Out-of-Sample Portfolio Variance ----
  out_of_sample_variance <- function(returns, portfolio_weights) {
    portfolio_returns <- returns %*% portfolio_weights
    return(var(portfolio_returns))
  }
  
  # Compute out-of-sample variances
  variance_glasso <- out_of_sample_variance(test_data, portfolio_weights_glasso)
  variance_linshrink <- out_of_sample_variance(test_data, portfolio_weights_linshrink)
  variance_nlshrink <- out_of_sample_variance(test_data, portfolio_weights_nlshrink)
  
  return(c(variance_glasso, variance_linshrink, variance_nlshrink))
}

# ---- Run the Experiment for Different Sample Sizes ----
results <- data.frame()
for (n_periods in sample_sizes) {
  variances_glasso <- numeric(n_simulations)
  variances_linshrink <- numeric(n_simulations)
  variances_nlshrink <- numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    variances <- run_experiment(n_periods)
    variances_glasso[i] <- variances[1]
    variances_linshrink[i] <- variances[2]
    variances_nlshrink[i] <- variances[3]
  }
  
  avg_variance_glasso <- mean(variances_glasso)
  avg_variance_linshrink <- mean(variances_linshrink)
  avg_variance_nlshrink <- mean(variances_nlshrink)
  
  results <- rbind(results, data.frame(SampleSize = n_periods,
                                       Glasso = avg_variance_glasso,
                                       LinearShrinkage = avg_variance_linshrink,
                                       NonLinearShrinkage = avg_variance_nlshrink))
}

# ---- Plot Average Variance Comparison ----
results_long <- reshape(results, 
                       varying = c("Glasso", "LinearShrinkage", "NonLinearShrinkage"), 
                       v.names = "Variance", 
                       timevar = "Method", 
                       times = c("Glasso", "LinearShrinkage", "NonLinearShrinkage"), 
                       direction = "long")

ggplot(results_long, aes(x = SampleSize, y = Variance, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Out-of-Sample Portfolio Variance Comparison",
       x = "Sample Size (Periods)", y = "Average Out-of-Sample Variance") +
  theme_minimal() +
  scale_color_manual(values = c("Glasso" = "blue", 
                                "LinearShrinkage" = "red", 
                                "NonLinearShrinkage" = "green"))
