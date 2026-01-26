# Comprehensive Comparison of Five Private Mean Estimation Methods
# Under GDP on One-Dimensional Data with Varying Sample Sizes
# 
# Methods compared:
# 0. Non-Private (Sample Mean) - Baseline
# 1. Shifted-Clipped-Mean (GDP) - Instance-optimal from Huang et al.
# 2. CoinPress (MVMRec_GDP) - Iterative confidence ball refinement
# 3. DP-MeanEst - Adaptive quantile-based method
# 4. Naive Data-Dependent - Simple baseline
#
# Test distributions: Gamma, Logistic, Gaussian
# Metric: Mean Squared Error (MSE)
# Sample sizes: 10^2, 10^2.5, 10^3, 10^3.5, 10^4

library(MASS)
library(ggplot2)

# ============================================================================
# LOAD ALL METHODS (same as before)
# ============================================================================

# --- Shifted-Clipped-Mean GDP (1D simplified version) ---
simple_clipped_mean_1d_gdp <- function(D, C, mu) {
  n <- length(D)
  D_clipped <- pmax(pmin(D, C), -C)
  empirical_mean <- mean(D_clipped)
  noise_scale <- (2 * C / n) / mu
  noise <- rnorm(1, mean = 0, sd = noise_scale)
  return(empirical_mean + noise)
}

simple_quantile_gdp <- function(D, q, mu) {
  n <- length(D)
  D_sorted <- sort(D)
  m <- ceiling(n * q)
  m <- max(1, min(m, n))
  empirical_q <- D_sorted[m]
  data_range <- diff(range(D))
  sensitivity <- data_range / sqrt(n)
  noise_scale <- sensitivity / mu
  noise <- rnorm(1, mean = 0, sd = noise_scale)
  return(empirical_q + noise)
}

shifted_clipped_mean_1d_simple <- function(D, mu, beta = 0.05) {
  n <- length(D)
  if (n < 10) stop("Need at least 10 samples")
  
  mu_shift <- mu / 3
  mu_clip <- mu / 3
  mu_mean <- mu / 3
  
  median_est <- simple_quantile_gdp(D, 0.5, mu_shift)
  D_shifted <- D - median_est
  abs_shifted <- abs(D_shifted)
  quantile_level <- 1 - 1/sqrt(n * mu_clip)
  quantile_level <- max(0.5, min(0.99, quantile_level))
  
  C <- simple_quantile_gdp(abs_shifted, quantile_level, mu_clip)
  C <- max(abs(C), 0.1)
  
  mean_shifted <- simple_clipped_mean_1d_gdp(D_shifted, C, mu_mean)
  result <- mean_shifted + median_est
  
  return(result)
}

shifted_clipped_mean_1d_gaussian_simple <- function(D, mu, sigma_guess = NULL, beta = 0.05) {
  n <- length(D)
  if (n < 10) stop("Need at least 10 samples")
  
  if (is.null(sigma_guess)) {
    sigma_guess <- mad(D, constant = 1.4826)
  }
  
  k <- sqrt(log(n))
  mu_center <- mu / 3
  mu_clip <- mu / 3
  mu_mean <- mu / 3
  
  center_est <- simple_quantile_gdp(D, 0.5, mu_center)
  D_shifted <- D - center_est
  C <- 3 * sigma_guess + k / sqrt(mu_clip)
  mean_shifted <- simple_clipped_mean_1d_gdp(D_shifted, C, mu_mean)
  result <- mean_shifted + center_est
  
  return(result)
}

# --- CoinPress GDP (Corrected Version) ---
# Helper: Gaussian mechanism
gaussian_mechanism <- function(value, sensitivity, epsilon) {
  # For GDP: noise std = sensitivity / epsilon
  noise_std <- sensitivity / epsilon
  noise <- rnorm(1, mean = 0, sd = noise_std)
  return(value + noise)
}

# Core iterative mean estimation
coinpress.multivariate_mean_iterative <- function(X, c, r, t, rho) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  n <- nrow(X)
  d <- ncol(X)
  current_mean <- c
  current_radius <- r
  
  for (i in 1:t) {
    # Center data around current mean
    X_centered <- sweep(X, 2, current_mean, "-")
    
    # Compute norms and clip to current radius
    norms <- sqrt(rowSums(X_centered^2))
    scale_factors <- pmin(1, current_radius / pmax(norms, 1e-10))
    X_clipped <- sweep(X_centered, 1, scale_factors, "*")
    X_clipped <- sweep(X_clipped, 2, current_mean, "+")
    
    # Add noise to mean
    sensitivity <- 2 * current_radius / n
    noisy_mean <- colMeans(X_clipped)
    for (j in 1:d) {
      noisy_mean[j] <- gaussian_mechanism(noisy_mean[j], sensitivity, rho[i])
    }
    
    current_mean <- noisy_mean
    
    # Update radius for next iteration
    if (i < t) {
      current_radius <- current_radius * sqrt(log(n) / n)
    }
  }
  
  return(current_mean)
}

# Main CoinPress mean estimation function
coinpress.meanEst <- function(D, epsilon, budget_split = c(0.1, 0.9), radius = NULL) {
  # CORRECTED: For GDP, epsilon_total^2 = epsilon_1^2 + epsilon_2^2 + ...
  # So if budget_split = c(0.25, 0.75), then:
  # epsilon_1 = sqrt(0.25) * epsilon
  # epsilon_2 = sqrt(0.75) * epsilon
  # Check: epsilon_1^2 + epsilon_2^2 = 0.25*epsilon^2 + 0.75*epsilon^2 = epsilon^2 ✓
  
  if (is.vector(D)) {
    D <- matrix(D, ncol = 1)
  }
  n <- nrow(D)
  d <- ncol(D)
  
  if (is.null(radius)) {
    radius <- 10 * sqrt(d)
  }
  center <- rep(0, d)
  
  # Normalize budget_split to sum to 1
  budget_split <- budget_split / sum(budget_split)
  
  # For GDP: rho[i] = sqrt(budget_split[i]) * epsilon
  rho <- sqrt(budget_split) * epsilon
  
  result <- coinpress.multivariate_mean_iterative(D, center, radius, length(rho), rho)
  if (d == 1) {
    return(as.numeric(result))
  }
  return(result)
}

# --- DP-MeanEst ---
dp.quantile <- function(D, q, epsilon, alpha, a, b, T) {
  if (epsilon <= 0) stop("Privacy budget epsilon must be positive")
  if (alpha <= 2) stop("order must be larger than 2")
  if (q <= 0 || q >= 1) stop("Quantile q must be numeric in (0,1)")
  
  x <- sort(D)
  n <- length(x)
  left <- a
  right <- b
  Z <- rnorm(T, mean = 0, sd = sqrt(T) / epsilon)
  mid <- (left + right) / 2
  
  for (t in 1:T) {
    true_count <- sum(x >= a & x <= mid)
    noisy_count <- true_count + Z[t]
    if (noisy_count < n * q) {
      left <- mid
    } else {
      right <- mid
    }
    mid <- (left + right) / 2
  }
  return(mid)
}

dp.quantile.forTesting <- function(D, q, epsilon, alpha, l_m, u_m, u_s) {
  x <- sort(D)
  n <- length(x)
  if (u_s <= 0) stop("scale must be positive")
  if (u_m - l_m <= 0) stop("search range must be positive")
  
  if (l_m == -Inf) {
    a <- -10*(log(n))^2 - 10
  } else {
    a <- l_m
  }
  
  if (u_m == Inf) {
    b <- 10*(log(n))^2 + 10
  } else {
    b <- u_m
  }
  
  T <- ceiling(log2((b - a) * n^alpha))
  
  if (is.character(q)) {
    q <- tolower(q)
    tau <- min(sqrt(2 * T * log(T / n^(2-alpha))) / epsilon, n/4 - 1)
    if (q == "lower") {
      q <- (tau + 2) / n
    } else if (q == "upper") {
      q <- 1 - ((tau + 1) / n)
    } else {
      stop("If q is character, it must be 'lower' or 'upper'")
    }
  }
  
  if (T < 1) stop("Number of steps T must be naturals")
  if (!is.numeric(q) || q <= 0 || q >= 1) {
    stop("Quantile q must be numeric in (0,1) or 'lower'/'upper'")
  }
  
  dp.quantile(D, q, epsilon, alpha, a, b, T)
}

dp.meanEst <- function(D, dp_info) {
  n <- length(D)
  ep_lowerQuantile <- dp_info$epsilon / log(n)^(1/4)
  ep_upperQuantile <- dp_info$epsilon / log(n)^(1/4)
  ep_meanEst <- dp_info$epsilon * sqrt(1 - 2 / log(n)^(1/2))
  
  dp_lowerQ = 0
  dp_upperQ = 0
  while (dp_lowerQ >= dp_upperQ) {
    dp_lowerQ <- dp.quantile.forTesting(D, dp_info$q[1], ep_lowerQuantile, dp_info$alpha, 
                                        dp_info$l_m, dp_info$u_m, dp_info$u_s)
    dp_upperQ <- dp.quantile.forTesting(D, dp_info$q[2], ep_upperQuantile, dp_info$alpha,
                                        dp_info$l_m, dp_info$u_m, dp_info$u_s)
  }
  
  data_clamped <- pmax(pmin(D, dp_upperQ), dp_lowerQ)
  z_noise <- rnorm(1, mean = 0, sd = (dp_upperQ - dp_lowerQ) / (n * ep_meanEst))
  return(mean(data_clamped) + z_noise)
}

# --- Naive Baseline ---
naive.datadep.mean <- function(D, epsilon) {
  n <- length(D)
  lower <- -10*(log(n))^2 - 10
  upper <- 10*(log(n))^2 + 10
  data_clamped <- pmax(pmin(D, upper), lower)
  sensitivity <- (upper - lower) / n
  z_noise <- rnorm(1, mean = 0, sd = sensitivity / epsilon)
  return(mean(data_clamped) + z_noise)
}

# ============================================================================
# WRAPPER FUNCTIONS
# ============================================================================

shifted_clipped_mean_1d <- function(D, mu, beta = 0.05) {
  result <- shifted_clipped_mean_1d_simple(D, mu, beta)
  return(as.numeric(result))
}

shifted_clipped_mean_1d_gaussian <- function(D, mu, R, sigma_min, sigma_max, beta = 0.05) {
  result <- shifted_clipped_mean_1d_gaussian_simple(D, mu, sigma_guess = sigma_min * 2, beta)
  return(as.numeric(result))
}

coinpress_1d <- function(D, mu_gdp, budget_split = c(0.1, 0.9), radius = NULL) {
  # D is 1D vector, mu_gdp is the GDP parameter (equivalent to epsilon in your notation)
  if (is.null(radius)) {
    # Use data-dependent radius
    data_range <- diff(range(D))
    radius <- max(10, data_range / 2)
  }
  
  result <- coinpress.meanEst(D, epsilon = mu_gdp, budget_split = budget_split, radius = radius)
  return(as.numeric(result))
}

dp_meanest_wrapper <- function(D, mu, alpha = 3) {
  epsilon <- mu 
  dp_info <- list(
    epsilon = epsilon,
    alpha = alpha,
    q = c("lower", "upper"),
    l_m = -Inf,
    u_m = Inf,
    u_s = 10
  )
  
  return(dp.meanEst(D, dp_info))
}

naive_wrapper <- function(D, mu) {
  epsilon <- mu
  return(naive.datadep.mean(D, epsilon))
}

# ============================================================================
# DATA GENERATION FUNCTIONS
# ============================================================================

generate_gamma_data <- function(n, shape = 2, rate = 0.5) {
  data <- rgamma(n, shape = shape, rate = rate)
  true_mean <- shape / rate
  return(list(data = data, true_mean = true_mean, 
              params = list(shape = shape, rate = rate)))
}

generate_logistic_data <- function(n, location = 5, scale = 2) {
  data <- rlogis(n, location = location, scale = scale)
  true_mean <- location
  return(list(data = data, true_mean = true_mean,
              params = list(location = location, scale = scale)))
}

generate_gaussian_data <- function(n, mean = 3, sd = 1) {
  data <- rnorm(n, mean = mean, sd = sd)
  true_mean <- mean
  return(list(data = data, true_mean = true_mean,
              params = list(mean = mean, sd = sd)))
}

# ============================================================================
# COMPARISON FRAMEWORK WITH VARYING SAMPLE SIZES
# ============================================================================

#' Run all four methods on a specific sample size and distribution
#' @param n sample size
#' @param params distribution parameters
#' @param mu_gdp GDP parameter
#' @param distribution_type "gamma", "logistic", or "gaussian"
#' @param n_trials number of independent trials
#' @return data frame with MSE results
compare_methods_single_n <- function(n, params, mu_gdp, distribution_type, n_trials = 50) {
  
  # Storage for squared errors
  se_shifted <- numeric(n_trials)
  se_coinpress <- numeric(n_trials)
  se_dpmeanest <- numeric(n_trials)
  se_naive <- numeric(n_trials)
  se_nonprivate <- numeric(n_trials)
  
  # Run trials
  for (trial in 1:n_trials) {
    
    # Generate fresh data for each trial
    if (distribution_type == "gamma") {
      data_obj <- generate_gamma_data(n, params$shape, params$rate)
    } else if (distribution_type == "logistic") {
      data_obj <- generate_logistic_data(n, params$location, params$scale)
    } else if (distribution_type == "gaussian") {
      data_obj <- generate_gaussian_data(n, params$mean, params$sd)
    }
    
    D <- data_obj$data
    true_mean <- data_obj$true_mean
    
    # Method 0: Non-private (baseline)
    tryCatch({
      est0 <- mean(D)
      se_nonprivate[trial] <- (est0 - true_mean)^2
    }, error = function(e) {
      se_nonprivate[trial] <- NA
    })
    
    # Method 1: Shifted-Clipped-Mean GDP
    tryCatch({
      if (distribution_type == "gaussian") {
        sigma <- params$sd
        R <- abs(true_mean) + 5 * sigma
        sigma_min <- sigma * 0.5
        sigma_max <- sigma * 2
        est1 <- shifted_clipped_mean_1d_gaussian(D, mu_gdp, R, sigma_min, sigma_max)
      } else {
        est1 <- shifted_clipped_mean_1d(D, mu_gdp)
      }
      se_shifted[trial] <- (est1 - true_mean)^2
    }, error = function(e) {
      se_shifted[trial] <- NA
    })
    
    # Method 2: CoinPress GDP
    tryCatch({
      # Use data-dependent radius
      data_range <- diff(range(D))
      radius <- max(10, data_range / 2)
      
      est2 <- coinpress_1d(D, mu_gdp, budget_split = c(0.1, 0.9), radius = radius)
      se_coinpress[trial] <- (est2 - true_mean)^2
    }, error = function(e) {
      se_coinpress[trial] <- NA
    })
    
    # Method 3: DP-MeanEst
    tryCatch({
      est3 <- dp_meanest_wrapper(D, mu_gdp)
      se_dpmeanest[trial] <- (est3 - true_mean)^2
    }, error = function(e) {
      se_dpmeanest[trial] <- NA
    })
    
    # Method 4: Naive Baseline
    tryCatch({
      est4 <- naive_wrapper(D, mu_gdp)
      se_naive[trial] <- (est4 - true_mean)^2
    }, error = function(e) {
      se_naive[trial] <- NA
    })
  }
  
  # Compute MSE for each method
  results <- data.frame(
    n = n,
    distribution = distribution_type,
    method = c("Non-Private", "Shifted-CM-GDP", "CoinPress-GDP", "DP-MeanEst", "Naive-DD"),
    mse = c(
      mean(se_nonprivate, na.rm = TRUE),
      mean(se_shifted, na.rm = TRUE),
      mean(se_coinpress, na.rm = TRUE),
      mean(se_dpmeanest, na.rm = TRUE),
      mean(se_naive, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# ============================================================================
# MAIN COMPARISON EXPERIMENT
# ============================================================================

run_full_comparison_varying_n <- function() {
  
  cat("\n")
  cat("========================================================================\n")
  cat("  COMPREHENSIVE COMPARISON WITH VARYING SAMPLE SIZES\n")
  cat("  Four Private Mean Estimation Methods Under GDP\n")
  cat("  Metric: Mean Squared Error (MSE)\n")
  cat("========================================================================\n")
  
  # Parameters
  mu_gdp <- 2.0        # GDP parameter
  n_trials <- 1000     # Number of trials per (n, distribution) combination
  
  # Sample sizes: 10^2, 10^2.5, 10^3, 10^3.5, 10^4
  n_values <- round(10^seq(2, 4, by = 0.5))
  
  cat("\nSample sizes to test:", paste(n_values, collapse = ", "), "\n")
  cat("Number of trials per configuration:", n_trials, "\n")
  cat("GDP parameter mu:", mu_gdp, "\n\n")
  
  # Set seed for reproducibility
  set.seed(234)
  
  # Define distributions
  distributions <- list(
    gamma = list(name = "Gamma", params = list(shape = 2, rate = 0.5)),
    logistic = list(name = "Logistic", params = list(location = 5, scale = 2)),
    gaussian = list(name = "Gaussian", params = list(mean = 3, sd = 1))
  )
  
  # Storage for all results
  all_results <- data.frame()
  
  # Run experiments
  for (dist_name in names(distributions)) {
    dist_config <- distributions[[dist_name]]
    
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("DISTRIBUTION:", toupper(dist_config$name), "\n")
    cat(rep("=", 70), "\n", sep = "")
    
    for (n in n_values) {
      cat(sprintf("  Running n = %d...", n))
      
      results <- compare_methods_single_n(
        n = n,
        params = dist_config$params,
        mu_gdp = mu_gdp,
        distribution_type = dist_name,
        n_trials = n_trials
      )
      
      all_results <- rbind(all_results, results)
      cat(" Done\n")
    }
    
    # Print summary table for this distribution
    cat("\n--- Summary Table (", dist_config$name, ") ---\n", sep = "")
    dist_results <- all_results[all_results$distribution == dist_name, ]
    
    # Reshape for better display
    summary_table <- reshape(dist_results[, c("n", "method", "mse")],
                             idvar = "n", timevar = "method", direction = "wide")
    colnames(summary_table) <- c("n", "Non-Private", "Shifted-CM-GDP", "CoinPress-GDP", 
                                 "DP-MeanEst", "Naive-DD")
    
    # Round for display
    summary_table[, -1] <- round(summary_table[, -1], 6)
    
    print(summary_table, row.names = FALSE)
    cat("\n")
  }
  
  return(all_results)
}

# ============================================================================
# VISUALIZATION
# ============================================================================

create_comparison_plot <- function(results) {
  
  # Rename distribution for better display
  results$Distribution <- factor(results$distribution,
                                 levels = c("gamma", "logistic", "gaussian"),
                                 labels = c("Gamma", "Logistic", "Gaussian"))
  
  # Rename methods for better display in legend
  results$method_display <- factor(results$method,
                                   levels = c("Non-Private", "Shifted-CM-GDP", "CoinPress-GDP", 
                                              "DP-MeanEst", "Naive-DD"),
                                   labels = c("Non-Private", "Shifted-CM", "CoinPress", 
                                              "Ours", "Naive-DD"))
  
  # Create the plot
  p <- ggplot(results, aes(x = n, y = mse, color = method_display, shape = method_display)) +
    geom_line(aes(linetype = method_display), size = 1) +
    geom_point(size = 3) +
    facet_wrap(~ Distribution, scales = "free_y", ncol = 3) +
    scale_x_log10(
      breaks = c(100, 1000, 10000),
      labels = c(expression(10^2), expression(10^3), expression(10^4)),
      limits = c(100, 10000)
    ) +
    scale_y_log10() +
    labs(
      x = "Sample Size (n)",
      y = "MSE (log scale)",
      color = "Method",
      shape = "Method",
      linetype = "Method",
      #title = "Comparison of Private Mean Estimation Methods (GDP)",
      #subtitle = "Mean Squared Error across different sample sizes and distributions"
    ) +
    scale_color_manual(
      values = c(
        "Non-Private" = "#000000",      # Black
        "Shifted-CM" = "#4DAF4A",       # Green
        "CoinPress" = "#377EB8",        # Blue
        "Ours" = "#E41A1C",             # Red
        "Naive-DD" = "#FF7F00"     # Orange
      )
    ) +
    scale_shape_manual(
      values = c(
        "Non-Private" = 8,        # Star/asterisk
        "Shifted-CM" = 16,        # Filled circle
        "CoinPress" = 17,         # Filled triangle
        "Ours" = 15,              # Filled square
        "Naive-DD" = 18      # Filled diamond
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Non-Private" = "dashed",
        "Shifted-CM" = "solid",
        "CoinPress" = "solid",
        "Ours" = "solid",
        "Naive-DD" = "solid"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(t = -3),
      legend.margin = margin(t = 0),
      legend.box.spacing = unit(0.5, "lines"),
      
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      strip.text = element_text(face = "bold", size = 18),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
      panel.grid.minor = element_line(color = "gray90", size = 0.3),
    )
  
  return(p)
}

# ============================================================================
# RUN EVERYTHING
# ============================================================================

# Run the full comparison
cat("\nStarting experiments...\n")
results <- run_full_comparison_varying_n()

cat("\n\n")
cat("========================================================================\n")
cat("  OVERALL RANKING ACROSS ALL CONFIGURATIONS\n")
cat("========================================================================\n\n")

# Compute average rank for each method across all (n, distribution) combinations
results$rank <- ave(results$mse, results$n, results$distribution, 
                    FUN = function(x) rank(x))

average_ranks <- aggregate(rank ~ method, data = results, FUN = mean)
average_ranks <- average_ranks[order(average_ranks$rank), ]

cat("Average Rank (1 = best):\n")
print(average_ranks, row.names = FALSE)

# Create and display the plot
cat("\n\nGenerating plot...\n")
p <- create_comparison_plot(results)
print(p)

# Save the plot
ggsave("mse_comparison_varying_n.png", p, width = 14, height = 5, dpi = 300)
cat("\nPlot saved as 'mse_comparison_varying_n.png'\n")

cat("\n========================================================================\n")
cat("  EXPERIMENT COMPLETE\n")
cat("========================================================================\n")

# Return results invisibly
invisible(results)