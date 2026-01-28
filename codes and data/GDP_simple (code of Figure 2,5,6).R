library(rmutil)  # for rlaplace, dlaplace
library(ggplot2)

# ========== EPSILON PRIME FINDING FUNCTIONS ==========

compute.integrals <- function(P_density, Q_density, epsilon_PQ, epsilon_QP, x_range) {
  integrand.PQ <- function(x) {
    p_val <- P_density(x)
    q_val <- Q_density(x)
    pmax(p_val - exp(epsilon_PQ) * q_val, 0)
  }
  
  integrand.QP <- function(x) {
    p_val <- P_density(x)
    q_val <- Q_density(x)
    pmax(q_val - exp(epsilon_QP) * p_val, 0)
  }
  
  integral_PQ <- integrate(integrand.PQ, lower = x_range[1], upper = x_range[2], 
                           subdivisions = 2000, rel.tol = 1e-10)$value
  integral_QP <- integrate(integrand.QP, lower = x_range[1], upper = x_range[2],
                           subdivisions = 2000, rel.tol = 1e-10)$value
  
  return(list(integral_PQ = integral_PQ, integral_QP = integral_QP))
}

compute.tau <- function(P_density, Q_density, epsilon, x_range) {
  result <- compute.integrals(P_density, Q_density, epsilon, epsilon, x_range)
  tau <- max(result$integral_PQ, result$integral_QP)
  return(list(tau = tau, integral_PQ = result$integral_PQ, integral_QP = result$integral_QP))
}

find.epsilon.prime <- function(P_density, Q_density, epsilon, x_range, tol = 1e-8, verbose = TRUE) {
  tau_result <- compute.tau(P_density, Q_density, epsilon, x_range)
  tau <- tau_result$tau
  
  if (verbose) {
    cat("Computing epsilon' for epsilon =", epsilon, "\n")
    cat("tau =", tau, "\n")
  }
  
  if (tau < 1e-10) {
    if (verbose) cat("tau is essentially zero\n")
    return(list(epsilon_prime = epsilon, tau = tau, 
                original_integrals = tau_result, final_integrals = tau_result))
  }
  
  if (tau_result$integral_PQ >= tau_result$integral_QP) {
    objective <- function(eps_prime) {
      result <- compute.integrals(P_density, Q_density, epsilon, eps_prime, x_range)
      abs(result$integral_QP - tau)
    }
    opt_result <- optimize(objective, interval = c(0, epsilon), tol = tol)
    epsilon_prime <- opt_result$minimum
    verify_result <- compute.integrals(P_density, Q_density, epsilon, epsilon_prime, x_range)
  } else {
    objective <- function(eps_prime) {
      result <- compute.integrals(P_density, Q_density, eps_prime, epsilon, x_range)
      abs(result$integral_PQ - tau)
    }
    opt_result <- optimize(objective, interval = c(0, epsilon), tol = tol)
    epsilon_prime <- opt_result$minimum
    verify_result <- compute.integrals(P_density, Q_density, epsilon_prime, epsilon, x_range)
  }
  
  if (verbose) {
    cat("Found epsilon' =", epsilon_prime, "\n")
    cat("Verification: PQ =", verify_result$integral_PQ, ", QP =", verify_result$integral_QP, "\n\n")
  }
  
  return(list(epsilon_prime = epsilon_prime, tau = tau, 
              original_integrals = tau_result, final_integrals = verify_result))
}

# ========== GDP Private Mean Estimation FUNCTIONS ==========

dp.quantile <- function(D, q, epsilon, alpha, a, b, T) {
  if (epsilon <= 0) stop("Privacy budget epsilon must be positive")
  if (alpha <= 2) stop("order must be larger than 2")
  if (q <= 0 || q >= 1) stop("Quantile q must be numeric in (0,1)")
  
  x <- sort(D)
  n <- length(x)
  
  left <- a
  right <- b
  
  # Fixed: Use T/epsilon^2 for variance, so sd = sqrt(T)/epsilon
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
    #a <- l_m - u_s * alpha * log(n) / 2
    a <- -u_s * log(n)^2  
  } else {
    a <- l_m
  }
  
  if (u_m == Inf) {
    #b <- u_m + u_s * alpha * log(n) / 2
    b <- u_s * log(n)^2  
  } else {
    b <- u_m
  }
  
  T <- ceiling(log2((b - a) * n^alpha))
  
  if (is.character(q)) {
    q <- tolower(q)
    tau <- min(sqrt(2 * T * log(T / n^(2-alpha))) / epsilon, n/4 - 1)  ### Avoid going too inward in the dataset
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
  
  # Privacy budget split
  ep_lowerQuantile <- dp_info$epsilon / log(n)^(1/4)
  ep_upperQuantile <- dp_info$epsilon / log(n)^(1/4)
  ep_meanEst <- dp_info$epsilon * sqrt(1 - 2 / log(n)^(1/2))  # Fixed: added sqrt
  
  # Calculate private quantiles
  dp_lowerQ = 0
  dp_upperQ = 0
  while (dp_lowerQ >= dp_upperQ) {
    dp_lowerQ <- dp.quantile.forTesting(D, dp_info$q[1], ep_lowerQuantile, dp_info$alpha, 
                                        dp_info$l_m, dp_info$u_m, dp_info$u_s)
    dp_upperQ <- dp.quantile.forTesting(D, dp_info$q[2], ep_upperQuantile, dp_info$alpha,
                                        dp_info$l_m, dp_info$u_m, dp_info$u_s)
  }
  
  # Clamp dataset
  data_clamped <- pmax(pmin(D, dp_upperQ), dp_lowerQ)
  
  # Add Gaussian noise
  z_noise <- rnorm(1, mean = 0, sd = (dp_upperQ - dp_lowerQ) / (n * ep_meanEst))
  
  return(mean(data_clamped) + z_noise)
}

# ========== LLR GDP TEST ==========

LLR.GDP.Test <- function(X, dp_info, P_density, Q_density) {
  n <- length(X)
  log_ratios <- log(P_density(X) / Q_density(X))
  
  test_statistic <- n*dp.meanEst(log_ratios, dp_info)
  
  return(list(
    test_statistic = test_statistic,
    log_ratios = log_ratios
  ))
}

# ========== ncLLR TEST ==========

clamp <- function(x, lower, upper) {
  pmax(lower, pmin(upper, x))
}

ncLLR <- function(X, P_density, Q_density, epsilon, epsilon_prime, add_noise = TRUE) {
  log_ratios <- log(P_density(X) / Q_density(X))
  clamped_log_ratios <- clamp(log_ratios, -epsilon_prime, epsilon)
  sum_clamped <- sum(clamped_log_ratios)
  
  if (add_noise) {
    gaussian_noise <- rnorm(1, 0, 2) # GDP guarantee for ncLLR
    test_statistic <- sum_clamped + gaussian_noise
  } else {
    test_statistic <- sum_clamped
  }
  
  return(list(
    test_statistic = test_statistic,
    sum_clamped = sum_clamped,
    log_ratios = log_ratios,
    clamped_log_ratios = clamped_log_ratios
  ))
}


# ========== Standard LLR ==========

LLR <- function(X, P_density, Q_density) {
  log_ratios <- log(P_density(X) / Q_density(X))
  test_statistic <- sum(log_ratios)
  
  return(list(
    test_statistic = test_statistic,
    log_ratios = log_ratios
  ))
}



# ========== Simulation ==========

# Define sample sizes and epsilon values to test
sample_sizes <- c(100, 200, 400, 800, 1600, 3200)
epsilon_values <- c(0.5, 1, 2)

# Initialize storage for results
results_all <- data.frame(
  epsilon = numeric(),
  n_samples = integer(),
  power_LLR = numeric(),
  power_ddcLLR = numeric(),
  power_ncLLR = numeric()
)

# Define densities
P_density <- function(x) dt(x, 1, 0)
Q_density <- function(x){ 0.5*dt(x, 1, 0) + 0.5*dt(x, 1.1, 0.1) }

x_range <- c(-100, 100)

# Loop over epsilon values
for (epsilon in epsilon_values) {
  cat("\n########################################\n")
  cat("####### EPSILON =", epsilon, "#######\n")
  cat("########################################\n")
  
  # Find epsilon prime for this epsilon value
  result <- find.epsilon.prime(P_density, Q_density, epsilon, x_range)
  cat("Epsilon prime found:", result$epsilon_prime, "\n")
  
  # Loop over sample sizes
  for (n_samples in sample_sizes) {
    cat("\n========================================\n")
    cat("Processing epsilon =", epsilon, ", n_samples =", n_samples, "\n")
    cat("========================================\n")
    
    # Setup for GDP test
    dp_info <- list(
      q = c("lower", "upper"),
      epsilon = epsilon,
      alpha = 3,
      l_m = -Inf,
      u_m = Inf,
      u_s = 10
    )
    
    # Step 1: Generate null distribution samples (1000 simulations)
    n_simulations_H0_samples <- 1000
    LLR_H0_samples <- numeric(n_simulations_H0_samples)
    GDP_H0_samples <- numeric(n_simulations_H0_samples)
    ncLLR_H0_samples <- numeric(n_simulations_H0_samples)
    
    cat("Generating null distribution with", n_simulations_H0_samples, "simulations...\n")
    pb <- txtProgressBar(min = 0, max = n_simulations_H0_samples, style = 3)
    set.seed(1234 + 100 * epsilon)
    
    for (i in 1:n_simulations_H0_samples) {
      # H0: Data from P
      X_from_P <- rt(n_samples, 1, 0)
      LLR_H0_samples[i] <- LLR(X_from_P, P_density, Q_density)$test_statistic
      GDP_H0_samples[i] <- LLR.GDP.Test(X_from_P, dp_info, P_density, Q_density)$test_statistic
      ncLLR_H0_samples[i] <- ncLLR(X_from_P, P_density, Q_density, epsilon, 
                                   result$epsilon_prime, add_noise = TRUE)$test_statistic
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Step 2: Calculate Power (1000 simulations under H1)
    n_simulations <- 1000
    LLR_H1 <- numeric(n_simulations)
    GDP_H1 <- numeric(n_simulations)
    ncLLR_H1 <- numeric(n_simulations)
    
    cat("\nCalculating Power with", n_simulations, "simulations...\n")
    pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
    set.seed(1234)
    
    for (i in 1:n_simulations) {
      # H1: Data from Q
      rQ <- function(n) {
        k <- rbinom(n, 1, 0.5) + 1
        rt(n, c(1.1, 1)[k], c(0.1, 0)[k])
      }
      X_from_Q <- rQ(n_samples)
      LLR_H1[i] <- LLR(X_from_Q, P_density, Q_density)$test_statistic
      GDP_H1[i] <- LLR.GDP.Test(X_from_Q, dp_info, P_density, Q_density)$test_statistic
      ncLLR_H1[i] <- ncLLR(X_from_Q, P_density, Q_density, epsilon, 
                           result$epsilon_prime, add_noise = TRUE)$test_statistic
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Calculate power using p-value formula: p = (1 + sum(T_H0 >= T_obs)) / (M + 1)
    # Reject if p-value <= 0.05
    power_LLR <- mean(sapply(LLR_H1, function(t_obs) {
      pvalue <- (1 + sum(LLR_H0_samples <= t_obs)) / (n_simulations_H0_samples + 1)
      pvalue <= 0.05
    }))
    
    power_ddcLLR <- mean(sapply(GDP_H1, function(t_obs) {
      pvalue <- (1 + sum(GDP_H0_samples <= t_obs)) / (n_simulations_H0_samples + 1)
      pvalue <= 0.05
    }))
    
    power_ncLLR <- mean(sapply(ncLLR_H1, function(t_obs) {
      pvalue <- (1 + sum(ncLLR_H0_samples <= t_obs)) / (n_simulations_H0_samples + 1)
      pvalue <= 0.05
    }))
    
    # Store results
    results_all <- rbind(results_all, data.frame(
      epsilon = epsilon,
      n_samples = n_samples,
      power_LLR = round(power_LLR, 3),
      power_ddcLLR = round(power_ddcLLR, 3),
      power_ncLLR = round(power_ncLLR, 3)
    ))
    
    # Print results for this combination
    cat("\n--- Results for epsilon =", epsilon, ", n_samples =", n_samples, "---\n")
    cat("LLR Power:", round(power_LLR, 3), "\n")
    cat("ddcLLR Power:", round(power_ddcLLR, 3), "\n")
    cat("ncLLR Power:", round(power_ncLLR, 3), "\n")
  }
}

# Print final summary table
cat("\n\n========================================\n")
cat("FINAL RESULTS SUMMARY\n")
cat("========================================\n")
print(results_all)

# Optional: Save results to CSV
write.csv(results_all, "simulation_results_all_methods.csv", row.names = FALSE)
cat("\nResults saved to 'simulation_results_all_methods.csv'\n")

# Create visualization comparing all three methods
library(tidyr)
library(dplyr)

# Reshape data for plotting
results_long <- results_all %>%
  pivot_longer(cols = starts_with("power_"), 
               names_to = "method", 
               values_to = "power") %>%
  mutate(method = case_when(
    method == "power_LLR" ~ "LLR (Non-private)",
    method == "power_ddcLLR" ~ "LLR.GDP.Test",
    method == "power_ncLLR" ~ "ncLLR"
  ))

# Create a combined plot with facets for all epsilon values
p_combined <- ggplot(results_long, 
                     aes(x = n_samples, y = power, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~epsilon, labeller = label_both) +
  scale_x_continuous(trans = "log2", breaks = sample_sizes) +
  labs(
    title = "Power Comparison Across Different Privacy Budgets",
    x = "Sample Size",
    y = "Power",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold")
  ) +
  ylim(0, 1)

print(p_combined)
ggsave("power_comparison_all_epsilon.png", p_combined, width = 14, height = 6)

# Optional: Create separate summary tables for each epsilon
for (eps in epsilon_values) {
  cat("\n--- Results for epsilon =", eps, "---\n")
  print(results_all[results_all$epsilon == eps, ])
}

###########################################################################
###########################################################################
###########################################################################

# CDF function for tulap
#   t - domain of CDF
#   median - median
#   lambda - Laplace parameter, as seen on Wikipedia
#   lcut - cut the leftmost lcut amount
#   rcut - cut the rightmost rcut amount
ptulap <- function (t, median = 0, lambda = 0, cut=0) {
  lcut=cut/2
  rcut=cut/2
  # Normalize
  t <- t - median
  
  # Split the positive and negative t calculations, and factor out stuff
  r <- round(t)
  g <- -log(lambda)
  l <- log(1 + lambda)
  k <- 1 - lambda
  negs <- exp((r * g) - l + log(lambda + ((t - r + (1/2)) * k)))
  poss <- 1 - exp((r * (-g)) - l + log(lambda + ((r - t + (1/2)) * k)))
  
  # Check for infinities
  negs[is.infinite(negs)] <- 0
  poss[is.infinite(poss)] <- 0
  
  # Truncate wrt the indicator on t's positivity
  is.leq0 <- t <= 0
  trunc <- (is.leq0 * negs) + ((1 - is.leq0) * poss)
  
  # Handle the cut adjustment and scaling
  cut <- lcut + rcut
  is.mid <- (lcut <= trunc) & (trunc <= (1 - rcut))
  is.rhs <- (1 - rcut) < trunc
  return (((trunc - lcut) / (1 - cut)) * is.mid + is.rhs)
}

library(rmutil)
library(dgof)

# One-sided KS test (greater alternative)
KS <- function(x, cdf, ep){
  n = length(x)
  # Use "greater" for one-sided test
  ks = ks.test(x, cdf, alternative = "less")$statistic
  
  # Add noise for differential privacy
  N = (1/n)*rnorm(n=1, mean = 0, sd = 1/ep)
  ks = ks + N
  return(ks)
}

Cramer <- function(x, cdf, ep) {
  U = sort(x)
  n = length(x)
  rank = seq_len(n)
  ###JA "U" should be replaced with cdf(U)
  omega2 = 1/(12 * n) + sum(((2*rank - 1)/(2*n)-cdf(U))^2)
  omega2 = sqrt(omega2/n)
  ### Add (1/n)*rnorm(1/ep)
  omega2 = omega2 + (1/n)*rnorm(1,0,1/ep)
  return(omega2)
}

cdf <- function(q){
  return (pnorm(q, m=0, s=1))
}

pdf <- function(x){
  return(dnorm(x, 0, 1))
}

# Reference distribution generation (no Kuiper for one-sided)
reference = function(n, ep, reps, type){
  # Generate data from laplacetic distribution
  x = rt(n*reps, 1, 0)
  x_mat = matrix(x, nrow=reps, ncol=n)
  
  if(type == "cramer")
    return(apply(X=x_mat, MARGIN=1, FUN=Cramer, cdf = pnorm, ep=ep))
  if(type == "ks")
    return(apply(X=x_mat, MARGIN=1, FUN=KS, cdf = pnorm, ep=ep))
}

# ========== Simulation ==========

set.seed(1234)

reps = 1000
al = 0.05
nVec = c(100, 200, 400, 800, 1600, 3200)
epVec = c(0.5, 1, 2)
mu1 = 0.1  # Fixed location shift to match first code

# Initialize results storage
results_KS = data.frame()

for (ep_idx in 1:length(epVec)) {
  ep = epVec[ep_idx]
  cat("\n########################################\n")
  cat("####### EPSILON =", ep, "#######\n")
  cat("########################################\n")
  
  p_GOF_KS = p_Cramer = rep(0, reps)
  power_GOF_KS = power_Cramer = rep(0, length(nVec))
  
  for (k in 1:length(nVec)){
    cat("\n========================================\n")
    cat("Processing epsilon =", ep, ", n_samples =", nVec[k], "\n")
    cat("========================================\n")
    
    n = nVec[k]
    
    # Generate null distributions
    cat("Generating null distribution with 1000 simulations...\n")
    null_cramer = ecdf(reference(n=n, ep=ep, reps=1000, type="cramer"))
    null_ks = ecdf(reference(n=n, ep=ep, reps=1000, type="ks"))
    
    cat("\nCalculating Power with", reps, "simulations...\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    
    for(i in 1:reps){
      # Generate data from laplacetic with location shift mu1 = 0.1
      # x = rlaplace(n, mu1, 1.1) 
      rQ <- function(n) {
        k <- rbinom(n, 1, 0.5) + 1
        rt(n, c(1.1, 1)[k], c(mu1, 0)[k])
      }
      x = rQ(n)
      
      # Calculate test statistics
      cramer = Cramer(x, cdf, ep) 
      ks = KS(x, cdf, ep)
      
      # Calculate p-values
      p_Cramer[i] = 1 - null_cramer(cramer)
      p_GOF_KS[i] = 1 - null_ks(ks)
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Calculate power
    power_GOF_KS[k] = mean(p_GOF_KS < al)
    power_Cramer[k] = mean(p_Cramer < al)
    
    # Store results for this combination
    results_KS = rbind(results_KS, data.frame(
      epsilon = ep,
      n_samples = n,
      power_KS = round(power_GOF_KS[k], 3),
      power_Cramer = round(power_Cramer[k], 3)
    ))
    
    # Print results for this combination
    cat("\n--- Results for epsilon =", ep, ", n_samples =", n, "---\n")
    cat("KS Power:", round(power_GOF_KS[k], 3), "\n")
    cat("Cramér Power:", round(power_Cramer[k], 3), "\n")
  }
}

# Print final summary table
cat("\n\n========================================\n")
cat("FINAL RESULTS SUMMARY\n")
cat("========================================\n")
print(results_KS)

# Save results
write.csv(results_KS, "ks_test_results_fixed_mu01.csv", row.names = FALSE)

# Optional: Create separate summary tables for each epsilon
for (eps in epVec) {
  cat("\n--- Results for epsilon =", eps, "---\n")
  print(results_KS[results_KS$epsilon == eps, ])
}

# ========================================
# VISUALIZATION
# ========================================
library(ggplot2)
library(tidyr)

# Reshape for ggplot
results_long = results_KS %>%
  pivot_longer(cols = starts_with("power_"), 
               names_to = "Test", 
               values_to = "Power",
               names_prefix = "power_")

# ========================================
# COMBINE WITH LLR/ddcLLR/ncLLR RESULTS 
# ========================================
if (exists("results_all")) {
  # Merge the dataframes
  results_combined = merge(results_all, results_KS, 
                           by = c("epsilon", "n_samples"),
                           all = TRUE)
  
  cat("\n\n========================================\n")
  cat("COMBINED RESULTS (LLR, ddcLLR, ncLLR, KS, Cramér)\n")
  cat("========================================\n")
  print(results_combined)
  
  write.csv(results_combined, "combined_all_methods_results.csv", row.names = FALSE)
  
  # Reshape for plotting all methods together
  results_all_long = results_combined %>%
    pivot_longer(cols = starts_with("power_"), 
                 names_to = "Test", 
                 values_to = "Power",
                 names_prefix = "power_")
  
  # Combined plot with all methods and epsilon values
  p_all_methods = ggplot(results_all_long, 
                         aes(x = n_samples, y = Power, color = Test, 
                             linetype = interaction(Test, factor(epsilon)), 
                             group = interaction(Test, epsilon))) +
    geom_line(size = 0.7) +
    geom_point(size = 1.3) +
    labs(title = paste0("All Methods Comparison Across All ε Values (μ₁ = ", mu1, ")"),
         x = "Sample Size",
         y = "Power") +
    theme_minimal() +
    scale_x_continuous(trans = "log10", breaks = nVec) +
    scale_color_manual(values = c("LLR" = "black",
                                  "ddcLLR" = "#E41A1C",
                                  "ncLLR" = "#FF7F00",
                                  "KS" = "#377EB8",
                                  "Cramer" = "#4DAF4A"),
                       labels = c("LLR" = "LLR (Non-private)",
                                  "ddcLLR" = "Ours",
                                  "ncLLR" = "ncLLR",
                                  "KS" = "KS Test",
                                  "Cramer" = "Cramér-von Mises"),
                       name = "Test") +
    scale_linetype_manual(values = c(
      # LLR - dashed for all epsilon values
      "LLR.0.5" = "dashed", "LLR.1" = "dashed", "LLR.2" = "dashed",
      # Others - solid
      "ddcLLR.0.5" = "solid", "ddcLLR.1" = "solid", "ddcLLR.2" = "solid",
      "ncLLR.0.5" = "solid", "ncLLR.1" = "solid", "ncLLR.2" = "solid",
      "KS.0.5" = "solid", "KS.1" = "solid", "KS.2" = "solid",
      "Cramer.0.5" = "solid", "Cramer.1" = "solid", "Cramer.2" = "solid"
    ),
    guide = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    ylim(0, 1)
  print(p_all_methods)
  
  # Faceted plot by epsilon
  p_faceted = ggplot(results_all_long, 
                     aes(x = n_samples, y = Power, color = Test, linetype = Test, group = Test)) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ epsilon, labeller = labeller(epsilon = function(x) paste0("ε = ", x)), ncol = 3) +
    labs(#title = paste0("All Methods Power Comparison by ε (θ1 = ", mu1, ")"),
      x = "Sample Size",
      y = "Power") +
    theme_minimal() +
    scale_x_continuous(trans = "log10", breaks = c(100, 200, 400, 800, 1600, 3200)) +
    scale_color_manual(values = c("LLR" = "black",
                                  "ddcLLR" = "#E41A1C",
                                  "ncLLR" = "#FF7F00",
                                  "KS" = "#377EB8",
                                  "Cramer" = "#4DAF4A"),
                       labels = c("LLR" = "LLR",
                                  "ddcLLR" = "Ours",
                                  "ncLLR" = "ncLLR",
                                  "KS" = "KS",
                                  "Cramer" = "CvM")) +
    scale_linetype_manual(values = c("LLR" = "dashed",
                                     "ddcLLR" = "solid",
                                     "ncLLR" = "solid",
                                     "KS" = "solid",
                                     "Cramer" = "solid"),
                          labels = c("LLR" = "LLR",
                                     "ddcLLR" = "Ours",
                                     "ncLLR" = "ncLLR",
                                     "KS" = "KS",
                                     "Cramer" = "CvM")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.box.margin = margin(t = -3),
          legend.margin = margin(t = 0),
          legend.box.spacing = unit(0.5, "lines"),
          legend.title = element_blank()) +
    ylim(0, 1) +
    theme(
      axis.text.x = element_text(angle = 39, hjust = 1, size = 15, face = "bold"),
      axis.text.y = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 18, face = "bold"),
      strip.text = element_text(size = 18, face = "bold"),
    )
  
  print(p_faceted)
  
  # Save the faceted plot
  ggsave("all_methods_faceted_comparison.png", p_faceted, width = 14, height = 6)
  cat("\nFaceted plot saved to 'all_methods_faceted_comparison.png'\n")
}

