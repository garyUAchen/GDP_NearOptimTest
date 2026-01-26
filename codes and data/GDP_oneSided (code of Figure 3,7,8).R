library(rmutil)
library(ggplot2)
library(reshape2)

# ========== GDP Private Mean Estimation FUNCTIONS ==========

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
    a <- -u_s*log(n)^2  
  } else {
    a <- l_m
  }
  
  if (u_m == Inf) {
    b <- u_s*log(n)^2  
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

# ========== LLR GDP TEST ==========

LLR.GDP.Test <- function(X, dp_info, P_density, Q_density) {
  n <- length(X)
  log_ratios <- log(P_density(X) / Q_density(X))
  
  test_statistic <- n * dp.meanEst(log_ratios, dp_info)
  
  return(list(
    test_statistic = test_statistic,
    log_ratios = log_ratios
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

# Define parameters
sample_sizes <- c(100, 200, 400, 800, 1600, 3200)
epsilon_values <- c(0.5, 1, 2)
mu1_values <- c(0, 0.03, 0.05, 0.1, 0.2, 0.4)

# Initialize storage for results
results_all <- data.frame(
  epsilon = numeric(),
  n_samples = integer(),
  mu1 = numeric(),
  power_LLR = numeric(),
  power_ddcLLR = numeric()
)

# Loop over epsilon values
for (epsilon in epsilon_values) {
  cat("\n########################################\n")
  cat("####### EPSILON =", epsilon, "#######\n")
  cat("########################################\n")
  
  # Loop over mu1 values (alternative hypothesis mean)
  for (mu1 in mu1_values) {
    cat("\n========================================\n")
    cat("Processing epsilon =", epsilon, ", mu1 =", mu1, "\n")
    cat("========================================\n")
    
    # Define densities for the log likelihood ratio test statistic
    P_density <- function(x) dnorm(x, 0, 1)
    Q_density <- function(x) dnorm(x, 0.02, 1)
    
    # Loop over sample sizes
    for (n_samples in sample_sizes) {
      cat("\nProcessing n_samples =", n_samples, "\n")
      
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
      
      cat("Generating null distribution with", n_simulations_H0_samples, "simulations...\n")
      pb <- txtProgressBar(min = 0, max = n_simulations_H0_samples, style = 3)
      set.seed(12345 + 1000 * epsilon + 100 * mu1 + n_samples)
      
      for (i in 1:n_simulations_H0_samples) {
        # H0: Data from P (mean = 0)
        X_from_P <- rnorm(n_samples, 0, 1)
        LLR_H0_samples[i] <- LLR(X_from_P, P_density, Q_density)$test_statistic
        GDP_H0_samples[i] <- LLR.GDP.Test(X_from_P, dp_info, P_density, Q_density)$test_statistic
        
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # Step 2: Calculate Power (1000 simulations under H1)
      n_simulations <- 1000
      LLR_H1 <- numeric(n_simulations)
      GDP_H1 <- numeric(n_simulations)
      
      cat("\nCalculating Power with", n_simulations, "simulations...\n")
      pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
      set.seed(54321 + 1000 * epsilon + 100 * mu1 + n_samples)
      
      for (i in 1:n_simulations) {
        # H1: Data from Q (mean = mu1)
        X_from_Q <- rnorm(n_samples, mu1, 1)
        LLR_H1[i] <- LLR(X_from_Q, P_density, Q_density)$test_statistic
        GDP_H1[i] <- LLR.GDP.Test(X_from_Q, dp_info, P_density, Q_density)$test_statistic
        
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # Calculate power using p-value formula
      power_LLR <- mean(sapply(LLR_H1, function(t_obs) {
        pvalue <- (1 + sum(LLR_H0_samples <= t_obs)) / (n_simulations_H0_samples + 1)
        pvalue <= 0.05
      }))
      
      power_ddcLLR <- mean(sapply(GDP_H1, function(t_obs) {
        pvalue <- (1 + sum(GDP_H0_samples <= t_obs)) / (n_simulations_H0_samples + 1)
        pvalue <= 0.05
      }))
      
      # Store results
      results_all <- rbind(results_all, data.frame(
        epsilon = epsilon,
        n_samples = n_samples,
        mu1 = mu1,
        power_LLR = round(power_LLR, 3),
        power_ddcLLR = round(power_ddcLLR, 3)
      ))
      
      cat("LLR Power:", round(power_LLR, 3), "\n")
      cat("ddcLLR Power:", round(power_ddcLLR, 3), "\n")
    }
  }
}

# Print final summary table
cat("\n\n========================================\n")
cat("FINAL RESULTS SUMMARY\n")
cat("========================================\n")
print(results_all)

# Save results
write.csv(results_all, "heatmap_results_all_methods.csv", row.names = FALSE)
cat("\nResults saved to 'heatmap_results_all_methods.csv'\n")


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
  # Generate data from logistic distribution
  x = rnorm(n*reps, 0, 1)
  x_mat = matrix(x, nrow=reps, ncol=n)
  
  if(type == "cramer")
    return(apply(X=x_mat, MARGIN=1, FUN=Cramer, cdf = pnorm, ep=ep))
  if(type == "ks")
    return(apply(X=x_mat, MARGIN=1, FUN=KS, cdf = pnorm, ep=ep))
}

# ========== Simulation ==========

reps = 1000
al = 0.05
nVec = c(100, 200, 400, 800, 1600, 3200)
epVec = c(0.5, 1, 2)
mu1_values = c(0, 0.03, 0.05, 0.1, 0.2, 0.4)  # Location shifts

# Initialize results storage
results_KS = data.frame()

for (ep_idx in 1:length(epVec)) {
  ep = epVec[ep_idx]
  print(paste("Running epsilon =", ep))
  
  for (mu1_idx in 1:length(mu1_values)) {
    mu1 = mu1_values[mu1_idx]
    print(paste("  Running mu1 =", mu1))
    
    p_GOF_KS = p_Cramer = rep(0, reps)
    power_GOF_KS = power_Cramer = rep(0, length(nVec))
    
    for (k in 1:length(nVec)){
      print(paste("    Progress: n =", nVec[k], "(", k, "/", length(nVec), ")"))
      n = nVec[k]
      
      # Generate null distributions (no Kuiper parameter)
      null_cramer = ecdf(reference(n=n, ep=ep, reps=1000, type="cramer"))
      null_ks = ecdf(reference(n=n, ep=ep, reps=1000, type="ks"))
      
      for(i in 1:reps){
        # Generate data from logistic with location shift mu1
        x = rnorm(n, mu1, 1) 
        
        # Calculate test statistics
        cramer = Cramer(x, cdf, ep) 
        ks = KS(x, cdf, ep)
        
        # Calculate p-values
        p_Cramer[i] = 1 - null_cramer(cramer)
        p_GOF_KS[i] = 1 - null_ks(ks)
      }
      
      # Calculate power
      power_GOF_KS[k] = mean(p_GOF_KS < al)
      power_Cramer[k] = mean(p_Cramer < al)
      
      # Store results for this combination
      results_KS = rbind(results_KS, data.frame(
        epsilon = ep,
        n_samples = n,
        mu1 = mu1,
        power_KS = round(power_GOF_KS[k], 3),
        power_Cramer = round(power_Cramer[k], 3)
      ))
    }
  }
}

# Print results
print(results_KS)

# Save results
write.csv(results_KS, "ks_test_results_one_sided_with_mu1.csv", row.names = FALSE)

# Merge with results_all (contains LLR and ddcLLR)
if (exists("results_all")) {
  results_combined = merge(results_all, results_KS, by = c("epsilon", "n_samples", "mu1"))
  
  cat("\n\n========================================\n")
  cat("COMBINED RESULTS (LLR, ddcLLR, KS, Cramér)\n")
  cat("========================================\n")
  print(results_combined)
  
  write.csv(results_combined, "combined_all_methods_heatmap_results.csv", row.names = FALSE)
} else {
  cat("\nWarning: results_all not found. Using only KS and Cramer results.\n")
  results_combined = results_KS
}

# ========================================
# VISUALIZATION
# ========================================

library(ggplot2)
library(reshape2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(cowplot)


if (exists("results_all")) {
  for (eps in epVec) {
    plots_list = list()
    plot_counter = 0
    
    for (test_type in c("power_LLR", "power_ddcLLR", "power_KS", "power_Cramer")) {
      plot_counter = plot_counter + 1
      
      # Use results_combined
      subset_data = results_combined[results_combined$epsilon == eps, c("mu1", "n_samples", test_type)]
      
      # Convert to data.frame
      subset_data = as.data.frame(subset_data)
      
      # Reshape data for heatmap - use reshape2 explicitly
      heatmap_data = reshape2::dcast(subset_data, mu1 ~ n_samples, value.var = test_type)
      rownames(heatmap_data) = heatmap_data$mu1
      heatmap_data$mu1 = NULL
      heatmap_matrix = as.matrix(heatmap_data)
      
      # Melt for ggplot - use reshape2 explicitly
      melted_data = reshape2::melt(heatmap_matrix)
      colnames(melted_data) = c("mu1", "n_samples", "power")
      
      # Create formatted labels
      melted_data$label = ifelse(melted_data$power == 1, 
                                 "1",
                                 sub("^0", "", sprintf("%.3f", melted_data$power)))
      
      # Create heatmap
      test_label = switch(test_type,
                          "power_LLR" = "LLR",
                          "power_ddcLLR" = "ddcLLR",
                          "power_KS" = "KS",
                          "power_Cramer" = "CvM")
      
      # Show y-axis label only for left column (plots 1 and 3)
      y_label = if(plot_counter %in% c(1, 3)) expression(theta[1]) else ""
      
      # Show x-axis label only for bottom row (plots 3 and 4)
      x_label = if(plot_counter %in% c(3, 4)) "Sample Size" else ""
      
      p = ggplot(melted_data, aes(x = factor(n_samples), y = factor(mu1), fill = power)) +
        geom_tile(color = "lightgray", size = 0.5) +
        geom_text(aes(label = label,
                      color = power > 0.4),
                  size = 5) +
        scale_color_manual(values = c("black", "white"), guide = "none") +
        scale_fill_gradient(low = "lightgray", high = "black",
                            limits = c(0, 1), name = "Power") +
        labs(title = test_label,
             x = x_label, y = y_label) +
        theme_minimal() +
        theme(
          axis.text.x = if(plot_counter %in% c(1, 2)) element_blank() else element_text(angle = 41, hjust = 1, size = 13, face = "bold"),
          axis.text.y = if(plot_counter %in% c(2, 4)) element_blank() else element_text(size = 13, face = "bold"),
          axis.title.x = if(plot_counter %in% c(1, 2)) element_blank() else element_text(size = 17, face = "bold"),
          axis.title.y = element_text(size = 17, face = "bold"),
          axis.ticks.x = if(plot_counter %in% c(1, 2)) element_blank() else element_line(),
          axis.ticks.y = if(plot_counter %in% c(2, 4)) element_blank() else element_line(),
          legend.position = "none",
          plot.title = element_text(size = 17, hjust = 0.5, face = "bold"),
          plot.margin = unit(c(0, 0, 0, 0), "pt"),
          panel.spacing = unit(0, "pt"),
          plot.background = element_blank()
        )
      
      plots_list[[length(plots_list) + 1]] = p
    }
    
    # Define row proportions (adjust these values as needed)
    # For example, 0.425 for top row and 0.575 for bottom row
    top_height = 0.425
    bottom_height = 0.575
    
    # Arrange four heatmaps in 2x2 grid with custom proportions
    plot_grid_2x2 = 
      ggdraw() +
      draw_plot(plots_list[[1]], 0,   bottom_height, 0.5, top_height) +  # Top-left
      draw_plot(plots_list[[2]], 0.5, bottom_height, 0.5, top_height) +  # Top-right
      draw_plot(plots_list[[3]], 0,   0,             0.5, bottom_height) +  # Bottom-left
      draw_plot(plots_list[[4]], 0.5, 0,             0.5, bottom_height)    # Bottom-right
    
    # Add title
    title = ggdraw() + 
      draw_label(paste0(""), 
                 fontface = 'bold', size = 14)
    
    final_plot_with_title = plot_grid(title, plot_grid_2x2, ncol = 1, 
                                      rel_heights = c(0.05, 1))
    
    print(final_plot_with_title)
    
    # Save the plot
    ggsave(filename = paste0("heatmap_4_methods_paper_epsilon_", eps, ".png"), 
           plot = final_plot_with_title, width = 10, height = 8, dpi = 300)
    cat("\nFour-method heatmap (for paper) saved for epsilon =", eps, "\n")
  }
}
