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

# ========== MLR GDP TEST ==========
# For distributions with Monotone Likelihood Ratio property,
# the sufficient statistic is X itself

MLR.GDP.Test <- function(X, dp_info) {
  n <- length(X)
  
  # For MLR, we directly use X as the test statistic
  test_statistic <- n * dp.meanEst(X, dp_info)
  
  return(list(
    test_statistic = test_statistic,
    data = X
  ))
}

# ========== Standard MLR Test ==========
# Standard test using sample mean of X

MLR.Test <- function(X) {
  test_statistic <- sum(X)
  
  return(list(
    test_statistic = test_statistic,
    data = X
  ))
}

# ========== Simulation ==========

# Define parameters for heatmap analysis
sample_sizes <- c(100, 200, 400, 800, 1600, 3200)
mu1_values <- c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2)
epsilon_values <- c(0.5, 1, 2)

# Initialize storage for results
results_heatmap <- data.frame(
  epsilon = numeric(),
  n_samples = integer(),
  mu1 = numeric(),
  power_MLR = numeric(),
  power_ddcMLR = numeric()
)

# Number of simulations
n_simulations_H0_samples <- 1000  # For null distribution
n_simulations_H1 <- 1000           # For power calculation

# Loop over epsilon values (each creates one heatmap)
for (epsilon in epsilon_values) {
  cat("\n########################################\n")
  cat("####### EPSILON =", epsilon, "#######\n")
  cat("########################################\n")
  
  # Loop over sample sizes (columns of heatmap)
  for (n_samples in sample_sizes) {
    
    # Setup for GDP test
    dp_info <- list(
      q = c("lower", "upper"),
      epsilon = epsilon,
      alpha = 3,
      l_m = -Inf,
      u_m = Inf,
      u_s = 10
    )
    
    # Step 1: Generate null distribution samples (H0: data from P, location = 0)
    MLR_H0_samples <- numeric(n_simulations_H0_samples)
    ddcMLR_H0_samples <- numeric(n_simulations_H0_samples)
    
    cat("Generating null distribution for n =", n_samples, "...\n")
    pb <- txtProgressBar(min = 0, max = n_simulations_H0_samples, style = 3)
    set.seed(12345 + 1000 * epsilon + n_samples)
    
    for (i in 1:n_simulations_H0_samples) {
      X_from_P <- rlogis(n_samples, 0, 1)
      MLR_H0_samples[i] <- MLR.Test(X_from_P)$test_statistic
      ddcMLR_H0_samples[i] <- MLR.GDP.Test(X_from_P, dp_info)$test_statistic
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Compute the expected value under H0 (center point)
    E0_T_MLR <- mean(MLR_H0_samples)
    E0_T_ddcMLR <- mean(ddcMLR_H0_samples)
    
    # Two-sided test statistic under H0: |T - E0[T]|
    MLR_H0_abs <- abs(MLR_H0_samples - E0_T_MLR)
    ddcMLR_H0_abs <- abs(ddcMLR_H0_samples - E0_T_ddcMLR)
    
    # Loop over mu1 values (rows of heatmap) for power calculation
    for (mu1 in mu1_values) {
      cat("Processing epsilon =", epsilon, ", mu1 =", mu1, ", n =", n_samples, "\n")
      
      # Step 2: Calculate Power under H1 (data from Q with location mu1)
      MLR_H1_samples <- numeric(n_simulations_H1)
      ddcMLR_H1_samples <- numeric(n_simulations_H1)
      
      cat("Calculating power...\n")
      pb <- txtProgressBar(min = 0, max = n_simulations_H1, style = 3)
      set.seed(54321 + 1000 * epsilon + 100 * abs(mu1) * 100 + n_samples)
      
      for (i in 1:n_simulations_H1) {
        # Generate data from the TRUE alternative (actual mu1)
        X_from_Q <- rlogis(n_samples, mu1, 1)
        
        # Calculate test statistics directly from X
        MLR_H1_samples[i] <- MLR.Test(X_from_Q)$test_statistic
        ddcMLR_H1_samples[i] <- MLR.GDP.Test(X_from_Q, dp_info)$test_statistic
        
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # Two-sided test statistic under H1: |T - E0[T]|
      # Use the SAME centering point E0_T from the null distribution
      MLR_H1_abs <- abs(MLR_H1_samples - E0_T_MLR)
      ddcMLR_H1_abs <- abs(ddcMLR_H1_samples - E0_T_ddcMLR)
      
      # Calculate power for two-sided test
      # Reject H0 if |T - E0[T]| is large (exceeds critical value from H0)
      power_MLR <- mean(sapply(MLR_H1_abs, function(t_obs) {
        pvalue <- (1 + sum(MLR_H0_abs >= t_obs)) / (n_simulations_H0_samples + 1)
        pvalue <= 0.05
      }))
      
      power_ddcMLR <- mean(sapply(ddcMLR_H1_abs, function(t_obs) {
        pvalue <- (1 + sum(ddcMLR_H0_abs >= t_obs)) / (n_simulations_H0_samples + 1)
        pvalue <= 0.05
      }))
      
      # Store results
      results_heatmap <- rbind(results_heatmap, data.frame(
        epsilon = epsilon,
        n_samples = n_samples,
        mu1 = mu1,
        power_MLR = round(power_MLR, 3),
        power_ddcMLR = round(power_ddcMLR, 3)
      ))
      
      cat("MLR Power:", round(power_MLR, 3), "\n")
      cat("ddcMLR Power:", round(power_ddcMLR, 3), "\n")
    }
  }
}

# Print final summary
cat("\n\n========================================\n")
cat("FINAL RESULTS SUMMARY\n")
cat("========================================\n")
print(results_heatmap)

# Save results
write.csv(results_heatmap, "heatmap_twosided_all_methods.csv", row.names = FALSE)
cat("\nResults saved to 'heatmap_twosided_all_methods.csv'\n")

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


KS <- function(x,cdf,ep,Kuiper = FALSE){
  n = length(x)
  ks = ks.test(x,cdf,alternative = "two.sided")$statistic
  if(Kuiper == "Kuiper"){ ### No Kuiper for one-sided ###
    Dplus = ks.test(x,cdf,alternative="greater")$statistic
    Dminus = ks.test(x,cdf,alternative="less")$statistic
    ks = Dplus + Dminus
  }
  #N = (1/n)*rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
  N = (1/n)*rnorm(n=1, mean = 0, sd = 1/ep)
  ks = ks+N
  return(ks)
}

Unknown_GOF<-function(x,ep,Kuiper = FALSE){
  #the minimum KS parameter estimates, estimate mean and sd of normal based on the sample x
  KSEstimator <- function(x){
    KSDis <- function(param){
      ks_min = ks.test(x,"pnorm",mean =param[1],sd = param[2],alternative="two.sided")$statistic
      if(Kuiper == "Kuiper"){
        Dplus = ks.test(x,"pnorm",mean =param[1],sd = param[2],alternative="greater")$statistic
        Dminus = ks.test(x,"pnorm",mean =param[1],sd = param[2],alternative="less")$statistic
        ks_min = Dplus+Dminus
      }
      return(ks_min)
    }
    res <- optim(c(0,1), fn = KSDis)$par
    return(list(mean = res[1], sd = res[2]))
  }
  ks_param = KSEstimator(x=x)
  fitted_cdf = function(t){
    return(pnorm(t,mean = ks_param$mean, sd = ks_param$sd))
  }
  #Then to evaluate the distance we can just use the above implementation of KS:
  ks = KS(x,fitted_cdf,ep = ep,Kuiper = Kuiper)
  return(ks)
}

Cramer <- function(x,cdf,ep) {
  U = sort(x)
  n = length(x)
  rank = seq_len(n)
  ###JA "U" should be replaced with cdf(U)
  omega2 = 1/(12 * n) + sum(((2*rank - 1)/(2*n)-cdf(U))^2)
  omega2 = sqrt(omega2/n)
  ### Add (1/n)*rlaplace(s=1/ep)
  #omega2 = omega2 + (1/n)*rlaplace(s=1/ep)
  ### Add (1/n)*rnorm(1/ep)
  omega2 = omega2 + (1/n)*rnorm(1,0,1/ep)
  return(omega2)
}

cdf <- function(q){
  return (pnorm(q,m=0,s=1))
}

pdf<- function(x){
  return(dnorm(x,0,1))
}

reference = function(n,ep,reps,type,Kuiper = FALSE){
  #x = rnorm(n*reps,0,1)
  #x = rlaplace(n*reps,0,1)
  x = rlogis(n*reps,0,1)
  x_mat = matrix(x,nrow=reps,ncol=n)
  if(type == "cramer")
    return(apply(X=x_mat,MARGIN=1, FUN=Cramer,cdf = pnorm,ep=ep))
  if(type == "ks")
    return(apply(X=x_mat,MARGIN=1, FUN=KS,cdf = pnorm,ep=ep,Kuiper = Kuiper))
  if(type == "Unknown")
    return(apply(X=x_mat,MARGIN=1, FUN=Unknown_GOF, ep=ep,Kuiper = Kuiper))
}


# Define parameters
reps = 1000
al = 0.05
nVec = c(100, 200, 400, 800, 1600, 3200)
epVec = c(0.5, 1, 2)
mu1_values = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2)  # Location shifts

# Initialize results storage
results_KS = data.frame()

for (ep_idx in 1:length(epVec)) {
  ep = epVec[ep_idx]
  print(paste("Running epsilon =", ep))
  
  for (mu1_idx in 1:length(mu1_values)) {
    mu1 = mu1_values[mu1_idx]
    print(paste("  Running mu1 =", mu1))
    
    p_GOF_KS = p_Cramer = p_KS_V = rep(0, reps)
    power_GOF_KS = power_Cramer = power_ksv = rep(0, length(nVec))
    
    for (k in 1:length(nVec)){
      print(paste("    Progress: n =", nVec[k], "(", k, "/", length(nVec), ")"))
      n = nVec[k]
      
      # Generate null distributions
      null_cramer = ecdf(reference(n=n, ep=ep, reps=1000, type="cramer"))
      null_ks = ecdf(reference(n=n, ep=ep, reps=1000, type="ks"))
      null_ksv = ecdf(reference(n=n, ep=ep, reps=1000, type="ks", Kuiper="Kuiper"))
      
      for(i in 1:reps){
        # Generate data from logistic with location shift mu1
        x = rlogis(n, mu1, 1) 
        
        # Calculate test statistics
        cramer = Cramer(x, cdf, ep) 
        ks = KS(x, cdf, ep, Kuiper=FALSE)
        ksv = KS(x, cdf, ep, Kuiper="Kuiper")
        
        # Calculate p-values
        p_Cramer[i] = 1 - null_cramer(cramer)
        p_GOF_KS[i] = 1 - null_ks(ks)
        p_KS_V[i] = 1 - null_ksv(ksv)
      }
      
      # Calculate power
      power_GOF_KS[k] = mean(p_GOF_KS < al)
      power_Cramer[k] = mean(p_Cramer < al)
      power_ksv[k] = mean(p_KS_V < al)
      
      # Store results for this combination
      results_KS = rbind(results_KS, data.frame(
        epsilon = ep,
        n_samples = n,
        mu1 = mu1,
        power_KS = round(power_GOF_KS[k], 3),
        power_Cramer = round(power_Cramer[k], 3),
        power_KS_Kuiper = round(power_ksv[k], 3)
      ))
    }
  }
}

# Print results
print(results_KS)

# Save results
write.csv(results_KS, "ks_test_results_with_mu1.csv", row.names = FALSE)

# ========================================
# MERGE WITH results_heatmap (MLR and ddcMLR)
# ========================================
if (exists("results_heatmap")) {
  results_combined = merge(results_heatmap, results_KS, 
                           by = c("epsilon", "n_samples", "mu1"),
                           all = TRUE)
  
  cat("\n\n========================================\n")
  cat("COMBINED RESULTS (MLR, ddcMLR, KS, Cramér, Kuiper)\n")
  cat("========================================\n")
  print(results_combined)
  
  write.csv(results_combined, "combined_twosided_all_methods_results.csv", row.names = FALSE)
} else {
  cat("\nWarning: results_heatmap not found. Using only KS results.\n")
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
    
    for (test_type in c("power_MLR", "power_ddcMLR", "power_KS", "power_Cramer")) {
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
                          "power_MLR" = "Non-Private",
                          "power_ddcMLR" = "Ours",
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
          axis.text.x = if(plot_counter %in% c(1, 2)) element_blank() else element_text(angle = 33, hjust = 1, size = 15, face = "bold"),
          axis.text.y = if(plot_counter %in% c(2, 4)) element_blank() else element_text(size = 15, face = "bold"),
          axis.title.x = if(plot_counter %in% c(1, 2)) element_blank() else element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.ticks.x = if(plot_counter %in% c(1, 2)) element_blank() else element_line(),
          axis.ticks.y = if(plot_counter %in% c(2, 4)) element_blank() else element_line(),
          legend.position = "none",
          plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
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
                 fontface = 'bold', size = 15)
    
    final_plot_with_title = plot_grid(title, plot_grid_2x2, ncol = 1, 
                                      rel_heights = c(0.05, 1))
    
    print(final_plot_with_title)
    
    # Save the plot
    ggsave(filename = paste0("heatmap_4_methods_paper_epsilon_", eps, ".png"), 
           plot = final_plot_with_title, width = 10, height = 8, dpi = 300)
    cat("\nFour-method heatmap (for paper) saved for epsilon =", eps, "\n")
  }
}
