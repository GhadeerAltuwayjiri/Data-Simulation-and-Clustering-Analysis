library(MASS)
library(flowClust)

# Define the cov_sim variations
cov_sim_values <- c(0, 0.1, 0.2, 0.5, 0.7)

# Initialize a list to hold the simulated data for each cov_sim value
simulated_data_lists <- list()

# Iterate over the cov_sim values
for (cov_index in seq_along(cov_sim_values)) {
  cov_sim <- cov_sim_values[cov_index]
  
  # Initialize lists for this cov_sim
  simulated_sample_list <- list()
  flow_sim_list <- list()
  MM_list <- list()
  Zhat_list1 <- list()
  Ztoy_list1 <- list()
  pi1 <- list()
  
  # Adjust the number of simulations if necessary
  num_simulations <- 5
  
  # Loop to generate data, perform flowClust, and store in lists
  for (j in 1:num_simulations) {
    set.seed(123 + j)
    
    # Calculate standard deviation for channels
    sd_channel <- apply(mu, 1, sd)
    
    # Generate modified mu for each sample
    mu_samp_list <- list()
    for (k in 1:5) {
      mu_samp_list[[k]] <- mu + matrix(rnorm(6 * 9, sd = cov_sim * sd_channel), nrow = 6, ncol = 9)
    }
    
    # Simulate datasets
    Ztoy <- sample(1:G, size = n, prob = pi, replace = TRUE)
    ytoy <- t(sapply(Ztoy, function(z) MASS::mvrnorm(n = 1, mu = mu_samp_list[[j]][, z], Sigma = sigma[, , z])))
    colnames(ytoy) <- Memory.mgd7.PS80.15.2[[9]]@varNames
    
    simulated_sample <- new("flowFrame", exprs = as.matrix(ytoy), parameters = parameters(Memory.mgd7.PS80.15.S1$lymphocyte))
    flow_sim <- flowClust(simulated_sample, K = 9, trans = 0, nu = Inf)
    
    if (!is.null(flow_sim@z)) {
      simulated_sample_list[[j]] <- simulated_sample
      flow_sim_list[[j]] <- flow_sim
      pi1[[j]] <- flow_sim@w
      # Store the simulated data and clustering results
      MM_list[[j]] <- ytoy
      Zhat <- apply(flow_sim@z, 1, which.max)
      Zhat_list1[[j]] <- Zhat
      Ztoy_list1[[j]] <- Ztoy
    } else {
      warning(paste("Clustering failed for simulation", j))
    }
  }
  
  # Store the lists for this cov_sim value
  simulated_data_lists[[paste("cov_sim", cov_sim, sep = "_")]] <- list(
    simulated_samples = simulated_sample_list, 
    flow_sims = flow_sim_list, 
    MM = MM_list, 
    Zhat = Zhat_list1,
    pi = pi1,
    Ztoy = Ztoy_list1
  )
}
