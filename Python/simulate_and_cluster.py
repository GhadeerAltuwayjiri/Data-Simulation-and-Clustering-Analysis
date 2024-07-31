import numpy as np
from sklearn.mixture import GaussianMixture
import warnings

# Define the cov_sim variations
cov_sim_values = [0, 0.1, 0.2, 0.5, 0.7]

# Initialize a list to hold the simulated data for each cov_sim value
simulated_data_lists = {}

# Assume mu and sigma are defined somewhere in the code
# Replace the following placeholders with actual values
mu = np.random.rand(6, 9)
sigma = np.array([np.eye(9) for _ in range(6)])  # Example covariance matrices
pi = np.random.dirichlet([1]*6)
G = 6
n = 100

# Iterate over the cov_sim values
for cov_sim in cov_sim_values:
    # Initialize lists for this cov_sim
    simulated_sample_list = []
    flow_sim_list = []
    MM_list = []
    Zhat_list1 = []
    Ztoy_list1 = []
    pi1 = []
    
    # Adjust the number of simulations if necessary
    num_simulations = 5
    
    # Loop to generate data, perform GaussianMixture, and store in lists
    for j in range(num_simulations):
        np.random.seed(123 + j)
        
        # Calculate standard deviation for channels
        sd_channel = np.std(mu, axis=1)
        
        # Generate modified mu for each sample
        mu_samp_list = []
        for k in range(5):
            mu_samp = mu + np.random.normal(scale=cov_sim * sd_channel[:, None], size=mu.shape)
            mu_samp_list.append(mu_samp)
        
        # Simulate datasets
        Ztoy = np.random.choice(range(G), size=n, p=pi)
        ytoy = np.array([np.random.multivariate_normal(mu_samp_list[j][:, z], sigma[z]) for z in Ztoy])
        
        try:
            flow_sim = GaussianMixture(n_components=G).fit(ytoy)
            Zhat = flow_sim.predict(ytoy)
            simulated_sample_list.append(ytoy)
            flow_sim_list.append(flow_sim)
            pi1.append(flow_sim.weights_)
            MM_list.append(ytoy)
            Zhat_list1.append(Zhat)
            Ztoy_list1.append(Ztoy)
        except Exception as e:
            warnings.warn(f"Clustering failed for simulation {j}: {e}")
    
    # Store the lists for this cov_sim value
    simulated_data_lists[f"cov_sim_{cov_sim}"] = {
        "simulated_samples": simulated_sample_list,
        "flow_sims": flow_sim_list,
        "MM": MM_list,
        "Zhat": Zhat_list1,
        "pi": pi1,
        "Ztoy": Ztoy_list1
    }
