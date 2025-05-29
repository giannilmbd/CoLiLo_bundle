# %%
import numpy as np
from scipy.stats import norm, gamma
#%%
# Define the parameters of the AR(1) process
phi = 0.8  # AR(1) coefficient
sigma = 1  # Innovation standard deviation
mu = 0  # Innovation mean
omega = 0.5  # Skewness parameter
gamma1 = 0.2  # Coefficient of skewness
gamma2 = 0.1  # Coefficient of kurtosis

# Define the parameters for discretizing the state space
n = 10  # Number of grid points
m = 3  # Number of standard deviations to include in grid
mu_grid = np.linspace(-m*np.sqrt(sigma**2 / (1-phi**2)), m*np.sqrt(sigma**2 / (1-phi**2)), n)  # Grid for state variable
#%%
# Compute the transition probabilities for each grid point
P = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if j == 0:
            cdf = norm.cdf((mu_grid[0] - phi*mu_grid[i] + mu) / sigma, loc=mu, scale=sigma)
            P[i, j] = cdf
        elif j == n-1:
            cdf = norm.cdf((mu_grid[n-1] - phi*mu_grid[i] + mu) / sigma, loc=mu, scale=sigma)
            P[i, j] = 1 - cdf
        else:
            cdf1 = norm.cdf((mu_grid[j] - phi*mu_grid[i] + mu) / sigma, loc=mu, scale=sigma)
            cdf2 = norm.cdf((mu_grid[j+1] - phi*mu_grid[i] + mu) / sigma, loc=mu, scale=sigma)
            P[i, j] = cdf2 - cdf1
#%%
# Compute the conditional mean and variance for each grid point
cond_mean = phi*mu_grid
cond_var = sigma**2 + omega**2
#%%
# Compute the conditional skewness and kurtosis for each grid point
cond_skew = gamma1 * np.sqrt(cond_var**3) / (cond_var - omega**2)**(3/2)
cond_kurt = gamma2 + 3*(gamma1**2)*(cond_var**2) / ((cond_var - omega**2)**2)
#%%
# Compute the conditional probabilities for each grid point
cond_prob = gamma.pdf(mu_grid, cond_skew, scale=np.sqrt(cond_var/cond_skew), loc=cond_mean)
#%%
# Normalize the conditional probabilities
cond_prob /= np.sum(cond_prob)
#%%
# Compute the discretized distribution of the innovation
innov_dist = [gamma(a=cond_var/cond_skew**2, scale=np.sqrt(cond_var/cond_skew)).pdf(mu_grid)]
#%%
# Compute the discretized distribution of the AR(1) process
for i in range(1, n):
    innov_dist.append(np.dot(P, innov_dist[i-1]) * cond_prob)
#%%
# Normalize the discretized distribution of the AR(1) process
innov_dist = np.array(innov_dist)
ar1_dist = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        ar1_dist[i, j] = np.sum(innov_dist[:i+1, j])
    ar1_dist[i, :] /= np.sum(ar1_dist[i, :])


# %%
