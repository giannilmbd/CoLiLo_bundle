#%%
import numpy as np
from scipy.stats import norm, gamma

# Define the parameters of the AR(1) process
T = 10000  # Length of time series
phi = 0.8  # AR(1) coefficient
sigma = 1  # Innovation standard deviation
mu = 0  # Innovation mean
omega = 0.5  # Skewness parameter
gamma1 = 0.2  # Coefficient of skewness
gamma2 = 0.1  # Coefficient of kurtosis

# Define the parameters for discretizing the state space
n = 100  # Number of grid points
m = 3  # Number of standard deviations to include in grid
mu_grid = np.linspace(-m*np.sqrt(sigma**2 / (1-phi**2)), m*np.sqrt(sigma**2 / (1-phi**2)), n)  # Grid for state variable

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
cond_var = sigma**2 / (1 - phi**2) * np.ones(n)
#%%
# Compute the conditional skewness and kurtosis for each grid point
cond_skew = gamma1 * np.sqrt(cond_var**3) / (cond_var - omega**2)**(3/2)
cond_skew = np.ones(n) * cond_skew # make it an array of length n

cond_kurt = gamma2 + 3*(gamma1**2)*(cond_var**2) / ((cond_var - omega**2)**2)
#%%
# Compute the discretized distribution of the AR(1) process
innov_dist = np.zeros((n, n))
for i in range(n):
    # Compute the conditional probabilities for each grid point
    cond_prob = gamma.pdf(mu_grid, cond_skew[i], scale=np.sqrt(cond_var[i]/cond_skew[i]), loc=cond_mean[i])

    # Normalize the conditional probabilities
    cond_prob /= np.sum(cond_prob)

    # Compute the discretized distribution of the innovation at this grid point
    innov_dist[i, :] = gamma.pdf(mu_grid, cond_skew[i]**2/cond_var[i], scale=np.sqrt(cond_var[i]/cond_skew[i]))
#%%
# Compute the discretized distribution of the AR(1) process
ar1_dist = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        ar1_dist[i, j] = np.sum(P[i, :]*innov_dist[j, :]*cond_prob[i])

# Normalize the discretized distribution of the AR(1) process
row_sums = np.sum(ar1_dist, axis=1, keepdims=True)
ar1_dist /= np.where(row_sums == 0, 1, row_sums) + 1e-10


# %%
import matplotlib.pyplot as plt

# Generate a sample from the true AR(1) process
ar1_sample = np.zeros(T)
ar1_sample[0] = norm.rvs(loc=mu, scale=sigma)
for t in range(1, T):
    ar1_sample[t] = phi*ar1_sample[t-1] + norm.rvs(loc=mu, scale=sigma)

# Compute the true density of the AR(1) process
true_density, true_bins = np.histogram(ar1_sample, bins=mu_grid)

# Plot the true density and the discretized density side by side
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

ax[0].bar(true_bins[:-1], true_density/T, width=true_bins[1]-true_bins[0])
ax[0].set_title('True Density')

ax[1].bar(mu_grid, ar1_dist[-1, :], width=mu_grid[1]-mu_grid[0])

ax[1].set_title('Discretized Density')
plt.show()
# %%
