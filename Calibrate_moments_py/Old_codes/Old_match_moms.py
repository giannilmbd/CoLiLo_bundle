# %%
import numpy as np
from scipy.optimize import minimize, LinearConstraint
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
# %%
def f1(param_vec, var1, skew1, kurt1):
    # Extract parameters
    p1, p2, mu2, mu3, s1, s2, s3 = param_vec

    # Calculate intermediate values
    # v1 = s1 ** 2
    # v2 = s2 ** 2
    # v3 = s3 ** 2
    p3 = 1 - p1 - p2
    mu1 = -p2 / p1 * mu2 - p3 / p1 * mu3
    m1=0.0
    m2= p1*(s1**2 + mu1**2) + p2*(s2**2 + mu2**2) + (1-p1-p2)*(s3**2 + mu3**2) - (p1*mu1 + p2*mu2 + (1-p1-p2)*mu3)**2
    m3 = (p1*(mu1**3 + 3*mu1*s1**2) + p2*(mu2**3 + 3*mu2*s2**2) + (1-p1-p2)*(mu3**3 + 3*mu3*s3**2) - 3*(p1*mu1 + p2*mu2 + (1-p1-p2)*mu3)*m2 - (p1*mu1**3 + p2*mu2**3 + (1-p1-p2)*mu3**3)) / m2**(3/2)
    m4 = (p1*(mu1**4 + 6*mu1**2*s1**2 + 3*s1**4) + p2*(mu2**4 + 6*mu2**2*s2**2 + 3*s2**4) + (1-p1-p2)*(mu3**4 + 6*mu3**2*s3**2 + 3*s3**4) - 4*(p1*mu1 + p2*mu2 + (1-p1-p2)*mu3)*(p1*mu1**3 + p2*mu2**3 + (1-p1-p2)*mu3**3) + 6*(p1*mu1 + p2*mu2 + (1-p1-p2)*mu3)**2*m2 - m2**2 - 3*(p1*mu1**4 + p2*mu2**4 + (1-p1-p2)*mu3**4)) / m2**2

    a1 = np.abs(np.exp(np.abs(m2 - var1) * 100) - 1)
    a2 = np.abs(m3 - skew1)
    a3 = np.abs(m4 - kurt1)
    out = np.linalg.norm([a1, a2, a3], 8)
    return out

# %%
# Set up vars vector
from_val = 0.009499 ** 2
to_val = 0.039910 ** 2
n_vals = 200
vars_vec = np.linspace(from_val, to_val, n_vals)

# %%
# Set up skewrng vector
from_val = -6.0349
to_val = 0.7388
n_valss = 199
step = (to_val - from_val) / (n_valss - 1)
skewrng_vec = np.concatenate(([0], np.arange(to_val, from_val - step, -step)))
# %%
# Set up kurt1 vector
from_val = 5.0
to_val = 40.0
n_valsk = 199
kurt1_vec = np.concatenate(([3.0], np.linspace(from_val, to_val, n_valsk) ))+ 3.0


 # %% 

param_grid = np.zeros((n_vals, 7))

# Define a function for parallel processing
def optimize_params(var1,skew1,kurt1):
    # Set initial parameter values
    p1, p2, mu2, mu3, s1, s2, s3 = 0.33, 0.33, 0.0, 0.0, 0.01, 0.01, 0.01
    init_params = [p1, p2, mu2, mu3, s1, s2, s3]
    constraint_matrix = [[1, 1, 0, 0, 0, 0, 0]]  # coefficient matrix for the constraint
    constraint_rhs = [1]  # right-hand side of the constraint
    constraint = LinearConstraint(constraint_matrix, -np.inf, constraint_rhs)

    # Minimize the objective function
    res = minimize(f1, init_params, args=(var1, skew1, kurt1),constraints=[constraint],
                   tol=1e-10, options={'verbose': 0,'disp': True,'maxiter': 2000},
                   bounds = [(0.0000001, .9999999), (0.0000001, .9999999), (None, None), (None, None), (0, None), (0, None), (0, None)])

    return res.x.squeeze()
# %% variance
num_jobs = -1  # Use all available CPU cores
skew1=0
kurt1=3
results = Parallel(n_jobs=num_jobs)(delayed(optimize_params)(var1,skew1,kurt1) for var1 in vars_vec)
newresults=np.reshape( results,(len(vars_vec),7))
f1_vals = [f1(newresults[i,:],vars_vec[i],skew1,kurt1) for i in range(0,len(vars_vec))]
# plot residuals
plt.plot(vars_vec,f1_vals)
plt.xlabel('Optimization run')
plt.ylabel('Residual')
plt.title('Residuals of f1')
plt.show()
# Save optimized parameter values to param_grid
for i in range(len(vars_vec)):
    param_grid[i, :] = newresults[i, :]
# Save results to CSV
np.savetxt("variance.csv", param_grid, delimiter=",")

# %% skewness

var1=vars_vec[0]
kurt1=3
num_jobs = -1  # Use all available CPU cores
results = Parallel(n_jobs=num_jobs)(delayed(optimize_params)(var1,skew1,kurt1) for skew1 in skewrng_vec)

newresults=np.reshape( results,(len(skewrng_vec ),7))
f1_vals = [f1(newresults[i,:],var1,skewrng_vec[i],kurt1) for i in range(0,len(skewrng_vec))]
# plot residuals
plt.plot(skewrng_vec,f1_vals)
plt.xlabel('moments')
plt.ylabel('Residual')
plt.title('Residuals of f1')
plt.show()
# Save optimized parameter values to param_grid
for i in range(len(skewrng_vec)):
    param_grid[i, :] = newresults[i, :]
# Save results to CSV
np.savetxt("skewness.csv", param_grid, delimiter=",")

# %% kurtosis
num_jobs = -1  # Use all available CPU cores
var1=vars_vec[0]
skew1=0
results = Parallel(n_jobs=num_jobs)(delayed(optimize_params)(var1,skew1,kurt1) for kurt1 in kurt1_vec)
newresults=np.reshape( results,(len(kurt1_vec),7))
f1_vals = [f1(newresults[i,:],var1,skew1,kurt1_vec[i]) for i in range(0,len(kurt1_vec))]

plt.plot(f1_vals)
plt.xlabel('Optimization run')
plt.ylabel('Residual')
plt.title('Residuals of f1')
plt.show()
# Save optimized parameter values to param_grid
for i in range(len(kurt1_vec)):
    param_grid[i, :] = newresults[i, :]
# Save results to CSV
np.savetxt("kurtosis.csv", param_grid, delimiter=",")


# %%
