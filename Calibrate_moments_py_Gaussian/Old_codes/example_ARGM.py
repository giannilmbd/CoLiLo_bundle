#%%
import numpy as np
import pandas as pd
from scipy.stats import norm
from warnings import warn
from discreteApproximation import discreteApproximation
from discreteARGM import discreteARGM
from sklearn.mixture import GaussianMixture
import importlib
#%%
# Load data
data = pd.read_csv('dGrowth.csv',header=None)
dGrowth = data.values
T = len(dGrowth)
dX = np.column_stack([np.ones(T-1), dGrowth[:-1]])
dY = dGrowth[1:]
#%%
# Estimate OLS parameters
betaHat = np.linalg.lstsq(dX, dY, rcond=None)[0]
rhoHat = betaHat[1]
muHat = betaHat[0] / (1 - betaHat[1])
regResid = dY - dX @ betaHat  # residuals
#%%
# Fit Gaussian mixture
nComp = 3  # number of components
gmObj = GaussianMixture(n_components=nComp, n_init=10).fit(regResid.reshape(-1, 1))
#%%
N = 15  # number of grid points
nMoments = 4  # number of moments to match
#%%
# Call discreteARGM function
P4, D4 = discreteARGM(muHat, rhoHat, gmObj, N, nMoments, 'even')
# 'method' can be 'even', 'gauss-hermite', 'gauss-legendre', or 'clenshaw-curtis'.

# %%
