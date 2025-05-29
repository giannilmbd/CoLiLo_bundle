# ----------------------------------------------------------------------------------------------------------------------------------------
# %% ------------------------------------  Load Packages ------------------------------------ #
import numpy as np
import pandas as pd
from scipy.stats import norm
import warnings
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import matplotlib 
from multiprocessing import Pool
import statistics 
from match_moments import match_moments
from calculate_moments import calculate_moments
from discreteGMAR import discrete_gmar
from discreteApproximation import discrete_approximation
warnings.filterwarnings('ignore') 

matplotlib.use("TkAgg") #
matplotlib.interactive(False)

#%% --------------------------------------------SET PARAMETERS -------------------------------------------------------- #
Nm = 25 #number of points for discretization

mu = np.array(.0555) # mean full DGP
A1 = np.array(0.5854) # Autorecursive parameter
# create a row vector of dim 2 for A2
A2 = np.array([0.8959, -0.3990]) # The original alorithm should work also for AR(2) but not tested

# Parameters of the Gaussian Mixture Model
# weights 
p1=.25
p2=.3
p3=1-p1-p2
pC = np.array([p1, p2, p3])
# means
muc1=.10
muc2=0.0
muc3=0.0
muC = np.array([muc1, muc2, muc3])
# Standard deviations
sigmaC1=.002;
sigmaC2=.004;
sigmaC3=.007;
sigmaC = np.array([sigmaC1, sigmaC2, sigmaC3])

# %% ----------------------------------------------------------------------------------------------------------------------------------------
T1, T2, T3, T4, TBar,var_, Sk, Kur =calculate_moments(pC, muC, sigmaC)

print(pd.DataFrame(data=[T1,var_,Sk,Kur],index=['mean(T1)','var','Sk','Kur']))
# How many moments to be used
n_moments = 4

# %% ----------------------------------------------------------------------------------------------------------------------------------------

p_even1, x_even1,match_order = discrete_gmar(mu, A1, pC, muC, sigmaC, Nm, n_moments, method='even')#'even''gauss-hermite'

PEven1_transpose = p_even1.T  # Transpose PEven1 matrix if necessary

eigenvalues, eigenvectors = np.linalg.eig(PEven1_transpose)
idx = np.argmax(np.real(eigenvalues))
pi = np.real(eigenvectors[:, idx]) / np.sum(np.real(eigenvectors[:, idx]))


XEven1_reshaped = np.reshape(x_even1, (1, -1))  # Reshape XEven1 to match pi dimensions

plt.figure()
plt.plot(XEven1_reshaped.flatten(), pi.flatten())
plt.xlabel('x')
plt.ylabel('Stationary distribution')
plt.show()
# %%
resultP=pd.DataFrame(data=p_even1)
resultP.to_csv('PiMx.csv',index=False,header=False)
resultX=pd.DataFrame(data=x_even1)
resultX.to_csv('XMx.csv',index=False,header=False)


# %%
