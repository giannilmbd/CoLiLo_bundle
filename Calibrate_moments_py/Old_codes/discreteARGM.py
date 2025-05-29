#%%
import numpy as np
from numpy.linalg import eigh
from scipy.stats import norm, multivariate_normal
from scipy.optimize import minimize
from itertools import product
import warnings
from sklearn.mixture import GaussianMixture

#%%

def gmm_pdf(x, gmm):
    return np.sum([w*multivariate_normal(mean=mu, cov=sigma).pdf(x)
                   for w, mu, sigma in zip(gmm.weights_, gmm.means_, gmm.covariances_)])

#%%
def discrete_gmar(mu, A, pC, muC, sigmaC, Nm, nMoments=2, method='even', nSigmas=None):
    # some error checking
    if any(p < 0 for p in pC):
        raise ValueError('mixture proportions must be positive')
    if any(s < 0 for s in sigmaC):
        raise ValueError('standard deviations must be positive')
    if sum(pC) != 1:
        raise ValueError('mixture proportions must add up to 1')

    pC, muC, sigmaC = map(np.atleast_1d, (pC, muC, sigmaC))
    A = np.atleast_1d(A)
    K = len(A)
    
    F = np.vstack([A, np.eye(K-1, K)])
    rho = abs(eigh(F, eigvals_only=True)[0])

    if rho >= 1:
        raise ValueError('spectral radius must be less than one')

    # compute conditional moments
    sigmaC2 = sigmaC**2
    T1 = pC @ muC  # mean
    T2 = pC @ (muC**2 + sigmaC2)  # uncentered second moment
    T3 = pC @ (muC**3 + 3*muC*sigmaC2)  # uncentered third moment
    T4 = pC @ (muC**4 + 6*(muC**2)*sigmaC2 + 3*sigmaC2**2)  # uncentered fourth moment

    TBar = np.array([T1, T2, T3, T4])

    gmObj = multivariate_normal(muC, sigmaC2, allow_singular=True)

    # Check that Nm is a valid number of grid points
    if not (isinstance(Nm, int) and Nm > 3):
        raise ValueError('Nm must be a positive integer greater than 3')

    # Check that nMoments is a valid number
    if not (isinstance(nMoments, int) and 1 <= nMoments <= 4):
        raise ValueError('nMoments must be either 1, 2, 3, 4')

    # set default nSigmas if not supplied
    if nSigmas is None:
        if rho <= 1-2/(Nm-1):
            nSigmas = np.sqrt(2*(Nm-1))
        else:
            nSigmas = np.sqrt(Nm-1)

    sigma = np.sqrt(T2 - T1**2)  # conditional standard deviation
    sigmaX = sigma * np.sqrt((np.eye(K**2) - np.kron(F, F)).sum())  # unconditional standard deviation

    # construct the one dimensional grid
    # (NOTE: Some quadrature methods are omitted as Python does not have built-in functions for them.)
    if method == 'even':
        X1 = np.linspace(mu - nSigmas*sigmaX, mu + nSigmas*sigmaX, Nm)
    else:
        raise ValueError(f"Unknown method: {method}")

    X = np.array(list(product(X1, repeat=K))).T

    P = np.full((Nm**K, Nm**K), np.nan)  # transition probability matrix
    P1 = np.full((Nm**K, Nm), np.nan)  # matrix to store transition probability
    P2 = np.kron(np.eye(Nm**(K-1)), np.ones(Nm))  # Nm^K * Nm^(K-1) matrix used to construct P
    scalingFactor = max(abs(X1))
    kappa = 1e-8

    for ii in range(Nm**K):
        condMean = mu*(1 - sum(A)) + A @ X[:, ii]
        xPDF = (X1 - condMean).T

    switcher = {
        'gauss-hermite': W * (gmObj_pdf(xPDF) / norm.pdf(xPDF, 0, sigma)),
        'GMQ': W * (gmObj_pdf(xPDF) / gmObj_pdf(X1.reshape(-1, 1))),
        'otherwise': W * gmObj_pdf(xPDF),
    }

    q = switcher.get(method, "Invalid method")

    if any(q < kappa):
        q[q < kappa] = kappa

        if nMoments == 1:  # match only 1 moment
            P1[ii, :] = discreteApproximation(X1, lambda x: (x - condMean) / scalingFactor, TBar[0] / scalingFactor, q, 0)
        else:  # match 2 moments first
            p, lambda_, momentError = discreteApproximation(X1, lambda x: [(x - condMean) / scalingFactor, ((x - condMean) / scalingFactor) ** 2], TBar[:2] / scalingFactor ** np.arange(1, 3), q, np.zeros(2))

            if np.linalg.norm(momentError) > 1e-5:  # if 2 moments fail, then just match 1 moment
                warnings.warn('Failed to match first 2 moments. Just matching 1.')
                P1[ii, :] = discreteApproximation(X1, lambda x: (x - condMean) / scalingFactor, TBar[0] / scalingFactor, q, 0)
            elif nMoments == 2:
                P1[ii, :] = p
            elif nMoments == 3:
                p_new, _, momentError = discreteApproximation(X1, lambda x: [(x - condMean) / scalingFactor, ((x - condMean) / scalingFactor) ** 2, ((x - condMean) / scalingFactor) ** 3], TBar[:3] / scalingFactor ** np.arange(1, 4), q, np.concatenate([lambda_, [0]]))

                if np.linalg.norm(momentError) > 1e-5:
                    warnings.warn('Failed to match first 3 moments.  Just matching 2.')
                    P1[ii, :] = p
                else:
                    P1[ii, :] = p_new
            else:  # 4 moments
                p_new, _, momentError = discreteApproximation(X1, lambda x: [(x - condMean) / scalingFactor, ((x - condMean) / scalingFactor) ** 2, ((x - condMean) / scalingFactor) ** 3, ((x - condMean) / scalingFactor) ** 4], TBar / scalingFactor ** np.arange(1, 5), q, np.concatenate([lambda_, [0, 0]]))

                if np.linalg.norm(momentError) > 1e-5:
                    p_new, _, momentError = discreteApproximation(X1, lambda x: [(x - condMean) / scalingFactor, ((x - condMean) / scalingFactor) ** 2, ((x - condMean) / scalingFactor) ** 3], TBar[:3] / scalingFactor ** np.arange(1, 4), q, np.concatenate([lambda_, [0]]))

                    if np.linalg.norm(momentError) > 1e-5:
                        warnings.warn('Failed to match first 3 moments.  Just matching 2.')
                        P1[ii, :] = p
                    else:
                        P1[ii, :] = p_new
                        warnings.warn('Failed to match first 4 moments.  Just matching 3.')
                else:
                    P1[ii, :] = p_new

            P[ii, :] = np.kron(P1[ii, :], P2[ii, :])






# %%
