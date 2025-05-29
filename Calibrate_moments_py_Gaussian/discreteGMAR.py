#%%
import numpy as np
from numpy.linalg import eigh
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.linalg import eigvals
from itertools import product
import warnings
from discreteApproximation import discrete_approximation
from gaussian_mixture_quadrature import gaussian_mixture_quadrature
from calculate_moments import calculate_moments


def discrete_gmar(mu, A, pC, muC, sigmaC, Nm, nMoments=2, method='even', nSigmas=None, country=None):
    # some error checking
    print(f"\033[1;34m \033[43m Starting Discretization \033[0m")

    if any(p < 0 for p in pC):
        raise ValueError('mixture proportions must be positive')
    if any(s < 0 for s in sigmaC):
        raise ValueError('standard deviations must be positive')
    if sum(pC) > 1 or sum(pC) < 0.999999:
        print(pC)
        print(sum(pC))
        raise ValueError('mixture proportions must sum to 1')
        
    K = np.size(A)  # number of states
    Nm = int(Nm)  # number of discrete values for each state
    if nSigmas is None:
        nSigmas = 4
    if nMoments < 1 or nMoments > 4:
        raise ValueError('nMoments must be between 1 and 4')

    # Create identity matrix of size K-1 K and stack it under A
    I = np.eye(K)
    if K == 1:
        F = A
        rho = A
    else:
        F = np.vstack((A, np.eye(K - 1, K)))
        rho = np.abs(eigvals(F))  # Compute the spectral radius of F
        rho = rho[0]

    if np.max(np.abs(rho)) >= 1:
        raise ValueError('AR(1) coefficient must be less than 1 in absolute value')

    # Compute conditional moments
    sigmaC2 = np.power(sigmaC, 2)
    T1, T2, T3, T4, TBar, _, _, _ = calculate_moments(pC, muC, sigmaC)
    TBar = TBar.flatten()

    # Default number of moments is 2
    if nMoments is None:
        nMoments = 2

    # Check that Nm is a valid number of grid points
    if not isinstance(Nm, int) or Nm < 3:
        raise ValueError('Nm must be a positive integer greater than 3')

    # Set default nSigmas if not supplied
    if nSigmas is None:
        if rho <= 1 - 2 / (Nm - 1):
            nSigmas = np.sqrt(2 * (Nm - 1))
        else:
            nSigmas = np.sqrt(Nm - 1)

    sigma = np.sqrt(T2 - T1 ** 2)
    temp = np.eye(K ** 2) - np.kron(F, F)
    temp = np.linalg.inv(temp)
    sigmaX = sigma * np.sqrt(temp[0, 0])  # unconditional standard deviation

    # Construct the one-dimensional grid
    if method == 'even':  # evenly-spaced grid
        X1 = np.linspace(mu - nSigmas * sigmaX, mu + nSigmas * sigmaX, Nm).flatten()
        W = np.ones(Nm)
    elif method == 'gauss-legendre':  # Gauss-Legendre quadrature
        X1, W = np.polynomial.legendre.leggauss(Nm)
        X1 = X1 * (nSigmas * sigmaX) + mu
    elif method == 'clenshaw-curtis':  # Clenshaw-Curtis quadrature
        X1, W = np.polynomial.chebyshev.chebgauss(Nm)
        X1 = X1 * (nSigmas * sigmaX) + mu
        X1 = np.flipud(X1)
        W = np.flipud(W)
    elif method == 'gauss-hermite':  # Gauss-Hermite quadrature
        if rho > 0.8:
            print('Model is persistent; even-spaced grid is recommended')
        X1, W = np.polynomial.hermite.hermgauss(Nm)
        X1 = mu + np.sqrt(2) * sigma * X1
        W = W / np.sqrt(np.pi)
    elif method == 'GMQ':  # Gaussian Mixture Quadrature
        if rho > 0.8:
            print('Model is persistent; even-spaced grid is recommended')
        X1, W = gaussian_mixture_quadrature(pC, muC, sigmaC, Nm)
        X1 = X1 + mu

    X = np.array(list(product(X1, repeat=K))).T  # K * Nm^K matrix of grid points

    P = np.empty((Nm ** K, Nm ** K))  # transition probability matrix
    P1 = np.empty((Nm ** K, Nm))  # matrix to store transition probability
    P2 = np.kron(np.eye(Nm ** (K - 1)), np.ones((Nm, 1)))  # Nm^K * Nm^(K-1) matrix used to construct P
    scalingFactor = np.max(np.abs(X1))
    kappa = 1e-8

    match_order = np.empty(Nm ** K)

    for ii in range(Nm ** K):
        condMean = mu * (1 - np.sum(A)) + np.dot(A, X[:, ii])
        xPDF = X1 - condMean  # Correctly compute xPDF

        # Replace GaussianMixture with direct Gaussian PDF evaluation
        tmp = norm.pdf(xPDF, loc=muC, scale=sigmaC)  # Evaluate Gaussian PDF
        q = W * tmp  # Scale by quadrature weights
        q[q < kappa] = kappa  # Numerical stability

        if nMoments == 1:  # Match only the first moment
            P1[ii, :], _, _ = discrete_approximation(
                X1, lambda x: np.array((x - condMean) / scalingFactor),
                TBar[:1] / scalingFactor,
                q=q, lambda0=np.zeros(1))
        else:  # Match 2 moments first
            p, lambda_val, momentError = discrete_approximation(
                X1, lambda x: np.vstack([(x - condMean) / scalingFactor,
                                          ((x - condMean) / scalingFactor) ** 2]),
                TBar[:2] / scalingFactor ** np.arange(1, 3),
                q=q, lambda0=np.zeros(2))
            match_order[ii] = 2.0
            if np.linalg.norm(momentError) > 1e-4:  # If 2 moments fail, match 1 moment
                warnings.warn('Failed to match first 2 moments. Just matching 1.')
                P1[ii, :], _, _ = discrete_approximation(
                    X1, lambda x: [(x - condMean) / scalingFactor],
                    TBar[:1] / scalingFactor,
                    q=q, lambda0=np.zeros(1))
                match_order[ii] = 1.0
            else:
                P1[ii, :] = p

        P[ii, :] = np.kron(P1[ii, :], P2[ii, :])

    eigs = np.max(np.abs(np.linalg.eigvals(P)))
    output_string = "\033[1;33m Max eig {:.5f}:\033[0m".format(eigs)
    print(output_string)
    print('%s %10.5f' % ("\033[1;33m Averaged Matched Order " + country + ":\033[0m", np.mean(match_order)))
    res = np.linalg.norm(momentError)
    return P, X, match_order, res
