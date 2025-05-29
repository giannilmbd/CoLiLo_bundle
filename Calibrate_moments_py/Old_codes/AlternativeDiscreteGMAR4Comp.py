import numpy as np
from itertools import product
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from discreteApproximation import discrete_approximation
from gaussian_mixture_quadrature import gaussian_mixture_quadrature

def discreteGMAR(mu, A, pC, muC, sigmaC, Nm, nMoments=2, method='even', nSigmas=None):
    # Error checking
    if any(pC < 0): 
        raise ValueError('mixture proportions must be positive')
    if any(sigmaC < 0):
        raise ValueError('standard deviations must be positive')
    if np.sum(pC) != 1:
        raise ValueError('mixture proportions must add up to 1')

    pC = np.squeeze(pC)
    muC = np.squeeze(muC)
    sigmaC = np.squeeze(sigmaC)

    K = len(A)
    A = np.squeeze(A)

    F = np.vstack((A, np.eye(K - 1, K)))  # matrix to represent AR(p) by VAR(1)
    rho = np.abs(np.linalg.eigvals(F)).max()  # spectral radius of F

    if rho >= 1:
        raise ValueError('spectral radius must be less than one')

    # Compute conditional moments
    sigmaC2 = sigmaC ** 2
    T1 = np.dot(pC, muC)  # mean
    T2 = np.dot(pC, muC ** 2 + sigmaC2)  # uncentered second moment
    T3 = np.dot(pC, muC ** 3 + 3 * muC * sigmaC2)  # uncentered third moment
    T4 = np.dot(pC, muC ** 4 + 6 * muC ** 2 * sigmaC2 + 3 * sigmaC2 ** 2)  # uncentered fourth moment

    TBar = np.array([T1, T2, T3, T4])

    nComp = len(pC)  # number of mixture components
    temp = np.zeros((1, 1, nComp))
    temp[0, 0, :] = sigmaC2
    gmObj = GaussianMixture(n_components=nComp, means_init=muC.reshape(-1, 1),
                            precisions_init=np.linalg.inv(temp.squeeze()), weights_init=pC)

    # Default number of moments is 2
    if nMoments is None:
        nMoments = 2

    # Check that Nm is a valid number of grid points
    if not isinstance(Nm, int) or Nm < 3:
        raise ValueError('Nm must be a positive integer greater than 3')

    # Check that nMoments is a valid number
    if not isinstance(nMoments, int) or nMoments < 1 or nMoments > 4 or not (nMoments == 1 or nMoments.is_integer()):
        raise ValueError('nMoments must be either 1, 2, 3, 4')

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
        X1 = np.linspace(mu - nSigmas * sigmaX, mu + nSigmas * sigmaX, Nm)
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
        X1, W = GaussianMixtureQuadrature(pC, muC, sigmaC, Nm)
        X1 = X1 + mu

    X = np.array(list(product(X1, repeat=K))).T  # K * Nm^K matrix of grid points

    P = np.empty((Nm ** K, Nm ** K))  # transition probability matrix
    P1 = np.empty((Nm ** K, Nm))  # matrix to store transition probability
    P2 = np.kron(np.eye(Nm ** (K - 1)), np.ones((Nm, 1)))  # Nm^K * Nm^(K-1) matrix used to construct P
    scalingFactor = np.max(np.abs(X1))
    kappa = 1e-8

    for ii in range(Nm ** K):
        condMean = mu * (1 - np.sum(A)) + np.dot(A, X[:, ii])
        xPDF = (X1 - condMean)
        
        if method == 'gauss-hermite':
            q = W * (gmObj.score_samples(xPDF.reshape(-1, 1)) / norm.pdf(xPDF, 0, sigma))
        elif method == 'GMQ':
            q = W * (gmObj.score_samples(xPDF.reshape(-1, 1)) / gmObj.score_samples(X1.reshape(-1, 1)))
        else:
            q = W * gmObj.score_samples(xPDF.reshape(-1, 1))

        q[q < kappa] = kappa

        if nMoments == 1:  # match only 1 moment
            P1[ii, :] = discreteApproximation(X1, lambda x: (x - condMean) / scalingFactor, TBar[0] / scalingFactor,
                                               q
        else:  # match 2 moments first 
            p, lambda_val, momentError =  discreteApproximation(X1, lambda x: [(x - condMean) / scalingFactor,
                                                                               ((x - condMean) / scalingFactor) ** 2],
                                                               TBar[:2] / scalingFactor ** np.arange(1, 3),
                                                               q, np.zeros(2))
            if np.linalg.norm(momentError) > 1e-5:  # if 2 moments fail, then just match 1 moment
                print('Failed to match first 2 moments. Just matching 1.')
                P1[ii, :] = discreteApproximation(X1, lambda x: (x - condMean) / scalingFactor,
                                                  TBar[0] / scalingFactor, q, 0)
            elif nMoments == 2:
                P1[ii, :] = p
            elif nMoments == 3:
                p_new, _, momentError = discreteApproximation(X1,
                                                              lambda x: [(x - condMean) / scalingFactor,
                                                                         ((x - condMean) / scalingFactor) ** 2,
                                                                         ((x - condMean) / scalingFactor) ** 3],
                                                              TBar[:3] / scalingFactor ** np.arange(1, 4),
                                                              q, np.array([lambda_val, 0]))
                if np.linalg.norm(momentError) > 1e-5:
                    print('Failed to match first 3 moments. Just matching 2.')
                    P1[ii, :] = p
                else:
                    P1[ii, :] = p_new
            else:  # 4 moments
                p_new, _, momentError = discreteApproximation(X1,
                                                              lambda x: [(x - condMean) / scalingFactor,
                                                                         ((x - condMean) / scalingFactor) ** 2,
                                                                         ((x - condMean) / scalingFactor) ** 3,
                                                                         ((x - condMean) / scalingFactor) ** 4],
                                                              TBar / scalingFactor ** np.arange(1, 5),
                                                              q, np.array([lambda_val, 0, 0]))
                if np.linalg.norm(momentError) > 1e-5:
                    print('Failed to match first 4 moments. Just matching 3.')
                    p_new, _, momentError = discreteApproximation(X1,
                                                                  lambda x: [(x - condMean) / scalingFactor,
                                                                             ((x - condMean) / scalingFactor) ** 2,
                                                                             ((x - condMean) / scalingFactor) ** 3],
                                                                  TBar[1:4] / scalingFactor ** np.arange(2, 5),
                                                                  q, np.array([lambda_val, 0]))
                    if np.linalg.norm(momentError) > 1e-5:
                        print('Failed to match first 3 moments. Just matching 2.')
                        P1[ii, :] = p
                    else:
                        P1[ii, :] = p_new
                        print('Failed to match first 4 moments. Just matching 3.')
                else:
                    P1[ii, :] = p_new

        P[ii, :] = np.kron(P1[ii, :], P2[ii, :])
