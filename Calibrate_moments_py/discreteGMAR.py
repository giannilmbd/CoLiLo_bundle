#%%
import numpy as np
from numpy.linalg import eigh
from scipy.stats import norm, multivariate_normal
from scipy.optimize import minimize
from scipy.linalg import eigvals
from itertools import product
import warnings
from sklearn.mixture import GaussianMixture
from discreteApproximation import discrete_approximation
from gaussian_mixture_quadrature import gaussian_mixture_quadrature
from calculate_moments import calculate_moments


def discrete_gmar(mu, A, pC, muC, sigmaC, Nm,nMoments=2, method='even', nSigmas=None,country=None):
    # some error checking
    print(f"\033[1;34m \033[43m Starting Discretization \033[0m")

    if any(p < 0 for p in pC):
        raise ValueError('mixture proportions must be positive')
    if any(s < 0 for s in sigmaC):
        raise ValueError('standard deviations must be positive')
    if sum(pC)>1 or sum(pC)<0.999999:
        print(pC)
        print(sum(pC))
        raise ValueError('mixture proportions must sum to 1')
        
    K = np.size(A)  # number of states
    Nm = int(Nm)  # number of discrete values for each state
    if nSigmas is None:
        nSigmas = 4
    if nMoments < 1 or nMoments > 4:
        raise ValueError('nMoments must be between 1 and 4')

    # Define a function to calculate the PDF
    def gm_pdf(x, gmm):
        return np.sum([w * multivariate_normal(mean=mu, cov=sigma).pdf(x)
                       for w, mu, sigma in zip(gmm.weights_init,
                                                gmm.means_init, gmm.covariances_)])

    # Create identity matrix of size K-1 K  and stack it under A
    I = np.eye(K)
    if K == 1:
        F = A
        rho = A
    else:
        # F=np.vstack((A,I[0:K-2,:]))
        # rho=np.linalg.eigvals(F)    
        F = np.vstack((A, np.eye(K - 1, K)))
        rho = np.abs(eigvals(F))  # Compute the spectral radius of F
        rho=rho[0]

    
    if np.max(np.abs(rho)) >= 1:
        raise ValueError('AR(1) coefficient must be less than 1 in absolute value')
    # compute conditional moments
    sigmaC2 = np.power(sigmaC,2)
    T1, T2, T3, T4, TBar,_,_,_ = calculate_moments(pC,muC,sigmaC)
    # T1 = np.inner(pC, muC)  # mean
    # T2 = np.inner(pC ,(np.power(muC,2) + sigmaC2))  # uncentered second moment
    # T3 = np.inner(pC , (np.power(muC,3) + 3*np.multiply(muC,sigmaC2)))  # uncentered third moment
    # T4 = np.inner(pC , (np.power(muC,4) + 6*(np.power(muC,2))*sigmaC2 + 3*np.power(sigmaC2,2)))  # uncentered fourth moment

    # TBar = np.array([T1, T2, T3, T4])

    nComp = len(pC)
    temp = np.zeros((1,1,nComp))
    temp[0,0,:] = sigmaC2
    # define the Gaussian mixture object
    gmObj = GaussianMixture(n_components=nComp, 
                            covariance_type='full',
                              means_init=muC.reshape(-1, 1),
                                weights_init=pC,
                                init_params='kmeans',
                                precisions_init=None)
    


    
    # Default number of moments is 2
    if nMoments is None:
        nMoments = 2

    # Check that Nm is a valid number of grid points
    if not isinstance(Nm, int) or Nm < 3:
        raise ValueError('Nm must be a positive integer greater than 3')

    # Check that nMoments is a valid number
    if nMoments < 1 or nMoments > 4:
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
        X1, W = gaussian_mixture_quadrature(pC, muC, sigmaC, Nm)
        X1 = X1 + mu

    X = np.array(list(product(X1, repeat=K))).T  # K * Nm^K matrix of grid points

    P = np.empty((Nm ** K, Nm ** K))  # transition probability matrix
    P1 = np.empty((Nm ** K, Nm))  # matrix to store transition probability
    P2 = np.kron(np.eye(Nm ** (K - 1)), np.ones((Nm, 1)))  # Nm^K * Nm^(K-1) matrix used to construct P
    scalingFactor = np.max(np.abs(X1))
    kappa = 1e-8

    match_order=np.empty(Nm ** K)
    for ii in range(Nm ** K):
        condMean = mu * (1 - np.sum(A)) + np.dot(A, X[:, ii])
        xPDF = (X1 - condMean)
        gmObj.fit(xPDF.reshape(-1, 1))

        
        if method == 'gauss-hermite':
            q = W * (gmObj.score_samples(xPDF.reshape(-1, 1)) / norm.pdf(xPDF, 0, sigma))
        elif method == 'GMQ':
            q = W * (gmObj.score_samples(xPDF.reshape(-1, 1)) / gmObj.score_samples(X1.reshape(-1, 1)))
        else:
            tmp=gmObj.score_samples(xPDF.reshape(-1, 1))
            q = W * gmObj.score_samples(xPDF.reshape(-1, 1))

        q[q < kappa] = kappa

        if nMoments == 1:  # match only 1 moment
            P1[ii, :], _, _ = discrete_approximation(X1, lambda x: [(x - condMean) / scalingFactor],
                                                    TBar[:1] / scalingFactor,
                                                    q=q, lambda0=np.zeros(1))
        else:  # match 2 moments first
            p, lambda_val, momentError = discrete_approximation(X1, lambda x: [(x - condMean) / scalingFactor,
                                                                            ((x - condMean) / scalingFactor) ** 2],
                                                            TBar[:2] / scalingFactor ** np.arange(1, 3),
                                                            q=q, lambda0=np.zeros(2))

            if np.linalg.norm(momentError) > 1e-5:  # if 2 moments fail, then just match 1 moment
                # warnings.warn('Failed to match first 2 moments. Just matching 1.')
                P1[ii, :], _, _ = discrete_approximation(X1, lambda x: [(x - condMean) / scalingFactor],
                                                        TBar[:1] / scalingFactor,
                                                        q=q, lambda0=np.zeros(1))
                        
                match_order[ii] =1.0
            else:
                P1[ii, :] = p

            if nMoments >= 3:
                p_new, _, momentError = discrete_approximation(X1, lambda x: [(x - condMean) / scalingFactor,
                                                                            ((x - condMean) / scalingFactor) ** 2,
                                                                            ((x - condMean) / scalingFactor) ** 3],
                                                            TBar[:3] / scalingFactor ** np.arange(1, 4),
                                                            q=q, lambda0=np.concatenate((lambda_val, [0])))

                if np.linalg.norm(momentError) > 1e-5:
                    # warnings.warn('Failed to match first 3 moments. Just matching 2.')
                    P1[ii, :] = p
                else:
                    P1[ii, :] = p_new

            if nMoments >= 4:
                p_new, _, momentError = discrete_approximation(X1, lambda x: [(x - condMean) / scalingFactor,
                                                                            ((x - condMean) / scalingFactor) ** 2,
                                                                            ((x - condMean) / scalingFactor) ** 3,
                                                                            ((x - condMean) / scalingFactor) ** 4],
                                                            TBar[:4] / scalingFactor ** np.arange(1, 5),
                                                            q=q, lambda0=np.concatenate((lambda_val, [0, 0])))

                match_order[ii] =4.0
                if np.linalg.norm(momentError) > 1e-5:
                    # warnings.warn('Failed to match first 4 moments. Just matching 3.')
                    p_new, _, momentError = discrete_approximation(X1, lambda x: [(x - condMean) / scalingFactor,
                                                                                ((x - condMean) / scalingFactor) ** 2,
                                                                                ((x - condMean) / scalingFactor) ** 3],
                                                                TBar[:3] / scalingFactor ** np.arange(1, 4),
                                                                q=q, lambda0=np.concatenate((lambda_val, [0])))

                    match_order[ii] =3.0
                    if np.linalg.norm(momentError) > 1e-5:
                        # warnings.warn('Failed to match first 3 moments. Just matching 2.')
                        P1[ii, :] = p
                        match_order[ii] =2.0
                    else:
                        P1[ii, :] = p_new
                else:
                    P1[ii, :] = p_new

        P[ii, :] = np.kron(P1[ii, :], P2[ii, :])
    eigs=np.max(np.abs(np.linalg.eigvals(P)  )  )
    output_string = "\033[1;33m Max eig {:.5f}:\033[0m".format(eigs)
    print(output_string)

    # print(f"\033[1;34m Averaged Matched Order"+ country+": {np.mean(match_order)}\033[0m")
    # print('%d %s cost $%.2f' % (6, 'bananas', 1.74))
    print('%s %10.5f' % ("\033[1;33m Averaged Matched Order "+ country+":\033[0m",np.mean(match_order)))
    res=np.linalg.norm(momentError)
    return P,X,match_order,res 






# %%
