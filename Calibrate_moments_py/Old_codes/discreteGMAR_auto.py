##########################################################################
# discreteGMAR
# (c) 2019 Alexis Akira Toda
#
# Purpose:
#       Discretize an AR(p) process with Gaussian mixture shocks
#
# Usage:
#       [P,X] = discreteGMAR(mu,A,pC,muC,sigmaC,Nm,nMoments,method,nSigmas)
#
# Inputs:
# mu        - unconditional mean
# A         - vector of coefficients of AR(p)
# pC        - vector of proportions of Gaussian mixtures
# muC       - vector of means of Gaussian mixtures
# sigmaC    - vector of standard deviations of Gaussian mixtures
# Nm        - number of grid points in one dimension
# Optional:
# nMoments  - number of moments to match (default = 2)
# method    - quadrature method (default = 'even')
# nSigmas	- grid spacing when using even-spaced grid

import numpy as np
import warnings
    
def discreteGMAR(mu = None,A = None,pC = None,muC = None,sigmaC = None,Nm = None,nMoments = None,method = None,nSigmas = None): 
    ## some error checking
    if np.any(pC < 0):
        raise Exception('mixture proportions must be positive')
    
    if np.any(sigmaC < 0):
        raise Exception('standard deviations must be positive')
    
    if sum(pC) != 1:
        raise Exception('mixture proportions must add up to 1')
    
    if pC.shape[1-1] < pC.shape[2-1]:
        pC = np.transpose(pC)
    
    if muC.shape[1-1] < muC.shape[2-1]:
        muC = np.transpose(muC)
    
    if sigmaC.shape[1-1] < sigmaC.shape[2-1]:
        sigmaC = np.transpose(sigmaC)
    
    K = len(A)
    if A.shape[1-1] > A.shape[2-1]:
        A = np.transpose(A)
    
    F = np.array([[A],[np.eye(K - 1,K)]])
    
    rho = np.abs(eigs(F,1))
    
    if rho >= 1:
        raise Exception('spectral radius must be less than one')
    
    # compute conditional moments
    sigmaC2 = sigmaC ** 2
    T1 = np.transpose(pC) * muC
    
    T2 = np.transpose(pC) * (muC ** 2 + sigmaC2)
    
    T3 = np.transpose(pC) * (muC ** 3 + np.multiply(3 * muC,sigmaC2))
    
    T4 = np.transpose(pC) * (muC ** 4 + np.multiply(6 * (muC ** 2),sigmaC2) + 3 * sigmaC2 ** 2)
    
    TBar = np.transpose(np.array([T1,T2,T3,T4]))
    nComp = len(pC)
    
    temp = np.zeros((1,1,nComp))
    temp[1,1,:] = sigmaC2
    gmObj = gmdistribution(muC,temp,pC)
    
    # Default number of moments is 2
    if len(varargin) == 6:
        nMoments = 2
    
    # Check that Nm is a valid number of grid points
    if not isnumeric(Nm)  or Nm < 3 or rem(Nm,1) != 0:
        raise Exception('Nm must be a positive integer greater than 3')
    
    # Check that nMoments is a valid number
    if not isnumeric(nMoments)  or nMoments < 1 or nMoments > 4 or not ((rem(nMoments,1) == 0) or (nMoments == 1)) :
        raise Exception('nMoments must be either 1, 2, 3, 4')
    
    # set default nSigmas if not supplied
    if len(varargin) < 9:
        if rho <= 1 - 2 / (Nm - 1):
            nSigmas = np.sqrt(2 * (Nm - 1))
        else:
            nSigmas = np.sqrt(Nm - 1)
    
    sigma = np.sqrt(T2 - T1 ** 2)
    
    temp = np.linalg.solve((np.eye(K ** 2) - kron(F,F)),np.eye(K ** 2))
    sigmaX = sigma * np.sqrt(temp(1,1))
    
    # construct the one dimensional grid
    if 'even' == method:
        X1 = np.linspace(mu - nSigmas * sigmaX,mu + nSigmas * sigmaX,Nm)
        W = np.ones((1,Nm))
    else:
        if 'gauss-legendre' == method:
            X1,W = legpts(Nm,np.array([mu - nSigmas * sigmaX,mu + nSigmas * sigmaX]))
            X1 = np.transpose(X1)
        else:
            if 'clenshaw-curtis' == method:
                X1,W = fclencurt(Nm,mu - nSigmas * sigmaX,mu + nSigmas * sigmaX)
                X1 = fliplr(np.transpose(X1))
                W = fliplr(np.transpose(W))
            else:
                if 'gauss-hermite' == method:
                    if rho > 0.8:
                        warnings.warn('Model is persistent; even-spaced grid is recommended')
                    X1,W = GaussHermite(Nm)
                    X1 = mu + np.sqrt(2) * sigma * np.transpose(X1)
                    W = np.transpose(W) / np.sqrt(np.pi)
                else:
                    if 'GMQ' == method:
                        if rho > 0.8:
                            warnings.warn('Model is persistent; even-spaced grid is recommended')
                        X1,W = GaussianMixtureQuadrature(pC,muC,sigmaC,Nm)
                        X1 = X1 + mu
    
    X = np.transpose(allcomb2(np.ones((K,1)) * X1))
    
    P = np.full([Nm ** K,Nm ** K],np.nan)
    
    P1 = np.full([Nm ** K,Nm],np.nan)
    
    P2 = kron(np.eye(Nm ** (K - 1)),np.ones((Nm,1)))
    
    scalingFactor = np.amax(np.abs(X1))
    kappa = 1e-08
    for ii in np.arange(1,Nm ** K+1).reshape(-1):
        condMean = mu * (1 - sum(A)) + A * X(:,ii)
        xPDF = np.transpose((X1 - condMean))
        if 'gauss-hermite' == method:
            q = np.multiply(W,np.transpose((pdf(gmObj,xPDF) / normpdf(xPDF,0,sigma))))
        else:
            if 'GMQ' == method:
                q = np.multiply(W,np.transpose((pdf(gmObj,xPDF) / pdf(gmObj,np.transpose(X1)))))
            else:
                q = np.multiply(W,np.transpose((pdf(gmObj,xPDF))))
        if np.any(q < kappa):
            q[q < kappa] = kappa
        if nMoments == 1:
            P1[ii,:] = discreteApproximation(X1,lambda x = None: (x - condMean) / scalingFactor,TBar(1) / scalingFactor,q,0)
        else:
            p,lambda_,momentError = discreteApproximation(X1,lambda x = None: np.array([[(x - condMean) / scalingFactor],[((x - condMean) / scalingFactor) ** 2]]),TBar(np.arange(1,2+1)) / (np.transpose(scalingFactor ** (np.arange(1,2+1)))),q,np.zeros((2,1)))
            if norm(momentError) > 1e-05:
                warnings.warn('Failed to match first 2 moments. Just matching 1.')
                P1[ii,:] = discreteApproximation(X1,lambda x = None: (x - condMean) / scalingFactor,TBar(1) / scalingFactor,q,0)
            else:
                if nMoments == 2:
                    P1[ii,:] = p
                else:
                    if nMoments == 3:
                        pnew,__,momentError = discreteApproximation(X1,lambda x = None: np.array([[(x - condMean) / scalingFactor],[((x - condMean) / scalingFactor) ** 2],[((x - condMean) / scalingFactor) ** 3]]),TBar(np.arange(1,3+1)) / (np.transpose(scalingFactor ** (np.arange(1,3+1)))),q,np.array([[lambda_],[0]]))
                        if norm(momentError) > 1e-05:
                            warnings.warn('Failed to match first 3 moments.  Just matching 2.')
                            P1[ii,:] = p
                        else:
                            P1[ii,:] = pnew
                    else:
                        pnew,__,momentError = discreteApproximation(X1,lambda x = None: np.array([[(x - condMean) / scalingFactor],[((x - condMean) / scalingFactor) ** 2],[((x - condMean) / scalingFactor) ** 3],[((x - condMean) / scalingFactor) ** 4]]),TBar / (np.transpose(scalingFactor ** (np.arange(1,4+1)))),q,np.array([[lambda_],[0],[0]]))
                        if norm(momentError) > 1e-05:
                            #warning('Failed to match first 4 moments.  Just matching 3.')
                            pnew,__,momentError = discreteApproximation(X1,lambda x = None: np.array([[(x - condMean) / scalingFactor],[((x - condMean) / scalingFactor) ** 2],[((x - condMean) / scalingFactor) ** 3]]),TBar(np.arange(1,3+1)) / (np.transpose(scalingFactor ** (np.arange(1,3+1)))),q,np.array([[lambda_],[0]]))
                            if norm(momentError) > 1e-05:
                                warnings.warn('Failed to match first 3 moments.  Just matching 2.')
                                P1[ii,:] = p
                            else:
                                P1[ii,:] = pnew
                                warnings.warn('Failed to match first 4 moments.  Just matching 3.')
                        else:
                            P1[ii,:] = pnew
            P[ii,:] = kron(P1(ii,:),P2(ii,:))
    
    return P,X