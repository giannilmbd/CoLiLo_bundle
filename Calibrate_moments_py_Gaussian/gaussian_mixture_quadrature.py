# %%
import numpy as np

def gaussian_mixture_quadrature(coeff, mu, sigma, n):
    if len(coeff) != len(mu) or len(coeff) != len(sigma) or len(mu) != len(sigma):
        raise ValueError("Coeff, Mu, Sigma must be of same length")

    coeff = np.asarray(coeff).flatten()
    # if coeff.shape[0] > coeff.shape[1]:
    #     coeff = coeff.T  # convert to row vector
    mu = np.asarray(mu).flatten()
    sigma = np.asarray(sigma).flatten()

    k = len(coeff)
    temp = np.zeros((2 * n + 1, k))
    sigma2 = sigma ** 2
    temp[0, :] = 1
    temp[1, :] = mu
    # for i in range(2, 2 * n):
    #     temp[i + 1, :] = mu * temp[i, :] + (i - 1) * sigma2 * temp[i - 1, :]
    for i in range(1, 2 * n):
        temp[i+1, :] = mu * temp[i, :] +(i) * sigma2 * temp[i - 1, :]
    poly_moments = temp @ coeff.T

    m = np.zeros((n + 1, n + 1))

    # In Matlab was
    # for n=1:N+1
    # M(n,:) = PolyMoments(n:N+n);
    # end
    for i in range(0, n+1):
        m[i , :] = poly_moments[i: n + i+1]
    if not np.allclose(m, m.T):
        print('Asymmetric matrix')
        print(np.max(m-m.T))
    r = np.linalg.cholesky(m)

    temp0 = np.diag(r)
    temp0 = np.delete(temp0, -1)
    beta = temp0[1:] / temp0[:-1]
    temp1 = np.diag(r, 1)
    temp2 = temp1 / temp0 
    alpha = temp2 - np.concatenate(([0], temp2[:-1]))

    t = np.diag(alpha) + np.diag(beta, -1) + np.diag(beta, 1)
    d, v = np.linalg.eig(t)

    x = np.sort(d)
    ind = np.argsort(d)

    w = np.zeros(n)
    for i in range(n):
        VV=v[:,i]
        w[i] = sum(coeff) * VV[1] ** 2 / np.dot(VV, VV)
    w = w[ind]

    return x, w

# %%
