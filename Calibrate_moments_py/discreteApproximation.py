# %%
import numpy as np
from scipy.optimize import minimize, basinhopping
# %%
# def entropy_objective(lambda_, Tx, TBar, q):
#     L, N = Tx.shape
#     if len(lambda_) != L or len(TBar) != L or len(q) != N:
#         raise ValueError('Dimensions of inputs are not compatible.')

#     Tdiff = Tx - TBar[:, None]
#     temp = q * np.exp(lambda_.T @ Tdiff)
#     obj = np.sum(temp)

#     temp2 = temp * Tdiff
#     gradObj = np.sum(temp2, axis=1)

#     return obj, gradObj


# def entropy_objective(lambda_, Tx, TBar, q):
#     L, N = Tx.shape
#     if len(lambda_) != L or len(TBar) != L or len(q) != N:
#         raise ValueError('Dimensions of inputs are not compatible.')

#     Tdiff = Tx - TBar[:, None]
#     log_temp = np.log(q) + lambda_.T @ Tdiff
#     max_log_temp = np.max(log_temp)
#     exp_log_temp = np.exp(log_temp - max_log_temp)

#     obj = np.sum(exp_log_temp) * np.exp(max_log_temp)

#     gradObj = np.sum(exp_log_temp * Tdiff, axis=1)

#     return obj, gradObj
# def entropy_objective(lambda_, Tx, TBar, q):
#     L, N = Tx.shape
#     if len(lambda_) != L or len(TBar) != L or len(q) != N:
#         raise ValueError('Dimensions of inputs are not compatible.')

#     Tdiff = Tx - TBar[:, None]
#     log_temp = np.log(q) + lambda_.T @ Tdiff
#     max_log_temp = np.max(log_temp)
#     shifted_log_temp = log_temp - max_log_temp
#     exp_shifted_log_temp = np.exp(shifted_log_temp)
#     sum_exp_shifted_log_temp = np.sum(exp_shifted_log_temp)

#     if sum_exp_shifted_log_temp == 0:
#         obj = max_log_temp
#     else:
#         obj = max_log_temp + np.log(sum_exp_shifted_log_temp)

#     gradObj = np.sum(exp_shifted_log_temp * Tdiff, axis=1)

#     if sum_exp_shifted_log_temp == 0:
#         hessianObj = np.zeros((len(lambda_), len(lambda_)))
#     else:
#         if len(lambda_) > 1:
#             hessianObj = exp_shifted_log_temp @ Tdiff.T @ Tdiff / N
#         else:
#             hessianObj = np.array([[np.sum(exp_shifted_log_temp * Tdiff * Tdiff)]]) / N

#     return obj, gradObj, hessianObj

def entropy_objective(lambda_, Tx, TBar, q):
    L, N = Tx.shape
    if len(lambda_) != L or len(TBar) != L or len(q) != N:
        raise ValueError('Dimensions of inputs are not compatible.')
    # Tdiff = Tx-repmat(TBar,1,N);
    TBar_repeated = np.expand_dims(TBar, axis=1)  # Expand dimensions to (L, 1)
    Tdiff= Tx- np.repeat(TBar_repeated, N, axis=1)  # Repeat along axis 1 N times
    # Transform to avoid overflow of exp
    # i.e. add to log of q the product of the lambda_ matrix and Tdiff
    temp_log = np.log(np.reshape(q, (1, N))) + lambda_.T @ Tdiff
    # max_log_val = np.log(np.finfo('d').max) - 10.0  # Maximum threshold for overflow in logarithmic scale
    # max_pos = np.argmax(temp_log, axis=1)

    # if temp_log[0, max_pos] >= max_log_val:
    #     temp_log[0, max_pos] = max_log_val

    # # Apply logarithmic sum trick
    # temp_max = np.max(temp_log)
    # temp_log_shifted = temp_log - temp_max
    # temp = np.exp(temp_log_shifted)
    # obj = np.log(np.sum(temp)) + temp_max
    pen=1.0
    if np.log(N)+np.max(temp_log)>=np.log(np.finfo('d').max): 
        temp_log[0,:]=5.0
        pen=100.0
    temp=np.exp(temp_log)*pen
    obj=np.sum(temp)

    #Compute gradient of objective function

   
    temp2 =np.repeat(temp,L,axis=0)*Tdiff
    gradObj = np.sum(temp2,axis=1)
    

    # % Compute hessian of objective function

    
    hessianObj = np.matmul(temp2,Tdiff.T)
   
    return obj, gradObj,hessianObj



#%%
def discrete_approximation(D, T, TBar, q=None, lambda0=None):
    # N=D.shape[0]
    N = D.shape[0] if D.ndim == 1 else D.shape[1]

    Tx = np.array(T(D))
    L = Tx.shape[0]

    if Tx.shape[1] != N or len(TBar) != L:
        raise ValueError('Dimension mismatch')

    if q is None:
        q = np.ones(N) / N

    if lambda0 is None:
        lambda0 = np.zeros(L)
    # trust-constr
    lower_bound = -10  # Adjust the lower bound as needed
    upper_bound = -lower_bound  # Adjust the upper bound as needed

    # Define the bounds as a list of tuples
    bounds = None # [(lower_bound, upper_bound)] * len(lambda0)
    # Nelder-Mead: NM
    # Powell: Powell
    # CG: CG (Conjugate Gradient)
    # BFGS: BFGS (Broyden-Fletcher-Goldfarb-Shanno)
    # L-BFGS-B: L-BFGS-B (Limited-memory BFGS with Bounds)
    # TNC: TNC (Truncated Newton Conjugate-Gradient)
    # COBYLA: COBYLA (Constrained Optimization BY Linear Approximations)
    # SLSQP: SLSQP (Sequential Least SQuares Programming)
    # trust-constr: trust-constr (Trust Region Constrained)
    # dogleg: dogleg
    method='COBYLA'#'TNC'#'Newton-CG'#'CG'#'L-BFGS-B'

    try:
        res = minimize(entropy_objective, lambda0, 
                       args=(Tx, TBar, q),
                       options={'gtol': 1e-8, 'maxiter': 1000},
                         method=method, jac=True,bounds=bounds)
        lambdaBar = res.x
    except:
        print('\033[91m Failed to find a solution from provided initial guess. Trying new initial guess.\033[0m')
        res = minimize(entropy_objective, 
                       np.random.random(len(lambda0))/100.0,
                        options={'gtol': 1e-8, 'maxiter': 1000},
                          args=(Tx, TBar, q), method=method, jac=True,bounds=bounds)
        # the initial guess was np.zeros_like(lambda0)
        lambdaBar = res.x
    
        
    # print(f'\033[33m Residual entropy problem:\033[0m {res.fun}')
    obj, gradObj,_ = entropy_objective(lambdaBar, Tx, TBar, q)
    Tdiff = Tx - TBar[:, None]
    p = (q * np.exp(lambdaBar.T @ Tdiff)) / obj
    # See equation 2.5 of QE paper
    moment_error = gradObj / obj

    return p, lambdaBar, moment_error


# %%
