from scipy.optimize import minimize,basinhopping
from calculate_moments import calculate_moments
import warnings
import numpy as np
warnings.filterwarnings('ignore') 
def match_moments(empirical_moments, muC_init, pC_init, sigmaC_init,method='an',country='NaC'):
    # print("\033[34m"+country+"\033[0m \n")
    def objective(params):
        muC = params[0]
        pC =1.0# np.abs(params[3])
        sigmaC = np.abs(params[1])
        # if np.sum(pC)>=1:
        #     res_pC=1-np.sum(pC)
        #     pC=pC+res_pC/2.0
        # pC=np.append(pC, 1-np.sum(pC))
        # pC=pC/np.sum(pC)

        # Compute the moments based on the current parameter values
        _, _, _, _, TBar,var_,Sk,Kur = calculate_moments(pC, muC, sigmaC)

        # Compute the L2 norm between the calculated moments and empirical moments
        d1=abs(TBar[0]-empirical_moments[0])
        
        d2=np.linalg.norm(np.array([var_,Kur])/empirical_moments[[1,3]]-1.0)
        

        distance = np.linalg.norm(np.array([d1,d2]))
        # print(distance)
        return distance
    
    def jacobian(params):
        muC = params[0]
        
        sigmaC = params[5]
        

        # Compute the moments based on the current parameter values
        _, _, _, _, TBar = calculate_moments(pC, muC, sigmaC)

        # Compute the Jacobian matrix
        dTBar_dmuC = np.array([
            [np.inner(pC, [1, 0, 0]), np.inner(pC, [0, 1, 0]), np.inner(pC, [0, 0, 1]), 0, 0, 0, 0, 0, 0],
            [2 * np.inner(pC, muC), 0, 0, np.inner([muC[0], muC[1], muC[2]], [1, 0, 0]), np.inner([muC[0], muC[1], muC[2]], [0, 1, 0]), np.inner([muC[0], muC[1], muC[2]], [0, 0, 1]), 0, 0, 0],
            [3 * np.inner(pC, np.power(muC, 2)), 0, 0, 3 * np.inner([muC[0], muC[1], muC[2]], np.power([muC[0], muC[1], muC[2]], 2)), 0, 0, 3 * np.inner([muC[0], muC[1], muC[2]], np.power(sigmaC, 2)), 0, 0],
            [4 * np.inner(pC, np.power(muC, 3)), 0, 0, 6 * np.inner([muC[0], muC[1], muC[2]], np.power([muC[0], muC[1], muC[2]], 2)), 0, 0, 3 * np.inner([muC[0], muC[1], muC[2]], np.power(sigmaC, 2)), 0, 0]
        ])

        dTBar_dpC = np.array([
            [muC[0], muC[1], muC[2], 0, 0, 0, 0, 0, 0],
            [(np.power(muC[0], 2) + np.power(sigmaC[0], 2)), (np.power(muC[1], 2) + np.power(sigmaC[1], 2)), (np.power(muC[2], 2) + np.power(sigmaC[2], 2)), muC[0], muC[1], muC[2], 0, 0, 0],
            [3 * muC[0] * np.power(sigmaC[0], 2), 3 * muC[1] * np.power(sigmaC[1], 2), 3 * muC[2] * np.power(sigmaC[2], 2), 3 * np.power(muC[0], 2) * sigmaC[0], 3 * np.power(muC[1], 2) * sigmaC[1], 3 * np.power(muC[2], 2) * sigmaC[2], 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, sigmaC[0], sigmaC[1], sigmaC[2]],
            [0, 0, 0, 0, 0, 0, 2 * muC[0] * sigmaC[0], 2 * muC[1] * sigmaC[1], 2 * muC[2] * sigmaC[2]]
        ])

        dTBar_dsigmaC = np.array([
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0]
        ])

        jacobian_matrix = np.concatenate((dTBar_dmuC, dTBar_dpC, dTBar_dsigmaC), axis=0)

        return jacobian_matrix






    # Define the initial parameter values
    initial_params = list(muC_init)+list(sigmaC_init)

    # Set the bounds for the optimization variables
    bounds = [(None, None)] * len(muC_init)  + [(1E-8, None)] * len(sigmaC_init)
    
    # Run the optimization to find the optimal parameter values
    # method:
    #     'Nelder-Mead': Nelder-Mead simplex algorithm
    # 'Powell': Powell's conjugate direction method
    # 'CG': Conjugate Gradient method
    # 'BFGS': Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
    # 'L-BFGS-B': Limited-memory BFGS with bounds
    # 'TNC': Truncated Newton algorithm
    # 'COBYLA': Constrained Optimization by Linear Approximation
    # 'SLSQP': Sequential Least SQuares Programming
    # 'trust-constr': Trust-region constrained optimization
    # Define the callback function

    # THIS IS NOT USED ANYMORE SINCE p3=1-sum(pC)

    def callback(xk):
        # Print the current iteration and parameter values
        print(f"Iteration: {callback.iteration}")
        print(f"Parameters: {xk}")
        print("-------------------------")
        callback.iteration += 1
    # constraints = [
    # {'type': 'eq', 'fun': sum_constraint}
    # ]
# Set the initial iteration count
    def constraint_mu(x):  # for probabilities
        return 1.0 if all(ele > 1e-8 for ele in x) else -1.0  # positivity


    
    # def constraint_pc(x): # for probabilities
    #     return 0.99999 - np.sum(np.abs(x[3:5]))# the upper limit is not included
    # const1=[{"fun": constraint_pc, "type": "ineq"},
    #         {"fun": constraint_mu, "type": "ineq"}]
    callback.iteration = 1
    # Optimize the objective function bounds=bounds,
    if method != 'an':
        result = minimize(objective, initial_params,
                       method='COBYLA',#L-BFGS-B
                    #    callback=callback,
                       tol=1e-14,
                       options={'gtol': 1e-14, 'disp': False}
                       )
        # ANNEALING/replacement basinhopping
    else:    
        # bounds = bounds,
        minimizer_kwargs = dict(bounds = bounds,
                                method = 'L-BFGS-B',tol= 1e-14)#,
                                #  constraints=(const1))


        result = basinhopping(objective, initial_params,
                              minimizer_kwargs=minimizer_kwargs,
                              niter=100,
                              T=5,
                              stepsize=1e-3,
                              interval=10,
                              disp=False,
                              niter_success=50,
                              seed=1984)

    # Retrieve the optimal parameter values
    muC_optimal, pC_optimal, sigmaC_optimal = result.x[:len(muC_init)], np.abs(result.x[len(muC_init):-len(sigmaC_init)]), np.abs(result.x[-len(sigmaC_init):])
    minimized_value = result.fun

    if np.sum(pC_optimal)>1:
        pC_optimal=pC_optimal/(np.sum(pC_optimal)+1E-8)
    pC=np.append(np.array(pC_optimal),1-np.sum(pC_optimal))
    pC_optimal=np.abs(pC/np.sum(pC))
    print("\033[44m \033[1;33m Country: " + country + "\033[0m")
    print(f"\033[91m Minimized Value:\033[0m {result.fun}")
    print(f"\033[91m Optimal muC:\033[0m {muC_optimal}")
    print(f"\033[91m Optimal pC:\033[0m {pC_optimal}")
    print(f"\033[91m Optimal sigmaC:\033[0m {sigmaC_optimal}")
    # print(f"\033[91m Sum of optimal pC:\033[0m {np.sum(pC_optimal)}")
    return muC_optimal, pC_optimal, sigmaC_optimal,minimized_value
