import numpy as np
import warnings

def calculate_moments(pC, muC, sigmaC):
    # warnings.filterwarnings("error", category=RuntimeWarning)
    pC=np.array(pC)
    muC=np.array(muC)
    sigmaC=np.array(sigmaC)
    try:
        sigmaC=np.abs(sigmaC)
        pC=np.abs(pC)
        sigmaC2 = np.power(sigmaC, 2.0)
        T1 = ( muC)  # mean
        T2 =  (muC** 2.0+sigmaC2)  # uncentered second moment
        T3 = (muC** 3.0 + 3.0 *(muC* sigmaC2)) # uncentered third moment
        T4 = (muC**4.0 + 6.0 * muC**2.0 * sigmaC2 + 3.0 * sigmaC2**2.0)  # uncentered fourth moment
        
        var_=(T2-T1**2)
        Sk=(T3-3.0*T1*var_-T1**3.0)/(var_**(3.0/2.0)+(1E-10))
        Kur=(T4+6.0*T1**2.0*T2-3.0*T1**4.0-4.0*T1*T3)/(var_**2.0+(1E-10))
        TBar = np.array([T1, T2, T3, T4])
    except RuntimeWarning as e:
        print("Warning encountered:", e)
        # do it anyway
        sigmaC=np.abs(sigmaC)
        pC=np.abs(pC)
        sigmaC2 = np.power(sigmaC, 2.0)
        T1 = ( muC)  # mean
        T2 =  (muC** 2.0+sigmaC2)  # uncentered second moment
        T3 = (muC** 3.0 + 3.0 *(muC* sigmaC2)) # uncentered third moment
        T4 = (muC**4.0 + 6.0 * muC**2.0 * sigmaC2 + 3.0 * sigmaC2**2.0)  # uncentered fourth moment
        
        var_=(T2-T1**2)
        Sk=(T3-3.0*T1*var_-T1**3.0)/(var_**(3.0/2.0)+(1E-10))
        Kur=(T4+6.0*T1**2.0*T2-3.0*T1**4.0-4.0*T1*T3)/(var_**2.0+(1E-10))
        TBar = np.array([T1, T2, T3, T4])


    return T1, T2, T3, T4, TBar,var_, Sk, Kur 
