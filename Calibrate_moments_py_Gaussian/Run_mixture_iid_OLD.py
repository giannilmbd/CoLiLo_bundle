#%%
# This code calibrates a Gaussian Mixture model for the iid case (see Run_mixture.py for the AR case)
# The iid is necessary for EZ preferences
# Then it generates a discrete approximation of the mixture model
# using the discretization method by  Leland E. Farmer and Alexis Akira Toda Quantitative Economics 2017
# Components of the program are (defined as .py files):
# 1. This master code
# 2. match_moments:
#   This function is used to match the moments of the DGP to those computed from the data. 
#   The latter are obtained by: 
#           a) firsyt quadratically detrending real GDP
#           b) fitting an AR(1) process on the cyclical component 
#           c) storing the fitted AR(1) auto-regressor coefficient and the first four moments in a csv file (uploaded here)
# 3. discreteGMAR: This function is used to generate a discrete approximation of the Gaussian Mixture Model
# 4. Subfunctions
#   4.1 discreteApproximation: This function is used to generate a discrete approximation of the Gaussian
#   4.2 entropyObjective: This function is used to compute the entropy of the discrete approximation
#   4.3 calculate_moments: This function is used to compute the moments of the Gaussian Mixture Model
#   4.4 other auxiliary function might be defined within the main ones

# ----------------------------------------------------------------------------------------------------------------------------------------
# %% ------------------------------------  Load Packages ------------------------------------ #
import os

# Set the working directory to the script's directory
current_directory = os.path.dirname(os.path.realpath(__file__))
os.chdir(current_directory)
import numpy as np
import pandas as pd
from scipy.stats import norm
import warnings
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import matplotlib 
from multiprocessing import Pool
import statistics 
import scipy.stats as stats
from match_moments import match_moments
from discreteGMAR import discrete_gmar
from discreteApproximation import discrete_approximation
import seaborn as sns
warnings.filterwarnings('ignore') 

matplotlib.use("TkAgg") #
matplotlib.interactive(False)

#%% --------------------------------------------SET PARAMETERS -------------------------------------------------------- #
Nm = 9 #number of points for discretization (9 and 25 is what we considered so far)
#%% --------------------------------------------------------------------------------------------------------------------#

mu = np.array(0) # mean full DGP
A1 = np.array(0.9) # Autorecursive parameter
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

# How many moments to be used
n_moments = 4

# #%% ------------------------------------  MATCH MOMENTS ------------------------------------ #
# %%
# Define the empirical moments: the file ../Rcode/ar_data.csv contains this data

# |cc | autocorr_coef|     stdev|   skewness|  kurtosis|
# |:--|-------------:|---------:|----------:|---------:|
# |AR |     0.9479129| 0.0284724| -1.9177994| 13.996373|
# |AU |     0.5569621| 0.0094989| -4.6836885| 38.945032|
# |BR |     0.9226162| 0.0165565| -1.4789713| 14.933441|
# |CA |     0.7597379| 0.0149151| -4.2676538| 39.643465|
# |CH |     0.7296003| 0.0113501| -2.6742009| 24.940978|
# |CL |     0.8436627| 0.0189830| -4.4540204| 35.773801|
# |CN |     0.8958501| 0.0181662| -1.6755076| 32.667846|
# |CO |     0.8762621| 0.0232624| -4.6197382| 40.634206|
# |CZ |     0.9209191| 0.0141897| -2.5272684| 23.240009|
# |DK |     0.9044048| 0.0118900| -1.1531083| 14.005966|
# |GB |     0.6728827| 0.0266076| -5.9497241| 55.225264|
# |HK |     0.9358785| 0.0165439| -0.1376921|  5.148525|
# |HU |     0.9307925| 0.0206868| -3.1244316| 34.774842|
# |ID |     0.9448044| 0.0186062| -3.1174243| 17.794093|
# |IL |     0.7673916| 0.0150605| -2.7237805| 25.070432|
# |IN |     0.6868980| 0.0363311| -4.5791379| 37.469463|
# |JP |     0.8017004| 0.0131341| -2.9950362| 18.999216|
# |KR |     0.8246345| 0.0124755| -3.2790864| 20.317529|
# |MX |     0.7231970| 0.0232187| -5.6567823| 49.700015|
# |MY |     0.7488990| 0.0268648| -2.6931036| 24.480780|
# |NO |     0.8294434| 0.0120659| -1.5954190| 13.350384|
# |NZ |     0.7452435| 0.0182962| -0.1351289| 20.320388|
# |PE |     0.7988687| 0.0399104| -3.9588792| 44.009195|
# |PH |     0.8517558| 0.0203565| -5.5237712| 45.586729|
# |PL |     0.7062489| 0.0159854| -1.8896419| 13.530424|
# |SA |     0.7728816| 0.0201638|  0.7162878|  6.614967|
# |SE |     0.8395650| 0.0132043| -2.2580895| 21.398271|
# |SG |     0.8777042| 0.0229075| -2.3469160| 17.029254|
# |TH |     0.8693735| 0.0214180| -1.0001150|  8.952373|
# |TR |     0.7518451| 0.0351631| -0.2274531|  7.269394|
# |US |     0.8486052| 0.0124904| -2.8195127| 32.742859|
# |XM |     0.7754802| 0.0170041| -2.5119007| 34.057660|
# |ZA |     0.7098734| 0.0216946| -5.5479243| 53.286816|

#%% load data
import pandas as pd # Use pandas for data manipulation
# emp_moms_csv=pd.read_csv('../Rcode/ar_data.csv') # THIS HAS CRAZY KURTOSIS
# OR PWT TABLE 
emp_moms_csv=pd.read_csv('../Rcode/ar_data_pwt.csv')

country_id=emp_moms_csv['cc'].values # countries name
all_probs=np.empty((len(country_id),Nm*Nm)) # initialize matrices of probabilities and discretized shocks for each country
all_shocks=np.empty((len(country_id),Nm)) # shocks
all_pC=np.empty((len(country_id),3)) # probabs
all_muC=np.empty((len(country_id),3)) # mean
all_sigmaC=np.empty((len(country_id),3)) # stdev
all_match=np.empty((len(country_id),Nm)) # match order
all_resids=np.empty((len(country_id),1)) # residuals
tab_results=pd.DataFrame(data={'cc':country_id,'Residual':None,'AverageMatch':None}) # results
#%% ------------------------------------------------ Container to Run the estimation/calibration for each country -------------------------------------------------


def country_level(cnt_c,country_id,muC,pC,sigmaC,Nm,emp_moms_csv,n_moments): 
    cntry=country_id[cnt_c]
    # print(cntry+"\n")
    # Load the empirical moments of the AR process
    A1=0.0# THIS IS THE IID CASE emp_moms_csv[emp_moms_csv['cc']==cntry]['autocorr_coef'].values[0]
    emp_m1=0.0
    emp_m2=emp_moms_csv[emp_moms_csv['cc']==cntry]['stdev'].values[0]**2
    emp_m3=emp_moms_csv[emp_moms_csv['cc']==cntry]['skewness'].values[0]
    emp_m4=emp_moms_csv[emp_moms_csv['cc']==cntry]['kurtosis'].values[0]
    empirical_moments = np.array([emp_m1,emp_m2,emp_m3,emp_m4]) # Mean, variance, skewness, kurtosis (not -3)
    
    # Note that the discretization algoritm coputes moments of the iid part, ie the AR part does not appear in the moments per se
    # Hence we should provide empirical moments of the iid part of the estimated AR process (see quick_stylized_facts_pwt_only.r)



    # Define the initial values of muC, pC, and sigmaC --- to be found optimally
    muC_init =muC.copy()
    pC_init = np.array([pC[0], pC[1]])
    sigmaC_init = sigmaC.copy()

    # Match the moments and retrieve the optimal values: This function uses a minimzation function (to be fine tuened internally to that function)

    muC_optimal, pC_optimal, sigmaC_optimal,resid =match_moments(empirical_moments,
                            muC_init, pC_init,
                              sigmaC_init,method='an',country=cntry) #'an' is basinhopping (proble: hops could be outside the constraint)
    
    # make sure that probabilities are in the range [0,1] and that the sum of the probabilities is 1
    # pC_res=1-np.sum(pC_optimal)
    # note that inside match moments pC is normalized, hence the optimum must be normalized too
  



    # ------------------------------------------------------ DISCRETE GMAR ------------------------------------------------------ #

    p_even1, x_even1,match_order,_ = discrete_gmar(mu, A1, pC, muC,
                                                  sigmaC, Nm, n_moments, method='even',country=cntry)#'even''gauss-hermite'
    pC=np.array(pC_optimal)
    # convert results into arrays to be used in the discretization
    muC=np.array(muC_optimal)
    
    sigmaC=np.array(sigmaC_optimal)
    return p_even1, x_even1,match_order,resid,pC,muC,sigmaC
#%%  --------------------- for debugging purposes ---------------------
# i=int(np.where(country_id=="VEN")[0])
# a,b,c,d=country_level(i,country_id,muC,pC,sigmaC,Nm,emp_moms_csv,n_moments)
#%%
# Parallelize process by country
# first contruct the list of arguments to be passed to the function
items = [(i,country_id,muC,pC,sigmaC,Nm,emp_moms_csv,n_moments) for i in range(len(country_id))]
# then parallelize the process by country
with Pool() as mp_pool: #processes=1processes=1
    results = mp_pool.starmap(country_level, items)
#%% COLLECT THE RESULTS
for cnt_c in range(len(country_id)):
    # stuck the results country by country
    all_probs[cnt_c,:]=np.reshape(results[cnt_c][0],(1,Nm*Nm))
    all_shocks[cnt_c,:]=np.reshape(results[cnt_c][1],(1,Nm))
    all_match[cnt_c,:]=np.reshape(results[cnt_c][2],(1,Nm))
    resid=results[cnt_c][3]
    
    all_pC[cnt_c,:]=np.reshape(results[cnt_c][4],(1,3))
    all_muC[cnt_c,:]=np.reshape(results[cnt_c][5],(1,3))
    all_sigmaC[cnt_c,:]=np.reshape(results[cnt_c][6],(1,3))
    
    tab_results.loc[tab_results['cc']==country_id[cnt_c],'Residual']=resid
    all_resids[cnt_c]=resid
    tab_results.loc[tab_results['cc']==country_id[cnt_c],'AverageMatch']=statistics.mean(all_match[cnt_c,:])
    print("\033[34m Done country "+country_id[cnt_c]+"\033[0m")
# end loop for each country



#%% SAVE PROBABILITIES AND GRID ---------------------------------------------------- #
df1=pd.DataFrame(all_probs)
df1.to_csv("Pmatrices_iid_"+str(Nm)+".csv",
           index=False,header=False)
df1=pd.DataFrame(all_shocks)
df1.to_csv("Xmatrices_iid_"+str(Nm)+".csv", 
           index=False,header=False)
data = {
    'Country': country_id,
    'p1':all_pC[:,0] ,
    'p2':all_pC[:,1] ,
    'p3':all_pC[:,2] ,
    'mu1': all_muC[:,0],
    'mu2': all_muC[:,1],
    'mu3': all_muC[:,2],
    'sig1': all_sigmaC[:,0],
    'sig2': all_sigmaC[:,1],
    'sig3': all_sigmaC[:,1],
}
df = pd.DataFrame(data)
df.to_csv('mixed_gaussian_params.csv', index=False)

print("CSV file saved as 'mixed_gaussian_params.csv'")




# %% ------------------------------  SIMULATE SERIES ---------------------------------------------------- #
plt.figure()  
def simulate_series(P, X, T):
    series = []
    indx = np.random.choice(range(len(X)))  # Choose a random initial state
    current_state = X[indx,0]
    for _ in range(T):
        series.append(current_state)
        probabilities = P[indx,:]  # Get transition probabilities for the current state
        indx = np.random.choice(range(len(X)), p=probabilities)  # Choose next state based on transition probabilities
        current_state = X[indx,0]

    return np.array(series)
# %% ----------------------- SIMULATE ALL COUNTRIES ---------------------------------------------------- #
n = 6
Ncntr = len(country_id)
m = int(np.ceil(Ncntr / n))
tab_results['Mean']=np.nan 
tab_results['Variance']=np.nan 
tab_results['Skewness']=np.nan 
tab_results['Kurtosis']=np.nan 
fig_match, axes_match = plt.subplots(n, m)
# fig_match.suptitle("Discretized distribution")
# Iterate over each case and plot in the respective subplot
T=10000
case_idx = 0
for i in range(n):
    for j in range(m):
        if case_idx < Ncntr:
            Prb=np.reshape(all_probs[case_idx, :], (Nm, Nm))
            if not np.all(np.isclose(np.sum(Prb,axis=1),1.0)):
                raise ValueError("probabilities don't add up to 1.")
            Xsh=np.reshape(all_shocks[case_idx, :], (Nm,1))
            
            series=simulate_series(Prb, Xsh, T)
            
            tab_results.loc[tab_results['cc']==country_id[case_idx],'Mean']=np.mean(series[int(.9*T):])
            tab_results.loc[tab_results['cc']==country_id[case_idx],'Variance']=np.var(series[int(.9*T):])
            tab_results.loc[tab_results['cc']==country_id[case_idx],'Skewness']=stats.skew(series[int(.9*T):])
            tab_results.loc[tab_results['cc']==country_id[case_idx],'Kurtosis']=stats.kurtosis(series[int(.1*T):])+3.0
            ax = axes_match[i, j]  # Get the current subplot
            ##
            # Density Plot and Histogram of all arrival delays
            sns.histplot(ax=ax,data=series, kde=True,stat='probability',
                        color = 'darkblue', 
                        edgecolor='black',
                        line_kws={'lw': 4})

            ax.set_title(country_id[case_idx])  # Set subplot title
            case_idx += 1
            ax.set(xlabel='', ylabel='')
# Adjust the layout
fig_match.subplots_adjust(wspace=0.1, hspace=0.1)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=.4)

plt.draw()
plt.gcf().set_size_inches(10, 7)

plt.savefig('Approx_distros_iid.png', bbox_inches='tight', pad_inches = 0) #dpi=100,



print(tab_results)
#%% --------------- COMPARE WITH EMPIRICAL MOMENTS ------------------------------- #
## --- NOTE THAT EMPIRICAL MOMENTS ARE OBTAINED BY SIMULATING THE AR(1) PROCESS IN ../Rcodes/quick_stylized_facts.r ---

# USE rpy2 TO AVAIL OF R FUNCTIONS
import rpy2.robjects.packages as rpackages
import rpy2 
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects.lib.dplyr as dplyr
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.ipython.ggplot import image_png
from rpy2 import robjects
utils = importr('utils')
base = importr('base')
knitr = importr('knitr')
kextra=importr('kableExtra')
ggrepel=importr('ggrepel')
readr=importr('readr')

#%%
emp_new=emp_moms_csv 
emp_new['Variance_emp']=emp_moms_csv['stdev'].values**2


emp_new=emp_new.rename(columns={'skewness': 'Skewness_emp', 'kurtosis': 'Kurtosis_emp'}) 
emp_new.head()
with (ro.default_converter + pandas2ri.converter).context():
  r_emp_new = ro.conversion.get_conversion().py2rpy(emp_new)
  r_tab_results = ro.conversion.get_conversion().py2rpy(tab_results)
tab_all=dplyr.full_join(r_tab_results,r_emp_new,by='cc')
print(base.names(tab_all))
utils.head(tab_all)




#%% VARIANCE
p2=ggplot2.ggplot(data=tab_all)+\
ggplot2.geom_point(ggplot2.aes_string(x='Variance', y='Variance_emp'))+\
ggplot2.labs(title='Variance',x='Discretized', y='Empirical')+\
ggplot2.geom_line(ggplot2.aes_string(x='Variance', y='Variance'),color='red',linetype='dashed') +\
 ggrepel.geom_text_repel(ggplot2.aes_string(x='Variance', y='Variance_emp',label = 'cc'))+\
ggplot2.geom_smooth(ggplot2.aes_string(x='Variance', y='Variance_emp'),method='lm')
p2.plot()
# display(image_png(p))
robjects.r.ggsave(filename="variance_emp_disc_iid.png", plot=p2, width=200, height=150, unit='mm')
#%% SKEWNESS
p3=ggplot2.ggplot(data=tab_all)+\
ggplot2.geom_point(ggplot2.aes_string(x='Skewness', y='Skewness_emp') )+\
ggplot2.labs(title='Skewness',x='Discretized', y='Empirical')+\
ggplot2.geom_line(ggplot2.aes_string(x='Skewness', y='Skewness'),color='red',linetype='dashed') +\
 ggrepel.geom_text_repel(ggplot2.aes_string(x='Skewness', y='Skewness_emp',label = 'cc'))+\
ggplot2.geom_smooth(ggplot2.aes_string(x='Skewness', y='Skewness_emp') ,method='lm')
p3.plot()
#  ggplot2.geom_text(ggplot2.aes_string(label = 'cc'), hjust = 0, vjust = 0)+ \
robjects.r.ggsave(filename="skewness_emp_disc_iid.png", plot=p3, width=200, height=150, unit='mm')
#%% KURTOSIS
p4=ggplot2.ggplot(tab_all)+\
ggplot2.geom_point(ggplot2.aes_string(x='Kurtosis', y='Kurtosis_emp'))+\
ggplot2.labs(title='Kurtosis',x='Discretized', y='Empirical')+\
ggplot2.geom_line(ggplot2.aes_string(x='Kurtosis', y='Kurtosis'),color='red',linetype='dashed') +\
 ggrepel.geom_text_repel(ggplot2.aes_string(x='Kurtosis', y='Kurtosis_emp',label = 'cc'))+\
ggplot2.geom_smooth(ggplot2.aes_string(x='Kurtosis', y='Kurtosis_emp'),method='lm')
p4.plot()
robjects.r.ggsave(filename="kurtosis_emp_disc_iid.png", plot=p4, width=200, height=150, unit='mm')
# display(image_png(p))
#%%    
tab_results.to_csv('calib_results_iid.csv',index=False,header=True)
tab_results.to_latex("calib_results_iid.tex",index=False,header=True,
                     float_format="%.5f",
                     caption='Discretization diagnostics',
                     label='tab:discret_diagn',
                     column_format='c c c c c c c',
                     longtable=True)
#%% ---------------------------------------------------  END ---------------------------------------------------- #
# Show the plots (MUST BE AT THE END OTHERWISE PLOTS WILL CLOSE)
plt.show()

# %%
