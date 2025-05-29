# Structure of Python codes

* Main code to calibrate moments is Run_mixture.py
* This code uses moments estimated from data with quick_stylized_facts.r from ../Rcodes
* It uses match_moments.py: this finds the parameters of the Gaussian Mixture that minimize the distance between its first fourth moments and the empirical counterparts. compute_moments.py is an auxiliary function used here and in the next. 
* It uses dicreteGMAR.py and dicreteApproximation.py. These discretize the AR(1) process
* Run_mixture.py saves the probability matrices Pmatrices.csv and Xmatrices.csv with the markow process for each country
* It also saves some figures and calibrario results both as csv and latex table
