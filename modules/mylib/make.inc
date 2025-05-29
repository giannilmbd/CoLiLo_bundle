#
SHELL = /bin/sh
#-funroll-loops -fopenmp -funroll-all-loops (not necessarily faster)-funswitch-loops 

 FC='gfortran-10'
# FOPTIMFLAGS='-O3  -ftree-vectorize -funroll-loops -funroll-all-loops -ffast-math   -ftree-vectorizer-verbose=4 -xf95-cpp-input  ' #gfortran

 FOPTIMFLAGS='-O3 -finline-functions -xf95-cpp-input -fPIC -lmwlapack -lmwblas' #gfortran-lmwlapack -lmwlapack -lmwblas
 LD='gfortran-10'
 INCLUDE='-I"/opt/matlab14/R2014a/extern/include/" -I"/opt/matlab14/R2014a/extern/examples/refbook/" -I"/opt/matlab14/R2014a/extern/include" '



#FC='ifort'
#FOPTIMFLAGS='-O3 -axP -unroll -funroll-loops -fpp -ipo -no-prec-div -static -fpp -Wl,-static-intel' #ifort
#LD='ifort'

# FC='f95'
# FOPTIMFLAGS='-O3 -finline-functions -ftree-vectorize -funroll-loops -funroll-all-loops -ffast-math -xf95-cpp-input ' #gfortran
# LD='f95'
