# Makefile for the interpolation test

# Compiler
FC=ifort #gfortran # -heap-arrays 5000000
FF=-O3 
# Compiler flags
FLIB= -recursive -assume bscc  -xHost -fp-model source -ftz -no-prec-div -no-prec-sqrt -ip -ipo   -extend-source  -free -diag-enable=warn -qopenmp -fPIC -warn all,noexternal,nounused -warn interface -check nooutput_conversion -traceback
LIBS= -L/urs/lib64 -llapack -lblas
INCLUDES= -I/usr/local 

all: test.out

test.out: test.f90 mypolynte_mod.o globals_4test_mod.o linespaced_mod.o 
	$(FC) $(FF) $(FLIB) \
	../OGPF/ogpf.f90 \
	../../chebyshev/chebyshev_polynomial.f90 \
	../mylib/solve_lapack.f90 \
	test.f90 mypolynte_mod.o globals_4test_mod.o linespaced_mod.o $(LIBS) -o test.out

mypolynte_mod.o: mypolynte_mod.f90 globals_4test_mod.o
	$(FC) $(FF) $(FLIB) $(INCLUDES) -c mypolynte_mod.f90

globals_4test_mod.o: globals_4test_mod.f90
	$(FC) $(FF) $(FLIB) $(INCLUDES) -c globals_4test_mod.f90

linespaced_mod.o: ../linespaced_mod.f90
	$(FC) $(FF) $(FLIB) $(INCLUDES) -c ../my_kinds_mod.f90 ../linespaced_mod.f90

clean:
	rm -f *.o *.mod test.out
