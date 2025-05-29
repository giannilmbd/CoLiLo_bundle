FC=ifort

FCO =-O3 -error-limit 1 -ipo -ip -xhost -qopenmp -recursive -free -xhost -flto  -assume bscc -ipo -r8   \
 -warn interfaces -Wl,-export-dynamic -align array64byte  -qmkl -no-prec-div  -qopt-prefetch=5 -mcmodel=medium -shared-intel
FLIB= -axCORE-AVX2 -fp-model precise -L/usr/local/lib64 -lnlopt
LIBS=-L/usr/local/lib -L/usr/local/lib64 -L/usr/lib64  -L/opt/intel/oneapi/vtune/2022.4.1/lib64/    #-lmkl_intel_thread -lmkl_core -liomp5 # -I/usr/local/include 
LINK_FNL_STATIC_HPC=-z muldefs -L/usr/local/lib -Wl,-rpath,/usr/local/lib64 -lnlopt \
 -L/opt/intel/oneapi/mkl/2023.1.0/lib/intel64 \
 -L~/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/nlopt-2.7.1 \
 -L/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin \
 -Bstatic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 \
 -Bdynamic -Wl,-rpath,/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin/

# DEBUGGING
FCODB =-g -O0 -recursive -r8 -free  -debug extended -mcmodel=medium -shared-intel \
-z muldefs -L/usr/local/lib -Wl,-rpath,/usr/local/lib64 -lnlopt\
-L/opt/intel/oneapi/mkl/2023.1.0/lib/intel64 \
-L~/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/nlopt-2.7.1 -lnlopt \
-L/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin -liomp5 \
-Bdynamic -Wl,-rpath,/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin


a.out: example_hybrd1.f90   # $(FLIB) $(LIBS) toolbox.mod  lapackhelperDouble.mod ogpf.mod #lapackhelperDouble.o ogpf.o -ldl -lstdc++ -lpthread 
	$(FC) $(FCO) $(FLIB)  -o  a.out \
	../src/minpack.f90 \
	example_hybrd1.f90 $(LINK_FNL_STATIC_HPC)

clean:
	rm *.o *.out *.mod