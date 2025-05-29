# include make.inc
# source /opt/intel/oneapi/setvars.sh
# source ~/Dropbox/fortran/IMSL/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/bin/fnlsetup.sh
# source ~/Dropbox/fortran/IMSL/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/bin/rdh800in212x64.sh
# source /opt/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/bin/fnlsetup.sh
# source /opt/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/bin/rdh800in212x64.sh
# FC=mpiifort # -heap-arrays 5000000: heap-arrays moves arrays from stack to heap in principle helping with large arrays: Any arrays known at compile-time to be larger than size are allocated on the heap instead of the stack. (never worked)
  FC=ifort
#FC=ifx 
# NOTE for DEBUGGING: Not all FLIB options are compatible with debugging; ip and ipo don't work (apparently only if at the end)
#######################################################################################
#############        NB  #########################3
# FOR PRECISION USE -fp-model keyword; extended should not be supported for fortran but appears to be the best (replaced by fast)
#-fp-model precise -fp-model except (or, more drastically, -fp-model strict) -axAVX,CORE-AVX2 
# The recommendation for obtaining floating-point values that are compliant with
#ANSI / IEEE standards for C++ and Fortran is:
#/fp:precise /fp:source (Windows)
#-fp-model precise –fp-model source (Linux or Mac OS X) -real_size 32
#######################################################################3
#### NB to get -fast working with IMSL (since it calls -static)  must use  $(LINK_FNL_STATIC_HPC)
# -no-prec-div  -fp-model fast=2  

#NB: -fast calls static libraries (-static) see below how to subset dynamic libraries or yum install glibc-static to get the static versions



debug=debug
FCO = -O3 -fast -qopenmp -recursive -free -xhost -flto  -assume bscc -ipo -r8  -Wl,-export-dynamic -align array64byte  -qmkl -assume nan-compares \
 -no-prec-div  -qopt-prefetch=5#  -qopt-report#-assume nan-compares -no-prec-div -fopenmp-target-do-concurrent -qopt-prefetch=5#-qmkl=parallel  -Ofast   # fast not working? -fopenmp-targets=spir64  THIS LOOKS FOR LIBOMPTARGET BUT CANNOT FIND IT
FLIB= $(FLIBA) -axCORE-AVX2 -fp-model precise  #-fopenmp-targets=nvptx64-nvidia-cuda #-fopenmp-targets=spir64 #The latter crashes as GPU set-up not correct
LIBS= -L/usr/lib64 -lm -lc -L/opt/intel/oneapi/mkl/latest/lib/intel64
    #-lmkl_intel_thread -lmkl_core -liomp5 # -I/usr/local/include 
LINK_FNL_STATIC_HPC=$(LIBS) -z muldefs -L/opt/intel/oneapi/mkl/latest/lib/intel64/ \
           -Bstatic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group \
		   -liomp5 -Bdynamic -Xlinker -rpath -Xlinker /opt/intel/oneapi/mkl/latest/lib/intel64

DIRECTORY=direct
FCO =-fp-stack-check -check bounds -g -O0 -recursive -r8 -free -traceback -debug extended # NOTE THAT IN THE launch.json MUST INCLUDE THE MISSING SHARED LIBRARIES -- see 
# LIBS=""
# # if using gfortran for debugging
# FC=gfortran
# FLIB=-frecursive -fopenmp -ffree-line-length-none -fcheck=all -Wall
# FCO=-g -O0 #-O3 # -g -O0 #-O3 
# LIBS=
# # end using fortran
# NOTE ON LIBRARIES: MAKE SURE LINK WITH /usr/lib OTHERWISE NO DEBUGGING  -L/opt/intel/oneapi/compiler/2022.0.2/linux/lib/ -lomptarget 
#-llapack -lblas
#LIBS=  -I/usr/local/include -I~/Dropbox/fortran/IMSL/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/include $(LINK_FNL_SHARED_IMSL) -L/mnt/OneTera/Dropbox/fortran/IMSL/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/lib/ -limsl#-L/usr/local/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/lapack/ #-lnlopt -I/usr/include/mkl #/usr/local/lib/libfortran_stdlib.a -I/usr/local/include/fortran_stdlib/GNU-10.3.0

#$(LINK_FNL_SHARED_IMSL) -L/mnt/OneTera/Dropbox/fortran/IMSL/RogueWave/imsl/fnl-2021.0.0/rdh800in212x64/lib/ -limsl

save1t.out:  # $(FLIB) $(F90FLAGS) $(LIBS)$(LINK_FNL_STATIC_HPC) toolbox.mod  lapackhelperDouble.mod ogpf.mod #lapackhelperDouble.o ogpf.o -ldl -lstdc++ -lpthread 
	$(FC)  $(FLIB) $(FCO) \
	../toolbox.f90 \
	../my_kinds_mod.f90 \
	../discretize.f90 \
	../mymoments.f90  \
	../my_I_O.f90  \
	../fortran-csv-module/src/csv_kinds.f90 \
	../fortran-csv-module/src/csv_parameters.f90 \
	 ../fortran-csv-module/src/csv_utilities.f90 \
	../fortran-csv-module/src/csv_module.F90 \
	../read_csv_mod.f90 \
	../linespaced_mod.f90 \
	$(direct)/globals.f90 \
        $(direct)/core_model.f90 \
	$(direct)/price.f90 \
	$(direct)/myfunc_update.f90 \
	save_coeffs_one_time_EZ_k.f90 $(LINK_FNL_STATIC_HPC) -o  save1t.out

# test.o: test.f90 lapackhelperDouble.mod ogpf.mod toolbox.mod
# 	$(FC) -c -I/usr/local/include/ "test.f90" -o "test.o" $(FLIB) $(LIBS) 


# toolbox.mod: ~/Dropbox/projects/GlobalSolutions/FortranBook/cefortran/installation/toolbox/toolbox.f90
# 	$(FC) -c ~/Dropbox/projects/GlobalSolutions/FortranBook/cefortran/installation/toolbox/toolbox.f90 $(FLIB) -o toolbox.o 




# lapackhelperDouble.mod: ~/Dropbox/fortran/LAPACK_helper/lapackhelperDouble.f90
# 	$(FC) -c ~/Dropbox/fortran/LAPACK_helper/lapackhelperDouble.f90  $(FLIB)

# ogpf.mod: ~/Dropbox/fortran/ogpf-master/src/ogpf.f90
# 	$(FC) -c ~/Dropbox/fortran/ogpf-master/src/ogpf.f90 $(FLIB)	

clean:
	unset OMP_NUM_THREADS
	rm save1t.out *.mod *.o
	

