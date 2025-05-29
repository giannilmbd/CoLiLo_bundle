FC=ifort
#-prof-use  THIS CAN BE USED ONLY AFTER CREATING THE PROFILE WITH -prof-gen 
FCO = -O3 -ipo -ip -xhost -qopenmp -recursive -free -xhost -flto  -assume bscc -ipo -r8   \
 -warn interfaces -Wl,-export-dynamic -align array64byte  -qmkl -no-prec-div  -qopt-prefetch=5
 

FLIB= -axCORE-AVX2 -fp-model precise  


LIBS=-L/usr/local/lib -L/usr/lib64      #-lmkl_intel_thread -lmkl_core -liomp5 # -I/usr/local/include 
LINK_FNL_STATIC_HPC=-z muldefs -L/opt/intel/oneapi/mkl/latest/lib/intel64/ \
           -Bstatic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
		   -Wl,--end-group -liomp5 -Bdynamic -Xlinker -rpath -Xlinker /opt/intel/oneapi/mkl/latest/lib/intel64/
# DEBUGGING
FCO =-g -O0 -recursive -r8 -assume bscc  -free  -debug extended # NOTE THAT IN THE launch.json MUST INCLUDE THE MISSING SHARED LIBRARIES -- see 
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
# ../toolbox.f90 \

test_simulate.out: markov_process_mod.f90
	$(FC) $(FCO) $(FLIB) -o test_simulate.out ../my_kinds_mod.f90 \
	../mymoments.f90  \
	../ansi_colors_mod.f90 \
	../my_I_O.f90  \
	../fortran-csv-module/src/csv_kinds.f90 \
	../fortran-csv-module/src/csv_parameters.f90 \
	../fortran-csv-module/src/csv_utilities.f90 \
	../fortran-csv-module/src/csv_module.F90 \
	../read_csv_mod.f90 \
	../../pyplot-fortran/src/pyplot_module.F90 \
	../gnufor2/gnufor2.f90 \
	 ../toolbox.f90 \
	simulate_series_mod.f90 \
	markov_process_mod.f90 $(LINK_FNL_STATIC_HPC) 


clean:
	rm *.o *.mod *.out
	
