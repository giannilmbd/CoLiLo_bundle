# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o: \
 ../src/SF_OPTIMIZE/adaptive_mix.f90 \
 ../src/SF_OPTIMIZE/brent.f90 \
 ../src/SF_OPTIMIZE/broyden1.f90 \
 ../src/SF_OPTIMIZE/broyden_mix.f90 \
 ../src/SF_OPTIMIZE/curvefit.f90 \
 ../src/SF_OPTIMIZE/fmin_Nelder_Mead.f90 \
 ../src/SF_OPTIMIZE/fmin_bfgs.f90 \
 ../src/SF_OPTIMIZE/fmin_cg.f90 \
 ../src/SF_OPTIMIZE/fmin_cg_cgplus.f90 \
 ../src/SF_OPTIMIZE/fmin_cg_minimize.f90 \
 ../src/SF_OPTIMIZE/froot_scalar.f90 \
 ../src/SF_OPTIMIZE/fsolve.f90 \
 ../src/SF_OPTIMIZE/fsolve_error.h90 \
 ../src/SF_OPTIMIZE/leastsq.f90 \
 ../src/SF_OPTIMIZE/leastsq_error.h90 \
 ../src/SF_OPTIMIZE/linear_mix.f90
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_routines.mod.stamp
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/cgfit_routines.mod.stamp
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o: include/sf_constants.mod
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o: include/sf_linalg.mod
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o.provides.build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/sf_optimize.mod.stamp
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/sf_optimize.mod.stamp: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/sf_optimize.mod src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/sf_optimize.mod.stamp Intel
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o.provides.build
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/SF_OPTIMIZE.f90.o.provides.build
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o.provides.build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_func_interface.mod.stamp
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_func_interface.mod.stamp: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/broyden_func_interface.mod src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_func_interface.mod.stamp Intel
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o.provides.build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_routines.mod.stamp
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_routines.mod.stamp: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/broyden_routines.mod src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/broyden_routines.mod.stamp Intel
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o.provides.build
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_broyden_routines.f90.o.provides.build
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_cgfit_routines.f90.o.provides.build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/cgfit_routines.mod.stamp
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/cgfit_routines.mod.stamp: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_cgfit_routines.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/cgfit_routines.mod src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/cgfit_routines.mod.stamp Intel
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_cgfit_routines.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_cgfit_routines.f90.o.provides.build
src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/build: src/SF_OPTIMIZE/CMakeFiles/SF_OPTIMIZELIB.dir/optimize_cgfit_routines.f90.o.provides.build
