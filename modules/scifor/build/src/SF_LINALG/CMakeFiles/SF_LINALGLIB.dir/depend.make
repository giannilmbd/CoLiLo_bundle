# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o: \
 ../src/SF_LINALG/linalg_auxiliary.f90 \
 ../src/SF_LINALG/linalg_blacs_aux.f90 \
 ../src/SF_LINALG/linalg_blas.f90 \
 ../src/SF_LINALG/linalg_build_tridiag.f90 \
 ../src/SF_LINALG/linalg_check_tridiag.f90 \
 ../src/SF_LINALG/linalg_eig.f90 \
 ../src/SF_LINALG/linalg_eigh.f90 \
 ../src/SF_LINALG/linalg_eigh_jacobi.f90 \
 ../src/SF_LINALG/linalg_eigvals.f90 \
 ../src/SF_LINALG/linalg_eigvalsh.f90 \
 ../src/SF_LINALG/linalg_external_products.f90 \
 ../src/SF_LINALG/linalg_get_tridiag.f90 \
 ../src/SF_LINALG/linalg_inv.f90 \
 ../src/SF_LINALG/linalg_inv_gj.f90 \
 ../src/SF_LINALG/linalg_inv_her.f90 \
 ../src/SF_LINALG/linalg_inv_sym.f90 \
 ../src/SF_LINALG/linalg_inv_triang.f90 \
 ../src/SF_LINALG/linalg_inv_tridiag.f90 \
 ../src/SF_LINALG/linalg_lstsq.f90 \
 ../src/SF_LINALG/linalg_p_blas.f90 \
 ../src/SF_LINALG/linalg_p_eigh.f90 \
 ../src/SF_LINALG/linalg_p_inv.f90 \
 ../src/SF_LINALG/linalg_solve.f90 \
 ../src/SF_LINALG/linalg_svd.f90 \
 ../src/SF_LINALG/linalg_svdvals.f90
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o: include/sf_blacs.mod
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o: include/sf_mpi.mod
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o.provides.build: src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/sf_linalg.mod.stamp
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/sf_linalg.mod.stamp: src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/sf_linalg.mod src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/sf_linalg.mod.stamp Intel
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o.provides.build
src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/build: src/SF_LINALG/CMakeFiles/SF_LINALGLIB.dir/SF_LINALG.f90.o.provides.build
