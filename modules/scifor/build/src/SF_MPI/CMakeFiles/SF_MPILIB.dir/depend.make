# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_BLACS.f90.o.provides.build: src/SF_MPI/CMakeFiles/SF_MPILIB.dir/sf_blacs.mod.stamp
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/sf_blacs.mod.stamp: src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_BLACS.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/sf_blacs.mod src/SF_MPI/CMakeFiles/SF_MPILIB.dir/sf_blacs.mod.stamp Intel
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_BLACS.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_BLACS.f90.o.provides.build
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/build: src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_BLACS.f90.o.provides.build
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_MPI.f90.o.provides.build: src/SF_MPI/CMakeFiles/SF_MPILIB.dir/sf_mpi.mod.stamp
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/sf_mpi.mod.stamp: src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_MPI.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/sf_mpi.mod src/SF_MPI/CMakeFiles/SF_MPILIB.dir/sf_mpi.mod.stamp Intel
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_MPI.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_MPI.f90.o.provides.build
src/SF_MPI/CMakeFiles/SF_MPILIB.dir/build: src/SF_MPI/CMakeFiles/SF_MPILIB.dir/SF_MPI.f90.o.provides.build
