# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build
src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/SF_DERIVATE.f90.o: \
 ../src/SF_DERIVATE/derivate_fjacobian_c.f90 \
 ../src/SF_DERIVATE/derivate_fjacobian_d.f90
src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/SF_DERIVATE.f90.o.provides.build: src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/sf_derivate.mod.stamp
src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/sf_derivate.mod.stamp: src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/SF_DERIVATE.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/sf_derivate.mod src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/sf_derivate.mod.stamp Intel
src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/SF_DERIVATE.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/SF_DERIVATE.f90.o.provides.build
src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/build: src/SF_DERIVATE/CMakeFiles/SF_DERIVATELIB.dir/SF_DERIVATE.f90.o.provides.build
