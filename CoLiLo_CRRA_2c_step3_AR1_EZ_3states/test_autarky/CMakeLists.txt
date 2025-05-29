cmake_minimum_required(VERSION 3.10)
project(MyFortranProject Fortran)

message("
##############################
For Release build (default):
    cmake ..
    make

For Debugging build:
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    make
##############################
")

# Use gfortran instead of ifx
set(CMAKE_Fortran_COMPILER gfortran)

# Set Fortran flags (gfortran-compatible)
# Set Fortran flags (gfortran-compatible)
set(CMAKE_Fortran_FLAGS_RELEASE
  "-O3 -fdefault-real-8 -fopenmp -ffree-form -flto -ffp-contract=off -fno-associative-math -fno-reciprocal-math -ffree-line-length-none -march=native"
)

set(CMAKE_Fortran_FLAGS_DEBUG
  "-g -O0 -fcheck=all  -fdefault-real-8 -fopenmp -ffree-form -fno-fast-math -ffp-contract=off -fno-associative-math -fno-reciprocal-math -flto -ffree-line-length-none \
    -march=native "
)

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")

# Set linker flags
set(LINK_FLAGS "-fopenmp -llapack -lblas -lm")

# Source files
set(SOURCE_FILES
    ../src/my_kinds_mod.f90
    ../src/my_toolbox.f90
    ../src/discretize.f90
    ../src/mymoments.f90
    ../src/my_I_O.f90
    ../src/csv_kinds.f90
    ../src/csv_parameters.f90
    ../src/csv_utilities.f90
    ../src/csv_module.F90
    ../src/read_csv_mod.f90
    ../src/write_csv_portable_module.f90
    ../src/linespaced_mod.f90
    src/globals.f90
    ../src/MarkovMoments_mod.f90
    ../src/pcu_CRRA_module.f90
    ../src/consumption1_autarky_module.f90
    ../src/optimization_interface_mod.f90
    ../src/minpack.f90
    ../src/initialize_and_save_mod_CRRA_2c.f90
    ../src/F_kappa.f90
    ../src/core_model.f90
    ../../modules/bspline-fortran/src/bspline_blas_module.F90
    ../../modules/bspline-fortran/src/bspline_kinds_module.F90
    ../../modules/bspline-fortran/src/bspline_oo_module.f90
    ../../modules/bspline-fortran/src/bspline_defc_module.F90
    ../../modules/bspline-fortran/src/bspline_module.f90
    ../../modules/bspline-fortran/src/bspline_sub_module.f90
    ../src/interp_mod.f90
    ../src/myfunc_update.f90
    ../src/solve_policy_function.f90
    ../src/solve_arrow_securities.f90
    src/core_model_autarky.f90
    src/bisection_price_autarky.f90
    src/value_autarky.f90
    src/value_autarky.f90
    ../src/welfare.f90
    ../src/find_kappa_mod2.f90
    src/main.f90
)

foreach(FILE ${SOURCE_FILES})
    if(FILE MATCHES "\\.F90$")
        set_source_files_properties(${FILE} PROPERTIES LANGUAGE Fortran)
    endif()
endforeach()

add_executable(my_program ${SOURCE_FILES})

# Include and link external libraries
# Correct NLOpt include + link paths
include_directories(${CMAKE_SOURCE_DIR}/../../nlopt/build/install/include)
target_link_libraries(my_program PRIVATE ${CMAKE_SOURCE_DIR}/../../nlopt/build/install/lib/libnlopt.a stdc++)

# LAPACK and BLAS
find_package(LAPACK REQUIRED)
target_link_libraries(my_program PRIVATE ${LAPACK_LIBRARIES})

# Optional diagnostics
message(STATUS "")
message(STATUS "=== Checking Required Libraries ===")

# LAPACK
message(STATUS "LAPACK        : ${LAPACK_FOUND}")

# BLAS
find_package(BLAS REQUIRED)
message(STATUS "BLAS          : ${BLAS_FOUND}")

# OpenMP
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP        : ${OpenMP_Fortran_FOUND}")
if(OpenMP_Fortran_FOUND)
  target_link_libraries(my_program PRIVATE OpenMP::OpenMP_Fortran)
endif()

# stdc++ check (often pulled in via nlopt)
find_library(STDCPP_LIB stdc++
  HINTS
    /usr/lib/gcc/x86_64-linux-gnu/13
    /usr/lib/gcc/x86_64-linux-gnu/12
    /usr/lib/x86_64-linux-gnu
    /usr/lib 
)
if(STDCPP_LIB)
  message(STATUS "libstdc++     : ${STDCPP_LIB}")
  target_link_libraries(my_program PRIVATE ${STDCPP_LIB})
else()
  message(WARNING "libstdc++ not found — required if using NLOpt with C++ components -- might not be needed")
endif()

# quadmath for gfortran
find_library(QUADMATH_LIB quadmath
  HINTS
    /usr/lib/gcc/x86_64-linux-gnu/13
    /usr/lib/gcc/x86_64-linux-gnu/12
    /usr/lib/x86_64-linux-gnu
    /usr/lib 
)
if(QUADMATH_LIB)
  message(STATUS "libquadmath   : ${QUADMATH_LIB}")
  target_link_libraries(my_program PRIVATE ${QUADMATH_LIB})
else()
  message(WARNING "libquadmath not found — may be required by gfortran")
endif()

message(STATUS "====================================")

add_custom_target(run
    COMMAND ulimit -s unlimited && ./my_program
)

add_custom_target(cleanup
    COMMAND unset OMP_NUM_THREADS && rm -f *.o *.mod *.out && ulimit -s unlimited
)
