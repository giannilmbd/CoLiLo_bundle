set(CMAKE_Fortran_COMPILER "/usr/bin/f95")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "12.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/bin/gcc-ar-12")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/bin/gcc-ranlib-12")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/usr/lib/gcc/x86_64-linux-gnu/12/finclude;/opt/intel/oneapi/umf/0.10/include;/opt/intel/oneapi/tbb/2022.1/include;/opt/intel/oneapi/mpi/2021.15/include;/opt/intel/oneapi/mkl/2025.1/include;/opt/intel/oneapi/ippcp/2025.1/include;/opt/intel/oneapi/ipp/2022.1/include;/opt/intel/oneapi/dpl/2022.8/include;/opt/intel/oneapi/dpcpp-ct/2025.1/include;/opt/intel/oneapi/dnnl/2025.1/include;/opt/intel/oneapi/dev-utilities/2025.1/include;/opt/intel/oneapi/dal/2025.5/include;/opt/intel/oneapi/dal/2025.5/include/dal;/opt/intel/oneapi/ccl/2021.15/include;/opt/intel/oneapi/pti/0.12/include;/usr/lib/gcc/x86_64-linux-gnu/12/include;/usr/local/include;/usr/include/x86_64-linux-gnu;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/intel/oneapi/tcm/1.3/lib;/opt/intel/oneapi/umf/0.10/lib;/opt/intel/oneapi/tbb/2022.1/lib/intel64/lib;/opt/intel/oneapi/pti/0.12/lib;/opt/intel/oneapi/mpi/2021.15/lib;/opt/intel/oneapi/mkl/2025.1/lib;/opt/intel/oneapi/ippcp/2025.1/lib;/opt/intel/oneapi/ipp/2022.1/lib;/opt/intel/oneapi/dnnl/2025.1/lib;/opt/intel/oneapi/dal/2025.5/lib;/opt/intel/oneapi/compiler/2025.1/lib;/opt/intel/oneapi/ccl/2021.15/lib;/usr/lib/gcc/x86_64-linux-gnu/12;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib;/opt/intel/oneapi/tbb/2022.1/lib/intel64/gcc4.8")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
