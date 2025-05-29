#!/bin/bash


set -e
current_dir=$(pwd)
# Step 0: Clean build directory
rm -rf build
mkdir build
rm -r ../nlopt/build 
# Step 1: Choose compiler
echo "Can provide compiler choice in line: 1=Intel, 2=gfortran"
if [ -z "$1" ]; then
    read -p "Enter your choice (1 or 2): " compiler_choice
else
    compiler_choice=$1
fi

if [ "$compiler_choice" -eq 1 ]; then
    echo "Using Intel Compiler"
    cp CMakeLists_intel.txt CMakeLists.txt
elif [ "$compiler_choice" -eq 2 ]; then
    echo "Using gfortran Compiler"
    cp CMakeLists_gfortran.txt CMakeLists.txt
else
    echo "Invalid choice. Exiting."
    exit 1
fi

# Step 2: Check for NLOpt static library and Fortran interface
NLOPT_TARGET_DIR="../nlopt/build/install"
if [[ -f "$NLOPT_TARGET_DIR/libnlopt.a" && -f "$NLOPT_TARGET_DIR/nlopt.f" ]]; then
    echo "NLOpt static library and interface already present, skipping build."
else
    echo "Building NLOpt statically..."
    cd ../nlopt
    mkdir -p build && cd build

    cmake \
    -DPython3_EXECUTABLE=$(which python3) \
    -DNLOPT_FORTRAN=ON \
    -DCMAKE_INSTALL_PREFIX=$PWD/install \
    -DBUILD_SHARED_LIBS=OFF \
  ..
#   -DCMAKE_C_FLAGS="-Ofast -march=native" \
#   -DCMAKE_CXX_FLAGS="-Ofast -march=native" \
#   -DCMAKE_Fortran_FLAGS="-Ofast -march=native" \
#   ..

    cmake --build . --parallel $(nproc)
    cmake --install .


    cd $current_dir
fi

# Step 3: Configure and build your project
cd build
cmake -G "Unix Makefiles" ..
cmake --build . --parallel $(nproc)

# Step 4: Run the compiled program
echo "Done Compilation... now running"
./my_program v 62
