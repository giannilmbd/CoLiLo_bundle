#!/bin/bash

echo "extension 1 for wrong stride; extension 2 for right stride"
if [ -z "$1" ]; then
    echo "No argument supplied, using default value 1"
    export ANNA_EXTENSION=1
else
    export ANNA_EXTENSION=$1
fi
echo "Using ANNA_EXTENSION=$ANNA_EXTENSION"

cd ~/Dropbox/projects/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/models/CoLiLoCodes4Circulation/CoLiLo_CRRA_2c_step3_AR1_EZ_3states/test_autarky

if [ "$ANNA_EXTENSION" == "1" ]; then
    cp src/value_autarky_wrong_stride.f90 src/value_autarky.f90
else
    cp src/value_autarky_right_stride.f90 src/value_autarky.f90
fi

if [ $? -ne 0 ]; then
    echo "Error: Copying value_autarky file failed."
    exit 1
fi

# Ensure build directory exists and is clean
if [ ! -d "build" ]; then
    echo "Build directory not found, creating it..."
    mkdir build
fi

cd build || { echo "Failed to cd to build directory"; exit 1; }
rm -rf *

# Run the build commands
cmake ..
cmake --build .
./my_program v 2






