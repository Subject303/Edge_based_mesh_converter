#!/bin/bash



COMPLIST=" ../src/src_shared/* "

COMPLIST=$COMPLIST" ../src/src_mpi/*"

LIB_LIST="./libs/lib "

INCLUDE_LIST="./libs/include "

#COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -g3 -fcheck=all -fbacktrace -ffpe-trap=zero,overflow  "
COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "

module purge
module load GCC

mkdir -p  "./modfiles"

mpifort $COMP_OPTIONS -o "../converter" $COMPLIST -I$INCLUDE_LIST -L$LIB_LIST -J"./modfiles" -lfmetis -lmetis

echo ' passed compilation '

echo ' running python preprocessor '

python3 "../python_scripts/case_preprocessor.py"

echo ' running main converter program parralel with mpi'

mpirun -np $PBS_NP ../converter

echo ' finished. '
