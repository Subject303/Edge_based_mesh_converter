#!/bin/bash



COMPLIST=" ../src/src_shared/* "

COMPLIST=$COMPLIST" ../src/src_serial/*"

COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "

module purge
module load GCC

gfortran $COMP_OPTIONS -o "../converter" $COMPLIST -J"./modfiles"

echo ' passed compilation '

echo ' running python preprocessor '

python3 "../python_scripts/case_preprocessor.py"

echo ' running main converter program serial '

../converter

echo ' finished. '
