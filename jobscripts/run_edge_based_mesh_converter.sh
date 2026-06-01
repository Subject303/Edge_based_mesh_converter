#!/bin/bash



COMPLIST=" ../src/src_shared/* "
echo $complist

# COMPLIST=$COMPLIST" ../src/src_serial/*"
# 
# LIB_LIST="./libs/lib "
# 
# INCLUDE_LIST="./libs/include "
# 
# #COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -g3 -fcheck=all -fbacktrace -ffpe-trap=zero,overflow  "
# COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "
# 
# module purge
# module load METIS/5.1.0-GCCcore-12.2.0 
# module load GCC
# 
# #IMPORT CONFIGS
# set -a
# source "./config/Testconfig.env"
# set +a
# 
# mkdir -p  "./modfiles"
# 
# mpifort $COMP_OPTIONS -o ./Flowsolver $COMPLIST -I$INCLUDE_LIST -L$LIB_LIST -J"./modfiles" -lfmetis -lmetis
# 
# echo ' passed compilation '
# 
# mpirun -np $PBS_NP ./Flowsolver
