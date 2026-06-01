#!/bin/bash
clear
set -e

MOD=$(ls ../modfiles/*.mod)
rm -f $MOD
rm -f "../converter"

LOG="./logfile.txt"
rm -f $LOG
touch $LOG
chmod 777 $LOG


COMPLIST=" ../src/src_shared/global_data.f90 ../src/src_shared/quicksort.f90 ../src/src_shared/utils.f90  ../src/src_shared/boundary_preprocessing_module.f90 ../src/src_shared/external_face_processing.f90 ../src/src_shared/internal_edge_processing.f90 ../src/src_shared/external_edge_processing.f90 ../src/src_shared/preprocessor_routine_module.f90 ../src/src_shared/setup_configuration_module.f90 "

COMPLIST=$COMPLIST" ../src/src_serial/read_data_serial.f90 ../src/src_serial/data_outputting_serial.f90 ../src/src_serial/data_preprocessing_serial.f90 ../src/src_serial/main_prog_serial.f90 "

COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -g3 -fcheck=all -fbacktrace -ffpe-trap=zero,overflow  "
# COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "

echo ' beginning compilation ' 2>&1 | tee $LOG

gfortran $COMP_OPTIONS -o "../converter" $COMPLIST -J"../modfiles" 2>&1 | tee -a $LOG

if [ -f "../converter"  ]; then
    echo " passed compilation " 2>&1 | tee -a $LOG
else
    echo "compilation failed" 2>&1 | tee -a $LOG
    exit
fi

echo ' running python preprocessor ' 2>&1 | tee -a $LOG

# python3 "../python_scripts/case_preprocessor.py" 2>&1 | tee -a $LOG

echo ' running main converter program serial ' 2>&1 | tee -a $LOG

../converter 2>&1 | tee -a $LOG

echo ' finished. ' 2>&1 | tee -a $LOG
