#!/bin/bash
#PBS -N Mesh_converter
#PBS -e logfolder/log.$PBS_JOBID.txt
#PBS -j eo
#PBS -m bf
#PBS -M j.thomas@lboro.ac.uk
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=40
#PBS -A Szmelter2025b

cd      ..

main_path=$(pwd)
MOD=' '$(ls $main_path'/modfiles')
MOD=${MOD//"dum"/}
cd "./modfiles"
rm -f $MOD
cd      ..
rm -f "./converter"

module purge
module load GCC
module load Python/3.10.8-GCCcore-12.2.0

case_file=$1

COMPLIST=" ./src/src_shared/global_data.f90 ./src/src_shared/quicksort.f90 ./src/src_shared/utils.f90  ./src/src_shared/boundary_preprocessing_module.f90 ./src/src_shared/volume_processing.f90 ./src/src_shared/preprocessor_routine_module.f90 ./src/src_shared/setup_configuration_module.f90 ./src/src_shared/manual_geom_generation.f90 ./src/src_shared/outputting_routines.f90 "

COMPLIST=$COMPLIST" ./src/src_serial/read_data_serial.f90 ./src/src_serial/data_outputting_serial.f90 ./src/src_serial/data_preprocessing_serial.f90 ./src/src_serial/main_prog_serial.f90 "

COMP_OPTIONS=" -o2 -fdec-format-defaults -ffree-line-length-none -g3 -fcheck=all -fbacktrace -ffpe-trap=zero,overflow  "
# COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "

echo ' beginning compilation '

gfortran $COMP_OPTIONS -o "./converter" $COMPLIST -J"./modfiles"

if [ -f "./converter"  ]; then
    echo " passed compilation "
else
    echo "compilation failed"
    exit
fi

echo ' running python preprocessor '

python3 "./python_scripts/case_preprocessor.py" -f $case_file

echo ' running main converter program serial '

# ./converter $case_file

echo ' finished. '

cd "./modfiles"
rm -f $MOD
cd      ..
rm -f "./converter"
