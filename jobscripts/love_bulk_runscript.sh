#!/bin/bash
#PBS -N Mesh_converter
#PBS -e logfolder/log.$PBS_JOBID.txt
#PBS -j eo
#PBS -m bf
#PBS -M j.thomas@lboro.ac.uk
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=40
#PBS -A Szmelter2025b

cd      $PBS_O_WORKDIR

module purge
module load GCC
module load Python/3.10.8-GCCcore-12.2.0

sh "./run_edge_based_mesh_converter.sh" "2_rad_sphere_poly.case"

sh "./run_edge_based_mesh_converter.sh" "2.case"

sh "./run_edge_based_mesh_converter.sh" "tet_cube.case"

wait
exit
