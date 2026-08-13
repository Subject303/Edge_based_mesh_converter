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

case_file="half_domain_slender.vtk"

sh "./run_edge_based_mesh_converter.sh" $case_file
