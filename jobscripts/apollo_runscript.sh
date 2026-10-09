#!/bin/bash
#SBATCH --job-name=edge-mesher
#SBATCH --output logfolder/log_mesher_.%j.txt
#SBATCH --ntasks-per-node=64
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --account=Szmelter2026b

cd      $SLURM_SUBMIT_DIR

module purge
module load GCC

# case_file="scaled-fine.vtk" 
case_file="t-birch-slender-visc-2nd2nd_1@45300.case" 

sh "./run_edge_based_mesh_converter.sh" $case_file
