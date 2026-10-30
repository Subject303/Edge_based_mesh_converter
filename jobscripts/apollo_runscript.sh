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

# case_file="404040.vtk" 
case_file="40-semi.case" 

sh "./run_edge_based_mesh_converter.sh" $case_file

# case_file="808080.vtk" 
# case_file="512k-tet.case" 
# 
# sh "./run_edge_based_mesh_converter.sh" $case_file

# case_file="160160160.vtk" 
# case_file="4000k-tet.case" 
# 
# sh "./run_edge_based_mesh_converter.sh" $case_file

# case_file="320320320.vtk" 
# 
# sh "./run_edge_based_mesh_converter.sh" $case_file
