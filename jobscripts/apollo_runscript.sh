#!/bin/bash
#SBATCH --job-name=BaseMpdataDEV
#SBATCH --output logfolder/log_mesher_.%j.txt
#SBATCH --ntasks-per-node=64
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --account=Szmelter2026b

cd      $SLURM_SUBMIT_DIR

module purge
module load GCC

# case_file="half_domain_slender.vtk"
# case_file="full_domain.vtk"
case_file="slender-fine.vtk"
# case_file="star1.case"

sh "./run_edge_based_mesh_converter.sh" $case_file
