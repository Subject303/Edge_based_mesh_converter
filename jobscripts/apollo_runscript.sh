#!/bin/bash
#SBATCH --job-name=BaseMpdataDEV
#SBATCH --output Logfolder/Apollotest.STD_OUT_LOG
#SBATCH --ntasks-per-node=64
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --account=Szmelter2026b

cd      $SLURM_SUBMIT_DIR

case_file="half_domain_slender.vtk"
# case_file="test_intake.vtk"
# case_file="star1.case"

sh "./run_edge_based_mesh_converter.sh" $case_file
