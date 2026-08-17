#!/bin/bash
#SBATCH --job-name=BaseMpdataDEV
#SBATCH --output Logfolder/mesh_log_job_$SLURM_JOB_ID.txt
#SBATCH --ntasks-per-node=64
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --account=Szmelter2026b

cd      $SLURM_SUBMIT_DIR

# pip list

module purge
# module load VTK/9.3.0-foss-2023a

pip list

# module load GCC
# module load Python/3.13.1-GCCcore-14.2.0

case_file="half_domain_slender.vtk"
# case_file="test_intake.vtk"
# case_file="star1.case"

sh "./run_edge_based_mesh_converter.sh" $case_file
