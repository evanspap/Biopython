#!/bin/bash

#SBATCH --nodes=1
#SBATCH -p short-40core
#SBATCH --mem=175G
#SBATCH --time 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40   # Adjust based on the number of cores available
#SBATCH --job-name=parallel_job

# Load the module containing `parallel` if it's not already in your path
#module load parallel
module load gnu-parallel/6.0


SLURM_CPUS_PER_TASK=40

# Run the hollow command in parallel
ls *.pdb | parallel -j $SLURM_CPUS_PER_TASK  hollow -o hollow{.}.pdb {}

#$WORKDIR/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r {} -A checkhydrogens -o {.}.pdbqt