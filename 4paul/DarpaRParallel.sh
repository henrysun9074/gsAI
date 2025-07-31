#!/bin/bash
#SBATCH --job-name=genomic-selection-parallel
#SBATCH --output=paul_parallel.out
#SBATCH --error=paul_parallel.err
#SBATCH --partition=common
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

cd /work/tfs3/gsAI/4paul
module load R/4.4.0
export R_LIBS=/work/tfs3/gsAI/4paul/rlib
export LD_LIBRARY_PATH=/work/tfs3/gsAI/4paul/rlib/BGLR:$LD_LIBRARY_PATH

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Starting R script at: $(date)"

echo "Checking LAPACK linkage:"
ldd /work/tfs3/gsAI/4paul/rlib/BGLR/libs/BGLR.so | grep lapack

Rscript DarpaAllParallel.R
# Rscript testR.R

echo "R script finished at: $(date)"
