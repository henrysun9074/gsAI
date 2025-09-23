#!/bin/bash
#SBATCH --job-name=genomic-selection
#SBATCH --output=paul.out
#SBATCH --error=paul.err
#SBATCH --partition=common
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

cd /work/tfs3/gsAI/4paul
module load R/4.4.0
export R_LIBS=/work/tfs3/gsAI/4paul/rlib

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Starting R script at: $(date)"

#rerunning noQC test cases
Rscript DarpaCV_NoQC.R

echo "R script finished at: $(date)"

