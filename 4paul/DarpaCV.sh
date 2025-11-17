#!/bin/bash
#SBATCH --job-name=genomic-selection
#SBATCH --output=Rmodels.out
#SBATCH --error=Rmodels.err
#SBATCH --partition=schultzlab 
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
export LD_LIBRARY_PATH=/hpc/group/schultzlab/hs325/miniconda3/envs/gsAI/lib:$LD_LIBRARY_PATH

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Starting R script at: $(date)"

# Rscript DarpaCV.R \
#     "phenTrainDarpa.xlsx" \
#     "/work/tfs3/gsAI/39kDarpaQCFiltered.csv" \
#     "allMAF05QC" \
#     0.2

Rscript DarpaCV_filtergen.R \
    "phenTrainDarpa.xlsx" \
    "/work/tfs3/gsAI/39kDarpaQCFiltered.csv" \
    "F2" \
    "F2MAF05QC" \
    0.2

echo "R script finished at: $(date)"

