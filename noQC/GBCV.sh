#!/bin/bash
#SBATCH --job-name=GB_10foldCV
#SBATCH --output=GBCV.out
#SBATCH --error=GBCV.err
#SBATCH --partition=common
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

cd /work/tfs3/gsAI
export LD_LIBRARY_PATH=/hpc/group/schultzlab/hs325/miniconda3/envs/gsAI/lib:$LD_LIBRARY_PATH

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"

# conda init 
conda activate gsAI
python GBCV.py

echo "Python script finished at: $(date)"