#!/bin/bash
#SBATCH --job-name=genomic_selection
#SBATCH --output=genomic_selection.out
#SBATCH --error=genomic_selection.err
#SBATCH --partition=gpu-common
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

cd /work/tfs3/gsAI
export LD_LIBRARY_PATH=/hpc/group/schultzlab/hs325/miniconda3/envs/gsAI/lib:$LD_LIBRARY_PATH

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"

module load CUDA
conda activate gsAI

python genomic_selection_pipeline.py
