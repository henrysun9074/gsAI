#!/bin/bash
#SBATCH --job-name=hpt_tuning_genomic_selection
#SBATCH --output=01_hpt.out
#SBATCH --error=01_hpt.err
#SBATCH --partition=gpu-common 
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

cd /work/tfs3/gsAI
export LD_LIBRARY_PATH=/hpc/group/schultzlab/hs325/miniconda3/envs/gsAI/lib:$LD_LIBRARY_PATH
module load CUDA

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate gsAI

python3 01_hpt.py \
    -o "Oct27_MAF01_F2" \
    -f "MAF01_DarpaQCFiltered.csv" \
    -g "F2" 

python3 01_hpt.py \
    -o "Oct27_MAF02_F2" \
    -f "MAF02_DarpaQCFiltered.csv" \
    -g "F2" 

