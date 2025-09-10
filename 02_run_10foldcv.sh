#!/bin/bash
#SBATCH --job-name=10foldCV_genomic_selection
#SBATCH --output=02_cv.out
#SBATCH --error=02_cv.err
#SBATCH --partition=common 
#SBATCH --ntasks=1
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

# PYTHON_EXEC=/hpc/group/schultzlab/hs325/miniconda3/envs/gsAI/bin/python
# # sanity check
# echo "Using Python: $($PYTHON_EXEC --version)"
# $PYTHON_EXEC -m pip list | grep -E "numpy|scikit-learn|xgboost"
# $PYTHON_EXEC genomic_selection_pipeline.py

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate gsAI
# python 01_hpt.py 
python 02_10foldcv.py
