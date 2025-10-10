#!/bin/bash
#SBATCH --job-name=10foldCV_genomic_selection
#SBATCH --output=02_cv.out
#SBATCH --error=02_cv.err
#SBATCH --partition=schultzlab 
#SBATCH --time=4-00:00:00
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
source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate gsAI

python3 02_10foldcv.py \
    -i "oct04" \
    -o "Oct09_allnewQC" \
    -f "39kDarpaQCFiltered.csv" \
    -g "all" 

python3 02_10foldcv.py \
    -i "sep22" \
    -o "Oct09_F0newQC" \
    -f "39kDarpaQCFiltered.csv" \
    -g "F0" 

python3 02_10foldcv.py \
    -i "sep22_F2QC" \
    -o "Oct09_F2newQC" \
    -f "39kDarpaQCFiltered.csv" \
    -g "F2" 

