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

## for -i searches within /MLmodels/models, -o saves within /MLmodels/gebvs, -f looks within /data

python3 MLmodels/02_crossval.py \
    -i "ModelDirName" \
    -o "OutDirName" \
    -f "TBA.csv" \
    -g "all" 

python3 MLmodels/02_crossval.py \
    -i "ModelDirName" \
    -o "OutDirName" \
    -f "TBA.csv" \
    -g "F2" 