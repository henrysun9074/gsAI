#!/bin/bash
#SBATCH --job-name=hpt_tuning_genomic_selection
#SBATCH --output=01_hpt.out
#SBATCH --error=01_hpt.err
#SBATCH --partition=schultzlab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

export LD_LIBRARY_PATH=/hpc/group/schultzlab/hs325/miniconda3/envs/gsAI/lib:$LD_LIBRARY_PATH
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh

cd /work/tfs3/gsAI
conda activate gsAI


## MAF 0.05
# python3 MLmodels/01_hpt.py \
#     -o "MAF0.05All" \
#     -f "MAF0.05.csv" \
#     -g "all" 

python3 MLmodels/01_hpt.py \
    -o "MAF0.05F2" \
    -f "MAF0.05.csv" \
    -g "F2" 

## MAF 0.01
# python3 MLmodels/01_hpt.py \
#     -o "MAF0.01All" \
#     -f "MAF0.01.csv" \
#     -g "all" 

python3 MLmodels/01_hpt.py \
    -o "MAF0.01F2" \
    -f "MAF0.01.csv" \
    -g "F2" 

## MAF 0.005
# python3 MLmodels/01_hpt.py \
#     -o "MAF0.005All" \
#     -f "MAF0.005.csv" \
#     -g "all" 

python3 MLmodels/01_hpt.py \
    -o "MAF0.005F2" \
    -f "MAF0.005.csv" \
    -g "F2" 
