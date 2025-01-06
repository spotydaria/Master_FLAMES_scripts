#!/bin/bash
#SBATCH --job-name=pops_split       # Job name
#SBATCH --time=05:00:00              # Time limit (3 hours)
#SBATCH --mem=32G

# Activate your Conda environment
# source /path/to/conda.sh  # This line loads Conda commands
source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate POPs

export POPS_TOOL="/home/dkramarenk/projects/tools/POPs"

python ${POPS_TOOL}/pops/munge_feature_directory.py \
  --gene_annot_path ${POPS_TOOL}/pops/example/data/utils/gene_annot_jun10.txt \
  --feature_dir ${POPS_TOOL}/POPs_f_Joel/to_munge/ \
  --save_prefix ${POPS_TOOL}/POPs_f_Joel_munged/munged \
  --max_cols 50

conda deactivate
