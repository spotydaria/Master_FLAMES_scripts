#!/bin/bash
#SBATCH --job-name=FLAMES_run     # Job name
#SBATCH --time=05:00:00              
#SBATCH --mem=32G
#SBATCH --output=FLAMES_run_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drkramarenko@gmail.com

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate FLAMES

export INDEX_FILE_NCLUDING_COLUMN="${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/FLAMES_index.ind"
export DESIRED_OUTPUT_DIRECTORY="${WD_PROJECT}/FLAMES_out/${FUMA_result_folder}/"

python ${WD_PROJECT}/FLAMES/FLAMES.py FLAMES \
-id $INDEX_FILE_NCLUDING_COLUMN \
-o $DESIRED_OUTPUT_DIRECTORY

conda deactivate
