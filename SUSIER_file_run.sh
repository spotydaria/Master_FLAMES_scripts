#!/bin/bash
#SBATCH --job-name=SUSIER_run       # Job name
#SBATCH --time=04:00:00              
#SBATCH --mem=32G
#SBATCH --output=SUSIER_run_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drkramarenko@gmail.com


# Activate your Conda environment
# source /path/to/conda.sh  # This line loads Conda commands
source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate Finemap_2024

if [ ! -d "${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/" ]; then
  echo "Creating directory ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/"
  mkdir -p ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/
fi

if [ ! -d "${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/" ]; then
  echo "Creating directory ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/"
  mkdir -p ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/
fi

if [ ! -d "${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/orig_loc/" ]; then
  echo "Creating directory ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/orig_loc/"
  mkdir -p ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/orig_loc/
fi


# Rscript in_file_path output_path N(susie) L(susie)
Rscript ${WD_PROJECT}/scripts/SUSIER_file_run.R ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/ ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/  35512 10

conda deactivate