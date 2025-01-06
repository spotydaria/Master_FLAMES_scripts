#!/bin/bash
#SBATCH --job-name=SUSIER_file_prep       # Job name
#SBATCH --time=01:00:00              
#SBATCH --mem=16G
#SBATCH --output=SUSIER_file_prep_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drkramarenko@gmail.com

# Check if the necessary folders exist and create them if they do not
if [ ! -d "${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/" ]; then
  echo "Creating directory ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/"
  mkdir -p ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/
fi

if [ ! -d "${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/" ]; then
  echo "Creating directory ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/"
  mkdir -p ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/
fi

if [ ! -d "${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/orig_loc" ]; then
  echo "Creating directory ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/orig_loc"
  mkdir -p ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/orig_loc
fi

# Activate your Conda environment
# source /path/to/conda.sh  # This line loads Conda commands
source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate Finemap_2024

Rscript ${WD_PROJECT}/scripts/SUSIER_file_prep.R ${WD_PROJECT}/data/sumst_processed/${SUMST_RF} ${WD_PROJECT}/data/sumst_processed/${SUSIE_OUT}_37_exclMYBPC3reg.txt /gpfs/work2/0/brugada/tools/LDmatrix/UKBB_LD_matrixes_SJ/ ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/  

conda deactivate
