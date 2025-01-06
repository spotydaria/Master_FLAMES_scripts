#!/bin/bash
#SBATCH --job-name=HCM_GWAS_process       # Job name
#SBATCH --time=01:00:00              
#SBATCH --mem=32G
#SBATCH --output=HCM_GWAS_process_%j.out

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate LAVA_2024

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

Rscript ${WD_PROJECT}/scripts/HCM_GWAS_process.R \
-f data/raw/hcm.gwama.250621.gel.nofin.txt \
-o data/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt

conda deactivate
