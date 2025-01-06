#!/bin/bash
#SBATCH --job-name=DCM_GWAS_process       # Job name
#SBATCH --time=01:00:00              
#SBATCH --mem=32G
#SBATCH --output=DCM_GWAS_process_%j.out

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate LAVA_2024

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

Rscript ${WD_PROJECT}/scripts/DCM_GWAS_process.R \
-f data/raw/Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv \
-o data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt

conda deactivate