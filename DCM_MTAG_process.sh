#!/bin/bash
#SBATCH --job-name=DCM_MTAG_process       # Job name
#SBATCH --time=01:00:00              
#SBATCH --mem=32G
#SBATCH --output=DCM_MTAG_process_%j.out

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate LAVA_2024


cd /home/dkramarenk/projects/LAVA/DCM_HCM

Rscript /home/dkramarenk/projects/LAVA/DCM_HCM/scripts/DCM_MTAG_process.R \
-f data/raw/DCM_MTAG_Ecc_global_LVESVi_processed.tsv \
-o data/sumst_processed/DCM_MTAG_37_exclMYBPC3reg.txt

conda deactivate
