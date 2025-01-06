#!/bin/bash
#SBATCH --job-name=compare_SNPs       # Job name
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=compare_SNPs_%j.out

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate LAVA_2024

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_sumst_proc="/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed"

Rscript ${WD_PROJECT}/scripts/Compare_SNPs.R --folder_path=${WD_sumst_proc}

conda deactivate
