#!/bin/bash
#SBATCH --job-name=FLAMES_annotate       # Job name
#SBATCH --time=04:00:00              
#SBATCH --mem=32G
#SBATCH --output=FLAMES_annotate_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drkramarenko@gmail.com
source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate FLAMES

export PATH_TO_DOWNLOADED_FEATURES="/home/dkramarenk/projects/tools/POPs/POPs_f_Joel_munged"
export PATH_TO_GENERATED_MAGMA_Z_SCORES="${WD_PROJECT}/data/FUMA/FUMA_${FUMA_result_folder}/magma"
export DESIRED_POPS_OUTPUT_PREFIX="${WD_PROJECT}/data/FUMA/FUMA_${FUMA_result_folder}/POPS/out"
export DESIRED_OUTPUT_DIRECTORY="${WD_PROJECT}/FLAMES_out/${FUMA_result_folder}/"
export PATH_TO_THE_DOWNLOADED_ANNOTATION_DATA_DIRECTORY="/home/dkramarenk/projects/tools/POPs/Annotation_data/"
export PATH_TO_INDEXFILE="${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/FLAMES_index.ind"

if [ ! -d "${WD_PROJECT}/FLAMES_out/${FUMA_result_folder}" ]; then
  echo "Creating directory ${WD_PROJECT}/FLAMES_out/${FUMA_result_folder}"
  mkdir -p ${WD_PROJECT}/FLAMES_out/${FUMA_result_folder}
fi

python ${WD_PROJECT}/FLAMES/FLAMES.py annotate \
-o $DESIRED_OUTPUT_DIRECTORY \
-a $PATH_TO_THE_DOWNLOADED_ANNOTATION_DATA_DIRECTORY \
-p ${DESIRED_POPS_OUTPUT_PREFIX}.preds \
-m ${PATH_TO_GENERATED_MAGMA_Z_SCORES}.genes.out \
-mt ${PATH_TO_GENERATED_MAGMA_Z_SCORES}.gsa.out \
-id $PATH_TO_INDEXFILE

conda deactivate
