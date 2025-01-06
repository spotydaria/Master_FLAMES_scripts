#!/bin/bash
#SBATCH --job-name=POPs_for_fuma      # Job name
#SBATCH --time=05:00:00              # Time limit (5 hours)
#SBATCH --mem=32G
#SBATCH --output=POPs_for_fuma_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drkramarenko@gmail.com

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate POPs

echo "WD_PROJECT: ${WD_PROJECT}"
echo "POPS_TOOL: ${POPS_TOOL}"
echo "FUMA_result_folder: ${FUMA_result_folder}"

if [ ! -d "${WD_PROJECT}/data/FUMA/${FUMA_result_folder}/POPS" ]; then
  echo "Creating directory ${WD_PROJECT}/data/FUMA/${FUMA_result_folder}/POPS"
  mkdir -p ${WD_PROJECT}/data/FUMA/${FUMA_result_folder}/POPS
fi
# Run the pops.py Python script with all parameters
python ${POPS_TOOL}/pops/pops.py \
--verbose \
--gene_annot_path ${POPS_TOOL}/pops/example/data/utils/gene_annot_jun10.txt \
--feature_mat_prefix ${POPS_TOOL}/POPs_f_Joel_munged/munged \
--num_feature_chunks 1155 \
--magma_prefix ${WD_PROJECT}/data/FUMA/${FUMA_result_folder}/magma \
--out_prefix ${WD_PROJECT}/data/FUMA/${FUMA_result_folder}/POPS/out



conda deactivate
