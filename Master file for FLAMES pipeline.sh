##########################################
#      Master file for FLAMES pipeline
# Date: 09/12/2024
##########################################

## Environments:
# LAVA
conda activate LAVA_2024 

# set veribles. to remove -> unset
export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

# /home/dkramarenk/projects/LAVA/DCM_HCM -> ${WD_PROJECT}
##########################################
#      0.  Downloading the data          #
##########################################
cd ${WD_PROJECT}

# HCM GWAS:  - build 37
# cp /archive/expcard/Lisa_Projects/HCM-GWAS/HCM-2022/hcm.gwama.250621.gel.nofin.txt ${WD_PROJECT}/data/raw/
# processed: ${WD_PROJECT}/data/HCM_stats_LAVA.txt

## HCM MTAG: 
# cp /nfs/archive01/expcard/Lisa_Projects/HCM-GWAS/HCM-2022/HCM-global_Francis-clusters_trait_1.txt ${WD_PROJECT}/data/raw/

## DCM GWAS: 
# cp /archive/expcard/DCM_GWAS_2023/data/summary_stats/MTAG/Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv.gz ${WD_PROJECT}/data/raw/
# processed: ${WD_PROJECT}/data/DCM_stats_LAVA.txt

## DCM MTAG: 
# cp /gpfs/work2/0/brugada/sjurgens/DCM_GWAS/data/META/V2/DCM_MTAG_Ecc_global_LVESVi_processed.tsv.gz ${WD_PROJECT}/data/raw/
# for DCM MTAG, you can use the mtag_Neff column for sample size filtering

## CC GWAS
# cp /gpfs/work2/0/brugada/sjurgens/DCM_HCM_shared_opposing/data/CC_GWAS__DCM__HCM_exclMYBPC3reg.txt.gz ${WD_PROJECT}/data/raw/

## CC_GWAS + MTAG with MRI GWAS
# raw & processed: /gpfs/work2/0/brugada/sjurgens/DCM_HCM_shared_opposing/data/CC_GWAS__DCM__HCM__MTAG_Ecc_LVESVi_LVconc.txt.gz
# cd /gpfs/work2/0/brugada/sjurgens/DCM_HCM_shared_opposing/data/CC_GWAS__DCM__HCM__MTAG_Ecc_LVESVi_LVconc.txt.gz ${WD_PROJECT}/data/raw/
# gunzip CC_GWAS__DCM__HCM__MTAG_Ecc_LVESVi_LVconc.txt.gz


#  1000G
cp /home/dkramarenk/projects/Project_DCM_GWAS/POP/1000_genomes_ALEX/1KGautos_fil_biallelic* ${WD_PROJECT}/data/1000G

## Files        Preprocessed        FUMA
## CC GWAS      + (exclMYBPC3reg)   +
## CC MTAG      + (exclMYBPC3reg)
## DCM GWAS     + full script
## DCM MTAG     + full script
## HCM GWAS     +
## HCM MTAG     +

##########################################
#      How file should look like        #
##########################################

# # SumStats
# CHR BP SNPID_UKB A1 A2 MAF BETA SE STAT P NMISS INFO_UKB
# 1 113418061 rs2360008 G A 0.2261 -0.01392 0.01282 -1.086 0.2775 379831 1
# 1 113418144 rs9429498 T C 0.001218 0.2492 0.1538 1.62 0.1051 379682 0.951417
# 1 113418415 rs1237670 G A 0.2242 -0.01563 0.01286 -1.215 0.2244 379470 0.998238

# # Input
# phenotype	cases	controls	filename
# asthma	44301	341521	vignettes/data/asthma.sumstats.txt
# bmi	NA	NA	vignettes/data/bmi.sumstats.txt
# depression	170756	329443	vignettes/data/depression.sumstats.txt
# diabetes	18483	366937	vignettes/data/diabetes.sumstats.txt
# hypothyroidism	13043	231847	vignettes/data/hypothyroidism.sumstats.txt
# neuro	NA	NA	vignettes/data/neuro.sumstats.txt
# rheuma	14361	43923	vignettes/data/rheuma.sumstats.txt

##########################################
#     1.1    Files processing            #
##########################################

sbatch ${WD_PROJECT}/scripts/DCM_MTAG_process.sh
sbatch ${WD_PROJECT}/scripts/HCM_MTAG_process.sh
sbatch ${WD_PROJECT}/scripts/CC_GWAS_process.sh
sbatch ${WD_PROJECT}/scripts/CC_MTAG_process.sh
sbatch ${WD_PROJECT}/scripts/DCM_GWAS_process.sh
sbatch ${WD_PROJECT}/scripts/HCM_GWAS_process.sh

####################################################################################

##########################################
# 1.2  Check so all SNPs on the same ref #
##########################################
sbatch ${WD_PROJECT}/scripts/compare_SNPs.sh

##########################################
#     1.3     Transfer to PC             #
##########################################

scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}/data/sumst_processed/* /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/sumst_processed/

##########################################
# 1.4   FUMA/MAGMA -> files to SNELLIUS  #
##########################################

scp /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_DCM/LAVA/FLAMES/FUMA_results/*.zip dkramarenk@snellius.surf.nl:/${WD_PROJECT}/data/FUMA
for zip in *.zip; do
    # Create a directory with the name of the zip file minus the .zip extension
    dir="${zip%.zip}"
    mkdir -p "$dir"
    # Unzip the contents into the newly created directory
    unzip "$zip" -d "$dir"
done

cd  ${WD_PROJECT}/data/FUMA/

mv FUMA_CCGWAS2 FUMA_CC_GWAS
mv FUMA_CC_MTAG1 FUMA_CC_MTAG
mv FUMA_DCM_GWAS1 FUMA_DCM_GWAS
mv FUMA_HCM_GWAS1 FUMA_HCM_GWAS
mv FUMA_DCM_MTAG1 FUMA_DCM_MTAG
mv FUMA_HCM_MTAG1 FUMA_HCM_MTAG

##########################################
# 2.1 Run PoPS on the generated MAGMA z-scores.
##########################################

conda activate POPs
# To deactivate an active environment, use
#     $ conda deactivate
export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export POPS_TOOL="/home/dkramarenk/projects/tools/POPs"

##########################################
# 2.2 Munge features
##########################################

# Features
wget https://www.dropbox.com/sh/o6t5jprvxb8b500/AABK2h7UIw3K85oHE85Ey2ZNa/data/PoPS.features.txt.gz?dl=1
gunzip PoPS.features.txt.gz

sbatch split_pops.sh
##Split 21GB file into smaller chunks

sbatch munge_joel.sh

##########################################
# 2.3 POPs 
##########################################

# # Convert tab-delimited to space-delimited
# awk 'BEGIN{FS=OFS="\t"} {$1=$1}1' ${WD_PROJECT}/data/FUMA/FUMA_CCGWAS2/magma.genes.out | column -t > ${WD_PROJECT}/data/FUMA/FUMA_CCGWAS2/magma_scores/magma.genes.out
# cp ${WD_PROJECT}/data/FUMA/FUMA_CCGWAS2/magma.genes.raw ${WD_PROJECT}/data/FUMA/FUMA_CCGWAS2/magma_scores/

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

sbatch --export=FUMA_result_folder="FUMA_CC_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",POPS_TOOL="/home/dkramarenk/projects/tools/POPs" ${WD_PROJECT}/scripts/POPs_for_fuma.sh
sbatch --export=FUMA_result_folder="FUMA_DCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",POPS_TOOL="/home/dkramarenk/projects/tools/POPs" ${WD_PROJECT}/scripts/POPs_for_fuma.sh
sbatch --export=FUMA_result_folder="FUMA_DCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",POPS_TOOL="/home/dkramarenk/projects/tools/POPs" ${WD_PROJECT}/scripts/POPs_for_fuma.sh
sbatch --export=FUMA_result_folder="FUMA_HCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",POPS_TOOL="/home/dkramarenk/projects/tools/POPs" ${WD_PROJECT}/scripts/POPs_for_fuma.sh
sbatch --export=FUMA_result_folder="FUMA_HCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",POPS_TOOL="/home/dkramarenk/projects/tools/POPs" ${WD_PROJECT}/scripts/POPs_for_fuma.sh
sbatch --export=FUMA_result_folder="FUMA_CC_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",POPS_TOOL="/home/dkramarenk/projects/tools/POPs" ${WD_PROJECT}/scripts/POPs_for_fuma.sh

##########################################
#     3. Format your credible sets.
##########################################

##########################################
#           3.1 Finemapping
##########################################
# conda create -n Finemap_2024
## Package Plan ##
#  environment location: /home/dkramarenk/.conda/envs/Finemap_2024
#
# To activate this environment, use
conda activate Finemap_2024
# To deactivate an active environment, use
#     $ conda deactivate

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export SUMST="DCM_GWAS_37_exclMYBPC3reg.txt"
export SUMST_RF="37_38_build.txt"
export SUSIE_OUT="DCM_GWAS"
cd ${WD_PROJECT}

sbatch  --export=SUSIE_OUT="DCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",SUMST_RF="37_38_build.txt" ${WD_PROJECT}/scripts/SUSIER_file_prep.sh
sbatch  --export=SUSIE_OUT="CC_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",SUMST_RF="37_38_build.txt" ${WD_PROJECT}/scripts/SUSIER_file_prep.sh
sbatch  --export=SUSIE_OUT="HCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",SUMST_RF="37_38_build.txt" ${WD_PROJECT}/scripts/SUSIER_file_prep.sh
sbatch  --export=SUSIE_OUT="DCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",SUMST_RF="37_38_build.txt" ${WD_PROJECT}/scripts/SUSIER_file_prep.sh
sbatch  --export=SUSIE_OUT="CC_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",SUMST_RF="37_38_build.txt" ${WD_PROJECT}/scripts/SUSIER_file_prep.sh
sbatch  --export=SUSIE_OUT="HCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM",SUMST_RF="37_38_build.txt" ${WD_PROJECT}/scripts/SUSIER_file_prep.sh


# Susie run
sbatch --export=SUSIE_OUT="DCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/SUSIER_file_run.sh
sbatch --export=SUSIE_OUT="CC_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/SUSIER_file_run.sh
sbatch --export=SUSIE_OUT="DCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/SUSIER_file_run.sh
sbatch --export=SUSIE_OUT="CC_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/SUSIER_file_run.sh
sbatch --export=SUSIE_OUT="HCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/SUSIER_file_run.sh
sbatch --export=SUSIE_OUT="HCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/SUSIER_file_run.sh

squeue -u dkramarenk

# scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}/data/SUSIE/CCGWAS/regions_to_analyze.txt /Users/drkramarenko/OneDrive/Computation/GWAS/Project_DCM/SUSIER/
##########################################
#        4.1 FLAMES annotate
##########################################

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
cd ${WD_PROJECT}

squeue -u dkramarenk

export FUMA_result_folder="CC_GWAS"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="CC_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_annotate.sh
# done

export FUMA_result_folder="DCM_GWAS"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="DCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_annotate.sh
# done

export FUMA_result_folder="HCM_GWAS"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="HCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_annotate.sh
# done

export FUMA_result_folder="DCM_MTAG"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="DCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_annotate.sh
# done

export FUMA_result_folder="HCM_MTAG"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="HCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_annotate.sh
# done

export FUMA_result_folder="CC_MTAG"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="CC_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_annotate.sh


# Check out
cat ${WD_PROJECT}/data/SUSIE/*/out/FLAMES_annotate_*

##########################################
#            4.2 FLAMES run
##########################################

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"    
cd ${WD_PROJECT}


export FUMA_result_folder="CC_GWAS"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="CC_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_run.sh

export FUMA_result_folder="DCM_GWAS"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="DCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_run.sh
# done

export FUMA_result_folder="HCM_GWAS"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="HCM_GWAS",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_run.sh

export FUMA_result_folder="DCM_MTAG"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="DCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_run.sh
# done

export FUMA_result_folder="HCM_MTAG"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="HCM_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_run.sh

export FUMA_result_folder="CC_MTAG"
cd ${WD_PROJECT}/data/SUSIE/${FUMA_result_folder}/out/
sbatch --export=FUMA_result_folder="CC_MTAG",WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM" ${WD_PROJECT}/scripts/FLAMES_run.sh

squeue -u dkramarenk

##########################################
#           5. Copy files to pc
##########################################

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}/FLAMES_out/* /Users/drkramarenko/OneDrive/Computation/GWAS/Project_SHARED_DCM_HCM/FLAMES_out/

for SUSIE_OUT in "DCM_GWAS" "HCM_GWAS" "CC_GWAS" "CC_MTAG" "DCM_MTAG" "HCM_MTAG"; do
  # Define the local target directory
  LOCAL_DIR="/Users/drkramarenko/OneDrive/Computation/GWAS/Project_SHARED_DCM_HCM/FLAMES_out_${SUSIE_OUT}/"
  
  # Check if the directory exists, create it if it doesn't
  if [ ! -d "$LOCAL_DIR" ]; then
    echo "Creating directory: $LOCAL_DIR"
    mkdir -p "$LOCAL_DIR"
  fi
  
  # Copy files from remote to local
  echo "Copying files from ${WD_PROJECT}/FLAMES_out/${SUSIE_OUT} to $LOCAL_DIR"
  scp dkramarenk@snellius.surf.nl:"${WD_PROJECT}/FLAMES_out/${SUSIE_OUT}/*" "$LOCAL_DIR/"
done


for SUSIE_OUT in "DCM_GWAS" "HCM_GWAS" "CC_GWAS" "CC_MTAG" "DCM_MTAG" "HCM_MTAG"; do
  # Define the local target directory
  LOCAL_DIR="/Users/drkramarenko/OneDrive/Computation/GWAS/Project_SHARED_DCM_HCM/SUSIER_${SUSIE_OUT}/out"
  
  # Check if the directory exists, create it if it doesn't
  if [ ! -d "$LOCAL_DIR" ]; then
    echo "Creating directory: $LOCAL_DIR"
    mkdir -p "$LOCAL_DIR"
  fi
  
  # Copy files from remote to local
  echo "Copying files from ${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out to $LOCAL_DIR"
  scp dkramarenk@snellius.surf.nl:"${WD_PROJECT}/data/SUSIE/${SUSIE_OUT}/out/*" "$LOCAL_DIR/"
done
# Define the project directory

ll /Users/drkramarenko/OneDrive/Computation/GWAS/Project_SHARED_DCM_HCM/SUSIER_

