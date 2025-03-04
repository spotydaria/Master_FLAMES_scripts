## Environments:
# LAVA
conda activate LAVA_2024

# set veribles. to remove -> unset
export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

scp /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/intersection_table_per_region_1KB.xlsx dkramarenk@snellius.surf.nl://home/dkramarenk/projects/LAVA/DCM_HCM/data/

# /home/dkramarenk/projects/LAVA/DCM_HCM -> ${WD_PROJECT}
##########################################
#      0.  Downloading the data          #
##########################################
cd ${WD_PROJECT}


R
##########################################
library(readr)   # For reading text files
library(dplyr)    # For data manipulation
library(stringr)  # For extracting strings
library(tidyr) 
library(data.table)  # for fast file reading (fread)
library(readxl)

flip_beta_if_needed <- function(df, dataset_name) {
  if (grepl("HCM", dataset_name)) {
    df <- df %>%
      mutate(
        # Flip HCM Betas to be negative if needed
        Beta_Flipped = ifelse(BETA > 0, -BETA, BETA),
        A1_Flipped = ifelse(BETA > 0, A2, A1),
        A2_Flipped = ifelse(BETA > 0, A1, A2)
      ) %>%
      select(-BETA, -A1, -A2) %>%
      rename(BETA = Beta_Flipped, A1 = A1_Flipped, A2 = A2_Flipped)
  } else if (grepl("DCM", dataset_name)) {
    df <- df %>%
      mutate(
        # Flip DCM Betas to be positive if needed
        Beta_Flipped = ifelse(BETA < 0, -BETA, BETA),
        A1_Flipped = ifelse(BETA < 0, A2, A1),
        A2_Flipped = ifelse(BETA < 0, A1, A2)
      ) %>%
      select(-BETA, -A1, -A2) %>%
      rename(BETA = Beta_Flipped, A1 = A1_Flipped, A2 = A2_Flipped)
  }
  return(df)
}

path_base <- "/home/dkramarenk/projects/LAVA/DCM_HCM/data"
## Upload GWAS/MTAG
dcm_file <- file.path(path_base, "sumst_processed", "DCM_GWAS_37_exclMYBPC3reg.txt")
hcm_file <- file.path(path_base, "sumst_processed", "HCM_GWAS_37_exclMYBPC3reg.txt")
cc_file  <- file.path(path_base, "sumst_processed",  "CC_GWAS_37_exclMYBPC3reg.txt")
dcm_mtag_file <- file.path(path_base, "sumst_processed", "DCM_MTAG_37_exclMYBPC3reg.txt")
hcm_mtag_file <- file.path(path_base, "sumst_processed", "HCM_MTAG_37_exclMYBPC3reg.txt")
cc_mtag_file  <- file.path(path_base, "sumst_processed",  "CC_MTAG_37_exclMYBPC3reg.txt")

# Read each. For large files, fread() is faster. If they're small, read.table() works too.
dcm_original <- fread(dcm_file) %>% mutate(Dataset = "DCM_GWAS")%>%
  mutate(chrpos = paste(CHR, BP, A1, A2, sep=":")) %>%
  flip_beta_if_needed("DCM_GWAS")

hcm_original <- fread(hcm_file) %>% mutate(Dataset = "HCM_GWAS")%>%
  mutate(chrpos = paste(CHR, BP, A1, A2, sep=":"))

hcm_original1 <- hcm_original %>%
  flip_beta_if_needed("HCM_GWAS")


cc_original  <- fread(cc_file)  %>% mutate(Dataset = "CC_GWAS")%>%
  mutate(chrpos = paste(CHR, BP, A1, A2, sep=":"))
dcm_mtag_original <- fread(dcm_mtag_file) %>% mutate(Dataset = "DCM_MTAG")%>%
  mutate(chrpos = paste(CHR, BP, A1, A2, sep=":"))
dcm_mtag_original1 <- dcm_mtag_original %>%
  mutate(Dataset = "DCM_MTAG", chrpos = paste(CHR, BP, A1, A2, sep=":")) %>%
  flip_beta_if_needed("DCM_MTAG")
hcm_mtag_original <- fread(hcm_mtag_file) %>% mutate(Dataset = "HCM_MTAG")%>%
  mutate(chrpos = paste(CHR, BP, A1, A2, sep=":"))
hcm_mtag_original1 <- hcm_mtag_original %>%
  mutate(Dataset = "HCM_MTAG", chrpos = paste(CHR, BP, A1, A2, sep=":")) %>%
  flip_beta_if_needed("DCM_MTAG")

cc_mtag_original  <- fread(cc_mtag_file)  %>% mutate(Dataset = "CC_MTAG")%>%
  mutate(chrpos = paste(CHR, BP, A1, A2, sep=":"))

# make a list of unique CHR:POS
chrposs <- unique(c(
  dcm_original$chrpos,
  hcm_original$chrpos,
  cc_original$chrpos,
  dcm_mtag_original$chrpos,
  hcm_mtag_original$chrpos,
  cc_mtag_original$chrpos
))

rsIDs <- unique(c(
  dcm_original$rsID,
  hcm_original$rsID,
  cc_original$rsID,
  dcm_mtag_original$rsID,
  hcm_mtag_original$rsID,
  cc_mtag_original$rsID
))

dataset_rsIDs <- data.frame(rsID = rsIDs)

joined_dataset <- dataset_rsIDs %>%
  left_join(dcm_original %>% transmute(rsID, BETA, SE, P, chrpos = paste(CHR, BP, A1, A2, sep=":"), chrpos_flip = paste(CHR, BP, A2, A1, sep=":")), by = "rsID") %>%
  dplyr::rename(BETA_DCM = BETA, SE_DCM = SE, P_DCM = P, CHRPOS_DCM = chrpos, CHRPOS_DCM_flip = chrpos_flip) %>%
  left_join(hcm_original %>% transmute(rsID, BETA, SE, P, chrpos = paste(CHR, BP, A1, A2, sep=":"), chrpos_flip = paste(CHR, BP, A2, A1, sep=":")), by = "rsID") %>%
  dplyr::rename(BETA_HCM = BETA, SE_HCM = SE, P_HCM = P, CHRPOS_HCM = chrpos, CHRPOS_HCM_flip = chrpos_flip) %>%
  left_join(cc_original %>% transmute(rsID, BETA, SE, P, chrpos = paste(CHR, BP, A1, A2, sep=":"), chrpos_flip = paste(CHR, BP, A2, A1, sep=":")), by = "rsID") %>%
  dplyr::rename(BETA_CC = BETA, SE_CC = SE, P_CC = P, CHRPOS_CC = chrpos, CHRPOS_CC_flip = chrpos_flip) %>%
  left_join(dcm_mtag_original %>% transmute(rsID, BETA, SE, P, chrpos = paste(CHR, BP, A1, A2, sep=":"), chrpos_flip = paste(CHR, BP, A2, A1, sep=":")), by = "rsID") %>%
  dplyr::rename(BETA_DCM_MTAG = BETA, SE_DCM_MTAG = SE, P_DCM_MTAG = P, CHRPOS_DCM_MTAG = chrpos, CHRPOS_DCM_MTAG_flip = chrpos_flip) %>%
  left_join(hcm_mtag_original %>% transmute(rsID, BETA, SE, P, chrpos = paste(CHR, BP, A1, A2, sep=":"), chrpos_flip = paste(CHR, BP, A2, A1, sep=":")), by = "rsID") %>%
  dplyr::rename(BETA_HCM_MTAG = BETA, SE_HCM_MTAG = SE, P_HCM_MTAG = P, CHRPOS_HCM_MTAG = chrpos, CHRPOS_HCM_MTAG_flip = chrpos_flip) %>%
  left_join(cc_mtag_original %>% transmute(rsID, BETA, SE, P, chrpos = paste(CHR, BP, A1, A2, sep=":"), chrpos_flip = paste(CHR, BP, A2, A1, sep=":")), by = "rsID") %>%
  dplyr::rename(BETA_CC_MTAG = BETA, SE_CC_MTAG = SE, P_CC_MTAG = P, CHRPOS_CC_MTAG = chrpos, CHRPOS_CC_MTAG_flip = chrpos_flip)

joined_dataset_2 <- joined_dataset %>%
   mutate(
    CHRPOS_REF = coalesce(CHRPOS_DCM_MTAG, CHRPOS_DCM, CHRPOS_HCM_MTAG, CHRPOS_HCM, CHRPOS_CC_MTAG, CHRPOS_CC),
    CHRPOS_DCM_final = ifelse(CHRPOS_DCM == CHRPOS_REF, CHRPOS_DCM, 
                              ifelse(CHRPOS_DCM_flip == CHRPOS_REF, CHRPOS_DCM_flip, NA)),
    BETA_DCM_final = ifelse(CHRPOS_DCM == CHRPOS_REF, BETA_DCM, 
                            ifelse(CHRPOS_DCM_flip == CHRPOS_REF, -BETA_DCM, NA)),
    CHRPOS_HCM_final = ifelse(CHRPOS_HCM == CHRPOS_REF, CHRPOS_HCM, 
                              ifelse(CHRPOS_HCM_flip == CHRPOS_REF, CHRPOS_HCM_flip, NA)),
    BETA_HCM_final = ifelse(CHRPOS_HCM == CHRPOS_REF, BETA_HCM, 
                            ifelse(CHRPOS_HCM_flip == CHRPOS_REF, -BETA_HCM, NA)),

    CHRPOS_CC_final = ifelse(CHRPOS_CC == CHRPOS_REF, CHRPOS_CC, 
                             ifelse(CHRPOS_CC_flip == CHRPOS_REF, CHRPOS_CC_flip, NA)),
    BETA_CC_final = ifelse(CHRPOS_CC == CHRPOS_REF, BETA_CC, 
                           ifelse(CHRPOS_CC_flip == CHRPOS_REF, -BETA_CC, NA)),

    CHRPOS_DCM_MTAG_final = ifelse(CHRPOS_DCM_MTAG == CHRPOS_REF, CHRPOS_DCM_MTAG, 
                                   ifelse(CHRPOS_DCM_MTAG_flip == CHRPOS_REF, CHRPOS_DCM_MTAG_flip, NA)),
    BETA_DCM_MTAG_final = ifelse(CHRPOS_DCM_MTAG == CHRPOS_REF, BETA_DCM_MTAG, 
                                 ifelse(CHRPOS_DCM_MTAG_flip == CHRPOS_REF, -BETA_DCM_MTAG, NA)),

    CHRPOS_HCM_MTAG_final = ifelse(CHRPOS_HCM_MTAG == CHRPOS_REF, CHRPOS_HCM_MTAG, 
                                   ifelse(CHRPOS_HCM_MTAG_flip == CHRPOS_REF, CHRPOS_HCM_MTAG_flip, NA)),
    BETA_HCM_MTAG_final = ifelse(CHRPOS_HCM_MTAG == CHRPOS_REF, BETA_HCM_MTAG, 
                                 ifelse(CHRPOS_HCM_MTAG_flip == CHRPOS_REF, -BETA_HCM_MTAG, NA)),

    CHRPOS_CC_MTAG_final = ifelse(CHRPOS_CC_MTAG == CHRPOS_REF, CHRPOS_CC_MTAG, 
                                  ifelse(CHRPOS_CC_MTAG_flip == CHRPOS_REF, CHRPOS_CC_MTAG_flip, NA)),
    BETA_CC_MTAG_final = ifelse(CHRPOS_CC_MTAG == CHRPOS_REF, BETA_CC_MTAG, 
                                ifelse(CHRPOS_CC_MTAG_flip == CHRPOS_REF, -BETA_CC_MTAG, NA))
  )

print("joined_dataset_2")
final_dataset_2 <- joined_dataset_2 %>%
  transmute(
    rsID,
    CHRPOS_REF,
    DCM_CHRPOS = CHRPOS_DCM_final, DCM_BETA = BETA_DCM_final, DCM_P = P_DCM, DCM_SE = SE_DCM,
    HCM_CHRPOS = CHRPOS_HCM_final, HCM_BETA = BETA_HCM_final, HCM_P = P_HCM, HCM_SE = SE_HCM,
    CC_CHRPOS = CHRPOS_CC_final, CC_BETA = BETA_CC_final, CC_P = P_CC, CC_SE = SE_CC,
    DCM_MTAG_CHRPOS = CHRPOS_DCM_MTAG_final, DCM_MTAG_BETA = BETA_DCM_MTAG_final, DCM_MTAG_P = P_DCM_MTAG, DCM_MTAG_SE = SE_DCM_MTAG,
    HCM_MTAG_CHRPOS = CHRPOS_HCM_MTAG_final, HCM_MTAG_BETA = BETA_HCM_MTAG_final, HCM_MTAG_P = P_HCM_MTAG, HCM_MTAG_SE = SE_HCM_MTAG,
    CC_MTAG_CHRPOS = CHRPOS_CC_MTAG_final, CC_MTAG_BETA = BETA_CC_MTAG_final, CC_MTAG_P = P_CC_MTAG, CC_MTAG_SE = SE_CC_MTAG
  ) %>%
    mutate(CHRPOS = coalesce(!!!select(., ends_with("_CHRPOS"))),     #loci_n = coalesce(!!!select(., ends_with("_loci"))),
         SNP = sub("^([^:]+:[^:]+):.*$", "\\1", CHRPOS)) %>%
 separate(
    SNP,
    into = c("CHR", "BP"),
    sep = ":",
    remove = FALSE
  ) %>%
  mutate(CHR = as.numeric(CHR))

print("final_dataset_2")

loci_file <- "/home/dkramarenk/projects/LAVA/DCM_HCM/data/intersection_table_per_region_1KB.xlsx"
my_loci <- read_excel(loci_file, sheet = "main") %>%
  dplyr::rename(loci_n = group_id)

print("my_loci")

final_dataset_all_loci <- final_dataset_2   %>%
  mutate(CHR= as.numeric(CHR), BP = as.numeric(BP)) %>%
  left_join(my_loci %>% transmute(CHR = as.numeric(CHR), POS_start, POS_end, loci_n), by = "CHR", relationship = "many-to-many") %>%   
    filter(
      BP >= POS_start,       # SNP within locus start
      BP <= POS_end         # SNP within locus end
    ) 


print("writing")

write.csv(final_dataset_all_loci, "data/DCM_HCM_CC_dataset_all_loci_flipped.csv", row.names = FALSE, na = "")
print("written")

##########################################

#!/bin/bash
#SBATCH --job-name=loci_inter       # Job name
#SBATCH --time=04:00:00              
#SBATCH --mem=32G
#SBATCH --output=loci_inter_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=drkramarenko@gmail.com

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate LAVA_2024

Rscript scripts/loci_inter.R

##########################################

sbatch scripts/loci_inter.sh
squeue -u dkramarenk

##########################################
#!/bin/bash

# Define the base path
path_base="/home/dkramarenk/projects/LAVA/DCM_HCM/data"

# Define the file paths
files=(
    "$path_base/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt"
    "$path_base/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt"
    "$path_base/sumst_processed/CC_GWAS_37_exclMYBPC3reg.txt"
    "$path_base/sumst_processed/DCM_MTAG_37_exclMYBPC3reg.txt"
    "$path_base/sumst_processed/HCM_MTAG_37_exclMYBPC3reg.txt"
    "$path_base/sumst_processed/CC_MTAG_37_exclMYBPC3reg.txt"
)

# Loop through each file and create a new file with first 100 lines
for file in "${files[@]}"; do
    if [[ -f "$file" ]]; then
        # Extract filename from path
        filename=$(basename "$file")
        # Define the new filename
        new_file="$path_base/sumst_processed/${filename}"
        # Extract the first 100 lines and save it
        head -n 10000 "$file" > "$new_file"
        echo "Created $new_file with the first 100 lines of $filename"
    else
        echo "File not found: $file"
    fi
done


chmod +x scripts/extract_head100.sh
./scripts/extract_head100.sh

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}/data/DCM_HCM_CC_dataset_all_loci_flipped.csv /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/sumst_processed/
