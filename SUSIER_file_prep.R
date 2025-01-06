# Main Script for Processing Summary Statistics to prepare for SuSiE
# Load necessary libraries
library(data.table)  # Efficient data manipulation
library(dplyr)       # Data manipulation and transformation
library(readr)       # Data input and output
library(tidyr)       # Data tidying
library(susieR)      # Implements SuSiE (Sum of Single Effects)
library(stringr)     # String manipulation
library(Rfast)       # Fast data manipulation


# Parse arguments for main script
args <- commandArgs(trailingOnly = TRUE)
# Example usage with paths for debugging (remove this line in production)
# args <- c("/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumstats/CC_GWAS__DCM__HCM_exclMYBPC3reg_hg38.txt", "/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumstats/CC_GWAS__DCM__HCM_exclMYBPC3reg.txt","/gpfs/work2/0/brugada/tools/LDmatrix/UKBB_LD_matrixes_SJ/", "/home/dkramarenk/projects/LAVA/DCM_HCM/data/SUSIE/CCGWAS/")
# args <- c("/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/37_38_build.txt", "/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt", "/gpfs/work2/0/brugada/tools/LDmatrix/UKBB_LD_matrixes_SJ/", "/home/dkramarenk/projects/LAVA/DCM_HCM/data/SUSIE/HCM_GWAS/")

# Check that correct number of arguments is provided
if (length(args) != 4) {
  stop("Usage: Rscript script_name.R sumstats_path sumstats_ref_path ld_matrix_path output_path")
}

# Assign arguments to variables
sumstats_ref_path <- args[1]
sumstats_path <- args[2]
ld_matrix_path <- args[3]
output_path <- args[4]
# N_CC <- as.numeric(args[5])

get_chr_start_end <- function(filename) {
  # Extract chromosome, start, and end positions from filename
  chr <- str_extract(filename, "(?<=chr)[0-9]+")
  start <- as.numeric(str_extract(filename, "(?<=snp\\.)[0-9]+(?=_)"))
  end <- as.numeric(str_extract(filename, "(?<=_)[0-9]+(?=\\.Rvar)"))
  return(list(chr = as.numeric(chr), start = start, end = end))
}

# Load summary statistics and LD matrix file paths
sumstats_ref <- fread(sumstats_ref_path)

# HCM <- fread("/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt")
# DCM <- fread("/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt")

sumst <- fread(sumstats_path) %>%
rename(CHR_37=CHR, BP_37=BP) %>%
left_join(sumstats_ref %>% select(!(rsID)), by = c("CHR_37", "BP_37"))

# Filter significant SNPs
sumst_filt <- sumst %>%
  filter(P < 5e-8)

# Extracting the start and end positions from the file names
ld_files <- list.files(ld_matrix_path, pattern = "UKB_b38_0.1_chr.*\\.Rvar", full.names = TRUE)

# Initialize an empty column in sumst_filt for the file names
sumst_filt$file <- NA

# Loop through each LD matrix file
for (file in ld_files) {
  # Print current file being processed
  # print(paste("Processing file:", file))
  
  # Get chromosome, start, and end positions
  file_info <- get_chr_start_end(file)
  
  # Filter sumst_filt for SNPs that are in the corresponding region
  matches <- with(sumst_filt, CHR_38 == file_info$chr & BP_38 >= file_info$start & BP_38 <= file_info$end)
  
  # Add the file name to the matched rows
  sumst_filt$file[matches] <- file
}


# Group by the file and summarize with min and max BP_38
summarized_files <- sumst_filt %>%
  group_by(file) %>%
  summarise(min_BP_38 = min(BP_38), max_BP_38 = max(BP_38), num_SNP = n()) %>%
  mutate(
    BP_38_diff = max_BP_38 - min_BP_38,
    chr_LD = str_extract(file, "(?<=chr)[0-9]+"),
    start_LD = as.numeric(str_extract(file, "(?<=snp\\.)[0-9]+(?=_)")),
    end_LD = as.numeric(str_extract(file, "(?<=_)[0-9]+(?=\\.Rvar)")),
    loci_n = row_number(),
    LD_edge_low = min_BP_38 - 5e5,
    LD_edge_high = max_BP_38 + 5e5
  ) %>%
  arrange(file)

# Write output for summarized files
write.table(summarized_files,
            file = file.path(output_path, "regions_to_analyze.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
            
summarized_files_noNA <- summarized_files %>%
filter(!is.na(file))

# Create and save final summary statistics
merged_data <- merge(sumst, bind_rows(lapply(summarized_files_noNA$file, fread)), by.x = c("CHR_38", "BP_38"), by.y = c("chrom", "pos")) %>%
filter(!is.na(CHR_38), !is.na(BP_38)) %>%
  mutate(
    flip = case_when(
      A1 == alt & A2 == ref ~ 0,
      A2 == alt & A1 == ref ~ 1,
      TRUE ~ NA
    ),
    BETA_final = case_when(
      flip == 0 ~ BETA,
      flip == 1 ~ -BETA,
      TRUE ~ NA
    ),
    A1_final = case_when(
      flip == 0 ~ A1,
      flip == 1 ~ A2,
      TRUE ~ NA
    ),
    A2_final = case_when(
      flip == 0 ~ A2,
      flip == 1 ~ A1,
      TRUE ~ NA
    )
  ) %>%
  transmute(CHR=CHR_38, BP=BP_38, rsID, SNP, A1 = A1_final, A2 = A2_final, BETA = BETA_final, SE, P,CHR_37, BP_37) %>%
  mutate(variant_id_38 = paste(CHR, BP, A2, A1, sep = ":")) %>%
  distinct(variant_id_38, .keep_all = TRUE) %>%
  mutate(variant_id_37 = paste(CHR_37, BP_37, A2, A1, sep = ":")) %>%
  mutate(credible_id_37 = paste(CHR_37, BP_37, A1, A2, sep = ":")) %>%
  select(-CHR_37, -BP_37)
  


write.table(merged_data,
            file = file.path(output_path, "full_sumstats_flipped.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# Print a message indicating that the files have been saved successfully
print(paste("All files have been saved successfully. Output folder:", output_path))
