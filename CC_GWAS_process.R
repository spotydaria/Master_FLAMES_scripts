suppressPackageStartupMessages(library(optparse))  # Parsing command-line options
suppressPackageStartupMessages(library(data.table))  # Efficient data reading and manipulation
suppressPackageStartupMessages(library(dplyr))  # Data manipulation and transformation
suppressPackageStartupMessages(library(readr))  # Reading and writing data
suppressPackageStartupMessages(library(tidyr))  # Data tidying

print("Libraries loaded successfully")

# Set up command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="data/raw/CC_GWAS__DCM__HCM_exclMYBPC3reg.txt", 
              help="Path to the input sumstats file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="data/sumst_processed/CC_GWAS_sumst_37_exclMYBPC3reg.txt",
              help="Path to the output file", metavar="file")
)

print("Command-line options set up")

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))
print(paste("Input file:", opt$file))
print(paste("Output file:", opt$output))
print(paste("Libraries loaded successfully at", Sys.time()))

# Read the TSV file
print("Reading data from input file")
options(scipen = 999)
# Read data
print("Reading data from input file")
dat1 <- fread(opt$file)
dat <- dat1
print(paste("Data read successfully. Rows:", nrow(dat1), ", Columns:", ncol(dat1)))

# 37 build
print("Processing data: selecting and renaming columns")
processed_data <- dat1 %>%
  transmute(
    rsID = SNP,
    CHR, # Extract CHR from CHRBP_B37
    BP, # Extract BP from CHRBP_B37
    SNP = paste0(CHR, ":", BP),
    A1 = EA,
    A2 = NEA,
    BETA = OLS_beta,
    SE = OLS_se,
    P = OLS_pval) %>%
  filter(!is.na(CHR))

print(paste0("Original sumst SNP count: ", nrow(dat1)))
print(paste0("Processed sumst SNP count: ", nrow(processed_data)))

# Write output
print("Writing processed data to output file")
write.table(processed_data, file=opt$output, col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
print(paste("Data written to output file successfully at", Sys.time()))