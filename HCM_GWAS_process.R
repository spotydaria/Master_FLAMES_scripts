options(scipen = 999)
suppressPackageStartupMessages(library(optparse))  # Parsing command-line options
suppressPackageStartupMessages(library(data.table))  # Efficient data reading and manipulation
suppressPackageStartupMessages(library(dplyr))  # Data manipulation and transformation
suppressPackageStartupMessages(library(readr))  # Reading and writing data
suppressPackageStartupMessages(library(tidyr))  # Data tidying

print(paste("Libraries loaded successfully at", Sys.time()))

# Set up command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="data/raw/Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv", 
              help="Path to the input sumstats file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="data/sumst_processed/DCM_stats_LAVA_37_exclMYBPC3reg.txt",
              help="Path to the output file", metavar="file")
)

print("Command-line options set up")

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))
print(paste("Input file:", opt$file))
print(paste("Output file:", opt$output))

# Read the TSV file
print("Reading data from input file")
dat1 <- fread(opt$file)
print(paste("Data read successfully. Rows:", nrow(dat1), ", Columns:", ncol(dat1)))

# Filter data based on sample size
print("Filtering data based on sample size (n_samples >= 0.96 * max(n_samples))")
dat <- dat1
dat <- dat[dat$n_samples >= 0.96 * max(dat$n_samples), ]
print(paste("Data filtered by sample size. Remaining Rows:", nrow(dat)))

# Select and rename columns
print("Selecting relevant columns")
dat <- dat[, c("rsid", "chr", "pos", "effect_allele", "noneffect_allele", "eaf", "beta", "se", "pvalue")]

# Calculate minor allele frequency (maf)
print("Calculating minor allele frequency (maf)")
dat$maf <- dat$eaf
dat[dat$eaf > 0.5, 'maf'] <- 1 - dat[dat$eaf > 0.5, 'eaf']

# Calculate effective sample size (n)
print("Calculating effective sample size (n)")
dat$n <- 1 / (2 * dat$maf * (1 - dat$maf) * (dat$se^2))

# Calculate Z-scores
print("Calculating Z-scores")
dat$z <- dat$beta / dat$se

# Filter out specific region in chromosome 11 (MYBPC3 region)
print("Filtering out specific region in chromosome 11 (MYBPC3 region)")
dat3 <- dat[-(which(dat$chr == 11 & dat$pos >= 30000000 & dat$pos <= 80000000)), ]
print(paste("Data filtered for specific region. Remaining Rows:", nrow(dat3)))

# Further processing and tidying data for output
print("Tidying and processing data for output")
processed_data <- dat3 %>%
  transmute(
    rsID = rsid,
    CHR = chr, # Extract CHR
    BP = pos, # Extract BP
    SNP = paste0(CHR, ":", BP),
    A1 = effect_allele,
    A2 = noneffect_allele,
    NMISS = n,
    BETA = beta,
    SE = se,
    P = pvalue) %>%
  filter(!is.na(CHR))

print(paste0("Original sumst SNP count: ", nrow(dat1)))
print(paste0("Processed sumst SNP count: ", nrow(processed_data)))

# Write output
print("Writing processed data to output file")
write.table(processed_data, file=opt$output, col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
print(paste("Data written to output file successfully at", Sys.time()))
