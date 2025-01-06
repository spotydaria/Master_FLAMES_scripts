suppressPackageStartupMessages(library(optparse))  # Parsing command-line options
suppressPackageStartupMessages(library(data.table))  # Efficient data reading and manipulation

# Set up command-line options
option_list <- list(
  make_option(c("-f", "--folder_path"), type="character", default=NULL, 
              help="Path to the folder containing the files", metavar="folder_path")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))
folder_path <- opt$folder_path
print(paste("Reading files from folder:", folder_path))

# Define file paths using the path provided in the bash script
files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# Read the files, extract rsID and SNP, and rename SNP column to specify file origin
data_tables <- lapply(files, function(file) {
  print(paste("Reading file:", file))
  dt <- fread(file, select = c("rsID", "SNP"))
  setnames(dt, "SNP", paste0("SNP_", gsub("data/sumst_processed/|_exclMYBPC3reg.*", "", basename(file))))
  return(dt)
})

print("Files read successfully")

# Perform iterative inner joins on the list of data.tables
print("Performing iterative inner joins on the data tables")
result_dt <- Reduce(function(x, y) merge(x, y, by = "rsID", all = FALSE), data_tables)
print("Inner joins completed")

# Additional cleaning to remove non-visible characters and ensure data types
print("Cleaning data to remove non-visible characters and ensure consistent data types")
result_dt[, (colnames(result_dt)[-1]) := lapply(.SD, function(x) gsub("[[:cntrl:]]", "", trimws(x))), .SDcols = patterns("^SNP_")]
result_dt[, (colnames(result_dt)[-1]) := lapply(.SD, as.character), .SDcols = patterns("^SNP_")]

# Recalculate SNP_equal with detailed mismatch information for debugging
print("Calculating SNP_equal and mismatch details")
result_dt[, SNP_equal := Reduce(`==`, .SD), .SDcols = patterns("^SNP_")]
result_dt[, SNP_mismatch_detail := do.call(paste, c(.SD, sep = " | ")), .SDcols = patterns("^SNP_")]

# Print results, including mismatch details for non-matching rows
print("Displaying first 10 rows of the result")
print(head(result_dt, 10))
print("Displaying mismatch details for non-matching rows")
print(head(result_dt[SNP_equal == FALSE, .(rsID, SNP_mismatch_detail)]))
