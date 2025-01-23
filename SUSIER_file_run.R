# Load necessary libraries
suppressPackageStartupMessages(library(data.table))  # Efficient data manipulation
suppressPackageStartupMessages(library(dplyr))       # Data manipulation and transformation
suppressPackageStartupMessages(library(readr))       # Data input and output
suppressPackageStartupMessages(library(tidyr))       # Data tidying
suppressPackageStartupMessages(library(susieR))      # Implements SuSiE (Sum of Single Effects)
suppressPackageStartupMessages(library(stringr))     # String manipulation
suppressPackageStartupMessages(library(Rfast))       # Fast data manipulation

# Parse arguments for main script
args <- commandArgs(trailingOnly = TRUE)
# args <-c( "/home/dkramarenk/projects/LAVA/DCM_HCM/data/SUSIE/HCM_MTAG/", "/home/dkramarenk/projects/LAVA/DCM_HCM/data/SUSIE/HCM_MTAG/out/", "35512", "10")

# Check that correct number of arguments is provided
if (length(args) != 4) {
  stop("Usage: Rscript in_file_path output_path N(susie) L(susie)")
}

# Assign arguments to variables

in_file_path <- args[1]
output_path <- args[2]
N_CC <- as.numeric(args[3])
susie_l <- as.numeric(args[4])

print(paste("Input file path:", in_file_path))
print(paste("Output path:", output_path))
print(paste("N_CC (sample size for SuSiE):", N_CC))
print(paste("L (number of effect components for SuSiE):", susie_l))

# Load summary statistics and summarized files
tryCatch({
  final_sumst <- fread(file.path(in_file_path, "full_sumstats_flipped.txt"))
  print("Successfully loaded full_sumstats_flipped.txt")
}, error = function(e) {
  stop("Error loading full_sumstats_flipped.txt: ", e$message)
})

tryCatch({
  summarized_files <- fread(file.path(in_file_path, "regions_to_analyze.txt")) %>%
    filter(!is.na(file))   %>%
    mutate(susier = NA, min_abs_corr = NA, cs_n = NA)
  print("Successfully loaded regions_to_analyze.txt")
}, error = function(e) {
  stop("Error loading regions_to_analyze.txt: ", e$message)
})

# Initialize the index file for FLAMES
index_file_path <- file.path(output_path, "FLAMES_index.ind")
tryCatch({
  write.table(data.frame(Filename = character(), Annotfiles = character()), 
              index_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  print(paste("Successfully initialized FLAMES index file at:", index_file_path))
}, error = function(e) {
  stop("Error initializing FLAMES index file: ", e$message)
})

for (i in 1:nrow(summarized_files)) {
  
  # Extract information for the current iteration
  chr_LD <- summarized_files$chr_LD[i]
  start_LD <- summarized_files$start_LD[i]
  end_LD <- summarized_files$end_LD[i]
  
  # Create the subset of summary statistics
  subset_final_sumst_1 <- final_sumst %>%
    filter(CHR == chr_LD & BP >= start_LD & BP <= end_LD) 
  
  # Read the Rvar file
  rvar_file <- summarized_files$file[i]
  rvar_data <- fread(rvar_file)
  
  # Replace the extension to read the RDS matrix file
  rds_file <- sub("\\.Rvar$", ".RDS", rvar_file)
  ld_matrix <- readRDS(rds_file)
    
  print(paste0("LD matrix file:", rds_file))

  rvar_data <- rvar_data %>%
    mutate(variant_id_38 = paste(chrom, pos, ref, alt, sep = ":"))
  
  # Set rownames and colnames in LD matrix
  rownames(ld_matrix) <- colnames(ld_matrix) <- rvar_data$variant_id_38
  
  # Find matching SNPs between sumstats and LD matrix
  common_variants <- intersect(subset_final_sumst_1$variant_id_38, rownames(ld_matrix))
  
  # Subset LD matrix and summary statistics to match common variants
  ld_matrix_subset <- ld_matrix[common_variants, common_variants]
  
  subset_final_sumst_1 <- subset_final_sumst_1 %>%
    filter(variant_id_38 %in% common_variants) %>%
    arrange(match(variant_id_38, rownames(ld_matrix_subset)))
  
  # Check for mismatches
  mismatches <- !(subset_final_sumst_1$variant_id_38 == rownames(ld_matrix_subset))
  
  if (any(mismatches)) {
    print(paste0("Mismatches found in iteration ", i, "/", nrow(summarized_files)))
    next
  } else {
    print(paste0("No mismatches, all variants align correctly in iteration ", i, "/", nrow(summarized_files)))
  }

  # Perform SuSiE fine-mapping
  min_abs_corr <- 0.5
  fitted_rss <- susie_rss(bhat = subset_final_sumst_1$BETA, 
                          shat = subset_final_sumst_1$SE, 
                          R = ld_matrix_subset, 
                          n = N_CC, L = susie_l,
                          min_abs_corr = min_abs_corr)
  fitted_summary <- summary(fitted_rss)

  credible_set_original <- as.data.frame(fitted_summary$vars) %>%
    mutate(variant_id_37= sapply(variable, function(x) subset_final_sumst_1[x, ]$variant_id_37)) %>%
    mutate(variant_id_38 = sapply(variable, function(x) subset_final_sumst_1[x, ]$variant_id_38)) %>%
    mutate(credible_id_37 = sapply(variable, function(x) subset_final_sumst_1[x, ]$credible_id_37)) %>%
    arrange(cs)


  credible_set_original_filt <- credible_set_original %>%
      filter(cs > 0)
       
  if (nrow(credible_set_original_filt) > 0) {
    rm(credible_set_original_filt)
    print(paste0("Success: min_abs_corr = 0.5"))
    credible_set_flames <- credible_set_original %>%
      filter(cs > 0) %>%
      select(cred1 = credible_id_37, prob1 = variable_prob)

    credible_set_cs <- as.data.frame(fitted_summary$cs) %>%
      arrange(cs)

    summarized_files$susier[i] <- 1
    summarized_files$min_abs_corr[i] <- min_abs_corr
    summarized_files$cs_n[i] <- max(as.numeric(credible_set_original$cs), na.rm = TRUE)


    unique_cs <- unique(credible_set_cs$cs)

    # Save results if there are credible sets
    # Save multiple files for credible sets if cs > 1
    for (j in seq_along(unique_cs)) {
      locus_name = paste0("locus_", i, ".", j, "_min_abs_corr_", min_abs_corr)
      print(paste0("Processing ", locus_name))
      cs_subset <- credible_set_original %>% filter(cs == unique_cs[j])
      output_filename <- paste0(output_path, "orig_loc/", locus_name, ".txt")
      annot_filename <- paste0(output_path, "FLAMES_annotated_", locus_name, ".txt")
      write.table(cs_subset,
                  file = output_filename,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
                                    print(paste("Saving credible_set_original file to:", output_filename))

      output_filename2 <- paste0(output_path, locus_name, ".txt")
      cs_subset_sus <- cs_subset %>%
      select(cred1 = credible_id_37, prob1 = variable_prob)
            write.table(cs_subset_sus,
                  file = output_filename2,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
                        print(paste("Saving cs_subset file to:", output_filename2))

      # Append each file to the index file for FLAMES
      index_entry <- data.frame(Filename = paste0(locus_name, ".txt"), Annotfiles = paste0("FLAMES_annotated_",locus_name,".txt"))
          write.table(index_entry, index_file_path, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE)
          write.table(credible_set_cs,
              file = paste0(output_path, "cs_stats_locus_", i, ".txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
          print(paste("Saving annotated file to:", paste0(output_path, "cs_stats_locus_", i, ".txt")))
    }

  } else {
    print(paste0("Failed: min_abs_corr = 0.5"))
    min_abs_corr <- 0.25
    fitted_rss <- susie_rss(bhat = subset_final_sumst_1$BETA, 
                            shat = subset_final_sumst_1$SE, 
                            R = ld_matrix_subset, 
                            n = N_CC, L = susie_l,
                            min_abs_corr = min_abs_corr)
    fitted_summary <- summary(fitted_rss)

    credible_set_original <- as.data.frame(fitted_summary$vars) %>%
      mutate(variant_id_37= sapply(variable, function(x) subset_final_sumst_1[x, ]$variant_id_37)) %>%
      mutate(variant_id_38 = sapply(variable, function(x) subset_final_sumst_1[x, ]$variant_id_38)) %>%
      mutate(credible_id_37 = sapply(variable, function(x) subset_final_sumst_1[x, ]$credible_id_37)) %>%
      arrange(cs)

    credible_set_original_filt <- credible_set_original %>%
    filter(cs > 0)

    if (nrow(credible_set_original_filt) > 0) {
      rm(credible_set_original_filt)
      print(paste0("Success: min_abs_corr = 0.25"))

      credible_set_flames <- credible_set_original %>%
        filter(cs > 0) %>%
        select(cred1 = credible_id_37, prob1 = variable_prob)

      credible_set_cs <- as.data.frame(fitted_summary$cs) %>%
        arrange(cs)

      summarized_files$susier[i] <- 1
      summarized_files$min_abs_corr[i] <- min_abs_corr
      summarized_files$cs_n[i] <- max(as.numeric(credible_set_original$cs), na.rm = TRUE)


      unique_cs <- unique(credible_set_cs$cs)

      # Save results if there are credible sets
      # Save multiple files for credible sets if cs > 1
      for (j in seq_along(unique_cs)) {
        locus_name = paste0("locus_", i, ".", j, "_min_abs_corr_", min_abs_corr)
        cs_subset <- credible_set_original %>% filter(cs == unique_cs[j])
        output_filename <- paste0(output_path, "/orig_loc/locus_", i, ".", j, ".txt")
        annot_filename <- paste0(output_path, "/FLAMES_annotated_locus_", i, ".", j, "_min_abs_corr", min_abs_corr, ".txt")
        write.table(cs_subset,
                    file = output_filename,
                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
          
        output_filename2 <- paste0(output_path, locus_name, ".txt")
        cs_subset_sus <- cs_subset %>%
        select(cred1 = credible_id_37, prob1 = variable_prob)
            write.table(cs_subset_sus,
                  file = output_filename2,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
                        print(paste("Saving cs_subset file to:", output_filename2))

        # Append each file to the index file for FLAMES
        index_entry <- data.frame(Filename = paste0(locus_name, ".txt"), Annotfiles = paste0("FLAMES_annotated_",locus_name,".txt"))
          write.table(index_entry, index_file_path, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE)
          write.table(credible_set_cs,
              file = paste0(output_path, "cs_stats_locus_", i, ".txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
          print(paste("Saving annotated file to:", paste0(output_path, "cs_stats_locus_", i, ".txt")))
      }
    
    } else {
    # make pseudo set
    subset_final_sumst_2 <- subset_final_sumst_1 %>%
    filter(P < 5e-8) %>%  # Filter rows where P is less than 5e-8
    arrange(P) 

      if (nrow(subset_final_sumst_2) > 0) {
        subset_final_sumst_3 <- subset_final_sumst_2 %>%        # Sort by P in ascending order
          slice_head(n = 1) 

        locus_name = paste0("locus_", i, ".","pseudo")
        print(paste0("Processing ", locus_name))

        output_filename2 <- paste0(output_path, locus_name, ".txt")
        cs_subset_sus <- subset_final_sumst_3 %>%
        select(cred1 = credible_id_37, prob1 = 1)
                write.table(cs_subset_sus,
                      file = output_filename2,
                      row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
                            print(paste("Saving cs_subset file to:", output_filename2))
        index_entry <- data.frame(Filename = paste0(locus_name, ".txt"), Annotfiles = paste0("FLAMES_annotated_",locus_name,".txt"))
        write.table(index_entry, index_file_path, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE)
                   
        summarized_files$susier[i] <- 2
        summarized_files$min_abs_corr[i] <- NA
        summarized_files$cs_n[i] <- 1
        print("Pseudo set")

      } else {

        summarized_files$susier[i] <- 0
        summarized_files$min_abs_corr[i] <- NA
        summarized_files$cs_n[i] <-NA

      print("No credible sets")
      }
      # Save credible set summary if available
    }
  }
}
# Print a message indicating that the files have been saved successfully
write.table(summarized_files, file = file.path(output_path, "regions_finemapped.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
print("All files have been saved successfully.")
