# Load necessary libraries
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

monopogen_directory <- args[1]
output_directory <- args[2]
sample_id <- args[3]

filter_snv_matrices <- function(monopogen_directory, output_directory, sample_id) {
  # Initialize list to store matrices and chromosome identifiers
  snv_matrices <- list()

  # 1. Load cellID.filter.csv
  common_cells_path <- file.path(output_directory, paste0(sample_id, ".cellID.filter.csv"))
  if (!file.exists(common_cells_path)) {
    stop(paste("Common cells file not found:", common_cells_path))
  }
  
  common_cells_df <- read.csv(common_cells_path, stringsAsFactors = FALSE)
  
  # Check if required columns exist
  required_columns <- c("cell", "index")
  if (!all(required_columns %in% colnames(common_cells_df))) {
    stop(paste("Common cells file must contain columns:", paste(required_columns, collapse=", ")))
  }
  
  # Create a named vector: names are cell barcodes, values are indices
  cell_info <- setNames(common_cells_df$index, common_cells_df$cell)
  
  message(paste("Loaded", length(cell_info), "common cell barcodes from", paste0(sample_id, ".cellID.filter.csv")))
  
  # 2. Load SNVs.filter.csv
  filtered_snvs_path <- file.path(output_directory, paste0(sample_id, ".SNVs.filter.csv"))
  if (!file.exists(filtered_snvs_path)) {
    stop(paste("Filtered SNVs file not found:", filtered_snvs_path))
  }
  
  filtered_snvs_df <- read.csv(filtered_snvs_path, stringsAsFactors = FALSE)
  
  # Create SNV IDs
  filtered_snvs_df <- filtered_snvs_df %>%
    mutate(SNV_ID = paste0(chr, ":", pos, ":", Ref_allele, ":", Alt_allele))
  
  snv_ids <- filtered_snvs_df$SNV_ID
  
  message(paste("Loaded", length(snv_ids), "filtered SNVs from", paste0(sample_id, ".SNVs.filter.csv")))

  # 3. Iterate through chromosomes 1 to 22
  for (chr_num in 1:22) {
    # Construct the SNV_mat.RDS filename
    snv_mat_filename <- paste0("chr", chr_num, ".SNV_mat.RDS")
    snv_mat_path <- file.path(monopogen_directory, snv_mat_filename)
    
    if (!file.exists(snv_mat_path)) {
      warning(paste("SNV matrix file not found for chromosome", chr_num, ":", snv_mat_path))
      next
    }
    
    # Load the SNV matrix
    snv_mat <- readRDS(snv_mat_path)
    
    # Ensure rownames and colnames are set
    if (is.null(rownames(snv_mat)) || is.null(colnames(snv_mat))) {
      warning(paste("Matrix in", snv_mat_filename, "does not have rownames or colnames. Skipping."))
      next
    }
    
    # Subset rows to SNVs in filtered_snvs_df
    matching_snvs <- intersect(rownames(snv_mat), snv_ids)
    if (length(matching_snvs) == 0) {
      warning(paste("No matching SNVs found in", snv_mat_filename))
      next
    }
    
    snv_mat_filtered <- snv_mat[matching_snvs, , drop = FALSE]
    
    # Subset columns to common cell barcodes
    matching_cells <- intersect(colnames(snv_mat_filtered), names(cell_info))
    if (length(matching_cells) == 0) {
      warning(paste("No matching cell barcodes found in", snv_mat_filename))
      next
    }
    
    snv_mat_filtered <- snv_mat_filtered[, matching_cells, drop = FALSE]
    
    # Reorder columns based on common_cells_df order
    snv_mat_filtered <- snv_mat_filtered[, common_cells_df$cell[common_cells_df$cell %in% matching_cells], drop = FALSE]
    
    # Define output filename
    #output_mat_filename <- paste0(sample_id, ".chr", chr_num, ".SNV_mat.filter.RDS")
    #output_mat_path <- file.path(output_directory, output_mat_filename)
    
    # Save the filtered matrix
    #saveRDS(snv_mat_filtered, file = output_mat_path)
    
    #message(paste("Filtered SNV matrix saved to", output_mat_path))

    # Check if it's a matrix
    if (!is.matrix(snv_mat_filtered)) {
	    #warning(paste("File does not contain a matrix:", output_mat_filename))
	    snv_mat_filtered <- as.matrix(snv_mat_filtered)
    }

    # Store in the list
    chr_label <- paste0("chr", chr_num)
    snv_matrices[[chr_label]] <- snv_mat_filtered
    #message(paste("Loaded SNV matrix for", chr_label, "with", nrow(snv_mat_filtered), "variants and", ncol(snv_mat_filtered), "cells."))
  }
  
  # 4. Optionally, save the cell barcode and index information as a list or another structure
  # For simplicity, we return them as a named vector
  #return(list(cell_info = cell_info, filtered_snvs = filtered_snvs_df))
  return(snv_matrices)
}

merge_snv_matrices <- function(snv_matrices, output_directory, sample_id) {
  # Check if any matrices were loaded
  if (length(snv_matrices) == 0) {
    stop("No valid SNV matrices were loaded.")
  }

  # Verify that all matrices have the same column names (cell barcodes)
  cell_barcodes <- lapply(snv_matrices, colnames)
  unique_cell_barcodes <- unique(sapply(cell_barcodes, paste, collapse = ","))

  if (length(unique_cell_barcodes) != 1) {
    stop("Column names (cell barcodes) are not consistent across all SNV matrices.")
  }

  # All matrices have the same cell barcodes
  common_cell_barcodes <- cell_barcodes[[1]]
  message(paste("All SNV matrices have", length(common_cell_barcodes), "consistent cell barcodes."))

  # Merge all matrices by row (variants)
  message("Merging all SNV matrices into a single matrix...")
  merged_snv_mat <- do.call(rbind, snv_matrices)

  # Remove duplicate rows if any (though variants are expected to be unique)
  merged_snv_mat <- unique(merged_snv_mat)

  message(paste("Merged SNV matrix has", nrow(merged_snv_mat), "unique variants and", ncol(merged_snv_mat), "cells."))

  # Optionally, set rownames and colnames explicitly (if needed)
  # rownames(merged_snv_mat) <- unique(rownames(merged_snv_mat))
  # colnames(merged_snv_mat) <- common_cell_barcodes

  # Save the merged matrix as RDS
  merged_rds_filename <- paste0(sample_id, ".SNV_mat.RDS")
  merged_rds_path <- file.path(output_directory, merged_rds_filename)
  saveRDS(merged_snv_mat, file = merged_rds_path)
  message(paste("Merged SNV matrix saved as RDS:", merged_rds_path))

  # Save the merged matrix as CSV
  #merged_csv_filename <- paste0(sample_id, ".SNV_mat.csv")
  #merged_csv_path <- file.path(output_directory, merged_csv_filename)

  # To save large matrices as CSV, consider writing in chunks or using data.table for efficiency
  # Here, we'll convert to data frame for simplicity
  #merged_df <- as.data.frame(merged_snv_mat, stringsAsFactors = FALSE)
  #merged_df$Variant_ID <- rownames(merged_snv_mat)

  # Reorder columns to have Variant_ID first
  #merged_df <- merged_df %>% select(Variant_ID, everything())

  # Write to CSV
  #write.csv(merged_df, file = merged_csv_path, row.names = FALSE)
  #message(paste("Merged SNV matrix saved as CSV:", merged_csv_path))

  message("Merging process completed successfully.")

  # Return the merged matrix invisibly
  invisible(merged_snv_mat)
}


# --------------------------------------
# Usage
# --------------------------------------
snv_matrices = filter_snv_matrices(monopogen_directory, output_directory, sample_id)
merge_snv_matrices(snv_matrices, output_directory, sample_id)
