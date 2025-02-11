# Load necessary libraries
library(progress)
library(parallel)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

snv_matrix <- args[1]
cell_barcodes <- args[2]
output_directory <- args[3]
sample_id <- args[4]

filter_snv_data <- function(
  snv_mat, 
  cell_barcodes
) {
  snv_mat <- readRDS(snv_mat)
  cell_barcodes <- read.delim(cell_barcodes)
  keep_barcodes <- cell_barcodes$x
  keep_barcodes_in_snv <- intersect(keep_barcodes, colnames(snv_mat))
  snv_mat_filtered <- snv_mat[, keep_barcodes_in_snv, drop = FALSE]
  return(snv_mat_filtered)
}

convert_snv_matrix <- function(
  snv_mat,
  nCores = 2,
  chunkSize = 1000
) {
  # 1) Parse rownames (e.g., "chr1:184413:C:A") into chromosome, position, REF, ALT
  row_splits <- strsplit(rownames(snv_mat), ":")

  chrom_vec <- sapply(row_splits, `[`, 1)
  pos_vec   <- as.integer(sapply(row_splits, `[`, 2))
  ref_vec   <- sapply(row_splits, `[`, 3)
  alt_vec   <- sapply(row_splits, `[`, 4)

  # Create Info dataframe
  Info <- data.frame(
    `#CHROM` = chrom_vec,
    POS      = pos_vec,
    REF      = ref_vec,
    ALT      = alt_vec,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # 2) Prepare empty matrices X, N, Z
  nVariants <- nrow(snv_mat)
  nCells    <- ncol(snv_mat)

  X <- matrix(0L, nrow = nVariants, ncol = nCells)
  N <- matrix(0L, nrow = nVariants, ncol = nCells)
  Z <- matrix(0,  nrow = nVariants, ncol = nCells)

  rownames(X) <- rownames(snv_mat)
  colnames(X) <- colnames(snv_mat)
  rownames(N) <- rownames(snv_mat)
  colnames(N) <- colnames(snv_mat)
  rownames(Z) <- rownames(snv_mat)
  colnames(Z) <- colnames(snv_mat)

  # 3) Fill in X (ALT counts), N (REF+ALT counts), Z (0 or 1) by parsing strings
  # Function to parse one row of character values
  parse_one_row <- function(row_values) {
    # row_values: character vector of length = nCells, e.g. c("14/0","0/0","4/0",...)
    # We'll return a small list with Xrow, Nrow, Zrow
    nCells_local <- length(row_values)
    
    Xrow <- integer(nCells_local)
    Nrow <- integer(nCells_local)
    Zrow <- numeric(nCells_local)
    
    for (j in seq_len(nCells_local)) {
      val <- row_values[j]
      if (grepl("/", val, fixed = TRUE)) {
        parts <- strsplit(val, "/")[[1]]
        ref_count <- suppressWarnings(as.numeric(parts[1]))
        alt_count <- suppressWarnings(as.numeric(parts[2]))
        
        if (!is.na(ref_count) && !is.na(alt_count)) {
          Xrow[j] <- alt_count
          Nrow[j] <- ref_count + alt_count
          Zrow[j] <- if (alt_count > 0) 1 else 0
        } else {
          # If parsing fails => 0
          Xrow[j] <- 0
          Nrow[j] <- 0
          Zrow[j] <- 0
        }
      } else {
        # Invalid pattern => 0
        Xrow[j] <- 0
        Nrow[j] <- 0
        Zrow[j] <- 0
      }
    }
    return(list(Xrow = Xrow, Nrow = Nrow, Zrow = Zrow))
  }
  
  # Create a progress bar for the total number of rows (nVariants)
  pb <- progress_bar$new(
    format = "Processing rows [:bar] :percent ETA: :eta",
    total = nVariants,
    clear = FALSE,
    width = 60
  )
  
  # Split row indices into chunks
  row_indices <- seq_len(nVariants)
  chunk_list  <- split(row_indices, ceiling(row_indices / chunkSize))
  
  for (chunk_id in seq_along(chunk_list)) {
    # The row indices for this chunk
    rows_in_chunk <- chunk_list[[chunk_id]]
    
    # Convert each row in this chunk to a character vector
    row_char_list <- lapply(rows_in_chunk, function(i) {
      as.character(snv_mat[i, ])
    })
    
    # Parallel apply for each row in the chunk
    result_list <- mclapply(
      X        = row_char_list,
      FUN      = parse_one_row,
      mc.cores = nCores
    )
    
    # Fill the corresponding positions in X, N, Z
    for (idx in seq_along(rows_in_chunk)) {
      i <- rows_in_chunk[idx]
      X[i, ] <- result_list[[idx]]$Xrow
      N[i, ] <- result_list[[idx]]$Nrow
      Z[i, ] <- result_list[[idx]]$Zrow
    }
    
    # Update the progress bar by the number of processed rows in this chunk
    pb$tick(length(rows_in_chunk))
  }

  # 4) Return a list in the desired format
  return(list(
	      X = X,
	      N = N,
	      Z = Z,
	      Info = Info
  ))
}

PlotCellMutationDist <- function(X, N, Z, output_directory) {
  # Open a PDF device
  pdf(paste0(output_directory, "CellMutationDist.pdf"), width = 13, height = 4)

  # Set up a layout for four plots side by side
  par(mfrow = c(1, 4))

  # 1) Histogram of rowSums(Z)
  hist(
    rowSums(Z, na.rm = TRUE),
    breaks = 30,
    main = "Cell distribution \n for each mutation call",
    xlab = "Number of cells"
  )

  # 2) Same rowSums(Z), but with x-axis limit
  hist(
    pmin(rowSums(Z, na.rm = TRUE), 50),
    breaks = seq(0, 50, by = 1),
    main = "Cell distribution \n for each mutation call",
    xlab = "Number of cells",
    xlim = c(0, 50)
  )

  # 3) Histogram of colSums(Z)
  hist(
    colSums(Z, na.rm = TRUE),
    breaks = 10,
    main = "Mutation distribution \n for each cell",
    xlab = "Number of mutations"
  )

  # 4) Histogram of colSums(Z)
  hist(
    pmin(colSums(Z, na.rm = TRUE), 30),
    breaks = seq(0, 30, by = 1),
    main = "Mutation distribution \n for each cell",
    xlab = "Number of mutations",
    xlim = c(0, 30)
  )

  # Close the PDF device
  dev.off()

  message("Histogram plots saved to: ", paste0(output_directory, "CellMutationDist.pdf"))
}

CheckFilterCutoffs <- function(X, N, Z, cut.off.mut = 5, cut.off.cell = 5) {
  # Step 1: Determine variant filtering criteria based on cut.off.mut
  # Create a logical vector indicating which variants meet the mutation cutoff
  filtered_mut <- rowSums(Z, na.rm = TRUE) >= cut.off.mut

  cat("If filtering were applied: ", sum(filtered_mut), " out of ", nrow(Z), " variants would be retained.\n")
  cat("Excluding ", nrow(Z) - sum(filtered_mut), " variants that do not meet the mutation cutoff.\n\n")

  # Subset the matrix to include only the variants that passed the mutation filtering
  Z_mut_filtered <- Z[filtered_mut, , drop = FALSE]

  # Step 2: Determine cell filtering criteria based on cut.off.cell using only the filtered variants
  filtered_cell <- colSums(Z_mut_filtered, na.rm = TRUE) >= cut.off.cell
  
  cat("If filtering were applied: ", sum(filtered_cell), " out of ", ncol(Z), " cells would be retained (using filtered variants only).\n")
  cat("Excluding ", ncol(Z) - sum(filtered_cell), " cells that do not meet the cell cutoff.\n")

  # Return the names of filtered variants (rows) and filtered cells (columns)
  return(list(filtered_variants = rownames(Z)[filtered_mut],
	      filtered_barcodes = colnames(Z)[filtered_cell])
  )
}

FilterCellMutation <- function(snv_mat, cell_barcodes, somatic_variants, output_directory, sample_id) {
  snv_mat_filtered <- snv_mat[, cell_barcodes, drop = FALSE]
  snv_mat_filtered <- snv_mat_filtered[somatic_variants, , drop = FALSE]
  
  # Save the filtered matrix as CSV
  filtered_csv_filename <- paste0(sample_id, ".SNV_mat.filter.csv")
  filtered_csv_path <- file.path(output_directory, filtered_csv_filename)

  # Here, we'll convert to data frame for simplicity
  filtered_df <- as.data.frame(snv_mat_filtered, stringsAsFactors = FALSE)
  filtered_df$Variant_ID <- rownames(snv_mat_filtered)

  # Reorder columns to have Variant_ID first
  filtered_df <- filtered_df %>% select(Variant_ID, everything())

  # Write to CSV
  write.csv(filtered_df, file = filtered_csv_path, row.names = FALSE)
  message(paste("Filtered SNV matrix saved as CSV:", filtered_csv_path))
}


# --------------------------------------
# Usage
# --------------------------------------
snv_matrix_filtered <- filter_snv_data(
  snv_mat = snv_matrix,			# SNV matrix
  cell_barcodes = cell_barcodes		# Cell barcodes
  )
result <- convert_snv_matrix(
  snv_mat = snv_matrix_filtered,
  nCores = 4,
  chunkSize = 1000
  )
PlotCellMutationDist(
  X = result$X,
  N = result$N,
  Z = result$Z,
  output_directory = output_directory
  )
output <- CheckFilterCutoffs(
  X = result$X,
  N = result$N,
  Z = result$Z,
  cut.off.mut = 5,
  cut.off.cell = 1
  )
FilterCellMutation(
  snv_mat = snv_matrix_filtered,
  cell_barcodes = output$filtered_barcodes,
  somatic_variants = output$filtered_variants,
  output_directory = output_directory,
  sample_id = sample_id
  )
