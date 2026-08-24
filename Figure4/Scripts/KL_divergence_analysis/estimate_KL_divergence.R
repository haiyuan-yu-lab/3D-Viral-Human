# Get the probability distribution for each amino acid for each position 
# Estimate KL divergence using the reference distribution
# Load libraries
library(Biostrings)
library(seqinr)
library(stringr)
library(philentropy)
# Read the protein alignment file (in FASTA format)
args <- commandArgs(trailingOnly = TRUE)
vir <- args[1]
parent_dir <- args[2]
alignment_dir <- file.path('path/to/MSAs') #MSAs can be requested from the authors
out_dir <-  file.path("/path/to/output") #output directory for KL divergence results
dir.create(out_dir, recursive = TRUE)
alignment_file_list <- list.dirs(alignment_dir, recursive = FALSE)
gene_list <- alignment_file_list %>%
  basename() %>%
  str_remove("_prot_sequences\\.fa$")
for (f in seq_along(alignment_file_list)) {
  alignment_file_dir <- alignment_file_list[f]
  print(paste0("Now processing ", gene_list[f]))
  alignment_files <- list.files(alignment_file_dir)
  for (file in seq_along(alignment_files)) {
    alignment <- seqinr::read.alignment(file = paste0(alignment_file_dir, "/",
                                                      alignment_files[file]),
                                        format = "fasta",
                                        forceToLower = FALSE)
    #create matrix of the alignment where each column in an amino acid
    alignment_matrix <- seqinr::as.matrix.alignment(alignment)
    #dropping sequences with more than 30% of gaps
    gap_percentages <- rowMeans(alignment_matrix == "-", na.rm = TRUE) * 100
    x_perc <- colMeans(alignment_matrix == "x", na.rm = TRUE) * 100
    rows_to_keep <- gap_percentages <= 30
    alignment_matrix <- alignment_matrix[rows_to_keep,]
    #get the frequency matrix for all positions combined
    amino_acids <- c("m", "w", "f", "y", "h", "q", "n", "k",
                     "d", "e", "c", "i", "t", "v", "p",
                     "a", "g", "s", "l", "r")
    flattened_vector <- as.vector(alignment_matrix)
    #not counting the gaps and 'x' in the total
    #count only the allowed amino-acids (levels) once
    #reference counts
    idx <- match(flattened_vector, amino_acids)
    counts <- tabulate(idx[!is.na(idx)], nbins = length(amino_acids))
    if (sum(counts) == 0) {
      expected_freq <- setNames(rep(0, length(counts)), amino_acids)
    } else {
      expected_freq <- counts / sum(counts)
      names(expected_freq) <- amino_acids
    }
    #per-position counts (rows = amino_acids, cols = positions)
    frequency_distribution <- apply(alignment_matrix, 2, function(col) {
      tabulate(match(col, amino_acids), nbins = length(amino_acids))
    })
    if (is.null(dim(frequency_distribution))) {
      frequency_distribution <- matrix(frequency_distribution,
                                       nrow = length(amino_acids))
    }
    rownames(frequency_distribution) <- amino_acids
    col_sums <- colSums(frequency_distribution)
    probability_distribution <- sweep(frequency_distribution, 2, col_sums, "/")
    #keep empty positions as NA
    probability_distribution[, col_sums == 0] <- NA
    kl_divergence <- data.frame("POS" = numeric(), "KL" = numeric())

    for (i in seq_len(ncol(probability_distribution))) {
      #create a matrix with the expected frequencies and observed probabilities
      prob_matrix <- rbind(probability_distribution[, i], expected_freq)
      if (any(is.na(probability_distribution[, i]))) {
        #if position has no observed amino acids, assign NA to KL divergence
        kl_divergence <- rbind(kl_divergence, data.frame("POS" = i,
                                                         "KL" = NA_real_))
      } else {
        # Calculate KL divergence for the current position and add to dataframe
        kl_divergence <- rbind(kl_divergence,
                               data.frame("POS" = i,
                                          "KL" = KL(prob_matrix,
                                                    test.na = TRUE,
                                                    unit = "log10",
                                                    est.prob = NULL,
                                                    epsilon = 1e-05)))
      }
    }
    out_file_name <- paste0(gene_list[f], "_", alignment_files[file], ".txt")
    write.table(kl_divergence, file = paste0(out_dir, out_file_name),
                col.names = TRUE, row.names = FALSE)
  }
}
