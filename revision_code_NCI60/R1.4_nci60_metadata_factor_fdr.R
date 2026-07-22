#!/usr/bin/env Rscript

# R1.4: Benjamini-Hochberg correction for NCI60 metadata-factor associations.
#
# This is the NCI60-only extraction of the revision analysis. It requires only
# the NCI60 MOFA factor matrix and metadata stored in the destination repository.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE))
script_dir <- if (length(script_file) > 0L) {
  dirname(normalizePath(script_file[[1]], mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
source(file.path(script_dir, "revision_helpers.R"))

repo_dir <- revision_repo_dir(
  package_dir = script_dir,
  env_name = "NCI60_REPO_DIR",
  required_paths = c(
    "results/mofa/z_matrix.RData",
    "support/all_metadata.RData"
  )
)
out_dir <- revision_output_dir(script_dir)

write_csv_base <- function(x, path) {
  write.csv(x, path, row.names = FALSE, quote = TRUE)
}

run_nci60_metadata_fdr <- function(repo_dir, out_dir) {
  objects <- new.env(parent = emptyenv())
  load(file.path(repo_dir, "results", "mofa", "z_matrix.RData"), envir = objects)
  load(file.path(repo_dir, "support", "all_metadata.RData"), envir = objects)

  if (!exists("Z_matrix", envir = objects, inherits = FALSE)) {
    stop("results/mofa/z_matrix.RData does not contain Z_matrix.")
  }
  if (!exists("all_metadata", envir = objects, inherits = FALSE)) {
    stop("support/all_metadata.RData does not contain all_metadata.")
  }

  Z_matrix <- objects$Z_matrix
  all_metadata <- objects$all_metadata
  if (!all(c("source", "target") %in% names(all_metadata))) {
    stop("all_metadata must contain source and target columns.")
  }

  factor_names <- colnames(Z_matrix)
  cell_lines <- rownames(Z_matrix)
  metadata_sources <- sort(unique(all_metadata$source))

  rows <- list()
  k <- 1L
  for (source_name in metadata_sources) {
    targets <- unique(all_metadata$target[all_metadata$source == source_name])
    x <- as.integer(cell_lines %in% targets)
    n_in_category <- sum(x == 1L)
    n_out_category <- sum(x == 0L)
    if (n_in_category < 3L || n_out_category < 3L) {
      next
    }

    for (factor_name in factor_names) {
      y <- Z_matrix[[factor_name]]
      fit <- lm(y ~ x)
      coefficient_table <- summary(fit)$coefficients
      rows[[k]] <- data.frame(
        metadata_category = source_name,
        factor = factor_name,
        n_in_category = n_in_category,
        n_out_category = n_out_category,
        score = unname(coefficient_table["x", "t value"]),
        p_value = unname(coefficient_table["x", "Pr(>|t|)"]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  if (length(rows) == 0L) {
    stop("No NCI60 metadata categories had at least three samples in and out of category.")
  }

  result <- do.call(rbind, rows)
  result$fdr_bh <- p.adjust(result$p_value, method = "BH")
  result$abs_score_gt_2 <- abs(result$score) > 2
  result$p_lt_0_05 <- result$p_value < 0.05
  result$fdr_lt_0_05 <- result$fdr_bh < 0.05
  result$fdr_lt_0_10 <- result$fdr_bh < 0.10
  result <- result[order(result$fdr_bh, result$p_value), , drop = FALSE]
  rownames(result) <- NULL

  write_csv_base(result, file.path(out_dir, "R1.4_nci60_metadata_factor_fdr.csv"))

  tissue_result <- result[grepl("^tissue:", result$metadata_category), , drop = FALSE]
  write_csv_base(tissue_result, file.path(out_dir, "R1.4_nci60_tissue_factor_fdr.csv"))

  result
}

nci60_result <- run_nci60_metadata_fdr(repo_dir, out_dir)
summary_lines <- c(
  "# R1.4 NCI60 Metadata-Factor Multiple-Testing Analysis",
  "",
  sprintf("- Metadata categories tested: %d", length(unique(nci60_result$metadata_category))),
  sprintf("- MOFA factors tested: %d", length(unique(nci60_result$factor))),
  sprintf("- Total category-factor tests: %d", nrow(nci60_result)),
  sprintf("- Tests with raw p < 0.05: %d", sum(nci60_result$p_lt_0_05, na.rm = TRUE)),
  sprintf("- Tests with BH FDR < 0.05: %d", sum(nci60_result$fdr_lt_0_05, na.rm = TRUE)),
  sprintf("- Tests with BH FDR < 0.10: %d", sum(nci60_result$fdr_lt_0_10, na.rm = TRUE)),
  sprintf("- Tests with abs(score) > 2: %d", sum(nci60_result$abs_score_gt_2, na.rm = TRUE)),
  "",
  "Top associations by BH FDR:",
  paste(
    utils::capture.output(print(head(
      nci60_result[, c("metadata_category", "factor", "n_in_category", "score", "p_value", "fdr_bh")],
      15
    ), row.names = FALSE)),
    collapse = "\n"
  )
)
writeLines(summary_lines, file.path(out_dir, "R1.4_nci60_multiple_testing_summary.md"))
cat(paste(summary_lines, collapse = "\n"), "\n")
