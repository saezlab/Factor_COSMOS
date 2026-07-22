#!/usr/bin/env Rscript

# R1.5: NCI60 MOON candidate-threshold sensitivity.
# The analysis reruns the TF-to-ligand branch at three candidate thresholds and
# recombines it with the stored receptor-to-TF/metabolite branch.

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
    "results/mofa/mofa_weights.RData",
    "data/cosmos/ligrec_TF_moon_inputs.Rdata",
    "support/dorothea_df.RData",
    "data/RNA/RNA_log2_FPKM_clean.csv",
    "results/cosmos/moon/moon_res_rec_to_TFmet.csv",
    "results/cosmos/moon/meta_network_filtered.csv",
    "support/c2.cp.v2022.1.Hs.symbols.gmt"
  )
)
out_dir <- revision_output_dir(script_dir)
revision_add_local_library(repo_dir)

suppressPackageStartupMessages({
  library(cosmosR)
  library(readr)
  library(dplyr)
})

thresholds <- c(1.5, 2, 2.5)
baseline_threshold <- 2

translate_column_hmdb <- function(my_column, hmdb_mapper_vec) {
  vapply(my_column, function(x) {
    x <- gsub("Metab__", "", x)
    x <- gsub("^Gene", "Enzyme", x)
    suffix <- stringr::str_extract(x, "_[a-z]$")
    x <- gsub("_[a-z]$", "", x)
    if (x %in% names(hmdb_mapper_vec)) {
      x <- paste("Metab__", hmdb_mapper_vec[x], sep = "")
    }
    if (!is.na(suffix)) {
      x <- paste(x, suffix, sep = "")
    }
    x
  }, character(1))
}

parse_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  pieces <- strsplit(lines, "\t", fixed = TRUE)
  rows <- lapply(pieces, function(x) {
    if (length(x) < 3) {
      return(NULL)
    }
    data.frame(source = x[1], target = x[-c(1, 2)], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

pathway_control_ora <- function(moon_res, meta_network, pathway_df,
                                top_score_cutoff = 2, background_cutoff = 1) {
  background_nodes <- moon_res$source[abs(moon_res$score) > background_cutoff]
  top_nodes <- moon_res$source[abs(moon_res$score) > top_score_cutoff]
  top_nodes <- top_nodes[top_nodes %in% c(meta_network$source, meta_network$target)]

  universe <- unique(background_nodes)
  out <- vector("list", length(top_nodes))
  i <- 1L

  for (node in top_nodes) {
    downstream_nodes <- unique(unlist(cosmosR:::keep_controllable_neighbours(
      meta_network,
      n_steps = 2,
      input_nodes = node
    )[, c("source", "target")]))
    downstream_nodes <- setdiff(downstream_nodes, node)
    downstream_nodes <- intersect(downstream_nodes, universe)

    if (length(downstream_nodes) == 0L) {
      next
    }

    pathway_terms <- split(pathway_df$target, pathway_df$source)
    res <- lapply(names(pathway_terms), function(term) {
      term_nodes <- intersect(unique(pathway_terms[[term]]), universe)
      if (length(term_nodes) == 0L) {
        return(NULL)
      }
      overlap <- intersect(downstream_nodes, term_nodes)
      k <- length(overlap)
      data.frame(
        node_of_interest = node,
        pathway = term,
        p_value = phyper(
          q = k - 1,
          m = length(term_nodes),
          n = length(universe) - length(term_nodes),
          k = length(downstream_nodes),
          lower.tail = FALSE
        ),
        overlap = k,
        downstream_n = length(downstream_nodes),
        pathway_n = length(term_nodes),
        stringsAsFactors = FALSE
      )
    })
    res <- do.call(rbind, res)
    if (!is.null(res) && nrow(res) > 0L) {
      out[[i]] <- res
      i <- i + 1L
    }
  }

  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(
      node_of_interest = character(),
      pathway = character(),
      p_value = numeric(),
      overlap = integer(),
      downstream_n = integer(),
      pathway_n = integer()
    )
  }
  out[order(out$p_value), , drop = FALSE]
}

run_tf_lig_threshold <- function(threshold, tf_weights, ligrec_inputs,
                                 dorothea_df, expressed_genes, rna_input) {
  dorothea_pkn <- dorothea_df[, c("source", "mor", "target")]
  names(dorothea_pkn) <- c("source", "interaction", "target")
  dorothea_pkn <- cosmosR:::filter_pkn_expressed_genes(
    names(expressed_genes),
    meta_pkn = dorothea_pkn
  )

  upstream_inputs <- setNames(tf_weights$feature_weights, tf_weights$Nodes)
  upstream_candidates_before_pkn <- upstream_inputs[abs(upstream_inputs) > threshold]

  lig_inputs <- ligrec_inputs$ligrec
  names(lig_inputs) <- gsub("___.+", "", names(lig_inputs))
  lig_inputs <- tapply(lig_inputs, names(lig_inputs), mean)
  lig_inputs <- as.numeric(lig_inputs)
  names(lig_inputs) <- names(tapply(
    ligrec_inputs$ligrec,
    gsub("___.+", "", names(ligrec_inputs$ligrec)),
    mean
  ))
  lig_inputs <- scale(lig_inputs, center = FALSE)
  downstream_inputs <- setNames(lig_inputs[, 1], row.names(lig_inputs))

  upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    upstream_candidates_before_pkn,
    dorothea_pkn
  )
  downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    downstream_inputs,
    dorothea_pkn
  )

  n_steps <- 1
  dorothea_pkn <- cosmosR:::keep_controllable_neighbours(
    dorothea_pkn,
    n_steps,
    names(upstream_inputs_filtered)
  )
  downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    downstream_inputs_filtered,
    dorothea_pkn
  )
  dorothea_pkn <- cosmosR:::keep_observable_neighbours(
    dorothea_pkn,
    n_steps,
    names(downstream_inputs_filtered)
  )
  upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    upstream_inputs_filtered,
    dorothea_pkn
  )

  meta_network_tf_lig <- dorothea_pkn
  before <- 1L
  after <- 0L
  i <- 1L
  moon_res <- NULL

  while (before != after && i < 10L) {
    before <- nrow(meta_network_tf_lig)
    moon_res <- cosmosR::moon(
      upstream_input = upstream_inputs_filtered,
      downstream_input = downstream_inputs_filtered,
      meta_network = meta_network_tf_lig,
      n_layers = n_steps,
      statistic = "ulm"
    )
    meta_network_tf_lig <- cosmosR:::filter_incohrent_TF_target(
      moon_res,
      dorothea_df,
      meta_network_tf_lig,
      rna_input
    )
    after <- nrow(meta_network_tf_lig)
    i <- i + 1L
  }

  moon_res <- as.data.frame(moon_res)
  moon_res$source <- as.character(moon_res$source)

  list(
    threshold = threshold,
    upstream_before_pkn = upstream_candidates_before_pkn,
    upstream_filtered = upstream_inputs_filtered,
    downstream_filtered = downstream_inputs_filtered,
    moon_res = moon_res,
    meta_network = meta_network_tf_lig,
    iterations = i - 1L
  )
}

objects <- new.env(parent = emptyenv())
load(file.path(repo_dir, "results", "mofa", "mofa_weights.RData"), envir = objects)
load(file.path(repo_dir, "data", "cosmos", "ligrec_TF_moon_inputs.Rdata"), envir = objects)
load(file.path(repo_dir, "support", "dorothea_df.RData"), envir = objects)
if (!all(c("weights", "ligrec_TF_moon_inputs", "dorothea_df") %in% ls(objects))) {
  stop("The NCI60 inputs did not contain weights, ligrec_TF_moon_inputs, and dorothea_df.")
}
weights <- objects$weights
ligrec_TF_moon_inputs <- objects$ligrec_TF_moon_inputs
dorothea_df <- objects$dorothea_df
data("HMDB_mapper_vec")

rna_input <- weights$RNA[, 4]
prot_input <- weights$proteo[, 4]
names(rna_input) <- gsub("_RNA", "", names(rna_input))
names(prot_input) <- gsub("_proteo", "", names(prot_input))
rna_log2_fpkm <- as.data.frame(read_csv(
  file.path(repo_dir, "data", "RNA", "RNA_log2_FPKM_clean.csv"),
  show_col_types = FALSE
))$Genes
expressed_gene_names <- rna_log2_fpkm

for (gene in names(rna_input)) {
  if (rna_input[gene] > -0.2 && rna_input[gene] < 0.2) {
    rna_input[gene] <- 0
  } else {
    rna_input[gene] <- sign(rna_input[gene]) * 10
  }
  if (gene %in% names(prot_input)) {
    if (prot_input[gene] > -0.05 && prot_input[gene] < 0.05) {
      rna_input[gene] <- 0
    } else {
      rna_input[gene] <- sign(rna_input[gene]) * 10
    }
  }
}

rna_log2_fpkm_missing <- rna_log2_fpkm[!(rna_log2_fpkm %in% names(rna_input))]
rna_log2_fpkm_vec <- rep(0, length(rna_log2_fpkm_missing))
names(rna_log2_fpkm_vec) <- rna_log2_fpkm_missing
rna_input <- c(rna_input, rna_log2_fpkm_vec)
expressed_genes <- setNames(rep(1, length(expressed_gene_names)), expressed_gene_names)

rec_weights <- ligrec_TF_moon_inputs$ligrec
names(rec_weights) <- gsub(".+___", "", names(rec_weights))
rec_weights <- tapply(rec_weights, names(rec_weights), mean)

lig_weights <- ligrec_TF_moon_inputs$ligrec
names(lig_weights) <- gsub("___.+", "", names(lig_weights))
lig_weights <- tapply(lig_weights, names(lig_weights), mean)

tf_weight_vec <- ligrec_TF_moon_inputs$TF
names(tf_weight_vec) <- gsub("_TF", "", names(tf_weight_vec))
tf_weight_vec <- tf_weight_vec[!(names(tf_weight_vec) %in% c(names(rec_weights), names(lig_weights)))]
tf_weights <- data.frame(
  Nodes = names(tf_weight_vec),
  feature_weights = as.numeric(tf_weight_vec)
)

baseline_rec_moon <- as.data.frame(read_csv(
  file.path(repo_dir, "results", "cosmos", "moon", "moon_res_rec_to_TFmet.csv"),
  show_col_types = FALSE
))
names(baseline_rec_moon)[1] <- "source"
baseline_rec_moon <- baseline_rec_moon[, c("source", "score", "level")]

baseline_rec_network <- as.data.frame(read_csv(
  file.path(repo_dir, "results", "cosmos", "moon", "meta_network_filtered.csv"),
  show_col_types = FALSE
))
baseline_rec_network$source <- translate_column_hmdb(
  baseline_rec_network$source,
  HMDB_mapper_vec
)
baseline_rec_network$target <- translate_column_hmdb(
  baseline_rec_network$target,
  HMDB_mapper_vec
)

pathway_df <- parse_gmt(file.path(repo_dir, "support", "c2.cp.v2022.1.Hs.symbols.gmt"))
pathway_df <- pathway_df[
  grepl("NABA_", pathway_df$source) | grepl("KEGG_", pathway_df$source),
  ,
  drop = FALSE
]

runs <- lapply(
  thresholds,
  run_tf_lig_threshold,
  tf_weights = tf_weights,
  ligrec_inputs = ligrec_TF_moon_inputs,
  dorothea_df = dorothea_df,
  expressed_genes = expressed_genes,
  rna_input = rna_input
)
names(runs) <- paste0("threshold_", thresholds)

combined_runs <- lapply(runs, function(run) {
  tf_moon <- run$moon_res[, c("source", "score", "level")]
  combined <- rbind(baseline_rec_moon, tf_moon) %>%
    group_by(source) %>%
    summarise(
      score = mean(score, na.rm = TRUE),
      level = mean(level, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  combined_network <- unique(rbind(
    baseline_rec_network[, c("source", "target", "interaction")],
    run$meta_network[, c("source", "target", "interaction")]
  ))

  list(
    combined_moon = combined,
    combined_network = combined_network,
    pathway_control = pathway_control_ora(combined, combined_network, pathway_df)
  )
})

candidate_census <- do.call(rbind, lapply(runs, function(run) {
  combined <- combined_runs[[paste0("threshold_", run$threshold)]]$combined_moon
  data.frame(
    threshold = run$threshold,
    upstream_candidates_before_pkn = length(run$upstream_before_pkn),
    upstream_candidates_after_pkn = length(run$upstream_filtered),
    downstream_ligands_after_pkn = length(run$downstream_filtered),
    tf_lig_scored_nodes = nrow(run$moon_res),
    combined_scored_nodes = nrow(combined),
    combined_abs_score_gt_1 = sum(abs(combined$score) > 1),
    combined_abs_score_gt_1_5 = sum(abs(combined$score) > 1.5),
    combined_abs_score_gt_2 = sum(abs(combined$score) > 2),
    iterations = run$iterations
  )
}))

baseline <- combined_runs[[paste0("threshold_", baseline_threshold)]]$combined_moon
baseline_top20 <- baseline$source[order(abs(baseline$score), decreasing = TRUE)][1:20]
baseline_top50 <- baseline$source[order(abs(baseline$score), decreasing = TRUE)][1:50]
baseline_gt2 <- baseline$source[abs(baseline$score) > 2]

score_stability <- do.call(rbind, lapply(thresholds, function(th) {
  current <- combined_runs[[paste0("threshold_", th)]]$combined_moon
  shared <- merge(
    baseline[, c("source", "score")],
    current[, c("source", "score")],
    by = "source",
    suffixes = c("_baseline", "_current")
  )
  current_top20 <- current$source[order(abs(current$score), decreasing = TRUE)][1:20]
  current_top50 <- current$source[order(abs(current$score), decreasing = TRUE)][1:50]
  current_gt2 <- current$source[abs(current$score) > 2]
  data.frame(
    threshold = th,
    shared_nodes_with_baseline = nrow(shared),
    spearman_score = suppressWarnings(cor(
      shared$score_baseline,
      shared$score_current,
      method = "spearman"
    )),
    pearson_score = suppressWarnings(cor(
      shared$score_baseline,
      shared$score_current,
      method = "pearson"
    )),
    sign_agreement = mean(sign(shared$score_baseline) == sign(shared$score_current)),
    top20_overlap_n = length(intersect(baseline_top20, current_top20)),
    top20_jaccard = length(intersect(baseline_top20, current_top20)) /
      length(union(baseline_top20, current_top20)),
    top50_overlap_n = length(intersect(baseline_top50, current_top50)),
    top50_jaccard = length(intersect(baseline_top50, current_top50)) /
      length(union(baseline_top50, current_top50)),
    gt2_overlap_n = length(intersect(baseline_gt2, current_gt2)),
    gt2_jaccard = length(intersect(baseline_gt2, current_gt2)) /
      length(union(baseline_gt2, current_gt2))
  )
}))

boundary_table <- data.frame(
  candidate = names(tf_weight_vec),
  feature_weight = as.numeric(tf_weight_vec),
  abs_feature_weight = abs(as.numeric(tf_weight_vec)),
  class = ifelse(
    abs(as.numeric(tf_weight_vec)) > 2.5,
    "retained_at_2.5",
    ifelse(
      abs(as.numeric(tf_weight_vec)) > 2,
      "lost_when_tightening_to_2.5",
      ifelse(
        abs(as.numeric(tf_weight_vec)) > 1.5,
        "gained_when_relaxing_to_1.5",
        "below_1.5"
      )
    )
  )
)
boundary_table <- boundary_table[boundary_table$class != "below_1.5", , drop = FALSE]
boundary_table <- boundary_table[order(boundary_table$class, -boundary_table$abs_feature_weight), , drop = FALSE]

top_node_preservation <- unique(unlist(lapply(thresholds, function(th) {
  current <- combined_runs[[paste0("threshold_", th)]]$combined_moon
  current$source[order(abs(current$score), decreasing = TRUE)][1:20]
})))
top_node_preservation <- do.call(rbind, lapply(top_node_preservation, function(node) {
  rows <- lapply(thresholds, function(th) {
    current <- combined_runs[[paste0("threshold_", th)]]$combined_moon
    current$rank_abs <- rank(-abs(current$score), ties.method = "min")
    row <- current[current$source == node, , drop = FALSE]
    if (nrow(row) == 0L) {
      return(data.frame(threshold = th, source = node, score = NA, rank_abs = NA, sign = NA))
    }
    data.frame(
      threshold = th,
      source = node,
      score = row$score[1],
      rank_abs = row$rank_abs[1],
      sign = sign(row$score[1])
    )
  })
  do.call(rbind, rows)
}))

baseline_pathway <- combined_runs[[paste0("threshold_", baseline_threshold)]]$pathway_control
baseline_pathway_top20 <- unique(baseline_pathway$pathway[order(baseline_pathway$p_value)])[1:min(
  20,
  length(unique(baseline_pathway$pathway))
)]
baseline_pair_top20 <- paste(
  baseline_pathway$node_of_interest,
  baseline_pathway$pathway,
  sep = "||"
)[order(baseline_pathway$p_value)][1:min(20, nrow(baseline_pathway))]

pathway_stability <- do.call(rbind, lapply(thresholds, function(th) {
  pathway_control <- combined_runs[[paste0("threshold_", th)]]$pathway_control
  top_pathways <- unique(pathway_control$pathway[order(pathway_control$p_value)])[1:min(
    20,
    length(unique(pathway_control$pathway))
  )]
  top_pairs <- paste(
    pathway_control$node_of_interest,
    pathway_control$pathway,
    sep = "||"
  )[order(pathway_control$p_value)][1:min(20, nrow(pathway_control))]
  data.frame(
    threshold = th,
    tested_node_pathway_pairs = nrow(pathway_control),
    unique_nodes_with_pathway_results = length(unique(pathway_control$node_of_interest)),
    unique_pathways_with_results = length(unique(pathway_control$pathway)),
    top20_pathway_overlap_n = length(intersect(baseline_pathway_top20, top_pathways)),
    top20_pathway_jaccard = length(intersect(baseline_pathway_top20, top_pathways)) /
      length(union(baseline_pathway_top20, top_pathways)),
    top20_pair_overlap_n = length(intersect(baseline_pair_top20, top_pairs)),
    top20_pair_jaccard = length(intersect(baseline_pair_top20, top_pairs)) /
      length(union(baseline_pair_top20, top_pairs))
  )
}))

pathway_top_table <- do.call(rbind, lapply(thresholds, function(th) {
  pathway_control <- combined_runs[[paste0("threshold_", th)]]$pathway_control
  pathway_control <- pathway_control[order(pathway_control$p_value), , drop = FALSE]
  pathway_control <- pathway_control[1:min(25, nrow(pathway_control)), , drop = FALSE]
  pathway_control$threshold <- th
  pathway_control[, c(
    "threshold", "node_of_interest", "pathway", "p_value", "overlap",
    "downstream_n", "pathway_n"
  )]
}))

write_csv(candidate_census, file.path(out_dir, "R1.5_nci60_candidate_census.csv"))
write_csv(score_stability, file.path(out_dir, "R1.5_nci60_score_stability.csv"))
write_csv(boundary_table, file.path(out_dir, "R1.5_nci60_boundary_gain_loss.csv"))
write_csv(top_node_preservation, file.path(out_dir, "R1.5_nci60_top_node_preservation.csv"))
write_csv(pathway_stability, file.path(out_dir, "R1.5_nci60_pathway_stability.csv"))
write_csv(pathway_top_table, file.path(out_dir, "R1.5_nci60_pathway_top_table.csv"))

summary_lines <- c(
  "# R1.5 NCI60 Threshold Sensitivity",
  "",
  "Thresholds tested: abs(upstream TF candidate score) > 1.5, 2, and 2.5.",
  "",
  "The comparison recombines each TF-to-ligand rerun with the stored receptor-to-TF/metabolite branch.",
  "It intentionally does not use results/cosmos/moon/full_moon_res_combined.csv.",
  "",
  "## Candidate Census",
  paste(capture.output(print(candidate_census, row.names = FALSE)), collapse = "\n"),
  "",
  "## MOON Score Stability Versus Threshold 2 Baseline",
  paste(capture.output(print(score_stability, row.names = FALSE)), collapse = "\n"),
  "",
  "## Pathway-Control Stability Versus Threshold 2 Baseline",
  paste(capture.output(print(pathway_stability, row.names = FALSE)), collapse = "\n"),
  "",
  "## Boundary Candidates",
  paste(capture.output(print(table(boundary_table$class))), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "R1.5_nci60_threshold_sensitivity_summary.md"))
cat(paste(summary_lines, collapse = "\n"), "\n")
