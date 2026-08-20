#!/usr/bin/env Rscript

# ==============================================================================
# Build RNA-seq Count Matrices from featureCounts and HTSeq Outputs
# ==============================================================================
# Author: Ankush Sharma
# Project: CLL RNA-seq treatment-response analysis
#
# Description:
#   This script builds sample metadata and gene-by-sample count matrices from
#   featureCounts and HTSeq count files. It is the R implementation of the
#   original Python count-matrix builder and is intended for GitHub/reproducible
#   analysis use.
#
# Expected input folder structure:
#   ../InputReadCounts/
#     featurecounts/<sample_id>/<sample_id>_featureCounts.txt
#     htseqReadCounts/<sample_id>/<sample_id>_htseq.txt
#
# Expected sample naming convention:
#   CLL149_01_DM
#     CLL149 = patient identifier
#     01     = treatment order
#     DM     = treatment code
#
# Outputs:
#   02_count_matrices/final_featurecounts/sample_metadata.tsv
#   02_count_matrices/final_featurecounts/featurecounts_gene_annotation.tsv
#   02_count_matrices/final_featurecounts/featurecounts_count_matrix.tsv
#   02_count_matrices/final_featurecounts/featurecounts_count_matrix_annotated.tsv
#   02_count_matrices/final_featurecounts/htseq_count_matrix.tsv
#   02_count_matrices/final_featurecounts/htseq_count_matrix_with_featurecounts_annotation.tsv
#   02_count_matrices/final_featurecounts/count_matrix_build_summary.txt
#
# Notes:
#   - featureCounts is used as the primary downstream count matrix in the
#     treatment-response analysis.
#   - HTSeq counts are parsed as a parallel count source and retained for
#     traceability/comparison.
#   - Run add_hgnc_gene_annotation.R after this script to add gene symbols,
#     gene names, Entrez IDs, and gene type annotations.
# ==============================================================================

options(stringsAsFactors = FALSE)

treatment_order <- c(
  "01" = "DM",
  "02" = "IB",
  "03" = "AC",
  "04" = "ZA",
  "05" = "PI",
  "06" = "NE",
  "07" = "RO",
  "08" = "NX",
  "09" = "BG"
)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

read_tsv <- function(path, comment.char = "") {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = comment.char
  )
}

write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_session_info <- function(output_path) {
  session_lines <- capture.output(sessionInfo())
  writeLines(session_lines, output_path)
  message("\nSession info:")
  message(paste(session_lines, collapse = "\n"))
}

stop_if_missing_dir <- function(path, label) {
  if (!dir.exists(path)) {
    stop(label, " directory not found: ", path, call. = FALSE)
  }
}

parse_sample_id <- function(sample_id) {
  match <- regexec("^(CLL[0-9]+)_([0-9]{2})_([A-Za-z0-9]+)$", sample_id)
  pieces <- regmatches(sample_id, match)[[1]]
  if (length(pieces) != 4) {
    stop("Unexpected sample name: ", sample_id, call. = FALSE)
  }

  patient <- pieces[2]
  order_code <- pieces[3]
  treatment <- pieces[4]
  expected_treatment <- unname(treatment_order[order_code])

  if (!is.na(expected_treatment) && treatment != expected_treatment) {
    stop(
      "Sample ", sample_id, " has treatment ", treatment,
      ", expected ", expected_treatment, " from order ", order_code,
      call. = FALSE
    )
  }

  data.frame(
    sample_id = sample_id,
    patient = patient,
    treatment_order = as.integer(order_code),
    treatment = treatment,
    is_reference_control = treatment == "DM",
    stringsAsFactors = FALSE
  )
}

discover_count_files <- function(parent_dir, suffix, label) {
  stop_if_missing_dir(parent_dir, label)
  files <- list.files(
    parent_dir,
    pattern = paste0(suffix, "$"),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No ", label, " files found under: ", parent_dir, call. = FALSE)
  }

  samples <- basename(dirname(files))
  names(files) <- samples
  files[order(names(files))]
}

make_sample_metadata <- function(samples) {
  metadata <- do.call(rbind, lapply(samples, parse_sample_id))
  metadata <- metadata[order(metadata$patient, metadata$treatment_order), ]
  rownames(metadata) <- NULL

  metadata$paired_design_block <- metadata$patient
  metadata$recommended_deseq2_design <- "~ patient + treatment"
  metadata$recommended_reference_level <- "DM"
  metadata
}

read_featurecounts_file <- function(path) {
  featurecounts <- read_tsv(path, comment.char = "#")
  required_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  missing_columns <- setdiff(required_columns, names(featurecounts))
  if (length(missing_columns) > 0) {
    stop(
      "featureCounts file is missing required columns: ",
      paste(missing_columns, collapse = ", "), "\nFile: ", path,
      call. = FALSE
    )
  }

  count_column <- tail(names(featurecounts), 1)
  list(
    annotation = featurecounts[, required_columns, drop = FALSE],
    counts = setNames(as.integer(featurecounts[[count_column]]), featurecounts$Geneid)
  )
}

build_featurecounts_outputs <- function(feature_files, samples) {
  first_sample <- samples[1]
  first_read <- read_featurecounts_file(feature_files[[first_sample]])

  gene_annotation <- first_read$annotation
  gene_annotation$ensembl_id <- sub("\\..*$", "", gene_annotation$Geneid)
  gene_annotation <- gene_annotation[
    ,
    c("Geneid", "ensembl_id", "Chr", "Start", "End", "Strand", "Length"),
    drop = FALSE
  ]

  count_list <- lapply(samples, function(sample_id) {
    read_featurecounts_file(feature_files[[sample_id]])$counts
  })
  names(count_list) <- samples

  all_gene_ids <- unique(unlist(lapply(count_list, names), use.names = FALSE))
  count_matrix <- data.frame(Geneid = all_gene_ids, stringsAsFactors = FALSE)
  for (sample_id in samples) {
    counts <- count_list[[sample_id]]
    count_matrix[[sample_id]] <- as.integer(counts[match(all_gene_ids, names(counts))])
  }
  count_matrix[is.na(count_matrix)] <- 0L

  list(annotation = gene_annotation, matrix = count_matrix)
}

read_htseq_file <- function(path) {
  htseq <- read.delim(
    path,
    sep = "\t",
    header = FALSE,
    col.names = c("Geneid", "count"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  htseq <- htseq[!grepl("^__", htseq$Geneid), , drop = FALSE]
  setNames(as.integer(htseq$count), htseq$Geneid)
}

build_htseq_matrix <- function(htseq_files, samples) {
  count_list <- lapply(samples, function(sample_id) read_htseq_file(htseq_files[[sample_id]]))
  names(count_list) <- samples

  all_gene_ids <- unique(unlist(lapply(count_list, names), use.names = FALSE))
  count_matrix <- data.frame(Geneid = all_gene_ids, stringsAsFactors = FALSE)
  for (sample_id in samples) {
    counts <- count_list[[sample_id]]
    count_matrix[[sample_id]] <- as.integer(counts[match(all_gene_ids, names(counts))])
  }
  count_matrix[is.na(count_matrix)] <- 0L
  count_matrix
}

merge_by_geneid_preserve_left <- function(annotation, count_matrix) {
  merged <- merge(annotation, count_matrix, by = "Geneid", all.x = TRUE, sort = FALSE)
  merged[match(annotation$Geneid, merged$Geneid), , drop = FALSE]
}

main <- function() {
  script_dir <- get_script_dir()
  analysis_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
  input_dir <- normalizePath(file.path(analysis_dir, "..", "InputReadCounts"), mustWork = TRUE)
  featurecounts_dir <- file.path(input_dir, "featurecounts")
  htseq_dir <- file.path(input_dir, "htseqReadCounts")
  output_dir <- file.path(analysis_dir, "02_count_matrices", "final_featurecounts")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  feature_files <- discover_count_files(featurecounts_dir, "_featureCounts.txt", "featureCounts")
  htseq_files <- discover_count_files(htseq_dir, "_htseq.txt", "HTSeq")
  metadata <- make_sample_metadata(names(feature_files))
  samples <- metadata$sample_id

  missing_htseq <- setdiff(samples, names(htseq_files))
  extra_htseq <- setdiff(names(htseq_files), samples)
  if (length(missing_htseq) > 0 || length(extra_htseq) > 0) {
    stop(
      "featureCounts and HTSeq sample sets do not match.\n",
      "Missing HTSeq: ", paste(missing_htseq, collapse = ", "), "\n",
      "Extra HTSeq: ", paste(extra_htseq, collapse = ", "),
      call. = FALSE
    )
  }

  feature_outputs <- build_featurecounts_outputs(feature_files, samples)
  htseq_matrix <- build_htseq_matrix(htseq_files, samples)

  feature_annotated <- merge_by_geneid_preserve_left(
    feature_outputs$annotation,
    feature_outputs$matrix
  )
  htseq_annotated <- merge(feature_outputs$annotation, htseq_matrix, by = "Geneid", all.y = TRUE, sort = FALSE)
  htseq_annotated$ensembl_id <- sub("\\..*$", "", htseq_annotated$Geneid)

  write_tsv(metadata, file.path(output_dir, "sample_metadata.tsv"))
  write_tsv(feature_outputs$annotation, file.path(output_dir, "featurecounts_gene_annotation.tsv"))
  write_tsv(feature_outputs$matrix, file.path(output_dir, "featurecounts_count_matrix.tsv"))
  write_tsv(feature_annotated, file.path(output_dir, "featurecounts_count_matrix_annotated.tsv"))
  write_tsv(htseq_matrix, file.path(output_dir, "htseq_count_matrix.tsv"))
  write_tsv(
    htseq_annotated,
    file.path(output_dir, "htseq_count_matrix_with_featurecounts_annotation.tsv")
  )

  summary_lines <- c(
    "CLL RNA-seq count matrix build summary",
    paste("Input directory:", input_dir),
    paste("Output directory:", output_dir),
    paste("Patients:", paste(sort(unique(metadata$patient)), collapse = ", ")),
    paste("Treatments:", paste(unique(metadata$treatment[order(metadata$treatment_order)]), collapse = ", ")),
    paste("Samples:", length(samples)),
    paste("featureCounts genes/features:", nrow(feature_outputs$matrix)),
    paste("HTSeq genes/features excluding __ summary rows:", nrow(htseq_matrix)),
    "Recommended paired differential expression design: ~ patient + treatment",
    "Recommended reference/control treatment level: DM"
  )
  writeLines(summary_lines, file.path(output_dir, "count_matrix_build_summary.txt"))
  message(paste(summary_lines, collapse = "\n"))

  write_session_info(file.path(output_dir, "session_info_count_matrix_build.txt"))
}

main()
