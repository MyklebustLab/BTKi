#!/usr/bin/env Rscript

# ==============================================================================
# Add HGNC / OrgDb Gene Annotation to RNA-seq Count Matrices
# ==============================================================================
# Author: Ankush Sharma
# Project: CLL RNA-seq treatment-response analysis
#
# Description:
#   This script annotates featureCounts and HTSeq count matrices with human gene
#   symbols, gene names, Entrez identifiers, and gene type information using
#   Bioconductor AnnotationDbi and org.Hs.eg.db.
#
# Expected input files:
#   02_count_matrices/final_featurecounts/featurecounts_gene_annotation.tsv
#   02_count_matrices/final_featurecounts/featurecounts_count_matrix.tsv
#   02_count_matrices/final_featurecounts/htseq_count_matrix.tsv
#
# Outputs:
#   02_count_matrices/final_featurecounts/featurecounts_gene_annotation_hgnc.tsv
#   02_count_matrices/final_featurecounts/featurecounts_count_matrix_annotated_hgnc.tsv
#   02_count_matrices/final_featurecounts/htseq_count_matrix_annotated_hgnc.tsv
#   02_count_matrices/final_featurecounts/hgnc_annotation_summary.tsv
#
# Notes:
#   - Run build_annotated_count_matrices.R before this script.
#   - Ensembl IDs are mapped through org.Hs.eg.db using keytype = "ENSEMBL".
#   - When an Ensembl ID maps to multiple entries, the first OrgDb match is used.
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path, call. = FALSE)
  }

  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
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

validate_required_columns <- function(data, required_columns, file_label) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      file_label, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

map_orgdb_column <- function(ensembl_ids, column_name) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    column = column_name,
    multiVals = "first"
  )
}

build_hgnc_annotation <- function(gene_annotation) {
  validate_required_columns(
    gene_annotation,
    c("Geneid", "ensembl_id", "Chr", "Start", "End", "Strand", "Length"),
    "featurecounts_gene_annotation.tsv"
  )

  ensembl_ids <- unique(gene_annotation$ensembl_id)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids) & ensembl_ids != ""]

  mapped <- data.frame(
    ensembl_id = ensembl_ids,
    gene_symbol = unname(map_orgdb_column(ensembl_ids, "SYMBOL")),
    gene_name = unname(map_orgdb_column(ensembl_ids, "GENENAME")),
    entrez_id = unname(map_orgdb_column(ensembl_ids, "ENTREZID")),
    gene_type = unname(map_orgdb_column(ensembl_ids, "GENETYPE")),
    stringsAsFactors = FALSE
  )

  annotated <- merge(gene_annotation, mapped, by = "ensembl_id", all.x = TRUE, sort = FALSE)
  annotated <- annotated[match(gene_annotation$Geneid, annotated$Geneid), , drop = FALSE]
  annotated[
    ,
    c(
      "Geneid", "ensembl_id", "gene_symbol", "gene_name", "entrez_id", "gene_type",
      "Chr", "Start", "End", "Strand", "Length"
    ),
    drop = FALSE
  ]
}

add_annotation_to_count_matrix <- function(count_matrix, annotation_hgnc, file_label, all_x = TRUE) {
  validate_required_columns(count_matrix, "Geneid", file_label)

  annotated <- merge(
    annotation_hgnc,
    count_matrix,
    by = "Geneid",
    all.x = all_x,
    all.y = !all_x,
    sort = FALSE
  )

  if (all_x) {
    annotated <- annotated[match(annotation_hgnc$Geneid, annotated$Geneid), , drop = FALSE]
  } else {
    annotated <- annotated[match(count_matrix$Geneid, annotated$Geneid), , drop = FALSE]
    annotated$ensembl_id <- sub("\\..*$", "", annotated$Geneid)
  }

  annotation_columns <- names(annotation_hgnc)
  sample_columns <- setdiff(names(annotated), annotation_columns)
  annotated[, c(annotation_columns, sample_columns), drop = FALSE]
}

build_annotation_summary <- function(annotation_hgnc) {
  data.frame(
    metric = c(
      "genes_total",
      "gene_symbols_mapped",
      "gene_names_mapped",
      "entrez_ids_mapped",
      "gene_types_mapped"
    ),
    value = c(
      nrow(annotation_hgnc),
      sum(!is.na(annotation_hgnc$gene_symbol) & annotation_hgnc$gene_symbol != ""),
      sum(!is.na(annotation_hgnc$gene_name) & annotation_hgnc$gene_name != ""),
      sum(!is.na(annotation_hgnc$entrez_id) & annotation_hgnc$entrez_id != ""),
      sum(!is.na(annotation_hgnc$gene_type) & annotation_hgnc$gene_type != "")
    ),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  script_dir <- get_script_dir()
  analysis_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
  count_dir <- file.path(analysis_dir, "02_count_matrices", "final_featurecounts")

  feature_annotation_file <- file.path(count_dir, "featurecounts_gene_annotation.tsv")
  feature_matrix_file <- file.path(count_dir, "featurecounts_count_matrix.tsv")
  htseq_matrix_file <- file.path(count_dir, "htseq_count_matrix.tsv")

  gene_annotation <- read_tsv(feature_annotation_file)
  feature_matrix <- read_tsv(feature_matrix_file)
  htseq_matrix <- read_tsv(htseq_matrix_file)

  gene_annotation_hgnc <- build_hgnc_annotation(gene_annotation)
  feature_annotated <- add_annotation_to_count_matrix(
    feature_matrix,
    gene_annotation_hgnc,
    "featurecounts_count_matrix.tsv",
    all_x = TRUE
  )
  htseq_annotated <- add_annotation_to_count_matrix(
    htseq_matrix,
    gene_annotation_hgnc,
    "htseq_count_matrix.tsv",
    all_x = FALSE
  )

  summary <- build_annotation_summary(gene_annotation_hgnc)

  write_tsv(gene_annotation_hgnc, file.path(count_dir, "featurecounts_gene_annotation_hgnc.tsv"))
  write_tsv(feature_annotated, file.path(count_dir, "featurecounts_count_matrix_annotated_hgnc.tsv"))
  write_tsv(htseq_annotated, file.path(count_dir, "htseq_count_matrix_annotated_hgnc.tsv"))
  write_tsv(summary, file.path(count_dir, "hgnc_annotation_summary.tsv"))

  message("Wrote HGNC/OrgDb annotations to: ", count_dir)
  print(summary)

  write_session_info(file.path(count_dir, "session_info_hgnc_annotation.txt"))
}

main()
