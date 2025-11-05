# Purpose: Checks statistics between the different datasets
# Author: Joseph Duraisingh
# Date: November 3 2025
# Version: 1.0
# Bugs and Issues:

#' Inspect Runes for Compatibility
#'
#' Performs a pre-harmonization scan across Seurat objects.
#' Summarizes gene counts, identifier types, assay structures, and metadata.
#'
#' @param datasets A list of Seurat objects.
#' @param sample_names Optional character vector naming each dataset.
#' @param verbose Logical; display progress messages.
#' @return A list containing:
#' \describe{
#'   \item{summary}{Tibble summarizing dataset-level metrics.}
#'   \item{overlap}{Matrix of feature overlaps between datasets.}
#'   \item{common_meta}{Vector of metadata fields shared by all datasets.}
#'   \item{notes}{Text summary.}
#' }
#' @examples
#' \dontrun{
#' runes <- loadRune("data/")
#' report <- runeInspection(runes)
#' }
#' @seealso [getSourceryReport()], [runeInspection()]
#' @export

# ------ Outward facing function -------------
runeInspection <- function(datasets, sample_names = NULL) {
  if (!all(purrr::map_lgl(datasets, ~ inherits(.x, "Seurat"))))
    stop("All elements in 'datasets' must be Seurat objects.")

  n <- length(datasets)
  if (is.null(sample_names)) sample_names <- paste0("Dataset_", seq_len(n))

  cli::cli_alert_info("Scanning {n} datasets for harmonization compatibility")

  # Basic gene summary
  gene_counts <- purrr::map_int(datasets, ~ nrow(.x))
  example_genes <- purrr::map_chr(datasets, ~ head(rownames(.x), 1))
  id_type <- purrr::map_chr(example_genes, guess_gene_type)

  # Gene overlap matrix
  overlap_mat <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j) {
    length(intersect(rownames(datasets[[i]]), rownames(datasets[[j]])))
  }))
  overlap_perc <- round(overlap_mat / matrix(gene_counts, n, n, byrow = TRUE) * 100, 1)

  # Assay structure summary
  assay_names <- purrr::map_chr(datasets, ~ Seurat::DefaultAssay(.x))
  assay_slots <- purrr::map_int(datasets, ~ length(Seurat::Assays(.x)))

  # Metadata summary
  metadata_fields <- purrr::map(datasets, ~ colnames(.x@meta.data))
  common_meta <- Reduce(intersect, metadata_fields)
  unique_meta <- purrr::map(metadata_fields, setdiff, y = common_meta)

  # Build summary table
  summary_tbl <- tibble::tibble(
    Dataset = sample_names,
    Genes = gene_counts,
    ID_Type = id_type,
    Default_Assay = assay_names,
    Num_Assays = assay_slots,
    Common_Metadata = lengths(common_meta),
    Unique_Metadata = purrr::map_chr(unique_meta, ~ paste(.x, collapse = ", "))
  )

  # Results object
  results <- list(
    summary = summary_tbl,
    overlap = overlap_perc,
    common_meta = common_meta,
    notes = cli::col_green("Scan complete: check summary and overlap for consistency.")
  )

  cli::cli_h1("Pre-Harmonization Summary")
  print(summary_tbl)
  cli::cli_h2("Gene Overlap (%)")
  print(round(overlap_perc, 1))
  cli::cli_alert_info("Common metadata fields: {paste(common_meta, collapse = ', ')}")

  return(results)
}

# --- Helper ----

#' Guess Gene Identifier Type
#' Simple heuristic to classify as Ensembl, Symbol, or Mixed.
#' @noRd
guess_gene_type <- function(gene_vector) {
  genes <- head(gene_vector, 50)
  if (all(grepl("^ENS", genes, ignore.case = TRUE))) return("ensembl")
  if (all(grepl("^[A-Za-z0-9-]+$", genes))) return("symbol")
  return("mixed/unknown")
}


