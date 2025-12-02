# Purpose: Checks statistics between the different datasets
# Author: Joseph Duraisingh
# Date: November 3 2025
# Version: 1.0
# Bugs and Issues:

#' Inspect Runes for Compatibility
#'
#' Performs a pre-harmonization scan across Seurat objects.
#' Summarizes gene counts, identifier types, assay structures, and metadata.
#' Compares datasets by quantifying overlaps.
#' Uses key lightweight diagnostics to help determine compatability between datsets for integration
#'
#' @param datasets A list of Seurat objects.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Tibble summarizing dataset-level metrics.}
#'   \item{overlap}{Matrix of feature overlaps between datasets.}
#'   \item{common_meta}{Vector of metadata fields shared by all datasets.}
#'   \item{unique_meta}{List of the unique metadata fields}
#'   \item{datasets}{The input datasets get passed out to be accessed downstream}
#' }
#'
#' @details
#' Does not modify the seurat objects, only extracts data out of them.
#' Checking for differences and whats going on under the hood for each dataset
#'
#' @examples
#' #using demo data
#' demo_dir <- system.file("extdata/demo_datasets", package = "SeuratSourcery")
#' if (dir.exists(demo_dir)) {
#'   runes <- summonData(demo_dir)
#'   insp <- runeInspection(runes)
#'   print(insp$summary)
#' }
#' @seealso [summonData()] for loading data,
#'  [activateRune()] for harmonization,
#'  [getSourceryReport()] for visualizations
#' @export

# ------ Outward facing function -------------

runeInspection <- function(datasets) {
  if (!is.list(datasets)) stop("'datasets' must be a list", call. = FALSE)
  if (is.null(names(datasets))) names(datasets) <- paste0("Dataset_", seq_along(datasets))
  if (!all(purrr::map_lgl(datasets, ~ inherits(.x, "Seurat"))))
    stop("All items in the list must be Seurat objects.", call. = FALSE)

  nm <- names(datasets)

  # take rna counts out
  counts_list <- purrr::map(datasets, get_counts)

  # gene + cell summary
  summary_tbl <- tibble::tibble(
    Dataset = nm,
    Genes = purrr::map_int(counts_list, nrow),
    Cells = purrr::map_int(counts_list, ncol),
    Mean_NonZero_Genes = purrr::map_dbl(
      counts_list,
      ~ mean(Matrix::colSums(.x > 0))
    ),
    Mean_UMI_Per_Cell = purrr::map_dbl(
      counts_list,
      ~ mean(Matrix::colSums(.x))
    ),
    Example_Gene = purrr::map_chr(
      counts_list,
      ~ rownames(.x)[1] %||% NA_character_
    )
  )

  # Gene type guess
  # summary_tbl$ID_Type <- purrr::map_chr(summary_tbl$Example_Gene, function(g) {
  #   if (is.na(g)) return("unknown")
  #   if (grepl("^ENS", g, ignore.case = TRUE)) return("ensembl")
  #   if (grepl("^[A-Za-z0-9-]+$", g)) return("symbol")
  #   return("mixed")
  # })
  summary_tbl$ID_Type <- purrr::map_chr(summary_tbl$Example_Gene, guess_gene_type)


  # Gene overlap
  gene_lists <- purrr::map(counts_list, rownames)

  overlap_mat <- outer(
    seq_along(gene_lists), seq_along(gene_lists),
    Vectorize(function(i, j) {
      length(intersect(gene_lists[[i]], gene_lists[[j]]))
    })
  )
  rownames(overlap_mat) <- colnames(overlap_mat) <- nm

  # metadata overlap
  meta_lists <- purrr::map(datasets, ~ colnames(.x[]))
  common_meta <- Reduce(intersect, meta_lists)
  unique_meta <- purrr::map(meta_lists, setdiff, y = common_meta)

  # return
  output <- list(
    summary = summary_tbl,
    overlap = overlap_mat,
    common_meta = common_meta,
    unique_meta = unique_meta,
    datasets = datasets
  )
  return(output)
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

get_counts <- function(obj) {
  assay <- obj[[Seurat::DefaultAssay(obj)]]

  # ---- v5 Assay5 structure ----
  if ("layers" %in% slotNames(assay)) {

    # gene and cell names stored separately
    featurenames <- rownames(assay)
    cellnames    <- colnames(assay)

    # 1. try raw counts
    if ("counts" %in% names(assay@layers)) {
      mat <- assay@layers$counts
      if (!is.null(mat)) {
        rownames(mat) <- featurenames
        colnames(mat) <- cellnames
        return(mat)
      }
    }

    # or...
    if ("data" %in% names(assay@layers)) {
      mat <- assay@layers$data
      if (!is.null(mat)) {
        rownames(mat) <- featurenames
        colnames(mat) <- cellnames
        return(mat)
      }
    }
  }

  # v3/v4 fallback
  if ("counts" %in% slotNames(assay)) {
    mat <- assay@counts
    featurenames <- rownames(mat)
    cellnames    <- colnames(mat)
    if (is.null(featurenames)) featurenames <- paste0("gene_", seq_len(nrow(mat)))
    if (is.null(cellnames)) cellnames <- paste0("cell_", seq_len(ncol(mat)))
    rownames(mat) <- featurenames
    colnames(mat) <- cellnames
    return(mat)
  }

  stop("No usable counts or data matrix found.")
}
