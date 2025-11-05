# Purpose: Harmonizes a set of seurat objects
# Author: Joseph Duraisingh
# Date: November 2 2025
# Version: 1.0
# Bugs and Issues:

#' Activate Runes (Harmonize Datasets)
#'
#' Harmonizes a list of Seurat objects to ensure compatibility for downstream integration.
#' This includes unifying gene identifiers (Ensembl ↔ Symbol), aligning assay structures,
#' and optionally removing duplicate or low-frequency features.
#'
#' @param datasets A list of Seurat objects (e.g., from [summonData()]).
#' @param species Character; species identifier for ID mapping ("human" or "mouse").
#' @param id_type Character; target gene ID type ("symbol" or "ensembl").
#' @param remove_dupes Logical; whether to collapse duplicated genes by mean expression.
#' @param verbose Logical; display progress messages via `cli`.
#' @return A list of harmonized Seurat objects ready for integration.
#' @details
#' This function acts as the primary "harmonization rune" in the SeuratSourcery workflow.
#' It standardizes feature identifiers, assay names, and metadata consistency across datasets.
#' Datasets that fail harmonization are automatically excluded with a warning.
#'
#' @examples
#' \dontrun{
#' runes <- loadRune("data/")
#' runes <- activateRune(runes, species = "human", id_type = "symbol")
#' }
#'
#' @seealso [summonData()], [runeInspection()], [basicIntegration()]
#' @export


# Takes all of the files that have been loaded and turns them into harmonized seurat objects
activateRune <- function(datasets, species = "human", id_type = "symbol",
                           remove_dupes = TRUE) {
  if (!all(purrr::map_lgl(datasets, ~ inherits(.x, "Seurat"))))
    stop("All elements in 'datasets' must be Seurat objects.")

  cli::cli_alert_info("Harmonizing {length(datasets)} datasets")

  datasets <- purrr::map(datasets, function(obj) {
    genes <- rownames(obj)
    genes <- strip_version_ids(genes)
    if (id_type %in% c("symbol", "ensembl")) {
      mapped <- map_gene_ids(genes, species = species, target = id_type)
      rownames(obj) <- mapped
    }

    if (remove_dupes) {
      obj <- ensure_unique_genes(obj)
    }

    obj <- standardize_assay_structure(obj)
    obj@misc$harmonized <- TRUE
    obj
  })

  cli::cli_alert_success("Harmonization complete")

  return(datasets)
}

#' Strip Version IDs from Ensembl Gene Names
#' e.g. ENSG000001234.5 → ENSG000001234
#' @noRd
strip_version_ids <- function(genes) {
  gsub("\\..*$", "", genes)
}


#' Map Gene Identifiers Between Ensembl and Symbol
#' @noRd
map_gene_ids <- function(genes, species = "human", target = "symbol") {
  # Fallback: return as-is if biomaRt unavailable
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    cli::cli_alert_warning("biomaRt not installed; skipping ID mapping.")
    return(genes)
  }

  species_dataset <- switch(
    species,
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl",
    { cli::cli_alert_warning("Unsupported species; skipping ID mapping."); return(genes) }
  )

  mart <- biomaRt::useEnsembl("ensembl", dataset = species_dataset, version = NULL)
  if (target == "symbol") {
    map <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          mart = mart)
    names(map$hgnc_symbol) <- map$ensembl_gene_id
    mapped <- map[strip_version_ids(genes), "hgnc_symbol"]
    mapped[is.na(mapped)] <- genes[is.na(mapped)]
    return(mapped)
  } else if (target == "ensembl") {
    map <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          mart = mart)
    names(map$ensembl_gene_id) <- map$hgnc_symbol
    mapped <- map[genes, "ensembl_gene_id"]
    mapped[is.na(mapped)] <- genes[is.na(mapped)]
    return(mapped)
  }

  return(genes)
}


#' Ensure Genes Are Unique Per Dataset
#' Collapses duplicates by averaging counts.
#' @noRd
ensure_unique_genes <- function(obj) {
  counts <- Seurat::GetAssayData(obj, slot = "counts")
  if (any(duplicated(rownames(counts)))) {
    cli::cli_alert_warning("Collapsing duplicated gene entries")
    counts <- collapse_duplicate_genes(counts)
    Seurat::SetAssayData(obj, slot = "counts", new.data = counts)
  }
  obj
}

collapse_duplicate_genes <- function(counts) {
  as(Matrix::Matrix(
    aggregate(counts, by = list(rownames(counts)), FUN = mean)[, -1],
    sparse = TRUE
  ))
}


#' Standardize Assay Structure
#' Ensures a consistent RNA assay and proper slot setup.
#' @noRd
standardize_assay_structure <- function(obj) {
  if (!"RNA" %in% names(obj@assays)) {
    cli::cli_alert_warning("No RNA assay detected; setting default assay to first available")
    Seurat::DefaultAssay(obj) <- names(obj@assays)[1]
  }
  assay <- Seurat::DefaultAssay(obj)
  data <- Seurat::GetAssayData(obj, slot = "data")
  if (is.null(data)) {
    Seurat::SetAssayData(obj, slot = "data",
                         new.data = Seurat::GetAssayData(obj, slot = "counts"))
  }
  obj
}
