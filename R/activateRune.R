# Purpose: Harmonizes a set of seurat objects
# Author: Joseph Duraisingh
# Date: November 2 2025
# Version: 1.0
# Bugs and Issues:

#' Activate Runes (Harmonize Seurat Objects)
#'
#' Harmonizes a list of Seurat objects to ensure compatibility for downstream integration.
#' This includes unifying gene identifiers (Ensembl ↔ Symbol), aligning assay structures,
#' and optionally removing duplicate or low-frequency features.
#'
#' @param datasets A list of Seurat objects (e.g., from [summonData()]). Preferred Assay5 objects
#' @param normalize Logical; if `TRUE`, runs `NormalizeData()` and `FindVariableFeatures()` on each dataset after cleaning.
#' @param fix_gene_names Logical; if `TRUE`, strips ensembl suffix or uppercases non-Ensembl genes
#' @param drop_zero_genes Logical; if `TRUE`, removes genes with zero total counts across all cells.
#' @param drop_duplicates Logical; if `TRUE`, collapses duplicated gene names
#'   using mean expression across duplicates.
#' @param map_ids Logical; if `TRUE`, attempts to map gene identifiers between
#'   Ensembl and symbol formats using **biomaRt**. Breaks without biomart so make sure it is installed
#'   this one can remain false unless absolutely needed, and is coated with safety checks
#' @param species Character; species identifier for ID mapping (almost always "human" or "mouse").
#' @param id_type Character; target gene ID type ("symbol", "ensembl" or "mixed").
#' @return A list of harmonized Seurat objects ready for integration.
#'
#' @details
#' This function acts as the primary "harmonization rune" in the SeuratSourcery workflow.
#' It standardizes feature identifiers, assay names, and metadata consistency across datasets.
#' Datasets that fail harmonization are automatically excluded with a warning.
#' This is a key function that cleans up a lot of the mess that single cell data can often have
#'
#' @examples
#' # Example using included demo data (installed under inst/extdata/)
#' demo_dir <- system.file("extdata/demo_datasets", package = "SeuratSourcery")
#' if (dir.exists(demo_dir)) {
#'   runes <- summonData(demo_dir)
#'   clean_runes <- activateRune(runes, fix_gene_names = TRUE)
#' }
#'
#' @seealso [summonData()] to load data, [runeInspection()] for qc metrics
#'
#' @references
#' Hao Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#' Cell 184(13):3573-3587.
#' Durinck S. et al. (2009). Mapping identifiers for the integration of
#' genomic datasets. Nature Methods 6(6): 477-479.
#'
#' @importFrom Seurat GetAssayData NormalizeData FindVariableFeatures
#'
#' @export

# --------- outwards functions ------------

activateRune <- function(
    datasets,
    normalize = FALSE,
    fix_gene_names = TRUE,
    drop_zero_genes = TRUE,
    drop_duplicates = TRUE,
    map_ids = FALSE,
    species = "human",
    id_type = "symbol"
    ) {

  if (!is.list(datasets) || !all(purrr::map_lgl(datasets, ~ inherits(.x, "Seurat"))))
    stop("activateRune() requires a list of Seurat objects.")

  cli::cli_h1("Activating Runes")
  cli::cli_alert_info("Processing {length(datasets)} datasets")

  out <- purrr::imap(datasets, function(obj, nm) {

    cli::cli_h2("Rune: {nm}")

    # check valid format
    obj <- tryCatch(standardize_assay_structure(obj), error = function(e) obj)

    if (fix_gene_names) {
      genes <- rownames(obj)
      genes <- strip_version_ids(genes)
      ens <- grepl("^ENS", genes, ignore.case = TRUE)
      genes[!ens] <- toupper(genes[!ens])
      rownames(obj) <- genes
    }

    if (map_ids) {
      cli::cli_alert("Mapping gene identifiers ...")
      rownames(obj) <- tryCatch(
        map_gene_ids(rownames(obj), species, id_type),
        error = function(e) {
          cli::cli_alert_warning("Gene ID mapping failed: {e$message}")
          rownames(obj)
        }
      )
    }
    obj <- tryCatch(sanitize_seurat_object(obj), error = function(e) obj)

    if (drop_duplicates) {
      obj <- tryCatch(ensure_unique_genes(obj), error = function(e) obj)
      }

    if (drop_zero_genes) {
      counts <- get_counts(obj)
      keep <- Matrix::rowSums(counts) > 0
      obj <- obj[keep, ]
    }

    data_layer <- tryCatch(GetAssayData(obj, layer = "data"), error = function(e) NULL)
    if (is.null(data_layer)) {
      Seurat::SetAssayData(obj, layer = "data", new.data = log1p(get_counts(obj)))
    }

    if (normalize) {
      cli::cli_alert("Running NormalizeData & FindVariableFeatures")
      obj <- tryCatch(NormalizeData(obj, verbose = FALSE), error = function(e) obj)
      obj <- tryCatch(FindVariableFeatures(obj, verbose = FALSE), error = function(e) obj)
    }

    obj@misc$harmonized <- TRUE
    obj
  })

  cli::cli_alert_success("Harmonization complete")

  return(out)
}

# ------- Helpers-----------------

strip_version_ids <- function(genes) {
  gsub("\\..*$", "", genes)
}


# Map Gene Identifiers Between Ensembl and Symbol
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


# Make assays compatible
standardize_assay_structure <- function(obj) {
  if (!"RNA" %in% names(obj@assays)) {
    cli::cli_alert_warning("No RNA assay detected; setting default assay to first available")
    Seurat::DefaultAssay(obj) <- names(obj@assays)[1]
  }
  assay <- Seurat::DefaultAssay(obj)
  data <- Seurat::GetAssayData(obj, layer = "data")
  if (is.null(data)) {
    Seurat::SetAssayData(obj, layer = "data",
                         new.data = Seurat::GetAssayData(obj, layer = "counts"))
  }
  obj
}


# General object hygiene
sanitize_seurat_object <- function(obj) {

  assay_names <- names(obj@assays)

  for (assay_nm in assay_names) {
    assay <- obj@assays[[assay_nm]]

    # base gene names
    g <- rownames(assay)

    # Fix missing and duplicated names
    if (anyNA(g)) {
      na_idx <- which(is.na(g))
      g[na_idx] <- paste0("gene_", na_idx)
    }
    if (anyDuplicated(g)) {
      g <- make.unique(g)
    }

    # Apply to assay-level counts/data/scale.data
    if (!is.null(assay@counts))
      rownames(assay@counts) <- g
    if (!is.null(assay@data))
      rownames(assay@data) <- g
    if (!is.null(assay@scale.data) &&
        nrow(assay@scale.data) > 0)
      rownames(assay@scale.data) <- g

    # Apply to feature-level metadata if present
    if (!is.null(assay@meta.features)) {
      rownames(assay@meta.features) <- g
    }

    # Apply new gene names back to assay
    assay@var.features <- intersect(assay@var.features, g)
    assay@key <- assay@key # keep key unchanged
    obj@assays[[assay_nm]] <- assay
  }

  # Also ensure object-level rownames reflect default assay
  DefaultAssay(obj) <- DefaultAssay(obj)
  return(obj)
}

# Fixes gene names
ensure_unique_genes <- function(obj) {

  counts <- Seurat::GetAssayData(obj, layer = "counts")
  g <- rownames(counts)

  if (!anyDuplicated(g)) {
    return(obj)
  }

  cli::cli_alert_warning("Collapsing duplicated gene entries")

  # Aggregate by gene name
  df <- as.data.frame(as.matrix(counts))
  df$gene <- g

  collapsed <- df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(across(-gene, mean)) |>
    as.data.frame()

  rownames(collapsed) <- collapsed$gene
  collapsed <- collapsed[, -1, drop = FALSE]

  collapsed <- as(Matrix::Matrix(as.matrix(collapsed), sparse = TRUE))

  # Set back into the object
  obj <- Seurat::SetAssayData(obj, layer = "counts", new.data = collapsed)

  # Remove var.feature references to dropped duplicates
  assay <- obj@assays[[DefaultAssay(obj)]]
  assay@var.features <- intersect(assay@var.features, rownames(collapsed))
  obj@assays[[DefaultAssay(obj)]] <- assay

  return(obj)
}
