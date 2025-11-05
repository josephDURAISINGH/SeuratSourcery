# Purpose: Reads files from a folder or one at a time to convert them into the initial seurat objects
# Author: Joseph Duraisingh
# Date: November 1 2025
# Version: 1.0
# Bugs and Issues:

#' Load Runes (Datasets)
#'
#' Loads single-cell datasets of multiple formats into Seurat objects.
#' Accepts folders or file paths and automatically detects supported types
#' (e.g. RDS, H5Seurat, H5AD, 10X, CSV/TSV).
#'
#' @param path A folder or vector of file paths to load.
#' @param pattern Optional regex to filter filenames.
#' @return A named list of Seurat objects, one per dataset loaded.
#' @examples
#' \dontrun{
#' data_list <- loadRune("data/")
#' }
#' @seealso [activateRune()], [runeInspection()]
#' @export


# ------ Outward Facing Functions --------------------

summonData <- function(path = ".", pattern = NULL, include_dirs = TRUE) {
  if (dir.exists(path)) {
    files <- list.files(path, full.names = TRUE, pattern = pattern)
    if (include_dirs) {
      dirs <- list.dirs(path, recursive = FALSE)
      dirs <- dirs[sapply(dirs, function(d) any(grepl("barcodes", list.files(d))))]
      files <- c(files, dirs)
    }
  } else if (all(file.exists(path))) {
    files <- path
  } else {
    stop("Invalid path")
  }
  if (length(files) == 0) stop("No files found at path.")

  cli::cli_alert_info("Detected {length(files)} files")

  # This is the important line that does the data loading
  datasets <- purrr::map(files, safely_load_dataset)
  datasets <- purrr::keep(datasets, ~ inherits(.x, "Seurat"))
  names(datasets) <- basename(files)
  # remove NULLs
  datasets <- purrr::compact(datasets)

  cli::cli_alert_success("Loaded {length(datasets)} datasets")

  return(datasets)
}

# for loading them in one at a time
#' Adds in a single dataset as a seurat object
#' @param data whatever object the other already added datasets are stored in
#' @param path filepath of data to be loaded
#' @return List of seurat objects with new one included
#' @examples
#' \dontrun{
#' data_list <- loadRune("")
#' }
#' @export
add_single_dataset <- function(data = NULL, path){
  return(c(data, safely_load_dataset(path)))
}

# ---------- Helper Functions -----------------------

safely_load_dataset <- function(file) {
  ext <- tolower(tools::file_ext(file))
  cli::cli_alert("Loading {.file {basename(file)}}")

  obj <- switch(
    ext,
    "rds" = readRDS(file),
    "h5seurat" = {
      if (requireNamespace("SeuratDisk", quietly = TRUE)) {
        SeuratDisk::LoadH5Seurat(file)
      } else {
        cli::cli_alert_warning("Package 'SeuratDisk' not available; skipping {.file {file}}")
        NULL
      }
    },
    "h5ad" = {
      if (requireNamespace("zellkonverter", quietly = TRUE)) {
        sce <- zellkonverter::readH5AD(file)
        Seurat::as.Seurat(sce)
      } else {
        cli::cli_alert_warning("Package 'zellkonverter' not available; skipping {.file {file}}")
        NULL
      }
    },
    "csv" = load_csv_as_seurat(file, sep = ","),
    "tsv" = load_csv_as_seurat(file, sep = "\t"),

    # If 10x folder included
    {
      if (dir.exists(file) && any(grepl("barcodes", list.files(file)))) {
        Seurat::Read10X(data.dir = file) |> Seurat::CreateSeuratObject()
      } else {
        cli::cli_alert_warning("Unsupported file type: {.file {file}}")
        NULL
      }
    }
  )

  if (!is.null(obj)) {
    obj@misc$source_file <- basename(file)
    obj@misc$file_type <- ext
  }
  return(obj)
}

# Helper function for CSV
load_csv_as_seurat <- function(file, sep = ",") {
  mat <- tryCatch({
    readr::read_delim(file, delim = sep, progress = FALSE)
  }, error = function(e) {
    cli::cli_alert_warning("Failed to read {.file {file}}: {e$message}")
    return(NULL)
  })
  if (is.null(mat)){
    return(NULL)
  }
  if (!"gene" %in% names(mat)){
    names(mat)[1] <- "gene"
  }
  rownames(mat) <- mat$gene
  mat <- mat[, -1, drop = FALSE]

  Seurat::CreateSeuratObject(counts = as.matrix(mat))
}
