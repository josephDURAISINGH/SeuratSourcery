# Purpose: Reads files from a folder or one at a time to convert them into the initial seurat objects
# Author: Joseph Duraisingh
# Date: November 1 2025
# Version: 1.0
# Bugs and Issues:

#' Load Runes (Datasets)
#'
#' Loads single-cell datasets of multiple formats into Seurat objects.
#' Accepts folders or file paths and automatically detects supported types
#' (e.g. RDS, H5AD, 10X, CSV/TSV).
#' Note: RDS and 10x much more thoroughly tested
#'
#' @param path A directory to folder or list of files to load.
#' @param pattern Optional regex to filter filenames.
#' @param include_dirs Logical; if `TRUE` (default), subdirectories are checked
#' Useful for loading in 10x data directories
#'
#' @return A named list of Seurat objects, one per dataset successfully loaded.
#'
#' @examples
#' # Loading in the demo data
#' demo_dir <- system.file("extdata/demo_datasets", package = "SeuratSourcery")
#' if (dir.exists(demo_dir)) {
#'   ds <- summonData(demo_dir)
#'   print(names(ds))
#' }
#'
#' @seealso [loadRune()] for loading in individual files,
#'  [activateRune()] for harmonization,
#'  [runeInspection()] for loading qc metrics
#'
#' @references
#' Hao Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#' Cell 184(13):3573-3587.
#'
#' @export


# ------ Outward Facing Functions --------------------

summonData <- function(path = ".", pattern = NULL, include_dirs = TRUE) {

  # Check it is a directory
  if (!(length(path) == 1 && dir.exists(path))) {
    stop("summonData() expects a directory. For single files, use add_single_dataset().")
  }
  # load files
  files <- list.files(path, full.names = TRUE, pattern = pattern)
  files <- files[!dir.exists(files)]

  #Check for 10x folders
  if (include_dirs) {

    subdirs <- list.dirs(path, recursive = FALSE)
    subdirs <- subdirs[subdirs != path]  # remove the main directory

    if (length(subdirs) > 0) {

      tenx_dirs <- subdirs[vapply(subdirs, is_tenx_dir, logical(1))]
      files <- c(files, tenx_dirs)
    }
  }

  if (length(files) == 0)
    stop("No supported files found in directory: ", path)

  cli::cli_alert_info("Detected {length(files)} files/folders")

  # load all the files
  datasets <- purrr::map(files, safely_load_dataset)
  names(datasets) <- make.unique(basename(files))

  # keep only things that turned into seurat objects
  datasets <- purrr::keep(datasets, ~ inherits(.x, "Seurat"))
  datasets <- purrr::compact(datasets)

  cli::cli_alert_success("Loaded {length(datasets)} Seurat objects")

  return(datasets)
}

# for loading a dataset (rune) into existing list one at a time
#' Adds in a single dataset as a seurat object
#' @param file filepath of data object to be loaded
#' @param data optional parameter for whatever object the other already added datasets are stored in,
#' to append them into a list
#'
#' @return List of seurat objects with new one appended
#'
#' @examples
#' # Adding a single dataset from extdata
#' demo_file <- system.file("extdata/demo_datasets/pbmc_demo1.rds",
#'                          package = "SeuratSourcery")
#' if (file.exists(demo_file)) {
#'   lst <- loadRune(demo_file)
#' }
#'
#' @export

loadRune <- function(file, data = NULL) {

  obj <- safely_load_dataset(file)

  if (is.null(obj)) {
    cli::cli_alert_warning("Could not load {.file {file}}; returned NULL.")
    return(data)
  }

  # Make sure we always return a list
  if (is.null(data)) {
    return(list(obj))
  } else if (is.list(data)) {
    data[[length(data) + 1]] <- obj
    return(data)
  } else {
    stop("loadRune(): 'data' must be a list or NULL.", call. = FALSE)
  }
}

# ---------- Helper Functions -----------------------

safely_load_dataset <- function(file) {
  # Normalize files
  raw_ext <- tools::file_ext(file)
  raw_ext <- sub("\\.gz$", "", raw_ext)
  ext <- tolower(raw_ext)
  cli::cli_alert("Loading {.file {basename(file)}}")

  obj <- tryCatch({
  out <- switch(
    ext,
    "rds" = {
      x <- readRDS(file)
      if (!inherits(x, "Seurat")) {
        cli::cli_alert_warning("{.file {file}} is an RDS but not a Seurat object")
        NULL  # skip
      } else x
    },
    "h5ad" = {
      if (requireNamespace("zellkonverter", quietly = TRUE)) {
        sce <- getNamespace("zellkonverter")[["readH5AD"]](file)
        Seurat::as.Seurat(sce)
      } else {
        cli::cli_alert_warning("Reading .h5ad requires zellkonverter. Install via BiocManager::install('zellkonverter')")
        NULL
      }
    },
    "csv" = load_csv_as_seurat(file, sep = ","),
    "tsv" = load_csv_as_seurat(file, sep = "\t"),

    # If 10x folder included
    {
      if (dir.exists(file) && is_tenx_dir(file)) {
        Seurat::Read10X(data.dir = file) |> Seurat::CreateSeuratObject()
      } else {
        cli::cli_alert_warning("Unsupported file type: {.file {file}}")
        NULL
      }
    }
  )
  out
  }, error = function(e) {
    cli::cli_alert_warning("Failed to load {.file {basename(file)}}: {e$message}")
    return(NULL)
  })

  if (!is.null(obj)) {

    # Sanitize gene names
    g <- rownames(obj)

    # Replace NA names
    if (anyNA(g)) {
      cli::cli_alert_warning("Missing gene names detected in {.file {basename(file)}}; replacing NA with placeholders.")
      na_idx <- which(is.na(g))
      g[na_idx] <- paste0("gene_", na_idx)
    }

    # Fix duplicates
    if (anyDuplicated(g)) {
      cli::cli_alert_warning("Duplicate gene names detected in {.file {basename(file)}}; making unique.")
      g <- make.unique(g)
    }

    # Apply sanitized names across all layers
    rownames(obj) <- g

    # Store metadata
    obj@misc$source_file <- basename(file)
    obj@misc$file_type <- ext
  }

  return(obj)
}


# Helper function for 10x
is_tenx_dir <- function(d) {
  f <- list.files(d)
  has_mtx  <- any(grepl("^matrix\\.mtx(\\.gz)?$", f))
  has_bar  <- any(grepl("^barcodes\\.tsv(\\.gz)?$", f))
  has_feat <- any(grepl("^features\\.tsv(\\.gz)?$", f))
  has_h5   <- any(grepl("\\.h5$", f))
  (has_mtx && has_bar && has_feat) || has_h5
}

# Helper function for CSV
load_csv_as_seurat <- function(file, sep = ",") {
  # Read file safely
  mat <- tryCatch(
    readr::read_delim(
      file,
      delim = sep,
      col_types = readr::cols(.default = readr::col_double(), gene = readr::col_character()),
      progress = FALSE
    ),
    error = function(e) {
      cli::cli_alert_warning("Failed to read {.file {file}}: {e$message}")
      return(NULL)
    }
  )

  if (is.null(mat) || nrow(mat) == 0) {
    cli::cli_alert_warning("Empty or invalid CSV for {.file {file}}")
    return(NULL)
  }

  # Identify gene column
  gene_col <- which(tolower(names(mat)) %in% c("gene", "genes", "symbol", "id", "ensembl", "gene_id"))

  if (length(gene_col) == 0) {
    cli::cli_alert_info("Assuming first column contains gene names in {.file {file}}")
    gene_col <- 1
  }

  gene_names <- mat[[gene_col]]

  # --- Step 3: Sanitize gene names ---
  if (anyNA(gene_names)) {
    cli::cli_alert_warning("Missing gene names in {.file {file}}; replacing NAs with placeholders")
    na_idx <- which(is.na(gene_names))
    gene_names[na_idx] <- paste0("gene_", na_idx)
  }

  if (anyDuplicated(gene_names)) {
    cli::cli_alert_warning("Duplicate gene names detected in {.file {file}}; making unique")
    gene_names <- make.unique(gene_names)
  }

  # Extract counts
  count_df <- mat[, -gene_col, drop = FALSE]

  # Validate columns are numeric
  if (!all(vapply(count_df, is.numeric, logical(1)))) {
    cli::cli_alert_warning("Non-numeric columns detected in {.file {file}}; attempting conversion")
    count_df <- data.frame(lapply(count_df, as.numeric))
  }

  count_mat <- as.matrix(count_df)
  rownames(count_mat) <- gene_names

  # Try to create seurat object
  obj <- tryCatch(
    Seurat::CreateSeuratObject(counts = count_mat),
    error = function(e) {
      cli::cli_alert_warning("Failed to create Seurat object from {.file {file}}: {e$message}")
      return(NULL)
    }
  )
  return(obj)
}
