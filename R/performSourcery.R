# Purpose: Wrapper for the entire preprocessing pipeline
# Author: Joseph Duraisingh
# Date: November 3 2025
# Version: 1.0
# Bugs and Issues:

#' Perform Full SeuratSourcery Workflow
#'
#' High-level wrapper that executes the full SeuratSourcery pipeline:
#' 1. Load datasets with [summonData()],
#' 2. Inspect with [runeInspection()],
#' 3. Optionally output graphs with [getSourceryReport()]
#' 4. Harmonize genes and metadata with [activateRune()]
#' 5. Inspect again with[runInspection()]
#' 6. Optionally visualize again with [getSourceryReport()]
#' 5. Return harmonized seurat objects
#'
#' @param path Optional folder path containing dataset files readable by \code{summonData()}.
#' @param datasets Optional list of Seurat objects. Used if \code{path} is not supplied.
#' @param visualize Logical; whether to generate summary plots from \code{getSourceryReport()},
#'  (both before and after harmonization).
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A harmonized list of seurat objects
#' @examples
#' @examples
#' \dontrun{
#'   # Load from folder:
#'   out <- performSourcery(path = "inst/extdata/demo_datasets")
#'   # Or supply datasets directly:
#'   dl <- summonData("inst/extdata/demo_datasets")
#'   out <- performSourcery(datasets = dl)
#' }
#'
#' @seealso [summonData()], [runeInspection()], [activateRune()], [getSourceryReport()]
#' @export

performSourcery <- function(path = NULL,
                            datasets = NULL,
                            visualize = TRUE,
                            verbose = TRUE) {
  # Summon
  if (!is.null(path)) {
    if (verbose) message("Summoning datasets from: ", path)
    raw <- summonData(path)
  } else {
    if (is.null(datasets))
      stop("Supply either `path` or `datasets`.", call. = FALSE)
    raw <- datasets
  }

  if (!is.list(raw))
    stop("`datasets` must be a list of Seurat objects.", call. = FALSE)


  # Scan (preharmonization)
  if (verbose) message("Running pre-harmonization inspection...")
  precheck <- runeInspection(raw)

  if (visualize) {
    if (verbose) message("Generating pre-harmonization report...")
    getSourceryReport(precheck)
  }

  # Harmonize
  if (verbose) message("Harmonizing datasets with activateRune()...")
  datasets <- activateRune(raw)


  # Scan (postharmonization)
  if (verbose) message("Running post-harmonization inspection...")
  postcheck <- runeInspection(datasets)

  if (visualize) {
    if (verbose) message("Generating post-harmonization report...")
    getSourceryReport(postcheck)
  }

  # return
  if (verbose) message("SeuratSourcery preprocessing complete.")
  return(datasets)
}
