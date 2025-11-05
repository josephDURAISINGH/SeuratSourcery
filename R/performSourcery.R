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
#' 3. Harmonize genes and metadata,
#' 4. Integrate data with [basicIntegration()],
#' 5. Return a unified Seurat object with provenance metadata.
#'
#' @param path Folder or vector of dataset file paths.
#' @param method Integration method ("SCT" or "standard").
#' @param dims Dimensions for integration.
#' @param visualize Logical; whether to generate plots via [getSourceryReport()].
#' @param verbose Logical; show progress.
#' @return A harmonized and integrated Seurat object.
#' @examples
#' \dontrun{
#' final_obj <- performSourcery("data/", method = "SCT", visualize = TRUE)
#' }
#' @seealso [summonData()], [runeInspection()], [basicIntegration()]
#' @export

performSourcery<- function(dataset){
  # Summon
  datasets <- summonData("raw_data/")

  # Scan (preharmonization)
  precheck <- runeInspection(datasets)
  getSourceryReport(precheck)

  # Harmonize
  datasets <- activateRune(datasets)

  # Scan (postharmonization)
  postcheck <- runeInspection(datasets)
  getSourceryReport(postcheck)

  #Simple integration pipeline
  return(datasets)
}

#' Create demo data for SeuratSourcery vignettes
#' @noRd
demo_rune_data <- function(n_genes = 200, n_cells = 50) {
  mat <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("Gene", seq_len(n_genes)),
                    paste0("Cell", seq_len(n_cells)))
  )
  Seurat::CreateSeuratObject(counts = mat)
}
