# Purpose: Wrapper for a basic downstream Seurat integration analysis
# Author: Joseph Duraisingh
# Date: November 4 2025
# Version: 1.0
# Bugs and Issues:

#' Perform Basic Seurat Integration
#'
#' Wraps the standard Seurat integration workflow for RNA datasets.
#' Automatically selects integration features, finds anchors, and integrates data.
#' Supports both standard and SCTransform-based normalization.
#'
#' @param datasets A list of Seurat objects (ideally harmonized).
#' @param method Integration method: `"SCT"` or `"standard"`.
#' @param dims Dimensions to use for integration.
#' @param verbose Logical; whether to show progress.
#' @return A Seurat object containing the integrated dataset.
#' @examples
#' \dontrun{
#' integrated <- basicIntegration(runes, method = "SCT")
#' }
#' @seealso [performSourcery()], [Seurat::FindIntegrationAnchors()], [Seurat::IntegrateData()]
#' @export

basicIntegration <- function(datasets, method = "SCT", dims = 1:30) {
  cli::cli_alert_info("Starting integration with {length(datasets)} datasets...")

  if (method == "SCT") {
    datasets <- lapply(datasets, Seurat::SCTransform, verbose = FALSE)
    features <- Seurat::SelectIntegrationFeatures(datasets)
    datasets <- Seurat::PrepSCTIntegration(datasets, anchor.features = features)
    anchors <- Seurat::FindIntegrationAnchors(datasets, normalization.method = "SCT",
                                              anchor.features = features)
    integrated <- Seurat::IntegrateData(anchors, normalization.method = "SCT")
  } else {
    features <- Seurat::SelectIntegrationFeatures(datasets)
    anchors <- Seurat::FindIntegrationAnchors(datasets, anchor.features = features)
    integrated <- Seurat::IntegrateData(anchors)
  }

  cli::cli_alert_success("Integration complete!")
  return(integrated)
}
