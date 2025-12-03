# Purpose: Run the Shiny app
# Author: Joseph Duraisingh
# Date: December 1 2025
# Version: 1.0
# Bugs and Issues:

#' Launch the SeuratSourcery Shiny Application
#'
#' Opens the interactive Shiny interface for the SeuratSourcery workflow.
#' The app provides a graphical environment for loading datasets, running the
#' harmonization rune, inspecting QC metrics, and visualizing plots generated
#' by \code{runeInspection()} and \code{getSourceryReport()}.
#'
#' This function simply launches the app located in the package's
#' \code{inst/shiny-scripts/} directory. The Shiny app is optional, and this
#' function requires the \pkg{shiny} package to be installed.
#'
#' @details
#' The Shiny application mirrors the core SeuratSourcery workflow:
#' \enumerate{
#'   \item Load single-cell datasets of various formats,
#'   \item Inspect dataset compatibility,
#'   \item Activate the harmonization rune,
#'   \item Visualize QC and gene-overlap plots.
#' }
#'
#' The Shiny app does not run automatically on package load and is entirely
#' optional, ensuring CRAN compliance.
#'
#' @return No return value; the function is called for its side effect
#'   of launching the Shiny application.
#'
#' @examples
#' \dontrun{
#'   runSeuratSourcery()
#' }
#'
#' @references
#' Chang W. et al. (2023). shiny: Web Application Framework for R. R package
#' version 1.8.0.
#'
#' @importFrom utils packageVersion
#'
#' @export


runSeuratSourcery <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "SeuratSourcery")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
