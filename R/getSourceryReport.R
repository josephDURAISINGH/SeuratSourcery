# Purpose: Outputs graphical representation of the datasets and their compatibility
# Author: Joseph Duraisingh
# Date: November 4 2025
# Version: 1.0
# Bugs and Issues:

#' Visualize Rune Metrics
#'
#' Produces graphical summaries of dataset characteristics prior to integration,
#' including gene counts per dataset and conserved gene proportions.
#'
#' @param datasets A list of Seurat objects or the output of [runeInspection()].
#' @param overlap_results Optional precomputed overlap matrix.
#' @param show_plot Logical; whether to display plots immediately.
#' @return A list containing ggplot2 objects:
#' \describe{
#'   \item{gene_bar}{Bar plot of total genes per dataset.}
#'   \item{conserved_bar}{Bar plot of conserved genes across datasets.}
#' }
#' @examples
#' \dontrun{
#' report <- runeInspection(runes)
#' plots <- getSourceryReport(report)
#' }
#' @seealso [runeInspection()], [performSourcery()]
#' @export

getSourceryReport <- function(datasets, overlap_results = NULL, show_plot = TRUE) {
  if (inherits(datasets, "list") && "summary" %in% names(datasets)) {
    summary_tbl <- datasets$summary
    overlap_mat <- datasets$overlap
  } else {
    # compute simple summary from Seurat list
    summary_tbl <- tibble::tibble(
      Dataset = names(datasets),
      Genes = purrr::map_int(datasets, ~ nrow(.x)),
      Cells = purrr::map_int(datasets, ~ ncol(.x))
    )
    overlap_mat <- outer(seq_along(datasets), seq_along(datasets),
                         Vectorize(function(i, j)
                           length(intersect(rownames(datasets[[i]]),
                                            rownames(datasets[[j]])))))
  }

  # Bar chart: total genes per dataset
  p_genes <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = Dataset, y = Genes, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = "Number of Genes per Dataset", y = "Genes", x = NULL) +
    ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.4, size = 3.5)

  # Shared / conserved genes
  all_genes <- unique(unlist(lapply(datasets, rownames)))
  presence <- sapply(all_genes, function(g)
    mean(sapply(datasets, function(obj) g %in% rownames(obj)))
  )
  conserved_count <- sum(presence == 1)
  shared_75 <- sum(presence >= 0.75)
  shared_50 <- sum(presence >= 0.5)

  df_conserved <- tibble::tibble(
    Category = c("Conserved (100%)", "Shared ≥75%", "Shared ≥50%"),
    Genes = c(conserved_count, shared_75, shared_50)
  )

  p_conserved <- ggplot2::ggplot(df_conserved, ggplot2::aes(x = Category, y = Genes, fill = Category)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = "Conserved Gene Proportions", x = NULL, y = "Gene Count") +
    ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.3, size = 3.5)

  if (show_plot) {
    print(p_genes)
    print(p_conserved)
  }

  return(list(gene_bar = p_genes, conserved_bar = p_conserved))
}
