#' Generate a Small Demo Seurat Object
#'
#' Creates a lightweight synthetic Seurat object for testing
#'
#' @param n_genes number of genes (rows).
#' @param n_cells number of cells (columns).
#'
#' @return A Seurat object containing a simulated sparse counts matrix.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Seurat CreateSeuratObject
demo_rune_data <- function(n_genes = 200, n_cells = 50) {
  mat <- matrix(
    stats::rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  Seurat::CreateSeuratObject(counts = mat)
}
