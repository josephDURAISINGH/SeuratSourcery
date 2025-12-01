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

getSourceryReport <- function(rune_obj, show_plots = TRUE) {

  if (!all(c("summary","overlap","datasets") %in% names(rune_obj))) {
    stop("getSourceryReport() requires the output of runeInspection().")
  }

  summary_tbl <- rune_obj$summary
  overlap_mat <- rune_obj$overlap
  seurat_list <- rune_obj$datasets

  # Other metrics
  qc_tbl <- purrr::map_dfr(
    names(seurat_list),
    function(nm) {
      obj <- seurat_list[[nm]]
      counts <- get_counts(obj)

      umi    <- Matrix::colSums(counts)
      feats  <- Matrix::colSums(counts > 0)

      mt_genes <- grep("^MT-", rownames(counts), ignore.case = TRUE, value = TRUE)
      pct_mito <- if (length(mt_genes) > 0) {
        Matrix::colSums(counts[mt_genes,,drop=FALSE]) / umi
      } else rep(0,length(umi))

      tibble::tibble(
        Dataset = nm,
        Mean_UMI = mean(umi),
        Mean_Features = mean(feats),
        Mean_Mito = mean(pct_mito),
        Sparsity = 1 - Matrix::nnzero(counts)/(nrow(counts)*ncol(counts)),
        Memory_MB = as.numeric(object.size(obj))/1e6
      )
    }
  )
# Conserved genes
  gene_lists <- lapply(seurat_list, function(obj) rownames(get_counts(obj)))
  all_genes  <- unique(unlist(gene_lists))

  presence <- sapply(all_genes, function(g)
    mean(sapply(gene_lists, function(lst) g %in% lst))
  )

  df_conserved <- tibble::tibble(
    Category = c("Conserved 100%", "Shared ≥75%", "Shared ≥50%"),
    Genes    = c(
      sum(presence == 1),
      sum(presence >= 0.75),
      sum(presence >= 0.5)
    )
  )

  # plots!
  p_genes <- ggplot2::ggplot(summary_tbl, aes(Dataset, Genes, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(aes(label = Genes), vjust = -0.3)

  p_mito <- ggplot2::ggplot(qc_tbl, aes(Dataset, Mean_Mito, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal()

  p_feats <- ggplot2::ggplot(qc_tbl, aes(Dataset, Mean_Features, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal()

  p_conserved <- ggplot2::ggplot(df_conserved, aes(Category, Genes, fill = Category)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal()

  #outputs
  if (show_plots) {
    print(p_genes)
    print(p_feats)
    print(p_mito)
    print(p_conserved)
  }

  list(
    summary = summary_tbl,
    qc = qc_tbl,
    overlap = overlap_mat,
    conserved = df_conserved,
    plots = list(
      genes = p_genes,
      features = p_feats,
      mito = p_mito,
      conserved = p_conserved
    )
  )

}
########### V2

# getSourceryReport <- function(datasets, show_plots = TRUE) {
#   # Extract metadata
#   summary_tbl <- purrr::map_dfr(
#     names(datasets),
#     function(nm) {
#       obj <- datasets[[nm]]
#       counts <- obj@assays$RNA@counts
#
#       # QC metrics
#       umi <- Matrix::colSums(counts)
#       feats <- Matrix::colSums(counts > 0)
#
#       # mitochondria content
#       mt_genes <- grep("^MT-", rownames(counts), value = TRUE, ignore.case = TRUE)
#       pct_mito <- if (length(mt_genes) > 0) {
#         Matrix::colSums(counts[mt_genes, , drop = FALSE]) / umi
#       } else rep(0, length(umi))
#
#       tibble::tibble(
#         Dataset = nm,
#         Cells = ncol(counts),
#         Genes = nrow(counts),
#         Memory_MB = as.numeric(object.size(obj)) / 1e6,
#         Sparsity = 1 - (Matrix::nnzero(counts) / (ncol(counts) * nrow(counts))),
#         Mean_UMI = mean(umi),
#         Median_UMI = median(umi),
#         Mean_Features = mean(feats),
#         Median_Features = median(feats),
#         Mean_Mito = mean(pct_mito),
#         Median_Mito = median(pct_mito)
#       )
#     }
#   )
#
#   # diagnostics for data
#   diagnostics <- purrr::map(
#     datasets,
#     function(obj) {
#       genes <- rownames(obj)
#
#       list(
#         duplicates = unique(genes[duplicated(genes)]),
#         has_duplicates = anyDuplicated(genes) > 0,
#         gene_type = dplyr::case_when(
#           all(grepl("^ENS", genes)) ~ "ENSEMBL",
#           all(grepl("^[A-Za-z]", genes)) ~ "SYMBOL",
#           TRUE ~ "mixed"
#         )
#       )
#     }
#   )
#   names(diagnostics) <- names(datasets)
#
#   # finding overlap of datasets
#   gene_lists <- lapply(datasets, rownames)
#   overlap_mat <- outer(
#     names(gene_lists), names(gene_lists),
#     Vectorize(function(a, b)
#       length(intersect(gene_lists[[a]], gene_lists[[b]])))
#   )
#   rownames(overlap_mat) <- colnames(overlap_mat) <- names(gene_lists)
#
#   # gene coservation categories (100%, 75%, 50%)
#   all_genes <- unique(unlist(gene_lists))
#   presence <- sapply(all_genes, function(g)
#     mean(sapply(gene_lists, function(lst) g %in% lst))
#   )
#
#   df_conserved <- tibble::tibble(
#     Category = c("Conserved 100%", "Shared ≥75%", "Shared ≥50%"),
#     Genes = c(sum(presence == 1),
#               sum(presence >= 0.75),
#               sum(presence >= 0.5))
#   )
#
#   # Plots
#   p_genes <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = Dataset, y = Genes, fill = Dataset)) +
#     ggplot2::geom_col(show.legend = FALSE) +
#     ggplot2::theme_minimal() +
#     ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.4) +
#     ggplot2::labs(title = "Gene Count per Dataset", y = "Genes", x = NULL)
#
#   # Features per dataset
#   p_features <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = Dataset, y = Median_Features, fill = Dataset)) +
#     ggplot2::geom_col(show.legend = FALSE) +
#     ggplot2::theme_minimal() +
#     ggplot2::labs(title = "Median Features per Cell", y = "nFeatures", x = NULL)
#
#   # Mito %
#   p_mito <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = Dataset, y = Mean_Mito, fill = Dataset)) +
#     ggplot2::geom_col(show.legend = FALSE) +
#     ggplot2::theme_minimal() +
#     ggplot2::labs(title = "Mean Mitochondrial Percentage", y = "% mito", x = NULL)
#
#   # Conserved gene
#   p_conserved <- ggplot2::ggplot(df_conserved, ggplot2::aes(x = Category, y = Genes, fill = Category)) +
#     ggplot2::geom_col(show.legend = FALSE) +
#     ggplot2::theme_minimal() +
#     ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.4) +
#     ggplot2::labs(title = "Conserved / Shared Genes", x = NULL, y = "Gene Count")
#
#   if (show_plots) {
#     print(p_genes)
#     print(p_features)
#     print(p_mito)
#     print(p_conserved)
#   }
#
#   # Return report
#   list(
#     summary = summary_tbl,
#     diagnostics = diagnostics,
#     overlap = overlap_mat,
#     plots = list(
#       gene_count = p_genes,
#       feature_count = p_features,
#       mito = p_mito,
#       conserved = p_conserved
#     )
#   )
#
#
#
# }
