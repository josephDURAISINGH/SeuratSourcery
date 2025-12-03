# Purpose: Outputs graphical representation of the datasets and their compatibility
# Author: Joseph Duraisingh
# Date: November 4 2025
# Version: 1.0
# Bugs and Issues:

#' Visualize Rune Metrics (QC report)
#'
#' Produces graphical summaries of dataset characteristics prior to integration,
#' including gene counts per dataset and conserved gene proportions.
#'
#' @param rune_obj The output of [runeInspection()].
#' @param show_plot Logical; if true prints out plots
#'
#' @return A list containing ggplot2 objects:
#' \describe{
#'   \item{summary}{Dataset-level gene/cell summary table.}
#'   \item{qc}{QC metrics table.}
#'   \item{overlap}{Gene overlap matrix.}
#'   \item{conserved}{Summary of conserved/shared genes.}
#'   \item{plots}{A list of ggplot2 objects: numgenes, numcells, numfeatures,
#'   mito content (genes that start with MT-), umi and conserved genes}
#' }
#'
#'
#' @examples
#' demo_dir <- system.file("extdata/demo_datasets", package = "SeuratSourcery")
#' if (dir.exists(demo_dir)) {
#'   runes <- summonData(demo_dir)
#'   insp <- runeInspection(runes)
#'   rep <- getSourceryReport(insp, show_plots = FALSE)
#'   rep$plots$genes
#' }
#'
#' @seealso
#' [runeInspection()], [performSourcery()], [activateRune()], [summonData()]
#'
#' @references
#' Hao Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#' Cell 184(13):3573-3587.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text theme_minimal
#'
#' @export

getSourceryReport <- function(rune_obj, show_plot = TRUE) {

  required <- c("summary", "overlap", "datasets")
  if (!all(required %in% names(rune_obj))) {
    stop("getSourceryReport() requires the output of runeInspection().",
         call. = FALSE)
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
      umi <- Matrix::colSums(counts)
      feats  <- Matrix::colSums(counts > 0)
      mt_genes <- grep("(^MT-|^mt-)", rownames(counts), ignore.case = TRUE, value = TRUE)
      pct_mito <- if (length(mt_genes) > 0) {
        Matrix::colSums(counts[mt_genes,,drop=FALSE]) / umi
      } else{ rep(0,length(umi))}

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
    mean(vapply(gene_lists, function(lst) g %in% lst, logical(1)))
  )

  df_conserved <- tibble::tibble(
    Category = c("AllGenes", "Shared50plus", "Shared75plus", "Shared100"),
    Genes = c(
      length(all_genes),
      sum(presence >= 0.5),
      sum(presence >= 0.75),
      sum(presence == 1)
    )
  )

  df_conserved$Category <- factor(
    df_conserved$Category,
    levels = c("AllGenes", "Shared50plus", "Shared75plus", "Shared100")
  )

  # plots!
  p_genes <- ggplot2::ggplot(summary_tbl, ggplot2::aes(Dataset, Genes, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.3) +
    ggplot2::theme_minimal()

  p_cells <- ggplot2::ggplot(summary_tbl,
                             ggplot2::aes(Dataset, Cells, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = Cells), vjust = -0.3) +
    ggplot2::theme_minimal()

  p_mito <- ggplot2::ggplot(qc_tbl, ggplot2::aes(Dataset, Mean_Mito, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal()

  p_feats <- ggplot2::ggplot(qc_tbl, ggplot2::aes(Dataset, Mean_Features, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal()

  p_conserved <- ggplot2::ggplot(df_conserved, ggplot2::aes(Category, Genes, fill = Category)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = Genes), vjust = -0.3) +
    ggplot2::theme_minimal()

  p_umi <- ggplot2::ggplot(qc_tbl,
                           ggplot2::aes(Dataset, Mean_UMI, fill = Dataset)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::theme_minimal()

  #outputs
  if (show_plots) {
    print(p_genes)
    print(p_cells)
    print(p_feats)
    print(p_mito)
    print(p_conserved)
    print(p_umi)
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
      conserved = p_conserved,
      cells = p_cells,
      umi = p_umi
    )
  )

}
