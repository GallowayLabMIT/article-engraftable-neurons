import(tidyverse)
import(pheatmap)
import(ggplot2)
import(tibble)

#' Normalizes short gene names to uniform format.
#'
#' The first letter of the gene short name is capitalized, with
#' all further letters in lower case.
#' @param x A list of gene names to normalize.
#' @return A list of normalized gene names
normalize_gene_name <- function(x) {
  substring(x, 1) <- toupper(substring(x, 1))
  substring(x, 2) <- tolower(substring(x, 2))
  x
}

#' Saves a gene clustering to a CSV for future reloading.
#' @param A CellDataSet containing full single-cell data.
#' @param A dataframe containing gene clustering information.
#' @param A filename to save the clustering information to
save_gene_clustering <- function(cds, gene_df, clustering_filename) {
  gene_lists <- gene_df[, c("id", "module")]
  gene_lists$short_name <- fData(cds)[gene_lists$id, ]$gene_short_name
  write.csv(gene_lists[order(gene_lists$module), ], clustering_filename)

}



normalize_heatmap_names <- function(mat) {
  result <- mat
  row.names(result) <- stringr::str_c("Cluster ", row.names(result))
  colnames(result) <- stringr::str_c(colnames(result))
  return(result)
}

#' Plots a heatmap of aggregate cell data
#' @param agg_expression A dataframe representing the aggregate expression values.
#' @param row_cluster An optional parameter specifying the gene cluster order.
#'    It can be returned by pheatmap$tree_row
#' @param ... Extra parameters passed to pheatmap::pheatmap. Common ones are
#'    'main' for the title and 'filename' to directly save the pheatmap
#' @return A pheatmap ggplot object
agg_heatmap <- function(agg_expression, row_cluster = TRUE, filename = NA, ...) {
  result <- pheatmap::pheatmap(normalize_heatmap_names(agg_expression),
                     cluster_rows = row_cluster, cluster_cols = FALSE,
                     scale = "column", clustering_method = "ward.D2",
                     fontsize = 10, ...)
  if (!is.na(filename)) {
    grDevices::dev.copy(grDevices::png, filename)
    grDevices::dev.off()
  }
  return(result)
}

#' Plots violin plots showing the distribution of gene cluster expression
#' across the entire cell data set.
violin_plots <- function(cds, gene_df, facet) {
  expression <- t(aggregate_gene_expression(
    cds, gene_df[, c("id", "module")]))
  cell_expression <- as.data.frame(pData(full_cds))
  cell_expression$cluster_expression <- expression[,1]
  ggplot(cell_expression) + aes_string(x=facet, y='cluster_expression') +
           geom_violin() + geom_boxplot(width=.1)
}
