#library(dplyr)
#library(flowCore)

##### Get dataframe for marker-cluster heatmap
#' Get dataframe for marker-cluster heatmap
#'
#' @param fcs_raw
#' @param use_markers
#' @param cell_clustering
#'
#' @return
#' @import flowCore dplyr
#'
#' @examples
get_mk_clust_heatmap_dataframe <- function(fcs_raw, use_markers, cell_clustering){
  cluster_expr_median <- data.frame(fsApply(fcs_raw[,use_markers], exprs), cell_clustering = cell_clustering, check.names = F) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  rownames(cluster_expr_median) <- cluster_expr_median$cell_clustering
  cluster_expr_median$cell_clustering <- NULL
  colnames(cluster_expr_median) <- names(use_markers)
  cluster_expr_median <- as.data.frame(cluster_expr_median)
  return(cluster_expr_median)
}

##### Get dataframe for deconvolution plot of markers to samples-clusters projection
#' Get dataframe for deconvolution plot of markers to samples-clusters projection
#'
#' @param fcs_raw
#' @param use_markers
#' @param cell_clustering
#'
#' @return
#' @import flowCore dplyr
#'
#' @examples
get_deconvol_dataframe <- function(fcs_raw, use_markers, cell_clustering){
  expr <- fsApply(fcs_raw[,use_markers], exprs)
  sample_ids <- rep(sampleNames(fcs_raw), fsApply(fcs_raw, nrow))
  deconv_expr <- data.frame(expr, cell_clustering = cell_clustering, sample_ids = sample_ids, check.names = F)
  sampl_cell_number <- deconv_expr %>% group_by(sample_ids) %>% summarise(cell_number = n())
  sampl_cell_number <- as.data.frame(sampl_cell_number)
  rownames(sampl_cell_number) <- sampl_cell_number$sample_ids
  cell_number <- deconv_expr %>% group_by(sample_ids, cell_clustering) %>%  summarise(cell_number = n())
  summarised_cluster <- deconv_expr %>% group_by(sample_ids, cell_clustering) %>%  summarise_all(median)
  summarised_cluster$cell_number <- cell_number$cell_number
  colnames(summarised_cluster)[match(use_markers, colnames(summarised_cluster))] <- names(use_markers)
  summarised_cluster$sample_cell_number <- sampl_cell_number[summarised_cluster$sample_ids, "cell_number"]
  summarised_cluster <- mutate(summarised_cluster, cell_rate = cell_number/sample_cell_number)
  summarised_cluster <- as.data.frame(summarised_cluster)
  return(summarised_cluster)
}
