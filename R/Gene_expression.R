
##### Get dataframe for marker-cluster heatmap
#' Get dataframe for marker-cluster heatmap
#'
#' @param fcs_raw
#' @param use_markers
#' @param cell_clustering
#'
#' @return
#' @importFrom flowCore fsApply exprs "exprs<-"
#' @importFrom dplyr group_by summarize_all funs
#' @importFrom stats median
#'
#' @examples
get_mk_clust_heatmap_dataframe <- function(fcs_raw, use_markers, cell_clustering, summarize_mode = 'median'){
  cluster_expr_median <- data.frame(flowCore::fsApply(fcs_raw[,use_markers], flowCore::exprs),
                                    cell_clustering = cell_clustering, check.names = F)
  cluster_expr_median <- dplyr::group_by(cluster_expr_median, cell_clustering)
  if(summarize_mode == 'median'){cluster_expr_median <- dplyr::summarize_all(cluster_expr_median, dplyr::funs(stats::median))}
  if(summarize_mode == 'mean'){cluster_expr_median <- dplyr::summarize_all(cluster_expr_median, dplyr::funs(mean))}
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
#' @importFrom magrittr "%>%"
#' @importFrom flowCore fsApply exprs "exprs<-" sampleNames "sampleNames<-"
#' @importFrom dplyr group_by summarise n
#' @importFrom stats median
#'
#' @examples
get_deconvol_dataframe <- function(fcs_raw, use_markers, cell_clustering, summarize_mode = 'median'){
  expr <- flowCore::fsApply(fcs_raw[,use_markers], flowCore::exprs)
  sample_ids <- rep(flowCore::sampleNames(fcs_raw), flowCore::fsApply(fcs_raw, nrow))
  deconv_expr <- data.frame(expr, cell_clustering = cell_clustering, sample_ids = sample_ids, check.names = F)
  sampl_cell_number <- deconv_expr %>% dplyr::group_by(sample_ids) %>% dplyr::summarise(cell_number = dplyr::n())
  sampl_cell_number <- as.data.frame(sampl_cell_number)
  rownames(sampl_cell_number) <- sampl_cell_number$sample_ids
  cell_number <- deconv_expr %>% dplyr::group_by(sample_ids, cell_clustering) %>%  dplyr::summarise(cell_number = dplyr::n())
  summarised_cluster <- deconv_expr %>% dplyr::group_by(sample_ids, cell_clustering)
  if(summarize_mode == 'median'){summarised_cluster <- dplyr::summarise_all(summarised_cluster, stats::median)}
  if(summarize_mode == 'mean'){summarised_cluster <- dplyr::summarise_all(summarised_cluster, mean)}
  summarised_cluster$cell_number <- cell_number$cell_number
  colnames(summarised_cluster)[match(use_markers, colnames(summarised_cluster))] <- names(use_markers)
  summarised_cluster$sample_cell_number <- sampl_cell_number[summarised_cluster$sample_ids, "cell_number"]
  summarised_cluster <- dplyr::mutate(summarised_cluster, cell_rate = cell_number/sample_cell_number)
  summarised_cluster <- as.data.frame(summarised_cluster)
  return(summarised_cluster)
}
