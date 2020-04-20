
##### Matrix for plotting the content of panels
#' Matrix for plotting the content of panels
#'
#' @param md_set
#'
#' @return
#' @importFrom flowCore sampleNames
#'
#' @examples
get_panel_content <- function(fcs_raw_set){
  samples <- unique(unlist(lapply(fcs_raw_set, function(x) flowCore::sampleNames(x))))
  #samples <- unique(unlist(lapply(md_set, function(x) x[,"file_name"])))
  panel_content <- matrix(data = "absent", ncol = length(fcs_raw_set), nrow = length(samples))
  colnames(panel_content) <- names(fcs_raw_set)
  rownames(panel_content) <- samples
  for (p in colnames(panel_content)){
    panel_content[rownames(panel_content) %in% flowCore::sampleNames(fcs_raw_set[[p]]), p] <- "presented"
  }
  #panel_content <- reshape2::melt(panel_content)
  #colnames(panel_content) <- c('samples', 'panels', 'status')
  return(panel_content)
}

#get_panel_content(md_set)

##### Extracting of common markers across panels and formed a comparison table
#' Extracting of common markers across panels and formed a comparison table
#'
#' @param use_marker_set
#'
#' @return
#'
#' @examples
get_common_markers <- function(use_marker_set){
  mk <- unique(unlist(lapply(use_marker_set, names)))
  mk_panel_content <- matrix(data = "absent", ncol = length(use_marker_set), nrow = length(mk))
  colnames(mk_panel_content) <- names(use_marker_set)
  rownames(mk_panel_content) <- mk
  for (p in colnames(mk_panel_content)){
    mk_panel_content[rownames(mk_panel_content) %in% names(use_marker_set[[p]]), p] <- "presented"
  }
  return(mk_panel_content)
}

#get_common_markers(use_marker_set)

##### Get expression data for one marker
get_mk_exp_data_cross_p <- function(fcs_raw_set, use_samples, use_marker_set, target_mk){
  mk_exp_data_cross_p <- unlist(lapply(names(fcs_raw_set), function(p) {
    lapply(use_samples, function(s) {
      flowCore::exprs(fcs_raw_set[[p]][[s]])[,use_marker_set[[p]][target_mk]]
    })
  }))
  return(mk_exp_data_cross_p)
}

##### Compiling the matrix with cluster abundance data for all panels
#' Compiling the matrix with cluster abundance data for all panels
#'
#' @param cell_clustering_list_set
#' @param use_samples
#'
#' @return
#'
#' @examples
get_abund_cross_panel <- function(cell_clustering_list_set, use_samples){
  abund_cross_p <- do.call(cbind, lapply(names(cell_clustering_list_set), function(p) {
    p_clusters <- unique(unlist(lapply(use_samples, function(s) unique(unlist(cell_clustering_list_set[[p]][s])))))
    tmp_abund <- matrix(data = 0, nrow = length(use_samples), ncol = length(p_clusters))
    colnames(tmp_abund) <- p_clusters
    rownames(tmp_abund) <- use_samples
    for(s in use_samples){
      portion_data <- table(cell_clustering_list_set[[p]][[s]])/length(cell_clustering_list_set[[p]][[s]])
      tmp_abund[s,] <- portion_data[p_clusters]
    }
    colnames(tmp_abund) <- paste(p, colnames(tmp_abund), sep = "-")
    return(tmp_abund)
  }))
  return(abund_cross_p)
}

#abund_cross_p <- get_abund_corr_cross_panel(cell_clustering_list_set, use_samples)

##### Produce three matrices with corr coefficient, p-value and adjusted p-value for abundance correlations within cross panel assay
#' Produce three matrices with corr coefficient, p-value and adjusted p-value for abundance correlations within cross panel assay
#'
#' @param abund_cross_p
#' @param method
#'
#' @return

#'
#' @examples
get_matrices_corr_info <- function(abund_cross_p, method_cortest = "spearman", method_padj = "BH"){
  corr_martix <- matrix(0, ncol = ncol(abund_cross_p), nrow = ncol(abund_cross_p))
  colnames(corr_martix) <- colnames(abund_cross_p)
  rownames(corr_martix) <- colnames(abund_cross_p)
  p_val_martix <- matrix(0, ncol = ncol(abund_cross_p), nrow = ncol(abund_cross_p))
  colnames(p_val_martix) <- colnames(abund_cross_p)
  rownames(p_val_martix) <- colnames(abund_cross_p)
  padj_martix <- matrix(0, ncol = ncol(abund_cross_p), nrow = ncol(abund_cross_p))
  colnames(padj_martix) <- colnames(abund_cross_p)
  rownames(padj_martix) <- colnames(abund_cross_p)
  for (i in colnames(abund_cross_p)){
    for (j in colnames(abund_cross_p)){
      tmp <- stats::cor.test(abund_cross_p[,i], abund_cross_p[,j], method = method_cortest)
      corr_martix[i,j] <- tmp$estimate
      p_val_martix[i,j] <- tmp$p.value
      if(is.na(corr_martix[i,j])){corr_martix[i,j] <- 0}
      if(is.na(p_val_martix[i,j])){p_val_martix[i,j] <- 1}
    }
  }
  padj_martix <- apply(p_val_martix, 2, function(x) stats::p.adjust(x,method = method_padj))
  corr_info <- list(corr_martix, p_val_martix, padj_martix)
  names(corr_info) <- c("corr_coef", "p_value", "padj")
  return(corr_info)
}

#corr_info <- get_matrices_corr_info(abund_cross_p)

##### Get expressing cell fraction data for target marker
#' Get dxpressing cell fraction data for target marker
#'
#' @param fcs_raw_set
#' @param cell_clustering_list_set
#' @param use_samples
#' @param use_marker_set
#' @param target_mk
#' @param min_cell_number
#'
#' @return
#' @importFrom flowCore exprs
#' @importFrom stats kmeans
#'
#' @examples
get_clustering_exp_cell_f_cross_panel <- function(fcs_raw_set, cell_clustering_list_set, use_samples, use_marker_set, target_mk,
                                       min_cell_number = 5){
  exp_cell_f_cross_p <- do.call(cbind, lapply(names(fcs_raw_set), function(p) {
    p_clusters <- unique(unlist(lapply(use_samples, function(s) unique(unlist(cell_clustering_list_set[[p]][s])))))
    tmp_exp_cell_f <- matrix(data = 0, nrow = length(use_samples), ncol = length(p_clusters))
    colnames(tmp_exp_cell_f) <- p_clusters
    rownames(tmp_exp_cell_f) <- use_samples
    for (s in use_samples) {
      for (cl in p_clusters){
        exp_data <- flowCore::exprs(fcs_raw_set[[p]][[s]])[cell_clustering_list_set[[p]][[s]] == cl , use_marker_set[[p]][target_mk]]
        #set.seed(123)
        if(sum(exp_data != 0) <= min_cell_number){next}
        tmp <- stats::kmeans(exp_data, centers = 2, nstart = 25)
        tmp_exp_cell_f[s,cl] <- sum(tmp$cluster == which.max(tmp$centers))/length(tmp$cluster)
      }
    }
    colnames(tmp_exp_cell_f) <- paste(p, colnames(tmp_exp_cell_f), sep = "-")
    return(tmp_exp_cell_f)
  }))
  return(exp_cell_f_cross_p)
}

#get_exp_cell_f_cross_panel(fcs_raw_set, cell_clustering_list_set, use_samples, use_marker_set, target_mk)

##### Extract expression data of one marker for coss-panel analysis as one vector
#' Extract expression data of one marker for coss-panel analysis as one vector
#'
#' @param fcs_raw_set
#' @param cell_clustering_list_set
#' @param use_samples
#' @param use_marker_set
#' @param target_mk
#' @param threshold
#' @param min_cell_number
#'
#' @return
#' @importFrom flowCore exprs
#'
#' @examples
get_threshold_exp_cell_f_cross_panel <- function(fcs_raw_set, cell_clustering_list_set, use_samples, use_marker_set, target_mk,
                                                  threshold = 20,  min_cell_number = 5){
  exp_cell_f_cross_p <- do.call(cbind, lapply(names(fcs_raw_set), function(p) {
    p_clusters <- unique(unlist(lapply(use_samples, function(s) unique(unlist(cell_clustering_list_set[[p]][s])))))
    tmp_exp_cell_f <- matrix(data = 0, nrow = length(use_samples), ncol = length(p_clusters))
    colnames(tmp_exp_cell_f) <- p_clusters
    rownames(tmp_exp_cell_f) <- use_samples
    for (s in use_samples) {
      for (cl in p_clusters){
        exp_data <- flowCore::exprs(fcs_raw_set[[p]][[s]])[cell_clustering_list_set[[p]][[s]] == cl , use_marker_set[[p]][target_mk]]
        if(sum(exp_data != 0) <= min_cell_number){next}
        tmp_exp_cell_f[s,cl] <- sum(exp_data >= threshold)/length(exp_data)
      }
    }
    colnames(tmp_exp_cell_f) <- paste(p, colnames(tmp_exp_cell_f), sep = "-")
    return(tmp_exp_cell_f)
  }))
  return(exp_cell_f_cross_p)
}

##### Preparing the table with cross-panel correlations for downloading

#' Preparing the table with cross-panel correlations for downloading
#'
#' @param abund_corr_info
#'
#' @return
#' @importFrom reshape2 melt
#'
#' @examples
get_write_corr_info_cross_p <- function(abund_corr_info){
  write_corr_info_cross_p <- reshape2::melt(abund_corr_info$corr_coef)
  rownames(write_corr_info_cross_p) <- paste0(as.character(write_corr_info_cross_p[,1]), "_vs_", as.character(write_corr_info_cross_p[,2]))
  tmp_pval <- reshape2::melt(abund_corr_info$p_value)
  rownames(tmp_pval) <- paste0(as.character(tmp_pval[,1]), "_vs_", as.character(tmp_pval[,2]))
  tmp_pval <- tmp_pval[rownames(write_corr_info_cross_p), 3]
  tmp_padj <- reshape2::melt(abund_corr_info$padj)
  rownames(tmp_padj) <- paste0(as.character(tmp_padj[,1]), "_vs_", as.character(tmp_padj[,2]))
  tmp_padj <- tmp_padj[rownames(write_corr_info_cross_p), 3]
  write_corr_info_cross_p <- cbind(write_corr_info_cross_p, tmp_pval, tmp_padj)
  colnames(write_corr_info_cross_p) <- c("cluster_1", "cluster_2", "corr_coef", "p_value", "padj")
  return(write_corr_info_cross_p)
}

#get_write_corr_info_cross_p(abund_corr_info)
