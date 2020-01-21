#library(flowCore)
#library(FlowSOM)
#library(ConsensusClusterPlus)
#library(dplyr)
#library(RColorBrewer)
#library(Rtsne)
#library(umap)

##### Get som object from flowSOM
#' Get som object from flowSOM
#'
#' @param fcs_raw
#' @param use_markers
#'
#' @return
#' @import flowCore FlowSOM
#'
#' @examples
get_som <- function(fcs_raw, use_markers){
  fsom <- ReadInput(fcs_raw, transform = FALSE, scale = FALSE)
  set.seed(1234)
  som <- BuildSOM(fsom, colsToUse = use_markers)
  return(som)

}

##### Get mc consensusCluster object from ConsensusClusterPlus
#' Get mc consensusCluster object from ConsensusClusterPlus
#'
#' @param som
#' @param maxK
#'
#' @return
#' @import flowCore FlowSOM ConsensusClusterPlus
#'
#' @examples
get_consensusClust <- function(som, maxK = 20){
  mc <- ConsensusClusterPlus(t(som$map$codes), maxK = maxK, reps = 100,
                             pItem = 0.9, pFeature = 1, title = "consensus_plots", plot = "png",
                             clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  return(mc)
}

##### Get a optimal number of clusters
#' Get a optimal number of clusters
#'
#' @param mc
#' @param rate_var_expl
#'
#' @return
#'
#' @examples
get_optimal_clusters <- function(mc, rate_var_expl = 0.9){
  l <- sapply(2:length(mc), function(x) mean(mc[[x]]$ml))
  l <- sapply(2:length(l), function(x) l[x-1]-l[x])
  l <- l/sum(l)
  l <- sapply(1:length(l), function(x) sum(l[1:x]))
  optimum <- length(mc)
  for(i in length(l):1){
    if(l[i] <= rate_var_expl){
      optimum <- length(mc) - (length(l)-i)
      break}}
  return(optimum)
}

##### Create the cluster annotation separated by samples
#' Create the cluster annotation separated by samples
#'
#' @param fcs_raw
#' @param som
#' @param mc
#' @param k
#'
#' @return
#' @import flowCore FlowSOM ConsensusClusterPlus
#'
#' @examples
get_cluster_annotation <- function(fcs_raw, som, mc, k){
  code_clustering <- mc[[k]]$consensusClass
  cell_clustering <- code_clustering[som$map$mapping[,1]]
  l <- fsApply(fcs_raw, nrow)
  l <- c(0, sapply(1:length(l), function(x) sum(l[1:x])))
  cell_clustering_list <- sapply(2:length(l), function(x) cell_clustering[(l[x-1]+1):l[x]])
  names(cell_clustering_list) <- sampleNames(fcs_raw)
  return(cell_clustering_list)
}

##### Create the cluster annotation as one vector
#' Create the cluster annotation as one vector
#'
#' @param som
#' @param mc
#' @param k
#'
#' @return
#' @import flowCore FlowSOM ConsensusClusterPlus
#'
#' @examples
get_cell_clustering_list <- function(som, mc, k){
  code_clustering <- mc[[k]]$consensusClass
  cell_clustering <- code_clustering[som$map$mapping[,1]]
  return(cell_clustering)
}

##### Get euclidean distance between clones
#' Get euclidean distance between clones
#'
#' @param fcs_raw
#' @param use_markers
#' @param cell_clustering
#'
#' @return
#' @import flowCore dplyr
#'
#' @examples
get_euclid_dist <- function(fcs_raw, use_markers, cell_clustering){
  expr_median <- data.frame(fsApply(fcs_raw[,use_markers], exprs), cell_clustering = cell_clustering, check.names = F) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  rownames(expr_median) <- expr_median$cell_clustering
  expr_median$cell_clustering <- NULL
  cluster_euclidean_distance <- data.frame(t(combn(rownames(expr_median),2)), dist=as.matrix(dist(expr_median))[lower.tri(as.matrix(dist(expr_median)))] )
  colnames(cluster_euclidean_distance) <- c("Cluster_1", "Cluster_2", "euclidean_distance")
  return(cluster_euclidean_distance)
}

#####  Get edges
#' Get edges
#'
#' @param cluster_euclidean_distance
#'
#' @return
#'
#' @examples
get_edges <- function(cluster_euclidean_distance){
  edges <- data.frame(from = cluster_euclidean_distance$Cluster_1,
                      to = cluster_euclidean_distance$Cluster_2,
                      width  = as.numeric(log2(cluster_euclidean_distance$euclidean_distance)/2),
                      smooth = T)
  edges <- cluster_euclidean_distance[,1:3]
  colnames(edges) <- c('from', 'to', 'width')
  edges$id <- rownames(edges)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$width <- as.numeric(edges$width)
  return(edges)
}

##### Filter out edges which have weight more than threshold
#' Filter out edges which have weight more than threshold
#'
#' @param edges
#' @param threshold
#'
#' @return
#'
#' @examples
filter_edges <- function(edges, threshold){
  cut_threshold <-  min(edges$width) + ((max(edges$width) - min(edges$width)) * threshold)
  filtered_edges <- edges[edges$width <= cut_threshold, ]
  return(filtered_edges)
}

#####  Get nodes
#' Get nodes
#'
#' @param edges
#' @param cell_clustering
#'
#' @return
#' @import RColorBrewer
#'
#' @examples
get_nodes <- function(edges, cell_clustering){
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  cluster_colour <- sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))), length(unique(cell_clustering)))
  nodes <- data.frame(id = unique(c(edges$to, edges$from)),
                      color = cluster_colour,
                      label = unique(c(edges$to, edges$from)),
                      value = as.numeric(log(table(cell_clustering))),
                      title = paste0("Cluster ",  unique(c(edges$to, edges$from)),
                                     "<br>number of cells: ", table(cell_clustering)))
  return(nodes)
}

##### Get indexs of rows with were sampled to subset
#' Get indexs of rows with were sampled to subset
#'
#' @param fcs_raw
#' @param plot_ncell
#'
#' @return
#' @import flowCore
#'
#' @examples
get_inds_subset <- function(fcs_raw, sampling_size = 0.5){
  sample_ids <- rep(sampleNames(fcs_raw), fsApply(fcs_raw, nrow))
  inds <- split(1:length(sample_ids), sample_ids)
  #tsne_ncells <- pmin(table(sample_ids), sampling_size)
  tsne_ncells <- as.integer((table(sample_ids) + 1) * sampling_size)
  names(tsne_ncells) <- names(table(sample_ids))
  tsne_inds <- lapply(names(inds), function(i){s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)})
  tsne_inds <- unlist(tsne_inds)
  return(tsne_inds)
}

##### Get subseted dataframe for UMAP ploting
#' Title
#'
#' @param fcs_raw
#' @param use_markers
#' @param clust_markers
#' @param tsne_inds
#' @param cell_clustering
#'
#' @return
#' @import flowCore umap Rtsne
#'
#' @examples
get_UMAP_dataframe <- function(fcs_raw, use_markers, clust_markers, tsne_inds, cell_clustering, method = "UMAP"){
  expr <- fsApply(fcs_raw[,use_markers], exprs)
  tsne_expr <- expr[tsne_inds, clust_markers]
  if(method == "UMAP"){
    umap_out <- umap(tsne_expr)
    umap_df <- data.frame(expr[tsne_inds, use_markers], umap_out$layout, cluster =  as.factor(cell_clustering)[tsne_inds])
  }
  if(method == "tSNE"){
    set.seed(1234)
    umap_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)
    umap_df <- data.frame(expr[tsne_inds, use_markers], umap_out$Y, cluster =  as.factor(cell_clustering)[tsne_inds])
  }
  colnames(umap_df) <- c(names(use_markers), "UMAP_1", "UMAP_2", "cluster")
  return(umap_df)
}

##### Create the data table to draw the abundance barplot
#' Create the data table to draw the abundance barplot
#'
#' @param fcs_raw
#' @param cell_clustering
#'
#' @return
#' @import flowCore
#'
#' @examples
get_abundance_dataframe <- function(fcs_raw, cell_clustering){
  sample_ids <- rep(sampleNames(fcs_raw), fsApply(fcs_raw, nrow))
  abundance_data <- table(cell_clustering, sample_ids)
  abundance_data <- t(t(abundance_data) / colSums(abundance_data)) * 100
  abundance_data <- as.data.frame(abundance_data)
  colnames(abundance_data) <- c('cluster', 'sample_ids', 'abundance')
  return(abundance_data)
}

##### Merging two or more clusters within cluster annotating vector
#' Merging two or more clusters within cluster annotating vector
#'
#' @param clusters
#' @param cluster_to_merge
#'
#' @return
#' @import dplyr
#'
#' @examples
cluster_merging <- function(clusters, cluster_to_merge){
  new_clusters <- clusters
  new_clusters[clusters %in% cluster_to_merge] <- cluster_to_merge[1]
  return(new_clusters)
}

