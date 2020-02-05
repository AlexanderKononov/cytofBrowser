#library(flowCore)

### Create test data
#test_fcs_raw <- fsApply(fcs_raw, function(x) x[1:300,])
#test_cell_clustering_list <- lapply(cell_clustering_list, function(x) x[1:100])
#table(test_cell_clustering_list[[1]])

###############################
####  Up-levele functions  ####
###############################

##### Get list_cell_ctDist object to correlation analysis
#' Get list_cell_ctDist object to correlation analysis
#'
#' @param cell_clustering_list
#'
#' @return
#'
#' @examples
get_list_cell_ctDist <- function(cell_clustering_list){
  cell_clustering_list <- as.list(cell_clustering_list)
  list_coldata_sce <- lapply(1:length(cell_clustering_list), function(x) data.frame(call = 1:length(cell_clustering_list[[x]]), cluster = cell_clustering_list[[x]]))
  list_cell_ctDist <- ct_sample_aggregator(list_coldata_sce, cell_type_column = "cluster")
  return(list_cell_ctDist)
}

###test for get_list_cell_ctDist
#test_list_cell_ctDist <- get_list_cell_ctDist(test_cell_clustering_list)
#table(test_cell_clustering_list[[1]])
#colSums(test_list_cell_ctDist[[1]])

### Functions to draw abundance correlation triangle plot
#########################################################

##### Get data of cluster abundace
#' Get data of cluster abundace
#'
#' @param list_cell_ctDist
#'
#' @return
#'
#' @examples
get_abundance <- function(list_cell_ctDist){
  cell_types <- unique(unlist(lapply(list_cell_ctDist, colnames)))
  abundence_data <- sapply(cell_types, function(x) abundence_meter(x, list_cell_ctDist))
  return(abundence_data)
}

### test
#test_abundence_data <- get_abundance(test_list_cell_ctDist)
#table(test_cell_clustering_list[[4]])/length(test_cell_clustering_list[[4]])
#test_abundence_data[4,]

##### Calculation of correlation between cluster by abundance data
#' Calculation of correlation between cluster by abundance data
#'
#' @param abundence_data
#' @param method
#'
#' @return
#'
#' @examples
get_abundance_correlation <- function(abundence_data, method = "spearman"){
  abundance_correlation <- do.call(rbind, lapply(colnames(abundence_data), function(i){
    do.call(rbind, lapply(colnames(abundence_data), function(j){
      res <- cor.test(abundence_data[,i], abundence_data[,j], method = method)
      data.frame(cluster_1 = i, cluster_2 = j, cor_significance = res$p.value, cor_coefficient = res$estimate)
    }))
  }))
  abundance_correlation <- abundance_correlation[abundance_correlation$cluster_1!=abundance_correlation$cluster_2,]
  rownames(abundance_correlation) <- NULL
  abundance_correlation$adj_pValue <- p.adjust(abundance_correlation$cor_significance, method = "BH")
  return(abundance_correlation)
}

###test
#test_abundance_correlation <- get_abundance_correlation(test_abundence_data)
#length(colnames(test_abundence_data))^2-length(colnames(test_abundence_data)) == nrow(test_abundance_correlation)

##### Extraction correlation coefficient from abundance correlation data
#' Extraction correlation coefficient from abundance correlation data
#'
#' @param abundance_correlation
#'
#' @return
#' @import reshape2
#'
#' @examples
get_corr_coef_matrix <- function(abundance_correlation){
  corr_coef_matrix <- dcast(abundance_correlation, cluster_1 ~ cluster_2, value.var = "cor_coefficient")
  rownames(corr_coef_matrix) <- corr_coef_matrix[,1]
  corr_coef_matrix[,1] <-NULL
  corr_coef_matrix <- as.matrix(corr_coef_matrix)
  corr_coef_matrix[is.na(corr_coef_matrix)] <- 1
  return(corr_coef_matrix)
}

##### Extraction p-value from abundance correlation data
#' Extraction p-value from abundance correlation data
#'
#' @param abundance_correlation
#'
#' @return
#' @import reshape2
#'
#' @examples
get_corr_pValue_matrix <- function(abundance_correlation){
  corr_pValue_matrix <- dcast(abundance_correlation, cluster_1 ~ cluster_2, value.var = "cor_significance")
  rownames(corr_pValue_matrix) <- corr_pValue_matrix[,1]
  corr_pValue_matrix[,1] <-NULL
  corr_pValue_matrix <- as.matrix(corr_pValue_matrix)
  corr_pValue_matrix[is.na(corr_pValue_matrix)] <- 0
  return(corr_pValue_matrix)
}

##### Extraction adjusted p-value from abundance correlation data
#' Extraction adjusted p-value from abundance correlation data
#'
#' @param abundance_correlation
#'
#' @return
#' @import reshape2
#'
#' @examples
get_adj_corr_pValue_matrix <- function(abundance_correlation){
  adj_corr_pValue_matrix <- dcast(abundance_correlation, cluster_1 ~ cluster_2, value.var = "adj_pValue")
  rownames(adj_corr_pValue_matrix) <- adj_corr_pValue_matrix[,1]
  adj_corr_pValue_matrix[,1] <-NULL
  adj_corr_pValue_matrix <- as.matrix(adj_corr_pValue_matrix)
  adj_corr_pValue_matrix[is.na(adj_corr_pValue_matrix)] <- 0
  return(adj_corr_pValue_matrix)
}

### Functions for complete correlation analysis
###############################################

##### Get list_expData object to correlation analysis
#' Get list_expData object to correlation analysis
#'
#' @param fcs_raw
#'
#' @return
#' @import flowCore
#'
#' @examples
get_list_expData <- function(fcs_raw){
  list_expData <- lapply(fsApply(fcs_raw, function(x) as.data.frame(exprs(x))),t)
  return(list_expData)
}

###test
#test_list_expData <- get_list_expData(test_fcs_raw)

##### Convert from signal dataframe to tidy dataframe with correlations between different clusters
#' Convert from signal dataframe to tidy dataframe with correlations between different clusters
#'
#' @param signals
#'
#' @return
#'
#' @examples
get_signals_between_clusters <- function(signals){
  signals_between_clusters <- do.call(rbind, lapply(names(signals), function(i){
    do.call(rbind, lapply(names(signals[[i]]), function(j){
      if(nrow(signals[[i]][[j]])==0){return(NULL)}
      data.frame(signaling_cluster = i, targetet_cluster = j,
                 gene_in_target_cluster = rownames(signals[[i]][[j]]),
                 signals[[i]][[j]][,c("cor_significance", "cor_coefficient")])
    }))
  }))
  rownames(signals_between_clusters) <- NULL
  signals_between_clusters$adj_pValue <- p.adjust(signals_between_clusters$cor_significance, method = "BH")
  return(signals_between_clusters)
}

###test
#test_signals_between_clusters <- get_signals_between_clusters(test_signals)

##### Convert from signal dataframe to tidy dataframe with correlations within one cluster
#' Convert from signal dataframe to tidy dataframe with correlations within one cluster
#'
#' @param signals
#'
#' @return
#'
#' @examples
get_signals_in_cluster <- function(signals){
  signals_in_cluster <- do.call(rbind, lapply(names(signals), function(i){
    if(nrow(signals[[i]][[i]])==0){return(NULL)}
    data.frame(cluster = i, gene = rownames(signals[[i]][[i]]),
               signals[[i]][[i]][,c("cor_significance", "cor_coefficient")])
  }))
  rownames(signals_in_cluster) <- NULL
  signals_in_cluster$adj_pValue <- p.adjust(signals_in_cluster$cor_significance, method = "BH")
  return(signals_in_cluster)
}

###test
#test_signals_in_cluster <- get_signals_in_cluster(test_signals)

### Low-level function for extraction of top significant corrlation signals
#' Low-level function for extraction of top significant corrlation signals
#'
#' @param cor_tdata
#'
#' @return
#'
#' @examples
take_top_correlation <- function(cor_tdata, adj_pValue = 0.01){
  a <- cor_tdata
  a <- a[!is.na(a$cor_coefficient),]
  a <- a[a$adj_pValue <= 0.01,]
  a <- a[order(a$cor_coefficient),]
  return(a)
}

##### Get the top of significant signals from signals_between_clusters object
#' Get the top of significant signals from signals_between_clusters object
#'
#' @param signals_between_clusters
#' @param use_markers
#'
#' @return
#' @import dplyr
#'
#' @examples
get_signals_between_clusters_top <- function(signals_between_clusters, use_markers){
  signals_between_clusters_top <- take_top_correlation(signals_between_clusters)
  signals_between_clusters_top <- signals_between_clusters_top[signals_between_clusters_top$gene_in_target_cluster %in% use_markers,]
  signals_between_clusters_top$gene_in_target_cluster <- names(use_markers)[match(signals_between_clusters_top$gene_in_target_cluster, use_markers)]
  signals_between_clusters_top$names <- paste0(signals_between_clusters_top$signaling_cluster, "-",
                                               signals_between_clusters_top$targetet_cluster, "_",
                                               signals_between_clusters_top$gene_in_target_cluster)
  rownames(signals_between_clusters_top) <- NULL
  return(signals_between_clusters_top)
}

##### Get the top of significant signals from signals_between_clusters object
#' Get the top of significant signals from signals_between_clusters object
#'
#' @param signals_in_cluster_top
#' @param use_markers
#'
#' @return
#' @import dplyr
#'
#' @examples
get_signals_in_cluster_top <- function(signals_in_cluster_top, use_markers){
  signals_in_cluster_top <- take_top_correlation(signals_in_cluster)
  signals_in_cluster_top <- signals_in_cluster_top[signals_in_cluster_top$gene %in% use_markers,]
  signals_in_cluster_top$gene <- names(use_markers)[match(signals_in_cluster_top$gene, use_markers)]
  signals_in_cluster_top$names <- paste0(signals_in_cluster_top$cluster, "_", signals_in_cluster_top$gene)
  rownames(signals_in_cluster_top) <- NULL
  return(signals_in_cluster_top)
}


###############################
#### Low-levele functions  ####
###############################

##### Functions to organise initial lists of data in a suitable form
####################################################################

### Low-level function to get cell_ctDist object for one sample
#' Low-level function to get cell_ctDist object for one sample
#'
#' @param coldata_sce
#' @param cell_id_column
#' @param cell_type_column
#'
#' @return
#'
#' @examples
cell_type_organiser <- function(coldata_sce, cell_id_column = "rownames", cell_type_column = "cluster"){
  if(cell_id_column == "rownames"){cell_id_column <- rownames(coldata_sce)}
  cell_types <- sort(unique(coldata_sce[,cell_type_column]))
  cell_ctDist <- matrix(data = F, nrow = length(cell_id_column), ncol = length(cell_types))
  rownames(cell_ctDist) <-as.character(cell_id_column)
  colnames(cell_ctDist) <- as.character(cell_types)
  for(i in cell_types){
    cell_ctDist[as.character(cell_id_column[coldata_sce[,cell_type_column]==i]),as.character(i)] <- T}
  return(cell_ctDist)
}

### Low-level function for unification of clusters set between samples within list_cell_ctDist object
#' Low-level function for unification of clusters set between samples within list_cell_ctDist object
#'
#' @param list_cell_ctDist
#'
#' @return
#'
#' @examples
cluster_unification <- function(list_cell_ctDist){
  cell_types <- unique(unlist(sapply(list_cell_ctDist, colnames)))
  new_list_cell_ctDist <- lapply(list_cell_ctDist, function(x){
    add_clusters <- setdiff(cell_types, colnames(x))
    add_ctDist <- matrix(data = F,  nrow = nrow(x), ncol = length(add_clusters))
    colnames(add_ctDist) <- add_clusters
    cbind(x, add_ctDist)
  })
  return(new_list_cell_ctDist)
}

##### Get list_cell_ctDist object for correlation analysis
#' Get list_cell_ctDist object for correlation analysis
#'
#' @param list_coldata_sce
#' @param cell_id_column
#' @param cell_type_column
#' @param do_cluster_unification
#'
#' @return
#'
#' @examples
ct_sample_aggregator <- function(list_coldata_sce, cell_id_column = "rownames", cell_type_column = "cluster", do_cluster_unification = TRUE){
  list_cell_ctDist <- lapply(list_coldata_sce, function(x) cell_type_organiser(x, cell_id_column = cell_id_column, cell_type_column = cell_type_column))
  if(do_cluster_unification){list_cell_ctDist <- cluster_unification(list_cell_ctDist)}
  return(list_cell_ctDist)
}

### Low-level function to calculate abundances for each cluster within each sample
#' Low-level function to calculate abundances for each cluster within each sample
#'
#' @param signalling_type
#' @param list_cell_ctDist
#'
#' @return
#'
#' @examples
abundence_meter <- function(signalling_type, list_cell_ctDist){
  abund_data <- sapply(list_cell_ctDist, function(x) sum(x[,signalling_type])/nrow(x))
  return(abund_data)
}

### Low-level function to get tt_expData object for one sample
#' Low-level function to get tt_expData object for one sample
#'
#' @param cell_ctDist
#' @param expData
#'
#' @return
#'
#' @examples
target_type_organiser <- function(cell_ctDist, expData){
  tt_expData <- lapply(colnames(cell_ctDist), function(x) as.matrix(expData[,cell_ctDist[,x]]))
  names(tt_expData) <- colnames(cell_ctDist)
  return(tt_expData)
}

##### Get list_tt_expData object for correlation analysis
#' Get list_tt_expData object for correlation analysis
#'
#' @param list_cell_ctDist
#' @param list_expData
#'
#' @return
#'
#' @examples
tt_sample_aggregator <- function(list_cell_ctDist, list_expData){
  list_tt_expData <- lapply(1:length(list_cell_ctDist), function(x) target_type_organiser(list_cell_ctDist[[x]], list_expData[[x]]))
  names(list_tt_expData) <- names(list_expData)
  return(list_tt_expData)
}

###test
#test_list_tt_expData <- tt_sample_aggregator(test_list_cell_ctDist, test_list_expData)

##### Functions to create background contrast objects used for correlation analysis
###################################################################################

##### Get bg_ctCor_data object with background correlation statistics
# Time-concuming (~30min)
#' Get bg_ctCor_data object with background correlation statistics
#'
#' @param list_cell_ctDist
#' @param list_expData
#' @param method
#'
#' @return
#'
#' @examples
backgraund_correlation <- function(list_cell_ctDist, list_expData, method = "spearman"){
  bg_htest <- matrix(data = NA, nrow = nrow(list_expData[[1]]), ncol = ncol(list_cell_ctDist[[1]]))
  rownames(bg_htest) <- rownames(list_expData[[1]])
  colnames(bg_htest) <- colnames(list_cell_ctDist[[1]])
  bg_coef <- matrix(data = NA, nrow = nrow(list_expData[[1]]), ncol = ncol(list_cell_ctDist[[1]]))
  rownames(bg_coef) <- rownames(list_expData[[1]])
  colnames(bg_coef) <- colnames(list_cell_ctDist[[1]])

  for(i in rownames(list_expData[[1]])){
    comparison_unit <- do.call(rbind, lapply(list_expData, function(x) data.frame(expression = x[i,])))
    comparison_unit <- cbind(comparison_unit, do.call(cbind, lapply(colnames(list_cell_ctDist[[1]]), function(x)
      rep(abundence_meter(x, list_cell_ctDist), time = sapply(list_expData, ncol)))))
    colnames(comparison_unit) <- c("expression",colnames(list_cell_ctDist[[1]]))
    htest <- sapply(comparison_unit[,-1], cor.test, y = comparison_unit$expression)
    bg_htest[i,] <- unlist(htest["p.value",])
    bg_coef[i,] <- unlist(htest["estimate",])
  }
  bg_padj <- apply(bg_htest, 2, function(x) p.adjust(x, method = "BH"))
  bg_ctCor_data <- list(bg_htest, bg_coef, bg_padj)
  names(bg_ctCor_data) <- c("bg_htest", "bg_coef", "bg_padj")
  return(bg_ctCor_data)
}

###test
#test_bg_ctCor_data <- backgraund_correlation(test_list_cell_ctDist, test_list_expData)

##### Get bg_anova object with background ANOVA statistics
# Time-concuming (~40min)
#' Get bg_anova object with background ANOVA statistics
#'
#' @param list_cell_ctDist
#' @param list_expData
#'
#' @return
#'
#' @examples
background_anova <- function(list_cell_ctDist ,list_expData){
  bg_anova <- do.call(cbind, lapply(colnames(list_cell_ctDist[[1]]), function(j){
    tt_anova_bg <- sapply(rownames(list_expData[[1]]), function(i){
      comparison_unit <- do.call(rbind, lapply(1:length(list_expData), function(x){
        if(length(list_expData[[x]][i,list_cell_ctDist[[x]][,j]])==0){return(NA)}
        data.frame(expression = list_expData[[x]][i,list_cell_ctDist[[x]][,j]], samples = x)}))
      comparison_unit$samples <- as.factor(comparison_unit$samples)
      if(length(levels(comparison_unit$sample))<2){return(NA)}
      aov_data <- aov(formula = expression ~ samples, data = comparison_unit)
      return(summary(aov_data)[[1]][["Pr(>F)"]][1])
    })
    return(tt_anova_bg)
  }))
  colnames(bg_anova) <- colnames(list_cell_ctDist[[1]])
  bg_anova_adj <- apply(bg_anova, 2, function(x){
    p_adj <- p.adjust(x, method = "BH")
    return(p_adj)})
  colnames(bg_anova_adj) <- colnames(list_cell_ctDist[[1]])
  bg_anova_result <- list(bg_anova, bg_anova_adj)
  names(bg_anova_result) <- c("bg_anova", "bg_anova_adj")
  return(bg_anova_result)
}

###test
#test_bg_anova <- background_anova(test_list_cell_ctDist, test_list_expData)

##### Functions to extract correlation signals different from the background
############################################################################

### Low-level function for correlation analysis of genes from one cluster with abundance of another cluster
# Time-concuming (~5min)
#' Low-level function for correlation analysis of genes from one cluster with abundance of another cluster
#'
#' @param target_type
#' @param abund_data
#' @param list_tt_expData
#' @param method
#'
#' @return
#'
#' @examples
target_type_cortest <- function(target_type, abund_data, list_tt_expData, method = "spearman"){
  tt_htest <- lapply(rownames(list_tt_expData[[1]][[1]]), function(i){
    comparison_unit <- lapply(1:length(abund_data), function(x) {
      if(length(list_tt_expData[[x]][[target_type]])==0){return(NULL)}
      tmp <- data.frame(expression = list_tt_expData[[x]][[target_type]][i,])
      tmp$abundance <- abund_data[x]
      return(as.data.frame(tmp))})
    comparison_unit <- do.call(rbind, comparison_unit)
    if(nrow(comparison_unit) < 3){return(NA)}
    return(cor.test(comparison_unit$expression, comparison_unit$abundance, method = method))
  })
  names(tt_htest) <- rownames(list_tt_expData[[1]][[1]])
  return(tt_htest)
}

###test
#target_type <- names(test_list_tt_expData[[1]])[1]
#test_abund_data <- abundence_meter( colnames(test_list_cell_ctDist[[1]])[1], test_list_cell_ctDist)
#test_tt_htest <- target_type_cortest(target_type = target_type , abund_data = test_abund_data, list_tt_expData = test_list_tt_expData)

### Low-level function for comparison of detected correlations between clusters with background correlations
#' Low-level function for comparison of detected correlations between clusters with background correlations
#'
#' @param gene_ctCor_data
#' @param bg_ctCor_data
#' @param bg_anova
#' @param bg_anova_alpha
#' @param bg_ctCor_alpha
#' @param gene_ctCor_alpha
#'
#' @return
#'
#' @examples
cell_filter <- function(gene_ctCor_data, bg_ctCor_data, bg_anova,
                        bg_anova_alpha = 0.01, bg_ctCor_alpha = 0.01, gene_ctCor_alpha = 0.01){
  filter_aggregator <- matrix(data = T, nrow = nrow(gene_ctCor_data[[3]]), ncol = ncol(gene_ctCor_data[[3]])  )
  rownames(filter_aggregator) <- rownames(gene_ctCor_data[[3]])
  colnames(filter_aggregator) <- colnames(gene_ctCor_data[[3]])
  for(i in colnames(bg_anova[[2]])){filter_aggregator[(bg_anova[[2]][,i]>=bg_anova_alpha),i] <- F}
  for(i in colnames(bg_ctCor_data[[3]])){
    filter_aggregator[(bg_ctCor_data[[3]][,i] <= bg_ctCor_alpha) & ((gene_ctCor_data[[2]][,i] * bg_ctCor_data[[2]][,i]) > 0),i] <- F
  }
  gene_ctDistr <- (gene_ctCor_data[[3]] <= gene_ctCor_alpha) & filter_aggregator
  rownames(gene_ctDistr) <- rownames(gene_ctCor_data[[3]])
  colnames(gene_ctDistr) <- colnames(gene_ctCor_data[[3]])
  return(gene_ctDistr)
}

##test
#test_gene_ctDistr <- cell_filter(test_gene_ctCor_data, test_bg_ctCor_data, test_bg_anova)

### Low-level function to format detected signals to result dataframe
#' Low-level function to format detected signals to result dataframe
#'
#' @param gene_ctDistr
#' @param bg_ctCor_data
#'
#' @return
#'
#' @examples
target_type_data_extractor <- function(gene_ctDistr, bg_ctCor_data){
  list_gene_ttStat <- lapply(colnames(gene_ctDistr), function(x){
    gene_tt_Stat <- as.data.frame(cbind(gene_ctDistr[,x], bg_ctCor_data[[1]][,x], bg_ctCor_data[[2]][,x]))
    colnames(gene_tt_Stat) <- c("filter_pass", "cor_significance", "cor_coefficient")
    return(gene_tt_Stat)
  })
  names(list_gene_ttStat) <- colnames(gene_ctDistr)
  return(list_gene_ttStat)
}


### Low-level function for correlation analysis of genes from all clusters with abundance of one cluster
# Time-concuming (~15min)
#' Low-level function for correlation analysis of genes from all clusters with abundance of one cluster
#'
#' @param abund_data
#' @param list_tt_expData
#' @param method
#'
#' @return
#'
#' @examples
htest_data_extractor <- function(abund_data, list_tt_expData, method = "spearman"){
  gene_ctSignif <- matrix(data = NA, nrow = nrow(list_tt_expData[[1]][[1]]), ncol = length(list_tt_expData[[1]]))
  rownames(gene_ctSignif) <- rownames(list_tt_expData[[1]][[1]])
  colnames(gene_ctSignif) <- names(list_tt_expData[[1]])
  ct_coef <- matrix(data = NA, nrow = nrow(list_tt_expData[[1]][[1]]), ncol = length(list_tt_expData[[1]]))
  rownames(ct_coef) <- rownames(list_tt_expData[[1]][[1]])
  colnames(ct_coef) <- names(list_tt_expData[[1]])
  for(i in names(list_tt_expData[[1]])){
    tt_htest <- target_type_cortest(target_type = i , abund_data = abund_data, list_tt_expData = list_tt_expData, method = method)
    if(all(is.na(tt_htest))){
      gene_ctSignif[,i] <- NA
      ct_coef[,i] <- NA
      next}
    gene_ctSignif[,i] <- sapply(names(tt_htest), function(x) tt_htest[[x]]$p.value)
    ct_coef[,i] <- sapply(names(tt_htest), function(x) tt_htest[[x]]$estimate)
  }
  gene_ctSignif_adj <- apply(gene_ctSignif, 2, function(x) p.adjust(x, method = "BH"))
  return(list(gene_ctSignif, ct_coef, gene_ctSignif_adj))
}

###test
#x <- colnames(test_list_cell_ctDist[[1]])[1]
#test_abund_data <- abundence_meter(x, test_list_cell_ctDist)
#test_gene_ctCor_data <- htest_data_extractor(test_abund_data, test_list_tt_expData)

##### Get signal_Stat object with complete signal information from correlation analysis
# Time-concuming (~30min per cell type)
#' Get signal_Stat object with complete signal information from correlation analysis
#'
#' @param list_cell_ctDist
#' @param list_tt_expData
#' @param bg_ctCor_data
#' @param bg_anova
#' @param method
#' @param bg_anova_alpha
#' @param bg_ctCor_alpha
#' @param gene_ctCor_alpha
#'
#' @return
#'
#' @examples
signal_extractor <- function(list_cell_ctDist, list_tt_expData, bg_ctCor_data, bg_anova, method = "spearman",
                             bg_anova_alpha = 0.01, bg_ctCor_alpha = 0.01, gene_ctCor_alpha = 0.01){
  signal_Stat <- lapply(colnames(list_cell_ctDist[[1]]), function(x){
    print(paste0("Detecting signals from ", x, " ..." ))
    abund_data <- abundence_meter(x, list_cell_ctDist)
    gene_ctCor_data <- htest_data_extractor(abund_data, list_tt_expData, method = method) ### Time-concuming (~15min)
    gene_ctDistr <- cell_filter(gene_ctCor_data, bg_ctCor_data, bg_anova,
                                bg_anova_alpha = bg_anova_alpha, bg_ctCor_alpha = bg_ctCor_alpha, gene_ctCor_alpha = gene_ctCor_alpha)
    list_gene_ttStat <- target_type_data_extractor(gene_ctDistr, bg_ctCor_data)
    return(list_gene_ttStat)
  })
  names(signal_Stat) <- colnames(list_cell_ctDist[[1]])
  return(signal_Stat)
}

###test
#test_signal_Stat <- signal_extractor(test_list_cell_ctDist, test_list_tt_expData, test_bg_ctCor_data, test_bg_anova)

##### The first function for filtration of signals
#' The first function for filtration of signals
#'
#' @param signal_Stat
#' @param pValue
#' @param threshold
#'
#' @return
#'
#' @examples
signal_filter <- function(signal_Stat, pValue = 0.01, threshold = 0.1){
  signals <- lapply(signal_Stat, function(x){
    tt_signal <- lapply(x, function(y){
      tmp <- y[y$filter_pass & (y$cor_significance <= pValue) & (abs(y$cor_coefficient) >= threshold),]
      tmp <- tmp[!is.na(tmp$filter_pass),]
      return(tmp)
    })
    names(tt_signal) <- names(x)
    return(tt_signal)
  })
  names(signals) <- names(signal_Stat)
  return(signals)
}

###test
#test_signals <- signal_filter(test_signal_Stat)


