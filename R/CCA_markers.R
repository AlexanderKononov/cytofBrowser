####################################################################
#### Functions of CCA Core (based complex correlation analysis) ####
####################################################################

##### Functions to create background contrast objects used for correlation analysis
###################################################################################

#### Calculate background correlation landscape between markers across the entire cell mass by naive without separation on samples
#' Background correlation landscape between markers
#'
#' @description Calculate background correlation landscape between markers across the entire cell mass by naive without separation on samples
#' @param list_expData
#' @param use_markers
#' @param method
#'
#' @return
#' @importFrom stats p.adjust
#'
#' @examples
get_bg_naive_corr_land_mk <- function(list_expData, use_markers = NULL, method = "spearman"){
  ifelse(is.null(use_markers), mk_names <- rownames(list_expData[[1]]), mk_names <- use_markers)
  bg_htest <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
  rownames(bg_htest) <- mk_names
  colnames(bg_htest) <- mk_names
  bg_coef <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
  rownames(bg_coef) <- mk_names
  colnames(bg_coef) <- mk_names

  comparison_unit <- as.data.frame(t(do.call(cbind, list_expData)))[,mk_names]

  for(i in mk_names){
    htest <- sapply(comparison_unit[,mk_names != i], cor.test, y = comparison_unit[,i], method = method)
    bg_htest[i,mk_names != i] <- unlist(htest["p.value",])
    bg_coef[i,mk_names != i] <- unlist(htest["estimate",])
  }
  bg_padj <- apply(bg_htest, 2, function(x) stats::p.adjust(x, method = "BH"))
  bg_naive_corr_mk <- list(bg_htest, bg_coef, bg_padj)
  names(bg_naive_corr_mk) <- c("bg_htest", "bg_coef", "bg_padj")
  return(bg_naive_corr_mk)
}

### test
#test_bg_naive_corr_mk <- get_bg_naive_corr_land_mk(test_list_expData)

#### Calculate background correlation landscape between markers across samples
#' Calculate background correlation landscape between markers across samples
#'
#' @param ist_expData
#' @param use_markers
#' @param method
#'
#' @return
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarise_all
#' @importFrom stats cor.test p.adjust
#'
#' @examples
get_bg_corr_land_mk <- function(list_expData, use_markers = NULL, method = "spearman"){
  ifelse(is.null(use_markers), mk_names <- rownames(list_expData[[1]]), mk_names <- use_markers)
  bg_htest <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
  rownames(bg_htest) <- mk_names
  colnames(bg_htest) <- mk_names
  bg_coef <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
  rownames(bg_coef) <- mk_names
  colnames(bg_coef) <- mk_names

  comparison_unit <- as.data.frame(t(do.call(cbind, list_expData)))[,mk_names]
  comparison_unit$samples <- rep(names(list_expData), lapply(list_expData, ncol))
  comparison_unit <- dplyr::group_by(comparison_unit, samples) %>% dplyr::summarise_all(median) %>% as.data.frame()
  comparison_unit$samples <-NULL
  for(i in mk_names){
    htest <- sapply(comparison_unit[,mk_names != i], stats::cor.test, y = comparison_unit[,i], method = method)
    bg_htest[i,mk_names != i] <- unlist(htest["p.value",])
    bg_coef[i,mk_names != i] <- unlist(htest["estimate",])
  }
  bg_padj <- apply(bg_htest, 2, function(x) stats::p.adjust(x, method = "BH"))
  bg_corr_mk <- list(bg_htest, bg_coef, bg_padj)
  names(bg_corr_mk) <- c("bg_htest", "bg_coef", "bg_padj")
  return(bg_corr_mk)
}

### test
#bg_corr_mk <- get_bg_corr_land_mk(ist_expData)

#### Get bg_anova_land_mk object with background ANOVA statistics for markers across samples
#' Get bg_anova_land_mk object with background ANOVA statistics for markers across samples
#'
#' @param list_expData
#' @param use_markers
#'
#' @return
#' @importFrom stats aov p.adjust
#'
#' @examples
get_bg_anova_land_mk <- function(list_expData, use_markers = NULL){
  ifelse(is.null(use_markers), mk_names <- rownames(list_expData[[1]]), mk_names <- use_markers)
  comparison_unit <- as.data.frame(t(do.call(cbind, list_expData)))[,mk_names]
  comparison_unit$samples <- rep(names(list_expData), lapply(list_expData, ncol))
  bg_anova <- sapply(mk_names, function(i){
    tmp <- as.data.frame(comparison_unit[,c('samples', i)])
    colnames(tmp) <- c('samples', 'expression')
    aov_data <- stats::aov(formula = expression ~ samples, data = tmp)
    return(summary(aov_data)[[1]][["Pr(>F)"]][1])
  })
  names(bg_anova) <- mk_names
  bg_anova_adj <-  stats::p.adjust(bg_anova, method = "BH")
  bg_anova_land_mk <- list(bg_anova, bg_anova_adj)
  names(bg_anova_land_mk) <- c("bg_anova", "bg_anova_adj")
  return(bg_anova_land_mk)
}

### test
#bg_anova_land_mk <- get_bg_anova_land_mk(list_expData)

####  Get anova_mk object with ANOVA statistics for markers within each cluster across samples
#' Get anova_mk object with ANOVA statistics for markers within each cluster across samples
#'
#' @param list_cell_ctDist
#' @param list_expData
#' @param use_markers
#'
#' @return
#' @importFrom stats aov p.adjust
#'
#' @examples
get_anova_mk <- function(list_cell_ctDist, list_expData, use_markers = NULL){
  ifelse(is.null(use_markers), mk_names <- rownames(list_expData[[1]]), mk_names <- use_markers)

  anova_mk <- lapply(colnames(list_cell_ctDist[[1]]), function(cl) {
    comparison_unit <- as.data.frame(t(do.call(cbind, lapply(1:length(list_expData), function(s)
      list_expData[[s]][,list_cell_ctDist[[s]][,cl]]) )))[,mk_names]
    comparison_unit$samples <- rep(names(list_expData), sapply(list_cell_ctDist, function(x) sum(x[,cl])))
    bg_anova <- sapply(mk_names, function(i){
      tmp <- as.data.frame(comparison_unit[,c('samples', i)])
      colnames(tmp) <- c('samples', 'expression')
      aov_data <- stats::aov(formula = expression ~ samples, data = tmp)
      return(summary(aov_data)[[1]][["Pr(>F)"]][1])
    })
    names(bg_anova) <- mk_names
    bg_anova_adj <-  stats::p.adjust(bg_anova, method = "BH")
    bg_anova_mk <- list(bg_anova, bg_anova_adj)
    names(bg_anova_mk) <- c("bg_anova", "bg_anova_adj")
    return(bg_anova_mk)
  })
  names(anova_mk) <- colnames(list_cell_ctDist[[1]])
  return(anova_mk)
}

### test
#anova_mk <- get_anova_mk(list_cell_ctDist, list_expData)

#### Get contrast correlation landscape for each clone each marker versus each marker of the entire cell mass
#' Contrast correlation landscape
#' @description Get contrast correlation landscape for each clone each marker versus each marker of the entire cell mass
#'
#' @param list_cell_ctDist
#' @param list_expData
#' @param use_markers
#' @param method
#'
#' @return
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarise_all
#' @importFrom stats cor.test p.adjust
#'
#' @examples
get_contrast_corr_land_mk <- function(list_cell_ctDist, list_expData, use_markers = NULL, method = "spearman"){
  ifelse(is.null(use_markers), mk_names <- rownames(list_expData[[1]]), mk_names <- use_markers)
  entire_expression <- as.data.frame(t(do.call(cbind, list_expData)))[,mk_names]
  entire_expression$samples <- rep(names(list_expData), lapply(list_expData, ncol))
  entire_expression <- dplyr::group_by(entire_expression, samples) %>% dplyr::summarise_all(median) %>% as.data.frame()
  rownames(entire_expression) <- entire_expression$samples
  entire_expression$samples <- NULL
  tt_expression <- lapply(colnames(list_cell_ctDist[[1]]), function(cl){
    tmp <- do.call(rbind, lapply(1:length(list_expData), function(s){
      sapply(mk_names, function(m){median(list_expData[[s]][m,list_cell_ctDist[[s]][,cl]])})
    }))
    colnames(tmp) <- mk_names
    rownames(tmp) <- names(list_expData)
    return(tmp)
  })
  names(tt_expression) <- colnames(list_cell_ctDist[[1]])

  contrast_corr_land_mk <- lapply(names(tt_expression), function(cl){
    bg_htest <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
    rownames(bg_htest) <- mk_names
    colnames(bg_htest) <- mk_names
    bg_coef <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
    rownames(bg_coef) <- mk_names
    colnames(bg_coef) <- mk_names
    for(i in mk_names){
      htest <- sapply(entire_expression, stats::cor.test, y = tt_expression[[cl]][,i], method = method)
      bg_htest[i,] <- unlist(htest["p.value",])
      bg_coef[i,] <- unlist(htest["estimate",])
    }
    bg_padj <- apply(bg_htest, 2, function(x) stats::p.adjust(x, method = "BH"))
    bg_ctCor_data <- list(bg_htest, bg_coef, bg_padj)
    names(bg_ctCor_data) <- c("bg_htest", "bg_coef", "bg_padj")
    return(bg_ctCor_data)
  })
  names(contrast_corr_land_mk) <- names(tt_expression)
  return(contrast_corr_land_mk)
}

### test
#contrast_corr_land_mk <- get_contrast_corr_land_mk(list_cell_ctDist, list_expData)

#### Get contrast naive correlation landscape for each clone each marker versus each marker of the entire cell mass
#' Get contrast naive correlation landscape for each clone each marker versus each marker of the entire cell mass
#'
#' @param list_cell_ctDist
#' @param list_expData
#' @param use_markers
#' @param method
#'
#' @return
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarise_all
#' @importFrom stats cor.test p.adjust
#'
#' @examples
get_contrast_naive_corr_land_mk <- function(list_cell_ctDist, list_expData, use_markers = NULL, method = "spearman"){
  ifelse(is.null(use_markers), mk_names <- rownames(list_expData[[1]]), mk_names <- use_markers)
  entire_expression <- as.data.frame(t(do.call(cbind, list_expData)))[,mk_names]
  entire_expression$samples <- rep(names(list_expData), lapply(list_expData, ncol))
  entire_expression <- dplyr::group_by(entire_expression, samples) %>% dplyr::summarise_all(median) %>% as.data.frame()
  rownames(entire_expression) <- entire_expression$samples
  entire_expression$samples <- NULL

  contrast_naive_corr_land_mk <- lapply(colnames(list_cell_ctDist[[1]]), function(cl){
    bg_htest <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
    rownames(bg_htest) <- mk_names
    colnames(bg_htest) <- mk_names
    bg_coef <- matrix(data = NA, nrow = length(mk_names), ncol = length(mk_names))
    rownames(bg_coef) <- mk_names
    colnames(bg_coef) <- mk_names
    for(m in mk_names){
      comparison_unit <- do.call(rbind, lapply(1:length(list_expData), function(s)
        data.frame(expression = list_expData[[s]][m,list_cell_ctDist[[s]][,cl]])))
      comparison_unit <- cbind(comparison_unit, do.call(cbind, lapply(mk_names, function(m2)
        rep(entire_expression[,m2], time = sapply(list_cell_ctDist, function(x) sum(x[,cl]))))))
      colnames(comparison_unit) <- c('tt_expression', mk_names)
      htest <- sapply(comparison_unit[,-1], stats::cor.test, y = comparison_unit$tt_expression, method = method)
      bg_htest[m,] <- unlist(htest["p.value",])
      bg_coef[m,] <- unlist(htest["estimate",])
    }
    bg_padj <- apply(bg_htest, 2, function(x) stats::p.adjust(x, method = "BH"))
    bg_ctCor_data <- list(bg_htest, bg_coef, bg_padj)
    names(bg_ctCor_data) <- c("bg_htest", "bg_coef", "bg_padj")
    return(bg_ctCor_data)
  })
  names(contrast_naive_corr_land_mk) <- colnames(list_cell_ctDist[[1]])
  return(contrast_naive_corr_land_mk)
}

### test
#contrast_naive_corr_land_mk <- get_contrast_naive_corr_land_mk(list_cell_ctDist, list_expData)


##### Functions to extract correlation signals ( not different from the background yet )
############################################################################

##### Low-level function for correlation analysis between markers of different clusters
#' Correlation analysis between markers of different clusters
#'
#' @param list_tt_expData
#' @param use_markers
#' @param method
#'
#' @return
#' @importFrom stats cor.test
#'
#' @examples
htest_data_extractor_mk <- function(list_tt_expData, use_markers = NULL, method = "spearman"){
  ifelse(is.null(use_markers), mk_names <- rownames(list_tt_expData[[1]][[1]]), mk_names <- use_markers)
  gene_to_gene_cor <- do.call(rbind, lapply(names(list_tt_expData[[1]]), function(i){
    do.call(rbind, lapply(names(list_tt_expData[[1]]), function(j){
      do.call(rbind, lapply(mk_names, function(g1){
        do.call(rbind, lapply(mk_names, function(g2){
          comparison_unit <- lapply(1:length(list_tt_expData), function(x) {
            if(length(list_tt_expData[[x]][[i]])==0){return(NULL)}
            if(length(list_tt_expData[[x]][[j]])==0){return(NULL)}
            tmp <- data.frame(signaling_gene = median(list_tt_expData[[x]][[i]][g1,]),
                              target_gene = median(list_tt_expData[[x]][[j]][g2,]))
            return(tmp)
          })
          comparison_unit <- do.call(rbind, comparison_unit)
          if(length(comparison_unit)<=1){return(NULL)}
          if(length(unique(comparison_unit$signaling_gene))<=1 & length(unique(comparison_unit$target_gene))<=1){return(NULL)}
          res <- stats::cor.test(comparison_unit$signaling_gene, comparison_unit$target_gene, method = method)
          return(data.frame(signaling_cluster = i, signaling_marker = g1, target_cluster = j, target_marker = g2,
                            cor_significance = res$p.value, cor_coefficient = res$estimate,
                            number_of_samples = length(comparison_unit$signaling_gene)))
        }))
      }))
    }))
  }))
  return(gene_to_gene_cor)
}

### test
#test_gene_to_gene_cor <- htest_data_extractor_mk(test_list_tt_expData)

#### Compare background correlation landscape versus target correlations
#' Compare background correlation landscape versus target correlations
#' and extract information whether signals pass-through or get the complete statistic information
#'
#' @param gene_to_gene_cor
#' @param anova_mk
#' @param bg_corr_mk
#' @param contrast_corr_land_mk
#' @param anova_mk_alpha
#' @param bg_corr_mk_alpha
#' @param contrast_corr_land_mk_alpha
#' @param statistic_output if TRUE, function return complete statistic information. Defaul default is FALSE, return pass or not inf
#'
#' @return
#' @importFrom magrittr "%>%"
#'
#' @examples
comparator_mk <- function(gene_to_gene_cor, anova_mk, bg_corr_mk, contrast_corr_land_mk,
                          anova_mk_alpha = 0.01, bg_corr_mk_alpha = 0.01, contrast_corr_land_mk_alpha = 0.01,
                          use_markers = NULL, statistic_output = F){
  ifelse(is.null(use_markers), mk_names <- names(anova_mk[[1]][['bg_anova']]), mk_names <- use_markers)
  sig_stat <- data.frame(anova_pVal_signal_mk = NA, anova_padj_signal_mk = NA, anova_pVal_targ_mk = NA, anova_padj_targ_mk = NA,
             bg_pVal = NA, bg_coef = NA, bg_padj = NA,
             contr_pVal_signal_mk = NA, contr_coef_signal_mk = NA, contr_padj_signal_mk = NA,
             contr_pVal_targ_mk = NA, contr_coef_targ_mk = NA, contr_padj_targ_mk = NA)
  for (i in 1:nrow(gene_to_gene_cor)){
    x <- gene_to_gene_cor[1,]
    if (!(x[['signaling_marker']] %in% mk_names)){next}
    if (!(x[['target_marker']] %in% mk_names)){next}
    tmp <- c(anova_mk[[x[['signaling_cluster']]]][['bg_anova']][x[['signaling_marker']]],
      anova_mk[[x[['signaling_cluster']]]][['bg_anova_adj']][x[['signaling_marker']]],
      anova_mk[[x[['target_cluster']]]][['bg_anova']][x[['target_marker']]],
      anova_mk[[x[['target_cluster']]]][['bg_anova_adj']][x[['target_marker']]],
      bg_corr_mk[['bg_htest']][x[['signaling_marker']], x[['target_marker']]],
      bg_corr_mk[['bg_coef']][x[['signaling_marker']], x[['target_marker']]],
      bg_corr_mk[['bg_padj']][x[['signaling_marker']], x[['target_marker']]],
      contrast_corr_land_mk[[x[['signaling_cluster']]]][['bg_htest']][x[['signaling_marker']], x[['target_marker']]],
      contrast_corr_land_mk[[x[['signaling_cluster']]]][['bg_coef']][x[['signaling_marker']], x[['target_marker']]],
      contrast_corr_land_mk[[x[['signaling_cluster']]]][['bg_padj']][x[['signaling_marker']], x[['target_marker']]],
      contrast_corr_land_mk[[x[['target_cluster']]]][['bg_htest']][x[['target_marker']], x[['signaling_marker']]],
      contrast_corr_land_mk[[x[['target_cluster']]]][['bg_coef']][x[['target_marker']], x[['signaling_marker']]],
      contrast_corr_land_mk[[x[['target_cluster']]]][['bg_padj']][x[['target_marker']], x[['signaling_marker']]]
    )
    names(tmp) <- c('anova_pVal_signal_mk', 'anova_padj_signal_mk', 'anova_pVal_targ_mk', 'anova_padj_targ_mk',
                    'bg_pVal', 'bg_coef', 'bg_padj',
                    'contr_pVal_signal_mk', 'contr_coef_signal_mk', 'contr_padj_signal_mk',
                    'contr_pVal_targ_mk', 'contr_coef_targ_mk', 'contr_padj_targ_mk')
    sig_stat <- rbind(sig_stat,tmp)
  }
  sig_stat <- sig_stat[-1,]

  if(statistic_output){return(sig_stat)}
  sig_summary <- data.frame(anova_pVal_signal_mk = (sig_stat$anova_pVal_signal_mk <= anova_mk_alpha),
                            anova_padj_signal_mk = (sig_stat$anova_padj_signal_mk <= anova_mk_alpha),
                            anova_pVal_targ_mk = (sig_stat$anova_pVal_targ_mk <= anova_mk_alpha),
                            anova_padj_targ_mk = (sig_stat$anova_padj_targ_mk <= anova_mk_alpha),
                            bg_coef = sig_stat$bg_coef,
                            bg_pVal = (sig_stat$bg_pVal > bg_corr_mk_alpha),
                            bg_padj = (sig_stat$bg_padj > bg_corr_mk_alpha),
                            contr_coef_signal_mk = sig_stat$contr_coef_signal_mk,
                            contr_pVal_signal_mk = (sig_stat$contr_pVal_signal_mk > contrast_corr_land_mk_alpha),
                            contr_padj_signal_mk = (sig_stat$contr_padj_signal_mk > contrast_corr_land_mk_alpha),
                            contr_coef_targ_mk = sig_stat$contr_coef_targ_mk,
                            contr_pVal_targ_mk = (sig_stat$contr_pVal_targ_mk > contrast_corr_land_mk_alpha),
                            contr_padj_targ_mk = (sig_stat$contr_padj_targ_mk > contrast_corr_land_mk_alpha))
  return(sig_summary)
}

### test
#sig_summary <- comparator_mk(gene_to_gene_cor, anova_mk, bg_corr_mk, contrast_corr_land_mk)

#### Filtering signals by background comparison information
#' Filtering signals by background comparison information
#'
#' @param gene_to_gene_cor
#' @param sig_summary
#' @param p_adjusted use adjusted p-value or not. Default TRUE, use adjusted p-value.
#' @param na_bg default TRUE - change NA in the background and contrast correlation to TRUE or "pass" value.
#'
#' @return
#'
#' @examples
filter_mk <- function(gene_to_gene_cor, sig_summary, threshold = 0.1, pValue = 0.01, p_adjusted = TRUE, na_bg = TRUE){
  if(p_adjusted){pass_list <- sig_summary[,c('anova_padj_signal_mk', 'anova_padj_targ_mk', 'bg_padj',
                                             'contr_padj_signal_mk', 'contr_padj_targ_mk')]}
  if(!p_adjusted){pass_list <- sig_summary[,c('anova_pVal_signal_mk', 'anova_pVal_targ_mk', 'bg_pVal',
                                              'contr_pVal_signal_mk', 'contr_pVal_targ_mk')]}
  colnames(pass_list) <- c('anova_signal_mk', 'anova_targ_mk', 'bg_corr', 'contr_signal_mk', 'contr_targ_mk')
  pass_list$anova_signal_mk[is.na(pass_list$anova_signal_mk)] <- na_bg
  pass_list$anova_targ_mk[is.na(pass_list$anova_targ_mk)] <- na_bg
  pass_list$bg_corr[is.na(pass_list$bg_corr)] <- na_bg
  pass_list$contr_signal_mk[is.na(pass_list$contr_signal_mk)] <- na_bg
  pass_list$contr_targ_mk[is.na(pass_list$contr_targ_mk)] <- na_bg
  pass_list$pass_decision <- pass_list$anova_signal_mk & pass_list$anova_targ_mk & pass_list$bg_corr
  pass_list$pass_type <- 'not_contrast'
  pass_list$pass_type[pass_list$contr_signal_mk] <- 'signal_mk_contrast'
  pass_list$pass_type[pass_list$contr_targ_mk] <- 'targ_mk_contrast'
  pass_list$pass_type[pass_list$contr_signal_mk & pass_list$contr_targ_mk] <- 'full_contrast'

  signals <- data.frame(gene_to_gene_cor, contrast = pass_list$pass_type)
  signals <- signals[pass_list$pass_decision,]
  rownames(signals) <- NULL
  signals <- signals[signals$cor_significance <= pValue,]
  signals <- signals[order(abs(signals$cor_coefficient)),]
  signals <- signals[abs(signals$cor_coefficient) > threshold,]
  signals <- signals[apply(signals, 1, function(x) !all(is.na(x))),]
  return(signals)
}

### test
#signals <- filter_mk(gene_to_gene_cor, sig_summary)

#### Extract signals which go from and to the same cluster
#' Extract signals which go from and to the same cluster
#'
#' @param signals
#'
#' @return
#'
#' @examples
get_signal_in_cluster_mk <- function(signals){
  if (nrow(signals)==0){return(NULL)}
  if (length(signals$signaling_cluster == signals$target_cluster)==0){return(NULL)}
  signal_in_cluster_mk <- signals[signals$signaling_cluster == signals$target_cluster,]
  signal_in_cluster_mk$target_cluster <- NULL
  colnames(signal_in_cluster_mk) <- c("cluster", "marker_1",  "marker_2", "cor_significance",
                                      "cor_coefficient", "number_of_samples", 'contrast')
  signal_in_cluster_mk$names <- paste0(signal_in_cluster_mk$cluster, "-", signal_in_cluster_mk$marker_1, "_",
                                       signal_in_cluster_mk$marker_2)
  return(signal_in_cluster_mk)
}

### test
#signal_in_cluster_mk <- get_signal_in_cluster_mk(signals)

#### Extract signals which go from one cluster and go to another
#' Extract signals which go from one cluster and go to another
#'
#' @param signals
#'
#' @return
#'
#' @examples
get_signal_between_cluster_mk <- function(signals){
  if (nrow(signals)==0){return(NULL)}
  if (length(signals$signaling_cluster != signals$target_cluster)==0){return(NULL)}
  signal_between_cluster_mk <- signals[signals$signaling_cluster != signals$target_cluster,]
  signal_between_cluster_mk$names <- paste0(signal_between_cluster_mk$signaling_cluster, "-",
                                            signal_between_cluster_mk$signaling_marker, "_",
                                            signal_between_cluster_mk$target_cluster, "_",
                                            signal_between_cluster_mk$target_marker)
  return(signal_between_cluster_mk)
}

### test
#signal_between_cluster_mk <- get_signal_between_cluster_mk(signals)
