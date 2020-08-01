
#' Take data to gating scatterplot
#'
#' @param fcs_raw
#' @param gating_subset
#' @param gating_mk1
#' @param gating_mk2
#'
#' @return
#' @export
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"
#'
#' @examples
get_data_for_gating <- function(fcs_raw, gating_subset, gating_mk1, gating_mk2){
  raw_df <- flowCore::fsApply(fcs_raw, function(x) flowCore::exprs(x))
  raw_df <- as.data.frame(raw_df)
  raw_df$original_cell_coordinates <- 1:nrow(raw_df)
  raw_df <- raw_df[gating_subset, c(gating_mk1, gating_mk2, "original_cell_coordinates")]
  colnames(raw_df) <- c("gating_mk1", "gating_mk2", "original_cell_coordinates")
  return(raw_df)
}

### test
#raw_df <- get_data_for_gating(fcs_raw, rep(TRUE, 13402), "I127Di", "Tm169Di")

#' Function adds set of gates to cell annotation object
#'
#' @param gates
#' @param gate_list
#' @param cell_annotation
#' @param method
#'
#' @return
#' @export
#'
#' @examples
get_cell_type_from_gates <- function(gates, gate_list, cell_annotation, method = 'pure'){
  target_gates <- as.matrix(gates[,as.character(unlist(gate_list))])
  populations <- apply(target_gates, 1, function(cell){paste(gate_list[cell], sep = "", collapse = "_&_")})
  populations[populations == ""] <- "untyped"
  new_cell_annotation <- cell_annotation
  if(method == 'pure'){
    new_name <- "pure_gated_cell_type"
  }
  if(method == 'squeeze'){
    purity <- rowSums(target_gates)
    populations[!(purity == 1)] <- "untyped"
    new_name <- "squeeze_gated_cell_type"
  }
  new_cell_annotation$new_gate <- populations
  if(any(grepl(new_name, colnames(new_cell_annotation)))){
    new_name <- paste0(new_name,"_",sum(grepl(new_name, colnames(new_cell_annotation)))+1)}
  colnames(new_cell_annotation)[colnames(new_cell_annotation) == "new_gate"] <- new_name
  return(new_cell_annotation)
}

#' Extract annotation fron fcs_raw object to cell annotation object
#'
#' @param entire_panel
#' @param fcs_raw
#' @param cell_annotation
#'
#' @return
#' @export
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"
#'
#' @examples
get_add_cell_annotation_from_data <- function(entire_panel, fcs_raw, cell_annotation){
  raw_df <- flowCore::fsApply(fcs_raw, function(x) flowCore::exprs(x)[,entire_panel])
  raw_df <- as.data.frame(raw_df)
  cell_annotation <- cbind(cell_annotation, raw_df)
  return(cell_annotation)
}


#### Adding of cell annatation to flowSet object
#' Adding of cluster-info to flowSet object
#'
#' @param fcs_raw
#' @param cell_clustering_list
#'
#' @return
#' @export
#' @importFrom flowCore sampleNames fr_append_cols
#' @importClassesFrom flowCore flowSet
#'
#' @examples
get_cell_annotation_fcs_files <- function(fcs_raw, cell_annotation, column_names = NULL){
  if(is.null(column_names)){column_names <- colnames(cell_annotation)[-1]}
  if(!all(flowCore::sampleNames(fcs_raw) %in% unique(cell_annotation$samples))){
    print("Cluster and data samples does not match")
    return(NULL)}
  clustered_fcs <- lapply(flowCore::sampleNames(fcs_raw), function(s) {
    addition_col <- as.matrix(cell_annotation[cell_annotation$samples == s, column_names])
    addition_col <- apply(addition_col, 2, function(x) as.integer(as.factor(x)))
    colnames(addition_col) <- column_names
    fs_d <- flowCore::fr_append_cols(fcs_raw[[s]], addition_col)
    return(fs_d)
  })
  names(clustered_fcs) <- flowCore::sampleNames(fcs_raw)
  clustered_fcs <- as(clustered_fcs, 'flowSet')
  return(clustered_fcs)
}

