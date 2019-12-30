#library(flowCore)
#library(matrixStats)
#library(Rtsne)


##### Create metadata
#' Create basic metadata for flowSet object
#' @description The function creates two-column dataframe with paths
#' and names of files by the list of paths to fcs files with CyTOF data.
#' this data frame use as base metadata or phenoData in flowSet object
#' from flowCore package.
#'
#' @param fcs_files list with paths
#'
#' @return The dataframe with paths to fcs files within one column and names of these files in another column
#' @export
#' @import flowCore
#'
#' @examples
get_fcs_metadata <- function(fcs_files){
  md <- data.frame(path = fcs_files)
  md$path <- as.character(md$path)
  md$file_name <- basename(md$path)
  return(md)
}

##### creat FlowSet object
#' Creat FlowSet object
#' @description The function takes dataframe with column contained  paths to
#' fcs files and creates flowSet object from flowCore package.
#'
#' @param md The input metadata should be dataframe format with
#' a column named path. The function assumed to process the result
#' of the function \dontrun{
#' get_fcs_metadata
#' }
#'
#' @return flowSet object from flowCore package
#' @export
#' @import flowCore
#'
#' @examples
get_fcs_raw <- function(md){
  pathes <- as.vector(md$path)
  fcs_raw <- read.flowSet(pathes, transformation = FALSE, truncate_max_range = FALSE)
  return(fcs_raw)
}

##### Create panel data
#' Create panel data to fowSet object
#' @description extract panel data about markers from flowSet
#' object from flowCore package. The function tries to assign
#' a description to markers and its type (technical or not)
#'
#' @param fcs_raw flowSet object from flowCore package
#'
#' @return data frame with information about the markers
#' @export
#' @import flowCore
#'
#' @examples
get_fcs_panel <- function(fcs_raw){
  tech_patterns <-list(computational_tech = c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE", "BCKG"),
                       marker_tech = c("_BC", "BCKG", "DNA", "Cisplatin"))
  panel <- as.data.frame(pData(parameters(fcs_raw[[1]]))[,c("name", "desc")])
  rownames(panel) <- NULL
  panel <- panel[sapply(panel$name, function(x) !any(sapply(tech_patterns$computational_tech, function(y) grepl(y,x)))),]
  panel$antigen <- sapply(strsplit(panel$desc, "_"), function(x) x[length(x)])
  panel$marker_class <- "type"
  panel$marker_class[sapply(panel$desc, function(x) any(sapply(tech_patterns$marker_tech, function(y) grepl(y,x))))] <- "state"
  panel <- as.data.frame(apply(panel, c(1,2), function(x) gsub(" ", "_", x)))
  return(panel)
}

##### Create use_marker
#' filtering out technical markers
#' @description the function formed list of markers which have class
#' mentioned as "type" in the panel. Other markers with class "state"
#' are considered as technical signals
#'
#' @param panel dataframe with marker panel data
#'
#' @return the list value of which are marker names and names of the
#' list element are names of antibody, relevant to the markers
#'
#'
#' @examples
get_use_marker <- function(panel){
  use_markers <- as.character(panel$name[panel$marker_class == "type"])
  names(use_markers) <- panel$antigen[panel$marker_class == "type"]
  names(use_markers) <- gsub("_(v)", "", names(use_markers), fixed = T)          ### It should be fixed (problem with "_(v)")
  return(use_markers)
}

##### Upload data from fcs files
upload_fcs_data <- function(fcs_files){
  md <- get_fcs_metadata(fcs_files)
  fcs_raw <- get_fcs_raw(md)
  panel <- get_fcs_panel(fcs_raw)
  use_markers <- get_use_marker(panel)
  return(list(fcs_raw,md, panel, use_markers))
}

#### Transformation from "count" to "asinh" data
#' Transformation from "count" to "asinh" data
#'
#' @param fcs_raw flowSet object
#' @param cofactor digit with used as denominator to transform data
#'
#' @return flowSet object with transform data by asinh function and divided
#' by cofactor
#' @import flowCore
#'
#' @examples
asinh_transformation <- function(fcs_raw, cofactor){
  computational_tech <- c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE")
  markers <- colnames(fcs_raw[[1]])[sapply(colnames(fcs_raw[[1]]), function(x) !any(sapply(computational_tech, function(y) grepl(y,x))))]
  fcs_asinh <- fsApply(fcs_raw, function(x, cf = cofactor, mk = markers){
    exprs(x)[,mk] <- asinh(exprs(x)[,mk] / cf)
    return(x)})
  return(fcs_asinh)
}

##### Transformation to a from 0 to 1 variable and removing outliers
#' Transformation to a from 0 to 1 variable and removing outliers
#'
#' @param fcs_raw
#' @param quantile
#'
#' @return
#' @export
#' @import flowCore matrixStats
#' @examples
outlier_by_quantile_transformation <- function(fcs_raw, quantile){
  computational_tech <- c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE")
  markers <- colnames(fcs_raw[[1]])[sapply(colnames(fcs_raw[[1]]), function(x) !any(sapply(computational_tech, function(y) grepl(y,x))))]
  fcs_outlier_by_quantile <- fsApply(fcs_raw, function(x, ql = quantile, mk = markers){
    rng <- colQuantiles(exprs(x)[,mk], probs = c(ql, 1-ql))
    expr_data <- t((t(exprs(x)[,mk]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    expr_data[expr_data < 0] <- 0
    expr_data[expr_data > 1] <- 1
    exprs(x)[,mk] <- expr_data
    return(x)
  })
  return(fcs_outlier_by_quantile)
}

#################
### tSNE plot ###
#################

##### Preparing data to tSNE
#' Preparing data to tSNE
#'
#' @param fcs_raw
#' @param use_markers
#' @param sampling_size
#'
#' @return
#' @export
#' @import Rtsne flowCore
#'
#' @examples
sampled_tSNE <- function(fcs_raw, use_markers, sampling_size = 1000){
  expr <- fsApply(fcs_raw[,use_markers], exprs)
  sample_ids <- rep(sampleNames(fcs_raw), fsApply(fcs_raw, nrow))

  ## Find and skip duplicates
  dups <- which(!duplicated(expr[, use_markers]))

  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids)

  ## How many cells to downsample per-sample
  tsne_ncells <- pmin(table(sample_ids), sampling_size)               ################ Number of cells to ploting

  ## Get subsampled indices
  set.seed(1234)
  tsne_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  tsne_inds <- unlist(tsne_inds)
  tsne_expr <- expr[tsne_inds, use_markers]

  ##### Run t-SNE
  set.seed(1234)
  tsne_result <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)
  #tsne_out <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2])
  tsne_out <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2], expr[tsne_inds, use_markers])
  colnames(tsne_out)[match(use_markers, colnames(tsne_out))] <- names(use_markers)
  return(tsne_out)
}

