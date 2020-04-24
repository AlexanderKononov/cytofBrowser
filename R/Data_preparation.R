
##### Create metadata
#' Create basic metadata for flowSet object
#' @description The function creates two-column dataframe with paths
#' and names of files by the list of paths to fcs files with CyTOF data.
#' this data frame use as base metadata or phenoData in flowSet object
#' from flowCore package.
#'
#' @param fcs_files list with paths
#' @export
#' @return The dataframe with paths to fcs files within one column and names of these files in another column
#'
#' @examples
get_fcs_metadata <- function(fcs_files){
  md <- data.frame(path = fcs_files)
  md$path <- as.character(md$path)
  md$file_name <- basename(md$path)

  samples <- gsub(".fcs", "", md$file_name)
  md$short_name <- unlist(lapply(1:length(samples), function(x){
    tmp <- unlist(strsplit(samples[x], ""))
    if(length(tmp) <= 12){return(samples[x])}
    return(paste0(paste0(tmp[1:7], collapse = ""), "_smpl_", as.character(x)))
  }))

  return(md)
}

### test
#md <- get_fcs_metadata(c("./test_data/c13_20190704_hnp_perf_11_0_Alex2.fcs","./test_data/c14_20190704_hnp_perf_11_0_Alex2.fcs",
#                        "./test_data/c13_20190704_hnp_perf_11_0_Alez1.fcs","./test_data/c14_20190704_hnp_perf_11_0_Alez1.fcs"))
#md <- get_fcs_metadata("./../../Toni_data/EC_200117_Freshly labelled PBMCs _1_0/Activation_Activation full panel unstim TILs_033.fcs")

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
#' @importFrom flowCore read.flowSet sampleNames "sampleNames<-"
#'
#' @examples
get_fcs_raw <- function(md){
  pathes <- as.vector(md$path)
  fcs_raw <- flowCore::read.flowSet(pathes, transformation = FALSE, truncate_max_range = FALSE)
  sampleNames(fcs_raw) <- unlist(gsub(".fcs", "", flowCore::sampleNames(fcs_raw)))
  return(fcs_raw)
}

### test
#fcs_raw <- get_fcs_raw(md)

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
#' @importFrom flowCore pData "pData<-" parameters "parameters<-"
#' @importClassesFrom flowCore flowSet
#'
#' @examples
get_fcs_panel <- function(fcs_raw){
  tech_patterns <-list(computational_tech = c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE", "<NA>"),
                       marker_tech = c("_BC", "BCKG", "DNA"))
  panel <- as.data.frame(flowCore::pData(flowCore::parameters(fcs_raw[[1]]))[,c("name", "desc")])
  rownames(panel) <- NULL
  panel <- panel[sapply(panel$name, function(x) !any(sapply(tech_patterns$computational_tech, function(y) grepl(y,x)))),]
  panel <- panel[sapply(panel$desc, function(x) !any(sapply(tech_patterns$computational_tech, function(y) grepl(y,x)))),]
  panel <- panel[!is.na(panel$desc),]
  panel$antigen <- sapply(strsplit(panel$desc, "_"), function(x) x[length(x)])
  #panel$antigen <- gsub(" \\(v)", "", panel$antigen)
  panel$marker_class <- "type"
  panel$marker_class[sapply(panel$desc, function(x) any(sapply(tech_patterns$marker_tech, function(y) grepl(y,x))))] <- "state"
  #panel <- as.data.frame(apply(panel, c(1,2), function(x) gsub(" ", "_", x)))
  return(panel)
}

### test
#panel <- get_fcs_panel(fcs_raw)

##### Create use_marker
#' filtering out technical markers
#' @description the function formed list of markers which have class
#' mentioned as "type" in the panel. Other markers with class "state"
#' are considered as technical signals
#'
#' @param panel dataframe with marker panel data
#' @export
#' @return the list value of which are marker names and names of the
#' list element are names of antibody, relevant to the markers
#'
#'
#' @examples
get_use_marker <- function(panel){
  use_markers <- as.character(panel$name[panel$marker_class == "type"])
  names(use_markers) <- panel$antigen[panel$marker_class == "type"]
  #names(use_markers) <- gsub("_(v)", "", names(use_markers), fixed = T)          ### It should be fixed (problem with "_(v)")
  return(use_markers)
}

### test
#use_markers <- get_use_marker(panel)

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
#' @param use_marker
#' @export
#' @return flowSet object with transform data by asinh function and divided
#' by cofactor
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"
#'
#' @examples
asinh_transformation <- function(fcs_raw, cofactor = 5, use_marker = NULL){
  if(is.null(use_marker)){
    computational_tech <- c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE")
    markers <- colnames(fcs_raw[[1]])[sapply(colnames(fcs_raw[[1]]), function(x) !any(sapply(computational_tech, function(y) grepl(y,x))))]
  }
  markers <- use_marker
  fcs_asinh <- flowCore::fsApply(fcs_raw, function(x, cf = cofactor, mk = markers){
    exprs(x)[,mk] <- asinh(exprs(x)[,mk] / cf)
    return(x)})
  return(fcs_asinh)
}

### test
#fcs_raw <- asinh_transformation(fcs_raw, 5)

##### Transformation to a from 0 to 1 variable and removing outliers
#' Transformation to a from 0 to 1 variable and removing outliers
#'
#' @param fcs_raw
#' @param quantile
#' @param use_marker
#'
#' @return
#' @export
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"
#' @importFrom matrixStats colQuantiles
#'
#' @examples
outlier_by_quantile_transformation <- function(fcs_raw, quantile = 0.01, use_marker = NULL){
  if(is.null(use_marker)){
    computational_tech <- c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE")
    markers <- colnames(fcs_raw[[1]])[sapply(colnames(fcs_raw[[1]]), function(x) !any(sapply(computational_tech, function(y) grepl(y,x))))]
  }
  markers <- use_marker
  fcs_outlier_by_quantile <- flowCore::fsApply(fcs_raw, function(x, ql = quantile, mk = markers){
    rng <- matrixStats::colQuantiles(exprs(x)[,mk], probs = c(ql, 1-ql))
    expr_data <- t((t(exprs(x)[,mk]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    expr_data[expr_data < 0] <- 0
    expr_data[expr_data > 1] <- 1
    exprs(x)[,mk] <- expr_data
    return(x)
  })
  return(fcs_outlier_by_quantile)
}

### test
#fcs_raw <- outlier_by_quantile_transformation(fcs_raw, 0.01)

##### Extract cell number
#' Extract cell number
#'
#' @param fcs_raw
#'
#' @return
#' @importFrom flowCore sampleNames "sampleNames<-" fsApply
#' @export
#' @examples
get_cell_number <- function(fcs_raw){
  cell_number <- data.frame(smpl = flowCore::sampleNames(fcs_raw), cell_nmbr = flowCore::fsApply(fcs_raw, nrow))
  colnames(cell_number) <- c('smpl', 'cell_nmbr')
  return(cell_number)
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
#' @importFrom flowCore fsApply sampleNames "sampleNames<-" exprs "exprs<-"
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#'
#' @examples
scatter_plot_data_prep <- function(fcs_raw, use_markers, sampling_size = 0.5, method = "tSNE",
                         perplexity = 30, theta = 0.5, max_iter = 1000, size_fuse = 5000){
  #sampling_size <- as.integer(sampling_size/length(fcs_raw))
  expr <- flowCore::fsApply(fcs_raw[,use_markers], flowCore::exprs)
  sample_ids <- rep(flowCore::sampleNames(fcs_raw), flowCore::fsApply(fcs_raw, nrow))

  ## Find and skip duplicates
  dups <- which(!duplicated(expr[, use_markers]))

  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids)

  ## How many cells to downsample per-sample
  #tsne_ncells <- pmin(table(sample_ids), sampling_size)             ################ Number of cells to ploting
  tsne_ncells <- as.integer((table(sample_ids) + 1) * sampling_size)
  if(!is.null(size_fuse) & (sum(tsne_ncells) > size_fuse)){tsne_ncells <- as.integer((tsne_ncells/sum(tsne_ncells))*size_fuse)}
  names(tsne_ncells) <- names(table(sample_ids))

  ## Get subsampled indices
  set.seed(1234)
  tsne_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
    intersect(s, dups)
  })

  tsne_inds <- unlist(tsne_inds)
  tsne_expr <- expr[tsne_inds, use_markers]

  if(method == "tSNE"){
    ##### Run t-SNE
    set.seed(1234)
    tsne_result <- Rtsne::Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE,
                         perplexity = perplexity, theta = theta, max_iter = max_iter)
    #tsne_out <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2])
    tsne_out <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2], expr[tsne_inds, use_markers])
    colnames(tsne_out) <- c("tSNE1", "tSNE2", use_markers)
  }

  if(method == "UMAP"){
    ##### Run UMAP
    umap_out <- umap::umap(tsne_expr)
    tsne_out <- data.frame(tSNE1 = umap_out$layout[, 1], tSNE2 = umap_out$layout[, 2], expr[tsne_inds, use_markers])
    colnames(tsne_out) <- c("tSNE1", "tSNE2", use_markers)
  }

  colnames(tsne_out)[match(use_markers, colnames(tsne_out))] <- names(use_markers)
  return(tsne_out)
}

#tSNE <- scatter_plot_data_prep(fcs_raw, use_markers, sampling_size = 0.1, method = "tSNE")
#ggplot(tSNE,  aes(x = tSNE1, y = tSNE2, color = tSNE[,names(use_markers)[10]])) +
#  geom_point(size = 0.2) +
#  labs(color = names(use_markers)[10])
