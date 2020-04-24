context("Checking the main pipeline from upload data to clustering")

test_that(desc = "Check the main pipeline from upload data to clustering", code = {
  test_files <- c("KPC1_stroma.fcs", "KPC2_stroma.fcs", "KPC3_stroma.fcs", "KPC4_stroma.fcs",
                  "KPC5_stroma.fcs", "KPC6_stroma.fcs", "KPC7_stroma.fcs")
  test_files <- system.file("extdata",test_files,package = "cytofBrowser")

  md <- cytofBrowser::get_fcs_metadata(test_files)
  testthat::expect_equal(dim(md), c(7,3), info = "Error in function: get_fcs_metadata")

  fcs_raw <- cytofBrowser::get_fcs_raw(md)
  testthat::expect_equal(length(fcs_raw), 7, info = "Error in function: get_fcs_raw")

  panel <- cytofBrowser::get_fcs_panel(fcs_raw)
  testthat::expect_equal(dim(panel), c(52,4), info = "Error in function: get_fcs_panel")
  testthat::expect_equal(all(panel$marker_class %in% c("type", "state")), TRUE,
               info = "Error in function: get_fcs_panel. problem with marker_class")

  use_markers <- cytofBrowser::get_use_marker(panel)
  testthat::expect_equal(all(c(use_markers %in% panel$name, names(use_markers) %in% panel$antigen)), TRUE,
               info = "Error in function: get_use_marker")

  fcs_asinh <- cytofBrowser::asinh_transformation(fcs_raw)
  testthat::expect_equal(length(fcs_asinh), length(fcs_raw), info = "Error in function: asinh_transformation")
  testthat::expect_equal(dim(fcs_asinh[[1]]), dim(fcs_raw[[1]]), info = "Error in function: asinh_transformation")
  fcs_raw <- fcs_asinh
  rm(fcs_asinh)

  fcs_outlier_by_quantile <- cytofBrowser::outlier_by_quantile_transformation(fcs_raw)
  testthat::expect_equal(length(fcs_outlier_by_quantile), length(fcs_raw),
               info = "Error in function: outlier_by_quantile_transformation")
  testthat::expect_equal(dim(fcs_outlier_by_quantile[[1]]), dim(fcs_raw[[1]]),
                         info = "Error in function: outlier_by_quantile_transformation")
  fcs_raw <- fcs_outlier_by_quantile
  rm(fcs_outlier_by_quantile)

  cell_number <- cytofBrowser::get_cell_number(fcs_raw)
  testthat::expect_equal(dim(cell_number), c(7,2), info = "Error in function: get_cell_number")

  tsne_out <- cytofBrowser::scatter_plot_data_prep(fcs_raw, use_markers, size_fuse = 300)
  testthat::expect_equal((dim(tsne_out)[1] <= 300) & (dim(tsne_out)[2] == length(use_markers)+2), TRUE,
                         info = "Error in function: scatter_plot_data_prep")
  tsne_out <- cytofBrowser::scatter_plot_data_prep(fcs_raw, use_markers, method = "UMAP", size_fuse = 300)
  testthat::expect_equal((dim(tsne_out)[1] <= 300) & (dim(tsne_out)[2] == length(use_markers)+2), TRUE,
                         info = "Error in function: scatter_plot_data_prep, UMAP")

  som <- cytofBrowser::get_som(fcs_raw, use_markers)
  testthat::expect_equal(dim(som$map$mapping)[[1]] >= 1, TRUE, info = "Error in function: get_som")

  mc <- cytofBrowser::get_consensusClust(som, maxK = 6)
  testthat::expect_equal(length(unique(mc[[6]]$consensusClass)) >= 2, TRUE,
                         info = "Error in function: get_consensusClust")

  cell_clustering_list <- cytofBrowser::get_cluster_annotation(fcs_raw, som, mc, 6)
  testthat::expect_equal(length(cell_clustering_list), 7, info = "Error in function: get_cluster_annotation")
  testthat::expect_equal(length(unique(unlist(cell_clustering_list))) , 6,
                         info = "Error in function: get_cluster_annotation")
  cell_clustering <- cytofBrowser::get_cell_clustering_vector(som, mc, 6)
  testthat::expect_equal(length(unique(cell_clustering)), 6,
                         info = "Error in function: get_cell_clustering_vector")

  cluster_euclidean_distance <- cytofBrowser::get_euclid_dist(fcs_raw, use_markers, cell_clustering)
  testthat::expect_equal(dim(cluster_euclidean_distance)[1] >= 1, TRUE,
                         info = "Error in function: get_euclid_dist")

  tsne_inds <- cytofBrowser::get_inds_subset(fcs_raw, size_fuse = 300)
  testthat::expect_equal(length(tsne_inds) <= 300, TRUE, info = "Error in function: get_inds_subset")

  umap_df <- cytofBrowser::get_UMAP_dataframe(fcs_raw, use_markers, use_markers, tsne_inds, cell_clustering)
  testthat::expect_equal(dim(umap_df)[1] <= 300, TRUE, info = "Error in function: get_UMAP_dataframe")

  umap_df <- cytofBrowser::get_UMAP_dataframe(fcs_raw, use_markers, use_markers,
                                              tsne_inds, cell_clustering, method = "tSNE")
  testthat::expect_equal(dim(umap_df)[1] <= 300, TRUE, info = "Error in function: get_UMAP_dataframe, tSNE")
})

