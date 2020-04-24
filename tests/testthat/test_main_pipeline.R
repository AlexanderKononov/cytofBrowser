context("Checking the main pipeline from upload data to clustering")

test_that(desc = "Check the main pipeline from upload data to clustering", code = {
  test_files <- c("KPC2_DNA_viSNE_CD45-_stroma.fcs", "KPC3_DNA_viSNE_CD45-_stroma.fcs",
                  "KPC4_DNA_viSNE_CD45-_stroma.fcs", "KPC5_DNA_viSNE_CD45-_stroma.fcs",
                  "KPC6_DNA_viSNE_CD45-_stroma.fcs", "KPC7_DNA_viSNE_CD45-_stroma.fcs",
                  "KPC8_DNA_viSNE_CD45-_stroma.fcs")
  test_files <- system.file("extdata",test_files,package = "cytofBrowser")

  md <- cytofBrowser::get_fcs_metadata(test_files)
  expect_equal(dim(md), c(7,3), info = "Error in function: get_fcs_metadata")

  fcs_raw <- cytofBrowser::get_fcs_raw(md)
  expect_equal(length(fcs_raw), 7, info = "Error in function: get_fcs_raw")

  panel <- cytofBrowser::get_fcs_panel(fcs_raw)
  expect_equal(dim(panel), c(52,4), info = "Error in function: get_fcs_panel")
  expect_equal(all(panel$marker_class %in% c("type", "state")), TRUE,
               info = "Error in function: get_fcs_panel. problem with marker_class")

  use_markers <- cytofBrowser::get_use_marker(panel)
  expect_equal(all(c(use_markers %in% panel$name, names(use_markers) %in% panel$antigen)), TRUE,
               info = "Error in function: get_use_marker")

  fcs_asinh <- cytofBrowser::asinh_transformation(fcs_raw)
  expect_equal(length(fcs_asinh), length(fcs_raw), info = "Error in function: asinh_transformation")
  expect_equal(dim(fcs_asinh[[1]]), dim(fcs_raw[[1]]), info = "Error in function: asinh_transformation")
  fcs_raw <- fcs_asinh
  rm(fcs_asinh)

  fcs_outlier_by_quantile <- cytofBrowser::outlier_by_quantile_transformation(fcs_raw)
  expect_equal(length(fcs_outlier_by_quantile), length(fcs_raw),
               info = "Error in function: outlier_by_quantile_transformation")

})

