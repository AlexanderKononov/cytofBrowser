context("Checking the main pipeline from upload data to clustering")

test_that(desc = "Check the main pipeline from upload data to clustering", code = {
  test_files <- c("KPC2_DNA_viSNE_CD45-_stroma.fcs", "KPC3_DNA_viSNE_CD45-_stroma.fcs",
                  "KPC4_DNA_viSNE_CD45-_stroma.fcs", "KPC5_DNA_viSNE_CD45-_stroma.fcs",
                  "KPC6_DNA_viSNE_CD45-_stroma.fcs", "KPC7_DNA_viSNE_CD45-_stroma.fcs",
                  "KPC8_DNA_viSNE_CD45-_stroma.fcs")
  test_files <- system.file("extdata",test_files,package = "cytofBrowser")
  #cytofBrowser::get_fcs_metadata(test_files)
})

