library("ideal")

context("Check that shiny app is generated")

dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
# vst <- vst(dds)
# cm <- counts(dds)
# cd <- colData(dds)

test_that("Shiny app is generated", {
  expect_is(ideal(), "shiny.appobj")
  expect_is(ideal(dds_obj = dds), "shiny.appobj")
})
