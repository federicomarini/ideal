#' wrapup_for_iSEE
#'
#' Combine data from a typical DESeq2 run
#'
#' Combines the DESeqDataSet input and DESeqResults into a SummarizedExperiment
#' object, which can be readily explored with iSEE.
#'
#' A typical usage would be after running the DESeq2 pipeline as specified in
#' one of the workflows which include this package, e.g. in the context of the
#' ideal package.
#'
#' @param dds A \code{\link{DESeqDataSet}} object.
#' @param res A \code{\link{DESeqResults}} object.
#'
#' @return A SummarizedExperiment object, with raw counts, normalized counts, and
#' variance-stabilizing transformed counts in the assay slots; and with colData
#' and rowData extracted from the corresponding input parameters
#'
#' @export
#'
#' @examples
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 10000, m = 8)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' se <- wrapup_for_iSEE(dds, res)
#' # library(iSEE)
#' # iSEE(se)
#' \dontrun{
#' # or with the well known airway package...
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'   colData = colData(airway),
#'   design = ~ cell + dex
#' )
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#' se_airway <- wrapup_for_iSEE(dds_airway, res_airway)
#' # iSEE(se_airway)
#' }
wrapup_for_iSEE <- function(dds, res) {
  # sanity checks on the objects
  stopifnot(all(rownames(dds) == rownames(res)))

  # dds to vst
  vst <- vst(dds)

  # initialize the container
  se <- SummarizedExperiment(
    assays = List(
      counts = counts(dds),
      normcounts = counts(dds, normalized = TRUE),
      vst_counts = assay(vst)
    )
  )

  # adding colData, taken directly from the DESeqDataSet object
  colData(se) <- colData(dds)

  # extract contrast info
  this_contrast <- sub(".*p-value: (.*)", "\\1", mcols(res, use.names = TRUE)["pvalue", "description"])

  # getting the rowData from the dds itself
  rdd <- rowData(dds)

  # modifying in advance the DESeqResults object
  res$log10_baseMean <- log10(res$baseMean)
  res$log10_pvalue <- -log10(res$pvalue)
  # and for the rowData
  rdd$log10_dispersion <- log10(rdd$dispersion)

  # adding rowData to se
  rowData(se)[[paste0("DESeq2_", gsub(" ", "_", this_contrast))]] <- res

  # merging in the existing rowData slot
  rowData(se) <- cbind(rowData(se), rdd)

  return(se)
}
