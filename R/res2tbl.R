#' Generate a tidy table with the results of DESeq
#'
#' Generate a tidy table with the results of DESeq
#'
#' @param deseqresult A [DESeq2::DESeqResults()] object
#'
#' @return A "tidy" data.frame with all genes
#' @export
#' 
#' @importFrom mosdef deresult_to_df
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 100, m = 8, betaSD = 1)
#' dds <- DESeq2::DESeq(dds)
#' res <- DESeq2::results(dds)
#' deseqresult2tbl(res)
deseqresult2tbl <- function(deseqresult) {
  .Deprecated(old = "deseqresult2tbl", new = "mosdef::deresult_to_df", 
              msg = paste0(
                "Please use `mosdef::deresult_to_df()` in replacement of the `deseqresult2tbl()` function, ",
                "originally located in the ideal package. \nCheck the manual page for ",
                "`?mosdef::deresult_to_df()` to see the details on how to use it, e.g. ",
                "refer to the new parameter definition and naming"))
  
  res_de <- mosdef::deresult_to_df(deseqresult)
  
  return(res_de)
}




#' Generate a tidy table with the DE genes from the results of DESeq
#'
#' Generate a tidy table with the DE genes from the results of DESeq
#'
#' @param deseqresult A [DESeq2::DESeqResults()] object
#' @param FDR Numeric value, the significance level for thresholding adjusted p-values
#'
#' @return A "tidy" data.frame with only genes marked as differentially expressed
#' @export
#' 
#' @importFrom mosdef deresult_to_df
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 100, m = 8, betaSD = 2)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' deseqresult2DEgenes(res)
deseqresult2DEgenes <- function(deseqresult,
                                FDR = 0.05) {
  .Deprecated(old = "deseqresult2DEgenes", new = "mosdef::deresult_to_df", 
              msg = paste0(
                "Please use `mosdef::deresult_to_df()` in replacement of the `deseqresult2DEgenes()` function, ",
                "originally located in the ideal package. \nCheck the manual page for ",
                "`?mosdef::deresult_to_df()` to see the details on how to use it, e.g. ",
                "refer to the new parameter definition and naming"))
  
  res_de <- mosdef::deresult_to_df(deseqresult, FDR = FDR)
  
  return(res_de)
}
