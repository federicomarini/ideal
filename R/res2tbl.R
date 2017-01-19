#' Generate a tidy table with the results of DESeq
#'
#' Generate a tidy table with the results of DESeq
#'
#' @param deseqresult A \code{\link{DESeqResults}} object
#'
#' @return A "tidy" data.frame with all genes
#' @export
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8, betaSD = 1)
#' dds <- DESeq2::DESeq(dds)
#' res <- DESeq2::results(dds)
#' deseqresult2tbl(res)
#'
deseqresult2tbl <- function(deseqresult) {
  # library("dplyr")
  if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)

  deseqresult <- cbind(rownames(deseqresult),deseqresult)
  names(deseqresult)[1] <- "id"

  dplyr::arrange_(deseqresult,"padj")
}




#' Generate a tidy table with the DE genes from the results of DESeq
#'
#' Generate a tidy table with the DE genes from the results of DESeq
#'
#' @param deseqresult A \code{\link{DESeqResults}} object
#' @param FDR Numeric value, the significance level for thresholding adjusted p-values
#'
#' @return A "tidy" data.frame with only genes marked as differentially expressed
#' @export
#'
#' @examples
#'
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8, betaSD = 2)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' deseqresult2DEgenes(res)
#'
deseqresult2DEgenes <- function(deseqresult,
                                FDR=0.05) {
  # library("dplyr")
  if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)

  deseqresult <- cbind(rownames(deseqresult),deseqresult)
  names(deseqresult)[1] <- "id"

  # deseqresult$id <- rownames(deseqresult)
  # rownames(deseqresult) <- NULL
  # deseqresult <- dplyr::tbl_df(deseqresult)
  # if("symbol" %in% names(deseqresult))
  #   deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:symbol)
  # else
  #   deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:padj)
  tmp <- dplyr::arrange_(deseqresult,"padj")
  res <- tmp[!(is.na(tmp$padj)) & tmp$padj <= FDR,]
  res
}

