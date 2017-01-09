


deseqresult2tbl <- function(deseqresult) {
  # library("dplyr")
  if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)

  deseqresult <- cbind(rownames(deseqresult),deseqresult)
  names(deseqresult)[1] <- "id"

  dplyr::arrange_(deseqresult,"padj")
}

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

