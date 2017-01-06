
# anno_df$entrez <- mapIds(org.Hs.eg.db,keys = anno_df$gene_id, keytype = "ENSEMBL", column = "ENTREZID")
# ggss <- rentrez::entrez_summary("gene", anno_df$entrez[1:150])

deseqresult2tbl <-
  function(deseqresult) {
    # library("dplyr")
    if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
    deseqresult <- as.data.frame(deseqresult)
    deseqresult$id <- rownames(deseqresult)
    rownames(deseqresult) <- NULL
    deseqresult <- dplyr::tbl_df(deseqresult)
    deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:padj)
    dplyr::arrange(deseqresult,padj)
  }

deseqresult2DEgenes <-
  function(deseqresult,FDR=0.05) {
    # library("dplyr")
    if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
    deseqresult <- as.data.frame(deseqresult)
    deseqresult$id <- rownames(deseqresult)
    rownames(deseqresult) <- NULL
    deseqresult <- dplyr::tbl_df(deseqresult)
    if("symbol" %in% names(deseqresult))
      deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:symbol)
    else
      deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:padj)
    tmp <- dplyr::arrange(deseqresult,padj)
    res <- tmp[!(is.na(tmp$padj)) & tmp$padj <= FDR,]
    res
  }

