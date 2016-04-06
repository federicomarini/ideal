# # example for running...
# library(airway)
# library(DESeq2)
# library(ggplot2)
# library(pheatmap)
# data(airway)
# dds_airway <- DESeqDataSet(airway, design= ~ cell + dex)
#
#
# rld_airway <- rlogTransformation(dds_airway)
# # constructing the annotation object
# anno_df <- data.frame(gene_id = rownames(dds_airway),
#                      stringsAsFactors=FALSE)
# library("AnnotationDbi")
# library("org.Hs.eg.db")
# anno_df$gene_name <- mapIds(org.Hs.eg.db,
#                            keys=anno_df$gene_id,
#                            column="SYMBOL",
#                            keytype="ENSEMBL",
#                            multiVals="first")
# rownames(anno_df) <- anno_df$gene_id
#
#
# bg_ids <- rownames(dds_airway)[rowSums(counts(dds_airway)) > 0]
# # library(topGO)
# # pca2go_airway <- pca2go(rld_airway,
# #                        annotation = anno_df,
# #                        organism = "Hs",
# #                        ensToGeneSymbol = TRUE,
# #                        background_genes = bg_ids)
#
# dds_airway <- DESeq(dds_airway)
# res_airway <- results(dds_airway)
# library("AnnotationDbi")
# library("org.Hs.eg.db")
# res_airway$symbol <- mapIds(org.Hs.eg.db,
#                             keys=row.names(res_airway),
#                             column="SYMBOL",
#                             keytype="ENSEMBL",
#                             multiVals="first")
# res_airway$entrez <- mapIds(org.Hs.eg.db,
#                             keys=row.names(res_airway),
#                             column="ENTREZID",
#                             keytype="ENSEMBL",
#                             multiVals="first")
# resOrdered <- as.data.frame(res_airway[order(res_airway$padj),])
# de_df <- resOrdered[resOrdered$padj < .05 & !is.na(resOrdered$padj),]
# de_symbols <- de_df$symbol
# bg_ids <- rownames(dds_airway)[rowSums(counts(dds_airway)) > 0]
# bg_symbols <- mapIds(org.Hs.eg.db,
#                      keys=bg_ids,
#                      column="SYMBOL",
#                      keytype="ENSEMBL",
#                      multiVals="first")
# # library(topGO)
# #
# # topgoDE_airway <- topGOtable(de_symbols, bg_symbols,
# #                              ontology = "BP",
# #                              mapping = "org.Hs.eg.db",
# #                              geneID = "symbol")
#
# ma_SUPALIVE(res_airway,dds_airway)
#
#
#
#
