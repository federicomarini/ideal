
plot_ma <- function(res_obj,
                    FDR = 0.05,
                    point_alpha = 0.2,
                    sig_color = 'red',
                    annotation_obj = NULL, # TODO: add a check, if not available skip this part
                    hlines = NULL,
                    title = NULL,
                    xlab = "mean of normalized counts - log10 scale",
                    ylim = NULL) {
  ma_df <- data.frame(
    mean = res_obj$baseMean,
    lfc = res_obj$log2FoldChange,
    padj = res_obj$padj,
    isDE = ifelse(is.na(res_obj$padj), FALSE, res_obj$padj < FDR),
    ID = rownames(res_obj))

  ma_df <- ma_df[ma_df$mean > 0,]

  if(!is.null(annotation_obj))
    ma_df$genename <- annotation_obj$gene_name[match(ma_df$ID,rownames(annotation_obj))]

  ma_df$logmean <- log10(ma_df$mean) # TO ALLOW FOR BRUSHING!!
  # ma_df$DE <- ifelse(ma_df$isDE,"yes","no")
  ma_df$DE <- ifelse(ma_df$isDE,"red","black")

  p <- ggplot(ma_df, aes_string(x = "logmean", y = "lfc", colour = "DE"))

  if(!is.null(hlines)) {
    p <- p + geom_hline(aes(yintercept = hlines), col = "lightblue", alpha = 0.4) +
      geom_hline(aes(yintercept = -hlines), col = "lightblue", alpha = 0.4)
  }
  p <- p + geom_hline(aes(yintercept = 0), col = "red", alpha = 0.4)

  p <- p + xlab(xlab) + ylab("log fold change")

  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(
    name = paste0("FDR = ",FDR),
    values = c("black", sig_color),
    labels = c("nonDE","DE"))

  if(is.null(ylim))
    p <- p + coord_cartesian(ylim = c(-2,2))
  else
    p <- p + coord_cartesian(ylim = ylim)

  if(!is.null(title))
    p <- p + ggtitle(title)

  p <- p + theme_bw()

  p
}



plot_ma_highlight <- function(res_obj,
                              FDR = 0.05,
                              point_alpha = 0.2,
                              sig_color = 'red',
                              annotation_obj = NULL, # TODO: add a check, if not available skip this part
                              hlines = NULL,
                              title = NULL,
                              xlab = "mean of normalized counts - log10 scale",
                              ylim = NULL,
                              intgenes = NULL,
                              intgenes_color = "steelblue",
                              labels_intgenes = TRUE) {
  ma_df <- data.frame(
    mean = res_obj$baseMean,
    lfc = res_obj$log2FoldChange,
    padj = res_obj$padj,
    isDE = ifelse(is.na(res_obj$padj), FALSE, res_obj$padj < FDR),
    ID = rownames(res_obj))

  ma_df <- ma_df[ma_df$mean > 0,]

  if(!is.null(annotation_obj))
    ma_df$genename <- annotation_obj$gene_name[match(ma_df$ID,rownames(annotation_obj))]

  ma_df$logmean <- log10(ma_df$mean) # TO ALLOW FOR BRUSHING!!
  # ma_df$DE <- ifelse(ma_df$isDE,"yes","no")
  ma_df$DE <- ifelse(ma_df$isDE,"red","black")

  p <- ggplot(ma_df, aes_string(x = "logmean", y = "lfc", colour = "DE"))

  if(!is.null(hlines)) {
    p <- p + geom_hline(aes(yintercept = hlines), col = "lightblue", alpha = 0.4) +
      geom_hline(aes(yintercept = -hlines), col = "lightblue", alpha = 0.4)
  }
  p <- p + geom_hline(aes(yintercept = 0), col = "red", alpha = 0.4)

  p <- p + xlab(xlab) + ylab("log fold change")

  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(
    name = paste0("FDR = ",FDR),
    values = c("black", sig_color),
    labels = c("nonDE","DE"))

  if(is.null(ylim))
    p <- p + coord_cartesian(ylim = c(-2,2))
  else
    p <- p + coord_cartesian(ylim = ylim)

  if(!is.null(title))
    p <- p + ggtitle(title)

  if(!is.null(intgenes)){

    # now here for the symbol
    res_df <- as.data.frame(res_obj)
    res_df$logmean <- log10(res_df$baseMean)

    if("symbol" %in% colnames(res_df)){
      # use the gene names
      df_intgenes <- res_df[res_df$symbol %in% intgenes, ]
      df_intgenes$myids <- df_intgenes$symbol

    } else {
      # use whatever is there as id
      df_intgenes <- res_df[rownames(res_df) %in% intgenes, ]
      df_intgenes$myids <- rownames(df_intgenes)
    }

    # df_intgenes <- res_df[res_df$symbol %in% intgenes,]

    p <- p + geom_point(data = df_intgenes,aes_string("logmean", "log2FoldChange"), color = intgenes_color, size = 4)

    if(labels_intgenes) {
      p <- p + geom_text(data = df_intgenes,aes_string("logmean", "log2FoldChange",label="myids"),
                         color = intgenes_color, size=5,hjust=0.25, vjust=-0.75)
    }

  }


  p <- p + theme_bw()

  p
}

