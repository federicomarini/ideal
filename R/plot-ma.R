
plot_ma <- function(object, which_beta, which_model = 'full',
                    sig_level = 0.10,
                    point_alpha = 0.2,
                    sig_color = 'red',
                    annotation_obj = NULL, # TODO: add a check, if not available skip this part
                    hlines = NULL

) {
  mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
  mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
  mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
  # mama$yesorno <- ifelse(mama$isDE,"yes","no")
  mama$yesorno <- ifelse(mama$isDE,"red","black")

  p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))

  if(!is.null(hlines)) {
    p <- p + geom_hline(aes(yintercept = hlines), col = "lightblue", alpha = 0.4) + geom_hline(aes(yintercept = -hlines), col = "lightblue", alpha = 0.4)
  }


  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-3,3)
  # p <- p + xlab('mean( log( counts + 0.5 ) )')
  # p <- p + ylab(paste0('beta: ', which_beta))

  #   if (!is.null(highlight)) {
  #     suppressWarnings({
  #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
  #     })
  #     if (nrow(highlight) > 0) {
  #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
  #     } else {
  #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
  #     }
  #   }


  p <- p + theme_bw()


  p
}



plot_ma_highlight <- function(object, which_beta, which_model = 'full',
                              sig_level = 0.10,
                              point_alpha = 0.2,
                              sig_color = "red",
                              intgenes = NULL,
                              intgenes_color = "steelblue",
                              labels_intgenes = TRUE,
                              annotation_obj = annotation_obj # SAME AS IN FUNC ABOVE

) {
  mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
  mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
  mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
  # mama$yesorno <- ifelse(mama$isDE,"yes","no")
  mama$yesorno <- ifelse(mama$isDE,"red","black")

  p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))
  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-4,4) # or even remove the y limit
  # p <- p + xlab('mean( log( counts + 0.5 ) )')
  # p <- p + ylab(paste0('beta: ', which_beta))

  #   if (!is.null(highlight)) {
  #     suppressWarnings({
  #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
  #     })
  #     if (nrow(highlight) > 0) {
  #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
  #     } else {
  #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
  #     }
  #   }

  if(!is.null(intgenes)){
    # check that genes are in the results object:
    # will check the id as rownames or an explicit symbol column


    # now here for the symbol
    res_df <- as.data.frame(object)
    res_df$logmean <- log10(res_df$baseMean)
    df_intgenes <- res_df[res_df$symbol %in% intgenes,]

    p <- p + geom_point(data = df_intgenes,aes(logmean, log2FoldChange), color = intgenes_color, size = 4)

    if(labels_intgenes) {
      p <- p + geom_text(data = df_intgenes,aes(logmean, log2FoldChange,label=symbol),
                         color = intgenes_color, size=5,hjust=0.25, vjust=-0.75)
    }

  }


  p <- p + theme_bw()

  p
}

