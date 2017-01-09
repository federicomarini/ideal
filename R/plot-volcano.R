
plot_volcano <- function(res_obj,
                         FDR = 0.05,
                         ylim_up = NULL,
                         vlines = NULL,
                         title = NULL,
                         ...) {

  mydf <- as.data.frame(res_obj)
  mydf$id <- rownames(mydf)
  mydf$isDE <- ifelse(is.na(res_obj$padj), FALSE, res_obj$padj < FDR)

  mydf <- mydf[mydf$baseMean > 0,]

  p <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(pvalue)")) + geom_point(aes_string(color = "isDE"), alpha = 0.4)

  if(!is.null(ylim_up))
    p <- p + coord_cartesian(ylim = c(0, ylim_up))
  else
    p <- p + coord_cartesian(ylim = c(0,20))

  if(!is.null(title))
    p <- p + ggtitle(title)

  p <- p + theme_bw() +
    scale_colour_manual(
      name = paste0("FDR = ",FDR),
      values = c("black", "red"),
      labels = c("nonDE","DE"))

  p <- p + geom_vline(aes(xintercept = 1), col = "lightblue", alpha = 0.6, size = 1.5) +
    geom_vline(aes(xintercept = -1), col = "lightblue", alpha = 0.6, size = 1.5)

  p
}

