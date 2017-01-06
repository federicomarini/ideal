
plot_volcano <- function(object,
                         FDR = 0.1,
                         ...) {

  mydf <- as.data.frame(object)
  mydf$id <- rownames(mydf)
  mydf$isDE <- ifelse(is.na(object$padj), FALSE, object$padj < FDR)

  p <- ggplot(mydf, aes(x = log2FoldChange, y = -log10(pvalue))) + geom_point(aes(color = isDE), alpha = 0.4)

  p <- p + ylim(0,20)
  p <- p + theme_bw() + scale_color_manual(values=c("black", "red"))

  p <- p + geom_vline(aes(xintercept = 1), col = "lightblue", alpha = 0.6, size = 1.5) +
    geom_vline(aes(xintercept = -1), col = "lightblue", alpha = 0.6, size = 1.5)

  p
}

