#' MA-plot from base means and log fold changes
#'
#' MA-plot from base means and log fold changes, in the ggplot2 framework, with
#' additional support to annotate genes if provided.
#'
#' The genes of interest are to be provided as gene symbols if a \code{symbol}
#' column is provided in \code{res_obj}, or else b< using  the identifiers specified
#' in the row names
#'
#' @param res_obj A \code{\link{DESeqResults}} object
#' @param FDR Numeric value, the significance level for thresholding adjusted p-values
#' @param point_alpha Alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param sig_color Color to use to mark differentially expressed genes. Defaults to red
#' @param annotation_obj A \code{data.frame} object, with row.names as gene
#' identifiers (e.g. ENSEMBL ids) and a column, \code{gene_name}, containing
#' e.g. HGNC-based gene symbols. Optional
#' @param hlines The y coordinate (in absolute value) where to draw horizontal lines,
#' optional
#' @param title A title for the plot, optional
#' @param xlab X axis label, defaults to "mean of normalized counts - log10 scale"
#' @param ylim Vector of two numeric values, Y axis limits to restrict the view
#' @param add_rug Logical, whether to add rug plots in the margins
#' @param intgenes Vector of genes of interest. Gene symbols if a \code{symbol}
#' column is provided in \code{res_obj}, or else the identifiers specified in the
#' row names
#' @param intgenes_color The color to use to mark the genes on the main plot.
#' @param labels_intgenes Logical, whether to add the gene identifiers/names close
#' to the marked plots
#' @param labels_repel Logical, whether to use \code{geom_text_repel} for placing the
#' labels on the features to mark
#'
#' @return An object created by \code{ggplot}
#' @export
#'
#' @examples
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design=~cell+dex)
#' # subsetting for quicker run, ignore the next two commands if regularly using the function 
#' gene_subset <- c(
#'   "ENSG00000103196",  # CRISPLD2
#'   "ENSG00000120129",  # DUSP1
#'   "ENSG00000163884",  # KLF15
#'   "ENSG00000179094",  # PER1
#'   rownames(dds_airway)[rep(c(rep(FALSE,99), TRUE), length.out=nrow(dds_airway))]) # 1% of ids
#' dds_airway <- dds_airway[gene_subset,]
#'                   
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#'
#' plot_ma(res_airway, FDR = 0.05, hlines = 1)
#'
#' plot_ma(res_airway, FDR = 0.1,
#'         intgenes = c("ENSG00000103196",  # CRISPLD2
#'                      "ENSG00000120129",  # DUSP1
#'                      "ENSG00000163884",  # KLF15
#'                      "ENSG00000179094")  # PER1
#'        )
#'
#'
plot_ma <- function(res_obj,
                    FDR = 0.05,
                    point_alpha = 0.2,
                    sig_color = 'red',
                    annotation_obj = NULL, # TODO: add a check, if not available skip this part
                    hlines = NULL,
                    title = NULL,
                    xlab = "mean of normalized counts - log10 scale",
                    ylim = NULL,
                    add_rug = TRUE,
                    intgenes = NULL,
                    intgenes_color = "steelblue",
                    labels_intgenes = TRUE,
                    labels_repel = TRUE) {
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

  if(!is.null(ylim))
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
      if(labels_repel) {
        p <- p + geom_text_repel(data = df_intgenes,aes_string("logmean", "log2FoldChange",label="myids"),
                           color = intgenes_color, size=5)
      } else {
        p <- p + geom_text(data = df_intgenes,aes_string("logmean", "log2FoldChange",label="myids"),
                         color = intgenes_color, size=5,hjust=0.25, vjust=-0.75)
      }    
    }
  }

  if(add_rug)
    p <- p + geom_rug(alpha = 0.3)

  p <- p + theme_bw()
  p
}

