#' Plot normalized counts for a gene
#'
#' Plot for normalized counts of a single gene, with jittered points superimposed
#' on the boxplot
#'
#' Note: this function relies on the [DESeq2::plotCounts()] function of DESeq2,
#' therefore pseudocounts of 0.5 are added to each point
#'
#' @param dds A [DESeq2::DESeqDataSet()] object.
#' @param gene A character, specifying the name of the gene to plot
#' @param intgroup Interesting groups: a character vector of
#' names in `colData(dds)` to use for grouping
#' @param annotation_obj A `data.frame` object, with `row.names` as gene
#' identifiers (e.g. ENSEMBL ids) and a column, `gene_name`, containing
#' e.g. HGNC-based gene symbols. Optional.
#' @param transform Logical value, corresponding whether to have log scale y-axis
#' or not. Defaults to TRUE.
#' @param labels_repel Logical value. Whether to use `ggrepel`'s functions to
#' place labels; defaults to TRUE.
#'
#' @return An object created by `ggplot`
#' 
#' @importFrom mosdef gene_plot
#' 
#' @export
#'
#' @examples
#' library("airway")
#' data("airway", package = "airway")
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'   colData = colData(airway),
#'   design = ~ cell + dex
#' )
#' ggplotCounts(dds_airway,
#'   gene = "ENSG00000103196", # CRISPLD2 in the original publication
#'   intgroup = "dex"
#' )
ggplotCounts <- function(dds, gene, intgroup = "condition", annotation_obj = NULL,
                         transform = TRUE, labels_repel = TRUE) {
  .Deprecated(old = "ggplotCounts", new = "mosdef::gene_plot", 
              msg = paste0(
                "Please use `mosdef::gene_plot()` in replacement of the `ggplotCounts()` function, ",
                "originally located in the ideal package. \nCheck the manual page for ",
                "`?mosdef::gene_plot()` to see the details on how to use it, e.g. ",
                "refer to the new parameter definition and naming"))
  
  p <- mosdef::gene_plot(
    de_container = dds,
    gene = gene,
    intgroup = intgroup,
    annotation_obj = annotation_obj,
    transform = transform,
    labels_repel = labels_repel
  )
    
  return(p)
}
