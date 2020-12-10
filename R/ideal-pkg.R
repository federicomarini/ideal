#' ideal: Interactive Differential Expression Analysis
#'
#' ideal makes differential expression analysis interactive, easy and reproducible.
#' The analysis of RNA-seq datasets is guided by the Shiny app as main component of
#' the package, which also provides a wide set of functions to efficiently extract
#' information from the existing data. The app can be also deployed on a Shiny
#' server, to allow its usage without any installation on the user's side.
#'
#' ideal makes differential expression analysis interactive, easy and reproducible.
#' The analysis of RNA-seq datasets is guided by the Shiny app as main component of
#' the package, which also provides a wide set of functions to efficiently extract
#' information from the existing data. The app can be also deployed on a Shiny
#' server, to allow its usage without any installation on the user's side.
#'
#' @import DESeq2
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame List
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @import shiny
#' @importFrom DT datatable
#' @import shinydashboard
#' @importFrom AnnotationDbi mapIds select keytypes
#' @importFrom shinyAce aceAutocomplete aceEditor getAceModes getAceThemes
#' updateAceEditor
#' @import BiocParallel
#' @import knitr
#' @import rmarkdown
#' @importFrom dplyr inner_join tbl_df filter mutate arrange last
#' @importMethodsFrom GOstats hyperGTest summary
#' @import GO.db
#' @importFrom UpSetR upset fromList
#' @importFrom goseq getgo goseq nullp
#' @import pcaExplorer
#' @importFrom gplots venn
#' @importFrom IHW ihw
#' @importFrom rentrez entrez_summary
#' @importFrom pheatmap pheatmap
#' @importFrom heatmaply heatmaply
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom limma goana topGO
#' @import topGO
#' @importFrom stringr str_count
#' @importFrom rintrojs introjs introjsUI
#' @importFrom shinyBS bsTooltip bsCollapse bsCollapsePanel
#' @importFrom ggrepel geom_text_repel
#' @importFrom base64enc dataURI
#' @importFrom grDevices dev.off pdf axisTicks
#' @import methods
#'
#' @author
#' Federico Marini \email{marinif@@uni-mainz.de}, 2016-2017
#'
#' Maintainer: Federico Marini \email{marinif@@uni-mainz.de}
#' @name ideal-pkg
#' @docType package
NULL

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription("ideal", fields = "Version")
  msg <- paste0("Welcome to ideal v", pkgVersion, "\n\n")
  citation <- paste0(
    "If you use ideal in your work, please cite:\n\n",
    "Federico Marini, Jan Linke, Harald Binder\n",
    "ideal: an R/Bioconductor package for Interactive Differential Expression Analysis\n",
    "BMC Bioinformatics, 2020 - https://doi.org/10.1186/s12859-020-03819-5\n"
  )
  packageStartupMessage(paste0(msg, citation))
}
