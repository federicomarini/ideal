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
#' @importFrom S4Vectors DataFrame
#' @import ggplot2
#' @import shiny
#' @importFrom DT datatable
#' @import shinydashboard
#' @importFrom AnnotationDbi mapIds select
#' @import shinyAce
#' @import BiocParallel
#' @import knitr
#' @import rmarkdown
#' @import plyr
#' @importFrom dplyr inner_join tbl_df
#' @importMethodsFrom GOstats hyperGTest summary
#' @import GO.db
#' @importFrom UpSetR upset fromList
#' @import goseq
#' @import pcaExplorer
#' @import rentrez
#' @importFrom pheatmap pheatmap
#' @import d3heatmap
#' @importFrom limma goana topGO
#' @import topGO
#' @import rintrojs
#' @importFrom shinyBS bsTooltip
#' @importFrom grDevices dev.off pdf
#' @import methods
#'
#' @author
#' Federico Marini \email{marinif@@uni-mainz.de}, 2016-2017
#'
#' Maintainer: Federico Marini \email{marinif@@uni-mainz.de}
#' @name ideal-pkg
#' @docType package
NULL
