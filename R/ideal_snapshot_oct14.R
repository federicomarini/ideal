# ideal.R

#' ideal: Interactive Differential Expression Analysis
#'
#' ideal makes differential expression analysis interactive, easy and reproducible.
#' This function launches the main application included in the package.
#'
#' @param dds_obj A \code{\link{DESeqDataSet}} object. If not provided, then a
#' \code{countmatrix} and a \code{expdesign} need to be provided. If none of
#' the above is provided, it is possible to upload the data during the
#' execution of the Shiny App
#' @param res_obj  A \code{\link{DESeqResults}} object. If not provided, it can
#' be computed during the execution of the application
#' @param annotation_obj A \code{data.frame} object, with row.names as gene
#' identifiers (e.g. ENSEMBL ids) and a column, \code{gene_name}, containing
#' e.g. HGNC-based gene symbols. If not provided, it can be constructed during
#' the execution via the org.eg.XX.db packages - these need to be installed
#' @param countmatrix A count matrix, with genes as rows and samples as columns.
#' If not provided, it is possible to upload the data during the execution of
#' the Shiny App
#' @param expdesign A \code{data.frame} containing the info on the covariates
#' of each sample. If not provided, it is possible to upload the data during the
#' execution of the Shiny App
#'
#' @return A Shiny App is launched for interactive data exploration and
#' differential expression analysis
#'
#' @export
#'
#' @examples
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
#' cm <- counts(dds)
#' cd <- colData(dds)
#'
#' # with the well known airway package...
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design=~cell+dex)
#' \dontrun{
#'
#'
#' ideal()
#' ideal(dds)
#' ideal(dds_airway)
#'
#' dds_airway <- DESeq(dds_airway)
#' res_airway <- results(dds_airway)
#' ideal(dds_airway, res_airway)
#' }
#'
ideal<- function(dds_obj = NULL,
                 res_obj = NULL,
                 annotation_obj = NULL,
                 countmatrix = NULL,
                 expdesign = NULL){

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("ideal requires 'shiny'. Please install it using
         install.packages('shiny')")
  }


  # create environment for storing inputs and values
  ## i need the assignment like this to export it up one level - i.e. "globally"
  ideal_env <<- new.env(parent = emptyenv())

  ## upload max 300mb files - can be changed if necessary
  options(shiny.maxRequestSize=300*1024^2)
  #
  options(shiny.launch.browser = TRUE)

  ## ------------------------------------------------------------------ ##
  ##                          Define UI                                 ##
  ## ------------------------------------------------------------------ ##

  source(file.path("R/ideal_ui.R"),local = TRUE)
  source(file.path("R/ideal_server.R"), local= TRUE)



  # components defined in separated .R files
  shinyApp(ui = ideal_ui, server = ideal_server)

}












############################# helper funcs #################################




# library(stringr)
sepguesser <- function(f) {
  separators_list = c(",", "\t", ";"," ")
  rl = readLines(f, warn = FALSE)
  rl = rl[rl != ""] # allow last line to be empty
  sephits_min = sapply(separators_list, function(x) min(stringr::str_count(rl, x))) #minimal number of separators on a line
  sep = separators_list[which.max(sephits_min)]
  sep
}










# combineTogether <- function(normCounts,resuTable,anns) {
#   combinedCountsAndRes <- inner_join(resuTable,normCounts,by="id")
#   anns2 <- anns[match(combinedCountsAndRes$id, anns[, 1]), ]
#   combinedCountsAndRes$Description <- anns2$description
#   return(combinedCountsAndRes)
# }
#
#
# combine_resucounts <- function(normCounts,resuTable) {
#   combinedCountsAndRes <- inner_join(resuTable,normCounts,by="id")
#   # anns2 <- anns[match(combinedCountsAndRes$id, anns[, 1]), ]
#   # combinedCountsAndRes$Description <- anns2$description
#   return(combinedCountsAndRes)
# }
#


# getGeneInfos <- function(obj, annopkg, idtype) {
#   # obj is a dds object...
#   ids <- rownames(obj)
#
#   mydf <- mapIds(eval(parse(text=annopkg)),keys=ids,column = "GENENAME",keytype = idtype)
#   mydf_2 <- AnnotationDbi::select(eval(parse(text=annopkg)),keys=ids,column = "GENENAME",keytype = idtype)
#
#   return(mydf_2)
# }
#
#






footer <- function(){
  tags$div(
    class = "panel-footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        # hr(),
        "ideal is a project developed by Federico Marini in the Bioinformatics division of the ",
        tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
        "- Institute for Medical Biostatistics, Epidemiology and Informatics",br(),
        "License: ",tags$a(href="https://opensource.org/licenses/MIT","MIT"), br(),

        "Development of the ideal package is on ",
        tags$a(href="https://github.com/federicomarini/ideal", "GitHub")
      )
    )
  )
}
