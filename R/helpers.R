# helpers.R





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

