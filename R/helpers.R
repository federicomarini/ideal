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

#' Make an educated guess on the separator character
#'
#' This function tries to guess which separator was used in a text delimited file
#'
#' @param file The name of the file which the data are to be read from
#' @param sep_list A vector containing the candidates for being identified as
#' separators. Defaults to \code{c(",", "\t", ";"," ")}
#'
#' @return A character value, corresponding to the guessed separator. One of ","
#' (comma), "\\t" (tab), ";" (semicolon)," " (whitespace)
#' @export
#'
#' @examples
#' sepguesser(system.file("extdata/design_commas.txt",package = "ideal"))
#' sepguesser(system.file("extdata/design_semicolons.txt",package = "ideal"))
#' sepguesser(system.file("extdata/design_spaces.txt",package = "ideal"))
#' mysep <- sepguesser(system.file("extdata/design_tabs.txt",package = "ideal"))
#'
#' # to be used for reading in the same file, without having to specify the sep
#'
sepguesser <- function(file, sep_list = c(",", "\t", ";"," ")) {
  separators_list = sep_list
  rl = readLines(file, warn = FALSE)
  rl = rl[rl != ""] # allow last line to be empty
  sephits_min = sapply(separators_list, function(x) min(stringr::str_count(rl, x))) #minimal number of separators on all lines
  sep = separators_list[which.max(sephits_min)]
  sep
}




# this function is a (long running) one-liner to install
# - all org.XX.eg.db packages
# and/or
# - all TxDb packages
setup_idealannotation <- function(orgdb_pkgs = TRUE,
                        TxDb_pkgs = TRUE) {

  all_orgdb <- c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
                 "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
                 "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
                 "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
                 "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
                 "org.Xl.eg.db")

  all_txdb <- c("TxDb.Athaliana.BioMart.plantsmart22",  "TxDb.Athaliana.BioMart.plantsmart25",
                "TxDb.Athaliana.BioMart.plantsmart28",  "TxDb.Btaurus.UCSC.bosTau8.refGene",
                "TxDb.Celegans.UCSC.ce11.refGene",  "TxDb.Celegans.UCSC.ce6.ensGene",
                "TxDb.Cfamiliaris.UCSC.canFam3.refGene",  "TxDb.Dmelanogaster.UCSC.dm3.ensGene",
                "TxDb.Dmelanogaster.UCSC.dm6.ensGene",  "TxDb.Drerio.UCSC.danRer10.refGene",
                "TxDb.Ggallus.UCSC.galGal4.refGene",  "TxDb.Hsapiens.UCSC.hg18.knownGene",
                "TxDb.Hsapiens.UCSC.hg19.knownGene",  "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts",
                "TxDb.Hsapiens.UCSC.hg38.knownGene",  "TxDb.Mmulatta.UCSC.rheMac3.refGene",
                "TxDb.Mmulatta.UCSC.rheMac8.refGene",  "TxDb.Mmusculus.UCSC.mm10.ensGene",
                "TxDb.Mmusculus.UCSC.mm10.knownGene",  "TxDb.Mmusculus.UCSC.mm9.knownGene",
                "TxDb.Ptroglodytes.UCSC.panTro4.refGene",  "TxDb.Rnorvegicus.UCSC.rn4.ensGene",
                "TxDb.Rnorvegicus.UCSC.rn5.refGene",  "TxDb.Rnorvegicus.UCSC.rn6.refGene",
                "TxDb.Scerevisiae.UCSC.sacCer2.sgdGene",  "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
                "TxDb.Sscrofa.UCSC.susScr3.refGen")

  if(orgdb_pkgs)
    BiocInstaller::biocLite(all_orgdb,suppressUpdates = TRUE)

  if(TxDb_pkgs)
    BiocInstaller::biocLite(all_txdb,suppressUpdates = TRUE)

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

