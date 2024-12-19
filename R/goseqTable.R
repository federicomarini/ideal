#' Extract functional terms enriched in the DE genes, based on goseq
#'
#' A wrapper for extracting functional GO terms enriched in a list of (DE) genes,
#' based on the algorithm and the implementation in the goseq package
#'
#' Note: the feature length retrieval is based on the [goseq::goseq()] function,
#' and requires that the corresponding TxDb packages are installed and available
#'
#' @param de.genes A vector of (differentially expressed) genes
#' @param assayed.genes A vector of background genes, e.g. all (expressed) genes
#' in the assays
#' @param genome A string identifying the genome that genes refer to, as in the
#' [goseq::goseq()] function
#' @param id A string identifying the gene identifier used by genes, as in the
#' [goseq::goseq()] function
#' @param testCats A vector specifying which categories to test for over representation amongst DE genes - can be any combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG"
#' @param FDR_GO_cutoff Numeric value for subsetting the results
#' @param nTop Number of categories to extract, and optionally process for adding
#' genes to the respective terms
#' @param orgDbPkg Character string, named as the `org.XX.eg.db`
#' package which should be available in Bioconductor
#' @param addGeneToTerms Logical, whether to add a column with all genes annotated
#' to each GO term
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#' 
#' @export
#' 
#' @importFrom mosdef run_goseq
#'
#' @examples
#'
#' library("airway")
#' data("airway", package = "airway")
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'   colData = colData(airway),
#'   design = ~ cell + dex
#' )
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#'
#' res_subset <- mosdef::deresult_to_df(res_airway)[1:100, ]
#' myde <- res_subset$id
#' myassayed <- rownames(res_airway)
#' \dontrun{
#' mygo <- goseqTable(myde,
#'   myassayed,
#'   testCats = "GO:BP",
#'   addGeneToTerms = FALSE
#' )
#' head(mygo)
#' }
#'
goseqTable <- function(de.genes, # Differentially expressed genes
                       assayed.genes, # background genes, normally = rownames(cds) or filtering to genes
                       #  with at least 1 read - could also be ls(org.Mm.egGO)
                       genome = "hg38",
                       id = "ensGene",
                       testCats = c("GO:BP", "GO:MF", "GO:CC"),
                       FDR_GO_cutoff = 1,
                       nTop = 200,
                       orgDbPkg = "org.Hs.eg.db",
                       # testKegg=TRUE,
                       # keggObject=mapPathwayToName("mmu"), # need the dedicated function!!
                       # writeOutput=FALSE,
                       addGeneToTerms = TRUE # ,
                       # outputFiles_goseq="",outputFiles_goseq_kegg=""
                       ## TODO TODO: bring back in action the function
                       ## add genes annotated to each term
                       ## do it by default only for bp?
                       ## tests at the beginning to see if the whole thing is feasible?
) {
  
  .Deprecated(old = "goseqTable", new = "mosdef::run_goseq", 
              msg = paste0(
                "Please use `mosdef::run_goseq()` in replacement of the `goseqTable()` function, ",
                "originally located in the ideal package. \nCheck the manual page for ",
                "`?mosdef::run_goseq()` to see the details on how to use it, e.g. ",
                "refer to the new parameter definition and naming"))
  
  res_enrich <- mosdef::run_goseq(
    de_genes = de.genes, 
    bg_genes = assayed.genes,
    genome = genome,
    id = id, 
    testCats = testCats, 
    mapping = orgDbPkg, 
    add_gene_to_terms = addGeneToTerms
  )
    
  return(res_enrich)
}
