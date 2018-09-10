#' Extract functional terms enriched in the DE genes, based on goseq
#'
#' A wrapper for extracting functional GO terms enriched in a list of (DE) genes,
#' based on the algorithm and the implementation in the goseq package
#'
#' Note: the feature length retrieval is based on the \code{\link{goseq}} function,
#' and requires that the corresponding TxDb packages are installed and available
#'
#' @param de.genes A vector of (differentially expressed) genes
#' @param assayed.genes A vector of background genes, e.g. all (expressed) genes
#' in the assays
#' @param genome A string identifying the genome that genes refer to, as in the
#' \code{\link{goseq}} function
#' @param id A string identifying the gene identifier used by genes, as in the
#' \code{\link{goseq}} function
#' @param testCats A vector specifying which categories to test for over representation amongst DE genes - can be any combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG"
#' @param FDR_GO_cutoff Numeric value for subsetting the results
#' @param nTop Number of categories to extract, and optionally process for adding
#' genes to the respective
#' @param orgDbPkg Character string, named as the \code{org.XX.eg.db}
#' package which should be available in Bioconductor
#' @param addGeneToTerms Logical, whether to add a column with all genes annotated
#' to each GO term
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#' @export
#'
#' @examples
#'
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design=~cell+dex)
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#'
#' res_subset <- deseqresult2DEgenes(res_airway)[1:100,]
#' myde <- res_subset$id
#' myassayed <- rownames(res_airway)
#'
#' \dontrun{
#' mygo <- goseqTable(myde,
#'                    myassayed,
#'                    testCats = "GO:BP",
#'                    addGeneToTerms = FALSE)
#' head(mygo)
#' }
#'
goseqTable <- function(de.genes,                  # Differentially expressed genes
                       assayed.genes,             # background genes, normally = rownames(cds) or filtering to genes
                       #  with at least 1 read - could also be ls(org.Mm.egGO)
                       genome = "hg38",
                       id= "ensGene",
                       testCats=c("GO:BP","GO:MF","GO:CC"),
                       FDR_GO_cutoff = 1,
                       nTop = 200,
                       orgDbPkg = "org.Hs.eg.db",
                       # testKegg=TRUE,
                       # keggObject=mapPathwayToName("mmu"), # need the dedicated function!!
                       # writeOutput=FALSE,
                       addGeneToTerms=TRUE # ,
                       # outputFiles_goseq="",outputFiles_goseq_kegg=""
                       ## TODO TODO: bring back in action the function
                       ## add genes annotated to each term
                       ## do it by default only for bp?
                       ## tests at the beginning to see if the whole thing is feasible?
)
{
  #  library(goseq)
  #  library(GO.db)
  gene.vector <- as.integer(assayed.genes %in% de.genes)
  names(gene.vector) <- assayed.genes
  fdr <- FDR_GO_cutoff

  pwf <- nullp(DEgenes=gene.vector,genome = genome, id= id,plot.fit=FALSE)

  goseq_out <-  goseq(pwf,genome=genome,id=id,test.cats=testCats)



  goseq_out$p.adj <- p.adjust(goseq_out$over_represented_pvalue,method="BH")

  # to reduce the load for adding the genes
  goseq_out <- goseq_out[seq_len(nTop),]

  if(addGeneToTerms) {
    # for adding the gene ids/names...
    gene2cat = getgo(de.genes,genome=genome,id=id,fetch.cats= testCats)
    names(gene2cat) = de.genes

    reversemap <- function(map) # as in goseq
    {
      tmp <- unlist(map, use.names = FALSE)
      names(tmp) <- rep(names(map), times = as.numeric(summary(map)[, 1]))
      return(split(names(tmp), as.vector(tmp)))
    }

    cat2gene = reversemap(gene2cat)
    # one list per GO term
    goseq_out$genes <- sapply(goseq_out$category, function(x) cat2gene[[x]])

    # TODO: replace identifiers/annotaions!!!
    ## and also TODO: do this only if genes are not already symbols
    goseq_out$genesymbols <- sapply(goseq_out$genes, function(x) sort(AnnotationDbi::mapIds(get(orgDbPkg),keys=x,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")))
    # coerce to char
    goseq_out$genes <- unlist(lapply(goseq_out$genes,function(arg) paste(arg,collapse=",")))
    # coerce to char
    goseq_out$genesymbols <- unlist(lapply(goseq_out$genesymbols,function(arg) paste(arg,collapse=",")))

  }

  return(goseq_out)
}

