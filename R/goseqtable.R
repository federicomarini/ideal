
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
  goseq_out <- goseq_out[1:nTop,]

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

