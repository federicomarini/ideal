### commented out on 24.11.2016

# lengthData <- fc_judith$annotation$Length
# names(lengthData) <- fc_judith$annotation$GeneID
# assayed.genes <- rownames(fc_judith$counts)
#
# goseqOutput <- function(de.genes,                  # Differentially expressed genes
#                         assayed.genes,             # background genes, normally = rownames(cds) or filtering to genes
#                         #  with at least 1 read - could also be ls(org.Mm.egGO)
#                         lengthData = lengthData,
#                         genome = "hg38",
#                         id= "ensGene",
#                         plotFit = TRUE,
#                         testCats=c("GO:BP","GO:MF","GO:CC"),
#                         FDR_GO_cutoff = 1,
#                         testKegg=TRUE,
#                         keggObject=mapPathwayToName("mmu"), # need the dedicated function!!
#                         writeOutput=FALSE,
#                         addGeneToTerms=FALSE,
#                         outputFiles_goseq="",outputFiles_goseq_kegg=""
#                         ## TODO TODO: bring back in action the function
#                         ## add genes annotated to each term
#                         ## do it by default only for bp?
#                         ## tests at the beginning to see if the whole thing is feasible?
# )
# {
#   library(goseq)
#   library(GO.db)
#   gene.vector <- as.integer(assayed.genes %in% de.genes)
#   names(gene.vector) <- assayed.genes
#   fdr <- FDR_GO_cutoff
#   #   head(gene.vector)
#   # selecting just the ones in the gene universe
#
#   # # lengthData, taken from the annotation
#   # lengthData <- fc_kleinert_paired$annotation$Length
#   # names(lengthData) <- fc_kleinert_paired$annotation$GeneID
#
#   mm10_lengthData_exp <- lengthData[names(lengthData) %in% assayed.genes]
#   pwf <- nullp(DEgenes=gene.vector,bias.data=mm10_lengthData_exp,plot.fit=plotFit)
#
#   #   GO.wallBP <- goseq(pwf,genome="mm10",id="geneSymbol",test.cats="GO:BP")
#   #   GO.wallMF <- goseq(pwf,genome="mm10",id="geneSymbol",test.cats="GO:MF")
#   #   GO.wallCC <- goseq(pwf,genome="mm10",id="geneSymbol",test.cats="GO:CC")
#   if(testKegg) {
#     mouseKEGGs <- keggObject
#     GO.kegg <- goseq(pwf,genome=genome,id=id,test.cats="KEGG") # test.cats="KEGG"?
#     GO.kegg$p.adj <- p.adjust(GO.kegg$over_represented_pvalue,method="BH")
#
#     GO.kegg <- GO.kegg[GO.kegg$p.adj < fdr,]
#     GO.kegg$pathway_name <- sapply(GO.kegg$category,function(keggID){mouseKEGGs$pathway_name[rownames(mouseKEGGs) %in% keggID]})
#
#
#     if(addGeneToTerms) {
#       gene2cat_kegg = getgo(de.genes,genome=genome,id="ensGene",fetch.cats= c("KEGG"))
#       names(gene2cat_kegg) = de.genes
#
#       cat2gene_kegg = goseq:::reversemapping(gene2cat_kegg)
#       GO.kegg$genes <- sapply(GO.kegg$category, function(x) cat2gene_kegg[[x]])
#       GO.kegg$genes <- sapply(GO.kegg$genes, function(x)
#       {
#         if(length(x)==0) return(NULL)
#         sort(AnnotationDbi::mapIds(org.Mm.eg.db,keys=x,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first"))
#       })
#
#       GO.kegg$genes <- unlist(lapply(GO.kegg$genes,function(arg) paste(arg,collapse=",")))
#     }
#
#     if(writeOutput) write.table(KEGGsign,file=outputFiles_goseq_kegg,sep="\t",quote=F,col.names=T,row.names=F)
#
#   }
#
#   GO.wall.all <- goseq(pwf,genome=genome,id=id,test.cats=testCats)
#
#   GO.wall.all <- GO.wall.all[1:200,]
#
#   # eventual selection for the DEgenes
#   #   gene2cat = getgo(de.genes),genome="mm10",id="geneSymbol",fetch.cats= c("GO:CC","GO:BP", "GO:MF"))
#   #   names(gene2cat) = de.genes
#   #   cat2gene = goseq:::reversemapping(gene2cat)
#
#   gene2cat = getgo(de.genes,genome=genome,id="ensGene",fetch.cats= c("GO:CC","GO:BP", "GO:MF"))
#   names(gene2cat) = de.genes
#   cat2gene = goseq:::reversemapping(gene2cat)
#
#   # cat2gene[["GO:1903047"]]
#
#   if(addGeneToTerms) {
#     GO.wall.all$genes <- sapply(GO.wall.all$category, function(x) cat2gene[[x]])
#     GO.wall.all$genesymbols <- sapply(GO.wall.all$genes, function(x) sort(AnnotationDbi::mapIds(org.Mm.eg.db,keys=x,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")))
#
#   }
#
#   GO.wall.all$p.adj <- p.adjust(GO.wall.all$over_represented_pvalue,method="BH")
#   GOsign <- GO.wall.all[GO.wall.all$p.adj < fdr,]
#
#
#   if(writeOutput) able(GOsign,file=outputFiles_goseq,sep="\t",quote=F,col.names=T,row.names=F)
#
#
#   if(testKegg) {
#     return(list(GO=GOsign, KEGG=GO.kegg))
#   } else {
#     return(list(GO=GOsign))
#   }
#
# }
