# ideal.R

#' Title
#'
#' @param dds
#' @param rlt
#' @param countmatrix
#' @param coldata
#' @param pca2go
#' @param annotation
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#'
ideal<- function(res_obj = NULL,
                 dds_obj = NULL,
                 annotation_obj = NULL,
                 countmatrix = NULL,
                 expdesign = NULL){

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("ideal requires 'shiny'. Please install it using
         install.packages('shiny')")
  }


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
  library(goseq)
  library(GO.db)
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
    cat2gene = goseq:::reversemapping(gene2cat)
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




library(stringr)
sepguesser <- function(f) {
  separators_list = c(",", "\t", ";"," ")
  rl = readLines(f, warn = FALSE)
  rl = rl[rl != ""] # allow last line to be empty
  sephits_min = sapply(separators_list, function(x) min(stringr::str_count(rl, x))) #minimal number of separators on a line
  sep = separators_list[which.max(sephits_min)]
  sep
}



# anno_df$entrez <- mapIds(org.Hs.eg.db,keys = anno_df$gene_id, keytype = "ENSEMBL", column = "ENTREZID")
# ggss <- rentrez::entrez_summary("gene", anno_df$entrez[1:150])

deseqresult2tbl <-
function(deseqresult) {
  library("dplyr")
  if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)
  deseqresult$id <- rownames(deseqresult)
  rownames(deseqresult) <- NULL
  deseqresult <- tbl_df(deseqresult)
  deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:padj)
  deseqresult %>% arrange(padj)
}

deseqresult2DEgenes <-
function(deseqresult,FDR=0.05) {
  library("dplyr")
  if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)
  deseqresult$id <- rownames(deseqresult)
  rownames(deseqresult) <- NULL
  deseqresult <- tbl_df(deseqresult)
  if("symbol" %in% names(deseqresult))
    deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:symbol)
  else
    deseqresult <- dplyr::select(deseqresult, id, baseMean, log2FoldChange:padj)
  tmp <- deseqresult %>% arrange(padj)
  res <- tmp[!(is.na(tmp$padj)) & tmp$padj <= FDR,]
  res
}



plot_ma <- function(object, which_beta, which_model = 'full',
                    sig_level = 0.10,
                    point_alpha = 0.2,
                    sig_color = 'red',
                    annotation_obj = NULL, # TODO: add a check, if not available skip this part
                    hlines = NULL

) {
  mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
  mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
  mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
  # mama$yesorno <- ifelse(mama$isDE,"yes","no")
  mama$yesorno <- ifelse(mama$isDE,"red","black")

  p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))

  if(!is.null(hlines)) {
    p <- p + geom_hline(aes(yintercept = hlines), col = "lightblue", alpha = 0.4) + geom_hline(aes(yintercept = -hlines), col = "lightblue", alpha = 0.4)
  }


  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-3,3)
  # p <- p + xlab('mean( log( counts + 0.5 ) )')
  # p <- p + ylab(paste0('beta: ', which_beta))

  #   if (!is.null(highlight)) {
  #     suppressWarnings({
  #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
  #     })
  #     if (nrow(highlight) > 0) {
  #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
  #     } else {
  #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
  #     }
  #   }


  p <- p + theme_bw()


  p
}



plot_ma_highlight <- function(object, which_beta, which_model = 'full',
                    sig_level = 0.10,
                    point_alpha = 0.2,
                    sig_color = "red",
                    intgenes = NULL,
                    intgenes_color = "steelblue",
                    labels_intgenes = TRUE,
                    annotation_obj = annotation_obj # SAME AS IN FUNC ABOVE

) {
  mama <- data.frame(mean=object$baseMean,lfc=object$log2FoldChange,padj = object$padj,isDE= ifelse(is.na(object$padj), FALSE, object$padj < sig_level),ID=rownames(object))
  mama$genename <- annotation_obj$gene_name[match(mama$ID,rownames(annotation_obj))]
  mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
  # mama$yesorno <- ifelse(mama$isDE,"yes","no")
  mama$yesorno <- ifelse(mama$isDE,"red","black")

  p <- ggplot(mama, aes(logmean, lfc,colour=yesorno))
  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color)) + ylim(-4,4) # or even remove the y limit
  # p <- p + xlab('mean( log( counts + 0.5 ) )')
  # p <- p + ylab(paste0('beta: ', which_beta))

  #   if (!is.null(highlight)) {
  #     suppressWarnings({
  #       highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
  #     })
  #     if (nrow(highlight) > 0) {
  #       p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
  #     } else {
  #       warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
  #     }
  #   }

  if(!is.null(intgenes)){
    # check that genes are in the results object:
    # will check the id as rownames or an explicit symbol column


    # now here for the symbol
    res_df <- as.data.frame(object)
    res_df$logmean <- log10(res_df$baseMean)
    df_intgenes <- res_df[res_df$symbol %in% intgenes,]

    p <- p + geom_point(data = df_intgenes,aes(logmean, log2FoldChange), color = intgenes_color, size = 4)

    if(labels_intgenes) {
      p <- p + geom_text(data = df_intgenes,aes(logmean, log2FoldChange,label=symbol),
                         color = intgenes_color, size=5,hjust=0.25, vjust=-0.75)
    }

  }


  p <- p + theme_bw()

  p
}



plot_volcano <- function(object,
                         FDR = 0.1,
                         ...) {

  mydf <- as.data.frame(object)
  mydf$id <- rownames(mydf)
  mydf$isDE <- ifelse(is.na(object$padj), FALSE, object$padj < FDR)

  p <- ggplot(mydf, aes(x = log2FoldChange, y = -log10(pvalue))) + geom_point(aes(color = isDE), alpha = 0.4)

  p <- p + ylim(0,20)
  p <- p + theme_bw() + scale_color_manual(values=c("black", "red"))

  p <- p + geom_vline(aes(xintercept = 1), col = "lightblue", alpha = 0.6, size = 1.5) +
    geom_vline(aes(xintercept = -1), col = "lightblue", alpha = 0.6, size = 1.5)

  p
}






ggplotCounts <- function(dds,gene,intgroup="condition",annotation_obj=annotation_obj,...){
  df <- plotCounts(dds,gene,intgroup,returnData = T,...)
  df$sampleID <- rownames(df)
  genesymbol <- annotation_obj$gene_name[match(gene,annotation_obj$gene_id)]
  jittered_df <- df
  # jittered_df$conditionj <- jitter(as.numeric(factor(jittered_df$condition)))
  jittered_df$countj <- jitter(jittered_df$count)

  onlyfactors <- df[,match(intgroup,colnames(df))]
  df$plotby <- interaction(onlyfactors)


  p <- df %>%
    ggplot(aes(x=plotby,y=count,col=plotby)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_text(data = jittered_df,aes(x=conditionj,y=countj,label=sampleID)) +
    geom_text(aes(label=sampleID),hjust=-.1,vjust=0) +
    labs(title=paste0("Normalized counts for ",genesymbol," - ",gene)) +
    scale_x_discrete(name="") +
    geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) +
    scale_color_discrete(name="Experimental\nconditions") +
    scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.1,NA)) +
    theme_bw()

  p

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


createLinkGO <- function(val) {
  sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank" class="btn btn-primary">%s</a>',val,val)
}

createLinkENS  <- function(val, species="Mus_musculus") {
  paste0('<a href="http://www.ensembl.org/',species,'/Gene/Summary?g=',val,'" target="_blank" class="btn btn-primary">',val,'</a>')
}

createLinkGeneSymbol <- function(val) {
  # possibilities:
  # ncbi
  # genecards
  paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',val,'[sym]" target="_blank" class="btn btn-primary">',val,'</a>')
}


getGeneInfos <- function(obj, annopkg, idtype) {
  # obj is a dds object...
  ids <- rownames(obj)

  mydf <- mapIds(eval(parse(text=annopkg)),keys=ids,column = "GENENAME",keytype = idtype)
  mydf_2 <- AnnotationDbi::select(eval(parse(text=annopkg)),keys=ids,column = "GENENAME",keytype = idtype)

  return(mydf_2)
}








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
