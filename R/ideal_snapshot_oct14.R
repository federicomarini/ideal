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
#'
#'
#'
#'
ideal<- function(
  res_obj = NULL,
  dds_obj = NULL,
  annotation_obj = NULL,
  countmatrix = NULL,
  expdesign = NULL){

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("ideal requires 'shiny'. Please install it using
         install.packages('shiny')")
  }

  # get modes and themes for the ace editor
  modes <- shinyAce::getAceModes()
  themes <- shinyAce::getAceThemes()

  ## upload max 300mb files - can be changed if necessary
  options(shiny.maxRequestSize=300*1024^2)

  ## ------------------------------------------------------------------ ##
  ##                          Define UI                                 ##
  ## ------------------------------------------------------------------ ##

  source(file.path("R/ideal_ui.R"),local = TRUE)
  source(file.path("R/ideal_server.R"), local= TRUE)





  # components defined in separated .R files
  shinyApp(ui = ideal_ui, server = ideal_server)

}












############################# helper funcs #################################

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
                         color = intgenes_color, size=5,hjust = 0, nudge_x = 0.2)
    }

  }


  p <- p + theme_bw()

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
    class = "footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        hr(),
        "ideal is a project developed by Federico Marini in the Bioinformatics division of the ",
        tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
        ". ",br(),
        "Development of the ideal package is on ",
        tags$a(href="https://github.com/federicomarini/ideal", "GitHub")
      )
    )
  )
}
